#include "detail.hpp"

namespace fxp_physics {
using namespace detail;
namespace {
bool concave_contact_shape(ShapeType type) {
    return type == ShapeType::mesh || type == ShapeType::heightfield;
}
bool rounded_contact_shape(ShapeType type) {
    return type == ShapeType::sphere || type == ShapeType::capsule;
}
} // namespace

// Execution scratch only: each job owns one collider's ordered candidate list
// and output. Retaining capacities avoids thousands of allocations per substep.
struct World::Impl::ContactData {
    struct Candidate {
        std::uint32_t index{};
        std::size_t cached_pair = std::numeric_limits<std::size_t>::max();
    };
    struct Job {
        std::vector<std::uint32_t> query;
        std::vector<Candidate> candidates;
        std::vector<Contact> contacts;
        std::vector<Pair> pairs;
    };
    std::vector<Job> jobs;
    std::vector<std::uint8_t> unchanged, matched, geometry_needed;
    std::vector<std::pair<BodyId, BodyId>> excluded;
    std::vector<Contact> previous;
    std::vector<Pair> previous_pairs;
};

void World::Impl::find_contacts() {
    rebuild();
    if (!contact_data)
        contact_data = std::make_shared<ContactData>();
    auto &data = *contact_data;
    const auto inactive = [](const Body &body) {
        return body.sleeping || body.desc.type == BodyType::static_body ||
               (body.desc.type == BodyType::kinematic && zero(body.desc.linear_velocity) &&
                zero(body.desc.angular_velocity));
    };
    data.unchanged.assign(colliders.size(), 0);
    // An empty pose cache also means invalidated history after removing the
    // final collider. Let the empty rebuild clear contacts and pair statistics.
    bool all_unchanged = !colliders.empty() && contact_poses.size() == colliders.size();
    for (std::size_t i = 0; i < colliders.size(); ++i) {
        data.unchanged[i] = contact_poses.size() == colliders.size() &&
                            contact_poses[i].position == proxies[i].pose.position &&
                            contact_poses[i].rotation.get_internal_vector() ==
                                proxies[i].pose.rotation.get_internal_vector() &&
                            inactive(bodies[proxies[i].body]);
        all_unchanged = all_unchanged && data.unchanged[i];
    }
    if (all_unchanged) {
        // Keep sensors and stay events, including the final-substep pair set.
        pairs = contact_pairs;
        for (auto &contact : contacts)
            contact.bounce = 0;
        return;
    }
    data.previous.swap(contacts);
    data.previous_pairs.swap(contact_pairs);
    const auto &previous = data.previous;
    const auto &previous_pairs = data.previous_pairs;
    contacts.clear();
    pairs.clear();
    data.matched.assign(previous.size(), 0);
    data.excluded.clear();
    for (const auto &spring : springs)
        if (!spring.desc.collide_connected && spring.desc.body_b)
            data.excluded.emplace_back((std::min)(spring.desc.body_a, spring.desc.body_b),
                                       (std::max)(spring.desc.body_a, spring.desc.body_b));
    std::sort(data.excluded.begin(), data.excluded.end());
    data.jobs.resize(colliders.size());

    // Every tree query is read-only and returns indices in the same order.
    // Resolve unchanged pairs here so their geometry need not be prepared.
    parallel_for(colliders.size(), colliders.size(), [&](std::size_t i) {
        auto &job = data.jobs[i];
        job.contacts.clear();
        job.pairs.clear();
        job.candidates.clear();
        const auto &ca = colliders[i];
        const auto &pa = proxies[i];
        const auto &ba = bodies[pa.body];
        tree.query(pa.bounds.expanded(config.contact_margin), job.query);
        for (const auto j : job.query) {
            if (j <= i)
                continue;
            const auto &cb = colliders[j];
            const auto &bb = bodies[proxies[j].body];
            if (ca.desc.body == cb.desc.body ||
                (ba.desc.type == BodyType::static_body && bb.desc.type == BodyType::static_body))
                continue;
            if (!(layers[ca.desc.layer] & (std::uint32_t(1) << cb.desc.layer)))
                continue;
            if (std::binary_search(data.excluded.begin(), data.excluded.end(),
                                   std::make_pair((std::min)(ba.id, bb.id), (std::max)(ba.id, bb.id))))
                continue;
            ContactData::Candidate candidate;
            candidate.index = j;
            if (data.unchanged[i] && data.unchanged[j]) {
                const auto key = std::make_pair(ca.id, cb.id);
                const auto old = std::lower_bound(
                    previous_pairs.begin(), previous_pairs.end(), key,
                    [](const Pair &pair, const auto &value) { return std::make_pair(pair.a, pair.b) < value; });
                if (old != previous_pairs.end() && old->a == ca.id && old->b == cb.id)
                    candidate.cached_pair = static_cast<std::size_t>(old - previous_pairs.begin());
            }
            job.candidates.push_back(candidate);
        }
    });

    // Gather only participating convex shapes. Sparse worlds therefore do not
    // pay to cook every collider. Lazy edge-direction caches must be filled
    // before shared geometry is read by the narrowphase workers.
    data.geometry_needed.assign(colliders.size(), 0);
    std::size_t pair_count = 0;
    for (std::size_t i = 0; i < data.jobs.size(); ++i) {
        pair_count += data.jobs[i].candidates.size();
        for (const auto candidate : data.jobs[i].candidates) {
            if (candidate.cached_pair != std::numeric_limits<std::size_t>::max())
                continue;
            const auto j = candidate.index;
            const auto type_a = colliders[i].desc.shape->type(), type_b = colliders[j].desc.shape->type();
            const bool directions = !rounded_contact_shape(type_a) && !rounded_contact_shape(type_b) &&
                                    !(type_a == ShapeType::box && type_b == ShapeType::box);
            if (!concave_contact_shape(type_a))
                data.geometry_needed[i] |= directions ? 3 : 1;
            if (!concave_contact_shape(type_b))
                data.geometry_needed[j] |= directions ? 3 : 1;
        }
    }
    parallel_for(colliders.size(), pair_count, [&](std::size_t i) {
        if (!data.geometry_needed[i])
            return;
        auto &cooked = proxies[i].geometry;
        if (!cooked)
            cooked = geometry(colliders[i].desc.shape->data(), proxies[i].pose);
        if (data.geometry_needed[i] & 2)
            (void)cooked->directions();
    });

    // Convex geometry is immutable throughout this phase. Concave triangle
    // geometry and its lazy directions are private to the current pair.
    const auto geometry_for = [&](std::size_t index, std::uint32_t triangle,
                                  Geometry &temporary) -> const Geometry & {
        const auto &shape = colliders[index].desc.shape->data();
        if (concave_contact_shape(shape.type)) {
            temporary = geometry(shape, proxies[index].pose, triangle);
            return temporary;
        }
        return *proxies[index].geometry;
    };
    parallel_for(colliders.size(), pair_count, [&](std::size_t i) {
        auto &job = data.jobs[i];
        const auto &ca = colliders[i];
        const auto &pa = proxies[i];
        const auto &ba = bodies[pa.body];
        for (const auto candidate : job.candidates) {
            const auto j = candidate.index;
            const auto &cb = colliders[j];
            const auto &pb = proxies[j];
            const auto &bb = bodies[pb.body];
            const auto key = std::make_pair(ca.id, cb.id);
            const auto begin = std::lower_bound(
                previous.begin(), previous.end(), key,
                [](const Contact &old, const auto &value) { return std::make_pair(old.a, old.b) < value; });
            if (candidate.cached_pair != std::numeric_limits<std::size_t>::max()) {
                job.pairs.push_back(previous_pairs[candidate.cached_pair]);
                for (auto old = begin; old != previous.end() && old->a == ca.id && old->b == cb.id; ++old) {
                    job.contacts.push_back(*old);
                    job.contacts.back().bounce = 0;
                }
                continue;
            }
            const auto ta =
                triangles(ca.desc.shape->data(), pa.pose, pb.bounds.expanded(config.contact_margin));
            const auto tb =
                triangles(cb.desc.shape->data(), pb.pose, pa.bounds.expanded(config.contact_margin));
            bool paired = false;
            for (auto ai : ta)
                for (auto bi : tb) {
                    Geometry temporary_a, temporary_b;
                    const auto &ga = geometry_for(i, ai, temporary_a);
                    const auto &gb = geometry_for(j, bi, temporary_b);
                    const bool sensor = ca.desc.sensor || cb.desc.sensor;
                    const Real margin = sensor ? 0_r : config.contact_margin;
                    const auto sep = separate(ga, gb, margin);
                    if (sep.distance > margin)
                        continue;
                    const auto points = manifold(ga, gb, sep, config.contact_margin);
                    if (!paired) {
                        job.pairs.push_back(
                            {ca.id, cb.id, ba.id, bb.id, sensor, points.front().point, sep.normal});
                        paired = true;
                    }
                    if (sensor ||
                        (ba.desc.type != BodyType::dynamic && bb.desc.type != BodyType::dynamic))
                        continue;
                    for (const auto &point : points) {
                        Contact c;
                        c.a = ca.id;
                        c.b = cb.id;
                        c.triangle_a = ai;
                        c.triangle_b = bi;
                        c.normal = sep.normal;
                        c.separation = point.separation;
                        c.local_a = untransform(ba.desc.pose, point.point);
                        c.local_b = untransform(bb.desc.pose, point.point);
                        const Real speed =
                            dot(point_velocity(bb, point.point) - point_velocity(ba, point.point), c.normal);
                        const Real restitution =
                            (std::max)(ca.desc.material.restitution, cb.desc.material.restitution);
                        c.bounce = speed < -config.restitution_threshold ? -speed * restitution : 0_r;
                        Real nearest = 0.03_r;
                        std::size_t selected = previous.size();
                        for (auto it = begin; it != previous.end() && it->a == c.a && it->b == c.b; ++it) {
                            const auto k = static_cast<std::size_t>(it - previous.begin());
                            if (data.matched[k] || it->triangle_a != ai || it->triangle_b != bi ||
                                (it->normal != c.normal && dot(it->normal, c.normal) < 0.95_r))
                                continue;
                            const Vec3 delta_a = it->local_a - c.local_a, delta_b = it->local_b - c.local_b;
                            // Both rounded lengths are nonnegative. Reject a
                            // distant corner before computing either square root.
                            if (length_greater_than(delta_a, nearest) || length_greater_than(delta_b, nearest))
                                continue;
                            const Real distance = length(delta_a) + length(delta_b);
                            if (distance < nearest) {
                                nearest = distance;
                                selected = k;
                            }
                        }
                        if (selected != previous.size()) {
                            // Each previous collider pair belongs exclusively
                            // to one job. Bytes avoid vector<bool>'s shared words.
                            data.matched[selected] = 1;
                            c.normal_impulse = previous[selected].normal_impulse;
                            c.tangent_impulse = previous[selected].tangent_impulse;
                            if (previous[selected].normal != c.normal)
                                c.tangent_impulse -= c.normal * dot(c.tangent_impulse, c.normal);
                        }
                        job.contacts.push_back(c);
                    }
                }
        }
    });

    // Scheduling never determines simulation order. The old nested i/j loop
    // is reconstructed exactly, including triangle and contact-point ordering.
    std::size_t contact_count = 0, touching_pairs = 0;
    for (const auto &job : data.jobs) {
        contact_count += job.contacts.size();
        touching_pairs += job.pairs.size();
    }
    contacts.reserve(contact_count);
    pairs.reserve(touching_pairs);
    for (const auto &job : data.jobs) {
        contacts.insert(contacts.end(), job.contacts.begin(), job.contacts.end());
        pairs.insert(pairs.end(), job.pairs.begin(), job.pairs.end());
    }
    stats.broadphase_pairs = static_cast<std::uint32_t>(pair_count);
    stats.contacts = static_cast<std::uint32_t>(contacts.size());
    contact_pairs = pairs;
    contact_poses.resize(proxies.size());
    for (std::size_t i = 0; i < proxies.size(); ++i)
        contact_poses[i] = proxies[i].pose;
}
} // namespace fxp_physics
