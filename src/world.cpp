#include "detail.hpp"
#include <numeric>
#include <optional>

namespace fxp_physics {
using namespace detail;
namespace {
bool immovable(BodyType type) {
    return type != BodyType::dynamic;
}
bool concave(ShapeType type) {
    return type == ShapeType::mesh || type == ShapeType::heightfield;
}
template <class Container> std::size_t index_of(const Container &values, std::uint64_t id) {
    // IDs are dense until removal. Keep the common lookup constant time, while
    // retaining ordered lookup for sparse IDs and loaded worlds.
    if (id && id <= values.size() && values[static_cast<std::size_t>(id - 1)].id == id)
        return static_cast<std::size_t>(id - 1);
    const auto it = std::lower_bound(values.begin(), values.end(), id,
                                     [](const auto &item, auto key) { return item.id < key; });
    if (it == values.end() || it->id != id)
        throw std::out_of_range("unknown physics handle");
    return static_cast<std::size_t>(it - values.begin());
}
template <class Container> std::uint64_t next_id(std::uint64_t &next, const Container &values) {
    require(values.size() < 1000000, "object capacity exceeded");
    if (next == UINT64_MAX)
        throw std::overflow_error("physics ID space exhausted");
    return next++;
}
struct Observer {
    unsigned &depth;
    explicit Observer(unsigned &d) : depth(d) { ++depth; }
    ~Observer() { --depth; }
};
Vec3 angular_velocity(Quat start, Quat end, Real dt) {
    Quat delta = end * inverse(start);
    delta.normalize();
    delta.minimize_rotation_angle();
    const auto v = delta.get_internal_vector();
    const Vec3 xyz{v.x, v.y, v.z};
    const Real len = length(xyz);
    return len > 0.0000001_r ? xyz * (2_r * fxp::atan2(len, v.w) / dt / len) : Vec3{};
}
Quat advance_rotation(Quat q, Vec3 velocity, Real dt) {
    if (zero(velocity))
        return q;
    const auto v = q.get_internal_vector();
    const Quat derivative = Quat(fxp::vec4<Real>{velocity.x, velocity.y, velocity.z, 0}) * q;
    Quat next(v + derivative.get_internal_vector() * (dt / 2_r));
    next.normalize();
    next.minimize_rotation_angle();
    return next;
}
} // namespace

World::World(Config config, std::uint64_t seed) {
    validate_config(config);
    impl_ = std::make_unique<Impl>(config, seed);
}
World::~World() = default;
World::World(World &&) noexcept = default;
World &World::operator=(World &&) noexcept = default;
void World::Impl::mutable_access() const {
    if (observer_depth || (stepping && !fixed_update))
        throw std::logic_error("world mutation outside fixed-update boundary");
}
std::size_t World::Impl::body_index(BodyId id) const {
    return index_of(bodies, id);
}
std::size_t World::Impl::collider_index(ColliderId id) const {
    return index_of(colliders, id);
}
Pose World::Impl::local_pose(const Collider &collider) const {
    Pose pose = collider.desc.local_pose;
    auto parent = collider.desc.parent;
    for (std::size_t depth = 0; parent; ++depth) {
        if (depth >= 64 || depth >= colliders.size())
            throw std::logic_error("collider hierarchy cycle or depth exceeds 64");
        const auto &p = colliders[collider_index(parent)];
        pose = compose(p.desc.local_pose, pose);
        parent = p.desc.parent;
    }
    return pose;
}
void World::Impl::rebuild() const {
    if (!dirty)
        return;
    bool topology_changed = proxies.size() != colliders.size();
    bool bounds_changed = topology_changed;
    proxies.resize(colliders.size());
    proxy_bounds.resize(colliders.size());
    for (std::size_t i = 0; i < colliders.size(); ++i) {
        const auto &c = colliders[i];
        auto &proxy = proxies[i];
        const auto b = body_index(c.desc.body);
        const auto pose = compose(bodies[b].desc.pose, local_pose(c));
        const bool new_collider = proxy.collider != c.id;
        topology_changed |= new_collider;
        if (new_collider || proxy.pose.position != pose.position ||
            proxy.pose.rotation.get_internal_vector() != pose.rotation.get_internal_vector()) {
            proxy.collider = c.id;
            proxy.pose = pose;
            proxy.bounds = bounds(c.desc.shape->data(), pose);
            proxy.geometry.reset();
            bounds_changed = true;
        }
        // Body indices may change even when this collider's pose does not.
        proxy.body = b;
        proxy_bounds[i] = proxy.bounds;
    }
    if (topology_changed || (bounds_changed && tree_refits >= 32)) {
        tree.build(proxy_bounds);
        tree_refits = 0;
    } else if (bounds_changed) {
        // Refit preserves all leaf bounds; query results are always ID-sorted,
        // so topology and query frequency cannot affect contact order/results.
        tree.refit(proxy_bounds);
        ++tree_refits;
    }
    dirty = false;
}
void World::Impl::recalculate_inertia(BodyId id) {
    auto &b = get(id);
    b.inverse_mass = immovable(b.desc.type) ? 0_r : 1_r / b.desc.mass;
    if (immovable(b.desc.type)) {
        b.inverse_inertia = {};
        return;
    }
    const auto m = b.desc.inertia;
    if (m.xx != 0 || m.yy != 0 || m.zz != 0 || m.xy != 0 || m.xz != 0 || m.yz != 0) {
        b.inverse_inertia = invert(m);
        return;
    }
    Bounds box;
    bool found = false;
    for (const auto &c : colliders)
        if (c.desc.body == id && !c.desc.sensor) {
            box.add(bounds(c.desc.shape->data(), local_pose(c)));
            found = true;
        }
    if (!found)
        box = {Vec3{-0.5_r, -0.5_r, -0.5_r}, Vec3{0.5_r, 0.5_r, 0.5_r}};
    // Bounding-box mass approximation, shifted to the requested center of mass.
    const Vec3 size = maximum(box.max - box.min, Vec3{0.01_r, 0.01_r, 0.01_r});
    const Vec3 offset = (box.max + box.min) / 2_r - b.desc.center_of_mass;
    const Real mass = b.desc.mass;
    Inertia tensor{
        mass * ((size.y * size.y + size.z * size.z) / 12_r + offset.y * offset.y + offset.z * offset.z),
        mass * ((size.x * size.x + size.z * size.z) / 12_r + offset.x * offset.x + offset.z * offset.z),
        mass * ((size.x * size.x + size.y * size.y) / 12_r + offset.x * offset.x + offset.y * offset.y),
        -mass * offset.x * offset.y,
        -mass * offset.x * offset.z,
        -mass * offset.y * offset.z};
    b.inverse_inertia = invert(tensor);
}
void World::Impl::wake_all() {
    for (auto &b : bodies)
        if (b.desc.type == BodyType::dynamic) {
            b.sleeping = false;
            b.quiet_ticks = 0;
        }
}
void World::Impl::wake_connected(BodyId id) {
    auto &changed = get(id);
    if (contacts.empty() && springs.empty() && changed.desc.type != BodyType::static_body) {
        changed.sleeping = false;
        changed.quiet_ticks = 0;
        return;
    }
    const auto roots = islands();
    const auto root = roots[body_index(id)];
    if (changed.desc.type == BodyType::static_body) {
        wake_all();
        return;
    }
    if (changed.desc.type == BodyType::kinematic) {
        std::vector<bool> affected(bodies.size(), false);
        const auto affect = [&](BodyId other) {
            if (other)
                affected[roots[body_index(other)]] = true;
        };
        for (const auto &contact : contacts) {
            const auto a = colliders[collider_index(contact.a)].desc.body;
            const auto b = colliders[collider_index(contact.b)].desc.body;
            if (a == id)
                affect(b);
            if (b == id)
                affect(a);
        }
        for (const auto &spring : springs) {
            if (spring.desc.body_a == id)
                affect(spring.desc.body_b);
            if (spring.desc.body_b == id)
                affect(spring.desc.body_a);
        }
        for (std::size_t i = 0; i < bodies.size(); ++i)
            if (affected[roots[i]]) {
                bodies[i].sleeping = false;
                bodies[i].quiet_ticks = 0;
            }
        return;
    }
    for (std::size_t i = 0; i < bodies.size(); ++i)
        if (roots[i] == root) {
            bodies[i].sleeping = false;
            bodies[i].quiet_ticks = 0;
        }
}
void World::Impl::invalidate(BodyId id) {
    contact_poses.clear();
    contact_pairs.clear();
    contacts.erase(std::remove_if(contacts.begin(), contacts.end(),
                                  [&](const auto &c) {
                                      return colliders[collider_index(c.a)].desc.body == id ||
                                             colliders[collider_index(c.b)].desc.body == id;
                                  }),
                   contacts.end());
    dirty = true;
}
BodyId World::create_body(const BodyDesc &input) {
    auto &w = *impl_;
    w.mutable_access();
    BodyDesc d = input;
    d.pose = validated_pose(d.pose);
    require(static_cast<unsigned>(d.type) <= static_cast<unsigned>(BodyType::dynamic), "invalid body type");
    require(d.mass >= 0.01_r && d.mass <= 10000, "mass must be 0.01..10000 kg");
    validate_vec(d.center_of_mass, 128, "center of mass too large");
    validate_vec(d.linear_velocity, w.config.max_linear_speed, "velocity too large");
    validate_vec(d.angular_velocity, w.config.max_angular_speed, "angular velocity too large");
    require(d.linear_damping >= 0 && d.linear_damping <= 100 && d.angular_damping >= 0 &&
                d.angular_damping <= 100 && fxp::abs(d.gravity_scale) <= 10 && d.rotation_locks <= 7,
            "invalid body settings");
    if (d.type == BodyType::static_body)
        require(zero(d.linear_velocity) && zero(d.angular_velocity), "static body cannot have velocity");
    const bool automatic = d.inertia.xx == 0 && d.inertia.yy == 0 && d.inertia.zz == 0 && d.inertia.xy == 0 &&
                           d.inertia.xz == 0 && d.inertia.yz == 0;
    if (!automatic)
        (void)invert(d.inertia);
    Impl::Body b;
    b.id = next_id(w.next_body, w.bodies);
    b.desc = d;
    b.previous = d.pose;
    b.target = d.pose;
    w.bodies.push_back(b);
    w.recalculate_inertia(b.id);
    w.dirty = true;
    return b.id;
}
ColliderId World::create_collider(const ColliderDesc &input) {
    auto &w = *impl_;
    w.mutable_access();
    ColliderDesc d = input;
    const auto &body = w.get(d.body);
    require(bool(d.shape), "collider needs a shape");
    d.local_pose = validated_pose(d.local_pose, 128_r);
    require(d.layer < 32, "layer must be in [0,31]");
    validate_material(d.material);
    require(!concave(d.shape->type()) || body.desc.type == BodyType::static_body,
            "concave/terrain colliders require static bodies");
    if (d.parent)
        require(w.colliders[w.collider_index(d.parent)].desc.body == d.body,
                "parent must belong to same body");
    const auto id = next_id(w.next_collider, w.colliders);
    w.colliders.push_back({id, d});
    try {
        validate_vec(w.local_pose(w.colliders.back()).position, 128, "compound offset too large");
        w.recalculate_inertia(d.body);
    } catch (...) {
        w.colliders.pop_back();
        --w.next_collider;
        throw;
    }
    w.wake_connected(d.body);
    w.invalidate(d.body);
    return id;
}
SpringId World::create_spring(const SpringDesc &d) {
    auto &w = *impl_;
    w.mutable_access();
    w.get(d.body_a);
    if (d.body_b)
        w.get(d.body_b);
    require(d.body_a != d.body_b, "spring bodies must differ");
    validate_vec(d.anchor_a, 128, "spring anchor too large");
    validate_vec(d.anchor_b, d.body_b ? 128_r : 8192_r, "spring anchor too large");
    require(d.rest_length >= 0 && d.rest_length <= 256 && d.frequency > 0 && d.frequency <= 30 &&
                d.damping_ratio >= 0 && d.damping_ratio <= 4 && d.max_force > 0 && d.max_force <= 1000000,
            "invalid spring");
    const auto id = next_id(w.next_spring, w.springs);
    w.springs.push_back({id, d, 0});
    w.contact_poses.clear();
    w.contact_pairs.clear();
    w.wake_connected(d.body_a);
    if (d.body_b)
        w.wake_connected(d.body_b);
    return id;
}
void World::remove_body(BodyId id) {
    auto &w = *impl_;
    w.mutable_access();
    const auto index = w.body_index(id);
    w.wake_connected(id);
    w.invalidate(id);
    w.colliders.erase(
        std::remove_if(w.colliders.begin(), w.colliders.end(), [&](auto &c) { return c.desc.body == id; }),
        w.colliders.end());
    w.springs.erase(std::remove_if(w.springs.begin(), w.springs.end(),
                                   [&](auto &s) { return s.desc.body_a == id || s.desc.body_b == id; }),
                    w.springs.end());
    w.bodies.erase(w.bodies.begin() + static_cast<std::ptrdiff_t>(index));
    w.dirty = true;
}
void World::remove_collider(ColliderId id) {
    auto &w = *impl_;
    w.mutable_access();
    const auto body = w.colliders[w.collider_index(id)].desc.body;
    w.wake_connected(body);
    w.invalidate(body);
    std::vector<ColliderId> removed{id};
    bool changed = true;
    while (changed) {
        changed = false;
        for (const auto &c : w.colliders)
            if (std::find(removed.begin(), removed.end(), c.desc.parent) != removed.end() &&
                std::find(removed.begin(), removed.end(), c.id) == removed.end()) {
                removed.push_back(c.id);
                changed = true;
            }
    }
    w.colliders.erase(std::remove_if(w.colliders.begin(), w.colliders.end(),
                                     [&](const auto &c) {
                                         return std::find(removed.begin(), removed.end(), c.id) !=
                                                removed.end();
                                     }),
                      w.colliders.end());
    w.recalculate_inertia(body);
    w.dirty = true;
}
void World::remove_spring(SpringId id) {
    auto &w = *impl_;
    w.mutable_access();
    const auto i = index_of(w.springs, id);
    w.wake_connected(w.springs[i].desc.body_a);
    if (w.springs[i].desc.body_b)
        w.wake_connected(w.springs[i].desc.body_b);
    w.springs.erase(w.springs.begin() + static_cast<std::ptrdiff_t>(i));
    w.contact_poses.clear();
    w.contact_pairs.clear();
}
BodyState World::body(BodyId id) const {
    const auto &b = impl_->get(id);
    return {id, b.desc, b.sleeping, b.quiet_ticks};
}
ColliderState World::collider(ColliderId id) const {
    const auto &c = impl_->colliders[impl_->collider_index(id)];
    return {id, c.desc};
}
std::vector<BodyId> World::body_ids() const {
    std::vector<BodyId> ids;
    for (const auto &b : impl_->bodies)
        ids.push_back(b.id);
    return ids;
}
void World::set_velocity(BodyId id, Vec3 linear, Vec3 angular) {
    auto &w = *impl_;
    w.mutable_access();
    auto &b = w.get(id);
    require(b.desc.type != BodyType::static_body, "static body velocity");
    validate_vec(linear, w.config.max_linear_speed, "velocity too large");
    validate_vec(angular, w.config.max_angular_speed, "angular velocity too large");
    if (linear != b.desc.linear_velocity || angular != b.desc.angular_velocity)
        w.wake_connected(id);
    b.desc.linear_velocity = limit(linear, w.config.max_linear_speed);
    b.desc.angular_velocity = limit(angular, w.config.max_angular_speed);
    b.has_target = false;
}
void World::add_force(BodyId id, Vec3 force, Vec3 torque) {
    auto &w = *impl_;
    w.mutable_access();
    auto &b = w.get(id);
    require(b.desc.type == BodyType::dynamic, "force requires dynamic body");
    validate_vec(b.force + force, 1000000, "force too large");
    validate_vec(b.torque + torque, 1000000, "torque too large");
    if (!zero(force) || !zero(torque))
        w.wake_connected(id);
    b.force += force;
    b.torque += torque;
}
void World::add_force_at(BodyId id, Vec3 force, Vec3 point) {
    validate_vec(point, 8192, "point too large");
    validate_vec(force, 1000000, "force too large");
    const Vec3 arm = point - impl_->com(impl_->get(id));
    validate_vec(arm, 128, "force arm too large");
    add_force(id, force, cross(arm, force));
}
void World::add_impulse(BodyId id, Vec3 impulse, Vec3 angular) {
    auto &w = *impl_;
    w.mutable_access();
    auto &b = w.get(id);
    require(b.desc.type == BodyType::dynamic, "impulse requires dynamic body");
    validate_vec(impulse, 1000000, "impulse too large");
    validate_vec(angular, 1000000, "angular impulse too large");
    if (!zero(impulse) || !zero(angular))
        w.wake_connected(id);
    b.desc.linear_velocity =
        limit(b.desc.linear_velocity + impulse * b.inverse_mass, w.config.max_linear_speed);
    b.desc.angular_velocity =
        limit(b.desc.angular_velocity + w.angular_mass(b, angular), w.config.max_angular_speed);
}
void World::add_impulse_at(BodyId id, Vec3 impulse, Vec3 point) {
    validate_vec(point, 8192, "point too large");
    validate_vec(impulse, 1000000, "impulse too large");
    const Vec3 arm = point - impl_->com(impl_->get(id));
    validate_vec(arm, 128, "impulse arm too large");
    add_impulse(id, impulse, cross(arm, impulse));
}
void World::teleport(BodyId id, Pose pose, bool clear_velocity) {
    auto &w = *impl_;
    w.mutable_access();
    pose = validated_pose(pose);
    auto &b = w.get(id);
    w.wake_connected(id);
    w.invalidate(id);
    b.desc.pose = pose;
    b.previous = pose;
    b.target = pose;
    b.has_target = false;
    if (clear_velocity) {
        b.desc.linear_velocity = {};
        b.desc.angular_velocity = {};
    }
}
void World::set_kinematic_target(BodyId id, Pose target) {
    auto &w = *impl_;
    w.mutable_access();
    target = validated_pose(target);
    auto &b = w.get(id);
    require(b.desc.type == BodyType::kinematic, "target requires kinematic body");
    require(length(target.position - b.desc.pose.position) / w.config.timestep <= w.config.max_linear_speed,
            "target exceeds linear speed limit; use teleport");
    require(length(angular_velocity(b.desc.pose.rotation, target.rotation, w.config.timestep)) <=
                w.config.max_angular_speed,
            "target exceeds angular speed limit");
    b.target = target;
    b.has_target = true;
}
void World::set_collider_pose(ColliderId id, Pose pose) {
    auto &w = *impl_;
    w.mutable_access();
    pose = validated_pose(pose, 128);
    auto &c = w.colliders[w.collider_index(id)];
    const auto old = c.desc.local_pose;
    c.desc.local_pose = pose;
    try {
        for (const auto &child : w.colliders)
            if (child.desc.body == c.desc.body)
                validate_vec(w.local_pose(child).position, 128, "compound offset too large");
        w.recalculate_inertia(c.desc.body);
    } catch (...) {
        c.desc.local_pose = old;
        throw;
    }
    w.wake_connected(c.desc.body);
    w.invalidate(c.desc.body);
}
void World::set_collider_parent(ColliderId id, ColliderId parent, bool preserve) {
    auto &w = *impl_;
    w.mutable_access();
    auto &c = w.colliders[w.collider_index(id)];
    if (parent)
        require(w.colliders[w.collider_index(parent)].desc.body == c.desc.body, "parent on different body");
    for (auto p = parent; p; p = w.colliders[w.collider_index(p)].desc.parent)
        require(p != id, "hierarchy cycle");
    const auto old = c.desc;
    const Pose absolute = w.local_pose(c);
    c.desc.parent = parent;
    if (preserve)
        c.desc.local_pose =
            parent ? relative(w.local_pose(w.colliders[w.collider_index(parent)]), absolute) : absolute;
    try {
        set_collider_pose(id, c.desc.local_pose);
    } catch (...) {
        c.desc = old;
        throw;
    }
}
void World::set_collider_filter(ColliderId id, std::uint8_t layer, bool sensor) {
    auto &w = *impl_;
    w.mutable_access();
    require(layer < 32, "layer out of range");
    auto &c = w.colliders[w.collider_index(id)];
    const auto old = c.desc;
    c.desc.layer = layer;
    c.desc.sensor = sensor;
    try {
        w.recalculate_inertia(c.desc.body);
    } catch (...) {
        c.desc = old;
        throw;
    }
    w.wake_connected(c.desc.body);
    w.invalidate(c.desc.body);
}
void World::set_material(ColliderId id, Material material) {
    auto &w = *impl_;
    w.mutable_access();
    validate_material(material);
    auto &c = w.colliders[w.collider_index(id)];
    c.desc.material = material;
    w.wake_connected(c.desc.body);
    w.invalidate(c.desc.body);
}
void World::set_layer_collision(std::uint8_t a, std::uint8_t b, bool enabled) {
    auto &w = *impl_;
    w.mutable_access();
    require(a < 32 && b < 32, "layer out of range");
    if (enabled) {
        w.layers[a] |= std::uint32_t(1) << b;
        w.layers[b] |= std::uint32_t(1) << a;
    } else {
        w.layers[a] &= ~(std::uint32_t(1) << b);
        w.layers[b] &= ~(std::uint32_t(1) << a);
    }
    w.contacts.clear();
    w.contact_poses.clear();
    w.contact_pairs.clear();
    w.wake_all();
}
bool World::layer_collision(std::uint8_t a, std::uint8_t b) const {
    require(a < 32 && b < 32, "layer out of range");
    return (impl_->layers[a] & (std::uint32_t(1) << b)) != 0;
}
void World::wake(BodyId id) {
    impl_->mutable_access();
    impl_->wake_connected(id);
}

Vec3 World::Impl::angular_mass(const Body &b, Vec3 v) const {
    if (immovable(b.desc.type) || zero(v) || b.desc.rotation_locks == 7)
        return {};
    b.rotation.set(b.desc.pose.rotation);
    Vec3 local = b.rotation.rotate(v, true);
    for (unsigned i = 0; i < 3; ++i)
        if (b.desc.rotation_locks & (1u << i))
            local[i] = 0;
    local = multiply(b.inverse_inertia, local);
    for (unsigned i = 0; i < 3; ++i)
        if (b.desc.rotation_locks & (1u << i))
            local[i] = 0;
    return b.rotation.rotate(local);
}
Vec3 World::Impl::point_velocity(const Body &b, Vec3 point) const {
    if (zero(b.desc.angular_velocity))
        return b.desc.linear_velocity;
    return b.desc.linear_velocity + cross(b.desc.angular_velocity, point - com(b));
}
std::vector<std::size_t> World::Impl::islands() const {
    std::vector<std::size_t> parent(bodies.size());
    std::iota(parent.begin(), parent.end(), std::size_t(0));
    const auto root = [&](std::size_t i) {
        while (parent[i] != i) {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    };
    const auto join = [&](BodyId aid, BodyId bid) {
        if (!aid || !bid)
            return;
        const auto a = body_index(aid), b = body_index(bid);
        if (immovable(bodies[a].desc.type) || immovable(bodies[b].desc.type))
            return;
        const auto ra = root(a), rb = root(b);
        parent[(std::max)(ra, rb)] = (std::min)(ra, rb);
    };
    for (const auto &c : contacts)
        join(colliders[collider_index(c.a)].desc.body, colliders[collider_index(c.b)].desc.body);
    for (const auto &s : springs)
        join(s.desc.body_a, s.desc.body_b);
    for (std::size_t i = 0; i < parent.size(); ++i)
        parent[i] = root(i);
    return parent;
}
void World::Impl::propagate_waking() {
    if (std::none_of(bodies.begin(), bodies.end(), [](const Body &b) {
            return (b.desc.type == BodyType::dynamic && !b.sleeping) ||
                   (b.desc.type == BodyType::kinematic &&
                    (!zero(b.desc.linear_velocity) || !zero(b.desc.angular_velocity)));
        }))
        return;
    const auto roots = islands();
    std::vector<bool> active(bodies.size(), false);
    for (std::size_t i = 0; i < bodies.size(); ++i)
        if (bodies[i].desc.type == BodyType::dynamic && !bodies[i].sleeping)
            active[roots[i]] = true;
    const auto moving = [](const Body &b) {
        return b.desc.type == BodyType::kinematic &&
               (!zero(b.desc.linear_velocity) || !zero(b.desc.angular_velocity));
    };
    for (const auto &c : contacts) {
        const auto a = body_index(colliders[collider_index(c.a)].desc.body),
                   b = body_index(colliders[collider_index(c.b)].desc.body);
        if (moving(bodies[a]))
            active[roots[b]] = true;
        if (moving(bodies[b]))
            active[roots[a]] = true;
    }
    for (const auto &s : springs)
        if (s.desc.body_b) {
            const auto a = body_index(s.desc.body_a), b = body_index(s.desc.body_b);
            if (moving(bodies[a]))
                active[roots[b]] = true;
            if (moving(bodies[b]))
                active[roots[a]] = true;
        }
    for (std::size_t i = 0; i < bodies.size(); ++i)
        if (active[roots[i]] && bodies[i].sleeping) {
            bodies[i].sleeping = false;
            bodies[i].quiet_ticks = 0;
        }
}
void World::Impl::sleep_islands() {
    if (std::none_of(bodies.begin(), bodies.end(), [](const Body &b) {
            return b.desc.type == BodyType::dynamic && !b.sleeping;
        }))
        return;
    const auto roots = islands();
    std::vector<bool> quiet(bodies.size(), true);
    std::vector<std::uint32_t> counters(bodies.size(), config.sleep_ticks);
    const auto moving = [](const Body &body) {
        return body.desc.type == BodyType::kinematic &&
               (!zero(body.desc.linear_velocity) || !zero(body.desc.angular_velocity));
    };
    for (const auto &contact : contacts) {
        const auto a = body_index(colliders[collider_index(contact.a)].desc.body);
        const auto b = body_index(colliders[collider_index(contact.b)].desc.body);
        if (moving(bodies[a]))
            quiet[roots[b]] = false;
        if (moving(bodies[b]))
            quiet[roots[a]] = false;
    }
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        const auto &b = bodies[i];
        if (immovable(b.desc.type))
            continue;
        if (!b.desc.allow_sleep || length_greater_than(b.desc.linear_velocity, config.sleep_linear_speed) ||
            length_greater_than(b.desc.angular_velocity, config.sleep_angular_speed) ||
            !zero(b.force) || !zero(b.torque))
            quiet[roots[i]] = false;
        counters[roots[i]] = (std::min)(counters[roots[i]], b.quiet_ticks);
    }
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        auto &b = bodies[i];
        if (immovable(b.desc.type))
            continue;
        b.quiet_ticks = quiet[roots[i]] ? (std::min)(config.sleep_ticks, counters[roots[i]] + 1) : 0;
        b.sleeping = b.quiet_ticks >= config.sleep_ticks;
        if (b.sleeping) {
            b.desc.linear_velocity = {};
            b.desc.angular_velocity = {};
        }
    }
}
void World::Impl::finish_events(std::vector<Pair> old) {
    events.clear();
    std::size_t a = 0, b = 0;
    const auto emit = [&](const Pair &p, EventType type) {
        events.push_back({tick, p.a, p.b, p.body_a, p.body_b, type, p.sensor, p.point, p.normal});
    };
    while (a < old.size() || b < pairs.size()) {
        if (b == pairs.size() ||
            (a < old.size() && std::make_pair(old[a].a, old[a].b) < std::make_pair(pairs[b].a, pairs[b].b)))
            emit(old[a++], EventType::end);
        else if (a == old.size() ||
                 std::make_pair(pairs[b].a, pairs[b].b) < std::make_pair(old[a].a, old[a].b))
            emit(pairs[b++], EventType::begin);
        else {
            if (old[a].sensor != pairs[b].sensor) {
                emit(old[a], EventType::end);
                emit(pairs[b], EventType::begin);
            } else
                emit(pairs[b], EventType::stay);
            ++a;
            ++b;
        }
    }
}
void World::step(const std::function<void(World &)> &fixed_update) {
    auto &w = *impl_;
    w.mutable_access();
    if (w.stepping)
        throw std::logic_error("recursive step");
    if (w.tick == UINT64_MAX)
        throw std::overflow_error("tick space exhausted");
    struct Guard {
        Impl &w;
        ~Guard() {
            w.stepping = false;
            w.fixed_update = false;
        }
    } guard{w};
    w.stepping = true;
    w.fixed_update = true;
    if (fixed_update)
        fixed_update(*this);
    w.fixed_update = false;
    const auto old_pairs = w.pairs;
    for (auto &b : w.bodies) {
        b.previous = b.desc.pose;
        if (b.desc.type == BodyType::kinematic && b.has_target) {
            b.desc.linear_velocity =
                (transform(b.target, b.desc.center_of_mass) - w.com(b)) / w.config.timestep;
            b.desc.angular_velocity =
                angular_velocity(b.desc.pose.rotation, b.target.rotation, w.config.timestep);
        }
    }
    const Real h = w.config.timestep / Real(w.config.substeps);
    std::vector<Impl::Pair> seen_pairs;
    for (std::uint32_t substep = 0; substep < w.config.substeps; ++substep) {
        for (auto &b : w.bodies) {
            if (b.desc.type == BodyType::static_body || b.sleeping)
                continue;
            if (b.desc.type == BodyType::dynamic) {
                b.desc.linear_velocity +=
                    (w.config.gravity * b.desc.gravity_scale + b.force * b.inverse_mass) * h;
                b.desc.angular_velocity += w.angular_mass(b, b.torque) * h;
                b.desc.linear_velocity = b.desc.linear_velocity / (1_r + b.desc.linear_damping * h);
                b.desc.angular_velocity = b.desc.angular_velocity / (1_r + b.desc.angular_damping * h);
                if (!zero(b.desc.angular_velocity)) {
                    Vec3 local = inverse(b.desc.pose.rotation).rotate(b.desc.angular_velocity);
                    for (unsigned i = 0; i < 3; ++i)
                        if (b.desc.rotation_locks & (1u << i))
                            local[i] = 0;
                    b.desc.angular_velocity = b.desc.pose.rotation.rotate(local);
                }
            }
            b.desc.linear_velocity = limit(b.desc.linear_velocity, w.config.max_linear_speed);
            b.desc.angular_velocity = limit(b.desc.angular_velocity, w.config.max_angular_speed);
            if (b.desc.type == BodyType::kinematic && b.has_target) {
                b.desc.pose = interpolate(b.previous, b.target, Real(substep + 1) / Real(w.config.substeps));
                w.dirty = true;
            } else if (b.desc.type == BodyType::kinematic &&
                       (!zero(b.desc.linear_velocity) || !zero(b.desc.angular_velocity))) {
                const Vec3 center = w.com(b) + b.desc.linear_velocity * h;
                b.desc.pose.rotation = advance_rotation(b.desc.pose.rotation, b.desc.angular_velocity, h);
                b.desc.pose.position = zero(b.desc.center_of_mass)
                                           ? center : center - b.desc.pose.rotation.rotate(b.desc.center_of_mass);
                w.dirty = true;
            }
            validate_vec(b.desc.pose.position, 8192, "body left supported world range");
        }
        w.find_contacts();
        w.propagate_waking();
        w.solve(h);
        // Observe contacts from every substep, including brief trigger crossings.
        std::vector<Impl::Pair> merged;
        std::size_t old_index = 0, new_index = 0;
        while (old_index < seen_pairs.size() || new_index < w.pairs.size()) {
            if (new_index == w.pairs.size() ||
                (old_index < seen_pairs.size() &&
                 std::make_pair(seen_pairs[old_index].a, seen_pairs[old_index].b) <
                     std::make_pair(w.pairs[new_index].a, w.pairs[new_index].b)))
                merged.push_back(seen_pairs[old_index++]);
            else {
                if (old_index < seen_pairs.size() && seen_pairs[old_index].a == w.pairs[new_index].a &&
                    seen_pairs[old_index].b == w.pairs[new_index].b)
                    ++old_index;
                merged.push_back(w.pairs[new_index++]);
            }
        }
        seen_pairs = std::move(merged);
        bool dynamic_advanced = false;
        for (auto &b : w.bodies)
            if (b.desc.type == BodyType::dynamic && !b.sleeping) {
                dynamic_advanced = true;
                b.desc.linear_velocity = limit(b.desc.linear_velocity, w.config.max_linear_speed);
                b.desc.angular_velocity = limit(b.desc.angular_velocity, w.config.max_angular_speed);
                const Vec3 center = w.com(b) + b.desc.linear_velocity * h;
                b.desc.pose.rotation = advance_rotation(b.desc.pose.rotation, b.desc.angular_velocity, h);
                b.desc.pose.position = zero(b.desc.center_of_mass)
                                           ? center : center - b.desc.pose.rotation.rotate(b.desc.center_of_mass);
                validate_vec(b.desc.pose.position, 8192, "body left supported world range");
            }
        if (dynamic_advanced)
            w.dirty = true;
    }
    w.sleep_islands();
    for (auto &b : w.bodies) {
        b.force = {};
        b.torque = {};
        b.desc.linear_velocity = limit(b.desc.linear_velocity, w.config.max_linear_speed);
        b.desc.angular_velocity = limit(b.desc.angular_velocity, w.config.max_angular_speed);
        if (b.has_target) {
            b.has_target = false;
            b.desc.linear_velocity = {};
            b.desc.angular_velocity = {};
        }
    }
    w.pairs = std::move(seen_pairs);
    ++w.tick;
    w.finish_events(old_pairs);
}
const std::vector<CollisionEvent> &World::events() const noexcept {
    return impl_->events;
}
void World::dispatch_events(const std::function<void(const CollisionEvent &)> &callback) const {
    if (!callback)
        return;
    Observer observer(impl_->observer_depth);
    for (const auto &event : impl_->events)
        callback(event);
}
Pose World::interpolated_pose(BodyId id, Real alpha) const {
    const auto &b = impl_->get(id);
    return interpolate(b.previous, b.desc.pose, alpha);
}
void World::render(Real alpha, const std::function<void(BodyId, const Pose &)> &callback) const {
    require(alpha >= 0 && alpha <= 1, "interpolation alpha out of range");
    if (!callback)
        return;
    Observer observer(impl_->observer_depth);
    for (const auto &b : impl_->bodies)
        callback(b.id, interpolate(b.previous, b.desc.pose, alpha));
}

std::vector<QueryHit> World::overlap(const Shape &shape, Pose pose, QueryFilter filter) const {
    require(!concave(shape.type()), "query shape must be convex");
    pose = validated_pose(pose);
    const auto &w = *impl_;
    w.rebuild();
    const auto area = bounds(shape.data(), pose);
    const auto source = geometry(shape.data(), pose);
    std::vector<QueryHit> hits;
    for (auto i : w.tree.query(area)) {
        const auto &c = w.colliders[i];
        if (c.desc.body == filter.ignore_body || !(filter.layers & (std::uint32_t(1) << c.desc.layer)) ||
            (c.desc.sensor && !filter.include_sensors))
            continue;
        const auto &p = w.proxies[i];
        for (auto triangle : triangles(c.desc.shape->data(), p.pose, area)) {
            const auto target = geometry(c.desc.shape->data(), p.pose, triangle);
            const auto s = separate(source, target, 0_r);
            if (s.distance <= 0) {
                const auto points = manifold(source, target, s, 0);
                hits.push_back({c.desc.body, c.id, triangle, 0, 0, points.front().point, -s.normal, true});
                break;
            }
        }
    }
    return hits;
}
CastResult World::cast(const Shape &shape, Pose pose, Vec3 delta, QueryFilter filter) const {
    require(!concave(shape.type()), "cast shape must be convex");
    pose = validated_pose(pose);
    validate_vec(delta, 16384, "cast displacement too large");
    const auto &w = *impl_;
    w.rebuild();
    const Real travel = length(delta);
    const auto source = geometry(shape.data(), pose);
    Bounds area = bounds(shape.data(), pose);
    Pose end = pose;
    end.position += delta;
    area.add(bounds(shape.data(), end));
    area = area.expanded(0.00001_r);
    CastResult result;
    for (auto i : w.tree.query(area)) {
        const auto &c = w.colliders[i];
        if (c.desc.body == filter.ignore_body || !(filter.layers & (std::uint32_t(1) << c.desc.layer)) ||
            (c.desc.sensor && !filter.include_sensors))
            continue;
        const auto &p = w.proxies[i];
        for (auto triangle : triangles(c.desc.shape->data(), p.pose, area)) {
            const auto target = geometry(c.desc.shape->data(), p.pose, triangle);
            Real fraction = 0;
            bool finished = false;
            for (unsigned iteration = 0; iteration < 64; ++iteration) {
                const auto current = translated(source, delta * fraction);
                const auto s = separate(current, target);
                if (s.distance <= 0.00001_r) {
                    const auto points = manifold(current, target, s, 0.00001_r);
                    result.hits.push_back({c.desc.body, c.id, triangle, fraction, travel * fraction,
                                           points.front().point, -s.normal,
                                           fraction == 0 && s.distance <= 0});
                    finished = true;
                    break;
                }
                const Real closing = dot(delta, s.normal);
                if (closing <= 0) {
                    finished = true;
                    break;
                }
                const Real advance = s.distance / closing;
                if (advance <= 0)
                    break;
                const Real next = fraction + advance;
                if (next > 1) {
                    finished = true;
                    break;
                }
                if (next <= fraction)
                    break;
                fraction = next;
            }
            if (!finished)
                result.complete = false;
        }
    }
    std::sort(result.hits.begin(), result.hits.end(), [](const auto &a, const auto &b) {
        return std::tie(a.fraction, a.collider, a.triangle) < std::tie(b.fraction, b.collider, b.triangle);
    });
    return result;
}
CastResult World::raycast(Vec3 origin, Vec3 direction, Real distance, QueryFilter filter) const {
    validate_vec(direction, 8192, "ray direction too large");
    require(!zero(direction) && distance >= 0 && distance <= 16384, "invalid ray");
    // A zero-radius point uses exactly the same convex sweep machinery.
    auto data = std::make_shared<Shape::Data>();
    data->type = ShapeType::sphere;
    data->bounds = {Vec3{}, Vec3{}};
    const Shape point(data);
    return cast(point, {origin, Quat{}}, unit(direction) * distance, filter);
}
CastResult World::spherecast(Real radius, Pose pose, Vec3 delta, QueryFilter filter) const {
    return cast(*Shape::sphere(radius), pose, delta, filter);
}
CastResult World::boxcast(Vec3 half, Pose pose, Vec3 delta, QueryFilter filter) const {
    return cast(*Shape::box(half), pose, delta, filter);
}
std::uint64_t World::tick() const noexcept {
    return impl_->tick;
}
const Config &World::config() const noexcept {
    return impl_->config;
}
Stats World::stats() const {
    auto result = impl_->stats;
    result.bodies = static_cast<std::uint32_t>(impl_->bodies.size());
    result.colliders = static_cast<std::uint32_t>(impl_->colliders.size());
    result.awake_dynamic_bodies = 0;
    for (const auto &b : impl_->bodies)
        if (b.desc.type == BodyType::dynamic && !b.sleeping)
            ++result.awake_dynamic_bodies;
    return result;
}
std::uint64_t World::random_u64() {
    impl_->mutable_access();
    return impl_->rng.next();
}
std::uint32_t World::random_bounded(std::uint32_t bound) {
    impl_->mutable_access();
    return impl_->rng.bounded(bound);
}
Real World::random_unit() {
    impl_->mutable_access();
    return impl_->rng.unit();
}
} // namespace fxp_physics
