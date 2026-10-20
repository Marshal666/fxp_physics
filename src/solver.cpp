#include "detail.hpp"
#include <atomic>
#include <condition_variable>
#include <exception>
#include <mutex>
#include <thread>

namespace fxp_physics {
namespace detail {
namespace {
// Jobs write disjoint output: contact generation owns per-job buffers; an
// island solver owns its dynamic bodies and accumulated impulses. Shared
// static and kinematic bodies are read-only. The atomic counter chooses where
// a job runs, never the order of constraints or the final merged contacts.
class JobExecutor {
  public:
    explicit JobExecutor(unsigned workers) {
        try {
            for (unsigned i = 1; i < workers; ++i)
                threads_.emplace_back([this] { worker(); });
        } catch (...) {
            stop();
            throw;
        }
    }
    ~JobExecutor() { stop(); }
    unsigned workers() const { return static_cast<unsigned>(threads_.size()) + 1; }
    void run(std::size_t count, std::function<void(std::size_t)> job) {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            job_ = std::move(job);
            count_ = count;
            next_.store(0, std::memory_order_relaxed);
            remaining_ = threads_.size();
            exception_ = {};
            ++generation_;
        }
        ready_.notify_all();
        consume();
        {
            std::unique_lock<std::mutex> lock(mutex_);
            done_.wait(lock, [this] { return remaining_ == 0; });
            job_ = {};
        }
        if (exception_)
            std::rethrow_exception(exception_);
    }

  private:
    std::vector<std::thread> threads_;
    std::mutex mutex_;
    std::condition_variable ready_, done_;
    std::function<void(std::size_t)> job_;
    std::atomic<std::size_t> next_{0};
    std::size_t count_{}, remaining_{};
    std::uint64_t generation_{};
    std::exception_ptr exception_;
    bool stopping_ = false;

    void stop() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stopping_ = true;
        }
        ready_.notify_all();
        for (auto &thread : threads_)
            thread.join();
    }
    void consume() {
        for (;;) {
            const auto index = next_.fetch_add(1, std::memory_order_relaxed);
            if (index >= count_)
                return;
            try {
                job_(index);
            } catch (...) {
                std::lock_guard<std::mutex> lock(mutex_);
                if (!exception_)
                    exception_ = std::current_exception();
                next_.store(count_, std::memory_order_relaxed);
                return;
            }
        }
    }
    void worker() {
        std::uint64_t seen = 0;
        for (;;) {
            std::unique_lock<std::mutex> lock(mutex_);
            ready_.wait(lock, [&] { return stopping_ || generation_ != seen; });
            if (stopping_)
                return;
            seen = generation_;
            lock.unlock();
            consume();
            lock.lock();
            if (--remaining_ == 0)
                done_.notify_one();
        }
    }
};

struct EffectiveMass {
    Real denominator{}, reciprocal{};
    explicit EffectiveMass(Real value = {}) : denominator((std::max)(value, 0_r)) {
        if (denominator > 0) {
            reciprocal = 1_r / denominator;
            // Huge or tiny inertias can give an unrepresentable reciprocal,
            // or too few significant bits. Retain division in those rare rows.
            if (reciprocal == Real::max() || reciprocal.raw_value() < (INT64_C(1) << 20))
                reciprocal = 0;
        }
    }
    Real apply(Real numerator) const {
        return reciprocal > 0 ? numerator * reciprocal
                              : denominator > 0 ? numerator / denominator : 0_r;
    }
};

struct SolverBody {
    Vec3 center{}, inertia_x{}, inertia_y{}, inertia_z{};
    Real inverse_mass{};
    bool active = false, prepared = false;
    Vec3 angular_response(Vec3 value) const {
        return inertia_x * value.x + inertia_y * value.y + inertia_z * value.z;
    }
};
struct SolverAxis {
    Vec3 direction{}, jacobian_a{}, jacobian_b{}, linear_a{}, linear_b{}, angular_a{}, angular_b{};
    EffectiveMass mass;
};
struct ContactRow {
    std::size_t a{}, b{}, contact{}, island{};
    SolverAxis normal, tangent1, tangent2;
    Real bias{}, mu_static{}, mu_dynamic{}, normal_impulse{}, tangent_impulse1{}, tangent_impulse2{};
};
struct SpringRow {
    std::size_t a{}, b{}, spring{}, island{};
    SolverAxis axis;
    EffectiveMass soft_mass;
    Real gamma{}, bias{}, max{}, impulse{};
};
struct SolverIsland {
    std::size_t contact_begin{}, contact_count{}, spring_begin{}, spring_count{};
};
unsigned execution_workers(std::uint32_t requested) {
    return requested ? requested : (std::max)(1u, (std::min)(8u, std::thread::hardware_concurrency()));
}
} // namespace

// This cache is execution state only. Neither its allocation history nor the
// number of workers is part of a save, checksum, or simulation decision.
struct SolverData {
    std::vector<SolverBody> bodies;
    std::vector<std::size_t> needed_bodies;
    std::vector<ContactRow> contacts;
    std::vector<SpringRow> springs;
    std::vector<SolverIsland> islands;
    std::vector<std::size_t> island_for_root, contact_order, spring_order, cursor;
    std::unique_ptr<JobExecutor> executor;
};
} // namespace detail

using namespace detail;

void World::set_worker_count(std::uint32_t count) {
    impl_->mutable_access();
    require(count <= 64, "worker count must be 0 (automatic) or 1..64");
    if (impl_->requested_worker_count != count) {
        impl_->requested_worker_count = count;
        if (impl_->solver_data)
            impl_->solver_data->executor.reset();
    }
}
std::uint32_t World::worker_count() const noexcept {
    return impl_->requested_worker_count;
}

void World::Impl::parallel_for(std::size_t count, std::size_t work,
                              std::function<void(std::size_t)> job) {
    const auto workers = execution_workers(requested_worker_count);
    if (workers > 1 && count > 1 && work >= 128) {
        if (!solver_data)
            solver_data = std::make_shared<SolverData>();
        auto &executor = solver_data->executor;
        if (!executor || executor->workers() != workers)
            executor = std::make_unique<JobExecutor>(workers);
        executor->run(count, std::move(job));
    } else {
        for (std::size_t i = 0; i < count; ++i)
            job(i);
    }
}

void World::Impl::solve(Real h) {
    if (contacts.empty() && springs.empty())
        return;
    if (std::none_of(bodies.begin(), bodies.end(), [](const Body &body) {
            return body.desc.type == BodyType::dynamic && !body.sleeping;
        }))
        return;
    if (!solver_data)
        solver_data = std::make_shared<SolverData>();
    auto &data = *solver_data;
    data.bodies.resize(bodies.size());
    for (auto &body : data.bodies)
        body.prepared = false;
    data.needed_bodies.clear();
    data.contacts.clear();
    data.springs.clear();
    data.contacts.reserve(contacts.size());
    data.springs.reserve(springs.size());
    data.islands.clear();
    data.island_for_root.assign(bodies.size(), bodies.size());
    const auto roots = islands();
    const auto active = [&](std::size_t index) {
        return index != bodies.size() && bodies[index].desc.type == BodyType::dynamic &&
               !bodies[index].sleeping;
    };
    const auto require_body = [&](std::size_t index) {
        auto &prepared = data.bodies[index];
        if (!prepared.prepared) {
            prepared.prepared = true;
            data.needed_bodies.push_back(index);
        }
    };
    const auto island_for = [&](std::size_t a, std::size_t b) {
        const auto root = roots[active(a) ? a : b];
        auto &index = data.island_for_root[root];
        if (index == bodies.size()) {
            index = data.islands.size();
            data.islands.emplace_back();
        }
        return index;
    };
    const auto axis_for = [&](std::size_t a, std::size_t b, Vec3 pa, Vec3 pb, Vec3 direction) {
        const auto &ba = data.bodies[a];
        SolverAxis axis;
        axis.direction = direction;
        axis.jacobian_a = cross(pa - ba.center, direction);
        axis.linear_a = direction * ba.inverse_mass;
        axis.angular_a = ba.angular_response(axis.jacobian_a);
        Real mass = ba.inverse_mass + dot(axis.jacobian_a, axis.angular_a);
        if (b != bodies.size()) {
            const auto &bb = data.bodies[b];
            axis.jacobian_b = cross(pb - bb.center, direction);
            axis.linear_b = direction * bb.inverse_mass;
            axis.angular_b = bb.angular_response(axis.jacobian_b);
            mass += bb.inverse_mass + dot(axis.jacobian_b, axis.angular_b);
        }
        axis.mass = EffectiveMass(mass);
        return axis;
    };
    for (std::size_t i = 0; i < contacts.size(); ++i) {
        const auto &contact = contacts[i];
        const auto &ca = colliders[collider_index(contact.a)];
        const auto &cb = colliders[collider_index(contact.b)];
        const auto a = body_index(ca.desc.body), b = body_index(cb.desc.body);
        if (!active(a) && !active(b))
            continue;
        ContactRow row;
        row.a = a;
        row.b = b;
        row.contact = i;
        row.island = island_for(a, b);
        require_body(a);
        require_body(b);
        data.contacts.push_back(row);
        ++data.islands[row.island].contact_count;
    }
    for (std::size_t i = 0; i < springs.size(); ++i) {
        auto &spring = springs[i];
        const auto a = body_index(spring.desc.body_a),
                   b = spring.desc.body_b ? body_index(spring.desc.body_b) : bodies.size();
        if (!active(a) && !active(b))
            continue;
        SpringRow row;
        row.a = a;
        row.b = b;
        row.spring = i;
        row.island = island_for(a, b);
        require_body(a);
        if (b != bodies.size())
            require_body(b);
        data.springs.push_back(row);
        ++data.islands[row.island].spring_count;
    }
    if (data.islands.empty())
        return;

    // Mark participating bodies serially, then give each body to exactly one
    // job. Its cached quaternion coefficients are mutable, so multiple contact
    // jobs must never initialize that body's inertia concurrently.
    parallel_for(data.needed_bodies.size(), data.needed_bodies.size(), [&](std::size_t i) {
        const auto index = data.needed_bodies[i];
        const auto &body = bodies[index];
        auto &prepared = data.bodies[index];
        prepared = {};
        prepared.prepared = true;
        prepared.center = com(body);
        prepared.active = active(index);
        if (prepared.active) {
            prepared.inverse_mass = body.inverse_mass;
            prepared.inertia_x = angular_mass(body, {1, 0, 0});
            prepared.inertia_y = angular_mass(body, {0, 1, 0});
            prepared.inertia_z = angular_mass(body, {0, 0, 1});
        }
    });

    // Counting sort preserves contact and spring order within every island,
    // including when removals leave gaps in handles and row/contact indices.
    std::size_t contact_offset = 0, spring_offset = 0;
    for (auto &island : data.islands) {
        island.contact_begin = contact_offset;
        island.spring_begin = spring_offset;
        contact_offset += island.contact_count;
        spring_offset += island.spring_count;
    }
    data.contact_order.resize(data.contacts.size());
    data.spring_order.resize(data.springs.size());
    data.cursor.resize(data.islands.size());
    for (std::size_t i = 0; i < data.islands.size(); ++i)
        data.cursor[i] = data.islands[i].contact_begin;
    for (std::size_t i = 0; i < data.contacts.size(); ++i)
        data.contact_order[data.cursor[data.contacts[i].island]++] = i;
    for (std::size_t i = 0; i < data.islands.size(); ++i)
        data.cursor[i] = data.islands[i].spring_begin;
    for (std::size_t i = 0; i < data.springs.size(); ++i)
        data.spring_order[data.cursor[data.springs[i].island]++] = i;

    const auto prepare_contact = [&](ContactRow &row) {
        const auto &contact = contacts[row.contact];
        const auto &ca = colliders[collider_index(contact.a)];
        const auto &cb = colliders[collider_index(contact.b)];
        const Vec3 point = (transform(bodies[row.a].desc.pose, contact.local_a) +
                            transform(bodies[row.b].desc.pose, contact.local_b)) /
                           2_r;
        const Vec3 normal = contact.normal, an = abs(normal);
        Vec3 basis{1, 0, 0};
        if (an.y < an.x && an.y <= an.z)
            basis = {0, 1, 0};
        else if (an.z < an.x && an.z < an.y)
            basis = {0, 0, 1};
        const Vec3 tangent1 = unit(cross(normal, basis)), tangent2 = cross(normal, tangent1);
        row.normal = axis_for(row.a, row.b, point, point, normal);
        row.tangent1 = axis_for(row.a, row.b, point, point, tangent1);
        row.tangent2 = axis_for(row.a, row.b, point, point, tangent2);
        row.bias = contact.separation > 0
                       ? -contact.separation / h
                       : (std::min)(config.max_correction_speed,
                                    config.stabilization *
                                        (std::max)(-contact.separation - config.contact_slop, 0_r) / h);
        if (contact.bounce > 0)
            row.bias = (std::max)(row.bias, contact.bounce);
        row.mu_static = fxp::sqrt(ca.desc.material.static_friction * cb.desc.material.static_friction);
        row.mu_dynamic = fxp::sqrt(ca.desc.material.dynamic_friction * cb.desc.material.dynamic_friction);
        row.normal_impulse = contact.normal_impulse;
        row.tangent_impulse1 = dot(contact.tangent_impulse, tangent1);
        row.tangent_impulse2 = dot(contact.tangent_impulse, tangent2);
    };
    const auto prepare_spring = [&](SpringRow &row) {
        const auto &spring = springs[row.spring];
        const Vec3 pa = transform(bodies[row.a].desc.pose, spring.desc.anchor_a);
        const Vec3 pb = row.b != bodies.size() ? transform(bodies[row.b].desc.pose, spring.desc.anchor_b)
                                               : spring.desc.anchor_b;
        const Real omega = 2_r * fxp::pi_v<Real> * spring.desc.frequency;
        const Real term = 2_r * spring.desc.damping_ratio * omega + h * omega * omega;
        // A positive frequency can round to zero stiffness in this timestep.
        // Such a spring has no representable response, rather than a division
        // by zero when computing its softness.
        if (term <= 0 || h * term <= 0)
            return;
        row.axis = axis_for(row.a, row.b, pa, pb, unit(pb - pa));
        row.gamma = row.axis.mass.denominator / (h * term);
        row.soft_mass = EffectiveMass(row.axis.mass.denominator + row.gamma);
        row.bias = (length(pb - pa) - spring.desc.rest_length) * omega * omega / term;
        row.max = spring.desc.max_force * h;
        row.impulse = clamp(spring.impulse, -row.max, row.max);
    };
    const auto velocity = [&](std::size_t a, std::size_t b, const SolverAxis &axis) {
        const auto &va = bodies[a].desc;
        if (b == bodies.size())
            return -dot(va.linear_velocity, axis.direction) - dot(va.angular_velocity, axis.jacobian_a);
        const auto &vb = bodies[b].desc;
        return dot(vb.linear_velocity - va.linear_velocity, axis.direction) +
               dot(vb.angular_velocity, axis.jacobian_b) - dot(va.angular_velocity, axis.jacobian_a);
    };
    const auto apply = [&](std::size_t a, std::size_t b, const SolverAxis &axis, Real change) {
        if (change == 0)
            return;
        if (data.bodies[a].active) {
            auto &body = bodies[a].desc;
            body.linear_velocity -= axis.linear_a * change;
            body.angular_velocity -= axis.angular_a * change;
        }
        if (b != bodies.size() && data.bodies[b].active) {
            auto &body = bodies[b].desc;
            body.linear_velocity += axis.linear_b * change;
            body.angular_velocity += axis.angular_b * change;
        }
    };
    const auto solve_island = [&](std::size_t index) {
        const auto &island = data.islands[index];
        const auto contact_end = island.contact_begin + island.contact_count;
        const auto spring_end = island.spring_begin + island.spring_count;
        for (auto i = island.contact_begin; i < contact_end; ++i) {
            auto &row = data.contacts[data.contact_order[i]];
            prepare_contact(row);
            apply(row.a, row.b, row.normal, row.normal_impulse);
            apply(row.a, row.b, row.tangent1, row.tangent_impulse1);
            apply(row.a, row.b, row.tangent2, row.tangent_impulse2);
        }
        for (auto i = island.spring_begin; i < spring_end; ++i) {
            auto &row = data.springs[data.spring_order[i]];
            prepare_spring(row);
            apply(row.a, row.b, row.axis, row.impulse);
        }
        for (std::uint32_t iteration = 0; iteration < config.velocity_iterations; ++iteration) {
            for (auto i = island.contact_begin; i < contact_end; ++i) {
                auto &row = data.contacts[data.contact_order[i]];
                if (row.normal.mass.denominator <= 0)
                    continue;
                const Real next = (std::max)(0_r, row.normal_impulse +
                                                     row.normal.mass.apply(row.bias - velocity(row.a, row.b, row.normal)));
                apply(row.a, row.b, row.normal, next - row.normal_impulse);
                row.normal_impulse = next;
                Real t1 = row.tangent_impulse1 -
                          row.tangent1.mass.apply(velocity(row.a, row.b, row.tangent1));
                Real t2 = row.tangent_impulse2 -
                          row.tangent2.mass.apply(velocity(row.a, row.b, row.tangent2));
                const Vec3 tangent{t1, t2, 0};
                if (length_greater_than(tangent, row.mu_static * next)) {
                    const Real magnitude = length(tangent);
                    if (magnitude > 0) {
                        const Real scale = (std::min)(1_r, row.mu_dynamic * next / magnitude);
                        t1 *= scale;
                        t2 *= scale;
                    }
                }
                apply(row.a, row.b, row.tangent1, t1 - row.tangent_impulse1);
                apply(row.a, row.b, row.tangent2, t2 - row.tangent_impulse2);
                row.tangent_impulse1 = t1;
                row.tangent_impulse2 = t2;
            }
            for (auto i = island.spring_begin; i < spring_end; ++i) {
                auto &row = data.springs[data.spring_order[i]];
                if (row.axis.mass.denominator <= 0)
                    continue;
                const Real next =
                    clamp(row.impulse - row.soft_mass.apply(velocity(row.a, row.b, row.axis) +
                                                           row.bias + row.gamma * row.impulse),
                          -row.max, row.max);
                apply(row.a, row.b, row.axis, next - row.impulse);
                row.impulse = next;
            }
        }
        for (auto i = island.contact_begin; i < contact_end; ++i) {
            const auto &row = data.contacts[data.contact_order[i]];
            auto &contact = contacts[row.contact];
            contact.normal_impulse = row.normal_impulse;
            contact.tangent_impulse = row.tangent1.direction * row.tangent_impulse1 +
                                      row.tangent2.direction * row.tangent_impulse2;
        }
        for (auto i = island.spring_begin; i < spring_end; ++i) {
            const auto &row = data.springs[data.spring_order[i]];
            springs[row.spring].impulse = row.impulse;
        }
    };
    parallel_for(data.islands.size(), data.contacts.size() + data.springs.size(), solve_island);
}
} // namespace fxp_physics
