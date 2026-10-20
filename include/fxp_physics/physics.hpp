#pragma once
#include <fxp/fxp.hpp>
#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <vector>

namespace fxp_physics {
using Real = fxp::real;
using Vec3 = fxp::vec3<Real>;
using Quat = fxp::quat<Real>;
using BodyId = std::uint64_t;
using ColliderId = std::uint64_t;
using SpringId = std::uint64_t;
using Bytes = std::vector<std::uint8_t>;
using namespace fxp::literals;

struct Pose {
    Vec3 position{};
    Quat rotation{};
};
Pose compose(const Pose &parent, const Pose &local);
Pose interpolate(const Pose &previous, const Pose &current, Real alpha);

// Symmetric tensor in body-local axes. An all-zero tensor requests automatic
// inertia from the compound's bounds. Explicit tensors must be positive definite.
struct Inertia {
    Real xx{}, yy{}, zz{}, xy{}, xz{}, yz{};
};
struct Material {
    Real static_friction = 0.6_r;
    Real dynamic_friction = 0.4_r;
    Real restitution = 0_r;
};

enum class ShapeType : std::uint8_t { box, sphere, capsule, convex, mesh, heightfield };
struct Triangle {
    std::uint32_t a{}, b{}, c{};
};

// Immutable cooked collision assets. Convex faces are outward-wound, planar,
// convex polygons. Mesh triangles are two-sided surfaces, not filled volumes.
class Shape {
  public:
    struct Data;
    static std::shared_ptr<const Shape> box(Vec3 half_extents);
    static std::shared_ptr<const Shape> sphere(Real radius);
    // Capsule axis is local Y; half_segment excludes the hemispherical ends.
    static std::shared_ptr<const Shape> capsule(Real radius, Real half_segment);
    static std::shared_ptr<const Shape> convex(std::vector<Vec3> vertices,
                                               std::vector<std::vector<std::uint32_t>> faces);
    static std::shared_ptr<const Shape> mesh(std::vector<Vec3> vertices, std::vector<Triangle> triangles);
    // Row-major samples: heights[z * width + x]. Origin is sample (0,0).
    static std::shared_ptr<const Shape> heightfield(std::uint32_t width, std::uint32_t depth, Real spacing,
                                                    std::vector<Real> heights);
    ShapeType type() const noexcept;
    const Data &data() const noexcept;

  private:
    friend class World;
    explicit Shape(std::shared_ptr<const Data> data);
    std::shared_ptr<const Data> data_;
};
using ShapePtr = std::shared_ptr<const Shape>;

enum class BodyType : std::uint8_t { static_body, kinematic, dynamic };
struct BodyDesc {
    BodyType type = BodyType::dynamic;
    Pose pose{};
    Vec3 linear_velocity{}, angular_velocity{};
    Real mass = 1_r;
    Vec3 center_of_mass{};
    Inertia inertia{};
    Real linear_damping = 0.05_r, angular_damping = 0.05_r;
    Real gravity_scale = 1_r;
    // Bit 0/1/2 locks local angular X/Y/Z. Linear motion remains unconstrained.
    std::uint8_t rotation_locks = 0;
    bool allow_sleep = true;
    std::uint64_t user_data = 0;
};
struct BodyState {
    BodyId id{};
    BodyDesc properties{};
    bool sleeping = false;
    std::uint32_t quiet_ticks = 0;
};
struct ColliderDesc {
    BodyId body{};
    ShapePtr shape;
    Pose local_pose{};
    // Optional collider parent on the same body. Local pose is relative to it.
    ColliderId parent{};
    Material material{};
    std::uint8_t layer = 0;
    bool sensor = false;
    std::uint64_t user_data = 0;
};
struct ColliderState {
    ColliderId id{};
    ColliderDesc properties{};
};
struct SpringDesc {
    BodyId body_a{}, body_b{}; // body_b == 0 attaches anchor_b to the world
    Vec3 anchor_a{}, anchor_b{};
    Real rest_length = 1_r;
    Real frequency = 2_r;
    Real damping_ratio = 0.7_r;
    Real max_force = 100000_r;
    bool collide_connected = false;
};

struct Config {
    Real timestep = 1_r / 30_r;
    std::uint32_t substeps = 2, velocity_iterations = 8;
    Vec3 gravity{0, -9.81_r, 0};
    Real contact_slop = 0.005_r;
    Real contact_margin = 0.01_r;
    Real stabilization = 0.2_r;
    Real max_correction_speed = 2_r;
    Real restitution_threshold = 1_r;
    Real sleep_linear_speed = 0.08_r, sleep_angular_speed = 0.08_r;
    std::uint32_t sleep_ticks = 30;
    Real max_linear_speed = 128_r, max_angular_speed = 32_r;
};

enum class EventType : std::uint8_t { begin, stay, end };
struct CollisionEvent {
    std::uint64_t tick{};
    ColliderId collider_a{}, collider_b{};
    BodyId body_a{}, body_b{};
    EventType type{};
    bool sensor{};
    Vec3 point{}, normal{}; // normal points from A toward B
};
struct QueryFilter {
    std::uint32_t layers = 0xffffffffu;
    BodyId ignore_body{};
    bool include_sensors = false;
};
struct QueryHit {
    BodyId body{};
    ColliderId collider{};
    std::uint32_t triangle = 0xffffffffu;
    Real fraction{}, distance{};
    Vec3 point{}, normal{}; // outward normal of target
    bool started_overlapping = false;
};
struct CastResult {
    std::vector<QueryHit> hits; // sorted by fraction, collider, triangle
    // False means at least one candidate exhausted the conservative advancement
    // budget. Never interpret an incomplete result as a guaranteed clear path.
    bool complete = true;
};
struct Stats {
    std::uint32_t bodies{}, awake_dynamic_bodies{}, colliders{}, broadphase_pairs{}, contacts{};
};

// SplitMix64, including fixed bounded-sampling rules; unsigned arithmetic only.
class Random {
  public:
    explicit Random(std::uint64_t seed = 0) : state_(seed) {}
    std::uint64_t next();
    std::uint32_t bounded(std::uint32_t exclusive_max);
    Real unit(); // [0,1), exactly 32 random fractional bits
    std::uint64_t state() const noexcept { return state_; }
    void restore(std::uint64_t value) noexcept { state_ = value; }

  private:
    std::uint64_t state_;
};

class World {
  public:
    explicit World(Config config = {}, std::uint64_t seed = 0);
    ~World();
    World(World &&) noexcept;
    World &operator=(World &&) noexcept;
    World(const World &) = delete;
    World &operator=(const World &) = delete;

    // Mutations are allowed between steps and in fixed_update. All callers must
    // submit the same sequence. Event/render callbacks may only inspect the world.
    BodyId create_body(const BodyDesc &desc = {});
    ColliderId create_collider(const ColliderDesc &desc);
    SpringId create_spring(const SpringDesc &desc);
    void remove_body(BodyId id);         // removes its colliders and connected springs
    void remove_collider(ColliderId id); // includes descendants
    void remove_spring(SpringId id);
    BodyState body(BodyId id) const;
    ColliderState collider(ColliderId id) const;
    std::vector<BodyId> body_ids() const;
    void set_velocity(BodyId id, Vec3 linear, Vec3 angular = {});
    void add_force(BodyId id, Vec3 force, Vec3 torque = {});
    void add_force_at(BodyId id, Vec3 force, Vec3 world_point);
    void add_impulse(BodyId id, Vec3 impulse, Vec3 angular_impulse = {});
    void add_impulse_at(BodyId id, Vec3 impulse, Vec3 world_point);
    void teleport(BodyId id, Pose pose, bool clear_velocity = false);
    void set_kinematic_target(BodyId id, Pose target);
    void set_collider_pose(ColliderId id, Pose local_pose);
    void set_collider_parent(ColliderId id, ColliderId parent, bool preserve_world_pose = true);
    void set_collider_filter(ColliderId id, std::uint8_t layer, bool sensor);
    void set_material(ColliderId id, Material material);
    void set_layer_collision(std::uint8_t a, std::uint8_t b, bool enabled);
    bool layer_collision(std::uint8_t a, std::uint8_t b) const;
    void wake(BodyId id);

    // Execution setting only; excluded from snapshots and state checksums.
    // Zero selects automatic parallelism (at most eight workers), one is serial.
    // The configured count is preserved by load; results are independent of it.
    void set_worker_count(std::uint32_t count);
    std::uint32_t worker_count() const noexcept;

    void step(const std::function<void(World &)> &fixed_update = {});
    const std::vector<CollisionEvent> &events() const noexcept;
    void dispatch_events(const std::function<void(const CollisionEvent &)> &callback) const;
    void render(Real alpha, const std::function<void(BodyId, const Pose &)> &callback) const;
    Pose interpolated_pose(BodyId id, Real alpha) const;

    std::vector<QueryHit> overlap(const Shape &shape, Pose pose, QueryFilter filter = {}) const;
    CastResult cast(const Shape &shape, Pose pose, Vec3 displacement, QueryFilter filter = {}) const;
    CastResult raycast(Vec3 origin, Vec3 direction, Real distance, QueryFilter filter = {}) const;
    CastResult spherecast(Real radius, Pose pose, Vec3 displacement, QueryFilter filter = {}) const;
    CastResult boxcast(Vec3 half_extents, Pose pose, Vec3 displacement, QueryFilter filter = {}) const;

    // Canonical little-endian, versioned, self-contained bytes; no pointers or
    // struct padding. Load is transactional. Callbacks belong to the host.
    Bytes save() const;
    void load(const Bytes &bytes);
    std::uint64_t checksum() const;
    std::uint64_t tick() const noexcept;
    const Config &config() const noexcept;
    Stats stats() const;
    std::uint64_t random_u64();
    std::uint32_t random_bounded(std::uint32_t exclusive_max);
    Real random_unit();

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
} // namespace fxp_physics
