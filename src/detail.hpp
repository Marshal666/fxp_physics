#pragma once
#include <fxp_physics/physics.hpp>
#include <algorithm>
#include <limits>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace fxp_physics {
namespace detail {
struct SolverData;
inline Real clamp(Real x, Real lo, Real hi) {
    return (std::max)(lo, (std::min)(x, hi));
}
inline Real dot(Vec3 a, Vec3 b) {
    return fxp::dot(a, b);
}
inline Vec3 cross(Vec3 a, Vec3 b) {
    return fxp::cross(a, b);
}
inline Real length(Vec3 a) {
    return fxp::length(a);
}
// Compare the rounded raw length without a square root. For an integer raw
// threshold t, round(sqrt(sum)) > t iff sum >= t*t + t + 1. Unlike dot(v,v),
// raw squares do not round each component before the comparison.
inline bool length_greater_than(Vec3 v, Real threshold) {
    if (threshold < 0)
        return true;
    if (threshold == Real::max())
        return false; // length saturates at Real::max().
    const auto t = static_cast<std::uint64_t>(threshold.raw_value());
    const auto boundary = fxp::detail::add(fxp::detail::mul(t, t), {0, t + 1});
    fxp::detail::u128 sum{};
    for (unsigned i = 0; i < 3; ++i) {
        const auto magnitude = fxp::detail::magnitude(v[i].raw_value());
        sum = fxp::detail::add(sum, fxp::detail::mul(magnitude, magnitude));
        if (sum >= boundary)
            return true;
    }
    return false;
}
inline Vec3 unit(Vec3 a, Vec3 fallback = {1, 0, 0}) {
    return fxp::normalize(&a) ? a : fallback;
}
inline Vec3 abs(Vec3 a) {
    return {fxp::abs(a.x), fxp::abs(a.y), fxp::abs(a.z)};
}
inline Vec3 minimum(Vec3 a, Vec3 b) {
    return {(std::min)(a.x, b.x), (std::min)(a.y, b.y), (std::min)(a.z, b.z)};
}
inline Vec3 maximum(Vec3 a, Vec3 b) {
    return {(std::max)(a.x, b.x), (std::max)(a.y, b.y), (std::max)(a.z, b.z)};
}
inline Quat inverse(Quat q) {
    q.unit_inverse();
    return q;
}
// Cache only quaternion coefficients, preserving the rounding/operation order
// of quat::rotate. A precomputed rotation matrix would give different bits.
struct Rotation {
    fxp::vec4<Real> quaternion{0, 0, 0, 1};
    Vec3 axis{};
    Real scale = 1_r, twice_w = 2_r;
    void set(Quat q) {
        const auto v = q.get_internal_vector();
        if (quaternion == v)
            return;
        quaternion = v;
        axis = {v.x, v.y, v.z};
        scale = v.w * v.w - dot(axis, axis);
        twice_w = 2_r * v.w;
    }
    Vec3 rotate(Vec3 v, bool inverted = false) const {
        if (zero_axis())
            return v * scale;
        const Vec3 a = inverted ? -axis : axis;
        return v * scale + twice_w * cross(a, v) + (2_r * dot(a, v)) * a;
    }
  private:
    bool zero_axis() const { return axis == Vec3{}; }
};
inline Vec3 transform(Pose p, Vec3 v) {
    return p.position + p.rotation.rotate(v);
}
inline Vec3 untransform(Pose p, Vec3 v) {
    return inverse(p.rotation).rotate(v - p.position);
}
inline Pose relative(Pose p, Pose w) {
    return {untransform(p, w.position), inverse(p.rotation) * w.rotation};
}
inline bool zero(Vec3 v) {
    return v == Vec3{};
}
inline Vec3 limit(Vec3 v, Real maximum_length) {
    if (dot(v, v) <= maximum_length * maximum_length)
        return v;
    const Real len = length(v);
    return len > maximum_length ? v * (maximum_length / len) : v;
}
void require(bool condition, const char *message);
void validate_vec(Vec3 value, Real limit, const char *message);
Pose validated_pose(Pose pose, Real position_limit = 8192_r);
void validate_material(const Material &material);
void validate_config(const Config &config);
Inertia invert(Inertia tensor);
Vec3 multiply(Inertia tensor, Vec3 value);

struct Bounds {
    Vec3 min{Real::max(), Real::max(), Real::max()};
    Vec3 max{Real::min(), Real::min(), Real::min()};
    void add(Vec3 p) {
        min = minimum(min, p);
        max = maximum(max, p);
    }
    void add(Bounds b) {
        add(b.min);
        add(b.max);
    }
    Bounds expanded(Real margin) const {
        const Vec3 v{margin, margin, margin};
        return {min - v, max + v};
    }
    bool intersects(Bounds b) const {
        return min.x <= b.max.x && max.x >= b.min.x && min.y <= b.max.y && max.y >= b.min.y &&
               min.z <= b.max.z && max.z >= b.min.z;
    }
};
Bounds transformed(Bounds b, Pose p);
struct Tree {
    struct Node {
        Bounds bounds;
        std::uint32_t left{}, right{}, item = 0xffffffffu;
    };
    std::vector<Node> nodes;
    void build(const std::vector<Bounds> &bounds);
    void refit(const std::vector<Bounds> &bounds);
    std::vector<std::uint32_t> query(Bounds bounds) const;
    void query(Bounds bounds, std::vector<std::uint32_t> &out) const;

  private:
    std::uint32_t build_node(std::vector<std::uint32_t> &ids, std::size_t begin, std::size_t end,
                             const std::vector<Bounds> &bounds);
};
struct Geometry {
    bool rounded = false;
    bool box = false;
    std::array<Vec3, 3> box_axes{}, box_half_axes{};
    Vec3 a{}, b{};
    Real radius{};
    std::vector<Vec3> vertices;
    std::vector<std::vector<std::uint32_t>> faces;
    std::vector<Vec3> normals;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> edges;
    Vec3 center{};
    std::vector<Vec3> centered_vertices;
    mutable std::vector<Vec3> edge_directions;
    Vec3 support(Vec3 direction) const;
    const std::vector<Vec3> &directions() const;
};
Geometry geometry(const Shape::Data &shape, Pose pose, std::uint32_t triangle = 0xffffffffu);
Bounds bounds(const Shape::Data &shape, Pose pose);
std::vector<std::uint32_t> triangles(const Shape::Data &shape, Pose pose, Bounds area);
struct ContactPoint {
    Vec3 point{};
    Real separation{};
};
struct Separation {
    Real distance = Real::min();
    Vec3 normal{1, 0, 0};
    Vec3 point_a{}, point_b{};
    std::int32_t box_axis = -1;
};
Separation separate(const Geometry &a, const Geometry &b, Real rejection_threshold = Real::max());
std::vector<ContactPoint> manifold(const Geometry &a, const Geometry &b, Separation separation, Real margin);
Geometry translated(Geometry value, Vec3 delta);
} // namespace detail

struct Shape::Data {
    ShapeType type{};
    Vec3 half_extents{};
    Real radius{}, half_segment{}, spacing{};
    std::uint32_t width{}, depth{};
    std::vector<Vec3> vertices;
    std::vector<std::vector<std::uint32_t>> faces;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> edges;
    std::vector<Triangle> mesh_triangles;
    std::vector<Real> heights;
    detail::Bounds bounds;
    detail::Tree tree;
};

struct World::Impl {
    struct Body {
        BodyId id{};
        BodyDesc desc;
        Pose previous{}, target{};
        bool has_target = false, sleeping = false;
        std::uint32_t quiet_ticks = 0;
        Vec3 force{}, torque{};
        Real inverse_mass{};
        Inertia inverse_inertia{};
        mutable detail::Rotation rotation;
    };
    struct Collider {
        ColliderId id{};
        ColliderDesc desc;
    };
    struct Spring {
        SpringId id{};
        SpringDesc desc;
        Real impulse{};
    };
    struct Contact {
        ColliderId a{}, b{};
        std::uint32_t triangle_a = 0xffffffffu, triangle_b = 0xffffffffu;
        Vec3 local_a{}, local_b{}, normal{}, tangent_impulse{};
        Real normal_impulse{}, separation{}, bounce{};
    };
    struct Pair {
        ColliderId a{}, b{};
        BodyId body_a{}, body_b{};
        bool sensor{};
        Vec3 point{}, normal{};
    };
    struct Proxy {
        ColliderId collider{};
        Pose pose{};
        detail::Bounds bounds{};
        std::size_t body{};
        std::optional<detail::Geometry> geometry;
    };
    Config config;
    Random rng;
    std::uint64_t tick = 0;
    BodyId next_body = 1;
    ColliderId next_collider = 1;
    SpringId next_spring = 1;
    std::array<std::uint32_t, 32> layers{};
    std::vector<Body> bodies;
    std::vector<Collider> colliders;
    std::vector<Spring> springs;
    std::vector<Contact> contacts;
    std::vector<Pair> pairs;
    std::vector<CollisionEvent> events;
    // Transient exact-pose contact cache. Queries never update this history.
    std::vector<Pose> contact_poses;
    std::vector<Pair> contact_pairs;
    std::shared_ptr<detail::SolverData> solver_data;
    struct ContactData;
    std::shared_ptr<ContactData> contact_data;
    std::uint32_t requested_worker_count = 0;
    mutable bool dirty = true;
    mutable detail::Tree tree;
    mutable std::vector<Proxy> proxies;
    mutable std::vector<detail::Bounds> proxy_bounds;
    mutable unsigned tree_refits = 0;
    mutable unsigned observer_depth = 0;
    bool stepping = false, fixed_update = false;
    Stats stats{};
    explicit Impl(Config cfg, std::uint64_t seed) : config(cfg), rng(seed) { layers.fill(0xffffffffu); }
    void mutable_access() const;
    std::size_t body_index(BodyId id) const;
    std::size_t collider_index(ColliderId id) const;
    Body &get(BodyId id) { return bodies[body_index(id)]; }
    const Body &get(BodyId id) const { return bodies[body_index(id)]; }
    Pose local_pose(const Collider &collider) const;
    void rebuild() const;
    void recalculate_inertia(BodyId id);
    void wake_all();
    void wake_connected(BodyId id);
    void invalidate(BodyId id);
    Vec3 com(const Body &b) const {
        // The usual origin-centred body needs no quaternion rotation here.
        return detail::zero(b.desc.center_of_mass) ? b.desc.pose.position
                                                  : detail::transform(b.desc.pose, b.desc.center_of_mass);
    }
    Vec3 angular_mass(const Body &b, Vec3 v) const;
    Vec3 point_velocity(const Body &b, Vec3 point) const;
    void find_contacts();
    void parallel_for(std::size_t count, std::size_t work,
                      std::function<void(std::size_t)> job);
    void solve(Real h);
    std::vector<std::size_t> islands() const;
    void propagate_waking();
    void sleep_islands();
    void finish_events(std::vector<Pair> old_pairs);
};
} // namespace fxp_physics
