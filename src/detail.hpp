#pragma once
#include <fxp_physics/physics.hpp>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace fxp_physics {
namespace detail {
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
    std::vector<std::uint32_t> query(Bounds bounds) const;

  private:
    std::uint32_t build_node(std::vector<std::uint32_t> &ids, std::size_t begin, std::size_t end,
                             const std::vector<Bounds> &bounds);
};
struct Geometry {
    bool rounded = false;
    Vec3 a{}, b{};
    Real radius{};
    std::vector<Vec3> vertices;
    std::vector<std::vector<std::uint32_t>> faces;
    std::vector<Vec3> normals;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> edges;
    Vec3 center{};
    Vec3 support(Vec3 direction) const;
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
};
Separation separate(const Geometry &a, const Geometry &b);
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
        Pose pose;
        detail::Bounds bounds;
        std::size_t body{};
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
    mutable bool dirty = true;
    mutable detail::Tree tree;
    mutable std::vector<Proxy> proxies;
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
    void impulse(Body &b, Vec3 value, Vec3 point);
    Real effective_mass(const Body &a, const Body *b, Vec3 pa, Vec3 pb, Vec3 axis) const;
    void find_contacts();
    void solve(Real h);
    std::vector<std::size_t> islands() const;
    void propagate_waking();
    void sleep_islands();
    void finish_events(std::vector<Pair> old_pairs);
};
} // namespace fxp_physics
