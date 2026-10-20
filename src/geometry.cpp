#include "detail.hpp"
#include <numeric>
#include <set>

namespace fxp_physics {
using namespace detail;
Pose compose(const Pose &a, const Pose &b) {
    if (b.rotation.get_internal_vector() == fxp::vec4<Real>{0, 0, 0, 1})
        return {zero(b.position) ? a.position : transform(a, b.position), a.rotation};
    Quat q = a.rotation * b.rotation;
    q.normalize();
    return {transform(a, b.position), q};
}
Pose interpolate(const Pose &a, const Pose &b, Real alpha) {
    require(alpha >= 0 && alpha <= 1, "interpolation alpha must be in [0,1]");
    Quat q;
    q.slerp(alpha, a.rotation, b.rotation);
    q.normalize();
    return {a.position + (b.position - a.position) * alpha, q};
}
std::uint64_t Random::next() {
    std::uint64_t z = (state_ += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}
std::uint32_t Random::bounded(std::uint32_t bound) {
    require(bound != 0, "random bound must be positive");
    const std::uint64_t threshold = (std::uint64_t(0) - bound) % bound;
    for (;;) {
        const auto x = next();
        if (x >= threshold)
            return static_cast<std::uint32_t>(x % bound);
    }
}
Real Random::unit() {
    return Real::from_raw(static_cast<std::int64_t>(next() >> 32));
}

namespace detail {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::invalid_argument(message);
}
void validate_vec(Vec3 v, Real max, const char *message) {
    require(fxp::abs(v.x) <= max && fxp::abs(v.y) <= max && fxp::abs(v.z) <= max, message);
}
Pose validated_pose(Pose p, Real max) {
    validate_vec(p.position, max, "position exceeds supported range");
    const auto q = p.rotation.get_internal_vector();
    require(fxp::abs(q.x) <= 2 && fxp::abs(q.y) <= 2 && fxp::abs(q.z) <= 2 && fxp::abs(q.w) <= 2,
            "invalid rotation");
    require(p.rotation.normalize(), "zero quaternion");
    p.rotation.minimize_rotation_angle();
    return p;
}
void validate_material(const Material &m) {
    require(m.static_friction >= 0 && m.static_friction <= 4 && m.dynamic_friction >= 0 &&
                m.dynamic_friction <= m.static_friction && m.restitution >= 0 && m.restitution <= 1,
            "invalid material");
}
void validate_config(const Config &c) {
    require(c.timestep >= 1_r / 1000_r && c.timestep <= 0.1_r, "timestep must be 0.001..0.1 seconds");
    require(c.substeps >= 1 && c.substeps <= 16 && c.velocity_iterations >= 1 && c.velocity_iterations <= 64,
            "invalid solver budget");
    validate_vec(c.gravity, 100_r, "gravity too large");
    require(c.contact_slop > 0 && c.contact_slop <= 0.1_r && c.contact_margin >= c.contact_slop &&
                c.contact_margin <= 0.2_r,
            "invalid contact tolerance");
    require(c.stabilization >= 0 && c.stabilization <= 1 && c.max_correction_speed > 0 &&
                c.max_correction_speed <= 32,
            "invalid stabilization");
    require(c.restitution_threshold >= 0 && c.restitution_threshold <= 128, "invalid bounce threshold");
    require(c.sleep_linear_speed >= 0 && c.sleep_linear_speed <= 1 && c.sleep_angular_speed >= 0 &&
                c.sleep_angular_speed <= 1 && c.sleep_ticks > 0 && c.sleep_ticks <= 1000000,
            "invalid sleep settings");
    require(c.max_linear_speed > 0 && c.max_linear_speed <= 128 && c.max_angular_speed > 0 &&
                c.max_angular_speed <= 32,
            "invalid speed limits");
}
Vec3 multiply(Inertia m, Vec3 v) {
    return {m.xx * v.x + m.xy * v.y + m.xz * v.z, m.xy * v.x + m.yy * v.y + m.yz * v.z,
            m.xz * v.x + m.yz * v.y + m.zz * v.z};
}
Inertia invert(Inertia m) {
    const Real scale = (std::max)({fxp::abs(m.xx), fxp::abs(m.yy), fxp::abs(m.zz), fxp::abs(m.xy),
                                   fxp::abs(m.xz), fxp::abs(m.yz)});
    require(scale > 0 && scale <= 100000000_r, "invalid inertia scale");
    m = {m.xx / scale, m.yy / scale, m.zz / scale, m.xy / scale, m.xz / scale, m.yz / scale};
    const Real xx = m.yy * m.zz - m.yz * m.yz, yy = m.xx * m.zz - m.xz * m.xz, zz = m.xx * m.yy - m.xy * m.xy;
    const Real xy = m.xz * m.yz - m.xy * m.zz, xz = m.xy * m.yz - m.xz * m.yy, yz = m.xy * m.xz - m.xx * m.yz;
    const Real det = m.xx * xx + m.xy * xy + m.xz * xz;
    require(m.xx > 0 && zz > 0 && det > 0.00000001_r,
            "inertia must be positive definite and well conditioned");
    Inertia out{xx / det / scale, yy / det / scale, zz / det / scale,
                xy / det / scale, xz / det / scale, yz / det / scale};
    require((std::max)({fxp::abs(out.xx), fxp::abs(out.yy), fxp::abs(out.zz)}) <= 1000000_r,
            "inertia too small");
    return out;
}
Bounds transformed(Bounds b, Pose p) {
    if (p.rotation.get_internal_vector() == fxp::vec4<Real>{0, 0, 0, 1})
        return {b.min + p.position, b.max + p.position};
    // Use the exact same quaternion operation as transformed geometry. A matrix
    // extent formula has a different fixed-point rounding path and can shrink a
    // boundary. Identity rotations and spheres have dedicated fast paths.
    Bounds out;
    for (unsigned i = 0; i < 8; ++i)
        out.add(transform(
            p, {(i & 1) ? b.max.x : b.min.x, (i & 2) ? b.max.y : b.min.y, (i & 4) ? b.max.z : b.min.z}));
    return out.expanded(Real::from_raw(16));
}
void Tree::build(const std::vector<Bounds> &bounds) {
    nodes.clear();
    if (bounds.empty())
        return;
    nodes.reserve(bounds.size() * 2);
    std::vector<std::uint32_t> ids(bounds.size());
    std::iota(ids.begin(), ids.end(), 0u);
    build_node(ids, 0, ids.size(), bounds);
}
std::uint32_t Tree::build_node(std::vector<std::uint32_t> &ids, std::size_t begin, std::size_t end,
                               const std::vector<Bounds> &boxes) {
    const auto index = static_cast<std::uint32_t>(nodes.size());
    nodes.emplace_back();
    Bounds total;
    for (auto i = begin; i < end; ++i)
        total.add(boxes[ids[i]]);
    nodes[index].bounds = total;
    if (end - begin == 1) {
        nodes[index].item = ids[begin];
        return index;
    }
    Vec3 extent = total.max - total.min;
    unsigned axis = 0;
    if (extent.y > extent.x)
        axis = 1;
    if (extent.z > extent[axis])
        axis = 2;
    const auto middle = begin + (end - begin) / 2;
    // The comparator has a unique final ID key. nth_element therefore selects
    // the same two sets on every STL; the recursive leaf layout is identical.
    std::nth_element(ids.begin() + static_cast<std::ptrdiff_t>(begin),
                     ids.begin() + static_cast<std::ptrdiff_t>(middle),
                     ids.begin() + static_cast<std::ptrdiff_t>(end), [&](auto a, auto b) {
                         const Real ca = boxes[a].min[axis] + boxes[a].max[axis],
                                    cb = boxes[b].min[axis] + boxes[b].max[axis];
                         return ca != cb ? ca < cb : a < b;
                     });
    const auto left = build_node(ids, begin, middle, boxes), right = build_node(ids, middle, end, boxes);
    nodes[index].left = left;
    nodes[index].right = right;
    return index;
}
void Tree::refit(const std::vector<Bounds> &bounds) {
    // Children follow their parent in build_node's preorder layout.
    for (auto i = nodes.size(); i != 0; --i) {
        auto &node = nodes[i - 1];
        if (node.item != 0xffffffffu)
            node.bounds = bounds[node.item];
        else {
            node.bounds = nodes[node.left].bounds;
            node.bounds.add(nodes[node.right].bounds);
        }
    }
}
void Tree::query(Bounds b, std::vector<std::uint32_t> &out) const {
    out.clear();
    if (nodes.empty())
        return;
    // Median splits bound the depth by 32 for 32-bit item indices. Traversal
    // holds at most one pending sibling per level, plus the current node.
    std::array<std::uint32_t, 33> stack{};
    std::size_t size = 1;
    while (size != 0) {
        const auto i = stack[--size];
        const auto &n = nodes[i];
        if (!n.bounds.intersects(b))
            continue;
        if (n.item != 0xffffffffu)
            out.push_back(n.item);
        else {
            stack[size++] = n.right;
            stack[size++] = n.left;
        }
    }
    std::sort(out.begin(), out.end());
}
std::vector<std::uint32_t> Tree::query(Bounds b) const {
    std::vector<std::uint32_t> out;
    query(b, out);
    return out;
}
static std::vector<std::pair<std::uint32_t, std::uint32_t>>
edge_indices(const std::vector<std::vector<std::uint32_t>> &faces) {
    std::vector<std::pair<std::uint32_t, std::uint32_t>> edges;
    for (const auto &f : faces)
        for (std::size_t i = 0; i < f.size(); ++i) {
            auto a = f[i], b = f[(i + 1) % f.size()];
            if (a > b)
                std::swap(a, b);
            edges.emplace_back(a, b);
        }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
}
static void finish_poly(Geometry &g) {
    g.center = {};
    for (const auto v : g.vertices)
        g.center += v / Real(static_cast<std::int64_t>(g.vertices.size()));
    g.centered_vertices.reserve(g.vertices.size());
    for (const auto v : g.vertices)
        g.centered_vertices.push_back(v - g.center);
    g.normals.reserve(g.faces.size());
    for (const auto &f : g.faces)
        g.normals.push_back(
            unit(cross(g.vertices[f[1]] - g.vertices[f[0]], g.vertices[f[2]] - g.vertices[f[0]])));
}
Geometry geometry(const Shape::Data &s, Pose p, std::uint32_t triangle) {
    Geometry g;
    if (s.type == ShapeType::sphere || s.type == ShapeType::capsule) {
        g.rounded = true;
        g.radius = s.radius;
        if (s.type == ShapeType::sphere)
            g.a = g.b = p.position;
        else {
            g.a = transform(p, {0, -s.half_segment, 0});
            g.b = transform(p, {0, s.half_segment, 0});
        }
        g.center = p.position;
        return g;
    }
    if (s.type == ShapeType::box) {
        // Keep the generic polyhedron representation for box/mesh and rounded
        // pairs, but cache the three OBB directions for the common box pair.
        // Building vertices from the same half axes makes the analytic support
        // calculation and clipped faces describe the same fixed-point box.
        g.box = true;
        g.center = p.position;
        for (unsigned i = 0; i < 3; ++i) {
            // Rotate the unit basis independently. Normalizing a very short
            // half axis would amplify its rounding error for thin, wide boxes.
            Vec3 basis{};
            basis[i] = 1;
            g.box_axes[i] = unit(p.rotation.rotate(basis));
            g.box_half_axes[i] = g.box_axes[i] * s.half_extents[i];
        }
        g.faces = s.faces;
        g.edges = s.edges;
        g.vertices.reserve(8);
        g.centered_vertices.reserve(8);
        for (unsigned i = 0; i < 8; ++i) {
            const Vec3 offset = ((i & 1) ? g.box_half_axes[0] : -g.box_half_axes[0]) +
                                ((i & 2) ? g.box_half_axes[1] : -g.box_half_axes[1]) +
                                ((i & 4) ? g.box_half_axes[2] : -g.box_half_axes[2]);
            g.vertices.push_back(p.position + offset);
            g.centered_vertices.push_back(g.vertices.back() - g.center);
        }
        g.normals.reserve(6);
        for (const auto axis : g.box_axes) {
            g.normals.push_back(-axis);
            g.normals.push_back(axis);
        }
        return g;
    }
    if (s.type == ShapeType::mesh) {
        const auto t = s.mesh_triangles.at(triangle);
        g.vertices = {transform(p, s.vertices[t.a]), transform(p, s.vertices[t.b]),
                      transform(p, s.vertices[t.c])};
        g.faces = {{0, 1, 2}, {2, 1, 0}};
        g.edges = {{0, 1}, {0, 2}, {1, 2}};
    } else if (s.type == ShapeType::heightfield) {
        const auto cell = triangle / 2, x = cell % (s.width - 1), z = cell / (s.width - 1);
        const auto vertex = [&](std::uint32_t dx, std::uint32_t dz) {
            return transform(p, {Real(x + dx) * s.spacing, s.heights[(z + dz) * s.width + x + dx],
                                 Real(z + dz) * s.spacing});
        };
        // Shared diagonal from (0,0) to (1,1), with upward-facing first faces.
        g.vertices = (triangle & 1) ? std::vector<Vec3>{vertex(0, 0), vertex(1, 1), vertex(1, 0)}
                                    : std::vector<Vec3>{vertex(0, 0), vertex(0, 1), vertex(1, 1)};
        g.faces = {{0, 1, 2}, {2, 1, 0}};
        g.edges = {{0, 1}, {0, 2}, {1, 2}};
    } else {
        g.faces = s.faces;
        g.edges = s.edges;
        g.vertices.reserve(s.vertices.size());
        for (auto v : s.vertices)
            g.vertices.push_back(transform(p, v));
    }
    finish_poly(g);
    return g;
}
Bounds bounds(const Shape::Data &s, Pose p) {
    if (s.type == ShapeType::sphere) {
        const Vec3 r{s.radius, s.radius, s.radius};
        return {p.position - r, p.position + r};
    }
    if (s.type == ShapeType::box) {
        // Use the same rounded basis as geometry(). Summed absolute half axes
        // enclose all eight corners exactly, including thin rotated boxes.
        Vec3 extent{};
        for (unsigned i = 0; i < 3; ++i) {
            Vec3 basis{};
            basis[i] = 1;
            extent += abs(unit(p.rotation.rotate(basis)) * s.half_extents[i]);
        }
        return {p.position - extent, p.position + extent};
    }
    return transformed(s.bounds, p);
}
std::vector<std::uint32_t> triangles(const Shape::Data &s, Pose p, Bounds area) {
    // Convex primitives have no triangle tree: avoid transforming its query
    // bounds for every broadphase pair.
    if (s.type != ShapeType::mesh && s.type != ShapeType::heightfield)
        return {0xffffffffu};
    const Bounds local = transformed(area, relative(p, Pose{}));
    if (s.type == ShapeType::mesh)
        return s.tree.query(local);
    if (!s.bounds.intersects(local))
        return {};
    const auto cell = [&](Real v, std::uint32_t maximum) {
        return static_cast<std::uint32_t>(
            clamp(fxp::floor(v / s.spacing), 0, Real(maximum)).to_integer<std::int64_t>());
    };
    const auto x0 = cell(local.min.x, s.width - 2), x1 = cell(local.max.x, s.width - 2);
    const auto z0 = cell(local.min.z, s.depth - 2), z1 = cell(local.max.z, s.depth - 2);
    std::vector<std::uint32_t> out;
    for (auto z = z0; z <= z1; ++z)
        for (auto x = x0; x <= x1; ++x) {
            const auto base = (z * (s.width - 1) + x) * 2;
            out.push_back(base);
            out.push_back(base + 1);
        }
    return out;
}
Vec3 Geometry::support(Vec3 d) const {
    if (rounded)
        return (dot(a, d) >= dot(b, d) ? a : b) + unit(d) * radius;
    Vec3 best = vertices.front();
    Real value = dot(centered_vertices.front(), d);
    for (std::size_t i = 1; i < vertices.size(); ++i) {
        const Real next = dot(centered_vertices[i], d);
        if (next > value) {
            value = next;
            best = vertices[i];
        }
    }
    return best;
}
const std::vector<Vec3> &Geometry::directions() const {
    if (edge_directions.empty()) {
        edge_directions.reserve(edges.size());
        for (const auto &edge : edges) {
            const Vec3 direction = unit(vertices[edge.second] - vertices[edge.first]);
            if (std::find(edge_directions.begin(), edge_directions.end(), direction) == edge_directions.end())
                edge_directions.push_back(direction);
        }
    }
    return edge_directions;
}
// Most primitive pairs need only a few distinct axes. Keep those on the stack,
// then retain logarithmic lookup for hulls with many edge directions.
class AxisSet {
    std::array<Vec3, 32> small_{};
    std::size_t size_ = 0;
    std::set<std::tuple<Real, Real, Real>> large_;

  public:
    bool insert(Vec3 axis) {
        if (!large_.empty())
            return large_.emplace(axis.x, axis.y, axis.z).second;
        const auto end = small_.begin() + static_cast<std::ptrdiff_t>(size_);
        if (std::find(small_.begin(), end, axis) != end)
            return false;
        if (size_ < small_.size()) {
            small_[size_++] = axis;
            return true;
        }
        for (const auto v : small_)
            large_.emplace(v.x, v.y, v.z);
        return large_.emplace(axis.x, axis.y, axis.z).second;
    }
};
static bool tiny_axis(Vec3 axis) {
    // For integer squared magnitude S, rounded_sqrt(S) < t exactly when
    // S < t*t - t + 1. Compare raw squares, preserving length()'s rounding
    // without taking a square root or rounding each component's square.
    constexpr auto t = static_cast<std::uint64_t>((0.0000001_r).raw_value());
    constexpr auto limit = t * t - t + 1;
    std::uint64_t sum = 0;
    for (unsigned i = 0; i < 3; ++i) {
        const auto v = fxp::detail::magnitude(axis[i].raw_value());
        if (v >= t)
            return false;
        sum += v * v;
    }
    return sum < limit;
}
static Vec3 closest_segment(Vec3 a, Vec3 b, Vec3 p) {
    const Vec3 d = b - a;
    const Real denom = dot(d, d);
    return denom > 0 ? a + d * clamp(dot(p - a, d) / denom, 0, 1) : a;
}
static std::pair<Vec3, Vec3> closest_segments(Vec3 a, Vec3 b, Vec3 c, Vec3 d) {
    // Normalize directions before the 2x2 solve to avoid fourth-power lengths.
    const Real la = length(b - a), lb = length(d - c);
    if (la <= 0.000001_r)
        return {a, closest_segment(c, d, a)};
    if (lb <= 0.000001_r)
        return {closest_segment(a, b, c), c};
    const Vec3 u = (b - a) / la, v = (d - c) / lb, r = a - c;
    const Real uv = dot(u, v), ur = dot(u, r), vr = dot(v, r), denom = 1_r - uv * uv;
    Real s = denom > 0.000001_r ? clamp((uv * vr - ur) / denom, 0, la) : 0_r;
    Real t = uv * s + vr;
    if (t < 0) {
        t = 0;
        s = clamp(-ur, 0, la);
    } else if (t > lb) {
        t = lb;
        s = clamp(uv * lb - ur, 0, la);
    }
    return {a + u * s, c + v * t};
}
static Vec3 box_support(const Geometry &box, Vec3 direction) {
    unsigned vertex = 0;
    for (unsigned i = 0; i < 3; ++i)
        if (dot(box.box_half_axes[i], direction) > 0)
            vertex |= 1u << i;
    return box.vertices[vertex];
}
static Separation separate_boxes(const Geometry &a, const Geometry &b, Real rejection_threshold) {
    Separation result;
    const Vec3 delta = b.center - a.center;
    const auto axis = [&](Vec3 normal, std::int32_t feature) {
        Real radius_a{}, radius_b{};
        for (unsigned i = 0; i < 3; ++i) {
            radius_a += fxp::abs(dot(a.box_half_axes[i], normal));
            radius_b += fxp::abs(dot(b.box_half_axes[i], normal));
        }
        const Real center_distance = dot(delta, normal);
        const Real distance = fxp::abs(center_distance) - radius_a - radius_b;
        if (distance > result.distance) {
            result.distance = distance;
            result.normal = center_distance < 0 ? -normal : normal;
            result.box_axis = feature;
        }
        return distance > rejection_threshold;
    };
    for (unsigned i = 0; i < 3; ++i)
        if (axis(a.box_axes[i], static_cast<std::int32_t>(i)))
            return result;
    // Identically oriented boxes need only three axes. This covers undeformed
    // stacks and avoids testing twelve exact duplicates.
    if (a.box_axes != b.box_axes) {
        for (unsigned i = 0; i < 3; ++i)
            if (axis(b.box_axes[i], static_cast<std::int32_t>(3 + i)))
                return result;
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j) {
                const Vec3 normal = cross(a.box_axes[i], b.box_axes[j]);
                if (!tiny_axis(normal) &&
                    axis(unit(normal), static_cast<std::int32_t>(6 + i * 3 + j)))
                    return result;
            }
    }
    result.point_a = box_support(a, result.normal);
    result.point_b = box_support(b, -result.normal);
    return result;
}
Separation separate(const Geometry &a, const Geometry &b, Real rejection_threshold) {
    if (a.box && b.box)
        return separate_boxes(a, b, rejection_threshold);
    Separation result;
    if (a.rounded && b.rounded) {
        const auto closest = closest_segments(a.a, a.b, b.a, b.b);
        const Vec3 delta = closest.second - closest.first;
        result.normal = delta;
        if (!fxp::normalize(&result.normal))
            result.normal = unit(b.center - a.center);
        result.distance = length(delta) - a.radius - b.radius;
        result.point_a = closest.first + result.normal * a.radius;
        result.point_b = closest.second - result.normal * b.radius;
        return result;
    }
    // Exact duplicate axes cannot improve the strict maximum below. Fold
    // opposite signs too, preserving the first tested orientation and all ties.
    AxisSet tested_axes;
    const auto axis = [&](Vec3 n) {
        if (zero(n))
            return false;
        const Vec3 key = (n.x < 0 || (n.x == 0 && (n.y < 0 || (n.y == 0 && n.z < 0)))) ? -n : n;
        if (!tested_axes.insert(key))
            return false;
        if (tiny_axis(n))
            return false;
        n = unit(n);
        const Vec3 ap = a.support(n), am = a.support(-n), bp = b.support(n), bm = b.support(-n);
        const Real positive = dot(bm - ap, n), negative = dot(am - bp, n);
        const bool flip = negative > positive;
        const Real separation = flip ? negative : positive;
        if (separation > result.distance) {
            result.distance = separation;
            result.normal = flip ? -n : n;
            result.point_a = flip ? am : ap;
            result.point_b = flip ? bp : bm;
        }
        return separation > rejection_threshold;
    };
    // Contact/overlap queries only need a separating witness for rejected
    // pairs. Accepted pairs still test every axis and retain identical ties.
    // Casts use the default threshold and require the full maximum distance.
    for (auto n : a.normals)
        if (axis(n))
            return result;
    for (auto n : b.normals)
        if (axis(n))
            return result;
    if (!a.rounded && !b.rounded) {
        // Preserve first-occurrence order and ties. Cache lazily because
        // rounded-vs-poly contacts never need these normalized directions.
        const auto &da = a.directions(), &db = b.directions();
        for (const auto ea : da)
            for (const auto eb : db)
                if (axis(cross(ea, eb)))
                    return result;
    } else {
        const auto &round = a.rounded ? a : b;
        const auto &poly = a.rounded ? b : a;
        for (auto v : poly.vertices)
            if (axis(v - closest_segment(round.a, round.b, v)))
                return result;
        for (auto edge : poly.edges) {
            const Vec3 p = poly.vertices[edge.first], q = poly.vertices[edge.second];
            const auto closest = closest_segments(round.a, round.b, p, q);
            if (axis(closest.second - closest.first) || axis(cross(round.b - round.a, q - p)) ||
                axis(closest_segment(p, q, round.a) - round.a) ||
                axis(closest_segment(p, q, round.b) - round.b))
                return result;
        }
    }
    return result;
}
static std::vector<Vec3> clip(std::vector<Vec3> polygon, Vec3 origin, Vec3 normal) {
    if (polygon.empty())
        return {};
    std::vector<Vec3> out;
    Vec3 previous = polygon.back();
    Real pd = dot(previous - origin, normal);
    for (auto current : polygon) {
        const Real cd = dot(current - origin, normal);
        if ((pd <= 0) != (cd <= 0))
            out.push_back(previous + (current - previous) * (pd / (pd - cd)));
        if (cd <= 0)
            out.push_back(current);
        previous = current;
        pd = cd;
    }
    return out;
}
static void finish_manifold(std::vector<ContactPoint> &out) {
    std::sort(out.begin(), out.end(), [](const auto &x, const auto &y) {
        return std::tie(x.point.x, x.point.y, x.point.z, x.separation) <
               std::tie(y.point.x, y.point.y, y.point.z, y.separation);
    });
    // Approximate equality is not transitive. Compact in the specified sorted
    // order, without allocating a second vector for every contact patch.
    std::size_t distinct = 0;
    for (const auto &point : out) {
        bool duplicate = false;
        for (std::size_t i = 0; i < distinct; ++i)
            if (!length_greater_than(point.point - out[i].point, 0.0001_r - Real::epsilon())) {
                duplicate = true;
                break;
            }
        if (!duplicate)
            out[distinct++] = point;
    }
    out.resize(distinct);
    if (out.size() > 4) {
        // Deepest point, then maximum distance from the chosen set. Ties retain
        // lexicographic order. This keeps the contact patch spatially distributed.
        std::array<ContactPoint, 4> reduced{};
        auto deepest = std::min_element(
            out.begin(), out.end(), [](const auto &x, const auto &y) { return x.separation < y.separation; });
        reduced[0] = *deepest;
        out.erase(deepest);
        for (unsigned count = 1; count < reduced.size(); ++count) {
            std::size_t best = 0;
            Real distance = -1_r;
            for (std::size_t i = 0; i < out.size(); ++i) {
                Real minimum_distance = Real::max();
                for (unsigned j = 0; j < count; ++j)
                    minimum_distance = (std::min)(minimum_distance, length(out[i].point - reduced[j].point));
                if (minimum_distance > distance) {
                    distance = minimum_distance;
                    best = i;
                }
            }
            reduced[count] = out[best];
            out.erase(out.begin() + static_cast<std::ptrdiff_t>(best));
        }
        out.assign(reduced.begin(), reduced.end());
    }
}
static std::pair<Vec3, Vec3> supporting_box_edge(const Geometry &box, unsigned axis, Vec3 direction) {
    Vec3 center = box.center;
    for (unsigned i = 0; i < 3; ++i)
        if (i != axis)
            center += dot(box.box_half_axes[i], direction) > 0 ? box.box_half_axes[i] : -box.box_half_axes[i];
    return {center - box.box_half_axes[axis], center + box.box_half_axes[axis]};
}
static std::vector<ContactPoint> generic_box_manifold(const Geometry &a, const Geometry &b, Separation s,
                                                     Real margin) {
    auto generic = a;
    generic.box = false;
    return manifold(generic, b, s, margin);
}
static std::vector<ContactPoint> box_manifold(const Geometry &a, const Geometry &b, Separation s,
                                             Real margin) {
    std::vector<ContactPoint> out;
    if (s.box_axis >= 6) {
        // An edge cross product can win by a few raw units while its normal
        // is almost parallel to a face. Keep the full SAT witness for casts,
        // but preserve a broad support patch for resting/slightly tilted
        // boxes instead of balancing a tower on one noisy edge point.
        Real alignment = 0.999_r;
        for (unsigned i = 0; i < 6; ++i) {
            const auto &box = i < 3 ? a : b;
            const Real next = fxp::abs(dot(box.box_axes[i % 3], s.normal));
            if (next > alignment) {
                alignment = next;
                s.box_axis = static_cast<std::int32_t>(i);
            }
        }
    }
    if (s.box_axis >= 6) {
        const unsigned axis_a = static_cast<unsigned>(s.box_axis - 6) / 3,
                       axis_b = static_cast<unsigned>(s.box_axis - 6) % 3;
        const auto edge_a = supporting_box_edge(a, axis_a, s.normal),
                   edge_b = supporting_box_edge(b, axis_b, -s.normal);
        const auto closest = closest_segments(edge_a.first, edge_a.second, edge_b.first, edge_b.second);
        out.push_back({(closest.first + closest.second) / 2_r, s.distance});
        return out;
    }
    const bool use_a = s.box_axis < 3;
    const auto &reference = use_a ? a : b;
    const auto &incident = use_a ? b : a;
    const unsigned ref_axis = static_cast<unsigned>(s.box_axis) % 3;
    const Vec3 toward_incident = use_a ? s.normal : -s.normal;
    const Vec3 normal = dot(reference.box_axes[ref_axis], toward_incident) > 0
                            ? reference.box_axes[ref_axis] : -reference.box_axes[ref_axis];
    const Vec3 origin = reference.center +
                        (dot(reference.box_axes[ref_axis], normal) > 0
                             ? reference.box_half_axes[ref_axis]
                             : -reference.box_half_axes[ref_axis]);
    unsigned inc_axis = 0;
    Real alignment = fxp::abs(dot(incident.box_axes[0], normal));
    for (unsigned i = 1; i < 3; ++i) {
        const Real next = fxp::abs(dot(incident.box_axes[i], normal));
        if (next > alignment) {
            inc_axis = i;
            alignment = next;
        }
    }
    const unsigned inc_face = inc_axis * 2 + (dot(incident.box_axes[inc_axis], normal) < 0 ? 1u : 0u);
    // A quadrilateral clipped against four planes has at most eight vertices.
    // Two fixed buffers avoid four allocations per face contact.
    std::array<Vec3, 8> polygon{}, temporary{};
    std::size_t count = 4;
    for (unsigned i = 0; i < 4; ++i)
        polygon[i] = incident.vertices[incident.faces[inc_face][i]];
    for (unsigned i = 0; i < 3 && count != 0; ++i) {
        if (i == ref_axis)
            continue;
        for (const bool positive : {false, true}) {
            const Vec3 side = positive ? reference.box_axes[i] : -reference.box_axes[i];
            const Vec3 plane = reference.center +
                               (positive ? reference.box_half_axes[i] : -reference.box_half_axes[i]);
            std::size_t next_count = 0;
            Vec3 previous = polygon[count - 1];
            Real pd = dot(previous - plane, side);
            for (std::size_t vertex = 0; vertex < count; ++vertex) {
                const Vec3 current = polygon[vertex];
                const Real cd = dot(current - plane, side);
                if ((pd <= 0) != (cd <= 0)) {
                    // Rounding at nearly coincident edges can make the
                    // polygon slightly nonconvex. Retain a safe generic
                    // fallback if that invalidates the eight-vertex bound.
                    if (next_count == temporary.size())
                        return generic_box_manifold(a, b, s, margin);
                    temporary[next_count++] = previous + (current - previous) * (pd / (pd - cd));
                }
                if (cd <= 0) {
                    if (next_count == temporary.size())
                        return generic_box_manifold(a, b, s, margin);
                    temporary[next_count++] = current;
                }
                previous = current;
                pd = cd;
            }
            count = next_count;
            polygon.swap(temporary);
            if (count == 0)
                break;
        }
    }
    out.reserve(count == 0 ? 1 : count);
    for (std::size_t i = 0; i < count; ++i) {
        const Real distance = dot(polygon[i] - origin, normal);
        if (distance <= margin)
            out.push_back({polygon[i] - normal * (distance / 2_r), distance});
    }
    if (out.empty())
        out.push_back({(s.point_a + s.point_b) / 2_r, s.distance});
    finish_manifold(out);
    return out;
}
std::vector<ContactPoint> manifold(const Geometry &a, const Geometry &b, Separation s, Real margin) {
    if (a.box && b.box && s.box_axis >= 0)
        return box_manifold(a, b, s, margin);
    std::vector<ContactPoint> out;
    if (a.rounded || b.rounded) {
        if (a.rounded && b.rounded)
            return {{(s.point_a + s.point_b) / 2_r, s.distance}};
        const auto &rounded = a.rounded ? a : b;
        const auto &poly = a.rounded ? b : a;
        const Vec3 n = a.rounded ? s.normal : -s.normal;
        // Project segment endpoints onto the opposing support plane. Keep both
        // only when the segment is tangent and the projected points are in it.
        const Vec3 plane = poly.support(-n);
        for (const auto end : {rounded.a, rounded.b}) {
            const Vec3 surface = end + n * rounded.radius;
            const Real sep = dot(plane - surface, n);
            if (fxp::abs(sep - s.distance) > 0.002_r)
                continue;
            const Vec3 on_plane = surface + n * sep;
            bool inside = true;
            for (std::size_t f = 0; f < poly.faces.size(); ++f)
                if (dot(on_plane - poly.vertices[poly.faces[f][0]], poly.normals[f]) > 0.002_r) {
                    inside = false;
                    break;
                }
            if (inside && (out.empty() || length(out.front().point - (surface + on_plane) / 2_r) > 0.001_r))
                out.push_back({(surface + on_plane) / 2_r, sep});
        }
        if (out.empty()) {
            // For rounded-edge contacts, recover closest points instead of an
            // arbitrary extremal vertex (which would create artificial torque).
            Real best = Real::max();
            Vec3 ra{}, pb{};
            for (auto e : poly.edges) {
                const auto c =
                    closest_segments(rounded.a, rounded.b, poly.vertices[e.first], poly.vertices[e.second]);
                const Real dist = length(c.second - c.first);
                if (dist < best) {
                    best = dist;
                    ra = c.first;
                    pb = c.second;
                }
            }
            out.push_back({(ra + n * rounded.radius + pb) / 2_r, s.distance});
        }
        return out;
    }
    std::size_t fa = 0, fb = 0;
    Real align_a = dot(a.normals[fa], s.normal), align_b = dot(b.normals[fb], -s.normal);
    for (std::size_t i = 1; i < a.faces.size(); ++i) {
        const Real alignment = dot(a.normals[i], s.normal);
        if (alignment > align_a) {
            fa = i;
            align_a = alignment;
        }
    }
    for (std::size_t i = 1; i < b.faces.size(); ++i) {
        const Real alignment = dot(b.normals[i], -s.normal);
        if (alignment > align_b) {
            fb = i;
            align_b = alignment;
        }
    }
    if ((std::max)(align_a, align_b) > 0.999_r) {
        const bool use_a = align_a >= align_b;
        const auto &reference = use_a ? a : b;
        const auto &incident = use_a ? b : a;
        const auto ref_face = use_a ? fa : fb;
        const Vec3 normal = reference.normals[ref_face];
        std::size_t inc_face = 0;
        Real inc_alignment = dot(incident.normals[inc_face], normal);
        for (std::size_t i = 1; i < incident.faces.size(); ++i) {
            const Real alignment = dot(incident.normals[i], normal);
            if (alignment < inc_alignment) {
                inc_face = i;
                inc_alignment = alignment;
            }
        }
        std::vector<Vec3> polygon;
        for (auto i : incident.faces[inc_face])
            polygon.push_back(incident.vertices[i]);
        const auto &face = reference.faces[ref_face];
        for (std::size_t i = 0; i < face.size(); ++i) {
            const Vec3 start = reference.vertices[face[i]],
                       end = reference.vertices[face[(i + 1) % face.size()]];
            polygon = clip(std::move(polygon), start, unit(cross(end - start, normal)));
        }
        const Vec3 origin = reference.vertices[face[0]];
        for (auto p : polygon) {
            const Real distance = dot(p - origin, normal);
            if (distance <= margin)
                out.push_back({p - normal * (distance / 2_r), distance});
        }
    }
    if (out.empty()) {
        // Edge-edge contacts use closest points on the supporting edges.
        Real best = Real::max();
        Vec3 point = (s.point_a + s.point_b) / 2_r;
        for (auto ea : a.edges)
            for (auto eb : b.edges) {
                const auto c = closest_segments(a.vertices[ea.first], a.vertices[ea.second],
                                                b.vertices[eb.first], b.vertices[eb.second]);
                const Real dist = length(c.second - c.first);
                if (dist < best) {
                    best = dist;
                    point = (c.first + c.second) / 2_r;
                }
            }
        out.push_back({point, s.distance});
    }
    finish_manifold(out);
    return out;
}
Geometry translated(Geometry g, Vec3 delta) {
    if (g.box) {
        // Saturation destroys the parallelepiped: use its actual vertices with
        // the generic path instead of applying an invalid analytic radius.
        for (const auto vertex : g.vertices)
            for (unsigned i = 0; i < 3; ++i)
                if ((delta[i] > 0 && vertex[i] > Real::max() - delta[i]) ||
                    (delta[i] < 0 && vertex[i] < Real::min() - delta[i]))
                    g.box = false;
    }
    g.a += delta;
    g.b += delta;
    g.center += delta;
    for (auto &v : g.vertices)
        v += delta;
    for (std::size_t i = 0; i < g.vertices.size(); ++i)
        g.centered_vertices[i] = g.vertices[i] - g.center;
    // Recompute after translation, including any saturated additions.
    g.edge_directions.clear();
    return g;
}
} // namespace detail

Shape::Shape(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
ShapeType Shape::type() const noexcept {
    return data_->type;
}
const Shape::Data &Shape::data() const noexcept {
    return *data_;
}
ShapePtr Shape::box(Vec3 half) {
    validate_vec(half, 128_r, "box too large");
    require(half.x >= 0.005_r && half.y >= 0.005_r && half.z >= 0.005_r, "box too small");
    auto d = std::make_shared<Data>();
    d->type = ShapeType::box;
    d->half_extents = half;
    d->bounds = {-half, half};
    for (unsigned i = 0; i < 8; ++i)
        d->vertices.push_back(
            {(i & 1) ? half.x : -half.x, (i & 2) ? half.y : -half.y, (i & 4) ? half.z : -half.z});
    d->faces = {{0, 4, 6, 2}, {1, 3, 7, 5}, {0, 1, 5, 4}, {2, 6, 7, 3}, {0, 2, 3, 1}, {4, 5, 7, 6}};
    d->edges = edge_indices(d->faces);
    return ShapePtr(new Shape(d));
}
ShapePtr Shape::sphere(Real radius) {
    return capsule(radius, 0);
}
ShapePtr Shape::capsule(Real radius, Real half_segment) {
    require(radius >= 0.005_r && radius <= 128 && half_segment >= 0 && half_segment + radius <= 128,
            "invalid capsule dimensions");
    auto d = std::make_shared<Data>();
    d->type = half_segment == 0 ? ShapeType::sphere : ShapeType::capsule;
    d->radius = radius;
    d->half_segment = half_segment;
    const Vec3 extent{radius, radius + half_segment, radius};
    d->bounds = {-extent, extent};
    return ShapePtr(new Shape(d));
}
ShapePtr Shape::convex(std::vector<Vec3> vertices, std::vector<std::vector<std::uint32_t>> faces) {
    require(vertices.size() >= 4 && vertices.size() <= 128 && faces.size() >= 4 && faces.size() <= 256,
            "invalid convex size");
    auto d = std::make_shared<Data>();
    d->type = ShapeType::convex;
    d->vertices = std::move(vertices);
    d->faces = std::move(faces);
    for (auto v : d->vertices) {
        validate_vec(v, 128, "convex vertex too large");
        d->bounds.add(v);
    }
    Vec3 center{};
    for (auto v : d->vertices)
        center += v / Real(static_cast<std::int64_t>(d->vertices.size()));
    std::vector<std::pair<std::uint32_t, std::uint32_t>> edges;
    std::vector<bool> used(d->vertices.size(), false);
    for (const auto &f : d->faces) {
        require(f.size() >= 3 && f.size() <= 128, "invalid face size");
        for (auto i : f) {
            require(i < d->vertices.size(), "convex index out of bounds");
            used[i] = true;
        }
        const Vec3 p = d->vertices[f[0]], raw = cross(d->vertices[f[1]] - p, d->vertices[f[2]] - p);
        require(length(raw) > 0.000001_r, "degenerate face");
        const Vec3 n = unit(raw);
        require(dot(center - p, n) < -0.000001_r, "convex must have outward faces and nonzero volume");
        for (auto v : d->vertices)
            require(dot(v - p, n) <= 0.00001_r, "nonconvex hull");
        for (std::size_t i = 0; i < f.size(); ++i) {
            const auto a = f[i], b = f[(i + 1) % f.size()];
            require(a != b && fxp::abs(dot(d->vertices[a] - p, n)) <= 0.00001_r, "invalid planar face");
            const Vec3 side = cross(d->vertices[b] - d->vertices[a], n);
            for (auto j : f)
                require(dot(d->vertices[j] - d->vertices[a], side) <= 0.00001_r, "nonconvex face");
            edges.emplace_back(a, b);
        }
    }
    require(std::all_of(used.begin(), used.end(), [](bool x) { return x; }), "unused hull vertex");
    std::sort(edges.begin(), edges.end());
    require(std::adjacent_find(edges.begin(), edges.end()) == edges.end(), "duplicate directed edge");
    for (auto e : edges)
        require(std::binary_search(edges.begin(), edges.end(), std::make_pair(e.second, e.first)),
                "hull must be closed");
    d->edges = edge_indices(d->faces);
    return ShapePtr(new Shape(d));
}
ShapePtr Shape::mesh(std::vector<Vec3> vertices, std::vector<Triangle> triangles) {
    require(!vertices.empty() && vertices.size() <= 1000000 && !triangles.empty() &&
                triangles.size() <= 1000000,
            "invalid mesh size");
    auto d = std::make_shared<Data>();
    d->type = ShapeType::mesh;
    d->vertices = std::move(vertices);
    d->mesh_triangles = std::move(triangles);
    for (auto v : d->vertices) {
        validate_vec(v, 8192, "mesh vertex too large");
        d->bounds.add(v);
    }
    std::vector<Bounds> boxes;
    for (auto t : d->mesh_triangles) {
        require(t.a < d->vertices.size() && t.b < d->vertices.size() && t.c < d->vertices.size(),
                "triangle index out of bounds");
        const Vec3 a = d->vertices[t.a], b = d->vertices[t.b], c = d->vertices[t.c];
        require(length(cross(b - a, c - a)) > 0.000001_r, "degenerate triangle");
        Bounds box;
        box.add(a);
        box.add(b);
        box.add(c);
        boxes.push_back(box);
    }
    d->tree.build(boxes);
    return ShapePtr(new Shape(d));
}
ShapePtr Shape::heightfield(std::uint32_t width, std::uint32_t depth, Real spacing,
                            std::vector<Real> heights) {
    require(width >= 2 && depth >= 2 && width <= 4096 && depth <= 4096 &&
                std::uint64_t(width) * depth <= 1000000,
            "invalid heightfield size");
    require(heights.size() == std::uint64_t(width) * depth && spacing >= 0.01_r && spacing <= 128,
            "invalid height samples/spacing");
    require(Real(width - 1) * spacing <= 8192 && Real(depth - 1) * spacing <= 8192,
            "heightfield exceeds world range");
    auto d = std::make_shared<Data>();
    d->type = ShapeType::heightfield;
    d->width = width;
    d->depth = depth;
    d->spacing = spacing;
    d->heights = std::move(heights);
    Real low = Real::max(), high = Real::min();
    for (auto h : d->heights) {
        require(fxp::abs(h) <= 8192, "height exceeds world range");
        low = (std::min)(low, h);
        high = (std::max)(high, h);
    }
    d->bounds = {{0, low, 0}, {Real(width - 1) * spacing, high, Real(depth - 1) * spacing}};
    return ShapePtr(new Shape(d));
}
} // namespace fxp_physics
