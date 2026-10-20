#include "../src/detail.hpp"
#include <iostream>

using namespace fxp_physics;
using namespace fxp_physics::detail;
using namespace fxp::literals;

namespace {
void check(bool value, const char *message) {
    if (!value)
        throw std::runtime_error(message);
}

// Original uncached support calculation: keep the subtraction at each query
// so this also checks translations whose fixed-point additions saturate.
Vec3 support_oracle(const Geometry &g, Vec3 direction) {
    std::size_t selected = 0;
    Real best = dot(g.vertices[0] - g.center, direction);
    for (std::size_t i = 1; i < g.vertices.size(); ++i) {
        const Real projection = dot(g.vertices[i] - g.center, direction);
        if (projection > best) {
            best = projection;
            selected = i;
        }
    }
    return g.vertices[selected];
}

void same_separation(const Geometry &a, const Geometry &b) {
    const auto cached = separate(a, b);
    auto fresh_a = a, fresh_b = b;
    fresh_a.edge_directions.clear();
    fresh_b.edge_directions.clear();
    const auto fresh = separate(fresh_a, fresh_b);
    check(cached.distance == fresh.distance && cached.normal == fresh.normal &&
              cached.point_a == fresh.point_a && cached.point_b == fresh.point_b,
          "geometry cache changed separation or tie ordering");
}

void geometry_caches() {
    Random random(918273645);
    const auto box = Shape::box({0.6_r, 0.8_r, 1.1_r});
    const auto sphere = Shape::sphere(0.5_r);
    const auto mesh = Shape::mesh({{-1, 0, -1}, {0, 0.2_r, 1}, {1, 0, -1}}, {{0, 1, 2}});
    for (unsigned i = 0; i < 80; ++i) {
        const Quat rotation(random.unit() * 6_r, unit({random.unit(), random.unit(), random.unit()}));
        const auto source = geometry(box->data(), {{1, -2, 3}, rotation});
        const auto other = geometry(box->data(), {{1.5_r, -2.4_r, 3}, Quat(0.7_r, {0, 1, 0})});
        check(source.edge_directions.empty(), "edge directions were eagerly computed");
        (void)separate(source, geometry(sphere->data(), {}));
        check(source.edge_directions.empty(), "rounded pair unnecessarily computed edge directions");
        same_separation(source, other);
        same_separation(source, geometry(mesh->data(), {}, 0));
        check(!source.edge_directions.empty(), "polyhedron pair did not retain edge directions");
        // Both ordinary motion and saturating arithmetic must invalidate inputs
        // to support and edge-direction queries consistently.
        for (const auto delta :
             {Vec3{0.3_r, -0.4_r, 0.7_r}, Vec3{Real::max(), Real::min(), Real::max()}}) {
            const auto moved = translated(source, delta);
            for (unsigned j = 0; j < 24; ++j) {
                const Vec3 direction{random.unit() - 0.5_r, random.unit() - 0.5_r,
                                     random.unit() - 0.5_r};
                check(moved.support(direction) == support_oracle(moved, direction),
                      "cached support disagreed with uncached projection");
            }
            same_separation(moved, other);
        }
    }
}

void tiny_axes() {
    Geometry a, b;
    a.vertices = a.centered_vertices = {{0, 0, 0}};
    b.vertices = b.centered_vertices = {{1, 0, 0}};
    const auto verify = [&](Vec3 axis) {
        a.normals = {axis};
        const bool rejected = separate(a, b).distance == Real::min();
        check(rejected == (length(axis) < 0.0000001_r),
              "axis rejection changed rounded length threshold");
    };
    // Squared magnitudes immediately below/above the rounded threshold.
    verify({Real::from_raw(3), Real::from_raw(151), Real::from_raw(401)});
    verify({Real::from_raw(3), Real::from_raw(198), Real::from_raw(380)});
    verify({});
    Random random(564738291);
    for (unsigned i = 0; i < 2000; ++i) {
        const auto component = [&]() {
            return Real::from_raw(static_cast<std::int64_t>(random.bounded(1001)) - 500);
        };
        verify({component(), component(), component()});
    }
}

void bounded_separation() {
    Random random(807060504);
    const std::vector<ShapePtr> shapes{
        Shape::box({0.7_r, 0.9_r, 1.1_r}), Shape::box({0.4_r, 0.5_r, 0.6_r}),
        Shape::sphere(0.75_r), Shape::capsule(0.4_r, 0.7_r),
        Shape::mesh({{-1, 0, -1}, {0, 0, 1}, {1, 0, -1}}, {{0, 1, 2}}),
        Shape::heightfield(2, 2, 2_r, {0, 0.3_r, -0.2_r, 0.1_r})};
    const auto random_vector = [&]() {
        return Vec3{random.unit() - 0.5_r, random.unit() - 0.5_r, random.unit() - 0.5_r};
    };
    for (unsigned i = 0; i < 300; ++i) {
        const Pose pa{random_vector() * 8_r, Quat(random.unit() * 6_r, unit(random_vector()))};
        const Pose pb{random_vector() * 8_r, Quat(random.unit() * 6_r, unit(random_vector()))};
        const auto a = geometry(shapes[i % shapes.size()]->data(), pa, 0);
        const auto b = geometry(shapes[(i / shapes.size()) % shapes.size()]->data(), pb, 0);
        const auto exact = separate(a, b);
        for (const auto threshold : {0_r, 0.02_r, 0.2_r, exact.distance,
                                     exact.distance - Real::epsilon(), exact.distance + Real::epsilon()}) {
            const auto bounded = separate(a, b, threshold);
            check((bounded.distance > threshold) == (exact.distance > threshold),
                  "bounded separation changed rejection");
            if (exact.distance <= threshold)
                check(bounded.distance == exact.distance && bounded.normal == exact.normal &&
                          bounded.point_a == exact.point_a && bounded.point_b == exact.point_b,
                      "bounded separation changed accepted contact data");
        }
    }
}

void box_contacts() {
    const auto cube = Shape::box({0.5_r, 0.5_r, 0.5_r});
    const auto a = geometry(cube->data(), {});
    for (const Real height : {0.98_r, 1_r, 1.01_r}) {
        const auto b = geometry(cube->data(), {{0, height, 0}, {}});
        const auto separation = separate(a, b);
        check(separation.normal == Vec3{0, 1, 0} && separation.distance == height - 1_r,
              "aligned box separation was incorrect");
        const auto points = manifold(a, b, separation, 0.02_r);
        check(points.size() == 4, "aligned boxes lost their four-point support patch");
        Vec3 average{};
        for (const auto point : points) {
            check(point.separation == separation.distance &&
                      fxp::abs(point.point.y - height / 2_r) <= Real::epsilon() &&
                      fxp::abs(point.point.x) == 0.5_r && fxp::abs(point.point.z) == 0.5_r,
                  "aligned box contact did not lie at a shared face corner");
            average += point.point / 4_r;
        }
        check(average.x == 0 && average.z == 0 && fxp::abs(average.y - height / 2_r) <= Real::from_raw(2),
              "aligned box support patch was unbalanced");
    }
    const auto rotated = geometry(cube->data(), {{0, 0.99_r, 0}, Quat(0.785398163_r, {0, 1, 0})});
    const auto rotated_separation = separate(a, rotated);
    const auto rotated_points = manifold(a, rotated, rotated_separation, 0.02_r);
    check(rotated_points.size() == 4, "rotated square intersection did not reduce to four contacts");
    for (const auto point : rotated_points)
        for (const auto *box : {&a, &rotated})
            for (const unsigned axis : {0u, 2u})
                check(fxp::abs(dot(point.point - box->center, box->box_axes[axis])) <= 0.500001_r,
                      "rotated box clipping emitted a point outside the shared face patch");

    const auto rod = Shape::box({1, 0.1_r, 0.1_r});
    const auto rod_a = geometry(rod->data(), {{}, Quat(0.7_r, {0, 1, 0})});
    auto rod_b = geometry(rod->data(), {{}, Quat(0.9_r, {0, 0, 1})});
    const auto edge_normal = unit(cross(rod_a.box_axes[0], rod_b.box_axes[0]));
    Real radius{};
    for (unsigned i = 0; i < 3; ++i)
        radius += fxp::abs(dot(rod_a.box_half_axes[i], edge_normal)) +
                  fxp::abs(dot(rod_b.box_half_axes[i], edge_normal));
    rod_b = translated(rod_b, edge_normal * (radius + 0.01_r));
    const auto rod_separation = separate(rod_a, rod_b);
    const auto rod_points = manifold(rod_a, rod_b, rod_separation, 0.02_r);
    check(rod_separation.box_axis == 6 && rod_points.size() == 1 &&
              fxp::abs(rod_separation.distance - 0.01_r) < 0.000001_r,
          "skew rod edges did not produce the expected separated edge contact");
    for (const Geometry *box : std::array<const Geometry *, 2>{&rod_a, &rod_b})
        for (unsigned i = 0; i < 3; ++i)
            check(fxp::abs(dot(rod_points.front().point - box->center, box->box_axes[i])) <=
                      rod->data().half_extents[i] + 0.005001_r,
                  "skew rod contact lies outside the two supporting edges");

    Random random(601204803);
    const auto vector = [&]() {
        return Vec3{random.unit() - 0.5_r, random.unit() - 0.5_r, random.unit() - 0.5_r};
    };
    unsigned edge_contacts = 0, clipped_contacts = 0;
    for (unsigned i = 0; i < 1500; ++i) {
        const auto shape_a = Shape::box({0.1_r + random.unit(), 0.1_r + random.unit(),
                                        0.1_r + random.unit()});
        const auto shape_b = Shape::box({0.1_r + random.unit(), 0.1_r + random.unit(),
                                        0.1_r + random.unit()});
        const auto ga = geometry(shape_a->data(), {vector() * 4_r, Quat(random.unit() * 6_r, unit(vector()))});
        const auto gb = geometry(shape_b->data(), {vector() * 4_r, Quat(random.unit() * 6_r, unit(vector()))});
        auto poly_a = ga, poly_b = gb;
        poly_a.box = poly_b.box = false;
        const auto actual = separate(ga, gb), oracle = separate(poly_a, poly_b);
        // Compare against the existing generic support-vertex SAT, which uses
        // independent projections and edge enumeration. Integer rounding of
        // analytic radii may differ by a few raw units.
        check(fxp::abs(actual.distance - oracle.distance) < 0.00001_r,
              "box SAT disagreed with generic convex SAT");
        const auto reversed = separate(gb, ga);
        check(fxp::abs(actual.distance - reversed.distance) < 0.00001_r,
              "box SAT separation was not symmetric");
        const auto moved_a = translated(ga, {100, -30, 25}), moved_b = translated(gb, {100, -30, 25});
        const auto moved = separate(moved_a, moved_b);
        check(actual.distance == moved.distance && actual.normal == moved.normal &&
                  actual.box_axis == moved.box_axis,
              "box SAT changed under a shared translation");
        if (actual.distance <= 0.02_r) {
            const auto points = manifold(ga, gb, actual, 0.02_r);
            check(!points.empty() && points.size() <= 4, "invalid box contact patch size");
            bool near_face = false;
            for (unsigned axis = 0; axis < 3; ++axis)
                near_face |= fxp::abs(dot(ga.box_axes[axis], actual.normal)) > 0.999_r ||
                             fxp::abs(dot(gb.box_axes[axis], actual.normal)) > 0.999_r;
            if (actual.box_axis >= 6 && !near_face) {
                ++edge_contacts;
                check(points.size() == 1, "edge/edge box contact generated a face patch");
                check(fxp::abs(points.front().separation - actual.distance) < 0.00001_r,
                      "edge contact depth disagreed with its separating axis");
            } else {
                ++clipped_contacts;
                for (const auto point : points)
                    check(point.separation <= 0.02_r, "clipped box contact exceeded margin");
            }
            const auto moved_points = manifold(moved_a, moved_b, moved, 0.02_r);
            check(points.size() == moved_points.size(), "box manifold size changed under translation");
            for (std::size_t j = 0; j < points.size(); ++j)
                check(points[j].point + Vec3{100, -30, 25} == moved_points[j].point &&
                          points[j].separation == moved_points[j].separation,
                      "box manifold changed under a shared translation");
        }
    }
    check(edge_contacts > 30 && clipped_contacts > 30, "random box cases missed contact feature types");

    // Very small relative angles exercise nearly parallel edge cross products.
    for (const Real angle : {Real::epsilon(), 0.000001_r, 0.0001_r, 0.01_r}) {
        const auto b = geometry(cube->data(), {{0, 0.999_r, 0}, Quat(angle, {0, 0, 1})});
        const auto separation = separate(a, b);
        check(separation.distance < 0 && dot(separation.normal, {0, 1, 0}) > 0.999_r,
              "nearly parallel boxes produced an unstable normal");
        const auto points = manifold(a, b, separation, 0.02_r);
        check(!points.empty() && points.size() <= 4, "nearly parallel box clipping failed");
    }
    // Any almost-face cross-axis witness must retain a support patch. This
    // occurs frequently as a resting stack acquires tiny angular corrections.
    const auto tilted_a = geometry(cube->data(), {{}, Quat(0.003_r, unit({1, 2, 3}))});
    const auto tilted_b = geometry(cube->data(), {{0, 0.999_r, 0}, Quat(-0.002_r, unit({3, 2, 1}))});
    auto witness = separate(tilted_a, tilted_b);
    // The X/Z edge cross product is nearly Y: exercise the manifold rule even
    // on backends where another nearly equal face axis wins this particular SAT.
    witness.box_axis = 6 + 0 * 3 + 2;
    witness.normal = -unit(cross(tilted_a.box_axes[0], tilted_b.box_axes[2]));
    check(manifold(tilted_a, tilted_b, witness, 0.02_r).size() == 4,
          "nearly face-aligned edge witness collapsed a resting support patch");
    const auto saturated = translated(a, {Real::max(), 0, 0});
    check(!saturated.box, "saturated box translation retained an invalid OBB cache");
}

void extreme_box_bounds() {
    Random random(8833);
    const auto wide = Shape::box({128, 0.005_r, 128});
    const auto tall = Shape::box({128, 128, 0.005_r});
    const auto rotation = [&]() {
        return Quat(random.unit() * 6_r,
                    unit({random.unit() - 0.5_r, random.unit() - 0.5_r, random.unit() - 0.5_r}));
    };
    for (unsigned i = 0; i < 200; ++i) {
        const Pose pose_a{{8192, -8192, 8192}, rotation()}, pose_b{{8192, -8191, 8192}, rotation()};
        auto a = geometry(wide->data(), pose_a), b = geometry(tall->data(), pose_b);
        const auto box_a = bounds(wide->data(), pose_a), box_b = bounds(tall->data(), pose_b);
        for (unsigned vertex = 0; vertex < 8; ++vertex)
            for (unsigned axis = 0; axis < 3; ++axis)
                check(a.vertices[vertex][axis] >= box_a.min[axis] &&
                          a.vertices[vertex][axis] <= box_a.max[axis] &&
                          b.vertices[vertex][axis] >= box_b.min[axis] &&
                          b.vertices[vertex][axis] <= box_b.max[axis],
                      "box broadphase bounds excluded a rounded vertex at supported range limits");
        const auto actual = separate(a, b);
        a.box = b.box = false;
        // Derive all face normals independently from transformed vertices,
        // including the thin edges, to bound dimension-dependent roundoff.
        for (auto *poly : {&a, &b})
            for (std::size_t face = 0; face < poly->faces.size(); ++face) {
                const auto &indices = poly->faces[face];
                poly->normals[face] = unit(cross(poly->vertices[indices[1]] - poly->vertices[indices[0]],
                                                poly->vertices[indices[2]] - poly->vertices[indices[0]]));
            }
        const auto oracle = separate(a, b);
        check(fxp::abs(actual.distance - oracle.distance) < 0.00001_r,
              "thin wide box axes amplified fixed-point rounding beyond tolerance");
    }
}

Bounds random_bounds(Random &random) {
    const Vec3 center{(random.unit() - 0.5_r) * 100_r, (random.unit() - 0.5_r) * 100_r,
                      (random.unit() - 0.5_r) * 100_r};
    const Vec3 extent{random.unit() * 4_r, random.unit() * 4_r, random.unit() * 4_r};
    return {center - extent, center + extent};
}

void tree_queries() {
    Random random(192837465);
    for (const auto count : {0u, 1u, 2u, 31u, 257u}) {
        std::vector<Bounds> boxes;
        for (unsigned i = 0; i < count; ++i)
            boxes.push_back(random_bounds(random));
        Tree tree;
        tree.build(boxes);
        std::vector<std::uint32_t> actual{999999};
        for (unsigned pass = 0; pass < 5; ++pass) {
            for (unsigned query = 0; query < 120; ++query) {
                const Bounds area = query == 0 ? Bounds{{-100, -100, -100}, {100, 100, 100}}
                                               : random_bounds(random).expanded(8_r);
                std::vector<std::uint32_t> expected;
                for (std::uint32_t i = 0; i < count; ++i)
                    if (boxes[i].intersects(area))
                        expected.push_back(i);
                tree.query(area, actual);
                check(actual == expected, "refitted tree query disagreed with brute force");
                check(tree.query(area) == expected, "allocating tree query disagreed with brute force");
            }
            // Randomizing every leaf makes the original spatial partition poor;
            // refit still has to return the complete, ordered candidate set.
            for (auto &box : boxes)
                box = random_bounds(random);
            tree.refit(boxes);
        }
    }
}
} // namespace

int main() {
    try {
        geometry_caches();
        tiny_axes();
        bounded_separation();
        box_contacts();
        extreme_box_bounds();
        tree_queries();
        std::cout << "Geometry caches and refitted tree queries matched their oracles\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
