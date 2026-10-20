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
        tree_queries();
        std::cout << "Geometry caches and refitted tree queries matched their oracles\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
