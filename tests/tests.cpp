#include <fxp_physics/physics.hpp>
#include "../src/detail.hpp"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace fxp_physics;
namespace {
void check(bool value, const std::string &message) {
    if (!value)
        throw std::runtime_error(message);
}
template <class F> void rejects(F function, const char *message) {
    bool caught = false;
    try {
        function();
    } catch (const std::exception &) {
        caught = true;
    }
    check(caught, message);
}
bool near(Real a, Real b, Real tolerance = 0.02_r) {
    return fxp::abs(a - b) <= tolerance;
}
BodyId add(World &w, ShapePtr shape, Vec3 position, BodyType type = BodyType::dynamic, bool sensor = false,
           std::uint8_t layer = 0) {
    BodyDesc b;
    b.type = type;
    b.pose.position = position;
    const auto id = w.create_body(b);
    ColliderDesc c;
    c.body = id;
    c.shape = std::move(shape);
    c.sensor = sensor;
    c.layer = layer;
    w.create_collider(c);
    return id;
}
void ground(World &w) {
    add(w, Shape::box({20, 0.5_r, 20}), {0, -0.5_r, 0}, BodyType::static_body);
}
void smoke() {
    World w;
    ground(w);
    const auto id = add(w, Shape::sphere(0.5_r), {0, 3, 0});
    for (int i = 0; i < 240; ++i)
        w.step();
    const auto b = w.body(id);
    check(near(b.properties.pose.position.y, 0.5_r, 0.04_r),
          "sphere did not rest on ground: " + b.properties.pose.position.y.to_string());
    check(b.sleeping, "sphere did not sleep");
    World copy;
    copy.load(w.save());
    check(copy.save() == w.save(), "save roundtrip differs");
}
ShapePtr hull() {
    return Shape::convex(
        {{-0.5_r, -0.5_r, -0.5_r},
         {0.5_r, -0.5_r, -0.5_r},
         {-0.5_r, 0.5_r, -0.5_r},
         {0.5_r, 0.5_r, -0.5_r},
         {-0.5_r, -0.5_r, 0.5_r},
         {0.5_r, -0.5_r, 0.5_r},
         {-0.5_r, 0.5_r, 0.5_r},
         {0.5_r, 0.5_r, 0.5_r}},
        {{0, 4, 6, 2}, {1, 3, 7, 5}, {0, 1, 5, 4}, {2, 6, 7, 3}, {0, 2, 3, 1}, {4, 5, 7, 6}});
}
ShapePtr plane_mesh() {
    return Shape::mesh({{-10, 0, -10}, {10, 0, -10}, {10, 0, 10}, {-10, 0, 10}}, {{0, 2, 1}, {0, 3, 2}});
}
void rng_and_validation() {
    Random rng(0);
    check(rng.next() == UINT64_C(0xe220a8397b1dcdaf), "RNG golden 1");
    check(rng.next() == UINT64_C(0x6e789e6aa1b965f4), "RNG golden 2");
    Random copy;
    copy.restore(rng.state());
    for (int i = 0; i < 1000; ++i) {
        check(rng.next() == copy.next(), "RNG restore");
        check(rng.bounded(7) < 7, "bounded RNG");
        copy.bounded(7);
    }
    rejects([&] { rng.bounded(0); }, "zero RNG bound accepted");
    rejects([] { Shape::box({0, 1, 1}); }, "zero box accepted");
    rejects(
        [] {
            Shape::convex({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
                          {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}});
        },
        "bad winding accepted");
    rejects([] { Shape::mesh({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}}, {{0, 1, 2}}); }, "degenerate mesh accepted");
    World w;
    const auto b = w.create_body();
    ColliderDesc c;
    c.body = b;
    c.shape = plane_mesh();
    rejects([&] { w.create_collider(c); }, "dynamic concave accepted");
    Config config;
    config.substeps = 0;
    rejects([&] { World bad(config); }, "zero substeps accepted");
    BodyDesc bad;
    bad.mass = 0;
    rejects([&] { w.create_body(bad); }, "zero dynamic mass accepted");
}
void forces_and_rotation() {
    Config cfg;
    cfg.gravity = {};
    World w(cfg);
    BodyDesc b;
    b.mass = 2;
    b.linear_damping = 0;
    b.angular_damping = 0;
    const auto id = w.create_body(b);
    ColliderDesc c;
    c.body = id;
    c.shape = Shape::box({0.5_r, 0.5_r, 0.5_r});
    w.create_collider(c);
    w.step([&](World &fixed) { fixed.add_force(id, {60, 0, 0}); });
    check(near(w.body(id).properties.linear_velocity.x, 1_r, 0.00001_r),
          "force not integrated across all substeps");
    w.step();
    check(near(w.body(id).properties.linear_velocity.x, 1_r, 0.00001_r), "force not cleared");
    w.add_impulse_at(id, {0, 1, 0}, {1, 0, 0});
    check(w.body(id).properties.angular_velocity.z > 0, "off-center impulse did not rotate");
    const auto before = w.body(id).properties.pose.rotation.get_internal_vector();
    w.step();
    check(w.body(id).properties.pose.rotation.get_internal_vector() != before, "orientation not integrated");
    b.rotation_locks = 7;
    const auto locked = w.create_body(b);
    w.add_impulse(locked, {}, {2, 3, 4});
    w.step();
    check(w.body(locked).properties.angular_velocity == Vec3{}, "rotation locks ignored");
}
void query_tests() {
    const std::vector<ShapePtr> shapes{Shape::box({0.5_r, 0.5_r, 0.5_r}), Shape::sphere(0.5_r),
                                       Shape::capsule(0.5_r, 0.5_r), hull()};
    for (const auto &target : shapes) {
        World w;
        add(w, target, {0, 0, 0}, BodyType::static_body);
        for (const auto &source : shapes) {
            check(w.overlap(*source, {{0.2_r, 0, 0}, {}}).size() == 1, "missing convex overlap");
            check(w.overlap(*source, {{3, 0, 0}, {}}).empty(), "false convex overlap");
            const auto hit = w.cast(*source, {{-5, 0, 0}, {}}, {10, 0, 0});
            check(hit.complete && hit.hits.size() == 1 && near(hit.hits.front().fraction, 0.4_r, 0.0001_r),
                  "convex cast wrong fraction");
        }
        const auto hit = w.raycast({-5, 0, 0}, {2, 0, 0}, 10);
        check(hit.complete && hit.hits.size() == 1 && near(hit.hits[0].distance, 4.5_r, 0.0001_r),
              "ray wrong distance");
        check(hit.hits[0].normal.x < -0.99_r, "ray normal sign");
        const auto inside = w.raycast({0, 0, 0}, {0, 1, 0}, 2);
        check(inside.complete && inside.hits[0].started_overlapping && inside.hits[0].fraction == 0,
              "inside cast contract");
        const auto away = w.raycast({-5, 0, 0}, {-1, 0, 0}, 10);
        check(away.complete && away.hits.empty(), "backward ray false hit");
    }
    World corners;
    add(corners, Shape::box({1, 1, 1}), {}, BodyType::static_body);
    check(corners.overlap(*Shape::sphere(0.5_r), {{1.4_r, 1.4_r, 0}, {}}).empty(),
          "sphere-box corner false positive");
    World capsule;
    add(capsule, Shape::capsule(0.5_r, 1), {}, BodyType::static_body);
    auto cap = capsule.raycast({0, 5, 0}, {0, -1, 0}, 10);
    check(cap.complete && near(cap.hits[0].distance, 3.5_r, 0.0001_r), "capsule end cap query");
    World rotated;
    BodyDesc b;
    b.type = BodyType::static_body;
    b.pose.rotation = Quat(fxp::pi_v<Real> / 4_r, Vec3{0, 1, 0});
    ColliderDesc c;
    c.body = rotated.create_body(b);
    c.shape = Shape::box({1, 1, 0.25_r});
    rotated.create_collider(c);
    auto obb = rotated.boxcast({0.2_r, 0.2_r, 0.2_r}, {{-5, 0, 0}, {}}, {10, 0, 0});
    check(obb.complete && !obb.hits.empty(), "rotated boxcast failed");
    World ties;
    add(ties, Shape::sphere(1), {}, BodyType::static_body);
    add(ties, Shape::sphere(1), {}, BodyType::static_body);
    auto hits = ties.raycast({-5, 0, 0}, {1, 0, 0}, 10);
    check(hits.hits.size() == 2 && hits.hits[0].collider < hits.hits[1].collider, "equal-hit tie ordering");
    const auto digest = ties.checksum();
    for (int i = 0; i < 5; ++i) {
        ties.overlap(*Shape::sphere(2), {});
        ties.raycast({-5, 0, 0}, {1, 0, 0}, 10);
    }
    check(ties.checksum() == digest, "queries mutated simulation");
}
void layers_events_parenting() {
    Config cfg;
    cfg.gravity = {};
    World w(cfg);
    const auto a = add(w, Shape::sphere(1), {}, BodyType::kinematic, true, 1);
    const auto b = add(w, Shape::sphere(1), {1, 0, 0}, BodyType::dynamic, false, 2);
    w.set_layer_collision(1, 2, false);
    w.step();
    check(w.events().empty(), "disabled layers generated event");
    w.set_layer_collision(1, 2, true);
    w.step();
    check(w.events().size() == 1 && w.events()[0].type == EventType::begin && w.events()[0].sensor,
          "sensor begin missing");
    w.dispatch_events([&](const CollisionEvent &) {
        rejects([&] { w.teleport(b, {}); }, "event callback mutation allowed");
        rejects([&] { w.random_u64(); }, "event callback RNG mutation allowed");
    });
    w.step();
    check(w.events()[0].type == EventType::stay, "stay missing");
    w.remove_body(b);
    World restored;
    restored.load(w.save());
    restored.step();
    w.step();
    check(w.save() == restored.save(), "removed-owner pair did not restore");
    check(w.events()[0].type == EventType::end, "removal end missing");
    check(w.overlap(*Shape::sphere(3), {}).empty(), "default query included sensor");
    QueryFilter filter;
    filter.include_sensors = true;
    filter.layers = 1u << 1;
    check(w.overlap(*Shape::sphere(3), {}, filter).size() == 1, "sensor query filtering");
    filter.ignore_body = a;
    check(w.overlap(*Shape::sphere(3), {}, filter).empty(), "ignore-body filter");
    World compound(cfg);
    const auto body = compound.create_body();
    ColliderDesc parent;
    parent.body = body;
    parent.shape = Shape::sphere(0.25_r);
    parent.local_pose.position = {1, 0, 0};
    const auto p = compound.create_collider(parent);
    ColliderDesc child = parent;
    child.parent = p;
    const auto ch = compound.create_collider(child);
    auto overlap = compound.overlap(*Shape::sphere(0.1_r), {{2, 0, 0}, {}});
    check(overlap.size() == 1 && overlap[0].collider == ch, "child transform incorrect");
    compound.set_collider_parent(ch, 0, true);
    check(near(compound.collider(ch).properties.local_pose.position.x, 2_r),
          "reparent did not preserve pose");
    compound.set_collider_parent(ch, p, false);
    rejects([&] { compound.set_collider_parent(p, ch); }, "parent cycle accepted");
    compound.step();
    check(compound.events().empty(), "compound self collision");
    compound.remove_collider(p);
    rejects([&] { compound.collider(ch); }, "child survived parent removal");
}
void kinematic_and_interpolation() {
    Config cfg;
    cfg.gravity = {};
    World w(cfg);
    const auto k = add(w, Shape::box({0.5_r, 0.5_r, 0.5_r}), {-1, 0, 0}, BodyType::kinematic);
    const auto d = add(w, Shape::sphere(0.5_r), {0, 0, 0});
    for (int i = 0; i < 15; ++i)
        w.step([&](World &world) {
            auto p = world.body(k).properties.pose;
            p.position.x += 0.04_r;
            world.set_kinematic_target(k, p);
        });
    check(w.body(d).properties.pose.position.x > 0.3_r, "kinematic body did not push dynamic body");
    const auto current = w.body(k).properties.pose;
    const auto old = w.interpolated_pose(k, 0);
    const auto half = w.interpolated_pose(k, 0.5_r);
    check(near(half.position.x, (current.position.x + old.position.x) / 2_r, 0.000001_r),
          "interpolation wrong");
    const auto hash = w.checksum();
    w.render(0.25_r, [&](BodyId, const Pose &) {
        rejects([&] { w.set_velocity(k, {}); }, "render callback mutation allowed");
    });
    check(w.checksum() == hash, "render changed state");
    w.teleport(k, {{5, 0, 0}, {}}, true);
    check(w.interpolated_pose(k, 0).position.x == 5, "teleport left stale interpolation");
    rejects([&] { w.step([&](World &x) { x.step(); }); }, "recursive step allowed");
    w.step();
}
void stacks_and_friction() {
    World w;
    ground(w);
    std::vector<BodyId> boxes;
    for (int i = 0; i < 4; ++i)
        boxes.push_back(add(w, Shape::box({0.5_r, 0.5_r, 0.5_r}), {0, 0.5_r + Real(i) * 1.01_r, 0}));
    for (int i = 0; i < 220; ++i)
        w.step();
    for (std::size_t i = 0; i < boxes.size(); ++i) {
        const auto b = w.body(boxes[i]);
        check(near(b.properties.pose.position.y, 0.5_r + Real(static_cast<std::int64_t>(i)), 0.15_r),
              "stack collapsed");
        check(fxp::abs(b.properties.pose.position.x) < 0.2_r, "stack slid sideways");
    }
    World restored;
    restored.load(w.save());
    w.add_impulse(boxes.back(), {1, 0, 0});
    restored.add_impulse(boxes.back(), {1, 0, 0});
    check(!w.body(boxes.front()).sleeping, "connected island did not wake");
    for (int i = 0; i < 10; ++i) {
        w.step();
        restored.step();
        check(w.save() == restored.save(), "stack continuation mismatch");
    }
    const auto slide = [](Real friction) {
        World world;
        ground(world);
        BodyDesc b;
        b.pose.position = {0, 0.5_r, 0};
        b.linear_velocity = {5, 0, 0};
        b.rotation_locks = 7;
        b.linear_damping = 0;
        const auto id = world.create_body(b);
        ColliderDesc c;
        c.body = id;
        c.shape = Shape::box({0.5_r, 0.5_r, 0.5_r});
        c.material.static_friction = friction;
        c.material.dynamic_friction = friction;
        world.create_collider(c);
        for (int i = 0; i < 45; ++i)
            world.step();
        return world.body(id).properties.linear_velocity.x;
    };
    check(slide(1_r) < slide(0_r) - 1_r, "friction did not slow sliding");
}
void springs_and_surfaces() {
    Config cfg;
    cfg.gravity = {};
    World spring(cfg);
    const auto body = add(spring, Shape::sphere(0.25_r), {3, 0, 0});
    SpringDesc s;
    s.body_a = body;
    s.anchor_b = {0, 0, 0};
    s.rest_length = 1;
    s.frequency = 2;
    s.damping_ratio = 1;
    const auto joint = spring.create_spring(s);
    for (int i = 0; i < 100; ++i)
        spring.step();
    check(near(spring.body(body).properties.pose.position.x, 1_r, 0.05_r), "spring did not converge");
    spring.remove_spring(joint);
    for (int surface = 0; surface < 2; ++surface) {
        World w;
        if (surface == 0)
            add(w, plane_mesh(), {}, BodyType::static_body);
        else
            add(w, Shape::heightfield(21, 21, 1, std::vector<Real>(441, 0_r)), {-10, 0, -10},
                BodyType::static_body);
        const auto cap = add(w, Shape::capsule(0.3_r, 0.5_r), {-2, 2, 0});
        const auto convex = add(w, hull(), {2, 2, 0});
        const auto ray = w.raycast({0, 4, 0}, {0, -1, 0}, 8);
        check(ray.complete && !ray.hits.empty() && near(ray.hits.front().distance, 4_r, 0.001_r),
              "surface ray failed");
        for (int i = 0; i < 150; ++i)
            w.step();
        check(near(w.body(cap).properties.pose.position.y, 0.8_r, 0.08_r), "capsule did not rest on surface");
        check(near(w.body(convex).properties.pose.position.y, 0.5_r, 0.08_r),
              "convex did not rest on surface");
        World copy;
        copy.load(w.save());
        check(copy.save() == w.save(), "surface save mismatch");
    }
}
World scenario() {
    Config cfg;
    cfg.substeps = 2;
    World w(cfg, 0x12345678);
    ground(w);
    const auto box = add(w, Shape::box({0.5_r, 0.5_r, 0.5_r}), {0, 2, 0});
    add(w, Shape::sphere(0.4_r), {1, 3, 0});
    add(w, Shape::capsule(0.3_r, 0.4_r), {-2, 2, 0});
    add(w, hull(), {3, 2, 0});
    add(w, Shape::box({0.5_r, 0.5_r, 0.5_r}), {-4, 0.5_r, 0}, BodyType::kinematic);
    add(w, Shape::sphere(1), {0, 1, 0}, BodyType::static_body, true, 3);
    ColliderDesc child;
    child.body = box;
    child.shape = Shape::sphere(0.2_r);
    child.local_pose.position = {0.6_r, 0, 0};
    w.create_collider(child);
    const auto spring_body = add(w, Shape::sphere(0.25_r), {5, 2, 0});
    SpringDesc spring;
    spring.body_a = spring_body;
    spring.anchor_b = {5, 3, 0};
    spring.rest_length = 1;
    spring.frequency = 3;
    w.create_spring(spring);
    add(w, Shape::heightfield(3, 3, 1, {0_r, 0.2_r, 0_r, 0.1_r, 0.3_r, 0.2_r, 0_r, 0.1_r, 0_r}), {8, 0, 0},
        BodyType::static_body);
    add(w, plane_mesh(), {0, -3, 0}, BodyType::static_body);
    return w;
}
void command(World &w) {
    const auto tick = w.tick();
    auto pose = w.body(6).properties.pose;
    pose.position.x += 0.02_r;
    w.set_kinematic_target(6, pose);
    if (tick == 15)
        w.add_impulse(2, {0.3_r, 1, 0.2_r});
    if (tick == 31)
        w.set_layer_collision(0, 3, false);
    if (tick == 40)
        w.set_layer_collision(0, 3, true);
    if (tick == 55) {
        w.remove_body(5);
        add(w, hull(), {3, 3, 0});
    }
    if (tick == 70)
        w.teleport(3, {{1, 2, 0}, {}}, true);
    if (tick % 13 == 0)
        w.add_force(2, {w.random_unit() - 0.5_r, 0, 0});
}
std::uint64_t query_digest(World &w) {
    std::uint64_t hash = UINT64_C(14695981039346656037);
    const auto mix = [&](std::uint64_t x) {
        for (unsigned i = 0; i < 8; ++i) {
            hash ^= (x >> (i * 8)) & 255u;
            hash *= UINT64_C(1099511628211);
        }
    };
    const auto cast = w.spherecast(0.1_r, {{-6, 0.5_r, 0}, {}}, {18, 0, 0}, {0xffffffffu, 0, true});
    mix(cast.complete ? 1 : 0);
    for (const auto &h : cast.hits) {
        mix(h.collider);
        mix(h.triangle);
        mix(static_cast<std::uint64_t>(h.fraction.raw_value()));
        mix(static_cast<std::uint64_t>(h.normal.x.raw_value()));
        mix(static_cast<std::uint64_t>(h.normal.y.raw_value()));
        mix(static_cast<std::uint64_t>(h.normal.z.raw_value()));
    }
    for (const auto &h : w.overlap(*Shape::box({2, 2, 2}), {{0, 1, 0}, {}}, {0xffffffffu, 0, true}))
        mix(h.collider);
    return hash;
}
void save_continuation() {
    World w = scenario();
    for (int i = 0; i < 37; ++i)
        w.step(command);
    w.add_force(3, {1, 2, 3});
    auto target = w.body(6).properties.pose;
    target.position.x += 0.02_r;
    w.set_kinematic_target(6, target);
    const auto bytes = w.save();
    World restored;
    restored.load(bytes);
    check(restored.save() == bytes, "snapshot byte roundtrip");
    for (int i = 0; i < 55; ++i) {
        w.step(command);
        restored.step(command);
        restored.render(0.3_r, [](BodyId, const Pose &) {});
        query_digest(restored);
        check(w.save() == restored.save(), "save/load diverged at tick " + std::to_string(w.tick()));
    }
    Bytes corrupt = bytes;
    corrupt[30] ^= 1;
    const auto before = restored.save();
    rejects([&] { restored.load(corrupt); }, "corruption accepted");
    check(before == restored.save(), "failed load mutated world");
    for (std::size_t cut : {std::size_t(0), std::size_t(8), bytes.size() / 2, bytes.size() - 1}) {
        Bytes truncated(bytes.begin(), bytes.begin() + static_cast<std::ptrdiff_t>(cut));
        rejects([&] { restored.load(truncated); }, "truncated save accepted");
    }
    // Recompute integrity hash to exercise semantic validation, not just corruption detection.
    corrupt = bytes;
    corrupt[8] = 99;
    std::uint64_t hash = UINT64_C(14695981039346656037);
    for (std::size_t i = 0; i < corrupt.size() - 8; ++i) {
        hash ^= corrupt[i];
        hash *= UINT64_C(1099511628211);
    }
    for (unsigned i = 0; i < 8; ++i)
        corrupt[corrupt.size() - 8 + i] = static_cast<std::uint8_t>(hash >> (i * 8));
    rejects([&] { restored.load(corrupt); }, "incompatible format accepted");
    check(before == restored.save(), "version rejection mutated world");
}
std::vector<std::string> reference_lines() {
    World w = scenario();
    std::vector<std::string> lines;
    for (int i = 0; i < 100; ++i) {
        w.step(command);
        std::ostringstream line;
        line << w.tick() << ' ' << w.checksum() << ' ' << query_digest(w);
        lines.push_back(line.str());
    }
    return lines;
}
void reference(const char *write_path = nullptr) {
    const auto actual = reference_lines();
    if (write_path) {
        std::ofstream out(write_path);
        check(bool(out), "cannot write reference");
        for (const auto &line : actual)
            out << line << '\n';
        return;
    }
    std::ifstream in(std::string(FXP_PHYSICS_TEST_DIR) + "/reference.txt");
    check(bool(in), "missing frozen reference.txt");
    for (const auto &line : actual) {
        std::string expected;
        check(bool(std::getline(in, expected)), "truncated reference");
        check(line == expected, "determinism mismatch: expected " + expected + ", actual " + line);
    }
    std::string extra;
    check(!std::getline(in, extra), "extra reference rows");
}
void edge_cases() {
    Config config;
    config.gravity = {};
    config.substeps = 4;
    World crossing(config);
    add(crossing, Shape::sphere(0.1_r), {}, BodyType::static_body, true);
    const auto moving = add(crossing, Shape::sphere(0.1_r), {-1, 0, 0}, BodyType::kinematic);
    crossing.set_kinematic_target(moving, {{1, 0, 0}, {}});
    crossing.step();
    check(crossing.events().size() == 1 && crossing.events()[0].type == EventType::begin,
          "transient substep sensor contact lost");
    crossing.step();
    check(crossing.events().size() == 1 && crossing.events()[0].type == EventType::end,
          "transient contact end missing");

    Config impact_config;
    impact_config.gravity = {};
    impact_config.restitution_threshold = 0;
    World impact(impact_config);
    BodyDesc body;
    body.linear_damping = 0;
    body.angular_damping = 0;
    body.pose.position = {-0.5_r, 0, 0};
    body.linear_velocity = {2, 0, 0};
    const auto a = impact.create_body(body);
    body.pose.position = {0.5_r, 0, 0};
    body.linear_velocity = {};
    const auto b = impact.create_body(body);
    ColliderDesc collider;
    collider.shape = Shape::sphere(0.5_r);
    collider.material.restitution = 1;
    collider.body = a;
    impact.create_collider(collider);
    collider.body = b;
    impact.create_collider(collider);
    impact.step();
    check(near(impact.body(a).properties.linear_velocity.x, 0_r, 0.001_r) &&
              near(impact.body(b).properties.linear_velocity.x, 2_r, 0.001_r),
          "elastic equal-mass sphere impact");

    Random rng(9876);
    World scatter(config);
    std::vector<Vec3> centers;
    const auto shape = Shape::sphere(0.3_r);
    for (int i = 0; i < 100; ++i) {
        Vec3 center{(Real(rng.bounded(81)) - 40_r) / 4_r, (Real(rng.bounded(81)) - 40_r) / 4_r,
                    (Real(rng.bounded(81)) - 40_r) / 4_r};
        centers.push_back(center);
        add(scatter, shape, center, BodyType::static_body);
    }
    for (int q = 0; q < 50; ++q) {
        const Vec3 center = centers[static_cast<std::size_t>(q)] + Vec3{0.25_r, 0, 0};
        std::vector<BodyId> expected;
        for (std::size_t i = 0; i < centers.size(); ++i)
            if (fxp::dot(centers[i] - center, centers[i] - center) <= 0.7_r * 0.7_r)
                expected.push_back(static_cast<BodyId>(i + 1));
        const auto hits = scatter.overlap(*Shape::sphere(0.4_r), {center, {}});
        std::vector<BodyId> actual;
        for (const auto &hit : hits)
            actual.push_back(hit.body);
        check(actual == expected, "broadphase/overlap differs from brute-force sphere oracle");
    }
    for (int i = 0; i < 8; ++i) {
        World rotated(config);
        BodyDesc box;
        box.type = BodyType::static_body;
        box.pose.position = {-100, 7, -80};
        box.pose.rotation = Quat(Real(i + 1) / 7_r, Vec3{1, 2, 3}, true);
        const auto id = rotated.create_body(box);
        ColliderDesc c;
        c.body = id;
        c.shape = Shape::box({2, 1, 3});
        rotated.create_collider(c);
        const auto p = rotated.body(id).properties.pose;
        const Vec3 corner = p.position + p.rotation.rotate({2, 1, 3});
        check(!rotated.overlap(*Shape::sphere(0.005_r), {corner, {}}).empty(),
              "rotated broadphase excluded a corner");
    }
    World shared(config), separate_assets(config);
    for (int i = 0; i < 4; ++i) {
        add(shared, shape, {Real(i), 0, 0}, BodyType::kinematic);
        add(separate_assets, Shape::sphere(0.3_r), {Real(i), 0, 0}, BodyType::kinematic);
    }
    check(shared.save() == separate_assets.save(), "asset allocation/sharing changed bytes");
}
void length_comparisons() {
    const auto verify = [](Vec3 v, Real threshold) {
        check(detail::length_greater_than(v, threshold) == (fxp::length(v) > threshold),
              "raw squared length comparison differs from rounded/saturated length");
    };
    const Real raw[] = {Real::min(), Real::max(), Real::from_raw(-1), 0_r,
                        Real::from_raw(1), 0.08_r, 1_r};
    for (auto x : raw)
        for (auto y : raw)
            for (auto z : raw)
                for (auto threshold : raw)
                    verify({x, y, z}, threshold);
    Random rng(81923);
    for (unsigned i = 0; i < 4000; ++i) {
        Vec3 v;
        for (unsigned axis = 0; axis < 3; ++axis) {
            auto value = static_cast<std::int64_t>(rng.next() >> 1);
            // Include tiny magnitudes where rounding individual fixed squares
            // would discard the length, and near-overflow magnitudes.
            value >>= (i % 63);
            v[axis] = Real::from_raw((rng.next() & 1) ? value : -value);
        }
        const auto length = fxp::length(v);
        verify(v, length);
        verify(v, Real::from_raw(length.raw_value() - 1));
        if (length != Real::max())
            verify(v, Real::from_raw(length.raw_value() + 1));
    }
}

void rotation_cache() {
    detail::Rotation cache;
    Random rng(17823);
    for (unsigned i = 0; i < 1200; ++i) {
        Quat q;
        if (i % 11) {
            fxp::vec4<Real> components;
            for (unsigned axis = 0; axis < 4; ++axis)
                components[axis] = Real(static_cast<int>(rng.bounded(2001)) - 1000) / 1000_r;
            q = Quat(components);
            check(q.normalize(), "generated zero rotation");
        }
        Vec3 v;
        for (unsigned axis = 0; axis < 3; ++axis) {
            auto raw = static_cast<std::int64_t>(rng.next() >> 1);
            raw >>= i % 63;
            v[axis] = Real::from_raw((rng.next() & 1) ? raw : -raw);
        }
        if (i % 13 == 0)
            v.x = Real::min();
        cache.set(q);
        check(cache.rotate(v) == q.rotate(v), "rotation cache changed fixed-point rounding");
        check(cache.rotate(v, true) == detail::inverse(q).rotate(v), "inverse rotation cache changed rounding");
        cache.set(q);
        check(cache.rotate(v) == q.rotate(v), "reused rotation cache changed result");
    }
}

void moving_broadphase() {
    Config config;
    config.gravity = {};
    World queried(config), unqueried(config);
    const auto sphere = Shape::sphere(0.3_r), probe = Shape::sphere(0.4_r);
    std::vector<BodyId> ids;
    for (int i = 0; i < 80; ++i) {
        const Vec3 p{Real(i % 10) * 2, 0, Real(i / 10) * 2};
        ids.push_back(add(queried, sphere, p, BodyType::kinematic));
        check(add(unqueried, sphere, p, BodyType::kinematic) == ids.back(), "body ID mismatch");
    }
    for (unsigned tick = 0; tick < 80; ++tick) {
        if (tick == 20) {
            // Remove/reinsert without changing population to exercise sparse
            // IDs and cache invalidation independently of vector size.
            queried.remove_body(ids[7]);
            unqueried.remove_body(ids[7]);
            ids.erase(ids.begin() + 7);
            ids.push_back(add(queried, sphere, {}, BodyType::kinematic));
            check(add(unqueried, sphere, {}, BodyType::kinematic) == ids.back(), "replacement ID mismatch");
        }
        for (std::size_t i = 0; i < ids.size(); ++i) {
            const Vec3 p{Real(static_cast<int>((i * 17 + tick * 3) % 41) - 20) / 2_r,
                         0, Real(static_cast<int>((i * 7 + tick) % 31) - 15) / 2_r};
            queried.teleport(ids[i], {p, {}});
            unqueried.teleport(ids[i], {p, {}});
        }
        // Force different rebuild/refit histories without different commands.
        for (std::size_t q = 0; q < 5; ++q) {
            const Vec3 center = queried.body(ids[q]).properties.pose.position + Vec3{0.25_r, 0, 0};
            std::vector<BodyId> expected, actual;
            for (const auto id : ids)
                if (fxp::length(queried.body(id).properties.pose.position - center) <= 0.7_r)
                    expected.push_back(id);
            for (const auto &hit : queried.overlap(*probe, {center, {}}))
                actual.push_back(hit.body);
            check(actual == expected, "moving broadphase differs from brute-force oracle");
        }
        queried.step();
        unqueried.step();
        check(queried.save() == unqueried.save(), "query/refit frequency changed simulation state");
        if (tick == 40) {
            World restored;
            restored.load(unqueried.save());
            unqueried = std::move(restored);
        }
    }
}

void benchmark() {
    for (int count : {100, 1000, 10000}) {
        Config cfg;
        cfg.gravity = {};
        World w(cfg);
        const auto shape = Shape::sphere(0.4_r);
        for (int i = 0; i < count; ++i) {
            const auto b = add(w, shape, {Real(i % 100) * 2_r, 0, Real(i / 100) * 2_r}, BodyType::kinematic);
            w.set_velocity(b, {0.5_r, 0, 0});
        }
        w.step();
        const auto begin = std::chrono::steady_clock::now();
        for (int i = 0; i < 10; ++i) {
            w.step();
            for (int q = 0; q < 100; ++q)
                w.overlap(*shape, {{Real(q) * 2_r, 0, 0}, {}});
        }
        const auto us =
            std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - begin)
                .count() /
            10;
        std::cout << count << " moving kinematic spheres, 100 overlaps/tick: " << us << " us/tick, "
                  << w.stats().broadphase_pairs << " candidate pairs\n";
    }
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc > 1 && std::string(argv[1]) == "--reference") {
            reference();
            std::cout << "100 frozen tick/state/query rows matched\n";
            return 0;
        }
        if (argc > 2 && std::string(argv[1]) == "--write-reference") {
            reference(argv[2]);
            return 0;
        }
        if (argc > 1 && std::string(argv[1]) == "--benchmark") {
            benchmark();
            return 0;
        }
        const std::pair<const char *, void (*)()> tests[] = {
            {"rng and validation", rng_and_validation},
            {"forces and rotation", forces_and_rotation},
            {"queries", query_tests},
            {"layers/events/parenting", layers_events_parenting},
            {"kinematic/interpolation", kinematic_and_interpolation},
            {"rest and sleep", smoke},
            {"stacks and friction", stacks_and_friction},
            {"springs and surfaces", springs_and_surfaces},
            {"save continuation", save_continuation},
            {"edge cases and broadphase oracle", edge_cases},
            {"exact length comparisons", length_comparisons},
            {"rotation coefficient cache", rotation_cache},
            {"moving broadphase and query independence", moving_broadphase}};
        for (const auto &test : tests) {
            test.second();
            std::cout << "PASS: " << test.first << std::endl;
        }
        std::cout << "physics tests passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
