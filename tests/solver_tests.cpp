#include <fxp_physics/physics.hpp>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace fxp_physics;

namespace {
void check(bool condition, const char *message) {
    if (!condition)
        throw std::runtime_error(message);
}
BodyId add(World &world, ShapePtr shape, BodyDesc body) {
    ColliderDesc collider;
    collider.body = world.create_body(body);
    collider.shape = std::move(shape);
    world.create_collider(collider);
    return collider.body;
}
ShapePtr hull(Vec3 half) {
    return Shape::convex(
        {{-half.x, -half.y, -half.z}, {half.x, -half.y, -half.z},
         {-half.x, half.y, -half.z}, {half.x, half.y, -half.z},
         {-half.x, -half.y, half.z}, {half.x, -half.y, half.z},
         {-half.x, half.y, half.z}, {half.x, half.y, half.z}},
        {{0, 4, 6, 2}, {1, 3, 7, 5}, {0, 1, 5, 4},
         {2, 6, 7, 3}, {0, 2, 3, 1}, {4, 5, 7, 6}});
}
World islands() {
    Config config;
    config.timestep = 1_r / 20_r;
    config.substeps = 4;
    config.velocity_iterations = 12;
    config.sleep_ticks = 15;
    World world(config, 12345);
    BodyDesc ground;
    ground.type = BodyType::kinematic;
    ground.pose.position = {0, -0.5_r, 0};
    add(world, Shape::box({80, 0.5_r, 10}), ground);
    const auto box = Shape::box({0.5_r, 0.5_r, 0.5_r});
    // Interleaved IDs make contact rows for different islands interleave too.
    for (unsigned layer = 0; layer < 3; ++layer) {
        for (unsigned tower = 0; tower < 44; ++tower) {
            BodyDesc body;
            body.pose.position = {Real(tower) * 3 - 65, Real(layer) + 0.5_r, 0};
            if (tower % 5 == 0)
                body.rotation_locks = 7;
            if (tower % 7 == 0) {
                body.mass = 2;
                body.center_of_mass = {0.05_r, 0, 0.02_r};
                body.inertia = {0.4_r, 0.3_r, 0.5_r, 0.01_r, -0.02_r, 0.01_r};
            }
            body.allow_sleep = tower != 19;
            const auto id = add(world, box, body);
            if (layer == 2 && tower % 4 == 1) {
                SpringDesc spring;
                spring.body_a = id;
                spring.anchor_a = {0.25_r, 0.2_r, 0};
                spring.anchor_b = body.pose.position + Vec3{0, 2, 0};
                spring.rest_length = 2;
                spring.frequency = 1;
                spring.max_force = 2;
                world.create_spring(spring);
            }
        }
    }
    // Springs between dynamic bodies join otherwise independent islands.
    SpringDesc spring;
    spring.body_a = 42;
    spring.body_b = 43;
    spring.rest_length = 3;
    spring.frequency = 0.5_r;
    spring.max_force = 1;
    spring.collide_connected = true;
    world.create_spring(spring);
    return world;
}
void commands(World &world) {
    const auto tick = world.tick();
    if (tick == 8)
        world.set_velocity(1, {0.1_r, 0, 0});
    if (tick == 15)
        world.set_velocity(1, {});
    if (tick == 24)
        world.add_impulse_at(52, {0.3_r, 0.2_r, 0.1_r}, {1, 3, 0});
    if (tick == 32)
        world.remove_body(48);
    if (tick == 48)
        world.wake(12);
}
void worker_determinism() {
    auto serial = islands();
    World two, four, automatic;
    const auto initial = serial.save();
    two.load(initial);
    four.load(initial);
    automatic.load(initial);
    serial.set_worker_count(1);
    two.set_worker_count(2);
    four.set_worker_count(4);
    automatic.set_worker_count(0);
    check(serial.save() == two.save() && two.save() == four.save() &&
              four.save() == automatic.save(),
          "worker configuration leaked into serialized state");
    bool enough_contacts = false, sleeping_contacts = false;
    for (unsigned tick = 0; tick < 80; ++tick) {
        serial.step(commands);
        two.step(commands);
        four.step(commands);
        automatic.step(commands);
        const auto state = serial.save();
        check(state == two.save() && state == four.save() && state == automatic.save(),
              "worker count changed the simulation");
        enough_contacts |= serial.stats().contacts >= 128;
        for (const auto id : serial.body_ids())
            sleeping_contacts |= serial.body(id).sleeping;
        if (tick == 20) {
            four.load(state);
            check(four.worker_count() == 4, "load changed worker configuration");
        }
        if (tick == 40) {
            two.set_worker_count(4);
            four.set_worker_count(2);
            automatic.set_worker_count(1);
        }
        if (tick == 60) {
            two.set_worker_count(2);
            four.set_worker_count(4);
            automatic.set_worker_count(0);
        }
    }
    check(enough_contacts, "parallel test did not reach the parallel workload threshold");
    check(sleeping_contacts, "worker test did not exercise omitted sleeping contact rows");
    // Destruction, moves and reloads must also stop/join existing workers.
    World moved(std::move(four));
    check(moved.save() == serial.save(), "moving a world changed its state");
    moved = std::move(two);
    moved.step(commands);
    serial.step(commands);
    check(moved.save() == serial.save(), "move assignment broke worker state");
    bool rejected = false;
    try {
        moved.set_worker_count(65);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    check(rejected, "unbounded worker count accepted");
}
World mixed_contacts(bool grounds_first) {
    Config config;
    config.timestep = 1_r / 20_r;
    config.substeps = 2;
    config.velocity_iterations = 4;
    World world(config);
    const auto grounds = [&] {
        BodyDesc ground;
        ground.type = BodyType::static_body;
        ground.pose.position = {0, -0.5_r, -8};
        add(world, hull({50, 0.5_r, 3}), ground);
        ground.pose.position = {};
        add(world, Shape::mesh({{-50, 0, -3}, {50, 0, -3}, {50, 0, 3}, {-50, 0, 3}},
                              {{0, 2, 1}, {0, 3, 2}}), ground);
        ground.pose.position = {-50, 0, 5};
        add(world, Shape::heightfield(51, 4, 2, std::vector<Real>(204, 0_r)), ground);
    };
    if (grounds_first)
        grounds();
    const ShapePtr shapes[] = {
        hull({0.45_r, 0.45_r, 0.45_r}), Shape::box({0.45_r, 0.45_r, 0.45_r}),
        Shape::capsule(0.25_r, 0.2_r), Shape::sphere(0.45_r)};
    for (unsigned surface = 0; surface < 3; ++surface) {
        for (unsigned layer = 0; layer < 2; ++layer) {
            for (unsigned column = 0; column < 48; ++column) {
                BodyDesc body;
                body.pose.position = {Real(column) * 2 - 47, 0.45_r + Real(layer) * 0.9_r,
                                      Real(surface) * 8 - 8};
                body.allow_sleep = false;
                body.angular_velocity = {0.02_r, 0.01_r, 0.03_r};
                ColliderDesc collider;
                collider.body = world.create_body(body);
                collider.shape = shapes[(column + layer) % 4];
                collider.sensor = layer == 1 && column % 8 == 0;
                world.create_collider(collider);
            }
        }
    }
    if (!grounds_first)
        grounds();
    return world;
}
void mixed_contact_workers() {
    for (const bool grounds_first : {false, true}) {
        auto serial = mixed_contacts(grounds_first);
        World parallel;
        parallel.load(serial.save());
        serial.set_worker_count(1);
        parallel.set_worker_count(4);
        bool enough_pairs = false;
        for (unsigned tick = 0; tick < 20; ++tick) {
            // Different query/cache histories must not change lazy geometry or
            // the canonical ordering of shared-ground and triangle contacts.
            (void)parallel.raycast({-47, 4, -8}, {0, -1, 0}, 8);
            if (tick == 8)
                parallel.load(serial.save());
            if (tick == 12)
                parallel.set_worker_count(2);
            serial.step();
            parallel.step();
            check(serial.save() == parallel.save(),
                  "parallel generic/rounded/concave contact generation changed state or events");
            enough_pairs |= serial.stats().broadphase_pairs >= 128;
        }
        check(enough_pairs, "mixed-shape test did not exercise parallel narrowphase");
    }
}
void spring_extremes() {
    Config config;
    config.timestep = 1_r / 20_r;
    config.gravity = {};
    World world(config);
    BodyDesc body;
    body.mass = 10000;
    body.pose.position = {3, 0, 0};
    body.rotation_locks = 7;
    const auto id = add(world, Shape::sphere(0.1_r), body);
    SpringDesc spring;
    spring.body_a = id;
    spring.frequency = Real::epsilon();
    spring.damping_ratio = 0;
    const auto tiny = world.create_spring(spring);
    world.step();
    check(world.body(id).properties.linear_velocity == Vec3{},
          "unrepresentable spring stiffness generated motion");
    world.remove_spring(tiny);
    spring.frequency = 3;
    spring.damping_ratio = 1;
    world.create_spring(spring);
    for (unsigned i = 0; i < 100; ++i)
        world.step();
    check(fxp::abs(world.body(id).properties.pose.position.x - 1_r) < 0.1_r,
          "heavy spring failed to converge with reciprocal effective mass");

    // Small inertia makes the effective mass so large that its reciprocal
    // needs the accurate division fallback rather than a rounded low-bit value.
    World tiny_body(config);
    body.mass = 0.01_r;
    body.rotation_locks = 0;
    body.pose.position = {2, 0, 0};
    body.inertia = {0.000001_r, 0.000001_r, 0.000001_r, 0, 0, 0};
    const auto small = add(tiny_body, Shape::sphere(0.1_r), body);
    spring.body_a = small;
    spring.anchor_a = {0, 0.1_r, 0};
    spring.rest_length = 2;
    tiny_body.create_spring(spring);
    for (unsigned i = 0; i < 10; ++i)
        tiny_body.step();
    const auto state = tiny_body.body(small).properties;
    check(fxp::length(state.linear_velocity) <= config.max_linear_speed &&
              fxp::length(state.angular_velocity) <= config.max_angular_speed + 0.001_r,
          "extreme effective mass produced an unbounded response");
}
void twenty_hz_stack() {
    Config config;
    config.timestep = 1_r / 20_r;
    config.substeps = 4;
    config.velocity_iterations = 12;
    World world(config);
    BodyDesc ground;
    ground.type = BodyType::static_body;
    ground.pose.position = {0, -0.5_r, 0};
    add(world, Shape::box({10, 0.5_r, 10}), ground);
    const auto box = Shape::box({0.5_r, 0.5_r, 0.5_r});
    for (unsigned i = 0; i < 8; ++i) {
        BodyDesc body;
        body.pose.position.y = Real(i) + 0.5_r;
        add(world, box, body);
    }
    for (unsigned i = 0; i < 180; ++i)
        world.step();
    for (unsigned i = 0; i < 8; ++i) {
        const auto body = world.body(i + 2);
        check(fxp::abs(body.properties.pose.position.y - (Real(i) + 0.5_r)) < 0.15_r,
              "20 Hz stack lost support");
        check(fxp::length(body.properties.linear_velocity) < 0.1_r,
              "20 Hz stack failed to settle");
    }
}
void local_rotation_lock() {
    Config config;
    config.timestep = 1_r / 20_r;
    config.substeps = 4;
    config.gravity = {};
    World world(config);
    BodyDesc body;
    body.pose.position = {1, 2, 3};
    body.pose.rotation = Quat(0.6_r, {0, 0, 1});
    body.rotation_locks = 1;
    body.inertia = {0.4_r, 0.3_r, 0.5_r, 0.03_r, 0.02_r, 0.01_r};
    body.center_of_mass = {0.1_r, 0.02_r, 0};
    const auto id = add(world, Shape::box({0.5_r, 0.6_r, 0.7_r}), body);
    SpringDesc spring;
    spring.body_a = id;
    spring.anchor_a = {0.4_r, 0.3_r, 0.5_r};
    spring.anchor_b = {4, 3, 1};
    spring.frequency = 1;
    world.create_spring(spring);
    bool rotated = false;
    for (unsigned i = 0; i < 60; ++i) {
        world.step();
        const auto state = world.body(id).properties;
        auto inverse = state.pose.rotation;
        inverse.unit_inverse();
        const auto local = inverse.rotate(state.angular_velocity);
        check(fxp::abs(local.x) < 0.00001_r,
              "cached angular response violated a rotated local rotation lock");
        rotated |= fxp::abs(local.y) + fxp::abs(local.z) > 0.01_r;
    }
    check(rotated, "rotation-lock test never exercised unlocked angular axes");
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--narrowphase") {
            mixed_contact_workers();
            std::cout << "narrowphase worker checks passed\n";
            return 0;
        }
        check(argc == 1, "usage: fxp_physics_solver_tests [--narrowphase]");
        worker_determinism();
        mixed_contact_workers();
        spring_extremes();
        twenty_hz_stack();
        local_rotation_lock();
        std::cout << "solver checks passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
