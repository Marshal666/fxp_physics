#include <fxp_physics/physics.hpp>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace fxp_physics;

namespace {
void check(bool condition, const std::string &message) {
    if (!condition)
        throw std::runtime_error(message);
}
struct Object { BodyId body; ColliderId collider; };
Object box(World &world, Vec3 position, Vec3 half = {0.5_r, 0.5_r, 0.5_r},
           BodyType type = BodyType::dynamic, bool sensor = false) {
    BodyDesc body;
    body.type = type;
    body.pose.position = position;
    ColliderDesc collider;
    collider.body = world.create_body(body);
    collider.shape = Shape::box(half);
    collider.sensor = sensor;
    return {collider.body, world.create_collider(collider)};
}
void compare_step(World &cached, World &fresh, unsigned count) {
    for (unsigned i = 0; i < count; ++i) {
        // Loading every tick forces rebuilding all transient caches, including
        // when only some islands sleep or a contact's normal changes.
        fresh.load(fresh.save());
        cached.step();
        fresh.step();
        check(cached.save() == fresh.save(),
              "contact cache changed continuation at tick " + std::to_string(cached.tick()));
    }
}
void sleeping_and_mutation() {
    Config config;
    config.timestep = 1_r / 20_r;
    config.substeps = 4;
    config.velocity_iterations = 12;
    World cached(config), fresh(config);
    cached.set_worker_count(1);
    fresh.set_worker_count(4);
    const auto floor = box(cached, {0, -0.5_r, 0}, {20, 0.5_r, 20}, BodyType::static_body);
    const auto bottom = box(cached, {0, 0.5_r, 0});
    const auto top = box(cached, {0, 1.51_r, 0});
    const auto separate = box(cached, {4, 0.5_r, 0});
    const auto trigger = box(cached, {0, 0.5_r, 0}, {0.6_r, 0.6_r, 0.6_r},
                             BodyType::static_body, true);
    fresh.load(cached.save());
    compare_step(cached, fresh, 180);
    check(cached.body(bottom.body).sleeping && cached.body(top.body).sleeping &&
              cached.body(separate.body).sleeping, "test stack did not sleep");
    compare_step(cached, fresh, 8);
    for (const auto &event : cached.events())
        check(event.type == EventType::stay, "stationary cached pair lost stay event");
    check(!cached.events().empty(), "sleeping contacts lost collision events");

    // Active visitors still discover a sleeping island through the broadphase.
    const auto visitor = box(cached, {-2, 0.5_r, 0});
    const auto visitor2 = box(fresh, {-2, 0.5_r, 0});
    check(visitor.body == visitor2.body, "creation IDs differ");
    cached.set_velocity(visitor.body, {8, 0, 0});
    fresh.set_velocity(visitor.body, {8, 0, 0});
    compare_step(cached, fresh, 8);
    check(!cached.body(bottom.body).sleeping, "incoming body failed to wake sleeping stack");
    check(cached.body(separate.body).sleeping, "unconnected stack woke unnecessarily");

    cached.set_collider_filter(trigger.collider, 3, true);
    fresh.set_collider_filter(trigger.collider, 3, true);
    cached.set_layer_collision(0, 3, false);
    fresh.set_layer_collision(0, 3, false);
    compare_step(cached, fresh, 2);
    cached.set_layer_collision(0, 3, true);
    fresh.set_layer_collision(0, 3, true);
    cached.set_material(floor.collider, {0.8_r, 0.5_r, 0});
    fresh.set_material(floor.collider, {0.8_r, 0.5_r, 0});
    compare_step(cached, fresh, 2);
    cached.remove_body(visitor.body);
    fresh.remove_body(visitor.body);
    cached.set_collider_pose(trigger.collider, {{5, 0, 0}, {}});
    fresh.set_collider_pose(trigger.collider, {{5, 0, 0}, {}});
    compare_step(cached, fresh, 3);
    cached.teleport(separate.body, {{6, 2, 0}, {}}, true);
    fresh.teleport(separate.body, {{6, 2, 0}, {}}, true);
    compare_step(cached, fresh, 3);
}
void stationary_kinematics_and_springs() {
    Config config;
    config.gravity = {};
    World cached(config), fresh(config);
    const auto a = box(cached, {}, {0.5_r, 0.5_r, 0.5_r}, BodyType::kinematic, true);
    const auto b = box(cached, {0.5_r, 0, 0}, {0.5_r, 0.5_r, 0.5_r}, BodyType::kinematic);
    fresh.load(cached.save());
    compare_step(cached, fresh, 3);
    check(cached.events().size() == 1, "stationary sensor pair missing");
    SpringDesc spring;
    spring.body_a = a.body;
    spring.body_b = b.body;
    const auto s = cached.create_spring(spring);
    check(fresh.create_spring(spring) == s, "spring IDs differ");
    compare_step(cached, fresh, 2);
    check(cached.events().empty(), "excluded stationary pair was retained");
    cached.remove_spring(s);
    fresh.remove_spring(s);
    compare_step(cached, fresh, 2);
    check(cached.events().size() == 1, "removed exclusion did not restore stationary pair");
    cached.set_kinematic_target(b.body, {{2, 0, 0}, {}});
    fresh.set_kinematic_target(b.body, {{2, 0, 0}, {}});
    compare_step(cached, fresh, 1);
    check(cached.events().size() == 1 && cached.events().front().type == EventType::end,
          "kinematic target did not invalidate stationary contact");
    compare_step(cached, fresh, 1);
    check(cached.events().empty(), "ended sensor contact persisted");
}
void removing_final_colliders() {
    for (const bool remove_bodies : {false, true}) {
        Config config;
        config.gravity = {};
        config.sleep_ticks = 1;
        World cached(config), fresh(config);
        cached.set_worker_count(1);
        fresh.set_worker_count(4);
        const auto floor = box(cached, {0, -0.5_r, 0}, {2, 0.5_r, 2}, BodyType::static_body);
        const auto resting = box(cached, {0, 0.5_r, 0});
        cached.step();
        cached.step();
        check(cached.stats().contacts > 0 && cached.stats().broadphase_pairs > 0 &&
                  cached.events().size() == 1 && cached.events().front().type == EventType::stay,
              "empty-world regression did not begin with cached contacts");
        if (remove_bodies) {
            cached.remove_body(resting.body);
            cached.remove_body(floor.body);
        } else {
            cached.remove_collider(resting.collider);
            cached.remove_collider(floor.collider);
        }
        fresh.load(cached.save());
        compare_step(cached, fresh, 1);
        check(cached.stats().colliders == 0 && cached.stats().contacts == 0 &&
                  cached.stats().broadphase_pairs == 0,
              "removing every collider retained stale contact statistics");
        check(cached.stats().bodies == (remove_bodies ? 0u : 2u),
              "removing colliders changed the remaining body count");
        check(cached.events().size() == 1 && cached.events().front().type == EventType::end &&
                  cached.events().front().collider_a == floor.collider &&
                  cached.events().front().collider_b == resting.collider,
              "removing every collider lost the final contact end event");
        compare_step(cached, fresh, 2);
        check(cached.events().empty() && cached.stats().contacts == 0 &&
                  cached.stats().broadphase_pairs == 0,
              "empty world retained ended contacts on later ticks");
    }
}
}
int main() {
    try {
        sleeping_and_mutation();
        stationary_kinematics_and_springs();
        removing_final_colliders();
        std::cout << "contact cache continuation, waking and mutation tests passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
