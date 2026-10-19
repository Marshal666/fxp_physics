#include <fxp_physics/physics.hpp>
#include <iostream>

int main() {
    using namespace fxp_physics;
    World world({}, 1234);
    BodyDesc ground;
    ground.type = BodyType::static_body;
    ground.pose.position = {0, -0.5_r, 0};
    ColliderDesc collider;
    collider.body = world.create_body(ground);
    collider.shape = Shape::box({20, 0.5_r, 20});
    world.create_collider(collider);
    BodyDesc debris;
    debris.pose.position = {0, 4, 0};
    const auto body = world.create_body(debris);
    collider.body = body;
    collider.shape = Shape::box({0.5_r, 0.5_r, 0.5_r});
    world.create_collider(collider);
    for (unsigned i = 0; i < 300; ++i)
        world.step();
    const auto bytes = world.save();
    World restored;
    restored.load(bytes);
    std::cout << "Tick " << world.tick() << ", height "
              << world.body(body).properties.pose.position.y.to_string() << ", asleep "
              << world.body(body).sleeping << ", bytes " << bytes.size() << ", checksum " << world.checksum()
              << '\n';
    return restored.checksum() == world.checksum() ? 0 : 1;
}
