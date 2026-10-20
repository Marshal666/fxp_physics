#include <fxp_physics/physics.hpp>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
using namespace fxp_physics;
using Clock = std::chrono::steady_clock;

struct Options {
    std::string scene = "playground", suite, trace, reference;
    unsigned bodies = 1000, ticks = 240, samples = 3, warmup = 0, queries = 0;
};

unsigned number(const std::string &text, const std::string &option) {
    if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
        throw std::invalid_argument(option + " requires a nonnegative integer");
    const auto value = std::stoull(text);
    if (value > std::numeric_limits<unsigned>::max())
        throw std::invalid_argument(option + " is too large");
    return static_cast<unsigned>(value);
}

Options options(int argc, char **argv) {
    Options o;
    bool samples_set = false, bodies_set = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help") {
            std::cout <<
                "fxp_physics_benchmark [--scene playground|kinematic|dynamic|stack]\n"
                "  --suite scaling  Run kinematic, dynamic and stack at 100, 500,\n"
                "                   1000, 5000 and 10000 moving bodies\n"
                "  --bodies N       Moving bodies (default 1000; stack default 128)\n"
                "  --ticks N        Measured ticks per sample (default 240)\n"
                "  --samples N      Independently recreated worlds (default 3)\n"
                "  --warmup N       Untimed initial ticks per sample (default 0)\n"
                "  --queries N      Timed sphere overlaps after each tick (default 0)\n"
                "  --trace FILE     Write per-tick state/query/event hashes; one sample\n"
                "  --check-reference FILE  Require an exact frozen trace; one sample\n"
                "CSV reports tick+query timings in microseconds. Setup, stats, hashes,\n"
                "reference checking and trace I/O are outside the timed region.\n"
                "All scenes use 60 Hz, four substeps and 12 solver iterations.\n";
            std::exit(0);
        }
        if (i + 1 == argc)
            throw std::invalid_argument("missing value for " + arg);
        const std::string value = argv[++i];
        if (arg == "--scene") o.scene = value;
        else if (arg == "--suite") o.suite = value;
        else if (arg == "--trace") o.trace = value;
        else if (arg == "--check-reference") o.reference = value;
        else if (arg == "--bodies") { o.bodies = number(value, arg); bodies_set = true; }
        else if (arg == "--ticks") o.ticks = number(value, arg);
        else if (arg == "--samples") { o.samples = number(value, arg); samples_set = true; }
        else if (arg == "--warmup") o.warmup = number(value, arg);
        else if (arg == "--queries") o.queries = number(value, arg);
        else throw std::invalid_argument("unknown option " + arg);
    }
    if (o.scene != "playground" && o.scene != "kinematic" &&
        o.scene != "dynamic" && o.scene != "stack")
        throw std::invalid_argument("unknown scene " + o.scene);
    if (!o.suite.empty() && o.suite != "scaling")
        throw std::invalid_argument("unknown suite " + o.suite);
    if (!o.suite.empty() && (bodies_set || o.scene != "playground" ||
                            !o.trace.empty() || !o.reference.empty()))
        throw std::invalid_argument("--suite scaling cannot combine with scene, bodies or trace options");
    if (o.scene == "stack" && !bodies_set) o.bodies = 128;
    if (o.scene == "playground" && bodies_set)
        throw std::invalid_argument("playground has a fixed population; omit --bodies");
    if (!o.ticks || !o.samples || !o.bodies || o.bodies > 10000 ||
        o.ticks > 100000 || o.warmup > 100000 || o.samples > 1000 || o.queries > 10000)
        throw std::invalid_argument("expected bodies 1..10000, ticks 1..100000, samples 1..1000, "
                                    "warmup 0..100000, queries 0..10000");
    if (!o.trace.empty() || !o.reference.empty()) {
        if (samples_set && o.samples != 1)
            throw std::invalid_argument("traces require --samples 1");
        o.samples = 1;
    }
    if (!o.trace.empty() && o.trace == o.reference)
        throw std::invalid_argument("trace output must not overwrite its reference");
    return o;
}

BodyId add(World &w, const ShapePtr &shape, Vec3 position,
           BodyType type = BodyType::dynamic, Material material = {},
           std::uint64_t user_data = 0) {
    BodyDesc body;
    body.type = type;
    body.pose.position = position;
    body.user_data = user_data;
    ColliderDesc collider;
    collider.body = w.create_body(body);
    collider.shape = shape;
    collider.material = material;
    w.create_collider(collider);
    return collider.body;
}

World playground() {
    Config config;
    config.timestep = 1_r / 60_r;
    config.substeps = 4;
    config.velocity_iterations = 12;
    World world(config, 1234);
    // Mirror examples/physics_playground.yaml, including viewer object IDs.
    // Every decimal is rounded once, exactly as the YAML decimal importer does.
    std::uint64_t object = 2; // Object 1 is the sun, with no physics body.
    add(world, Shape::box({15, 0.5_r, 13}), {0, -0.5_r, 0},
        BodyType::static_body, {0.7_r, 0.5_r, 0}, object++);
    BodyDesc projectile;
    projectile.pose.position = {-4, 1, 8};
    projectile.mass = 12;
    projectile.linear_velocity = {0, 0, -16};
    projectile.linear_damping = 0.01_r;
    projectile.user_data = object++;
    ColliderDesc ball;
    ball.body = world.create_body(projectile);
    ball.shape = Shape::sphere(0.75_r);
    ball.material = {0.4_r, 0.25_r, 0.1_r};
    world.create_collider(ball);
    const auto cube = Shape::box({0.5_r, 0.5_r, 0.5_r});
    const Real wall_x[] = {-6.04_r, -5.02_r, -4_r, -2.98_r, -1.96_r};
    const Real height[] = {0.50_r, 1.52_r, 2.54_r, 3.56_r, 4.58_r, 5.60_r, 6.62_r, 7.64_r};
    for (unsigned y = 0; y < 4; ++y)
        for (const auto x : wall_x)
            add(world, cube, {x, height[y], 0}, BodyType::dynamic, {}, object++);
    const Real tower_x[] = {3.50_r, 3.72_r, 3.94_r, 4.16_r, 4.38_r, 4.60_r, 4.82_r, 5.04_r};
    for (unsigned y = 0; y < 8; ++y)
        add(world, cube, {tower_x[y], height[y], -1}, BodyType::dynamic,
            {0.65_r, 0.45_r, 0}, object++);
    add(world, Shape::box({15, 0.75_r, 0.25_r}), {0, 0.75_r, -12.75_r},
        BodyType::static_body, {}, object++);
    add(world, Shape::box({15, 0.75_r, 0.25_r}), {0, 0.75_r, 12.75_r},
        BodyType::static_body, {}, object++);
    add(world, Shape::box({0.25_r, 0.75_r, 12.5_r}), {-14.75_r, 0.75_r, 0},
        BodyType::static_body, {}, object++);
    add(world, Shape::box({0.25_r, 0.75_r, 12.5_r}), {14.75_r, 0.75_r, 0},
        BodyType::static_body, {}, object++);
    return world;
}

World scenario(const Options &o) {
    if (o.scene == "playground") return playground();
    Config config;
    config.timestep = 1_r / 60_r;
    config.substeps = 4;
    config.velocity_iterations = 12;
    if (o.scene != "stack") config.gravity = {};
    World world(config, 1234);
    const auto shape = o.scene == "stack" ? Shape::box({0.5_r, 0.5_r, 0.5_r})
                                         : Shape::sphere(0.4_r);
    if (o.scene == "stack")
        add(world, Shape::box({100, 0.5_r, 100}), {0, -0.5_r, 0}, BodyType::static_body);
    for (unsigned i = 0; i < o.bodies; ++i) {
        BodyDesc body;
        body.type = o.scene == "kinematic" ? BodyType::kinematic : BodyType::dynamic;
        if (o.scene == "stack") {
            const auto column = i / 8;
            body.pose.position = {Real(column % 50) * 2 - 49,
                                  Real(50 + 101 * (i % 8)) / 100,
                                  Real(column / 50) * 2 - 49};
        } else {
            body.pose.position = {Real(i % 100) * 2, 0, Real(i / 100) * 2};
            body.linear_velocity = {0.5_r, 0, 0};
            body.allow_sleep = false;
            body.linear_damping = 0;
        }
        ColliderDesc collider;
        collider.body = world.create_body(body);
        collider.shape = shape;
        world.create_collider(collider);
    }
    return world;
}

struct Hash {
    std::uint64_t value = UINT64_C(14695981039346656037);
    void mix(std::uint64_t x) {
        for (unsigned i = 0; i < 8; ++i) {
            value ^= (x >> (i * 8)) & 255u;
            value *= UINT64_C(1099511628211);
        }
    }
    void mix(Real x) { mix(static_cast<std::uint64_t>(x.raw_value())); }
    void mix(Vec3 v) { mix(v.x); mix(v.y); mix(v.z); }
    void mix(const std::vector<QueryHit> &hits) {
        mix(hits.size());
        for (const auto &hit : hits) {
            mix(hit.body); mix(hit.collider); mix(hit.triangle);
            mix(hit.fraction); mix(hit.distance); mix(hit.point); mix(hit.normal);
            mix(hit.started_overlapping ? 1u : 0u);
        }
    }
};

std::uint64_t queries(const World &world, const Shape &probe) {
    Hash hash;
    hash.mix(world.overlap(probe, {{0, 1, 0}, {}}, {0xffffffffu, 0, true}));
    const auto ray = world.raycast({-4, 10, 0}, {0, -1, 0}, 20, {0xffffffffu, 0, true});
    hash.mix(ray.complete ? 1u : 0u);
    hash.mix(ray.hits);
    return hash.value;
}

std::uint64_t events(const World &world) {
    Hash hash;
    hash.mix(world.events().size());
    for (const auto &event : world.events()) {
        hash.mix(event.tick); hash.mix(event.collider_a); hash.mix(event.collider_b);
        hash.mix(event.body_a); hash.mix(event.body_b);
        hash.mix(static_cast<std::uint64_t>(event.type));
        hash.mix(event.sensor ? 1u : 0u);
        hash.mix(event.point); hash.mix(event.normal);
    }
    return hash.value;
}

void benchmark(const Options &o, bool header = true) {
    std::ofstream trace;
    std::ifstream reference;
    if (!o.trace.empty()) {
        trace.open(o.trace);
        if (!trace) throw std::runtime_error("cannot open trace output: " + o.trace);
    }
    if (!o.reference.empty()) {
        reference.open(o.reference);
        if (!reference) throw std::runtime_error("cannot open reference: " + o.reference);
    }
    const auto probe = Shape::box({2, 2, 2});
    const auto timed_probe = Shape::sphere(0.4_r);
    std::uint64_t expected_checksum = 0, expected_queries = 0;
    if (header)
        std::cout << "scene,bodies,ticks,warmup,queries,sample,mean_us,p50_us,p95_us,max_us,"
                     "peak_pairs,peak_contacts,peak_awake,checksum,query_checksum,query_hits\n";
    for (unsigned sample = 0; sample < o.samples; ++sample) {
        World world = scenario(o);
        for (unsigned tick = 0; tick < o.warmup; ++tick) world.step();
        std::vector<double> times;
        times.reserve(o.ticks);
        std::uint32_t peak_pairs = 0, peak_contacts = 0, peak_awake = 0;
        std::uint64_t query_hits = 0;
        for (unsigned tick = 0; tick < o.ticks; ++tick) {
            const auto begin = Clock::now();
            world.step();
            for (unsigned q = 0; q < o.queries; ++q)
                query_hits += world.overlap(*timed_probe,
                    {{Real(q % 100) * 2, 0, Real(q / 100) * 2}, {}}).size();
            times.push_back(std::chrono::duration<double, std::micro>(Clock::now() - begin).count());
            const auto stats = world.stats();
            peak_pairs = std::max(peak_pairs, stats.broadphase_pairs);
            peak_contacts = std::max(peak_contacts, stats.contacts);
            peak_awake = std::max(peak_awake, stats.awake_dynamic_bodies);
            if (!o.trace.empty() || !o.reference.empty()) {
                std::ostringstream row;
                row << o.scene << ',' << world.tick() << ',' << world.checksum() << ','
                    << queries(world, *probe) << ',' << events(world);
                if (!o.trace.empty()) trace << row.str() << '\n';
                if (!o.reference.empty()) {
                    std::string expected;
                    if (!std::getline(reference, expected) || expected != row.str())
                        throw std::runtime_error("reference mismatch at tick " + std::to_string(world.tick()) +
                                                 ": expected " + expected + "; actual " + row.str());
                }
            }
        }
        const auto checksum = world.checksum(), query_checksum = queries(world, *probe);
        if (sample && (checksum != expected_checksum || query_checksum != expected_queries))
            throw std::runtime_error("independent benchmark samples diverged");
        expected_checksum = checksum;
        expected_queries = query_checksum;
        const auto mean = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
        std::sort(times.begin(), times.end());
        const auto percentile = [&](unsigned p) { return times[(times.size() * p + 99) / 100 - 1]; };
        std::cout << o.scene << ',' << world.stats().bodies << ',' << o.ticks << ',' << o.warmup
                  << ',' << o.queries << ',' << sample + 1 << ',' << std::fixed << std::setprecision(3)
                  << mean << ',' << percentile(50) << ',' << percentile(95) << ',' << times.back()
                  << ',' << peak_pairs << ',' << peak_contacts << ',' << peak_awake << ','
                  << checksum << ',' << query_checksum << ',' << query_hits << '\n';
    }
    if (!o.reference.empty()) {
        std::string extra;
        if (std::getline(reference, extra)) throw std::runtime_error("reference has extra rows");
    }
    if (!o.trace.empty() && !trace) throw std::runtime_error("failed to write trace");
}
} // namespace

int main(int argc, char **argv) {
    try {
        auto o = options(argc, argv);
        if (o.suite.empty()) benchmark(o);
        else {
            bool header = true;
            for (const auto *scene : {"kinematic", "dynamic", "stack"})
                for (const unsigned count : {100u, 500u, 1000u, 5000u, 10000u}) {
                    o.scene = scene;
                    o.bodies = count;
                    benchmark(o, header);
                    header = false;
                }
        }
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "benchmark: " << e.what() << '\n';
        return 1;
    }
}
