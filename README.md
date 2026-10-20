# fxp_physics

A standalone C++17 deterministic rigidbody engine for RTS simulation. Physics
scalars and math use `fxp::real` and `fxp_math`; the engine currently uses no
`fxp::wide` and no floating-point arithmetic. Integer types hold IDs, indices,
ticks, random state, and serialized bytes. The viewer is not a dependency.

## Included in 0.1

- Static, kinematic, and dynamic bodies; linear/angular velocity, forces,
  impulses, torque, mass, center of mass, inertia tensors, damping and rotation locks.
- Box, sphere, Y-axis capsule, validated convex hull, static triangle mesh,
  and static heightfield collision shapes.
- Compound colliders with local transforms and parenting on the same body.
- 32 collision layers with a symmetric matrix; independent query masks;
  sensors and deterministic begin/stay/end events.
- Convex overlap and translational shape casts, plus raycast, spherecast and
  oriented boxcast conveniences. Query shape may also be a capsule or convex hull.
- Static/dynamic Coulomb friction, restitution, damped distance springs,
  sleeping islands, and render-pose interpolation.
- Fixed ticks and fixed substeps, stable IDs, deterministic RNG, checksums,
  and versioned self-contained save/load as `std::vector<std::uint8_t>`.
- Headless behavior tests and a frozen per-tick state/query determinism trace.

This is an initial implementation with deliberately bounded inputs, a single
simulation thread, and no viewer integration. It is intended to be measured
against actual game workloads before production use.

## Build

Clone with `--recurse-submodules`, or run `git submodule update --init --recursive`.
The nested dependency is pinned to fxp_math commit
`f18e16b2f6c404da51075ce013315369096cbf98`. Dependency tests are not rebuilt by default.

Windows, including VS 2026, from an ordinary PowerShell:

```powershell
powershell -ExecutionPolicy Bypass -File tools/build-windows.ps1 -Configuration Release
powershell -ExecutionPolicy Bypass -File tools/build-windows.ps1 -Configuration Debug
powershell -ExecutionPolicy Bypass -File tools/build-windows.ps1 -Configuration Release -Portable
.\build\windows-release\fxp_physics_example.exe
.\build\windows-release\fxp_physics_tests.exe --benchmark
```

The helper locates Visual Studio, initializes its x64 developer environment, and
uses its bundled CMake/Ninja. Alternatively use the Windows presets from a VS
developer shell. CMake 3.20 supports manual configure; presets require 3.21+.

Linux with GCC or Clang, CMake and Ninja:

```sh
cmake --preset linux-release
cmake --build --preset linux-release
ctest --preset linux-release
./build/linux-release/fxp_physics_example
./build/linux-release/fxp_physics_tests --benchmark
```

Use `-DCMAKE_CXX_COMPILER=clang++` in a separate build directory to select Clang.
`FXP_PHYSICS_SANITIZERS=ON` enables address/undefined-behavior sanitizers on GCC/Clang.
CI checks Windows/MSVC and Linux/GCC/Clang, Debug/Release, and native/portable
integer arithmetic against **the same** reference file. Linux CI is provided;
this workspace's implementation was exercised locally on Windows.

Embed in a game:

```cmake
add_subdirectory(path/to/fxp_physics)
target_link_libraries(my_game PRIVATE fxp::physics)
```

Tests/examples default off when embedded. An existing `fxp::fxp` target is reused;
the application must provide exactly the pinned math version and consistent math
backend definitions across translation units.

## Basic API

```cpp
#include <fxp_physics/physics.hpp>
#include <cassert>
using namespace fxp_physics;

World world(Config{}, 1234); // fixed 30 Hz, two substeps, deterministic seed
BodyDesc body;
body.type = BodyType::kinematic;
const BodyId unit = world.create_body(body);

ColliderDesc collider;
collider.body = unit;
collider.shape = Shape::capsule(0.4_r, 0.5_r);
world.create_collider(collider);

world.step([&](World& fixed) {
    Pose target = fixed.body(unit).properties.pose;
    target.position.x += 0.05_r;
    fixed.set_kinematic_target(unit, target);
});
world.dispatch_events([](const CollisionEvent& event) {
    // Record game commands for a later fixed update here.
});
world.render(0.5_r, [](BodyId id, const Pose& pose) {
    // Convert poses for the renderer at this boundary only.
});

Bytes bytes = world.save();
World restored;
restored.load(bytes);
assert(restored.checksum() == world.checksum());
```

All units are metres, seconds, kilograms and radians. Coordinates are right
handed with Y up. `Pose` is the body's frame; `center_of_mass` is relative to it.
The engine keeps integration about the center of mass when rotation changes.
`BodyState`/`ColliderState` are copies, not mutable references into the world.

`step()` advances exactly one configured tick. There is no wall-clock input or
automatic catch-up policy. The host owns scheduling and its ordered command log.
Mutation is permitted between steps or inside the fixed-update callback. Event
and render callbacks are read-only and reject mutations, including RNG consumption.
Do not move/destroy the world or retain references to callback pose arguments
while a callback is executing. The world and its queries are not thread-safe.

`set_kinematic_target` reaches the target over the next tick, derives contact
velocities from that motion, and stops afterward. `set_velocity` instead persists
kinematic velocity. Kinematic bodies push dynamic bodies; they are not themselves
stopped or separated by collisions. An RTS movement controller should use casts
and overlaps to decide valid kinematic targets.

Forces last one tick, across all its substeps. Impulses act immediately at the
fixed-update boundary. Teleportation invalidates contact impulses and resets both
interpolation poses. Interpolation uses the two latest completed poses and thus
adds one tick of display latency; render calls never alter physics state.

## Colliders, parenting and mass

A body owns any number of colliders. A collider may parent another collider on
that same body, with a maximum hierarchy depth of 64. Collider removal removes
its descendants; body removal removes its colliders and connected springs.
`set_collider_parent` preserves the body-space pose by default. Changes to collider
poses happen at fixed boundaries, recalculate automatic inertia, invalidate old
contact impulses, and wake affected bodies. They do not create independent child
rigidbodies or physically simulated joints. Animated collider changes do not
impart a separately computed child velocity.

By default, the center of mass is the body origin and inertia uses the bounding
box of all nonsensor colliders, including a parallel-axis shift to the specified
center of mass. This is an approximation, especially for capsules and sparse
compounds. Supply a positive-definite `Inertia` tensor and an explicit center of
mass for authored mass properties. Collider changes retain the chosen center of
mass. Gyroscopic torque and automatic exact volume integration are not included.

Convex hulls require a closed manifold of outward-wound convex planar faces;
the factory validates this. It does not compute a hull from a point cloud.
Triangle meshes are two-sided surfaces, not filled volumes. Heightfields use a
fixed diagonal per cell and support rotation/translation as static colliders.
Mesh/heightfield assets are immutable; replace their collider to change terrain.
Shared-edge smoothing, terrain holes, deformable terrain, and mesh cooking tools
are not part of this version. Low-poly terrain and modest hulls are recommended.

## Contacts, queries and springs

The broadphase uses a deterministic median AABB tree. Changed bounds refit the
tree; structural changes and every 32 refits rebuild its spatial partition.
Query candidates remain sorted by collider ID regardless of tree topology.
Unchanged colliders retain their bounds and cooked convex geometry.
Static meshes have a cooked local BVH; terrain queries traverse affected cells.
Primitive and hull contact tests use analytic segment distances and separating
axes. Face clipping creates up to four contact points per convex pair; triangle
contacts may produce several manifolds. The solver uses warm-started sequential
impulses with two friction tangents and fixed iteration counts. Friction combines
by geometric mean; restitution combines by maximum. Springs use a soft distance
constraint parameterized by frequency and damping ratio, with a force cap.

Collision events are ordered by collider pair. A pair is observed if it occurs
in any substep during the tick; its last observed point/normal is reported. A
pair absent for the next complete tick generates `end`. Removal and filtering
changes produce their events on the following step. Sensor-mode changes emit
`end` then `begin`. Static/static and same-body pairs do not generate events.
Physical contact events include the configured contact margin; sensor events
require actual touching/overlap. Events contain IDs even if an owner was removed.

Overlap returns one hit per collider, in collider-ID order; the lowest matching
triangle is reported for meshes/terrain. Touching counts as overlap. Cast results
are sorted by fraction, collider ID, then triangle ID, and may contain multiple
triangle hits for one collider. `hits.front()` is the closest hit when nonempty;
`!hits.empty()` provides an any-hit answer. Cast orientation remains fixed.
Normals point outward from the target. Initial overlap returns fraction zero.
Rays normalize a nonzero direction; distance and sweep displacement must be bounded.
The query mask is independent of the simulation's collision matrix.

Casts use conservative advancement with a 64-iteration limit and a 0.00001 m
contact tolerance. Always inspect `CastResult::complete`: false reports a
candidate that did not converge, and must not be interpreted as a clear path.
These queries are not general continuous collision detection for moving bodies.
There is no rotational sweep. Concave query shapes are not supported; convex
queries against concave targets are supported.

## Determinism and bytes

For the same engine/math version, assets, configuration, seed, initial state,
and ordered fixed-update calls, supported platforms must produce identical state,
query outputs, and events. Different creation orders or tickrates are different
inputs. Object IDs increase monotonically and are never reused. Containers are
iterated in ID order; broadphase ties use input IDs; contact and island ordering
are explicit. Sleeping counts completed ticks and propagates across connected
dynamic bodies. Thread scheduling, pointer values and wall clocks never decide
simulation results. RNG is SplitMix64 with specified rejection sampling and a
32-bit fractional conversion for `random_unit()`.

Save/load is permitted at completed tick boundaries. The format stores raw fixed
values in explicit little-endian fields, with a magic number, format/math version
and a 64-bit FNV-1a integrity hash. It includes shapes, body properties, previous
render poses, kinematic targets, pending forces, cached contact impulses, springs,
sleep counters, active event pairs, recent events, RNG state, layers and next IDs.
Host callbacks and host-owned game/command state are not included. Asset sharing
does not affect bytes: shapes are serialized by value for each collider.
Snapshots have a 256 MiB size limit. This initial encoding favors a simple contract
over compression or deduplicating large assets.

Load validates into a temporary world, including IDs, references, shapes and
numeric invariants, then swaps it into place. Failed loads leave the original
world intact. Incompatible versions are rejected, with no implicit migration.
`checksum()` is the hash of the canonical snapshot payload, not raw struct memory.
It includes cached state needed for exact continuation. It is a desync/integrity
check, not a cryptographic authenticity mechanism. Calling save/checksum does not
consume RNG or change state. Rendering/query frequency does not affect the hash.

## Numeric and practical limits

Use `fxp::real` literals, integers or explicit deterministic import conversions.
There is no float/double path in physics. Math saturation is not a substitute for
a supported physical domain. Current checked input limits include:

| Input | Limit |
| --- | --- |
| Body position | each component within +/-8192 m |
| Primitive half extents / convex vertex coordinates | up to 128 m |
| Primitive radius / box half extent | at least 0.005 m |
| Collider offsets and center-of-mass offsets | each component within +/-128 m |
| Dynamic mass | 0.01 to 10,000 kg |
| Linear / angular speed | configured cap, at most 128 m/s / 32 rad/s |
| Applied force / torque / impulse | each component within +/-1,000,000 |
| Tick timestep / substeps | 0.001 to 0.1 seconds / 1 to 16 |
| Convex hull | 128 vertices, 256 faces |
| Mesh / heightfield | 1,000,000 vertices or samples; mesh 1,000,000 triangles |

These are input bounds, not a promise that every extreme combination is stable.
Ill-conditioned or excessive inertia is rejected. Favor ordinary RTS scales,
moderate mass ratios, and small convex hulls. Use casts for fast projectiles until
CCD is added. Bounds/geometry failures during a step indicate an unsupported
simulation; restore a prior save rather than continuing a partially executed tick.
Host exceptions in fixed-update callbacks similarly leave their already-applied
commands in place; callbacks should normally not throw.

Cloth, general joints/ragdolls, continuous collision detection, rotating casts,
runtime fracture, multithreading, and rollback APIs are not implemented.

## Validation and performance

The tests cover analytic query distances and normals, convex pair combinations,
rotated casts, layers/sensors, parenting/cycle rejection, kinematic pushing,
interpolation isolation, forces/torques, resting/sleeping, stacked boxes, friction,
springs, terrain/mesh support, mutation guards, corrupted/truncated saves, and
byte-exact save/load continuation with different render/query schedules. Additional
regressions cover transient substep sensor crossings, elastic equal-mass impacts,
rotated bounds, asset-sharing independence, and broadphase results against a
brute-force sphere-overlap oracle at negative and positive coordinates.

`tests/reference.txt` contains 100 frozen rows of tick, full-state checksum and
query-result checksum for a mixed scene with creation/removal, forces, RNG,
kinematic targets, materials, triggers, terrain and springs. Normal tests only
read it. It was generated explicitly with the portable scalar backend. Changes
to this file require reviewing the simulation change and verifying all backends;
never regenerate it as part of configure/build/test.

```powershell
# Only for deliberate reference maintenance after reviewing a simulation change:
.\build\windows-release-portable\fxp_physics_tests.exe --write-reference tests/reference.txt
```

`fxp_physics_benchmark` runs without the viewer, window system or downloaded
dependencies. It reproduces the viewer's physics playground and provides sparse
kinematic/dynamic spheres and eight-box stacks at 100, 500, 1,000, 5,000 and
10,000 moving bodies. All use 60 Hz, four substeps and 12 velocity iterations.
Build it with `FXP_PHYSICS_BUILD_BENCHMARKS=ON`; tests also build it for replay
validation.

```sh
./build/linux-release/fxp_physics_benchmark --scene playground --ticks 240 --samples 3
./build/linux-release/fxp_physics_benchmark --scene dynamic --bodies 10000 --ticks 120 --samples 3 --queries 100
./build/linux-release/fxp_physics_benchmark --suite scaling --ticks 10 --samples 2 > scaling.csv
```

CSV output includes mean, median, p95 and maximum tick time, contact/candidate
counts and full-state/query checksums. Construction, statistics, hashing and
trace I/O are outside the timed region. Optional overlap queries are timed.
Each sample recreates the same scene and must end with identical checksums.
Use `--warmup N` to measure a later simulation interval; timings do not control
the number of simulated ticks.

`tests/playground_reference.csv` adds 240 frozen state/query/event rows captured
from the original engine before these optimizations. CI checks both reference
files on GCC, Clang, MSVC and macOS ARM, Debug/Release, native/portable arithmetic.
Release artifacts also compare all five populations across platforms and
backends. Timing numbers are reported, not used as flaky CI pass/fail limits.
The viewer repository has its own workflow because nested submodule workflows
do not run on a viewer push.

See [measured performance and reproduction details](docs/performance.md).
The arithmetic optimizations also modify the nested `fxp_math` dependency; carry
those changes with this engine when updating submodule revisions. They preserve
exact rounding, saturation and frozen results without floating-point arithmetic.
The original `fxp_physics_tests --benchmark` command remains available.
