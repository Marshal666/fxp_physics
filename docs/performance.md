# Physics performance measurements

## Current: stack and parallel solver optimization

Measured on 2026-10-20, Linux x86-64, Intel Core Ultra 7 165H, GCC 13.3
Release (`-O3 -DNDEBUG`). The baseline is physics commit `f944ceb` (the previous
optimization pass); both builds use math commit
`028198c02fc88852c22d484120ee7da12a62e6b6`. The after results use the implementation
accompanying this report. No viewer scheduling or rendering changes are included.

**10,000 awake boxes in separated towers improved from 2,620.126 to 371.209 ms
per tick with eight workers, about 7.06x faster. This still misses the 50 ms budget
for 20 Hz physics.** Once all those boxes sleep, the cost is 3.868 ms/tick. Disturbing
one eight-box tower costs 33.817 ms/tick on average, with occasional ticks above
50 ms. A touching pile remains considerably more expensive because its connected
constraint island is solved serially.

These are headless physics measurements, not a 60 FPS viewer claim. Four substeps
and 12 velocity iterations are retained throughout. At 20 Hz, each measured tick
advances 50 ms of simulation regardless of its execution time. Existing viewer
stress YAML files still specify 60 Hz; set `physics.tick_hz: 20` when reproducing
the 20 Hz configuration interactively.

### Awake scaling at 20 Hz

Times below are milliseconds per complete tick. Before is one measured run;
after values are medians of the means from three independently recreated worlds.
Each run measures the first 12 ticks, with no warmup. These are initial contact
and settling costs, not sustained long-running averages. All boxes stay awake
through this interval. Each scene adds one static ground to the dynamic count.

| Dynamic boxes | Before, 1 worker | After, 1 worker | After, 4 workers | After, 8 workers | Before / 8 workers |
| --- | ---: | ---: | ---: | ---: | ---: |
| 100 | 32.163 | 11.539 | 7.618 | 7.396 | 4.35x |
| 500 | 143.212 | 51.974 | 18.145 | 21.840 | 6.56x |
| 1,000 | 272.374 | 102.069 | 34.345 | 49.083 | 5.55x |
| 5,000 | 1305.906 | 501.521 | 173.679 | 212.737 | 6.14x |
| 10,000 | 2620.126 | 1017.260 | 399.059 | 371.209 | 7.06x |

The `stack` scene uses separate eight-box towers, 1,250 towers at 10,000 boxes.
It reaches 40,000 contact points and permits parallel island solving. Both
one-worker builds were pinned to logical CPU 1. Parallel runs were unpinned, and
all timing commands ran sequentially without concurrent builds or test jobs.
This CPU has heterogeneous cores; scheduling and clocks affect results, and more
workers do not always improve performance. Four workers were faster than eight
at some smaller populations. These measurements are local, not Windows or macOS
performance claims.

### Sleeping and waking 10,000 boxes

The settled and disturbed phases run 600 untimed warmup ticks, then measure 120
ticks per sample with eight workers and three independent samples. The disturbed
phase applies a horizontal impulse of 5 to one central bottom box. Exactly eight
bodies remain awake throughout this measured phase; adjacent towers stay asleep.
All three scenarios retain the same 10,000 boxes and ground.

Mean and p95 are medians of the per-sample statistics. Maximum is the slowest
individual tick across all three samples. The awake row uses the 12-tick window
above, so these rows intentionally describe different simulation intervals.

| State | Mean ms | p95 ms | Maximum ms | Peak awake bodies |
| --- | ---: | ---: | ---: | ---: |
| All awake | 371.209 | 422.909 | 450.705 | 10000 |
| All asleep | 3.868 | 4.502 | 5.125 | 0 |
| One tower disturbed | 33.817 | 42.451 | 58.825 | 8 |

A low average does not guarantee a deadline: the disturbed case has a 58.825 ms
worst tick, above the 50 ms physics budget. Blocking ticks also consume viewer
frame time; rendering at 60 FPS needs a separate 16.667 ms frame budget.

Matching before/after sleeping comparisons use one worker, 600 warmup ticks and
12 measured ticks (one baseline sample, three after samples):

| Dynamic boxes | Phase | Before mean ms | After mean ms | Speedup |
| --- | --- | ---: | ---: | ---: |
| 100 | settled | 12.701 | 0.012 | 1020.4x |
| 1,000 | settled | 114.734 | 0.118 | 970.6x |
| 100 | disturbed | 16.420 | 1.534 | 10.7x |

No 10,000-body sleeping baseline was measured; do not extrapolate one from these
smaller populations. The raw results also include 1,000 disturbed boxes after
optimization, for which a clean matching baseline was not retained.

### Connected piles at 20 Hz

`pile` places the same eight-box towers directly beside one another. Contacts
connect them into a single dynamic island. Contact generation can still run in
parallel, but its sequential impulse solve cannot. These after measurements use
eight workers, the first 12 ticks, and three samples. Every body remains awake.
There is no before measurement for this scene.

| Dynamic boxes | Mean ms | p95 ms | Maximum ms | Peak contact points |
| --- | ---: | ---: | ---: | ---: |
| 1,000 | 693.637 | 851.207 | 914.195 | 31116 |
| 10,000 | 7961.698 | 10174.861 | 10388.921 | 354064 |

Independent-island parallelism cannot make this connected workload meet 20 Hz.
Further work would need to address solving within a large island, for example a
deterministic constraint partition with stable ordering. Any such change needs
its own stack stability and replay validation.

### Playground and sparse scenes at 60 Hz

These are freshly measured with the same benchmark driver on both libraries.
The playground runs 240 ticks. Sparse scenes have 10,000 moving spheres, run
30 ticks and include 100 timed overlap queries per tick. There is no warmup;
all table entries are medians of three sample means. Substeps and iterations
remain four and 12. These rows use 60 Hz, unlike the stack tables above.

| Scene | Before, 1 worker ms | After, 1 worker ms | After, 8 workers ms |
| --- | ---: | ---: | ---: |
| Playground | 7.944 | 2.351 | 2.337 |
| Sparse kinematic | 9.630 | 10.392 | 8.302 |
| Sparse dynamic | 12.721 | 13.578 | 11.389 |

There is a tradeoff in sparse serial workloads: their means increased by about
7–8% with the new contact pipeline. Eight workers improved both sparse cases
relative to the baseline. The playground improved by about 3.4x.

Raw second-pass results: [before](stack-performance-before.csv),
[after](stack-performance-after.csv). CSV rows record the frequency, phase,
workers, warmup, sample length and awake/contact counts. Timing excludes scene
construction, warmup, statistics, state/query/event hashing and trace I/O. The
same driver performs untimed verification queries after every measured tick on
both libraries. The earlier report below used an older driver; compare each
report's own paired measurements.

### Implementation and determinism

- Cache world inertia, contact Jacobians, linear/angular impulse responses and
  effective masses before solver iteration. Use scalar tangent impulses and
  guarded fixed-point reciprocals, with division for unrepresentable or poorly
  resolved reciprocals. Skip inactive constraints.
- Retain contacts for unchanged inactive pairs and skip collision detection and
  solver setup when the entire world is stationary. Preserve sensors and
  begin/stay/end events, and invalidate caches after supported mutations.
- Specialize box pairs with oriented-box separating axes, bounded face clipping
  and supporting-edge contacts. Preserve face contacts for nearly aligned axes,
  and keep generic fallbacks for clipping degeneracies and saturated casts.
- Use a persistent worker pool for broadphase queries, geometry preparation,
  contact generation, inertia preparation and independent-island solving. Merge
  outputs in collider order and keep the constraint order inside every island.
  Cache storage is reused across steps. Shared geometry is prepared before reads.

`World::set_worker_count(1)` selects serial execution; zero, the default, uses
hardware concurrency capped at eight workers. Worker count changes scheduling,
not simulation inputs, and is excluded from snapshots/checksums. Loading preserves
the destination world's worker setting. Internal workers finish before `step()`
returns; concurrent host access to a world is still unsupported. All simulation
arithmetic remains fixed-point integer arithmetic. This pass adds no SSE/AVX
requirement and makes no further math-library changes.

Solver rounding and box contacts intentionally change trajectories. Simulation/
save version is now **2**; version-1 saves are rejected rather than silently
continuing under different rules. Frozen 100-tick mixed-scene and 240-tick
playground references were deliberately regenerated with portable arithmetic and
checked against native output. Normal tests never regenerate them.

All nine physics tests pass on native and portable Release builds and under
portable AddressSanitizer/UBSan. The viewer builds and its two scene/playback
tests pass. Regression coverage includes randomized box geometry, thin/extreme
shapes, resting stacks at 20 Hz, springs and locked rotations, cache invalidation,
sensors and incoming wake-ups, removing every collider, mixed-shape parallel
contacts and byte-exact continuation across worker settings. Twenty benchmark
workloads have identical per-tick state/query/event aggregate hashes and counts
across native/portable arithmetic and one/four workers. Timing samples also agree
across one/four/eight workers where their scenarios overlap.

Local ASan uses `ASAN_OPTIONS=detect_leaks=0` because LeakSanitizer cannot run in
the sandbox. Local ThreadSanitizer cannot start (`unexpected memory mapping`), so
it has no local pass result. A focused Helgrind run found no data-race or lock-order
diagnostics, but did report its condition-variable notification heuristic for
notifications made after unlocking. CI adds a separate TSan job, plus GCC, Clang,
MSVC and macOS ARM, native/portable, Debug/Release checks. Release artifacts compare
20 scenarios across 16 platform/backend/worker outputs. Those remote jobs must
run before claiming Windows/macOS or TSan validation for this change.

### Reproduce this pass

From the physics repository:

```sh
cmake -S . -B build/performance -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DFXP_PHYSICS_BUILD_TESTS=ON -DFXP_PHYSICS_BUILD_BENCHMARKS=ON
cmake --build build/performance --parallel
ctest --test-dir build/performance --output-on-failure

./build/performance/fxp_physics_benchmark --scene stack --bodies 10000 --hz 20 --substeps 4 --iterations 12 --workers 8 --ticks 12 --samples 3
./build/performance/fxp_physics_benchmark --scene stack --bodies 10000 --hz 20 --workers 8 --phase settled --ticks 120 --samples 3
./build/performance/fxp_physics_benchmark --scene stack --bodies 10000 --hz 20 --workers 8 --phase disturbed --ticks 120 --samples 3
./build/performance/fxp_physics_benchmark --scene pile --bodies 10000 --hz 20 --workers 8 --ticks 12 --samples 3
python3 benchmarks/run_validation.py build/performance/fxp_physics_benchmark validation-results
python3 benchmarks/compare_results.py validation-results --expected-files 2
```

Repeat the awake stack command with 100, 500, 1,000 and 5,000 boxes and workers
1/4/8. Prefix one-worker timing commands with `taskset -c 1` on Linux. For the
small sleeping comparisons, use workers 1, ticks 12 and bodies 100/1,000. The
settled/disturbed default warmup is 600 ticks; inspect the awake counts instead
of assuming a specified warmup guarantees sleep. Windows uses `.exe` filenames.
To reproduce the baseline, compile this same benchmark source against the older
library with `FXP_PHYSICS_BENCHMARK_LEGACY_WORKERS` defined and workers 1.

## Historical: first optimization pass at 60 Hz

The report below describes the previous implementation and its original paired
measurements. Its statements about unchanged references and snapshot encoding
apply to that first pass only; the current version-2 changes are documented above.

Measured on 2026-10-20 with GCC 13.3, Linux x86-64, Intel Core Ultra 7 165H,
Release (`-O3 -DNDEBUG`), pinned to logical CPU 1. The baseline was physics
`b30e2fde07ea6bd31a9c3cddcb25f7105fd0d068` with math
`f18e16b2f6c404da51075ce013315369096cbf98`. The optimized column describes the
first pass, subsequently committed as physics `f944ceb`. These are local measurements,
not Windows or macOS performance claims.

The playground averages **8.185 ms/tick instead of 14.244 ms** (1.74x faster).
Sparse moving populations of 10,000 bodies average **9.784 ms kinematic** and
**12.913 ms dynamic**, including 100 overlap queries per tick. Dense stacks remain
expensive: **10,000 awake boxes with 40,000 contacts do not meet 60 Hz**.

All scenes retain 60 Hz, four substeps and 12 velocity iterations. No solver
quality settings were lowered. These are physics timings; rendering and host
game logic still need their own share of a 16.667 ms frame budget.

### Results

Times are milliseconds per complete tick, including optional queries. Before,
after and after-p95 are medians of the respective statistics from three
independent runs. The playground runs 240 ticks; sparse scenes run 30 ticks with
100 overlap queries each; stacks run the first 12 ticks with no extra queries.
There are no warmup ticks. Stack measurements cover initial settling, not a
steady-state or sleeping population. Stacks have the requested number of
dynamic boxes plus one static ground body.

| Scene | Total bodies | Before mean | After mean | Speedup | After p95 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Playground | 34 | 14.244 | 8.185 | 1.74x | 11.761 |
| Sparse kinematic | 100 | 0.241 | 0.078 | 3.08x | 0.090 |
| Sparse kinematic | 500 | 1.138 | 0.389 | 2.93x | 0.488 |
| Sparse kinematic | 1,000 | 2.495 | 0.842 | 2.96x | 1.095 |
| Sparse kinematic | 5,000 | 13.677 | 4.797 | 2.85x | 6.184 |
| Sparse kinematic | 10,000 | 28.711 | 9.784 | 2.93x | 12.714 |
| Sparse dynamic | 100 | 0.389 | 0.106 | 3.69x | 0.125 |
| Sparse dynamic | 500 | 1.936 | 0.529 | 3.66x | 0.628 |
| Sparse dynamic | 1,000 | 4.117 | 1.129 | 3.64x | 1.375 |
| Sparse dynamic | 5,000 | 22.481 | 6.249 | 3.60x | 7.708 |
| Sparse dynamic | 10,000 | 46.804 | 12.913 | 3.62x | 15.771 |
| Stack | 101 | 44.202 | 24.461 | 1.81x | 33.898 |
| Stack | 501 | 207.909 | 108.845 | 1.91x | 150.642 |
| Stack | 1,001 | 409.529 | 209.584 | 1.95x | 286.937 |
| Stack | 5,001 | 2054.922 | 1020.861 | 2.01x | 1390.430 |
| Stack | 10,001 | 4173.208 | 2037.034 | 2.05x | 2764.709 |

Raw runs: [before](performance-before.csv), [after](performance-after.csv).
All 48 corresponding rows have identical full-state hashes, query hashes,
candidate/contact counts, awake counts and query hit counts. The playground
peaks at 70 candidate pairs and 130 contacts; sparse scenes have zero contacts;
stacks reach four contacts per dynamic box. Counts are sampled at completed
tick boundaries, so intermediate substeps may contain more contacts.

These measurements do not guarantee a frame deadline. The largest individual
optimized playground tick in these samples was 16.716 ms; the largest sparse
10,000-dynamic-body tick was 18.923 ms. Rebuilding the tree periodically costs
more than ordinary refits. Larger contact-heavy populations still need further
solver work, and a sparse benchmark cannot establish their frame budget.

### Changes

- Refit existing broadphase nodes instead of rebuilding every substep. Rebuild
  after structural changes or 32 refits; return candidates in collider-ID order.
  Reuse query storage and use the balanced tree's bounded traversal stack.
- Retain bounds and cooked convex geometry until a collider's exact world pose
  changes. Cache edge topology, centered support vertices and normalized edge
  directions. This trades some resident memory for less repeated work.
- Reject separated contact/overlap pairs as soon as an axis proves they cannot
  touch. Accepted contacts preserve the original axis order and tie handling;
  casts still compute the full separation distance.
- Cache quaternion coefficients for angular inertia without changing the order
  of fixed-point operations. Skip exact-zero impulses/rotations, use direct
  lookups for dense IDs, and compress island paths.
- Compare raw squared lengths against the exact rounded-length boundary for
  friction, sleeping and tiny-axis rejection. This avoids square roots without
  substituting rounded `dot(v,v)` comparisons.
- In `fxp_math`, use exact integer Newton square roots on native wide arithmetic;
  portable builds seed a restoring calculation with a 64-bit integer root.
  Specialize fixed multiplication's rounding/saturation and directly normalize
  vectors with one nonzero component. No floating-point seed, approximate
  reciprocal, SSE/AVX requirement or platform-dependent simulation branch was
  introduced. Native and portable paths implement the same integer result.

The public API, snapshot encoding, original frozen reference, solver iteration
order and timestep settings remain unchanged.

### Reproduce

From the physics repository, in a compiler environment with CMake and Ninja:

```sh
cmake -S . -B build/performance -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DFXP_PHYSICS_BUILD_BENCHMARKS=ON
cmake --build build/performance --parallel
ctest --test-dir build/performance --output-on-failure

./build/performance/fxp_physics_benchmark --scene playground --ticks 240 --samples 3
./build/performance/fxp_physics_benchmark --scene kinematic --bodies 10000 --ticks 30 --samples 3 --queries 100
./build/performance/fxp_physics_benchmark --scene dynamic --bodies 10000 --ticks 30 --samples 3 --queries 100
./build/performance/fxp_physics_benchmark --scene stack --bodies 10000 --ticks 12 --samples 3
```

Repeat the last three commands at 100, 500, 1,000 and 5,000 moving bodies to
reproduce the table. On Linux, prefix timed commands with `taskset -c 1` for the
same affinity policy. On Windows, run from an x64 developer shell and use the
`.exe` filenames. `--suite scaling --ticks 10 --samples 2` runs all 15 scalable
cases with one uniform measurement interval for CI.

Construction, statistics, checksums and trace I/O are excluded from timing.
`--queries` adds timed sphere-overlap queries. `--warmup` changes the simulated
interval, so comparisons must use the same value. The deterministic scene
generator in the viewer repository, `tools/generate_physics_stress.py`, produces
matching YAML populations for interactive testing.

### Determinism and validation

The original 100-tick state/query reference is unchanged. A new 240-tick
playground state/query/event reference was captured using the original library
before optimization. All 68 math tests and seven physics tests passed in native
and portable Release builds. All 16 benchmark scene/population combinations also
matched across the two backends, including every requested size. The viewer's
scene/playback tests and portable AddressSanitizer/UBSan checks passed locally.
LeakSanitizer is unsupported by the local sandbox's tracing
environment; local sanitizer runs set `ASAN_OPTIONS=detect_leaks=0`. CI keeps
normal sanitizer settings.

Regression tests cover cache invalidation after motion/removal/load, query
frequency independence, broadphase refits against brute-force bounds, saturated
support queries, accepted/rejected separation thresholds and exact rounding
boundaries. Arithmetic development checks compared millions of cases with the
original algorithms and independent Python integers, including all existing
51,004 math and 5,530 geometry reference rows on native and portable backends.

The root viewer workflow and standalone physics workflow test GCC, Clang, MSVC
and macOS ARM in Debug/Release with native/portable math. Eight Release
platform/backend artifacts compare hashes and counts for all 15 scalable
scenes; wall-clock timings are excluded from equality checks. Windows/macOS
runs must complete on their CI hosts before claiming those platforms were
validated for this change.
