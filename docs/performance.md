# Physics performance measurements

Measured on 2026-10-20 with GCC 13.3, Linux x86-64, Intel Core Ultra 7 165H,
Release (`-O3 -DNDEBUG`), pinned to logical CPU 1. The baseline was physics
`b30e2fde07ea6bd31a9c3cddcb25f7105fd0d068` with math
`f18e16b2f6c404da51075ce013315369096cbf98`. The optimized column uses the
working-tree changes accompanying this report. These are local measurements,
not Windows or macOS performance claims.

The playground averages **8.185 ms/tick instead of 14.244 ms** (1.74x faster).
Sparse moving populations of 10,000 bodies average **9.784 ms kinematic** and
**12.913 ms dynamic**, including 100 overlap queries per tick. Dense stacks remain
expensive: **10,000 awake boxes with 40,000 contacts do not meet 60 Hz**.

All scenes retain 60 Hz, four substeps and 12 velocity iterations. No solver
quality settings were lowered. These are physics timings; rendering and host
game logic still need their own share of a 16.667 ms frame budget.

## Results

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

## Changes

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

## Reproduce

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

## Determinism and validation

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
