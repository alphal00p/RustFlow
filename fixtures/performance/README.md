# Differential-equation benchmark provenance

These two files are verbatim copies from AMFlow 2.0 commit
[`26005517a288086c4cb4d1b26d829691bc088485`](https://gitlab.com/multiloop-pku/amflow/-/tree/26005517a288086c4cb4d1b26d829691bc088485):

| Local file | Upstream path | SHA-256 |
|---|---|---|
| `upstream12-diffeq.wl` | `examples/differential_equation_solver_cpp/input_files/diffeq` | `828dd1ae41291007fd3b37b664e695d67e0daebde3e9330f0e178f3d39cc40e3` |
| `upstream12-boundary.wl` | `examples/differential_equation_solver_cpp/RegularExpansion/bc` | `fcbe018f60ca9cf1bf909dec36d3608993d22fe5c8bb97a86df432870c1f92c6` |

The matrix has 12 masters. The supplied boundary is at eta=1/2 and
epsilon=101/9999900. Its component precision annotations are 180 decimal digits;
these describe the upstream input, not an independent accuracy verification.
The benchmark continues to eta=1/10+i/5. This endpoint is chosen for this benchmark,
not taken from an upstream timing claim.

The fixtures are used only as supplied boundaries for the standalone DE comparison.
They are **not** used by the automatic integral evaluator or the mandatory four-target
acceptance test. This comparison does not validate recursive boundary generation.

The upstream MIT license is preserved in
[`../amflow-2.0/LICENSE.md`](../amflow-2.0/LICENSE.md); see also the root `NOTICE`.
`scripts/benchmark_de.py` checks both hashes and parses only lists and numeric
literals. It never executes Wolfram code.
