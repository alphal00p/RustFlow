# Native dependent cut inventories

This milestone accepts an overcomplete native HEPKit physical denominator
inventory before constructing independent RustRed families. It preserves the
original positive-energy measure and sums all exact partial-fraction and
numerator weights at the same nonzero epsilon before one Laurent fit.

## Domain and ownership

The admitted measure is a complete positive-energy final state, with exact real
rational quadratic coefficients and Gram entries, independent cut shells,
nonnegative squared masses, and uniform `+i0` virtual directions. Every real-only
uncut pole must pass the existing exact bound away from zero over supported
compact phase space. Admission precedes all partial fractions, cut-zero pruning,
and numerator cancellation. Unsupported prescriptions and interior phase-space
poles retain typed failures even for zero requested integrals.

The common regulator acts on virtual diagonal scalar products. Original cut
forms have rate zero; every retained virtual quadratic has a strictly positive
rate. The adapter checks these rates on original exact forms. It rebuilds the
native inventory from admitted canonical rational values before decomposition,
so an earlier symbolic rank estimate cannot change the physical specialization.
The exact tiny-coefficient regression distinguishes the derivative contribution
`-10^-1000/4` from zero without using numerical tolerances.

HEPKit owns affine partial fractions, power transfer, sector extraction, basis
completion, numerator rewriting, graph cut partitions and momentum routing.
Linnet/Spenso/Idenso remain the graph and tensor owners. RustRed owns cut-aware
IBPs; existing ordinary/cut recursive boundary providers and phase-space
terminals evaluate the independent children. Symbolica owns exact and MPFR
arithmetic. No process-specific numerical terminal, graph parser, reduction
algorithm, fit algorithm, or numerical cache format is added.

Raised cuts use the existing normalized identity `D*C_n=C_(n-1)`, `C_0=0`.
Original oriented cut momenta and positive normalizations survive slot remapping.
Original represented numerator denominator guards survive rational cancellation,
removed required cuts, and empty groups. Exact forbidden epsilon samples still
fail for a zero result. Backend conditions retain their existing owners.

The public grouped entrypoints are `PreparedCutCombination::from_hepkit` and
`GraphIntegral::prepare_cut_combination`. Existing descriptive Python cut methods
and CLI graph requests share them. A bounded native partial-fraction budget is
exposed through the existing CLI option and Python
`max_partial_fraction_states`. Explicit mass slots remain original physical
slots; a nonterminal child that loses every selected slot is rejected. Finite
samples report working precision only; independently refined fits report their
measured accuracy evidence.

## Independent comparisons

At `P²=1`, the massless two-body cut inventory includes the smooth dependent pole
`r²-2`. Its unit value is `-Phi_2/2`; raising the first cut gives
`(1/4-epsilon)*Phi_2`. Tests include reversed energy orientation, denominator
permutations, exact complex numerator weights, removed cuts, and actual native
Diagram/Model/Kinematics objects through Rust, CLI DOT reload and Python.

The nonfactorized example has three massless cuts `r²`, `l²`, `(P-r-l)²` and
virtual poles `k²-2`, `k²-3`, `(k-r-l)²-5`. The virtual bubble depends on the
variable invariant `(r+l)²`; it cannot be replaced by a constant times phase
space. Tests use an independent double-beta/Feynman-parameter series, All and
partial mass placements, additional epsilon samples and increased order and
precision. Equal-mass specialization checks the homogeneous partial-fraction
identity and a raised virtual pole. A separate Laurent test verifies that the
individual bubble UV poles cancel before fitting the combined value.

The optional `upstream_dependent_cut_oracle.py` driver ran the pinned original
AMFlow 2.0 computational sources with Kira on the two independent complete
child families. The archived inputs retain the exact measure, epsilon, physical
point, actual top-level auxiliary positions, requested precision, source hashes
and runtime configuration. This is an upstream oracle for the independent
children; it does not claim upstream acceptance of a raw overcomplete family.

| Independent child, epsilon = 1/13 | AMFlow recorded value | SolveIntegrals seconds |
| --- | --- | ---: |
| `k²-2`, `(k-r-l)²-5` | `0.00318444953599229637406421702536975702` | 129.144518 |
| `k²-3`, `(k-r-l)²-5` | `0.00314557426656282219605916631513760248` | 128.614578 |

Both reference values carry **24 decimal digits of recorded precision**. Their
difference is `-0.00003887526942947417800505071023215454`; cancellation loses
about 2.2 relative digits. `upstream-combination.json` records the nominal
uncertainty propagated from that precision and the independent series value.
Small raw residuals do not upgrade the recorded reference accuracy. Numerical
agreement is asserted at 20 relative digits.

The explicit oracle assertions scale nonzero errors by the expected value;
zero coefficients use an absolute threshold. The fitter's `verified_digits`
retains the existing library convergence convention and is reported separately.
Neither working precision nor a small raw residual is relabeled as a stronger
verified precision claim.

## Native timing and validation

The fresh native process uses one CPU, 20 requested digits, 70 working decimal
digits and expansion order 80. The independent refinement raises these to 90
digits and order 112. Both numerical and persistent exact-system caches are
disabled. The already prepared second sample is not a numerical cache hit.

| Work | Native All | Native explicit slot 5 | Original AMFlow All |
| --- | ---: | ---: | ---: |
| Cold preparation | 1.1126 s | 19.2283 s | Included below |
| Combined sample at epsilon = 1/13 | 1.9430 s | 7.2745 s | Included below |
| Preparation plus first sample | 3.0556 s | 26.5028 s | 257.7591 s, sum of two child jobs |
| Prepared sample at epsilon = 1/17 | 1.8999 s | 7.2298 s | Not measured |

These are qualified single-case timings. Native evaluation starts from the
overcomplete inventory; upstream evaluates its two independent children.
Physical point, normalization and epsilon agree, while reducers, requested
precision and numerical profiles differ. This does not establish general
performance parity. `single-cpu-samples.json` and `numerical-comparison.json`
preserve the complete profile and exact values.

The separate UV-cancelling Laurent reconstruction through epsilon zero takes
175.84 s in the two-thread regression suite, using 33 fit samples, 29 validation
samples and one refinement. It passes both the shared fitter's 20-digit
verification and the independent finite-coefficient comparison. That timing
excludes preparation and is not compared with the upstream finite-sample jobs.

All **68 scoped release tests pass**, covering the exact native inventory unit,
cut schema and IBPs, ten dependent-inventory regressions, native Rust/CLI/Python
objects, partial placements, phase-space terminals, and ordinary HEPKit and
weighted-projection regressions. Strict release all-target `python_stubgen`
Clippy, formatting, diff checks and the optional oracle driver's Python syntax
check also pass. The source-verifying runner records every executable hash and
test count in `verified-tests.json`; final logs are in `logs/`.

The frozen port source digest is
`412f4afb7bc20c57bcf34de7547673b0a6e834e00817cd85b71d01d3f477b9da`.
`source-files.json`, `validation.json` and `frozen-artifacts.json` retain the
source, dependency, toolchain and immutable binary/library provenance. This
gate uses locked Symbolica `c3408e4`, HEPKit `a3d1c8a` and RustRed `7c1ed03`.
Newer owner and Symbolica revisions require their separate integration gate;
they do not retroactively change this scientific evidence.

Development failures are retained separately in `development/`: the initial
fixture omitted native momentum declarations, and the next frozen gate exposed
the corrected scalar-product/dimension ordering issue. Their results are not
included in the 68-test final gate.
