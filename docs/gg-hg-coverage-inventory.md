# Full mixed QCD–EW Higgs-plus-jet coverage inventory

The pinned application at `mg5_higgs_ew_plugin` commit
`89f64b93d0bdbd8ee90eb85737021189229b034a` provides enough scientific data to
validate the complete supplied-system amplitude. The existing production check
covers only its closed first two nonplanar masters. The following inventory is
not a claim that the full amplitude has passed.

| Input | Planar | Nonplanar |
| --- | ---: | ---: |
| Complete differential-system dimension | 48 | 61 |
| Logarithmic letters | 36 | 75 |
| Nonzero constant letter-matrix entries | 1,835 | 2,754 |
| Distinct expanded square-root radicands | 2 | 8 |
| Epsilon orders required by form factors | 0–4 | 1–4 |
| Distinct master indices referenced by form factors | 42 | 57 |
| Invariant permutations referenced by form factors | 4 | 3 |

The original driver evaluates four permutations for each family. Testing every
master through epsilon order four therefore means 2,180 complex coefficients per
mass choice before forming the four complex form factors. Both W and Z choices
must be covered. The 8-root count identifies exact polynomial duplicates; it
does not assert algebraic independence modulo squares.

Every connection entry is an exact constant linear combination of logarithms,
so the native `CanonicalAlgebraicSystem` owner can retain this representation.
No application-specific differential-equation solver or cache is needed. Sparse
scientific matrices, exact letters and precision metadata have been exported to
`target/gg-hg-inventory`. The [machine-readable inventory](../reports/validation/2026-10-04-gg-hg-input-inventory.json)
records source hashes and extraction scripts.

Root germs require explicit conversion from the original convention. Its
`ExpCan4_Reg.wl` lines 361–378 and 608–611 flip a threshold root above the
corresponding threshold. For example, `sqrt(-b)` is opposite to the principal
square root at positive physical `b`. Nonthreshold roots retain their original
convention. A blanket principal-root seed is incorrect. The supplied canonical
boundary also carries its logarithmic continuation history.

The independent dimensionless invariants are `s/MV²`, `t/MV²`, and `b=MH²/MV²`,
with `u=b-s-t`. The application fixes `b=13074/5399` for W and `14631/7775` for Z.
The four permutations are `(s,u,b)`, `(s,t,b)`, `(u,t,b)`, and `(t,s,b)`.
`expew.wls` lines 46–58 and 87–126 apply `MV^-4` and `-1/(4*pi)^4` to the
form factors. The effective UFO couplings and Lorentz tensors are applied after
this normalization; phases and prefactors must each appear exactly once.

The selected W source has inherited integration-error estimates of
`1.6272e-27` and `3.5312e-28`. Its long mantissas and Mathematica component
Accuracy values above 238 digits do not replace those estimates. A conservative
`1e-26` per-component envelope and a 24-digit source cap are appropriate for
the first regular-point check. Final form-factor errors must sum the transported
component errors multiplied by their exact coefficient magnitudes. At the
selected point, a deliberately conservative bound formed from fully expanded
terms has coefficient norms between `5.62e3` and `5.66e4`; cancellations and
other points can change the usable relative accuracy. This diagnostic assumes
the same `1e-26` envelope for every permutation and does not establish their
actual errors.

Native HEPKit already owns the amplitude machinery. Load one shared
`Arc<Model>` from JSON, use `ModelProcessExt` and `Process::generate_amplitude`,
then `Amplitude::squared`, `SquaredAmplitude::sum_spins` and `sum_colors` with
Idenso/Spenso contraction. This supports coherent W/Z tensor interference,
Ward checks, and spin/color averaging without new tensor algebra. Raw UFO
loading belongs to HEPKit's `UfoLoader` and may export JSON beforehand; the
RustFlow native library remains Python-free. The plugin's effective `gggH`
vertices take already evaluated form factors. Generating them does not derive
the physical two-loop numerator or its IBP reduction.

The grid archives contain differential matrices, precanonical coefficient maps,
boundary files, indexes and diagnostic logs. They contain no explicit physical
propagator lists or integral power-vector definitions. Those are available from
separate primary publications: the [planar paper](https://arxiv.org/pdf/1810.05138)
gives its 9-slot topology, normalization, basis transformations and 48
precanonical routings; the [independent full-amplitude paper](https://arxiv.org/pdf/2007.09813)
gives planar/nonplanar 9-slot families, master lists and ancillary helicity
amplitudes. Its basis differs from the plugin's, so indexwise matching would be
invalid. These nine-slot families fit within RustRed’s compiled arity registry. The
[application paper](https://arxiv.org/pdf/2010.09451) supplies matrix-element
benchmarks and describes its effective-UFO construction.

The next acceptance steps are complete 48/61 transport at a non-grid point,
all four permutations for W and Z, uncertainty-aware assembly of all four form
factors, and native tensor/helicity or matrix-element comparison. Each numerical
stage needs an unchanged original oracle and independent native precision/order
refinement. Exact published-to-plugin basis matching and automatic native
two-loop graph reduction follow supplied-amplitude validation.

The plugin has no root license file, and neither grid archive includes one.
Its computational implementation remains an externally loaded oracle rather
than MIT source. Durable fixtures contain independently serialized mathematical equations,
newly computed numerical outputs and their provenance. They include no original
solver implementation or verbatim paper text. The absence of a source-code
license does not create a pending approval requirement for these scientific data.

The first complete-system milestone now passes: at the recorded non-grid W
point, all 545 coefficients of the 48/61 systems through epsilon order four
agree with the unchanged original solver. Native precision/order refinement and
binary cache restart also pass, with 21 propagated digits for a 20-digit request under the component criterion
`error <= 10^-digits * max(1, |coefficient|)`. This is an absolute error
criterion below unit magnitude, not 21 significant digits for every small
nonzero coefficient.
The [validation report](../reports/validation/2026-10-04-gg-hg-full-systems.json)
retains every native output, source accuracy policy, exact coordinates, root
germs and hashes. Smaller observed comparison residuals do not increase the
reported accuracy. Full crossed W/Z form factors remain the next stage.

The complete crossed supplied-system comparison also passes: sixteen W/Z and
permutation cases cover all 4,360 epsilon coefficients at two precision/order
settings, with exact cache restarts. All eight form factors and all 4,376
supplied precanonical diagnostic rows agree with the original mathematical
expressions. Exact reconstruction checks exclude nonlinear or misidentified
master terms, and the original full form-factor expressions are evaluated
independently of the serialized aggregation.

The [crossed validation report](../reports/validation/2026-10-04-gg-hg-crossed-form-factors.json)
records component-wise propagated errors. These support relative form-factor
digits W `[20, 19, 20, 19]` and Z `[20, 19, 19, 19]`; uniform 20 significant
form-factor digits have **not** been established. The largest relative error
estimate is `9.75e-20`, despite much smaller observed solver differences.
The original stored `delta` is an inherited sum of dropped-tail estimates,
not independent evidence that permits increasing the declared 24-digit source
cap. Stronger boundary evidence is required for that remaining precision gap.

Offline mathematical fixtures and tests are staged under
`target/gg-hg-offline`; their fresh production-library validation remains a
separate gate. The next physics check maps both masses to one exact
`MH²=1` phase-space point and contracts the W/Z contributions using native
HEPKit tensors. The sixteen cases above use separate W and Z contexts and do
not yet establish coherent interference.
