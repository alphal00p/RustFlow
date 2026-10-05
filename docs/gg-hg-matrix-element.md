# Native HEPKit validation of the coherent gg → gH matrix element

The full native HEPKit calculation reproduces the independent MG5/ALOHA binary128 result at the common physical W/Z point. The native result is

```text
|A_W + A_Z|² = 4.330896620239807810685905609324107407099586592501332217072507265531229985795405567426467250200459763661630547134293e-8
```

The observable sums final spins/colors and averages the two incoming gluons, with `IDEN = 256`. It is the pure electroweak squared contribution, not interference with a HEFT amplitude. The [validation report](../reports/validation/2026-10-04-gg-hg-squared-matrix-element.json) records all eight complex form factors, their source uncertainties, the exact point, native outputs, oracle values, source/binary hashes, and the generated-oracle scaffold provenance.

The existing form-factor evidence supports **19 relative digits** in this matrix element. Its propagated absolute allowance is `1.03696527954088948e-27`, corresponding to `2.39434318218261064e-20` relative. The uniform 20-digit goal remains open. The much smaller native/oracle difference, `1.30970995865925013e-40` absolute or approximately `3.03e-33` relative, checks the independent tensor and arithmetic implementations; it does not strengthen the input form-factor evidence.

## Exact inputs and native owners

Use `MH² = 1`, `MW² = 5399/13074`, `MZ² = 7775/14631`, `α = 1/128`, `αs = 118/1000`, and

```text
s = 7173070292440521/111284741846000
t = -12058167788971/339319588980
u = 1 - s - t
GF = π α MZ² / (sqrt(2) MW² (MZ² - MW²)).
```

Momentum labels are `G0,G1` incoming and `G2,H3` outgoing, with `t=(P0-P2)²` and metric `+---`. Form factors already include `-MV^(-4)/(4π)^4`; native model couplings are applied once. Native model `value` fields containing binary64 defaults are never used in the evaluation. Original analytic parameter definitions and coupling definitions are expanded, and every supplied decimal form factor is parsed as an exact rational before MPFR evaluation.

The calculation uses native `Process::generate_amplitude`, `Amplitude::squared`, spin/color sums, Idenso `SymbolicTensor::expanded` and `simplify_algebra`, native SU(3) Casimir formulas, and `Kinematics::apply`. Symbolica owns scalar differentiation, polynomial conversion, and the 256/384-bit evaluations. No tensor or color contraction algorithm is reimplemented here. The native 256/384-bit change is `4.055674264672502e-84` absolute.

Native polynomial grouping certifies that the final expression is homogeneous quadratic in the 16 real form-factor components. The conditional source-error estimate uses

```text
Σ_i |∂f/∂x_i| δ_i + (1/2) Σ_ij |∂²f/∂x_i∂x_j| δ_i δ_j.
```

Each real/imaginary component conservatively receives the full corresponding complex form-factor allowance. An additional `1e-65` arithmetic reserve exceeds the observed precision-refinement change. These are propagated numerical allowances, not outward-rounded interval certificates.

## Original-UFO scope compatibility

The pinned plugin's three Lorentz tensors contain reciprocals of closed indexed momentum products. Direct scalar parsing distributes those reciprocals into inverse vectors, losing the contraction scope; native Idenso correctly rejects that open-index inverse.

The private importer preserves the original expression before normalization using Symbolica's `Token` AST. It wraps the three reciprocal momentum products in native `spenso::bracket`. Idenso then validates their scalar interfaces and contracts them to native `spenso::dot`. This adds scope metadata only: all original Lorentz expression strings were checked against the pinned UFO, and the external model and shared HEPKit sources remain unchanged. Generic regressions cover closed denominator contraction, retained free numerator ports, negative powers, nested division, idempotence, unchanged scalar expressions, and continued rejection of an open denominator.

The separately recorded Python 3/metadata compatibility fixes retain the original tensor expressions. The [optional validation workspace](../tools/hepkit-amplitude/README.md) provides runnable preparation, native import, amplitude generation, and scalar comparison commands. Its dependencies are separate from the root RustFlow library; Python is used only by the optional native UFO importer.

The [fresh workspace acceptance](../reports/validation/2026-10-05-native-amplitude-tool.json) reproduces the scalar result after preparing and importing a new private UFO copy. Its five release tests and strict all-feature Clippy check pass; a deliberately incorrect oracle value is rejected without a success report. This reruns the native contraction using the recorded form factors, not their numerical integration.

## Independent oracle

The oracle uses native generated MG5 3.7.0/ALOHA helicity, color, and tensor expressions, supplied with the same coherent form factors. The binary128 run uses native MP kernels/model functions and explicit precision-kind/scaffold changes, all recorded in the report. Four linker aliases connect two native MP naming conventions without modifying the generated tensor bodies. Direct MP COMMON input avoids a legacy parameter-card string truncation. Its result is `4.33089662023980781068590560932409431e-8`.

This validates one complete supplied-amplitude point. It does not claim automatic generation or reduction of the two-loop master families, automatic AMFlow boundary construction for the full amplitude, or a performance comparison.

## HEFT contribution and interference

At the same exact phase point, the optional native tool now also reproduces the plugin's infinite-top HEFT contribution and its interference with the two-loop W/Z amplitude:

```text
|A_HEFT|²                 = 1.7683949285604083825551314387453478289821e-3
2 Re(A_HEFT conj(A_W+A_Z)) = 1.6658019366610435271343417452342751052535e-5
```

The [HEFT validation report](../reports/validation/2026-10-05-gg-hg-heft-interference.json) records independent original MG5/ALOHA binary128 evaluations, native 256/384-bit refinement, exact scalar coefficients from the plugin's HEFT bridge, and the coupling-order selections. Relative arithmetic differences are approximately `1.08e-34` for the HEFT square and `1.58e-33` for the interference. This is the infinite-top model; finite-top QCD amplitudes have not been validated.

The native HEPKit contraction is shared with the pure EW calculation. Symbolica separates its scalar squared amplitude into the three coupling orders and checks exact reconstruction. The interference is exactly linear in the 16 real EW form-factor components, so native differentiation propagates the supplied allowances without an omitted quadratic remainder. Including the arithmetic reserve gives an absolute allowance of approximately `1.38618e-25`, or `8.32140e-21` relative: **20 conditional relative digits**. The separate pure EW square retains its 19-digit bound.

The [production-tool acceptance](../reports/validation/2026-10-05-heft-tool.json) reruns the argument-based input preparation and native contraction and checks that the refactored pure EW result is byte-for-byte unchanged. These checks use recorded numerical form factors; they do not establish automatic AMFlow initialization of the complete amplitude or full-workflow performance parity.
