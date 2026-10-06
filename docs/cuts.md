# Cut integrals and phase-space terminals

`cuts::CutFamily` pairs an ordinary algebraic family with its cut measure. It records the cut slots, exact oriented cut momenta and per-loop causal conventions. There is no implicit conversion to an ordinary integral family: `ReductionBackend::reduce_cut` must explicitly support the measure. The RustRed adapter uses native `SectorConfig::deltas`, cut restrictions and cut-aware ordering. A required cut with nonpositive power is zero by its distributional definition, independently of scaleless-sector analysis. Raised cuts are supported by exact IBPs. Both the native factorized solver and RustRed’s sparse runtime bridge receive the explicit cut mask; neither treats a cut family as an ordinary uncut family.

The massive two-body terminal is `phase_space::PreparedTwoBodyPhaseSpace`. It accepts one integration loop, one external channel, exactly two normalized quadratic physical denominators, both cut, and no extra denominator slots. The specialized invariant and squared masses must be exact real rationals, with nonnegative masses and a timelike channel. `LoopPrescription::Insensitive` identifies the pure phase-space integration loop. This is a two-body terminal; general mixed real/virtual cut AMF recursion is not implemented by this entrypoint.

## Orientation, support and normalization

A squared denominator cannot distinguish `q` from `-q`. Each `CutDefinition::PositiveEnergy` therefore retains the exact momentum used in `delta_+(q²-m²)`. The caller additionally supplies `FutureTimelikeChannel`, an explicit declaration that the specified external total momentum is future directed. Scalar invariants alone cannot establish that declaration.

For `s=P²`, `lambda=(s-m1²-m2²)²-4*m1²*m2²`, the terminal checks the physical support exactly. Reversing the total cut energy gives zero. Below-threshold and pseudothreshold regions give zero, even where `lambda` is positive. The physical threshold `lambda=0` returns an unsupported-endpoint error because raised cuts require a separate endpoint treatment. A routing `q=c*k+external` retains its `|c|^(-D)` integration Jacobian.

Native HEPKit `Kinematics::two_body_phase_space` supplies the exact four-dimensional density per solid angle. RustFlow multiplies by `4*pi` and the dimensional angular ratio

```text
2^(8-2D) * pi^((4-D)/2) * lambda^((D-4)/2) * s^((4-D)/2)
    * Gamma(3/2) / Gamma((D-1)/2).
```

There is no initial-state flux or identical-particle factor in this integral. The result at `D=4` is `sqrt(lambda)/(8*pi*s)`. Required gamma factors and powers use Symbolica MPFR values through `Precision`; exact inputs never pass through `f64`.

A positive cut power `n` denotes `(-1)^(n-1)/(n-1)!` times the `(n-1)`th derivative of the on-shell delta. Differentiating with respect to its squared mass raises the power with coefficient `n`. Reduction coefficients are applied at nonzero epsilon before Laurent reconstruction.

For the supported all-cut two-body component with ordinary `+i0`, the independent uncut normalization check is

```text
Phi_n = 2 * (4*pi)^(-D/2) * Im(I_uncut_n).
```

This relation is tested using the existing auxiliary-mass differential-equation solver for the uncut unequal-mass bubble, rather than a second call to the phase-space formula. It must not be applied to a general cut graph as the imaginary part of the whole graph.

## Massless N-body terminal

`phase_space::PreparedMasslessPhaseSpace` and `solve_massless_phase_space` accept a complete massless final state with `N=L+1 >= 2` positive-energy cut momenta. All physical denominators must be cut and normalized, every cut mass must vanish exactly, and the oriented momenta must sum to the declared future-timelike channel or its negative. Several external coordinates and affine loop shifts are allowed. Symbolica computes the exact determinant of the independent loop-routing matrix; the measure retains `|det|^(-D)`. Reversed total energy and pinched required cuts give exact zero.

Let `a=(D-2)/2` and let `C2` be the existing HEPKit-normalized massless two-body volume at `s=1`. Repeated phase-space factorization with measure `dt/(2*pi)` gives

```text
Phi_N(s) = C2^(N-1)/(2*pi)^(N-2)
    * Gamma(a)*Gamma(2*a)^(N-1)/(Gamma((N-1)*a)*Gamma(N*a))
    * s^((N-1)*a-1) * |det|^(-D).
```

The convolution converges for `Re(a)>0`; this gamma expression supplies its meromorphic dimensional continuation. The N=2 gamma quotient cancels exactly, retaining finite cases such as D=2. Other exceptional gamma samples that cannot be evaluated finitely return a typed numerical error. No initial flux or identical-particle factor is inserted. The native three-body Dalitz density, integrated over the exact massless triangle, independently checks the four-dimensional normalization.

The unit volume is certified by this geometry. Raised cuts and polynomial numerators still require the caller's cut-aware reduction backend, and every nonzero residual must reduce onto that unit volume. Missing reductions and residuals outside this class are errors. Reduction nonzero conditions are checked at every exact epsilon sample. An uncut denominator in an ISP slot is invalid input.

Both terminal APIs share one sample-evaluation and independent-fit owner. The new massless fitter uses a conservative leading pole allowance `-2*L`, extended by any additional exact epsilon poles of the reduction weights; it never assumes all raised massless cuts are finite. `PreparedMasslessPhaseSpace::leading_power` exposes this bound. The existing two-body API retains its off-threshold leading-zero convention.

This class does not implement massive N-body phase space, mixed real/virtual cut recursion, or singular-threshold endpoint limits. The independent uncut check used for its all-cut `N=L+1` components is `2*(-1)^N*(4*pi)^(-L*D/2)*Im(I_uncut)` with the ordinary `+i0` prescription. It is not a prescription for the imaginary part of an arbitrary cut graph.

## Mixed cut auxiliary-mass flow

`cut_flow::PreparedCutFlow` prepares a cut-aware differential system and supplies
its infinity boundaries automatically. `evaluate`, `evaluate_samples` and `solve`
reuse ordinary finite-epsilon transport, endpoint projection and independently
refined Laurent fitting. Pass the same `CutFamily`, declared future channel,
exact kinematics and cut-aware reduction backend used by the terminal APIs.

The current domain contains one complete positive-energy final state in the
prescription-insensitive loop subspace, with all virtual directions using
`PlusI0`. Inputs specialize to exact real coefficients and nonnegative squared
masses. Every uncut physical denominator receives the auxiliary mass; `Auto` and
`All` select this placement, while an explicit list must contain the same slots.
Cuts and ISP slots are never deformed. Partial placements, algebraic distribution
cuts and opposite virtual causal signs return typed unsupported errors.

The boundary construction uses compact physical support. In the rest frame of
the future total momentum, all on-shell final-state energies are bounded and
nonnegative. No cut momentum can become hard. With all uncut lines deformed,
each remaining soft virtual integration appears only polynomially at every
asymptotic order and is scaleless. The surviving region has all virtual
directions hard and all real directions soft. Existing region expansion and
HEPKit's native tensor projector factor it into ordinary vacuum integrals and
unchanged phase-space shells. Ordinary vacuum factors recursively use the
existing AMF/FT boundary engine. Current phase-space leaves are the massive
two-body and massless N-body terminals above. This proof does not cover a
partial auxiliary-mass placement, an incomplete final-state cut, or conflicting
causal assignments.

Each cut carries its original oriented momentum and unscaled denominator through
boundary conversion. If native normalization produces `d = c D`, the adapter
restores `D` and multiplies the coefficient by `c^(-n)` for cut power `n`, requiring
exact positive rational `c`. It proves the cut shells independent before partial
fractions and interprets numerator cancellation of a required cut as its exact
cut-zero identity. Unknown or ambiguous surviving cut identities and unresolved
reductions remain errors. Routing Jacobians are retained by the exact hard
vacuum and soft phase-space factors. No imaginary part of the mixed graph is
used to obtain its cut value.

An uncut line depending only on real momenta needs a separate support check:
`Insensitive` does not supply a causal prescription for an interior pole. For a
single future external channel, the adapter rewrites that denominator in the
final-state scalar products, then uses `q_i.q_j >= 0` and
`sum_(i<j) q_i.q_j = (s - sum_i m_i²)/2` to bound it on the full physical region.
Bounds containing zero are rejected, including singular endpoints. The bound
can be conservative; rejection does not prove that a physical singularity exists.

Symbolic system caches use the existing source-sensitive keys and the explicit
cut backend identity, including oriented momenta and loop prescriptions. Physical
support and domain checks run again on preparation. This adapter does not add a
persistent numerical cut-value cache or turn working precision into an achieved
accuracy claim. General partial cut recursion, massive N-body leaves, internal
real-phase-space poles and mixed virtual signs remain part of the broader parity
goal.

The automatic connected benchmark has cut lines `r²-1`, `(P-r)²-4`, virtual
lines `k²`, `(k-r)²`, and `P²=25`. Its independent result is
`Phi_2^D(25;1,4) Gamma(eps) Gamma(1-eps)^2/Gamma(2-2eps) exp(i*pi*eps)`.
Raising the first cut differentiates both the volume and its attached bubble,
giving the relative factor `-7*(1-2eps)/96-eps`. Tests also cover the second
raised cut, virtual pinches, rank-two scalar numerators, positive-energy reversal,
nonzero reduction conditions and a triangular routing with determinant six.

The second class has massless three-body cuts and an uncut denominator
`(q1+q2)²-M²`, with `M²>s`. Its result is
`-Phi_3^D(s)/M² * 2F1(1,a;3*a;s/M²)`, where `a=D/2-1`.
Tests compare a convergent independent beta-integral series, raised uncut powers
and the 20-digit independently refined finite coefficient
`-(1-log(2))/(128*pi^3)` at `s=1`, `M²=2`.

`scripts/upstream_cut_oracle.py` optionally regenerates both examples using the
unchanged pinned AMFlow 2.0 computational sources and Kira. The recorded results
and normalization conversion are in
[`reports/validation/2026-10-06-mixed-cut-flow`](../reports/validation/2026-10-06-mixed-cut-flow).
These examples exercise the generic adapter; their analytic formulas occur only
in tests and the reference driver.

## Native HEPKit graph input

Use the application's existing `Arc<Model>` and a finalized native `FeynmanDiagram` containing `DiagramCut` metadata. No separate model registry or DOT cut parser is introduced.

```rust,ignore
let graph = GraphIntegral::new(diagram, &kinematics)?;
let (family, weights) = graph.cut_integral_group(
    0,                         // selected native DiagramCut
    &point,
    epsilon,
    4,
    vec![LoopPrescription::Insensitive],
    &context,
)?;
let targets = weights.keys().cloned().collect::<Vec<_>>();
let prepared = PreparedTwoBodyPhaseSpace::new(
    &family,
    &FutureTimelikeChannel { external: vec![Rational::from(1)] },
    &targets,
    &KinematicPoint::default(), // the graph adapter already specialized the point
    &backend,
    &options,
    &context,
)?;
let samples = prepared.evaluate_samples(&epsilon_samples, &options, &context)?;
```

The returned `weights` include the contracted numerator, projector and native overall factors exactly once. Apply them to each target's value at the corresponding exact epsilon before fitting a weighted answer. An empty weight map is an exact zero; it does not require a numerical terminal. `PreparedTwoBodyPhaseSpace::solve` performs independent-grid and increased-precision/order Laurent validation for its individual integral targets.

The adapter selects exactly one cut by index. Native `DiagramCut::cut` contains the left half-edge of each crossing. A source half-edge gives the native routed momentum; a target half-edge gives its negative. HEPKit's `LoopMomentumBasis` supplies the routing, including the original external-coordinate labels and dependent-momentum elimination. `CutFamily::new` checks that every resulting oriented momentum squared agrees with its retained normalized denominator.

Cut conversion keeps every original denominator slot and preserves signed powers. It bypasses ordinary sector deletion and partial fractions. Dependent denominators are rejected until a cut-aware decomposition can retain their measure identity. Numerator terms that cancel a required cut are removed by the exact cut-zero identity. Other cuts in the graph's inventory are not combined with the selected cut.

Ordinary `GraphIntegral::integral_groups` continues to reject cut metadata. Native compact DOT `is_cut` tags apply to matched dangling initial-state legs that HEPKit sews; they are not arbitrary internal-edge Boolean attributes. The graph-cut tests construct their cuts through native builder/partition APIs.

## Ownership and current checks

| Operation | Reused owner |
|---|---|
| Cut partitions, half-edge identity and momentum routing | HEPKit `DiagramCut`, `LoopMomentumBasis`, Linnet |
| Model denominator formulas and numerator contraction | HEPKit, Spenso, Idenso |
| Four-dimensional two-body measure and three-body normalization check | HEPKit `Kinematics::two_body_phase_space`, `three_body_phase_space` |
| Exact N-body affine routing determinant | Symbolica `Matrix<Q>::det` |
| Cut IBPs and sector restrictions | RustRed native solver and `Restrictions` |
| Arbitrary-precision algebra, powers and gamma arithmetic | Symbolica and the existing `Precision` wrapper |
| Epsilon fitting and independent refinement | Existing RustFlow fitter |

The focused tests cover native cut IBPs, missing cuts, cache restart and measure identity; the unequal-mass volume and raised cuts at `s=25, m1²=1, m2²=4`; the finite value `sqrt(6)/(25*pi)` and mass-derivative factors `-7/96` and `-11/192`; independent uncut discontinuities at epsilon `1/13` and `1/17`; threshold orientation, routing Jacobians and precision refinement. Native graph tests cover the future-channel sum, numerator cancellation, reversed partitions, edge reversal and the massless volume.

Both reduction adapters follow RustRed’s compiled runtime registry, defaulting to 1–16 scalar-product slots. RustFlow delegates to the upstream dispatcher and preserves local development in that checkout. Finite search and representation limits remain. General cut boundary recursion, cut-preserving partial fractions, mixed-sheet physical normalization, nonzero widths, algebraic/complex phase-space masses, and singular threshold endpoint evaluation remain outside this numerical terminal. Algebraic `Distribution` cuts can be reduced, but do not automatically define a positive-energy phase-space volume.
