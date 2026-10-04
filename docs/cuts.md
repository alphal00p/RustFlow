# Cut integrals and two-body phase space

`cuts::CutFamily` pairs an ordinary algebraic family with its cut measure. It records the cut slots, exact oriented cut momenta and per-loop causal conventions. There is no implicit conversion to an ordinary integral family: `ReductionBackend::reduce_cut` must explicitly support the measure. The RustRed adapter uses native `SectorConfig::deltas`, cut restrictions and cut-aware ordering. A required cut with nonpositive power is zero by its distributional definition, independently of scaleless-sector analysis. Raised cuts are supported by exact IBPs.

The current numerical terminal is `phase_space::PreparedTwoBodyPhaseSpace`. It accepts one integration loop, one external channel, exactly two normalized quadratic physical denominators, both cut, and no extra denominator slots. The specialized invariant and squared masses must be exact real rationals, with nonnegative masses and a timelike channel. `LoopPrescription::Insensitive` identifies the pure phase-space integration loop. This is a two-body terminal; general mixed real/virtual cut AMF recursion is not implemented by this entrypoint.

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
| Four-dimensional two-body measure | HEPKit `Kinematics::two_body_phase_space` |
| Cut IBPs and sector restrictions | RustRed native solver and `Restrictions` |
| Arbitrary-precision algebra, powers and gamma arithmetic | Symbolica and the existing `Precision` wrapper |
| Epsilon fitting and independent refinement | Existing RustFlow fitter |

The focused tests cover native cut IBPs, missing cuts, cache restart and measure identity; the unequal-mass volume and raised cuts at `s=25, m1²=1, m2²=4`; the finite value `sqrt(6)/(25*pi)` and mass-derivative factors `-7/96` and `-11/192`; independent uncut discontinuities at epsilon `1/13` and `1/17`; threshold orientation, routing Jacobians and precision refinement. Native graph tests cover the future-channel sum, numerator cancellation, reversed partitions, edge reversal and the massless volume.

The native runtime limit remains 1–12 denominators and the RustRed checkout is unchanged. General cut boundary recursion, cut-preserving partial fractions, mixed-sheet physical normalization, nonzero widths, algebraic/complex phase-space masses, and singular threshold endpoint evaluation remain outside this numerical terminal. Algebraic `Distribution` cuts can be reduced, but do not automatically define a positive-energy phase-space volume.
