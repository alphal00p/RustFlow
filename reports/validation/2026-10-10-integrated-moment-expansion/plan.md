# Integrated large-eta coefficient exploration

Status: isolated design/census; no weighted closure, production replacement, or numerical claim.

The input is `examples/finite_density/massless_three_loop_chain.json`, cut `[0,3]`, both original targets. Fixed-shell all-uncut placement is `[1,2,4]`. Retain exact prepared inverse-propagator indices and coefficients before deformation. Independent physical mass coefficient derivatives precede physical mass assignments and eta; do not silently keep a reconstructed momentum numerator fixed when shifted physical factors change it.

## Source owners and coefficient API

1. `DensityInput::prepare`, `occupied_cut`, `at_physical_masses`, `region_family` provide exact target sums, routing, factor maps and roles. Save both independent-mass and physical coefficient maps. Raw raised C_n/H indices retain all physical surface terms through the distribution moment; do not differentiate the on-shell numerator a second time.
2. Use unchanged `cut_regions::enumerate_compact` through a source-bound isolated copy because it is crate-private. Save every region, exact determinant and classification; never select only the expected all-hard region by topology.
3. `regions::expand_region` generates symbolic eta exponents and t=eta^(-1/2) coefficients. `integrand::projected_factor_region` uses native HEPKit tensor projection and hard/soft separation. Preserve residual denominators. A virtual-soft polynomial requires ordinary native zero evidence and never erases compact moments.
4. Implement a reusable symbolic compact adapter following `IntegratedOccupiedBoundary::integrate_polynomial`: Minkowski Gram-to-energy/spatial conversion, rotational projection in d=D-1, then independent radial/distribution moments. Return exact rational coefficients times explicit n-independent angular/Gamma normalization and chemical-potential powers. Retain raw C_n signs, upper derivatives, joint high-D lower contacts and denominator conditions.
5. Hard factors use `integrand::to_integrals` then ordinary RustRed only, bubble_subloops=false. Keep independent vacuum-master channels formal. No reference values or special bubble formula.
6. Multiply exact target coefficients and Jacobians before combining channels. Preserve eta exponent offsets, logarithms and odd half-grades; collapse to integer grades only after a structural parity proof and native checks.

Schema: sum over channels of `eta^alpha(epsilon) log(eta)^ell M_channel(D,mu) sum_n c_n(epsilon) t^n`. Each M binds the full ordinary family/master plus compact normalization. The first fitting field is exact D=13/2; symbolic Q(epsilon) coefficients remain preferable when affordable. Known prefix means every included grade is present; omitted later grades are unknown. Different hard periods are not presumed independent; coefficientwise annihilation is sufficient.

## Bounds and material limitations

Census and initial coefficient probes each stop at 180 seconds. Initial request <=24 coefficients, <=200,000 retained monomials and <=128 MiB retained serialized algebra. Begin with four integer orders and exact half-grade output, then measure growth. These do not bound internal CAS scratch. Only a feasible pilot justifies total coefficient work <=900 seconds and maximum64 coefficients (48 fitting+16 untouched holdout).

Existing region half-order cap is100; existing tensor rank cap is32 (general HEPKit rank20). A naive24 integer-coefficient series can reach rank46, and64 cannot silently bypass either cap. Report this limitation. Investigate a generic Gaussian/Symanzik-before-expansion route, whose tensor rank is the original numerator rank and whose parameter monomials map to ordinary vacuum integrals, potentially in shifted dimensions. Current Gaussian terminal is numerical/one-positive-denominator; a reusable parameter-monomial/dimension-shift lowering owner has not yet been located. No benchmark-specific formula is an acceptable substitute.

All contributing regions and any fractional/log branches must be explicit. At generic D, Taylor-expanded all-mass integrands generate no logs, but resonance/dimensional continuation remains separate. Save predictions before independent finite-eta comparisons. Finite recurrence fitting is not proof; require unused coefficients and precision/truncation/continuation checks. Existing occupation-aware closure remains paused.
