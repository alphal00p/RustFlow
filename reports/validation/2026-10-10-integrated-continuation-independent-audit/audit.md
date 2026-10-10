# Independent continuation audit

This is a read-only review of the frozen source-tied coefficient producer, same-series direct continuation, and conditional scalar ODE. No producer coefficients, fit, reference values or production code are changed.

## Same-series transformation and scalar bound

The exact lower-triangular composition in predict.py is correct:

    H(w)=(1-w)^alpha G(w/[R(1-w)]),
    h_k=sum_(n<=k) g_n R^-n (-1)^(k-n) binom(alpha-n,k-n).

Its inverse is the implemented R^k sum_(m<=k)h_m binom(alpha-m,k-m). The physical stripped function is R^alpha w^-alpha H(w), with alpha=1/4, R=4. This retains the original fractional eta prefactor; evaluating H(1) is an endpoint diagnostic, not a newly proven convergent Taylor sum.

For the actual scalar source, the normalized transformed kernel is Phi=(1-a w)^nu/(1-b w), nu=5/4, 0<=a,b<=1. The native ordinary ratio and the simplex/compact scalar measure are positive at D=13/2. Consequently H(w)/h0 is a positive normalized average of Phi, and the proposed bounds are sound:

    M(r)/h0 <= (1+r)^nu/(1-r),
    H(w)/h0 >= (1-w)^nu,  0<=w<r<1.

The relative Taylor tail is bounded by their ratio times (w/r)^N/(1-w/r). The exact rational upward fourth-root enclosure proposed in bound-plan.json is a suitable implementation. No admissible r exists at w=1, so this does not bound eta0. Positivity of the scalar measure does not justify using abs(h0) for the raised target.

The prospective truncations, precisions, eta points and Pade degrees are explicitly fixed. Pade residuals and denominators sampled at a few points are diagnostics, not an all-order equation or a pole-free-path certificate. The code states this limitation. Exact composition uses the same frozen64 coefficients; no coefficients are inferred from a reference.

## Physical raised target and signs

The actual converted target numerator is the sum of the two unshifted completion monomials, so its eta-deformed numerator remains the original momentum polynomial for this input. This conclusion is source-specific, not a general permission to hold every original momentum numerator fixed along eta. Its physical mass derivative was performed before the eta convention.

For the repository's raw cuts, C_n=(-1)^(n-1) delta^(n-1)/(n-1)! and q^2 C_n=C_(n-1). Thus C2 contains a full mass derivative of the shell kernel, including the energy-residue derivative and the moving upper surface. Neither q^2 C2 nor the upper surface can be dropped. The saved generator retains off-shell compact diagonal terms until this product identity is applied. Its exact native target signs already include the Wick/cut convention; the physical normalization adapter correctly adds no further target or cut sign.

An absolute raised-kernel bound must retain three bulk pieces and a surface piece. The separate raised-norm.md derives one conservative bound directly from the emitted Gaussian mean and complete original numerator. It is a possible proof input for a separately fixed bound computation, not a retrospective claim that the current raised endpoint prediction is certified.

## Scalar candidate equation and ordinary transport

The recurrence convention is sum_j p_j(n)c_(n+j)=0. Multiplication by t^(n+r), with g=t^beta F and beta=-1/4, gives

    sum_j t^(r-j) p_j(theta-beta-j) g = t^beta B(t).

The startup polynomial B=-14/117 and the saved64-coefficient conversion check are consistent. With theta_t=-eta*d/deta, the eta companion signs in transport.rs are correct, including the auxiliary state t^beta=eta^(1/4). The hard-master and compact normalizations remain factored, so the first transported component is precisely the stripped physical scalar function.

The finite transport starts in the source-proved convergent region and varies boundary truncation, precision and starting eta. Its results are explicitly candidate-internal refinements. Endpoint matching through the existing Frobenius owner preserves the same normalization. The first endpoint component is the candidate stripped scalar limit. The saved exponent set is positive except0, so the mode classification is consistent with a finite candidate limit; numerical matching and physical_limit inherit the existing tolerance policy. The extra1/16 matching control failed its step cap and must remain a failed control.

Neither finite-heldout agreement nor successful ordinary transport proves that this ODE annihilates the physical integral. The exact K^-4 pole witness correctly obstructs only the fixed direct-eta certificate ansatz with pole order<=2; it is not a proof that no identity exists. The report's candidate-only and physical_endpoint_accepted=false labels are essential.

## Spending and scope

Coefficient generation, fitting, numerical transport and mathematical checks have distinct recorded caps and timing. The reported2.0393seconds is cumulative coefficient execution, including the preserved compact-index failure; it is not total compilation or total project time. The same-series fallback uses no new coefficients or equation input. No ten-digit physical eta0 claim follows from the present finite-eta bounds, Pade profiles or conditional ODE alone. The source report's resource limitation for transient CAS scratch and simultaneous Atom heap remains explicit.

## Generic machinery and source-specific inputs

The region enumeration, Gaussian numerator map, ordinary proof-backed vacuum reduction, independent-sphere moment library, exact series composition and ordinary transport are reusable methods. The actual physical factor maps, independent mass jets, exponent 1/4, simplex mean, kernel bounds 1 and 4, native master normalization, and the raised component norm are necessarily derived from this input. Their explicit specialization is legitimate proof data; it does not authorize a topology-specific production formula or extrapolation to a different graph. Higher virtual tensor rank and several ordinary hard-master channels still require their corresponding generic owners.

The separately derived raised absolute norm is validated algebraically by exact rational checks and independent review. Evaluating it on the already fixed bound grid is a separate report operation; this audit neither introduces another profile nor reads reference values to set its constants.
