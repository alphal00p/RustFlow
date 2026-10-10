# What the same-series transformation establishes

The fixed transformation is R = 4, w = R/(eta+R), alpha = 1/4 and H(w) = (1-w)^alpha G(w/[R(1-w)]). Its coefficients are obtained by an exact lower-triangular rational binomial map from the saved 64 coefficients. Applying the inverse map recovers all 64 coefficients for each target exactly. No unknown tail coefficient is inferred.

The source audit has projective U = 1, virtual compact quadratic 0 <= Q <= 1, and pure compact transfer 0 <= h <= 4. In the scalar Gaussian representation, after the common hard and compact normalization is removed, the transformed kernel is

    (1-a*w)^(5/4) / (1-b*w),
    a = 1-Q/4 in [3/4,1], b = 1-h/4 in [0,1].

The original scalar weight is a positive parameter/compact measure. Its total mass equals the absolute first scalar coefficient after the common prefactor is included. These facts come from the routed source and Gaussian representation, not from the independent finite-eta reference. They establish analyticity for |w| < 1. The same open disk holds for the finite C2 jet, with extra finite denominator powers and integrable origin factors.

For any circle radius r < 1, the scalar normalized kernel obeys the absolute bound M(r)/|h0| <= (1+r)^(5/4)/(1-r). For a real sample 0 <= w < 1 it also obeys H(w)/h0 >= (1-w)^(5/4): each numerator is at least this positive quantity and the denominator is at most 1. Thus the relative truncation error of the first N exact Taylor coefficients is at most

    ((1+r)/(1-w))^(5/4) * (w/r)^N / [(1-r)*(1-w/r)],  w < r < 1.

`bound-plan.json` fixes the five radii before evaluation. `scalar-bound.py` computes rational upward enclosures of the fourth roots on a denominator 10^60 grid and checks the fourth-power inequalities exactly. Its minimum over the fixed admissible radii is a valid bound. The 64-term bounds are 3.135e-28 at eta = 8 and 2.839e-8 at eta = 2. These bound truncation of the exact mathematical series; they do not enclose rounding from decimal master/normalization evaluation. The numerical 30/50-digit check is separate.

At eta = 0, w = 1, no radius satisfies w < r < 1. Neither this Cauchy bound nor finite-prefix analyticity supplies an endpoint bound. Sixty-four terms fail the prescribed endpoint numerical criterion, and all endpoint Padé profiles fail as well. Four endpoint Padé precision checks also fail. No extra order, precision, or fitted profile was introduced after those results.

The raised target is a signed C2 bulk/surface functional. Its first coefficient can contain cancellations and is not an absolute norm. A valid raised circle bound must sum absolute finite radial, numerator-jet, kernel-jet, and surface component norms. The separate audited calculation is now bound in raised-bound-plan.json. It gives M_raised(r) <= (1987/2772) B0 + (476/1573) B1, with B0 and B1 specified in the source derivation. Exact rational radial moments and the signed leading coefficient 43/264 are independently replayed. The same fixed radius grid gives 64-term relative tail bounds 6.609e-27 at eta = 8 and 3.959e-7 at eta = 2, using T/(|S_N|-T) only when the exact rational partial sum satisfies |S_N| > T. No scalar |h0| positivity bound is reused. This again excludes physical normalization rounding and does not bound the endpoint.

The Padé diagnostics use only the same saved coefficient prefix. Their prefix residuals and denominator values are recorded; no pole-free complex path, all-order ODE, or endpoint error theorem is asserted. The method is therefore useful for the bounded fallback assessment, while the endpoint objective remains unresolved.
