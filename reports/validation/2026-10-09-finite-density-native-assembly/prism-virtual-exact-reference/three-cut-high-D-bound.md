# Raised prism triangle after UV meromorphic continuation

This is a validation-only sufficient endpoint bound for the three-cut,
one-virtual-loop sector. It is not a complete prism thermal prescription or a
computed Laurent reference, and it is not production admission.

Write the three future null compact vectors as q_i=(r_i,r_i n_i),
0<r_i<mu, and h_ij=2 r_i r_j(1-n_i.n_j)>=0. The triangle simplex has
A=x0+x1+x2=1 and F_eta=eta+x0*x1*h12+x0*x2*h13+x1*x2*h23. The original
supplemental numerator is linear in the remaining virtual loop, so Gaussian
integration introduces only its mean and no covariance. Its exact numerator
is N(a) from finite-density-prism-reference-construction.md; a is the independent
squared mass of the raised occupied line, differentiated at fixed original
spatial momentum and original momentum numerator.

After taking the Gaussian virtual integral as its UV-meromorphic continuation,
the base kernel is Gamma(3-D/2) F_eta^(D/2-3) times its polynomial parameter
numerator. A finite mass jet J=1 lowers the F power by at most one. On a fixed
compact domain, Re D>8 makes both F powers nonnegative and bounds them uniformly
for 0<=eta<=eta0. Every parameter factor x_i is bounded, so F=0 parameter
faces create no new inverse factors in this sufficient high-D region. Genuine
Gamma poles are excluded by a generic noninteger D witness, then retained in
the dimensional continuation. This continuation is performed before the compact
integral; there is no assertion of an unregulated loop-momentum UV convergence
strip at these values of D.

The transfer denominators outside the triangle are h12*h13*h23. One raised
mass derivative can increase one transfer power to two. Near one pairwise
angular collision, the transverse space has real dimension D-2, so the worst
factor |x|^-4 is locally integrable for Re D>6. Near a simultaneous collision
of all three directions, take two independent transverse vectors x,y. The
worst product is |x|^-2a |y|^-2b |x-y|^-2c, with a+b+c<=4, max(a,b,c)<=2.
The total scaling condition 2(Re D-2)>2(a+b+c), together with the three pairwise
conditions, again follows from Re D>6. These conditions also control the
cluster subregions; one must not check only a generic isolated pair.

For radial origins, expand the finite original mass derivative before bounding.
The compact measure is product r_i^(D-3) dr_i. Inverse transfer factors each
contribute at most one inverse power to their incident radii; the Jacobian
mass derivative contributes r_1^-2, and h_1j'= -1+r_j/r_1. Derivatives of the
Gaussian numerator and F likewise have only finite polynomial factors and at
most the explicit r_1^-1 from E_1'. A conservative product majorant is obtained
by allowing up to five inverse powers of each r_i. Re D>7 suffices for that
majorant. The moving upper support at r_1=mu>0 is a separate term with its
original negative sign and lower singular degree; it is included, not discarded.
Thus Re D>8 is a common sufficient compact domination domain for this one mass
jet, including angular clusters and all radial origins.

At fixed eta>0, the independent mass derivative is taken inside the admitted
common regulated neighborhood. The explicit derivative kernels at a=0 then
have the uniform high-D compact majorant above, so dominated convergence allows
their eta->0 limit in this sector after the stated virtual UV continuation.
This does not justify starting with a real positive-mass eta=0 parameter
integral: h_1j(a) can be negative near collinearity for arbitrarily small a>0.
It also does not establish the common finite-temperature cut identity on its
own, nor interchange arbitrary line-regulator, temperature, analytic-index
and eta limits.

The remaining prism obstacle is substantial. The two-cut sectors have two
coupled virtual directions with U=A*B-f^2 and Phi/U projective degenerations;
they are outside the currently proved independent-rank-one block class. Their
five original virtual coefficient identities can be checked by exact native
IBP replay, but that algebraic fact alone does not prove their off-shell mass
jets commute with the common thermal endpoint. After that proof, one still
needs the complete three-cut angular/radial/simplex subtraction forest and a
refined Laurent quadrature. The existing local UV/collinear subtraction formulas
remain useful for this final numerical construction.
