# Coupled prism virtual endpoint: exact joint eta-sector candidate proof

This is a reference-construction derivation for the five ordinary virtual
prism coefficients. It is not production admission, a full thermal cut-sum
proof, or a numerical reference. Independent audit is still required.

For the six-denominator two-loop family of the prism construction document,
write alpha=(a,b,c,d,e,f), S=sum(alpha),

```
U=(a+c+e+f)*(b+d+f)-f^2,
Phi=abc+acd+acf+bcf+abd+bcd+bde+bdf+adf,
F_eta=Phi+eta*U*S                 (external invariant h=1).
```

All coefficients are positive. Put the largest alpha equal to one and resolve
the remaining five alpha ratios in all 720 Hepp charts. Each such chart has
alpha_order[k]=product(t_1,...,t_k), and projective Jacobian
`t_1^4*t_2^3*t_3^2*t_4`. Include eta in (0,1) as a sixth variable and perform
only exact positive coordinate blowups. A blowup chooses a minimal set S of
variables meeting the support of every nonconstant monomial, splits by its
largest coordinate p, and substitutes z_i=z_p*y_i for i in S\{p}. This is an
exact cover of the cube up to measure-zero boundaries, with Jacobian
z_p^(|S|-1). At every step preserve the complete exponent vectors and the
monomial map for the original eta.

The bounded exact scripts `/tmp/prism_virtual_sector_probe.py` and
`/tmp/prism_joint_eta_sector_probe.py` give 3008 final sectors, maximum eight
blowups after the primary Hepp chart. In every sector

```
U=z^u*U0(z), F_eta=z^f*F0(z), eta_original=z^m,
```

where U0 and F0 have nonnegative coefficients and a strictly positive constant
term. They are bounded away from zero on the closed cube. The check retains
all original monomials, so positivity is exact, not a sampled rank or numerical
positivity test. Numerical exponent samples are not used.

Take the Mellin transform integral_0^1 eta^(s-1) I(eta) d eta. For the scalar
unit-power family (P=6, two loops) a sector has monomial exponents

```
E_i(s,D)=j_i+6*(u_i-f_i)+D*(f_i-3*u_i/2)+s*m_i,
```

where j includes the transformed original eta^-1 measure and all Jacobians.
Finite coordinate Taylor subtraction implements the meromorphic continuation
of this compact parameter integral; U0,F0 and their derivatives are smooth on
the closed cube at every generic fixed D. Its possible eta Mellin pole families
are therefore

```
s= -lambda_i(D)-n/m_i, n>=0, m_i>0,
lambda_i(D)=[j_i+1+6*(u_i-f_i)+D*(f_i-3*u_i/2)]/m_i.
```

Axes with m_i=0 give ordinary UV continuation poles independent of s. The
exact full-sector enumeration finds D slopes
`{0,1/2,3/4,5/6,1}` and no negative D slope. Every zero-slope axis has
u_i=f_i=0, j_i+1=0 and m_i=1, with at most one such axis in each sector.
Thus the only zero-slope leading term is the simple regular eta^0 branch;
there is no multiple eta^0 Mellin pole/logarithm from these charts.
For the scalar family every other lambda_i is positive at Re D>4.

The original five tensor/jet coefficients have P<=7 and total loop-momentum
polynomial rank at most six. A Gaussian mean contributes adj(A)*B/U; a
covariance contributes adj(A)/(lambda*U). For c covariances and remaining
means, the worst inverse U power from the Gaussian numerator is at most
U^(-6+c). The lambda integration shifts the usual U exponent by -c and the
F exponent by +c. Positive parameter numerators (including alpha_4 for the
squared D4 slot) cannot worsen a boundary exponent. A common conservative
integer bound for all five original coefficients is consequently

```
U^(1-3D/2) * F_eta^(D-7)
```

up to bounded polynomial factors and genuine meromorphic Gamma prefactors.
For lower P/rank terms the difference only adds nonnegative U/F monomial
powers. This keeps exactly the same D slopes. Checking all 3008 charts gives

```
max_{m_i>0, f_i-3u_i/2>0}
 (7*f_i-u_i-j_i-1)/(f_i-3u_i/2) = 18.
```

Re D>18 therefore puts every nonregular candidate Mellin exponent strictly
positive even for this generous tensor/jet bound. The regular axes have u=f=0,
so changing these integer U/F powers cannot create a negative regular-branch
power. Choose a generic noninteger D away from the finite relevant Gamma/UV
pole hyperplanes; do not replace genuine poles by zero.

The intended final inference is to deform the Mellin inversion contour after
finite UV continuation: only the simple s=0 physical branch survives eta->0,
and every nonanalytic branch has a strictly positive real exponent in this
high-D region. Meromorphic continuation in D then defines the original
coefficient family. The argument must explicitly justify the finite Taylor
remainder and Mellin-contour growth, and identify the s=0 residue with the
continued eta=0 original-integrand mass jet, before it is promoted from a
candidate proof. It must not be replaced by the much weaker observation that
individual coordinate axes of Phi/U have positive powers.

Even a completed virtual proof leaves the compact/thermal steps: the three
occupied-pair compact Beta formulas must share the original line/temperature
regulator with the triple-cut amplitude, all mass derivatives must hold the
original numerator fixed and retain moving upper support, and the complete
three-cut overlapping endpoint subtraction forest must be evaluated and
refined. The full prism Laurent reference remains ungenerated.
