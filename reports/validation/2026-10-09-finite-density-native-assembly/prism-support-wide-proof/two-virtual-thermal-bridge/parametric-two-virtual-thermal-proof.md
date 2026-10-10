# Two virtual loops: a uniform UV and thermal continuation bound

Draft for the additive physical permit. Existing production admission is unchanged.
This is a generic rank-two virtual argument, not a prism value or topology test.

## Exact UV charts

Let the actual active virtual routing vectors be nonzero rational vectors r_j in
R^2 of full rank. The native first Symanzik polynomial is

    U = sum_{i<j} det(r_i,r_j)^2 alpha_i alpha_j.

In a maximum-parameter primary chart set alpha_0=1. Let S contain the rows not
parallel to r_0. U vanishes exactly when all alpha_j with j in S vanish. Split
that coordinate cube by a largest S coordinate t and write alpha_j=t y_j for
j in S, with the chosen y_j=1. This is a complete finite subdivision up to
measure-zero interfaces. Then U=t U0, where U0 is polynomial, smooth on the
closed chart, and bounded below by the exact positive coefficient of
alpha_0 alpha_j. There are no other UV zero sets in this chart.

Every full kinematic coefficient V_ij, W_i, and U alpha_j contains at least one
S parameter. This follows from the bordered Gram identities and can also be
checked directly on the retained exact native monomials. Consequently F/U has
smooth bounded coefficient functions on this resolved chart, with the same
property for all finite parameter derivatives. Signed W_i coefficients cause
no division problem: the exact divisibility, not their signs, is used here.

## Finite UV continuation

For a finite label let P be the sum of positive virtual powers, R the actual
polynomial momentum degree, J the total required shell/upper jet order, and
P_S the sum of positive powers in S. Gaussian tensor terms can be bounded by

    t^(A-1) f(t),   A=P_S-D/2-R,
    f(t) containing (F/U)^(lambda-j),   lambda=D-P+c, c>=0, j<=J.

All omitted powers of t are nonnegative, and U0 and its inverse have bounded
parameter derivatives. Positive parameter weights are integrable on the
remaining compact variables. Choose an integer N with Re(A)+N>0 and subtract
N Taylor coefficients:

    integral t^(A-1) f(t) dt
      = integral t^(A-1) [f(t)-sum_{k<N} f^(k)(0)t^k/k!] dt
        + sum_{k<N} f^(k)(0)/(k!(A+k)).

This is the ordinary meromorphic continuation; avoid its pole hyperplanes and
the overall Gaussian Gamma poles. If Re(lambda)-J-N>0, all required Taylor
coefficients and the remainder depend continuously on the right-half-plane
value of F/U, including its zeros. Their derivatives are bounded by positive
powers of |F/U| times smooth polynomial factors. A fixed N works throughout a
small generic complex-D neighborhood. Because P_S>=1, the conservative strict
condition

    Re(D)/2 > P+R+J

allows both inequalities for N. Thus arbitrary finite R and J do not prevent a
common high-D proof domain. This argument is specific to two virtual loops;
the arbitrary-virtual query does not by itself supply the corresponding UV
uniformity theorem for higher virtual rank.

## Common Feynman regulators and compact domination

For real spatial momenta and q_i^0=sqrt(r_i^2-i epsilon_i), epsilon_i>0, write
q_i^0=omega_i-i gamma_i. The square-root identity gives

    omega_i omega_j-gamma_i gamma_j >= r_i r_j,
    Re(2q_i.q_j) >= 2r_i r_j(1-cos theta_ij)=h_ij.

The exact native off-null decomposition therefore yields

    Re F >= eta U sum(alpha) + sum_ij V_ij h_ij.

Diagonal shell masses and internal Feynman masses contribute only imaginary
terms. The Gaussian power uses the branch analytic in this right half-plane.
The preceding UV Taylor continuation is consequently continuous uniformly as
eta and these line regulators tend to zero, on bounded compact energies. At
large energies its finite derivatives grow polynomially; fixed-temperature
Fermi factors decay exponentially. Keep the regulator inside the Fermi strip
and remove it at fixed positive T. No arbitrary joint complex-regulator/T
limit is asserted.

Pure transfers satisfy Re[eta-c(q_i-q_j)^2]>=eta+c h_ij. They may be dominated
separately by their physical h powers. For two compact momenta, sufficient
additional bounds are

    Re(D)/2-1 > P0+J,
    Re(D)-2-2J_i-P0-J > 0  for each occupied loop i,

where P0 is the actual positive transfer power sum. Together with the explicit
inverse-energy shell-jet bounds, these give an integrable compact majorant
even before any virtual h power is extracted. Higher upper occupation
contacts are handled by energy integration by parts, not a false pointwise
bound on derivatives of the Fermi function.

Positive one-sided physical mass jets are obtained at fixed eta from a small
nonnegative real mass box, as in the bridge design. Negative or arbitrary
complex physical mass directions are not included. The jet formula at zero
mass is then controlled by the above uniform right-half-plane Gaussian bound.

Full-rank V=0 supports are included: their finite-eta coefficients are nonzero,
but every continued Gaussian Taylor term tends to zero in this high-D domain.
They do not become native finite-eta source-zero sectors. Rank-deficient
virtual supports continue to use their separate unrestricted-polynomial zero
proof.

The complete native joint-eta Mellin tree is still required for the endpoint
power/log expansion and the surviving regular branch. This UV argument does
not replace that tree. It supplies the missing uniformity with the actual
thermal line regulators; all identities remain bound to the complete cut
prescription and fixed-original numerator jets.
