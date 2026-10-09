# Sufficient finite-density contour domains

These are sufficient certificates for energy and auxiliary-mass continuation.
They do not establish ultraviolet/infrared convergence, weighted reduction
closure, integrated boundaries, physical endpoint projection, or numerical
accuracy. A failed certificate does not establish a pinch: another routing or a
more precise Landau analysis may admit the same physical input.

The separate [massless graph-channel derivation](finite-density-massless-contours.md)
and [compact-distribution construction](finite-density-massless-distributions.md)
describe possible positive-eta admissions. They do not enable massless numerical
evaluation or establish its eta=0 endpoint.

## A routing-dependent heavy-uncut-edge certificate

Fix a connected occupied cut set and its exact real momentum routing. Use future
occupied momenta q_j=(E_j,**q**_j), with 0<=E_j<=mu_j, and ordinary native virtual
momenta k_l. For every uncut physical edge write

```text
p_e = sum_l A_el*k_l + sum_j B_ej*q_j,
b_e = sum_j abs(B_ej)*mu_j,
delta_e = m_e²-b_e².
```

The proposed certificate requires delta_e>0 for every uncut edge, including an
edge with A_e=0. Mass squares, coefficients and chemical-potential magnitudes
must be real, with exact comparison used when they are rational. The minimum
positive delta_e gives a uniform margin over the entire occupied support.
An absent or numerator-only line need not create a pole, but checking all
physical uncut slots is a conservative sufficient implementation.

The cut shells retain their physical masses. Deform only physical uncut
quadratics by D_e(eta)=p_e²-m_e²-eta; their Euclidean inverse factors acquire
+eta. The initial proof is for real eta>=0.

### Simultaneous virtual Wick rotation

For fixed real spatial momenta and compact energies define
x_e=sum_l A_el*k_l0, b=sum_j B_ej*E_j and
omega_e²=|A_e*k_spatial+B_e*q_spatial|²+m_e²+eta.
Then |b|<=b_e<m_e<=omega_e. Rotate all virtual energy variables together by the
native +i0 Wick homotopy k_l0=exp(i*theta)*x_l, 0<theta<=pi/2. Along it,

```text
D_e = [exp(i*theta)*x_e+b]²-omega_e²,
Im D_e = 2*sin(theta)*x_e*[cos(theta)*x_e+b].
```

If the imaginary part vanishes, either x_e=0, giving D_e=b²-omega_e²<0, or
cos(theta)*x_e+b=0, giving D_e=-sin(theta)²*x_e²-omega_e²<0. Thus no denominator
vanishes anywhere on the common energy homotopy. At its real starting boundary,
the usual +i0 prescription fixes the orientation. This is a simultaneous
continuation of the vacuum amplitude; no separate principal-value assignment is
made for an individual denominator. The globally reversed native convention
used for future occupied shells is equivalent by reversing virtual dummy
integration variables.

After the Wick rotation, the real part of every Euclidean inverse propagator is

```text
Re rho_e = |A_e*k_spatial+B_e*q_spatial|² + (A_e*k_0)²
           + m_e²-(B_e*E)² + Re eta
         >= delta_e+Re eta > 0.
```

Therefore the completed Euclidean representation continues analytically through
the right half-plane Re eta>=0, initially in an integration convergence domain.
For complex eta this statement is analytic continuation of the complete contour
from real eta, not a claim that an undeformed real Minkowski energy contour is
valid for every complex mass. A numerical transport path may detour around
apparent reduction poles while staying inside this certified half-plane. A
path leaving it needs additional continuation evidence.

### Common external-energy continuation and raised lines

A common path to the occupied external energies is q_j0=lambda*E_j,
0<=lambda<=1, with all spatial momenta fixed. At lambda=0 the external energies
are Euclidean/spacelike; the same bound holds for all lambda because
|sum B_ej*lambda*E_j|<=b_e. This provides one inherited continuation for the full
uncut amplitude to its simultaneous on-shell limits, rather than independent
prescriptions for each cut denominator.

The strict margin persists in a sufficiently small neighborhood of the external
energies and independent uncut masses. It consequently permits the derivatives
of the continued test amplitude required by finite-order raised shell
distributions. It does not settle products of distributions at a threshold
mu_j²=m_j², massless lower endpoints, or infrared divergences of the occupied
measure. Those require separate support and dimensional-limit admissions.
For smooth physical mass differentiation away from support changes, record an
open above-threshold domain or a strict below-threshold zero for each occupied
line, and differentiate before identifying line masses.

Ultraviolet decay and the vanishing of Wick-rotation arcs must first be justified
in a common convergence or analytic-regularization domain. Dimensional or index
continuation then applies to the complete identities. The positive denominator
bound alone supplies neither that convergence proof nor a weighted IBP closure
certificate.

### Scope and provenance of an implementation

A certificate should store the selected cut slots, future orientations, exact
full routing/inverse routing and determinant, A and B for each uncut edge, mu_j,
m_e², b_e and delta_e. It should identify the +i0/common-external-energy path,
fixed-shell deformation, and allowed auxiliary half-plane. Keep this data in
measure/cache identity together with any separate endpoint/convergence evidence.

The certificate depends on the chosen virtual complement: mixing a compact
momentum into a virtual coordinate changes B and can turn this sufficient test
on or off without changing the underlying amplitude. It is not a graph-rank or
name-based admission. To assemble a full finite-density integral, every
contributing cut set needs this or another valid continuation certificate, or
an independently justified exact support-zero result. The ordinary zero-cut
term remains part of the assembly.

## The less restrictive massive sunset certificate

For the selected two-loop three-edge graph with positive charged masses m1,m2
and neutral mass M, |m1-m2|<M is sufficient even when the heavy-edge test fails.
The one-cut Feynman polynomial and fully occupied denominator bounds are

```text
F(x) = x*m_uncut²+(1-x)*M²-x*(1-x)*m_cut²
     = [x*m_uncut-(1-x)*M]²
       +x*(1-x)*[(m_uncut+M)²-m_cut²] > 0,
rho_neutral >= M²-(m1-m2)² > 0.
```

All-uncut deformation adds eta to the one-cut Feynman polynomial and to the
fully occupied neutral inverse propagator. Their right-half-plane continuation
is therefore admitted by the same positive-real-part argument. The physical
point m1=m2=1/2,M=1,mu=1 is strictly inside this mass domain and away from the
occupied thresholds. The proof and its independent-mass neighborhood are also
recorded in [finite-density.md](finite-density.md).
