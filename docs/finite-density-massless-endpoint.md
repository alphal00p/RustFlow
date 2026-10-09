# Narrow massless endpoint evidence and the next native gate

This document accompanies the implemented
[`massless_endpoint.rs`](../src/finite_density/massless_endpoint.rs). Its certificates
read only admitted input definitions, exact routed momenta, target powers and
original polynomials. It does not use any supplied answer or an evaluated virtual period.
The general degree-evidence queries retain false numerical-admission flags.
A separate sealed `MasslessFlowEvidence` covers only the two proved two-loop
classes and binds the actual boundary/source origin prescription. The native
numerical gates remain separate from these algebraic certificates. The
one-virtual-loop proof below is stronger than the general two-cut degree
evidence; the two should not be silently identified.

## One occupied null vector and one virtual loop

Let every active uncut momentum be `a_e k+b_e q`, with exact real rational
`a_e!=0`, all physical masses zero, and the same positive auxiliary squared mass
eta in every uncut physical factor. An uncut nonpositive index is a finite
polynomial insertion. At least one positive uncut power must remain. No inverse
energy completion is allowed. Set `P=sum_e max(n_e,0)` over uncut physical slots.

For the active positive powers, ordinary Feynman parameters on `sum x_e=1` give

```
A=sum a_e^2 x_e,
B=sum a_e b_e x_e,
Q=sum b_e^2 x_e-B^2/A,
sum x_e [(a_e k+b_e q)^2+eta]
 = A (k+B q/A)^2 + eta+Q q^2.
```

The exact bounds are

```
A >= min_active a_e^2 > 0,
0 <= Q <= max_active b_e^2.
```

The second bound is weighted Cauchy--Schwarz. They hold also at every simplex
face; this is why an arbitrary multi-virtual-loop graph cannot inherit this
proof. A polynomial numerator is evaluated by a finite Gaussian mean/covariance
expansion. A term with c covariance pairs has the parameter factor

```
constant(D,indices) * A^(-D/2-c) * (eta+Q q^2)^(D/2-P+c)
```

times a bounded rational function of `A,B,x`, polynomial external contractions,
and possibly nonnegative explicit powers of eta from numerator insertions. The
constant contains the Gaussian Gamma factor. This representation is first an
identity on a UV-convergent open set of complex D (one may start with independent
indices), and then its exact meromorphic continuation defines it outside that
set. It is not a claim that the original unregulated k integral converges at
large D. Gaussian Gamma poles must be retained, and isolated sampling dimensions
at those poles are not admitted by this bound.

For a real contour displacement delta, a future shell of squared mass a has
`q_E=(delta+i sqrt(r^2+a),r_vector)`, hence

```
q_E^2=delta^2+2i delta sqrt(r^2+a)-a.
```

At a=0, `Re(q_E^2)=delta^2>=0`. Therefore `F=eta+Q q_E^2` lies in a common closed
right-half-plane cone. This is the continuation of the complete virtual
amplitude, not a prescription for separate propagator poles. Both signs of delta
are covered. For `0<=E<=M`,

```
|F| <= eta+Q_max(delta^2+2|delta|M).
```

Take the independent original cut-mass derivatives before a=0. At fixed eta>0
they are valid in a one-sided mass neighborhood; its size may shrink with eta.
After setting a=0, any total jet order J lowers the F exponent by at most J.
Derivatives of shell energies or Jacobians contribute inverse powers of E of
order at most 2J relative to the massless radial shell measure. No coefficient
acquires a simplex singularity because A has the positive lower bound above.
Thus the sufficient strict inequalities are

```
Re D/2-P-J > 0,
Re D-2-2J > 0.
```

On every compact interval away from zero, the first gives uniform vanishing as
`(eta,delta)->(0,0)` of all required jets. At the compact origin the second makes
`E^(Re D-3-2J)` integrable and supplies the zero boundary jets needed for the
specified massless origin extension. For an original H_0 occupation,
`J=n_cut-1`. With a separate upper H_s, add `max(s-1,0)` before applying this bound,
because integration by parts transfers those energy derivatives to the complete
kernel. Positive lower H indices need the explicit dimensional zero-jet proof;
the conservative bound there is enlarged by their finite contact order.

The line/energy regulator must first be removed at fixed positive temperature
(and fixed eta), so that the remaining Fermi smoothing has real energy argument.
Equivalently a separately proved controlled complex-energy path would need a
uniform Matsubara-pole margin, such as |delta|/T<pi; that alternative is not
assumed here. An unrestricted joint complex-Fermi delta,T limit is not justified:
if the occupation argument is E-i*delta-mu, delta=(2n+1)*pi*T can hit a thermal
pole at E=mu. The virtual-amplitude cone bound by itself does not exclude it.
This ordered regulator removal is part of the prescription and its identity.

For the resulting real Fermi smoothing, on an interval about the positive upper
endpoint mu transfer its finitely many derivatives to the kernel; the enlarged
J bound controls them. The noncompact finite-temperature
tail is bounded by a fixed exponentially decaying function times a polynomial
for temperatures in `(0,T0]`. At E=0 the radial bound removes the integration-by-
parts boundary terms. Consequently, after that fixed-T removal of the energy regulator, the real
thermal pairing has zero joint (eta,T) limit in this high-D domain after UV
meromorphic continuation. Continuing that
zero family in D fixes the single-cut endpoint. This conclusion concerns the
same common thermal continuation, rather than a separately assigned scaleless
value or a homogeneity-only argument.

The implemented flow computes the finite-eta native AMF solution and applies the
shared physical projector. The proof authorizes which nonanalytic endpoint modes
vanish; it does not supply the answer to the differential equation. The sealed
permit binds the actual deformation, cut set and source/origin convention. A
basis label outside the finite polynomial/jet scope cannot inherit the target
proof.

## Two compact null vectors, one spacelike denominator direction

With two occupied momenta q1,q2, solve exact affine shifts of every remaining
virtual loop so that every uncut momentum is `R_e k+b_e(q1-q2)`. The solution is a
rational rank test, invariant under changes of loop routing and edge orientation.
If no virtual loops remain, it simply requires every uncut factor to be a
nonzero multiple of the transfer squared. This pure compact case needs no UV
continuation and includes the double cut of the massless two-loop sunset.

For virtual loops, every nonempty consistent intersection of affine forms
`R_e k+b_e p=0`, `p!=0`, is enumerated and closed under linear dependence. If its
normal rank is r, scalar local infrared integrability requires

```
r*Re D-2 sum_(closed intersection) n_e > 0.
```

The general degree-evidence query records half of this degree as an eta-power
margin. Polynomial tensor
numerators cannot worsen local fixed-p infrared behavior. This finite test does
not on its own prove the uniform eta/h remainder of an arbitrary UV-subtracted
multi-loop amplitude; that obligation remains explicit.

A formal Gaussian Wick expansion of the *actual original numerator* keeps
independent mean/covariance parameters, the two independent cut masses, external
energies and `h=-(q1-q2)_M^2`. Covariance parameters carry one power of h; virtual
integration supplies `h^(L_virtual*D/2-P_uncut)`. Independent parameters may retain
terms that vanish in the real period, yielding a conservative envelope; no
period cancellation can be used to improve the bound. Original cut-mass
operators are then applied, including shell Jacobians and moving upper support,
before either cut mass is set to zero.

For each resulting envelope monomial `E1^r E2^s h^lambda`, write
`h=4 E1 E2 t`, `t=(1-z)/2`. Required strict compact bounds are

```
Re D-2+r+Re lambda > 0  (unless E1 is fixed at a positive upper surface),
Re D-2+s+Re lambda > 0  (unless E2 is fixed at a positive upper surface),
(Re D-2)/2+Re lambda > 0,
Re D>2.
```

Higher upper derivatives differentiate smooth positive-radius factors and leave
these remaining origin/collinear exponents unchanged. Joint radial scaling is
controlled by the two radial bounds. Lower-contact zero extensions remain
separate from this endpoint degree calculation.

For the seven-edge supplied definition and cut [0,4], the generated worst bounds
are `2D-7>0`, `3D/2-5>0`, and `D/2-1>0`. At D=15/4 their margins are 1/2,5/8,7/8.
The q1^2 covariance term of the original numerator is retained through the raised
mass derivative; dropping it before differentiation would give incomplete
structural evidence. No bubble coefficient enters these calculations.

For the prism supplied definition and cut [1,5], the exact affine rank condition
fails: the virtual denominator momenta retain separate dependence on the two
compact vectors. This certificate rejects that case without calling it a pinch.
A broader parametric sector proof is still required.

## Implemented native integration

`MasslessFlowEvidence` has private variants for the one-virtual null cone and
pure-compact transfer classes. Its constructor derives the proof from the
prepared input and actual deformation, and consumers recheck the bound family.
Public degree-evidence data cannot construct this permit. The joint high-D
massless origin prescription is part of source and cache identity. Generated
lower-contact zero domains and every finite weighted basis index use that same
prescription; massive empty-support rules remain separate.

The integrated boundary owner accepts dimensionally continued massless
polynomial compact seeds together with recursive hard vacuum factors only under
this sealed evidence. Inverse compact energies remain excluded. The flow
performs native weighted closure, constructs occupied region boundaries,
transports at finite eta and applies the common physical endpoint projector.
The assembly includes vacuum, both single cuts and the double cut of the
massless two-loop sunset, including the original raised cut and medium numerator.
The whole-graph numerical gates are separate from the certificate and boundary
unit tests listed below. General multiloop degree evidence, including the E7
example above, still does not authorize four-loop numerical admission.


## Audit of the entire finite native label set

`MasslessFlowEvidence::validate_labels` checks each physical label in the union
of the native basis, reconstructed targets, and retained candidate keys and
values. It first rejects positive completion indices, negative occupation
indices, wrong arity and mismatched family/deformation. These checks also apply
to labels which subsequently vanish, so a physical zero cannot mask an invalid
index domain. Zero required-cut powers give the explicit required-cut zero.
Positive lower contacts use the dimensional zero-jet origin prescription.

For a singleton occupied loop, set `P=sum max(n_e,0)` over ordinary physical
slots and `J=n_cut-1+max(H_upper-1,0)`. Nonpositive ordinary indices and completion
numerators are finite polynomial Gaussian insertions. If P=0, the unconstrained
virtual polynomial integral is scaleless for every eta; this classification is
allowed for a generated label and does not replace the sector's native AMF. If
P>0, record `D/2-P-J>0` and `D-2-2J>0`.

For two compact loops, put `J_i=n_cut,i-1+max(H_upper,i-1,0)` and `J=J_1+J_2`.
A conservative label bound discards numerator improvements and requires

```
D-2-2J_i-P-J > 0  for each compact radius,
D/2-1-P-J > 0    for the collinear transfer.
```

These follow by bounding every original mass/upper-energy derivative by the
worst inverse energy power and transfer power. Negative ordinary indices supply
finite polynomial insertions and do not worsen the bound. For all-compact P=0
labels the compact moment is retained. It is never classified as a scaleless
unrestricted virtual integral.

Every inequality has a strictly positive D coefficient, so every finite label
set has a nonempty common high-D proof domain. The recorded witness is a half
integer strictly above all its bounds; its fractional part avoids the
integer-shifted Gaussian Gamma poles. It is a witness before meromorphic
continuation, not a restriction that later sample dimensions remain large.
True dimensional poles of the actual finite-eta amplitudes and native source
conditions are still checked by their numerical/condition owners.

For lower contact index ell, the sufficient zero-origin bound is enlarged to
`D>2*n_cut+ell+2*max(H_upper-1,0)`. This is interpreted pointwise for each finite
native proof or application. No assertion is made that a single finite dimension
covers the source program's entire unbounded integer domain.

The permit stores private fields, binds the complete assigned factor/role and
routing arrays, shell masses and chemical magnitudes, exact shifts, source
options and input identity, and is not deserializable. The source factory and
boundary wrapper both recheck that binding. General four-loop degree evidence
cannot be converted into this permit. The boundary moment cache distinguishes
this flow-origin prescription from the existing terminal-only origin mode.

The coordinated native build passed all nine endpoint-certificate unit tests,
including the sealed permit, genuine loop shear, E7 degree bounds, prism rank
rejection and generated-label audit. The full finite-density unit filter passed
47 tests; occupied flow-boundary and terminal suites passed six and eight tests,
respectively. The boundary regression integrates a continued raised massless
seed outside its bare origin convergence range together with a hard Gaussian
factor, and checks upper contacts and lower-origin zero jets. Exact binary
hashes, source snapshot, logs, runtimes and memory are retained in
[`massless-native-boundary-gates.json`](../reports/validation/2026-10-09-finite-density-native-assembly/massless-native-boundary-gates.json).
These 61 passing tests do not by themselves establish a complete massless
numerical flow or any of the required four-loop results.

The audit retains the compact-origin inequality for every positive required-cut
loop even when the complete label is classified zero. Thus generated P=0
virtual polynomials and a lower contact on one of two compact loops still carry
a concrete common high-D factorization witness for the other compact factors.
Only an already-vanishing required cut needs no such domain witness.
