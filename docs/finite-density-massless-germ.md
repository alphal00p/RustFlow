# Massless singleton germs and independent virtual blocks

This extends the *proof* for a singleton occupied massless loop while keeping
all virtual physical denominators at the same positive auxiliary mass eta. It
is not an evaluated vacuum/bubble formula, and it does not cover multicut
external kinematics. The thermal order removes the energy/line regulator first
at fixed positive T and eta. No arbitrary joint eta/delta limit is asserted.

## Exact structural conditions

There is one future compact vector q with mu>0, h>=1 virtual loop vectors k,
zero physical masses, and uncut momenta `R_e k+b_e q`. Every uncut physical
factor is shifted by the same positive eta. Numerators and nonpositive ordinary
indices are finite polynomials, including polynomial medium insertions. Positive
completion indices are excluded. The exact active positive-power row matrix R
must have full virtual rank h. If it does not, a nonsingular rational loop
change separates at least one unconstrained polynomial virtual integration;
that finite label is scaleless rather than a massive-germ case.

For positive Schwinger parameters,

```
A=sum_e alpha_e R_e^T R_e,
B=sum_e alpha_e R_e^T b_e,
C=sum_e alpha_e b_e^2,
U=det A,
V=U*(C-B^T A^-1 B),
phi0=U*sum_e alpha_e.
```

By Cauchy--Binet, U is a sum of squared rational h-row minors times
positive parameter monomials. The bordered Gram determinant
`det([[A,B],[B^T,C]])` equals V and is likewise a sum of squared
(h+1)-row minors. Hence U and V have nonnegative coefficients (V may
vanish identically), as required by the positive-polynomial sector argument.
Full rank gives U>0 in the interior. Gaussian minimization gives
`0<=C-B^T A^-1 B<=C`, hence

```
0 <= V/phi0 <= Qmax := max_e b_e^2.
```

The implemented singleton certificate conservatively rejects external-only
uncut factors R_e=0. Although the displayed ratio bound would still apply to
such factors, their on-cone contour admission has not been folded into this
certificate. All active quadratics and the actual auxiliary deformation are
bound to the certificate. Rational
routing provides an exact Qmax with no numerical minimization.

## Uniform UV continuation of the external-invariant germ

Put z=q_E^2/eta and rescale the overall Schwinger parameter. The relevant
polynomial is `phi(z)=phi0+z*V`. For `|z|<1/(2 Qmax)` (any finite disk if Qmax=0),
`phi(z)/phi0` is bounded away from zero on the entire positive projective domain.
This does **not** bound U away from zero; high-D ultraviolet faces remain.

The bound implies `Newton(V) subset Newton(phi0)`: a monomial outside that
polytope would dominate phi0 along a real positive toric ray and violate the
ratio bound. Resolve the finitely many positive-polynomial sectors for U and
phi0. The pulled-back V has no smaller monomial valuation on any sector face.
Thus every leading unit of phi(z) stays uniformly nonzero on a sufficiently
small z disk, including the sector boundaries. Tensor Gaussian coefficients
introduce finitely many extra U denominators and polynomial numerators, but do
not introduce a new kinematic singular polynomial.

Parametric partial integration gives a constructive meromorphic continuation:
for each sector scaling direction with a nonintegrable degree at the desired
D-neighborhood, a finite number of degree-raising operators produces a convergent
integrand plus rational degree denominators. The same operators work throughout
the small z disk because the leading units are uniformly nonzero. The resulting
amplitude is therefore meromorphic in D and holomorphic in z near zero. This
uses the finite continuation theorem, not bare simplex dominated convergence.
See Panzer, *Feynman integrals and hyperlogarithms*, Corollary 2.2.26 and Remark
2.2.27, printed pp.36–37 (PDF pp.43–44):
[official thesis PDF](https://people.maths.ox.ac.uk/panzer/paper/phd.pdf). The same source explains why
identities proved in an open convergence domain extend under these exact
parametric continuation operators.

An equivalent concrete implementation uses finite Taylor subtraction in each
resolved sector coordinate: subtract enough Taylor coefficients at zero to
make the remainder integrable, and retain the analytically integrated Taylor
terms with their exact degree denominators. The analytic z disk and D-neighborhood
are common to this finite operation. No arbitrary renormalization scale or
finite subtraction prescription is introduced; all added-back terms preserve
the unrenormalized meromorphic integral.

## Fixed-temperature regulator removal on the unbounded energy tail

The local germ disk alone does not control the finite-temperature energy
integral, whose support is unbounded. Fix positive T and eta first. On every
resolved positive-parameter sector write

```
phi(z)=phi0*(1+z*r),  0<=r<=Qmax,  Re z>=0.
```

The leading unit of phi0 is bounded away from zero on the closed sector cube.
The Newton support inclusion makes r and the finitely many sector derivatives
needed for continuation bounded there. Since `|1+z*r|>=1`, every derivative
introduced by the finite parametric Taylor subtraction or degree-raising
operators costs at most a finite power of z. Complex powers have bounded phase
on compact D-neighborhoods away from meromorphic poles. Thus the same convergent
sector representation gives a bound of the form

```
|continued virtual amplitude and each required finite jet|
  <= C(D,eta)*(1+|z|)^M
```

for some finite M on the right half-plane. The constant may depend on the fixed
positive eta; this step does not exchange eta with regulator removal. Gaussian
tensor means and covariances add only finitely many such factors.

For the regulated null shell `z=(delta^2+2i*delta*E)/eta`, this is at most
polynomial growth in E uniformly for sufficiently small delta at fixed eta.
At fixed T, take `|delta|<pi*T/2`. Then the complex Fermi denominator is uniformly
separated from its Matsubara zeros and its modulus is bounded by a constant
times the real Fermi weight. Each of its finitely many derivatives has the same
exponential large-E falloff, with constants permitted to depend on T. The
polynomial virtual bound is therefore integrable on the full large-energy tail.
At E=0 use the separately retained high-D origin/jet inequalities.

Dominated convergence now removes the line/energy regulator at fixed positive
T and eta over the entire thermal pairing. Only afterward is the local
holomorphic germ used at q_E^2=0, followed by the real-Fermi T->0 and physical
eta endpoint in the stated high-D domain. This supplies the tail justification
for the ordered thermal prescription; arbitrary joint complex delta/T limits
remain outside the certificate.

## Exact scaling of the finite original mass/energy jets

Work directly with Gaussian tensor means and covariances. A numerator component
with c covariance pairs has an external tensor polynomial times a scalar germ
whose exact scaling is

```
eta^(h*D/2-P+c) * f_c(q_E^2/eta,D),
P=sum positive active virtual/ordinary physical powers, c>=0.
```

Any explicit eta powers from negative physical indices are nonnegative.
There are no artificial inverse q^2 factors from a singular tensor projector.
Every coefficient of the j-th invariant jet at q_E^2=0 therefore scales as
`eta^(h*D/2-P+c-j)` times a meromorphic coefficient. The coefficients need not be
evaluated to certify the exponent.

At fixed T,eta remove the energy regulator on the common positive-eta branch,
then take the independent original mass derivatives and the specified massless
origin continuation. A finite cut/upper-energy jet budget J gives sufficient
strict bounds

```
h*Re D/2-P-J > 0,
Re D>2J+2,
```

with the usual enlarged lower-contact zero-jet bound when needed. These provide
a nonempty high-D open domain away from meromorphic poles. A generic exact
witness may be chosen as D=n+1/(2h+1), with the integer n strictly above every
bound. Then kD/2 is nonintegral for 1<=k<=h, avoiding the integer-shifted
Gaussian/subgraph degree poles. A half-integer witness alone would not suffice
for arbitrary h (h=4 is a counterexample). Compact radial
pairings and positive upper surfaces are then controlled exactly as in the
one-virtual-loop proof. Every required coefficient tends to zero with eta on
that already fixed branch. Continue the resulting zero family in D only after
this endpoint construction.

This argument proves more than homogeneity: it establishes the finite analytic
invariant jets by uniform ultraviolet continuation, identifies the same
regulated thermal branch, and keeps every original mass/support derivative.
It still requires native weighted closure, integrated large-eta boundary data,
transport and the shared physical projector; no single-cut answer is injected.

## Implemented sealed structural evidence

The sealed permit in `massless_endpoint.rs` retains:

- the exact singleton future routing, active rows, full-rank witness or an
  explicit free polynomial-loop witness;
- all physical masses zero and precisely all uncut quadratics shifted;
- U, V, phi0 identities and the exact Gaussian-minimization bound Qmax;
- finite polynomial degree, cut/upper/lower jet indices, P and eta degree;
- the fixed-T regulator order, joint origin convention, and exact continuation
  theorem version/identity;
- per-finite-label proof domains, with no assertion that one finite D covers all
  integer labels in the native source program.

The current permit implements this singleton class and the independent-block
class below. Its source and boundary identities retain the common origin
prescription plus the structural theorem version. E7 satisfies these two
structural classes; numerical source closure and AMF validation remain required.
The prism and E8 nontrivial multicuts need stronger independent work.


# Two compact vectors with independent rank-one virtual blocks

## Exact reusable routing test

First solve the exact affine virtual-loop shift already used by the general
degree query, so every uncut momentum is `R_e k+b_e p`, `p=q1-q2`. Separate rows
with R_e=0; their momentum must have b_e!=0 and they are pure transfer factors.
For nonzero virtual rows, normalize each by its first nonzero entry and group
parallel rows. Require exactly h distinct groups spanning the h-dimensional
virtual space. Their representative rows form an invertible rational matrix.
The resulting virtual coordinates K_b make every uncut virtual denominator
`a_e K_b+b_e p`, with a_e!=0, and no denominator couples distinct blocks.

This is invariant under a nonsingular virtual basis change and uses no graph
name. Polynomial numerators may couple blocks and the compact vectors; they do
not invalidate the quadratic separation. Negative ordinary indices are finite
polynomial insertions. For each finite generated label, recheck the active
positive rows: a missing block leaves an unrestricted polynomial virtual
integration and the entire label is a dimensional scaleless zero.

The E7 double cut [0,4] passes with two independent virtual row classes. The
prism cut [1,5] already fails the earlier single-spacelike-direction test, so it
cannot inherit this proof. Other graphs with the same exact structural property
can pass regardless of topology name.

## Independent Gaussian parameter blocks

For each active block, normalize its positive parameters separately. With
P_b the sum of its positive denominator powers,

```
A_b=sum x_e a_e^2 >= min_e a_e^2 > 0,
Q_b=sum x_e b_e^2-(sum x_e a_e b_e)^2/A_b,
0<=Q_b<=max_e b_e^2,
F_b=eta+Q_b*h,  h=-(q1-q2)_M^2>=0.
```

Finite Gaussian tensor expansion of the original polynomial produces means
proportional to p and covariance contractions internal to each block. A term
with c_b covariance pairs carries
`F_b^(D/2-P_b+c_b)`, c_b>=0, times a meromorphic Gaussian Gamma coefficient,
bounded rational functions of the positive A_b, and external polynomial
contractions. No coefficient is evaluated. The polynomial coupling of blocks
requires a finite sum of such products, not a factorization assumption about
the numerator itself.

Each Gaussian identity is first valid in an open UV convergence domain, then
continued meromorphically. At the high-D proof point choose a generic dimension
away from all integer-shifted Gamma poles. The positive-parameter integrals
below converge after that Gaussian UV continuation; no general multi-loop
forest subtraction is needed for this special quadratic direct sum.

## Conservative uniform endpoint bounds

Let J_i be the cut-mass plus upper-energy jet budget on compact loop i, and
J=J_1+J_2. Let P0 sum the positive powers of the pure transfer factors. Require

```
Re D/2-P_b-J > 0      for every active virtual block b,
Re D-2-2J_i-P0-J > 0 for each compact radius i,
Re D/2-1-P0-J > 0    for the collinear angle.
```

All slopes are positive, so every finite index set has a nonempty common high-D
open domain. These conservative bounds do not use numerator cancellations.
They may be much stronger than the eta=0 degree envelope's optimal strip.

After at most J derivatives, every virtual F_b exponent remains strictly
positive. Its Schwinger-parameter factor is therefore bounded uniformly for
compact h and 0<eta<=1, including Q_b=0 faces. The positive index weights are
integrable and A_b is uniformly bounded below. Pure transfer derivatives are
bounded by a constant times `h^(-P0-J)`. Shell/energy derivatives contribute at
worst `E_i^(-2J_i)`. Multiplying the compact shell/angle measure gives exactly the
last two strict inequalities above. This supplies an integrable dominating
function for the eta endpoint, including all original mass and moving-upper
terms. Positive lower contacts retain their separately checked high-D zero jets.

At h>0 the virtual parameter integrals converge pointwise to their eta=0 values,
including a block with Q_b identically zero (whose positive eta power vanishes).
Dominated convergence therefore identifies the hard endpoint with the original
massless common-contour amplitude in this proof domain. Meromorphic continuation
in D is taken afterward. No zero or nonzero endpoint value is supplied to AMF.

The same bound extends to a common closed right-half-plane eta path: for u>=0,
`|eta+u| >= (|eta|+u)/sqrt(2)` when Re eta>=0, and the argument of each F_b stays
in a fixed half-plane sector. Complex powers have uniformly bounded phase on
compact D-neighborhoods. Physical masses are differentiated at fixed eta before
setting them to zero; the existing positive-eta channel proof provides that
independent mass neighborhood.

As for the singleton proof, the thermal energy regulator is removed at fixed
positive T and eta before the real Fermi T->0 limit. Arbitrary complex-Fermi
joint delta/T paths are not covered. Upper distribution derivatives transfer a
finite number of energy derivatives to the complete kernel; those are included
in J_i. The complete native sector still needs source closure, occupied region
boundaries, transport, and the shared symbolic-epsilon endpoint projector.

## Implemented block permit

The sealed block variant retains the rational virtual affine shift, independent row
representatives, each physical slot's block/transfer assignment, and the exact
A_b lower/Q_b upper bounds. For each finite label record active ranks, P_b, P0,
J_i and the above strict witness. A rank-deficient label is explicitly
scaleless; invalid completion/occupation/tail indices must fail before that zero
classification. The source-origin and boundary seed identities remain the same
ordered joint dimensional prescription, additionally bound to this structural
proof version.

Combined with the singleton massive-germ proof above, this covers the
structural endpoint admission of every occupied cut sector of the E7 input. It does not cover
the prism multicuts or the E8 common-pole complete-cut prescription. The structural
certificate supplies no numerical integral value and does not bypass native
weighted reduction, boundary integration or AMF transport.


The optional [free-virtual-loop zero domains](finite-density-free-virtual-zeros.md)
can expose the same active-rank polynomial zero to native guarded discovery.
That source optimization is false by default and requires this bound sealed
proof; it does not change the endpoint prescription or replace AMF transport.
