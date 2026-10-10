# Moving cut shells with fixed occupied spatial support

This is an exploratory deformation assessment, not production admission or
three-/four-loop acceptance. The existing fixed-shell evaluator is unchanged.
The user's proposed simple cut is

\[
 C_\eta(q)=\theta(q^0)\theta(B-\mathbf q^2)
                 \delta(q^2-a-\eta),\qquad
 a=m^2,\quad B=\mu^2-a.
\]

For \(\mu\geq0\) the spatial theta agrees with the stated
\(\theta(\mu-\sqrt{\mathbf q^2+a})\). Initially take \(a>0,B>0\)
and real \(\eta\geq0\). After integrating the positive-energy shell,

\[
 E_\eta(r)=\sqrt{r^2+a+\eta},\qquad
 dW_\eta=\frac{\theta(B-r^2)}{2E_\eta(r)}\,d^d\mathbf q,
 \qquad d=D-1.
\]

The occupied ball does remain fixed. Keeping \(\theta(\mu-q^0)\) instead
would produce the different ball \(r^2<B-\eta\), which eventually vanishes.
Here the large-mass boundary is nonzero and decreases as a power of eta.

The regulated thermal definition and its physical endpoint remain those in
[the cutting-rules paper](https://arxiv.org/html/1609.04339v2). The deformed
interpolation is an additional construction; neither the physical cutting rules
nor [distributional finite-density IBPs](https://arxiv.org/html/2304.05427v2)
alone prove its generic complex continuation. The calculations below establish
specific small examples and identify the extra boundary work.

## Surface terms and physical raised targets

Write \(f=q^2-a-\eta\), \(h=B-\mathbf q^2\), and use the existing native
normalizations \(C_n(f)=(-1)^{n-1}\delta^{(n-1)}(f)/(n-1)!\),
\(H_0(h)=\theta(h)\), \(H_j(h)=C_j(h)\) for \(j\geq1\).
Eta changes only the shell, so \(\partial_\eta C_n=nC_{n+1}\).
Physical mass differentiation also changes the occupied radius:

\[
 \partial_a[C_1(f)H_0(h)]
   = C_2(f)H_0(h)-C_1(f)H_1(h).
\]

More generally its normalized \((n-1)\)-th mass derivative is

\[
 \Gamma_n(f,h)=C_n(f)H_0(h)
       -\sum_{j=1}^{n-1}\frac1j C_{n-j}(f)H_j(h).
\]

These are identities before energy integration. The usual Euclidean raising
operator \((-\partial_a)^{n-1}/(n-1)!\) and the existing simple-cut minus
sign must still be applied to the complete cut amplitude. In particular, a bare
\(C_nH_0\) is not the physical raised line for \(n>1\).
For simple cuts one can equivalently write the upper energy endpoint as
\(\sqrt{\mu^2+\eta}\), but differentiating this representation in eta
produces its moving-energy surface term. These are alternative representations
of the same effect; adding both would double count it.

For a spatial vector field \(v\), differentiation of the compactly supported
integrand retains

\[
 0=\int d^d\mathbf q\left\{
 \frac{\theta(h)}{2E_\eta}
 \left[(\nabla\!\cdot v)F+v\!\cdot\nabla F
             -\frac{v\!\cdot\mathbf q}{E_\eta^2}F\right]
 -\frac{\delta(h)}{E_\eta}(v\!\cdot\mathbf q)F\right\}.
\]

Lower positive-energy contacts also belong to the original distribution system.
They vanish in the massive example by disjoint support, not by dropping the
lower theta. The massless origin requires its separate dimensional argument.

## Explicit one-loop connection and boundary

Omit the common angular factor \(\Omega_{d-1}/(2\pi)^d\), put
\(A=a+\eta\), and define

\[
 I=\frac12\int_0^{\sqrt B}\frac{r^{d-1}\,dr}{\sqrt{A+r^2}},\qquad
 S=\frac{B^{(d-2)/2}}{4\sqrt{A+B}}.
\]

Here \(S\) is the radial integral with the upper \(H_1(B-r^2)\),
including its Jacobian. Spatial radial IBP and direct surface differentiation
give a closed two-state connection,

\[
 \frac{d}{d\eta}\binom I S=
 \begin{pmatrix}(d-1)/(2A)&-B/A\\0&-1/[2(A+B)]\end{pmatrix}
 \binom I S.
\]

This displayed derivation is independent of native rule discovery. The native
guarded-source experiment is a separate validation artifact and must pass its
own closure and replay checks before being used as evidence of automatic closure.

For \(|B/A|<1\) a uniformly convergent expansion supplies the boundary:

\[
 I=\frac{B^{d/2}}{2\sqrt A}
   \sum_{k\geq0}\frac{(-1)^k\binom{2k}{k}}{4^k(d+2k)}
       \left(\frac BA\right)^k,
 \qquad
 S=\frac{B^{(d-2)/2}}{4\sqrt A}
   \sum_{k\geq0}\frac{(-1)^k\binom{2k}{k}}{4^k}
       \left(\frac BA\right)^k.
\]

These coefficients are ordinary compact polynomial moments. On the positive
real axis the decreasing alternating tail bounds the omitted terms. No endpoint
value or numerical oracle is needed to initialize transport. With the positive
square root at real eta, the shell-integrated functions are holomorphic on
\(\mathbb C\setminus(-\infty,-a]\). The delta and theta symbols themselves
are not evaluated at complex arguments; this analytic continuation defines the
complex flow. The connection's finite singularities are \(-a\) and
\(-a-B=-\mu^2\). The tested positive-real and upper-half-plane paths are
homotopic within that branch domain.

The physical mass derivative satisfies \(\partial_a I=I'-S\).
For the original Minkowski polynomial \(N=q^2+(q^0)^2\), define
\(K=\frac12\int r^{d-1}\sqrt{A+r^2}\,dr\). Radial IBP gives
\(K=[AI+2B(A+B)S]/(d+1)\). Holding the physical-point converted
coefficients fixed along auxiliary flow gives the extension \(aI+K\),
not \(AI+K\). Both recover the original numerator at eta zero. Physical
mass differentiation of that endpoint, at fixed original numerator and mu, is

\[
 \left.\partial_a(aI+K)\right|_{\eta=0}
   =\tfrac32 I+aI'-(a+\mu^2)S.
\]

The \(I\) term from differentiating the original input coefficient must be
included. Eta differentiation alone supplies neither that term nor the physical
radius derivative.

Validation uses \(a=1/4,\mu=1\), spatial dimensions \(7/5,3,9/2\),
60/90-digit independent radial quadrature, complex eta points, and separate
surface/physical-mass checks. Thirty RustFlow transport runs vary working
precision, Taylor order, boundary order, starting eta and path. Saved predictions
pass 150 independent comparisons and 135 refinement comparisons at a
\(10^{-28}\) relative criterion; the largest observed relative difference is
\(4.88\times10^{-35}\). This tests the displayed radial connection through
the existing transport owner, not an automatically prepared generic AMF family.
See [the comparison report](../reports/validation/2026-10-10-moving-shell-exploration/transport-comparison.json).

At \(a=0\), the direct simple endpoint first converges for
\(\operatorname{Re}d>1\); one physical mass derivative needs
\(\operatorname{Re}d>3\). In that initial strip, at mu one,
\(I(0)=1/[2(d-1)]\) and \(\partial_a I(0)=-(d-2)/[4(d-3)]\).
The independent radial check verifies these before dimension continuation.
This does not validate arbitrary higher mass derivatives or threshold contacts.

## What changes at two loops

For a single occupied line in the massive sunset, shift the cut shell and both
virtual squared masses by eta. If the physical shell is below the virtual
two-particle threshold, its gap stays positive on the positive eta ray:

\[
 (\sqrt{b+\eta}+\sqrt{c+\eta})^2-(a+\eta)>0.
\]

The gap increases with derivative at least three. With two future cut momenta,
\((q_1-q_2)^2\leq(\sqrt{a+\eta}-\sqrt{b+\eta})^2\).
Thus a strict physical transfer-mass gap
\(M>|\sqrt a-\sqrt b|\) protects the double-cut virtual denominator
after shifting it too. The positive-mass sunset used here satisfies both
conditions. These are example-specific continuation bounds, not a generic
multiloop contour theorem.

The double-cut large-eta boundary expands into polynomial spatial moments.
The mixed boundary changes more substantially: \(q^0_{\rm cut}/\sqrt\eta\to1\).
The leading hard subgraph has unit virtual masses and external invariant one,
rather than zero external momentum. For the scalar bubble its Feynman-parameter
polynomial is \(1-x(1-x)\geq3/4\); this also proves its leading coefficient
has the specified below-threshold branch. The diagnostic quadrature verifies
the expected scaling \(\eta^{(D-5)/2}\) of the complete single-cut radial
term. At \(D=12/5\), the leading hard coefficient is about
0.06497221859; setting the external momentum to zero gives about 0.05584532753,
which is incorrect. These values are independent validation data, never
production boundary inputs.

Consequently the current vacuum-only hard-boundary recursion cannot be reused
unchanged. A native ordinary AMF child with the surviving external kinematics,
and a region argument for the shifted shell energies, would be needed. The
simple one-loop success does not establish a multiloop closure improvement.
The bounded two-loop native closure experiment is still required before choosing
between this interpolation and the current fixed-shell approach.
