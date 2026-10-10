# Moving cut shells with fixed physical spatial occupations

This is a mathematical assessment of the proposed auxiliary deformation, not a production implementation or an amplitude evaluation. No oracle values or evaluated bubble formulas are used. Let `d=D-1`, `a_i=m_i^2`, `R_i^2=mu_i^2-a_i`, and initially assume `mu_i>m_i>=0`.

For each occupied line use

\[
q_i=(E_i(\eta),\mathbf r_i),\quad E_i(\eta)=\sqrt{\eta+a_i+r_i^2},\quad r_i\le R_i,
\]

with the principal square root continued from positive eta. Its simple-shell measure is `d^d r/(2 E_i(eta))`. Shift every uncut physical denominator by the same `-eta`. The spatial occupation stays fixed as eta varies. The equivalent polynomial upper factor is

\[
g_i=\mu_i^2-a_i-r_i^2=\mu_i^2-a_i-q_{i0}^2+q_i^2.
\]

It is independent of eta before restricting the shell. On the shell it equals `mu_i^2+eta-q_i0^2`; replacing it by that shell expression before taking distribution derivatives is not valid without the corresponding distribution identities.

## Small sunset continuation

For the mixed sector, let the occupied external line have physical squared mass `a`, and the two virtual lines have squared masses `b,c`. The normalized one-loop Feynman denominator is

\[
F_\eta(x)=\eta[1-x(1-x)]+xb+(1-x)c-x(1-x)a,\qquad 0\le x\le1.
\]

The eta coefficient lies in `[3/4,1]`. If `sqrt(a)<sqrt(b)+sqrt(c)`, the physical term is positive on the simplex (with the familiar endpoint qualifications if a virtual mass is zero). Equivalently the normal threshold gap

\[
(\sqrt{b+\eta}+\sqrt{c+\eta})^2-(a+\eta)
\]

is positive at eta zero and increases with derivative at least three for positive eta. Hence the small massive sunset `a=b=1/4,c=1` has no crossing on positive eta. The all-massless case has `F_eta>=3 eta/4` for eta positive. More strongly, `Re eta>0` gives a positive real parametric denominator for the massless case; with a positive physical gap this also covers the closed right half-plane. Ordinary UV poles in D remain meromorphic conditions, not contour singularities to discard.

For the double-cut sector the remaining neutral denominator is

\[
D_{12}(\eta)=(E_1(\eta)-E_2(\eta))^2-|\mathbf r_1-\mathbf r_2|^2-M^2-\eta.
\]

On the real positive path, future timelike Cauchy gives

\[
D_{12}(\eta)\le(\sqrt{a_1+\eta}-\sqrt{a_2+\eta})^2-M^2-\eta.
\]

Thus `M>|m1-m2|` is a sufficient uniform physical gap, preserved for every nonnegative eta. The small equal-mass sunset has `D12<=-1-eta`. The all-massless case has `D12<=-eta`, with the physical collinear/origin limit still requiring dimensional continuation at eta zero.

There is also a useful complex bound. Set `x_i=a_i+r_i^2`. For `Re eta>=0`, principal roots satisfy `Re sqrt(x_i+eta)>=sqrt(x_i)`. Using

\[
E_1-E_2=\frac{x_1-x_2}{E_1+E_2}
\]

shows `|E1-E2|<=|sqrt(x1)-sqrt(x2)|`, and therefore

\[
\Re D_{12}(\eta)\le D_{12}(0)-\Re\eta.
\]

The simultaneous `x1=x2=0` case has zero energy difference directly. This proves a right-half-plane continuation domain for these sunset sectors. It is not a certificate for arbitrary graph channels or arbitrary complex shell deformations.

## Large-eta boundary changes

The cut energy is hard even though its spatial momentum is bounded:

\[
E_i=\sqrt\eta\left(1+\frac{x_i}{2\eta}+O(\eta^{-2})\right),\qquad
\frac1{2E_i}=\frac1{2\sqrt\eta}\left(1-\frac{x_i}{2\eta}+O(\eta^{-2})\right).
\]

For the double cut,

\[
D_{12}=-\eta-M^2-|\mathbf r_1-\mathbf r_2|^2+
\frac{(x_1-x_2)^2}{4\eta}+O(\eta^{-2}).
\]

Consequently its expansion consists of compact polynomial spatial moments. Two simple cuts and one simple neutral propagator scale as `eta^-2` before numerator powers. This is a useful simple boundary for this sector.

For the mixed sector the virtual loop scales as `k=sqrt(eta) K`, while the occupied external momentum tends to `sqrt(eta) u`. The leading virtual integral has denominators

\[
K^2-1,\qquad (K-u)^2-1,\qquad u^2=1.
\]

It is a massive two-point period at external invariant one, below threshold four. It is not a vacuum tadpole. The uniform Feynman bound above shows that apparent soft-energy cancellations of one denominator do not create an extra pinched soft region in this one-loop case: the full parametric denominator stays proportional to eta.

For a scalar mixed sunset the overall power is `eta^(D/2-5/2)`: one occupied factor `eta^-1/2` and a virtual bubble factor `eta^(D/2-2)`. Its coefficient is the compact spatial volume times the unevaluated unit-scale two-point period. A native implementation would need a generic external-momentum ordinary boundary owner, for example a nested auxiliary flow for that period whose own large auxiliary mass has vacuum boundaries. The current vacuum-only recursive boundary cannot be reused unchanged. No evaluated analytic bubble seed is proposed.

For higher graphs, leading occupied energies can generate nontrivial hard external invariants and potentially hard thresholds. The sunset proof alone does not establish a simpler general boundary theory.

## Raised original propagators require upper-surface terms

Define the raw normalized distributions

\[
C_n(f)=\frac{(-1)^{n-1}}{(n-1)!}\delta^{(n-1)}(f),\qquad
H_0(g)=\theta(g),\quad H_j(g)=\frac{(-1)^{j-1}}{(j-1)!}\delta^{(j-1)}(g)\ (j\ge1).
\]

For one line, `f=q^2-a-eta`, `g=mu^2-a-r^2`. With independent physical line mass a, fixed mu, and the original polynomial numerator held fixed before shell evaluation,

\[
\frac1{k!}\partial_a^k[C_1(f)H_0(g)]
=C_{k+1}(f)H_0(g)-\sum_{j=1}^{k}\frac1j C_{k+1-j}(f)H_j(g).
\]

This follows from `partial_a^l C1=l! C_(l+1)` and `partial_a^j H0=-(j-1)! Hj` for positive j. The original Euclidean occupied correction at power `n=k+1` has the usual external phase `(-1)^n` multiplying this complete combination. In particular, power two requires `C2 H0-C1 H1`.

Using only `C_n(f) H_0(g)` for a raised original line gives the wrong physical eta-zero target. Although the two simple-shell occupation definitions agree on shell, their off-shell distribution derivatives do not. This is a concrete mandatory design constraint for the proposal.

An exact way to see the same identity is to set `mu_eff(eta)=sqrt(mu^2+eta)`. On the future simple shell, fixed physical spatial occupation equals `theta(mu_eff(eta)-q0)`. Differentiating that equality in the independent physical mass a, at fixed eta and mu, gives

\[
C_n(q^2-a-\eta)\theta(\mu_{\rm eff}(\eta)-q_0)
=C_n(f)H_0(g)-\sum_{j=1}^{n-1}\frac1j C_{n-j}(f)H_j(g).
\]

Thus the complete raised family can equivalently use an energy occupation with moving upper endpoint `mu_eff`. In that representation eta derivatives do produce upper-surface terms through `partial_eta mu_eff=1/(2 mu_eff)`. The fixed-spatial representation removes explicit eta dependence from H by retaining the finite initial C/H combinations instead; it does not remove these distributional effects from the mathematics. This is a zero-temperature distribution identity away from degenerate threshold; it does not identify arbitrary finite-T regulators.

Eta differentiation is simpler: g is eta independent and `partial_eta C_n=n C_(n+1)`. The initial physical target is the finite C/H combination above; a native weighted family may then differentiate each C/H term with the ordinary moving-shell rule. It must not mistake the complete physical mass derivative for a pure cut-index raise.

Upper-surface pieces can dominate large eta for raised original propagators: differentiating the fixed physical radius changes the compact volume already at order `eta^-1/2`. Thus increasing the original occupied propagator power need not improve the large-eta power by one. The boundary planner must retain the actual C/H target combination and permit cancellations.

An elementary one-line check at `D=3` makes the distinction explicit. With spatial measure `d^2r/(2pi)^2` and `mu>sqrt(a)`, the simple Euclidean occupied correction is

\[
I_{1,\mathrm{occ}}(\eta)=-\frac{\sqrt{\mu^2+\eta}-\sqrt{a+\eta}}{4\pi}.
\]

Its original raised counterpart is `-partial_a I1=-1/(8 pi sqrt(a+eta))`. By contrast the bare fixed-ball `C2 H0` seed is

\[
\frac1{8\pi}\left[\frac1{\sqrt{\mu^2+\eta}}-\frac1{\sqrt{a+\eta}}\right].
\]

Subtracting `C1 H1`, whose surface integral is `1/(8 pi sqrt(mu^2+eta))`, restores the correct mass derivative. The complete raised quantity scales as `eta^-1/2`, while the bare C2 term scales as `eta^-3/2`. This is a direct radial check of the distribution identity, not a proposed replacement for native flow evaluation.

For `N=g12+u1*u2`, the fixed Wick-mapped numerator on two future shells is `r1.r2-2 E1 E2`, so it carries a leading eta power. Mass derivatives act on its shell-evaluated energies through C_n, while explicit original numerator coefficients are held fixed. No original numerator may be replaced prematurely by a mass-shell identity before differentiating.

At `mu=m`, the radial cutoff degenerates. Its C/H products must inherit the fixed-T Fermi prescription and its limiting surface factors; assigning a bare value to theta at the endpoint is insufficient. Massless eta zero likewise retains the existing dimensional-origin prescription. Away from threshold, `mu<m` gives an empty occupied region throughout this deformation.

## Assessment

The proposed deformation has a concrete safe common path and elementary fully occupied boundary for the small sunset. It also removes eta derivatives of occupation factors, provided the upper factor is represented off shell correctly. The main cost is a new hard two-point boundary even for a mixed two-loop sunset, plus finite upper-surface target combinations for raised original lines. A useful next comparison is a native one-loop density family and a native mixed-sunset boundary period, with exact eta-zero target matching and equal accuracy/refinement criteria. This assessment does not yet prefer it as a replacement for the existing fixed-shell scheme on arbitrary graphs.
