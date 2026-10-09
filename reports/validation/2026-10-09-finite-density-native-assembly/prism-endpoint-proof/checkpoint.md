# Coupled prism endpoint proof checkpoint

Status: mathematical validation evidence only. No production admission, no integral value, and no complete four-loop numerical result. Work paused before support-wide implementation and final proof review.

The exact routed two-cut prism diagnostics produce positive Symanzik polynomials U and V and F_eta = eta U sum(alpha) + h V, where h is the nonnegative invariant of two future null occupied momenta. The bounded-normal probe is explicitly incomplete and is not an admission certificate. The joint sector construction is a separate complete diagnostic for the displayed six-positive-denominator scalar support.

## Exact joint resolution and independent replay

`prism_joint_eta_sector_probe.py` resolves the positive polynomial F_eta over all 720 primary Hepp charts, with eta included as a sector coordinate. It produces 3008 positive-unit charts of depth at most eight. Each chart factors U and F_eta as monomials times positive units; eta_old is a monomial with nonnegative exponents. The independent audit replays every substitution, factor valuation, Jacobian and eta map. All 2068 branching nodes contain every required pivot. Integrating a unit test density with the exact transformed Jacobian gives volume 1/120 in each original Hepp chart, and six in the complete projective maximum gauge.

For virtual loop count two and scalar power sum P=6, localized eta Mellin poles have exponents [j+1+P(u-f)+D(f-3u/2)+n]/m for m>0, with n nonnegative. All D slopes are nonnegative. All 940 zero-D axes have u=f=0, j+1=0, m=1, and every original alpha valuation zero. At most one occurs in each chart. Thus they give a simple regular eta^0 term, not a negative integer term or eta^0 logarithm, in a nonresonant high-dimensional neighborhood.

The fixed original reference terms obey: C0A/B P=6,R<=4; C1A/B P=7,R<=6; CC P=7,R<=4. Gaussian numerator decomposition therefore admits the conservative monomial bound U^(1-3D/2) F_eta^(D-7), with remaining parameter factors polynomial. The independent audit finds the strict sufficient bound Re D>18 for every nonregular eta exponent. C1 already includes the single original mass derivative; no second mass-jet penalty should be charged to this bound.

## UV continuation and the local remainder

Use a smooth eta cutoff equal to one near zero and supported below one. Continue the axes with m=0 by finite Taylor subtraction before taking the high-D endpoint. These axes produce meromorphic dimension/index poles, not new eta pole slopes. Their derivatives do not lower valuations in other coordinates. This is UV meromorphic continuation, not a claim that the bare high-D integral converges.

For vertical Mellin decay, apply global eta Euler derivatives before partitioning into charts. Every nonconstant derivative of an F power contains a factor R_eta=eta U sum(alpha)/F_eta. Its monomial valuations m+u-f are nonnegative; on each regular axis it has positive valuation one. Repeated derivatives give bounded polynomials in R_eta with that factor, plus explicitly positive eta numerator terms. Thus they kill the regular Mellin pole and preserve a positive strip. Integration by parts is first justified in a convergence/large-Re(s) domain, then continued; applying it before sector partition avoids artificial interface terms. Repeated integration by parts gives arbitrary vertical decay. The converse Mellin mapping theorem then yields a constant plus O(eta^gamma) for a strictly positive gap gamma. Reference: Flajolet, Gourdon and Dumas, Theorem 4, https://specfun.inria.fr/dumas/Publications/FlGoDu95.pdf . Genuine dimension poles remain explicit.

## Fixed-term compact and thermal completion

Keep Gaussian tensor structures as polynomials in the external vectors; do not introduce inverse Gram projectors. Each scalar coefficient is homogeneous, K_eta(h)=h^lambda G(eta/h). The completed local Mellin endpoint bounds G near zero. For the opposite ratio, positivity of the Schur complement gives |Schur_12|<=sqrt(C11 C22), hence V/(U sum(alpha)) has a finite uniform bound. The existing positive massive germ argument applies in h/eta near zero after UV meromorphic continuation. Smoothness on the intermediate positive ratio interval then gives |K_eta(h)|<=C(h+eta)^Re(lambda).

For the displayed fixed original terms, the worst virtual coefficient has lambda=D-7 and the pure transfer factor has power at most two. At Re D>18 their product has a positive residual h+eta power, so it is uniformly bounded on the compact energy domain. The single shell-mass derivative costs at worst E_i^-2; its radial endpoint requires only Re D>4. Angular measures require Re D>2. Upper occupation surfaces are at strictly positive mu and create no new energy-origin singularity. These estimates support dominated compact eta limits term by term, after the virtual UV continuation.

The thermal prescription remains ordered: remove the line/energy regulator at fixed positive T and eta, within a strip that avoids complex Fermi poles; then take the real-Fermi T->0 limit, the certified endpoint limit, and meromorphically continue D. Positive channel polynomials preserve the finite-eta gap. No arbitrary joint complex regulator/T limit is claimed.

## Work remaining at pause

This is a fixed-support and fixed-original-term proof draft. Every generated finite basis, candidate, and target label must receive its own active-positive-support certificate and finite numerator/jet bound. Deleting denominator factors can change the Symanzik support and must not inherit the full six-factor proof. Rank-deficient virtual supports need the separate polynomial-scaleless proof. A reusable owner should obtain exact Symanzik data from the pinned HEPKit implementation, verify complete positive-unit sector trees with explicit budgets, and retain exact support/routing/continuation identities. The current scripts are proof diagnostics, not that owner. The complete prism common thermal continuation and support-wide proof must be reviewed before production admission. E8 dependent-cut/contact issues remain separate and unresolved by this argument.
