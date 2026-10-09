# Massless compact distributions with a dimensional prescription

This is a mathematical distribution construction, not a production admission.
The required common-contour prerequisite is discussed separately in
[the massless channel certificate](finite-density-massless-contours.md).
Production still requires positive occupied shell masses. The construction
applies after angular projection of the original polynomial. Keep the
independent shell variable until every cut derivative is taken.

Let `d=D-1`, `t=|q|²`, `a=d/2+j`, and take the original real-energy monomial
`E^r t^j`, with integer `r,j>=0`. Literal Euclidean `P0^r` adds its existing
Wick factor `i^r`. Write `A_d=Omega_d/(2*pi)^d` and
`C_n(s)=(-1)^(n-1) delta^(n-1)(s)/(n-1)!`.

For sufficiently large Re(D), integrate t before the energy distribution:

```text
integral_0^infinity dt t^(a-1) C_n(E²-t)
  = (-1)^(n-1) binom(a-1,n-1) |E|^(2a-2n).
```

The spatially integrated density, including the original energy numerator, is

```text
B_n(E) = A_d/2 * (-1)^(n-1) binom(a-1,n-1)
         * E^r |E|^(2a-2n).
lambda = d+2j+r-2n; beta=lambda+1.
```

For `mu>0`, `H_lower=theta(E)` and `H_upper=theta(mu-E)`, the raw cut value is

```text
M_n = A_d/2 * (-1)^(n-1) binom(a-1,n-1) * mu^beta/beta.
```

The original Euclidean occupied correction is `I_n=(-1)^n M_n`, hence

```text
I_n = -A_d/2 * binom(d/2+j-1,n-1)
      * mu^(d+r+2j+1-2n)/(d+r+2j+1-2n).
```

This agrees independently with taking `(-d/dm²)^(n-1)/(n-1)!` of the complete
simple-line expression
`-A_d/2 integral_m^mu E^r(E²-m²)^(d/2+j-1)dE`, including moving support,
before mass removal. `Re(beta)>0` is a sufficient dominated mass-removal
condition for this monomial. The energy-variable representation already
contains the Fermi surface contributions; omitting them changes the answer.

For an upper occupation distribution `H_s(mu-E)`, `s>=1`, put `ell=s-1`.
The raw compact value is

```text
A_d/2 * (-1)^(n+s-2)
  * binom(a-1,n-1) * falling(lambda,ell)/ell!
  * mu^(lambda-ell).
```

The upper endpoint is away from the massless origin when `mu>0`. All native
normalization factors remain external: multiply these spatial-measure raw
values by `(2*pi)^d/pi^(D/2)` per occupied loop. Do not apply an extra cut sign
after the existing target Wick map and whole-amplitude measure conversion.

## Lower occupation distributions

For `H_l(E)=C_l(E)` with `l>=1`, choose the sufficient open domain
`Re(D)>2n+l` for any fixed finite n,l and nonnegative r,j. In this domain the
radial shell distribution and its needed derivatives are defined on the
half-line, `B_n(E)` is `C^(l-1)` near E=0, and all its first l-1 jets vanish.
Therefore

```text
integral dE B_n(E) C_l(E) = B_n^(l-1)(0)/(l-1)! = 0.
```

For `mu>0`, the upper theta is identically one near zero. A positive-index
upper distribution instead has support at E=mu, disjoint from E=0. Both cases
give zero. Since this finite-index product is identically zero on an open
high-Re(D) domain, its common meromorphic dimensional continuation is the zero
family. This is an origin-product prescription, **not real empty support**:
the massless cone and E=0 do meet. A measure identity must retain that distinction.
The argument is finite-index; each finite source certificate can use a common
large-Re(D) domain. It does not choose one finite D for unbounded indices.

For a mixed graph at positive auxiliary mass, this local argument can also
pair with virtual factors smooth in the compact variables near the origin.
Separate analytic powers on virtual propagators must first make UV integrations
and the finitely many required derivatives converge. A common regulated energy
contour is an additional prerequisite. Meromorphic continuation is performed
only after those identities are proved together.

For a graph satisfying the positive-eta channel certificate, the strict complete
F-polynomial bound persists in a neighborhood of a compact origin. After the
necessary virtual UV regulation, the amputated test amplitude and each finite
number of required derivatives are therefore smooth there. Angular integration
of its local Taylor series leaves nonnegative radial powers, so the same
high-Re(D) zero-jet argument applies term by term. This assertion depends on the
jointly regulated complete amplitude; it does not apply to a separately assigned
singular virtual denominator.

For positive shell mass the lower-contact support is disjoint. In the large-D
domain above, removal of this mass agrees with the zero-jet continuation: the
radial densities and the finitely many needed energy derivatives have vanishing
limits at the origin. Outside that domain the prescribed continuation of the
joint product is used, not a new fixed-D product. A future measure/cache identity
must record that dimensional prescription and its common-contour certificate.

## Required boundaries of the claim

- `mu=0` makes the upper and lower endpoints coincide and needs a separate
  prescribed product/limit. The argument above does not admit it.
- An uncancelled `beta=0` denominator is a meromorphic pole, but the binomial
  coefficient can cancel it. Simplify the complete exact ratio before sampling
  D. For example `n=2,r=1,j=0` gives `(d/2-1)/(d-2)=1/2`, hence the finite
  seed `-A_d/4` at d=2. Values at `Re(beta)<=0` are not generally convergent
  radial integrals or a finite massive limit. Continuing the complete formula
  must be an explicit dimensional prescription.
- Extra inverse-energy factors are not admitted by the positive-mass smoothness
  certificate at mass zero. They need their own origin analysis.
- A virtual pole or collinear singularity at the origin can invalidate the
  smooth-test-function argument. Positive auxiliary mass and a justified
  common contour must precede the local product proof.
- Multiplying separately renormalized distributions at a fixed singular D can
  introduce other contact schemes. The zero above is fixed by continuation of
  the jointly defined high-Re(D) product, not by such an independent choice.
- None of these local identities justify the eta->0 physical endpoint of an
  entire graph. Virtual soft or collinear regions require a compatible common
  regulator and an independent endpoint argument before Laurent reconstruction.

An essential algebraic regression is `g*C2=C1` at zero shell mass, with
`g=E²-t`. Taking the two monomials through the formula yields this identity
exactly. Replacing g by its shell value zero before differentiating incorrectly
loses C1. This illustrates why the original numerator must be retained.

## Isolated native moment API

`CompactShell::dimensionally_continued_massless_raised_moment(d,power,r,j,p)`
evaluates the complete Euclidean occupied formula above using native precision.
It requires exactly zero shell mass, positive chemical potential, d>0 and
1<=power<=32. The binomial/beta ratio is cancelled exactly before assigning d,
so a removable zero is distinguished from a genuine meromorphic pole. The real
energy convention remains `E^r=(P0/i)^r`; a literal Euclidean `P0^r` requires its
existing Wick factor.

The helper also permits signed integer r,j as a jointly dimensionally continued
radial moment family. For negative powers, start sufficiently far in the
large-Re(D) convergence domain for those fixed indices before continuing the
complete formula. This is not an authorization to multiply separately continued
inverse-energy and shell distributions at a singular fixed dimension. In
particular, the nonnegative-power lower-contact proof above and the numerical
graph admission are not widened by this isolated scalar-moment API. Its values
carry neither a common graph contour certificate nor an eta=0 endpoint proof.
