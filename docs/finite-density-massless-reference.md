# Independent massless two-loop validation plan

Independent derivation and validation-only generator. No answers in production admission, boundary selection,
closure or transport. The complete result uses the common thermal/massless prescription established
by the one-virtual-loop cone and pure-compact certificates; the independent
reference does not replace their native admission checks.
The regulator must be removed at fixed positive temperature (or inside an
explicit Matsubara-pole-free hierarchy) before T tends to zero.

## Exact physical input

Clone only the definition of
`examples/finite_density/massive_two_loop_sunset.json`. Set all three independent
physical mass squares to zero, retaining its exact incidence, routing and
chemical assignments. Use mu=1 initially and optionally mu=3/2 as a homogeneity
check. Keep these targets:

```
scalar: powers [1,1,1], numerator 1,
raised original numerator: powers [2,1,1], numerator g1_2+u1*u2,
numerator_convention: shifted_euclidean,
Laurent range [-2,0].
```

The neutral routing is P1-P2 and both occupied future momenta are q1,q2.
Their Wick convention is P_E0=i E, so on the two-cut sector

```
rho_neutral=h=2*r1*r2*(1-z),
N=-2*E1*E2+r1*r2*z.
```

Use the original independent mass derivative to generate the raised target;
do not freeze this on-shell N before differentiation. The input labels and
graph name must not select these formulas in production.

## Direct independent radial/angular result

Let d=D-1 and

```
A_d=area(S^(d-1))/(2*pi)^d=2^(1-d)*pi^(-d/2)/Gamma(d/2),
dW=A_d*r^(D-3)*dr/2,
t=(1-z)/2, alpha=(D-2)/2,
K_b=<t^b>=B(alpha+b,alpha)/B(alpha,alpha).
```

All normalizations below are the full unscaled Euclidean loop normalization,
with d^D P/(2*pi)^D per loop; there is no MSbar or EulerGamma prefactor.
The vacuum is scaleless. A one-cut virtual amplitude has only a null external
vector and vanishes in the prescribed dimensionally continued endpoint, once
the common one-virtual-loop proof and finite original mass/support jets apply.
That zero is not an independent-per-propagator prescription.

Thus the two-cut scalar gives

```
I_scalar = integral dW1 dW2 / h
         = A_d^2 * mu^(2D-6) / [8*(D-3)*(D-4)].
```

For the raised numerator introduce a=m1^2, E1=sqrt(r1^2+a),
h(a)=-a+2*E1*r2-2*r1*r2*z and N(a)=-2*E1*r2+r1*r2*z.
Then h_a=-1+r2/E1 and N_a=-r2/E1. Applying -d/da to the complete simple
two-cut term gives at a=0

```
bulk = dW1 dW2 * [N/(2*r1^2*h)+r2/(r1*h)+N*h_a/h^2],
upper = +A_d*mu^(D-4)/4 * integral dW2 * N/h |r1=mu.
```

The upper term is negative because N<0. With n=D-2, the three scalar factors
of the bulk and upper contribution, divided by A_d^2*mu^(2n-2)/64, are

```
2*(K_-1-2)/[n*(n-2)],
[1/(n-1)^2-1/(n*(n-2))]*(K_-2+2*K_-1),
-2*(K_-1+2)/n.
```

Using K_-1=2*(n-1)/(n-2), K_-2=4*(n-1)/(n-4), their sum simplifies to

```
I_raised_N = -A_d^2*mu^(2D-6)
  * (2*D^3-25*D^2+98*D-119)
  / [16*(D-2)*(D-3)*(D-4)*(D-6)].
```

For the two-cut calculation, all three pieces are ordinary convergent compact
integrals for real D>6 (and mu>0). This is a direct independent reference strip
for the complete raised numerator and moving surface. No virtual period or
finite-density reduction coefficient is used in these formulas. Continue the
complete expression meromorphically only after deriving it in that strip.

At D=7 (epsilon=-3/2), both targets have especially simple nonzero references:

```
I_scalar   = mu^8/(393216*pi^6),
I_raised_N = -7*mu^8/(983040*pi^6).
```

An additional sample at D=13/2 avoids a special integer dimension. The native
AMF calculation must independently construct every vacuum/cut contribution,
its physical source closure and region boundary, then use the shared endpoint
owner. These formulas may be consumed only after predictions are written.

## Laurent coefficients in the same raw normalization

Set D=4-2*epsilon. Define

```
L_mu=4-2*EulerGamma+2*log(pi)-4*log(mu).
```

The exact formulas imply

```
I_scalar = -mu^2/(64*pi^4)
             * [1/epsilon + L_mu+2 + O(epsilon)],
I_raised_N = -mu^2/(512*pi^4)
             * [1/epsilon + L_mu+14 + O(epsilon)].
```

Both epsilon^-2 coefficients vanish. In particular the raised expression's
individual bulk terms contain stronger-looking poles which cancel against
each other and the upper surface. The relation

```
I_raised_N/I_scalar
 = -(2*D^3-25*D^2+98*D-119)/[2*(D-2)*(D-6)]
 = 1/8+(3/2)*epsilon+O(epsilon^2)
```

provides an algebraic check on the finite coefficient. Keep A_d and the mu
power unexpanded in reference generation until the Laurent owner acts.

## Actual validation sequence

1. Independently audit the polynomial simplification, raw measure conversion,
   vacuum/single-cut common prescription, and fixed-original-numerator mass
   derivative including the upper surface. No numerical gate is claimed yet.
2. In a validation-only Rust test, evaluate the exact Gamma/Beta reference at
   two arbitrary precisions and verify the D=7 rational-pi special values.
   An optional direct convergent angular/radial quadrature at D=7 separately
   checks the three raised bulk/surface pieces before their sum.
3. Save native complete predictions at D=7 and D=13/2, with scalar and raised
   targets, before reading references. Refine digits, boundary series order,
   start scale independently, preserving guard precision and source policy.
   Compare vacuum, both single cuts, the double cut and total, including zeros.
4. Save five native full Laurent profiles, adding an independently changed
   epsilon grid. Compare all coefficients [-2,0], checking the epsilon^-2 zero
   with the absolute threshold and each nonzero coefficient relatively.
5. Record the source/binary/input hashes, exact regulator order, nonzero
   conditions, projected branches, source replay, actual runtime and peak RSS.

Suggested criteria match the existing massive validation: relative 1e-12 for
nonzero coefficients above 1e-20, absolute 1e-25 below it, and imaginary-zero
absolute 1e-25. Precision must be established by completed independent changes,
not inferred from working digits. This narrow test does not establish the
required four-loop massless endpoint or resolve the raised prism obstruction.

## Completed independent validation

The validation-only implementation is
[`tests/finite_density_massless_reference.rs`](../tests/finite_density_massless_reference.rs).
All four ordinary checks pass, and the explicit reference generator passes.
Its [reference report](../reports/validation/2026-10-09-finite-density-native-assembly/independent-massless-reference/README.md)
saves all raised components, fixed-dimension values, Laurent coefficients and
precision/node refinements, together with actual source and binary hashes.
The finest direct D7 quadrature differs from the exact Beta values by at most
3.995e-31 relatively. These are independent reference results; the subsequent native AMF comparisons are recorded separately and pass at
D=13/2, D=15/4 and through Laurent order zero. Native D7 remains blocked by an
ordinary-boundary indicial resonance; it does not invalidate the independent
convergent D7 quadrature. No four-loop result follows from this calculation.
