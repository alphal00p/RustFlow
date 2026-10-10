# Independent three-loop chain reference

The definition is `examples/finite_density/massless_three_loop_chain.json`.
This is a new intermediate validation case; it does not replace any four-loop
acceptance target. The graph has three vertices and five edges, with two
parallel pairs and their closing edge. It has no articulation vertex and is
not a product of vacuum graphs. Its exact routings are
`[P1,P2,P3,P1-P3,P2-P3]`; only slots 0 and 3 carry the chemical potential.
The scalar target has all powers one. The second target raises physical slot 0
and keeps the original shifted-Euclidean numerator `g1_2+u1*u2` fixed.
All masses vanish and mu=1 in the fixture.

## Complete cut sum and original mass derivative

The vacuum is dimensionally scaleless. Each singleton has two virtual loops,
no external-only denominator, and one massless on-shell external momentum.
The raised singleton has virtual power sum P=4 and mass-jet budget J=1;
the other singleton has P=5, J=0. The existing UV-meromorphic singleton germ
has positive margin `D-P-J>0` for D>5; compact origins require D>4 for
one mass jet. These sectors are zero under the common regulated thermal and
dimensional continuation, not by a real principal-value convention.

The sole nonzero sector cuts slots `[0,3]`. Set
`qA=P1`, `qB=P1-P3`, `p=P3=qA-qB`, and `h=p_E^2`.
The uncut virtual bubble has loop momentum P2 and denominators
`P2^2*(P2-p)^2`; the remaining neutral propagator is `1/h`.
For reference generation only, its scalar coefficient is

```
C(D)=(4*pi)^(-D/2) * Gamma(2-D/2)*Gamma(D/2-1)^2/Gamma(D-2),
b=D/2-3, alpha=(D-2)/2, n=3*D/2-5,
Ad=area(S^(D-2))/(2*pi)^(D-1),
K_s=B(alpha+s,alpha)/B(alpha,alpha).
```

The bubble first moment is exactly P2=p/2 by the change of variables
P2 -> p-P2. The production evaluator must still use its native AMF flow;
this Gaussian/Beta representation belongs only to the independent reference.

Introduce the independent squared mass `a=m0^2` before differentiating.
With fixed original spatial vectors,

```
E1=sqrt(r1^2+a), E2=r2,
h=-a+2*(E1*E2-r1*r2*z), h_a=-1+E2/E1,
u1=i*E1, u3=i*(E1-E2),
Nbar=(h-a)/4-E1^2/2+E1*E2/2,
Nbar_a=-1+E2/(2*E1).
```

The first term of Nbar comes from g1_2 and the second pair from u1*u2.
This is the original numerator after an exact virtual first moment. Its
implicit energy dependence is differentiated along with h. At a=0 write
`t=(1-z)/2`, so `h=4*r1*r2*t` and
`Nbar=-r1^2/2+r1*r2*(1/2+t)`.

The raised occupied target is the negative mass derivative of the unraised
two-cut amplitude. With `Wi=Ad*r_i^(D-2)/(2*Ei) dr_i`, its bulk is

```
integral W1 W2 C*h^b * [Nbar/(2*E1^2)-Nbar_a-b*Nbar*h_a/h],
```

and its **positive moving upper-surface term** is

```
Ad*mu^(D-4)/4 * integral W2 C*h^b*Nbar | E1=mu, a=0.
```

The surface is kept as a separate nonzero reference component. There is no
step that differentiates a precomputed massless answer while freezing its
shell substitution.

## Independent Beta evaluation

Let `F=Ad^2*C*4^b*mu^(2*n)/4` and `T=n^2-1`. The scalar is
`F*K_b/n^2`. The raised pieces, in units of F, are

```
measure:   -K_b/(4*n^2) + (K_b/4+K_(b+1)/2)/T,
numerator:  K_b/n^2 - K_b/(2*T),
transfer:   b*K_b/(4*n^2)-b*K_b/(4*T)
           +b*K_(b-1)/(4*n^2)-b*K_(b-1)/(4*T),
upper:     -K_b/(4*n)+K_b/(4*(n+1))+K_(b+1)/(2*(n+1)).
```

Thus their sum divided by the scalar simplifies exactly to

```
R(D)=(6*D^3-67*D^2+244*D-294)/(4*(D-5)*(3*D-8)).
R(4-2*epsilon)=3/8-11*epsilon/16+7*epsilon^2/32+O(epsilon^3).
```

The unsimplified Beta pieces also define an independent check of this rational
cancellation. No physical dimension is substituted before cancellation.
For direct compact quadrature, D=13/2 and D=31/3 satisfy the strict original
jet bounds: virtual bubble jets require D/2-3>0; pure-transfer angular and
radial factors with one mass jet require D>6. Avoid the discrete virtual Gamma
poles. These dimensions are selected by the proof, before native results.

## Laurent coefficients and normalization

The raw Euclidean scalar reference is exactly

```
S(epsilon,mu)=-mu^2/(2048*pi^6*epsilon^2) * (pi/mu^2)^(3*epsilon)
 * Gamma(1+epsilon)*Gamma(1-epsilon)^3
 /[(1-2*epsilon)^2*(1-3*epsilon)^2
   *Gamma(1-2*epsilon)*Gamma(1-3*epsilon)].
```

Put `A=-mu^2/(2048*pi^6)` and
`L=10-3*EulerGamma+3*log(pi)-6*log(mu)`.
The raw scalar coefficients at orders [-3,-2,-1,0] are

```
[0, A, A*L, A*(L^2/2+13-3*pi^2/4)].
```

The raised coefficients are

```
[0, 3*A/8, A*(3*L/8-11/16),
 A*(3*(L^2/2+13-3*pi^2/4)/8-11*L/16+7/32)].
```

As an independent normalization check, multiplying by
`(4*pi)^6*(exp(EulerGamma)/pi)^(3*epsilon)` at mu=1 gives scalar
`[0,-2,-20,-126+3*pi^2/2]` and raised
`[0,-3/4,-49/8,-543/16+9*pi^2/16]`.
Native comparisons use the raw Euclidean values. No supplied oracle
coefficients or native predictions are read by this reference derivation.
Numerical generation, direct compact quadrature and native comparisons remain
separate gates whose actual results must be recorded before claiming acceptance.
