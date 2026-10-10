# Independent connected hard-boundary check

The ordinary native integral has denominators K²-1, L²-1 and
(K-L)²-1, each to power one, and measures d^D K/(i pi^(D/2)) and
d^D L/(i pi^(D/2)). The dimension is D=12/5 (epsilon=4/5).
All one-loop ultraviolet subgraphs require D<4 and the whole graph requires
D<3; positive unit masses remove infrared divergences. Thus the real
Euclidean Schwinger representation is absolutely convergent here.

After integrating both Gaussian momenta and the overall Schwinger scale,
the positive Euclidean pi-normalized value is

```
Gamma(3-D) integral_{a+b+c=1; a,b,c>=0} (ab+ac+bc)^(-D/2) da db.
```

The native Minkowski value is its negative: each of the three denominators
contributes a minus sign under Wick rotation, while the two specified native
loop measures remove the two rotation factors. No (4 pi)^(-D) conversion is
applied to this already pi-normalized ordinary coefficient.

Partition the simplex into its three maximum-parameter sectors. In the
sector a is maximal, put (a,b,c)=(1,x,y)/(1+x+y), with x,y in [0,1].
The Jacobian is (1+x+y)^(-3), hence the integrand is
(1+x+y)^(D-3) (x+y+xy)^(-D/2). Split the square into x>=y and y>=x;
symmetry gives a factor two, and x=t,y=tv gives

```
6 integral_0^1 dt integral_0^1 dv
  t^(1-D/2) (1+t+tv)^(D-3) (1+v+tv)^(-D/2).
```

At D=12/5 the remaining endpoint weight is t^(-1/5). The substitution
t=s^5 gives the smooth bounded integrand

```
30 s^3 (1+s^5(1+v))^(-3/5) (1+v+s^5 v)^(-6/5).
```

The separate t=s^10 map gives an independent coordinate-map refinement.
`reference.py` evaluates these positive integrals by tensor Gauss-Legendre
quadrature at independently increased node counts and decimal precisions.
Its Gamma factor is computed with exact Bernoulli coefficients and a shifted
positive-real Stirling expansion; the first omitted term bounds the
log-Gamma truncation. Quadrature convergence is empirical, not interval
certification. The requested comparison is ten decimal digits.

`native.rs` reads no reference file and calls only the existing
`RecursiveBoundary` with `TadpolesOnly` and `bubble_subloops=false`. Separate
native precision/order profiles are saved before any comparison. This is a
validation-only ordinary hard coefficient, not a formula supplied to the
finite-density flow or an alternative source of its boundary constants.
