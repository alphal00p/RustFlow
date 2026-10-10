# Generic same-series remainder assessment

This is a conditional analytic design, independent of the finite-eta reference formula. It supplies neither producer coefficients nor a production endpoint value. No generic 64-term endpoint accuracy certificate follows from the available moment prefix.

## A usable large-eta bound

Suppose an exact Gaussian/Schwinger owner first proves, for every contributing region and finite numerator/mass-jet term, a representation of the form

```
F(eta) = sum_j eta^lambda_j L_j[phi_j(z,1/eta)],
phi_j(z,t) = P_j(z,t) product_l (1+t Q_jl(z))^gamma_jl.
```

Here `z` is a compact parameter/radial/angular domain after a complete UV resolution; `L_j` is the explicitly continued finite-order parameter distribution including compact surface functionals. The statement must identify the actual physical integral and all branches, rather than only its formal asymptotic expansion. For all-mass shifts, projective normalization makes the virtual mass term `eta U`, and an exact Schur/channel proof can provide `|Q_jl|<=M`. Real nonnegative `Q` is an additional channel property, not a generic consequence of finite density or mass shifting. Pure compact transfer denominators must be included too. The chosen graph's bounded compact momenta alone do not substitute for the complete parameter-face proof.

If the resolved `Q`, coefficient units and their required derivatives are bounded, the circle `|t|=R<1/M` has the uniform gap `|1+tQ|>=1-MR`. Finite differentiation produces only finitely many bounded factors and explicitly shifted exponents. Thus a certified finite-order norm for every `L_j` gives a computable circle bound `C_j(R)`. This also proves a convergent normalized germ there; an unproved region expansion cannot claim that conclusion automatically.

A concrete UV norm owner can use Taylor subtraction. With exactly `N` subtracted terms of degrees `0,...,N-1`,

```
T_a(phi) = sum_(k=0)^(N-1) phi^(k)(0)/(k! (a+k))
         + integral_0^1 x^(a-1) [phi(x)-Taylor_(N-1)(phi,x)] dx,
Re(a)+N > 0.
```

Away from all retained poles `a+k=0`, Taylor's integral remainder gives

```
|T_a(phi)| <= sum_(k=0)^(N-1) ||phi^(k)||_infinity/(k! |a+k|)
           + ||phi^(N)||_infinity/(N! (Re(a)+N)).
```

Sequential finite UV axes give an iterated bound, including their mixed derivatives. This norm must include physical mass jets, upper Fermi derivatives and the genuine joint high-D lower-origin prescription. It cannot be replaced by the absolute value of a signed integrated coefficient. At fixed D it is feasible in principle; near dimensional poles its explicit denominators can be large.

For a normalized analytic channel `G(t)=sum a_n t^n`, a proved circle bound `C(R)` yields

```
|G(t)-sum_(n=0)^(N-1) a_n t^n|
 <= C(R) (|t|/R)^N / (1-|t|/R),   |t|<R.
```

This is the ordinary Cauchy coefficient estimate followed by a geometric sum; the coefficient contour formula and convergence domain are given in [DLMF 1.10.7](https://dlmf.nist.gov/1.10.E7) and [DLMF 1.10(i)](https://dlmf.nist.gov/1.10.i). It is a useful certificate at sufficiently large positive eta once the missing explicit norm is supplied. All region prefactors and normalization periods multiply the absolute error bound. Relative accuracy additionally needs a lower bound on the nonzero complete target, which signed cancellations do not provide.

## Why this does not yet reach the endpoint with 64 terms

Even if the exact representation proves a sole finite negative cut `eta in [-M,0]`, its conformal variable

```
w = (sqrt(eta+M)-sqrt(eta))/(sqrt(eta+M)+sqrt(eta)),
eta = M (1-w)^2/(4w)
```

maps infinity to zero and the physical endpoint to `w=1`. The finite-eta branch is fixed by continuation from positive eta. A Taylor bound on `|w|<1` does not contract at `w=1`. The complete normalized target may also require explicit endpoint powers before a bounded holomorphic function is obtained; that normalization needs proof and cannot be inferred from fitted exponents.

Even under the optimistic additional assumptions `|H(w)|<=B` in the unit disk and a known radial endpoint modulus `|H(1)-H(w)|<=K(1-w)^sigma`, `0<sigma<=1`, the first N coefficients need not determine the endpoint accurately. The two analytic functions `0` and `c w^N` have identical coefficients through order `N-1`. They satisfy these bounds whenever

```
|c| <= min(B, K/N^sigma),
```

because `1-w^N <= min(1,N(1-w)) <= N^sigma(1-w)^sigma` on the real radius. Their endpoints differ by `c`. This is an explicit information bound, not a numerical failure of the actual graph. A uniform norm and ordinary Holder continuity alone therefore cannot certify twelve digits from N=64 unless their constants are extraordinarily restrictive. The reference target's actual constants have not been computed.

A valid finite-N fallback could combine an interior truncation bound and an independent endpoint estimate:

```
|H(1)-S_N(w)| <= K(1-w)^sigma
                  + C(R)(w/R)^N/(1-w/R),  0<w<R<1.
```

Optimizing this certified expression would decide whether 64 terms suffice. It does not create missing coefficients or require an ODE. At present no explicit `K`, `C(R)`, complete branchwise normalization or all-region analytic identity is supplied by the coefficient generator, so this remains a design rather than an implemented bound. Higher smoothness or a proved endpoint expansion can improve the rate, but each term and remainder must be independently certified; fitting them from the same 64 coefficients does not close the argument.

## Four-loop and raised-target obligations

The construction is independent of a one-invariant bubble formula. For a generic multi-loop case it requires all compact Gram invariants, exact off-null Gaussian channels, complete UV/soft faces, a finite-order distribution norm, and every generated numerator/physical-mass jet. General sign-changing channels may have additional cuts; then the displayed conformal map is not justified. Raised targets and Fermi surfaces are signed functionals, so positive Stieltjes-measure or monotone Padé bounds cannot be assumed. Existing high-D endpoint existence proofs are qualitative and do not by themselves return the constants needed above.

The next useful bounded artifact is therefore a graph-independent norm certificate for the exact producer representation at a large-eta circle. It can certify the first series evaluation without solving the full endpoint problem. A 64-term endpoint or Laurent claim remains unsupported until the additional constants/continuation authority are available. Independent finite-eta quadrature checks remain valuable empirical validation and are deliberately separate from this certificate design.
