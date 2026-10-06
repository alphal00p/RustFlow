# Linear propagators and FT recursion

`solve_integrals` and `solve_integral_projections` honor `RecursionMode::Ft`
for families containing linear denominators. They use the same recursive
Feynman-parameter evaluator as quadratic families: combine two active
denominators, derive and close their differential equations with the configured
reducer, recursively determine a regular interior boundary with fewer active
denominators, integrate the weighted power/log series, and finish at the existing
Gaussian terminal. Exact target weights are applied before epsilon fitting.
No eikonal formula is used as a numerical terminal or family dispatch rule.

The existing FT domain check remains mandatory. RustRed constructs the exact
Symanzik polynomials from the full affine scalar-product forms; coefficients of
`U` and `-F` in the active sector must be real and nonnegative. This is a
sufficient parameter-domain criterion, not a classification of every physical
Euclidean momentum configuration. It also admits mixed quadratic/linear forms.
Inputs that fail the criterion return `Error::Unsupported`; resource exhaustion
and incomplete reduction retain their existing typed errors.

AMF retains its separate `D(x)=D+x(u·l)²` deformation followed by a Frobenius
endpoint projection. Each linear denominator must have one real rational loop
branch `u`. The pinned original AMFlow implementation uses the same factorized
branch requirement in `RawBranchOf`/`GenerateSquare`. FT combines the complete
denominators directly and consequently does not require this rank-one condition.
The change introduces no graph algorithm or tensor machinery: exact algebra,
Symanzik construction, IBP reduction, Gaussian integration, and the series engine
remain with their existing Symbolica, RustRed, and RustFlow owners; Gaussian
tensor connectivity already uses native Linnet.

The integration regression uses two independent eikonal directions,

```text
I(eps) = integral d^D q/(i*pi^(D/2))
         1 / [(q²+i0)(2q·v-2+i0)(2q·w-2+i0)],
D = 4-2eps, v² = w² = 0, v·w = 1.
```

Positive Feynman parameters `t,u` give the completed-square mass
`M²=2t+2u+2tu`. Integrating `t` and then `u` yields
`I(eps)=-2^(-1-eps) Gamma(eps)^2 Gamma(1-eps)` in the convergent strip
`0<Re(eps)<1`. Raising the first linear power multiplies this by `-eps/2`.
These expressions are independent test oracles. The tests use `eps=1/3,2/5`,
increased working precision/order, the public Laurent API through `eps^0`, and
an epsilon-dependent projection. Changing `v·w` to `-1` makes
`M²=2t+2u-2tu`; the production Euclidean guard rejects it.

A two-loop regression uses `l1²`, `l2²`, and
`2l1·v+2l2·w-2`, with `v²=w²=1`, `v·w=0`. Its linear form has rank two
in loop/external space and is deliberately rejected by AMF's rank-one
deformation. The active parameter polynomials are
`U=a1*a2`, `-F=2*a1*a2*b+(a1+a2)*b²`, so FT admits it directly.
Independent Schwinger integrations give
`-2^(3-4eps) Gamma(1-eps)^2 Gamma(4eps-3)` in the strip
`3/4<Re(eps)<1`. The regression uses `eps=7/8`, precision/order refinement,
and a mixed loop routing with determinant `-2`; it checks the associated
`2^D` integration Jacobian explicitly.

The live reference is original AMFlow 2.0
`26005517a288086c4cb4d1b26d829691bc088485` with its FT recursion, Wolfram
15.0.1, Kira 3.1 and Fermat. At `eps=1/3` it returns
`-3.85663610547731999372773380013725693079` with **24 recorded digits**.
Its scaled residual against the gamma formula is `5.36e-27`; the acceptance
comparison asks for 20 digits and never treats the longer printed mantissa as
additional precision. The original automatic solve took 46.833677 seconds
with one worker, 48 working digits and order 96 plus extra order 50. This is
an oracle timing, not a matched native performance claim. Source hashes,
settings and the kernel log are retained in
[the reference report](../reports/validation/2026-10-06-linear-ft/oracle-metadata.json).

Regenerate the reference with the existing optional runner and
`--case eikonal --mode sample --recursion FT`; see
[the runtime setup](upstream-oracle.md). Ordinary Rust regressions read the
retained output and require neither Mathematica nor an external reducer.

This extension does not implement arbitrary mixed-sheet cuts, general cut
boundary recursion, physical-region FT, or a new prescription inference rule.
