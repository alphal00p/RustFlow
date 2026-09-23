# Massless two-point subloops

The native adapter can integrate a massless subloop with denominators
`K²` and `(K+q)²` before asking RustRed to search IBPs. It expresses the result
in terms of the same subloop with powers `(1,1)`, so the noninteger power
of `q²` remains inside an ordinary integral in the original family.

For positive integers `a,b`, consider a numerator
`(K·v₁)…(K·vᵣ)`. Let `Wⱼ` sum all choices of `j` disjoint pairs of labelled
vectors: every pair contributes `vᵢ·vₖ`, and every unpaired vector contributes
`q·vᵢ`. Relative to the scalar `(1,1)` subloop, the tensor integral is

```text
Σⱼ Cⱼ (q²)^(j+2-a-b) Wⱼ,  0 ≤ j ≤ floor(r/2)

Cⱼ = (-1)^r / 2^j
     × P(2-D/2, a+b-j-2)
     × P(D/2-1, 1-a+r-j)
     × P(D/2-1, 1-b+j)
     / [(a-1)! (b-1)! P(D-2, 2-a-b+r)]

P(z,n) = Γ(z+n)/Γ(z).
```

This follows by combining the two propagators with a Feynman parameter,
shifting `K` to complete the square, contracting the even shifted tensors,
and integrating the remaining beta function. Since every shift in `P` is
integer, the implementation uses finite products and reciprocals, preserving
exact rational coefficients. Examples are `I₂₁/I₁₁ = -(D-3)/q²`, the vector
projection `-q·v/2`, and the rank-two projection
`[D(q·v)²-q²v²]/[4(D-1)]`.

The implementation reconstructs the momentum routing from the quadratic
forms, shifts numerator factors, removes powers of `K²`, and converts the
remaining scalar products back to existing denominator coordinates. Terms
with a nonpositive bubble power are scaleless one-denominator integrals.
The transfer square must be external or equal an already active propagator;
otherwise this optimization declines the sector. A propagator with a different
mass would require a different family and is not silently introduced.

Pure remaining-loop coordinates precede bubble-dependent coordinates in the
native integral order. Every proposed analytic rule must strictly decrease
that exact order. A rule that fails this check is left to native IBPs, as are
unsupported routings and exceptional rational dimensions. This order is also
used by RustRed when it searches the sector, preventing cycles between the two
sources of rules. Cache migration rejects previously searched sectors whose
order would change.

Tests compare analytic tensor identities and the integrated adapter against
independent native IBP reductions, including a two-loop insertion. Additional
tests check the ordering of raised powers and rank-three numerators in the
original paper's deformed family.
