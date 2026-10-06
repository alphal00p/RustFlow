# Registered-root singular endpoints

`AlgebraicSystem::prepare_frobenius(max_dimension, context)` prepares a local
power/logarithm expansion at zero, including branch points of its registered
square roots. `pullback_coordinate(z, c+z)` centers an endpoint at `x=c`;
`pullback_coordinate(z, 1/z)` treats infinity. Both substitutions are exact and
include the path Jacobian. The existing regular transport/cache APIs continue
to require regular points.

The adapter uses restriction of scalars. If `r_i²=R_i(x)` and `Y'=A Y`, adjoin
the vectors `m_S Y`, where `m_S` is a product of distinct registered roots. Then

```
(m_S Y)' = (sum_{i in S} R_i'/(2 R_i)) m_S Y + m_S A Y.
m_S m_T = (product_{i in S intersection T} R_i) m_(S symmetric_difference T).
```

This is a rational differential system, even when no simultaneous change of
coordinate rationalizes the roots. Only the subgroup of root monomials generated
by the exact normalized connection is included. Unused radicals do not enlarge
the system; a product of several roots can require only two monomials. The
matrix dimension is `components * epsilon_orders * generated_monomials`, at
most `components * epsilon_orders * 2^registered_roots`. The caller's dimension
limit is checked before dense allocation and as the subgroup grows. The current
bit-mask representation allows fewer than `usize::BITS` registered roots.

`rational_lift` exposes this exact connection for reuse with existing solvers.
Endpoint preparation uses the existing SCC decomposition, meromorphic Fuchsian
normalization, exact indicial exponents and coupled resonant logarithmic
recurrence. It does not add a separate series implementation. Exact polynomial
and quotient normalization belong to Symbolica; no graph or tensor operations
are introduced. This independently derived product rule implements the class
of root/logarithm expansions discussed in the DiffExp paper, without translating
the GPL reference implementation.

The prepared object's `match_boundary` accepts an `EpsilonBoundary`, the same
`RootSeed` map as regular algebraic transport, working precision, expansion
order, local logarithm winding and `RunContext`. For an ordinary system, use one
epsilon coefficient row. The common Laurent leading power is retained. At the
regular matching point, exact radicand substitution and the existing root seed
logic select every `m_S Y`, fixing the physical solution inside the larger
rational solution space. This also works with dependent radicals when the
requested denominator is invertible in the existing formal root algebra.
Nonunit denominators keep their existing typed rejection.

The result, `AlgebraicEndpointExpansion`, evaluates original integral components
and extracts finite endpoint limits of the retained epsilon coefficients through
`coefficient_limits()`. Divergent auxiliary root-monomial components are not
projected into that limit. A divergent original coefficient returns an error.
The endpoint is extracted from exponents and coefficients; no small numerical
mass is substituted.

The order of limits matters: `EpsilonSystem` has already expanded in epsilon,
so this adapter cannot identify symbolic dimensional sectors such as `x^epsilon`
afterwards. For `Y'=epsilon/x Y`, its logarithmic epsilon coefficient therefore
causes `coefficient_limits()` to fail, rather than being silently discarded.
Dimensional-sector selection requires retaining symbolic epsilon exponents before
specialization, as supported by the rational `PreparedFrobenius` interface; this
algebraic coefficient-hierarchy adapter does not infer that missing provenance.

Both matching and evaluation use `log(x)+2*pi*i*winding` for **all** fractional
powers and explicit logarithms. This is winding around the local expansion
center. Root seeds at the matching point and subsequent winding must describe
the caller's intended path. The adapter does not infer a physical prescription
or a global homotopy. Matching/evaluation must remain inside the local
convergence domain, excluding other singularities.

An expansion carries working-precision values, not verified uncertainty or
cache evidence. Repeat preparation evaluation/matching with fresh boundary
data at higher precision and increased order to establish numerical stability.
Endpoint cancellation checks inherit the rational solver's numerical thresholds;
they are not interval proofs. A low-accuracy input boundary is not promoted by
increasing working precision. Irregular exponential sectors, indicial roots
outside the rational solver's support, and systems beyond its Fuchsian
normalization ability still return typed failures.

Analytic tests cover `exp(±2 sqrt(x))` and monodromy, coupled
`sqrt(x) log(x)` sectors, finite epsilon hierarchies, dependent roots with
different signs, a multi-root product with an unused generator, infinity with
divergent auxiliary rows, retained source holes, resource limits and
cancellation. These local tests do not establish general large-system
performance parity or automatic singular-point insertion into the physical
boundary cache.
