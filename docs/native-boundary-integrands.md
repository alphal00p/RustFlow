# Native boundary-integrand decomposition

`integrand::to_integrals` delegates Gaussian-rational partial fractions,
independent sectors, scalar-product completion and numerator rewriting to
HEPKit's native `IntegralFamily`. This is the same family owner used for native
diagram input. The adapter changes notation and returns the existing RustFlow
`IntegralTerm` objects; it does not add a second native decomposition algorithm.

Input scalar-product coordinates must be distinct symbols spanning the declared
loop/external space. Their grouped loop-loop then loop-external order is mapped
by exact product identities to HEPKit's interleaved order. Completed families
preserve the caller's loop/external names, Gram matrix, dimension, epsilon symbol
and normalization. Original denominator scaling stays in the numerator weight.
Internal momentum, dimension and denominator symbols are chosen to avoid every
input coordinate and scalar coefficient, including successive suffix collisions.

The native path admits rational functions of scalar parameters with exact
Gaussian-rational number coefficients. HEPKit uses Symbolica's native
`AlgebraicExtension::complex(Q)` field, enforcing `i^2 = -1` during rank,
partial fractions, completion and rewriting. Gaussian numerators also work:
native rewriting uses exact substitution, and the existing exact coefficient
adapter extracts the final integral powers. No free imaginary parameter is
introduced: rows proportional over `Q(i)` need not be proportional over `Q(t)`.

Functions, noninteger powers and number coefficients outside the admitted
Gaussian-rational field retain the existing exact Symbolica `AtomField` path.
The reserved formal imaginary parameter also remains outside native admission.
This conservative structural boundary does not claim every coefficient
expression supported internally by the owner. Regressions distinguish both
routes through their limit errors and verify exact reconstruction.

Admission is structural and happens before native decomposition. A native error
or exhausted state budget never activates the other path. Native partial-fraction
state and power limits return typed limit errors; malformed native families
retain typed input errors. A returned partial-fraction weight must be independent
of the loop scalar products. No partially accumulated result escapes an error.

The shared HEPKit revision is
[`b96600b`](https://github.com/ValentinHirschi/gammaloop/commit/b96600b0085d9ddfa9e6acbc11fa72ec6163253c).
It extends exact affine extraction and native matrix operations to Gaussian
rationals, preserving original row denominators. The earlier extraction fix
also covers momentum and Symanzik coefficient grouping. The previous owner could
lose an expanded symbolic coefficient such as `(a+b)/10^1000` during statistical
zero testing and return a loop-dependent weight. The direct regression checks
both scalar weights and exact reconstruction; reconstruction alone would miss
that failure. All native crates share one pinned owner revision; Symbolica and
RustRed revisions are unchanged. Source-sensitive cache fingerprints therefore
change, while the binary schema does not.

This is an ownership and exactness improvement within the boundary recursion.
It does not establish automatic evaluation of every topology or remove the
remaining linear/cut, coefficient-field or endpoint limitations. The isolated
validation report records its measured source and numerical scope separately
from the community publication wheel.
