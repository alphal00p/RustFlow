# Native boundary-integrand decomposition

`integrand::to_integrals` delegates rational-coefficient partial fractions,
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

The native path currently admits real rational functions of scalar parameters
in denominator coefficients. Gaussian numerators also work: native rewriting
uses exact substitution, and the existing exact coefficient adapter extracts
the final integral powers. Denominators outside the native rational coefficient
field retain the previous exact Symbolica `AtomField` path. In particular, the
imaginary unit cannot be replaced by a free parameter when computing rank:
rows proportional over `Q(i)` need not be proportional over `Q(t)`.

Admission is structural and happens before native decomposition. A native error
or exhausted state budget never activates the other path. Native partial-fraction
state and power limits return typed limit errors; malformed native families
retain typed input errors. A returned partial-fraction weight must be independent
of the loop scalar products. No partially accumulated result escapes an error.

The shared HEPKit revision is
[`8bfd027`](https://github.com/ValentinHirschi/gammaloop/commit/8bfd027a8df276ab640adf1f238833bbcf9bb9c1).
It fixes affine extraction and partial fractions using native exact rational
polynomials, preserving original row denominators. The previous owner could
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
