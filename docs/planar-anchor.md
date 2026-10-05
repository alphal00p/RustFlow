# Planar singular-seed regeneration

The complete 48-component planar boundary has been regenerated at the exact
Euclidean point `s=-1/10, t=-1/25, b=-1/50`, through epsilon order four. The
240 coefficients pass a conservative 40-digit absolute-error target. This is
an exact upstream singular-seed calculation with explicit GPL input evidence;
it does not use a regular grid point as the boundary.

The native calculation factors the even root valuations on
`s=-x²/10, t=-x²/25, b=-x²/50` and solves the known-leading canonical Frobenius
column before ordinary transport to `x=1`. Its root values are
`sqrt((200+x²)/50)` and `x*sqrt(1/50)`, both positive for `x>0`. All 1,200
residue equations pass symbolically with the 68 GPL constants left as independent
symbols. The 80-digit/order-80 calculation starting at `1/8` and the
110-digit/order-112 calculation starting at `1/10` differ by at most
`2.944e-82` across the serialized coefficients (the in-process check reports
`2.661e-82`).

An unchanged original `ExpCan4_Reg.wl` singular-core calculation on
`lambda=x²` independently propagates the same singular formulas to the anchor.
It requests 45 digits, reports a `3.407e-48` tail estimate, and agrees with the
high native profile to `1.050e-100`. These observed differences do not increase
the 40-digit accuracy cap.

The GPL endpoint constants are native values refined from
100 digits/order 100 to 130 digits/order 140 with different matching points.
All 68 were independently checked against the original GPL evaluator at the
regular point `31/32`; the strongest original comparison differs by at most
`7.77e-59`. Its direct endpoint helper fails on divergent suffix components;
that failed output is retained and excluded as a reference. The original planar
anchor comparison receives the native GPL endpoint constants, so it checks
propagation independently but is not an independent oracle for those constants.

A declared conservative empirical allowance of `1e-50` per complex GPL constant
is propagated explicitly. GPL dependence occurs only in component 25, through
three constants at epsilon powers two, three and four. Exact native polynomial
grouping bounds their source errors using the gradient and Hessian; a separately
refined three-order sensitivity column propagates those errors to all 240
coefficients. The largest resulting source allowance is `3.487e-48`. Adding
16 times the original reported tail, 16 times each native refinement change,
and a `1e-65` arithmetic reserve yields a maximum total allowance below
`5.81e-47`.

These are conservative numerical evidence estimates conditional on the GPL
input allowance, not outward-rounded interval certificates. The reusable
[planar-regenerated-anchor.json](../fixtures/mg5/planar-regenerated-anchor.json) carries every value and its absolute allowance with a
40-digit cap. No physical form-factor or matrix-element precision is upgraded
until those errors have been propagated to the requested observable.

The [complete planar report](../reports/validation/2026-10-04-planar-exact-seed.json)
and [GPL evidence](../reports/validation/2026-10-04-planar-gpl-constants.json)
retain all coefficients, uncertainty terms, original metadata and source hashes.
The fixture explicitly identifies the canonical component order and both
principal root germs at the Euclidean anchor.

The independent [nonplanar anchor](../fixtures/mg5/nonplanar-regenerated-anchor.json)
contains all 305 coefficients at `(-1/10,-1/5,-1)`, with eight principal germs
and a conservative `1e-40` component allowance. Its exact rational/log/pi/zeta
singular boundary needs no GPL constants. The [nonplanar report](../reports/validation/2026-10-04-nonplanar-exact-seed.json)
records two native precision/order/start-point profiles and an independent
original calculation; their largest differences are `2.041e-72` and `5.451e-51`.
Both anchors remain supplied analytic-boundary calculations, with subsequent
physical continuation and error propagation checked separately.
