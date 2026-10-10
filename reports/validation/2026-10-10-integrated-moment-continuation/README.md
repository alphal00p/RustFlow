# Same-series continuation fallback

The prescribed fallback improves finite-eta evaluation but does not resolve the endpoint. It uses the same frozen 64 source-derived coefficients for both original targets, with R = 4, w = R/(eta+R), and an exact binomial composition. All 128 coefficient inverse checks pass. The 72 partial-sum/Padé predictions were saved before reference comparison; generation took 1.62 s within the 120 s cap.

At eta = 2, the 64-term transformed sums have relative errors 2.41e-16 (scalar) and 5.20e-17 (raised). At eta = 8 both also pass. The fixed source-based circle bounds independently control exact-series truncation: at eta = 8 the scalar/raised relative bounds are 3.135e-28 / 6.609e-27; at eta = 2 they are only 2.839e-8 / 3.959e-7. The stronger observed eta-2 agreement is empirical. These bounds use exact rational upward root enclosures and exclude decimal master/normalization rounding.

The endpoint fails: 64-term errors are 5.22e-3 (scalar) and 2.19e-4 (raised). The prescribed 31/32 Padé endpoint errors are about 5.24e-5 for both. All 24 endpoint reference checks fail the prospective 12-digit criterion. Four endpoint Padé precision checks fail as well, at orders 24/24 and 31/32. No added order, precision or fitted profile replaces these failures.

The scalar norm uses positive parameter/compact weights. The raised norm is separately derived from absolute C2 measure, original-numerator, kernel and upper-surface components; it never substitutes the absolute signed first coefficient as a norm. Both geometric bounds require w < r < 1. At eta = 0, w = 1, no such radius exists. Details and exact inequalities are in `bounds.md`, `scalar-bound.json`, and `raised-bound.json`.

Padé denominator values and unused-prefix residuals are recorded, but no pole-free path or all-order identity is asserted. This report makes no certified endpoint, Laurent, full three-loop amplitude, or four-loop claim. The separate scalar guessed-ODE experiment is not an input to this fallback.

`summary.json` records all final profiles and failures. The first comparator rejected the historical reference's integer-zero input counters because it expected boolean false; its original source and failure record are retained. Only this metadata check was corrected, and predictions were unchanged.
