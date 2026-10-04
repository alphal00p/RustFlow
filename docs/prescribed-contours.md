# Certified polynomial-prescription contours

Implementation: `src/contour.rs`; tests: `tests/prescribed_contours.rs`. It accepts exact real univariate singularity polynomials and exact real polynomial i0 declarations. All roots are isolated through native Symbolica; root identities are compared exactly. Shared simple real roots require compatible prescription sides, computed from the certified sign of the declaration's derivative. Unprescribed real roots use the explicit default side. Multiple prescribed roots are unsupported; off-axis roots remain obstacles.

Every input passes the repository's exact scalar validator before conversion, so floating coefficients, functions, undeclared symbols and complex coefficients are rejected. Polynomials must have constant rational denominators after conversion. A zero polynomial or singular endpoint is an input error. `plan_with_context` checks cancellation before/between/after native operations; an ongoing native isolation/refinement call cannot be preempted by the token.

The planner retains each complex root's original certified disk. It only refines a real root when its interval cannot resolve endpoint membership, derivative sign, or required detour width. No high-precision refinement of all complex obstacles is performed.

## Geometric argument

For the actual rounded real center c and an obstacle disk with rational center z and radius r,

`max(abs(Re(z)-c), abs(Im(z))) - r`

is an exact lower bound on distance from c to its root. The minimum D of these bounds and endpoint distances determines a conservative detour radius `(D - margin)/8`. Nonpositive/unresolved clearance returns a typed accuracy error. The chosen radius must exceed the selected real-root interval uncertainty plus the rounding margin; this selected certificate can be narrowed if needed.

Each *actual rounded triangle vertex* is then checked exactly to have L1 distance less than D from c. Convexity proves the entire triangle excludes every other root. The real entry/exit must bracket the whole certified interval of the chosen real zero, and the apex has the requested nonzero imaginary side; the two open nonreal segments therefore avoid that zero as well. Finally, detour real projections must not overlap. Connecting segments lie on the real axis between certified real crossings. Reverse travel uses the identical geometric path.

The disk and vertex checks carry the proof; the coarse precision margin alone is not treated as an enclosure certificate. Conservative disk overlap or insufficient precision can return `Error::Accuracy` even if a different contour would exist. The real-endpoint/sign/width refinement attempts are bounded, but native exact root calls currently have no runtime cancellation interface.

## Validation and upstream limits

The five previous geometry regressions pass with the disk-based planner, including the close polynomial `(x-1/2)^2+10^-100` at120digits, complex roots of a declared polynomial, shared exact factors, conflicting signs, endpoints and reversed paths. Direct tests of exact analytic disks exercise the rounding-margin collapse and subtraction of nonzero obstacle radii independently of native root isolation.

The complete polynomial test for `delta=16/10^43+1/10^68`, `Q=(x-1/2)^2+delta^2` remains explicitly ignored in the default suite because native `isolate_roots` itself exceeds a60-second bound. It is preserved unchanged apart from that explicit annotation, and a separate Symbolica-only bounded reproducer is preserved in `repros/symbolica-close-root-isolation`. This is a distinct issue from the prior high-precision `refined` slow case. Nothing merges its roots or relaxes the geometry assertion.

The production planner calls the existing `family::scalar_symbols` validator. Before integration, its standalone harness linked a validated RustFlow library and Symbolica75f8350; `reports/validation/2026-10-04-prescribed-contours.json` records those actual results and hashes. Repository-wide checks are recorded separately by the integration milestone.

## Analytic continuation acceptance

The final standalone suite reports9passed,0failed,1explicitly ignored native isolation case (0.38seconds reported by the Rust test runner). For

`y'=[1/(2(x−1/5))+1/(4(x−4/5))] y`, `y(0)=1`,

with declaration `(x−1/5)(x−4/5)+i0`, the planner supplies a lower first detour and upper second detour. Existing rational `CompiledSystem::transport` returns the analytic value `y(1)=1+i`; reverse transport along the planned reverse path restores1. This test uses actual adaptive numerical transport, not only geometric assertions.

At60decimal working digits/order64, forward analytic error is `7.2266455081721967473902489042803014490513080865204152460954861459e-56` and reverse error `7.7016623278192015288109854161468849416394547201344689802358137825e-56`, with186steps each. At80digits/order96, forward error is `2.110871649630927717126635286049126356435629632788602355446738354384838804013052230126e-78` and reverse error `1.115742755424132759974818468872841764796955067113857603062496773854658200397575427366e-78`, with100steps each. Cross-precision difference is `7.226645508172196747390058703763194269045603198329177453364529619248138368226812408016e-56`. Requested error controls were24and36digits; extra agreement is this analytic test's observed result, not a general working-precision guarantee.
