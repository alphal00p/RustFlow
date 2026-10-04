# MG5 mixed QCD–EW plugin: transport cache audit

Pinned clone: `/common/dev/mg5_higgs_ew_plugin`, master commit
`89f64b93d0bdbd8ee90eb85737021189229b034a`. No root LICENSE/COPYING/NOTICE file
was found; the bundled ALOHA templates name an ALOHA license which is not present
as a separate license file. No computational source is being copied into the
Rust library. The two scientific boundary grids are Git LFS objects; download
sizes and verified SHA256 hashes are in `grid-downloads.json`. Both archives
were downloaded without modifying the source checkout (3,038,906,949 bytes total).

## What the application actually does

The EW real-emission backend evaluates four complex form factors for Higgs plus
jet production, assembling a planar and a nonplanar family at four permutations
of the Mandelstam invariants. Its master-integral bank stores *epsilon
coefficients*, not evaluations at a small epsilon: the driver sets epsilon order
4 and substitutes selected coefficients into `FFlist`. That formula references
planar MI orders 0–4 and nonplanar orders 1–4. Its highest referenced MI indices
are 47 and 61 respectively; these are not a proof of full DE dimensions.

The independent coordinates are dimensionless `(s/MV²,t/MV²,mH²/MV²)` with
`u=b-s-t`. The bridge selects W or Z and snaps the third coordinate to
`13074/5399` or `14631/7775` within a tolerance of `1e-4`. It rejects a negative
Gram determinant and uses four invariant permutations. This is specific physics
normalization, not an interchangeable convention for a general integral bank.

`ExpCan4_Reg.wl:641–651` first imports the boundary-point index associated with
the family and W/Z label, falling back to a general family index. Lines 657–681
prefer candidates with the same coordinate-sign pattern as the destination;
if none exist they use all indexed points. Mathematica `Nearest` ranks raw
Euclidean distance in these dimensionless coordinates. `iBClist=1` selects the
first candidate. There is no comparison of predicted segment count, inherited
error, or requested accuracy at selection time. `FindContourMV` at lines 573–575
builds a straight line, after which singularity and branch checks run. Thus a
nearby starting point is a heuristic for cost, not a guarantee of a cheap or
valid continuation.

Boundary files are named by family plus Mathematica `Hash[coordinates]`. Each
contains `{coordinate rules, full MI epsilon-coefficient matrix, delta}`.
`SetGlobalV` (594–615) verifies exact coordinate equality, loads the inherited
error `delta0`, sets root/threshold sign information and prepares the contour.
`ContourExpandMV` adds per-segment estimates to that inherited error and rejects
when the total exceeds the requested threshold (528–531). This is an estimated
error, not an interval certificate. The filenames do not encode a matrix/basis
hash, normalization, precision policy, epsilon coverage, or a branch-path ID.
The folder and driver settings implicitly provide part of that identity.

The library's saving routines support uniform samples `x=k/nsp` along completed
transport, selected from the computed local series and avoiding singular
expansion centers. They save complete MI coefficient vectors and error estimates
and optionally append coordinates to the grid index (684–725). **The distributed
application defaults are different from a growing intermediate-point cache:**
`expewmi.m:16` and its FIFO counterpart set `$NSave=1`, so only the destination is
saved; both top-level scripts set `$grid=False` and delete newly computed endpoint
files after assembling form factors. Existing points are reused. Setting grid
mode and a larger `nsp` supports persistent densification, but saving does not
update the currently loaded in-memory candidate list inside that batch.

There is also a separate Fortran form-factor cache: 1000 entries keyed by the
three gluon momenta and the two masses, using a fixed absolute `1e-13` comparison.
It caches the four complex output coefficients, not master integrals or
intermediate transport states. The insertion code saturates its size counter,
so after filling the table it repeatedly replaces slot 1 rather than performing
an actual circular eviction. This should not be copied into a new cache design.

FIFO mode keeps Mathematica processes alive to avoid repeated launches and
multi-kernel lockups; it does not change the default endpoint-deletion policy.
The grid filler preserves input-row order, skips previously filled rows, and
leaves failed entries as NaN for later retries. These are batch/MC scheduling
semantics, separate from nearest-boundary transport.

## Native architecture to implement

Make a long-lived physical-transport session own a point bank by default.
Prepared symbolic systems and the point bank are separate cache layers: the
former avoids algebraic setup; the latter avoids transporting repeatedly from
one distant origin. Save accepted regular checkpoints generated during actual
transport, with retention budgets; do not require a separate grid-generation
phase. Keep local series too when useful, together with their certified or
estimated convergence domain and achieved error. Only promote entries after
required numerical validation succeeds, or mark provisional entries explicitly.

Compatibility must include canonical DE/basis and variable identities,
normalization, dimension and epsilon conventions/coverage, physical prescription,
algebraic-root sheet and continuation class, and format/algorithm versions.
Entries hold exact coordinates or exact encodings of declared input precision,
full basis coefficient vectors, per-component errors and validation status,
precision/order provenance, and parent-anchor lineage. Requested accuracy must
participate in numerical reuse; a higher working precision alone does not prove
that a previous result satisfies a stricter request.

Use a spatial index to shortlist compatible points, then rank predicted work:
segment count from distances to singularities, necessary precision uplift,
threshold crossings, coupled-block cost and reuse of an existing local series.
Carry inherited error through the new transport, and reject seeds with
insufficient accuracy; fall back to another compatible seed if a path is
invalid. A Euclidean nearest seed may be more expensive than a farther seed in
a better-conditioned region. Deterministic tie-breaking and diagnostics should
record the selected anchor and why candidates were rejected.

An exact-coordinate, compatible, sufficiently accurate hit can return directly.
A nearby point supplies a *boundary*, followed by transport to the exact requested
destination. Do not substitute a nearest-neighbor value or silently interpolate
in a Monte Carlo integrand: that changes the sampled function and can bias an
observable. Keep MI caching independent of couplings that enter only later
prefactors; amplitude-result caching must include every coupling/model parameter
that affects that amplitude. Numerical failures remain typed failures, suitable
for explicit retry, rather than accepted zero or NaN cache entries.

Persist point blobs atomically under content-based identities, validate their
checksums, and rebuild/update a compatible point index transactionally. A
process-safe writer protocol must prevent half-written boundaries and index
entries. Merge duplicate-coordinate results by validated accuracy and epsilon
coverage; keep differently normalized or different-sheet entries separate.
Eviction should reflect cost saved and spatial coverage, not only recency.

## Meaningful cache regression

`cache-selection-fixture.json` uses the application's exact W/Z coordinate
normalization and same-sign physical domain, but an explicitly independent
analytic scalar DE rather than falsely labelling it an extracted EW master.
First transport from s=4 to s=12 and retain an accepted checkpoint at s=8; then
request s=79/10. Closer wrong-basis, wrong-sheet, wrong-mass and inadequate-error
candidates must be rejected, the s=8 anchor should be selected over the original
s=4 seed, and the final answer must agree with a fresh calculation. A stricter
accuracy request and a changed normalization must not reuse an incompatible
numerical entry. The grids now provide actual EW-family subblocks, as described below.


## Verified scientific data and actual closed subsystems

Both LFS archives match their expected SHA256 hashes. They contain 32,699 and
32,703 entries (nonplanar and planar), with approximately 3.95 GB and 2.92 GB of
uncompressed data. `*.inventory.json` records every archive member and size.
`selected-grid-files.json` records hashes of the selectively extracted matrices,
precanonical functions, W/Z indexes, and three representative boundary files per
family. No license/readme file is present inside either archive.

The full connection matrices have dimensions 61 (nonplanar) and 48 (planar).
`actual-rational-subsystems.json` exports two closed blocks and proves their rows
have no dependency on omitted masters. Planar master 7 satisfies
`dJ/ds = -epsilon*J/s`, with zero t and b derivatives. Nonplanar masters 1 and 2
satisfy `dJ/ds = epsilon*A(s)*J` with
`A=[[2/(1-s)+1/s,1/(1-s)],[4/(1-s)+4/s,2/(1-s)]]`; their t and b derivatives
also vanish. These are scientific matrix entries, not copied solver code.

Both exports include the same exact positive-s W-grid point, coefficients for
epsilon orders 0–4, source hashes, and inherited grid error estimates. Those
estimates are about 1.63e-27 (planar) and 3.53e-28 (nonplanar); stored numerical
mantissas have over 200 digits, which is not an accuracy certificate. The planar
scalar's physical branch has `Log[-s]=Log[s]-i*pi`. Its normalization includes
`-8*zeta(3)/3` at order 3 and `-pi^4/30` at order 4, so a pure
`exp(-epsilon*Log[-s])` expression alone is not the full normalized integral.
`actual-scalar-validation.json` compares all five coefficients against the
supplied analytic-origin expressions with a 20-digit assertion. The initial
normalization-only exploratory failure is preserved separately.

These blocks allow actual-data point-bank tests without loading a full EW
amplitude: for positive s, transport a coefficient vector using the closed
rational system, save intermediate regular points, then choose a compatible
checkpoint for a new destination and transport from that exact point. The
scalar has the exact propagation multiplier `exp(-epsilon*log(s1/s0))` along a
positive real interval, giving an independent coefficient-convolution check
that preserves the source normalization constants. This validates transport and
cache semantics; full EW form factors additionally require the full basis,
precanonical maps, four invariant permutations, and the amplitude prefactors.


## Original DiffExp check on the actual NP grid block

The runnable oracle and machine-readable fixture are under
`../mg5-np-pointbank/`. The fixture uses independent coordinates s,t,b with
explicitly dependent u=b-s-t, plain decimal complex coefficients, and exact
rational point coordinates. Starting from a physical W-grid point, original
DiffExp transports s0→s0+1→s0+2 and directly s0→s0+2, at working precision/order
80/80 and 100/112. All ten epsilon coefficients agree between routes and between
settings at the asserted 24 digits. The original grid error estimate is carried
as an upward-rounded 1e-27 estimate per coefficient; all input/output accuracy
claims remain capped at 24 digits. Exact returned coordinates and the analytic
first-epsilon increment `[1,2]*log((s1-1)/(s0-1))` were checked separately.
`provenance.json` preserves commands, hashes, raw results and individual phase
timings. This is a physical point-bank transport fixture, not full 61-master
transport, form-factor evaluation, or full amplitude acceptance.
