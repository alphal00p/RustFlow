# pflow: massless propagator integrals at fixed d0 from a one-mass flow (route B)

Module `src/pflow` (feature `native`). Front ends: `examples/pflow.rs` (flags), and `examples/dispersive.rs` with a class config
whose `"kind"` is `"propagator"`. Both call `pflow::cli::run`.

## Method
- Give one line e of a massless two-point family a mass: 1/(k_e^2 + z)^(a_e) with z = m^2/p^2. Then J(z) obeys the one-mass DE of
  RustRed (e-block: the masters with a_e >= 1).
- The discontinuity rho(s) = -(1/pi) Im J(-s + i0) obeys the homogeneous e-block. The z-free masters are z-independent sources and
  drop out.
- Where rho comes from:
  - on (0, 1): the non-integer Frobenius families at z = 0 (classes k h, h = d0/2, k = 1..L);
  - on (1, inf): the families at infinity (classes -k h);
  - their coefficients are fitted to the leading region rows (soft sets at 0, hard subgraphs at infinity; graph.rs).
- Two half-line patch-kernel transports take rho to the threshold neighbourhoods. There, local families at z = -1 (classes
  frac(j h), j = 0..2L) give the termwise integrals to s = 1.
- Boundary directions that the rows leave free are fixed by the continuity of the integer-class projection at z = -1.

## What the outputs are (E-normalisation)
E is Euclidean, measure d^dk / pi^(d/2), p^2 = 1. RustRed's I = (-1)^(sum a) E, with D_e = k^2 - z (Minkowski). For a master a with
a_e = index of the massive line, `value(m)` = FP int rho s^-m ds:
- m = 1: E[a] at z = 0, the massless value.
- m = n + 1 >= 2: (a_e)_n / n! x E[a + n 1_e]. The CLI also writes `integral` and `integral_value` = E[a + n 1_e].
- m = -k <= 0: (-1)^(a_e - 1) (a_e)_n / n! x E[a_e -> -n] with n = k + 1 - a_e. This is the e-pinched integral with line e as a
  numerator of power n; it vanishes when n < 0. `integral` / `integral_value` as above.
- `--extra a;b` gives non-master targets through their RustRed rule and the moments. The z-free masters in a rule need
  `--efree FILE`.

## Acceptance (every run)
All checks are gates with threshold 10^-T, T = `--target-digits` or by default min(w, digits of the used provided entries) - gate_margin
(gate_margin = 20, `--gate-margin`):
- boundary fits: residual of every over-determined fit;
- order budget: wf - w - 10 must cover the digits lost per order at a spurious point, and the budget must not cost rank (counted
  from the first nonzero row);
- every nonzero row representable (n >= n_lo; a leading row within the family's orders). n_lo = min(--nlo, lowest nonzero
  row - 2) per class; the threshold families use min(--nlo, -10);
- threshold families span the e-block; the threshold decompositions; the integer-class continuity at z = -1. This is an
  accuracy estimate of the threshold join: measured at four loops it tracks the precision of the outputs that depend on the
  join (no1 at w45 vs w55). At small eps the near-degenerate classes frac(j h) cost digits (~25 at |eps| = 1.3e-3 at four
  loops with a census-top basis); the decomposition's condition number is reported but is a far too pessimistic bound;
- the continuity-fit residual. Directions that continuity cannot see are UNDETERMINED: their effect on every requested output
  must be below the gate, else the run is not accepted;
- every master has boundary rows (`--allow-skipped` to override). Every irreducible lower-loop graph is provided;
- every boundary row is computable: a Gamma pole in a region or exact-series value (a log region at a non-generic d0) or a NaN
  row fails the run; it is never read as "no region" (`--allow-skipped` does not override this);
- entries of the provided table agree with the exact evaluator wherever it can reduce them (each bad entry is reported once, with
  its number of uses);
- independent checks: exactly vanishing masters are numerical zeros; masters with a closed form (graph evaluator) agree with it.

A run that fails a gate still writes its output (`"accepted": false`, failures in `diag.acceptance`), and the front ends exit
with code 3. `--accept-diagnostics` makes the exit code 0. `digits_cap` = min(w, digits of the used provided entries) bounds the digits a run may
claim. `digits_estimate` is the minimum of `digits_cap` and of -log10 of every gate quantity that measures an error (boundary
fits, continuity-fit residual, threshold decompositions, threshold continuity, undetermined directions, vanishing and closed-form
masters). `digits_by_check` lists each of them, and the front ends print `accepted`, `digits_estimate` and the limiting check on
one line. `--two-precision K` estimates the digits independently.

Transport diagnostics (information): `transport_m<m>_h<half>` per final transport and `transports_worst` (largest overlap
discrepancy, series tail and dropped imaginary part, log10 relative); `tail_heals` counts the healed patches of the run.

## Lower-loop boundary values (`--provided FILE`)
```
{"schema": "pflow-provided.v2", "normalisation": "E, P^2=1", "d0": "p/q", "master_definitions_sha256": HEX,
 "sources": [{"file", "sha256", "role"}, ...], "unverified": [ids],
 "values": [{"edges": [[u, v, a], ...], "ext": [s, t], "value": "...", "abs_err": "...", "id": ...}]}
```
- Written by `tools/pflow/answer2prov.py` from p-integral's d-exact reductions (rows from accepted lower-loop samples:
  `sample2prov.py`). v1 tables (a nominal `digits` stamp, no provenance) are refused. The table is sealed (`content_sha256`); its
  source digests are recomputed where the files are reachable (see "Content seals" below).
- Indices are integers or "p/q" strings; d0 must be the run's, exactly.
- Every value carries its ACHIEVED absolute uncertainty `abs_err`, propagated by the producer through the reduction
  (sum |c_i| |B_i| err_i plus rounding). Rows built from a source without a measured uncertainty are listed in `unverified`; a run
  that USES one fails the gate "provided values verified" unless `--allow-unverified-provided` (recorded in the output).
- Graphs are keyed by `graph::canonical`: index-0 lines and p-carrying bridges contracted, series lines merged, the minimum over
  relabellings.
- The exact evaluator runs first. An entry is used only for a graph it cannot reduce; `digits_cap` is min(w, the relative digits of
  the USED entries), and the gate "provided values meet the target" requires them to reach the target. Entries it can reduce are
  cross-checked within their own `abs_err`.

## Samples (the output, schema `pflow-route-b.v2`)
A sample authenticates itself for every consumer:
- `d0` (exact), `frame` and `frame_sha256`: the normalisation, the graph, the massive line, the master list and the boundary
  identity (`boundary`: the provided table's `master_definitions_sha256`, `answer_meaning_sha256` and normalisation, or null).
  d0-independent: all samples of one grid share it; samples built from different master definitions never do;
- `inputs`: every file read (DE blocks, manifest, graph, provided table with its sources and definitions digest, efree, class
  config) with its sha256;
- `target_digits` (the effective gate), `accepted` and `acceptance.failures` (both precision passes; a requested `--extra` that
  cannot be computed is a failure), `accept_diagnostics`, `allow_unverified_provided`, `allow_skipped`;
- every value with `abs_err` = max(10^-digits_estimate max(|v|, 10^-30 x the moment's scale), the two-precision difference), and
  `integral_abs_err`; every extra with its propagated `abs_err` (`--efree` files: schema `pflow-efree.v1` with d0, abs_err and a
  content seal);
- `downgrades`: explicit reasons why the sample is weaker than production quality (no two-precision pass: abs_err is then one
  run-wide estimate; provided sources whose digest could not be verified; `--allow-unverified-provided`, `--allow-skipped`,
  `--accept-diagnostics`);
- `content_sha256`: the producer's seal over the whole document.

## Root of trust
- The content hash (`content_sha256`) catches ACCIDENTAL changes only: a stale, mis-copied or hand-edited file. Anyone who edits a
  file can re-seal it, so a seal does not bind a result to its inputs.
- What binds a result to its inputs is the CONTRACT:
  - the fit contract pins every node's sample file by sha256, requires the frame digest (graph, massive line, masters, boundary
    definitions) and is itself identified by its `spec_sha256` in the fit output;
  - the code's legacy list (`pflowio.LEGACY_4L_NODES`) is the only route for the 12 archived four-loop v1 samples;
  - re-sealed edits fail the pin; re-sealed inputs with unverifiable sources carry downgrades, which production contracts refuse.
- Five-loop contracts (and any production contract) must be COMMITTED before their samples are used: in the campaign directory,
  with the contract file's sha256 recorded in the progress log at commit time. A fit is valid only with a contract whose recorded
  sha256 predates the fit.

## Uncertainty of a sample value (`abs_err`, `abs_err_kind`)
`abs_err` bounds the value that is WRITTEN, the one of the highest working precision w. It is the larger of:
- the run-wide estimate 10^-digits_estimate max(|v|, 10^-30 x the moment's scale);
- with `--two-precision K`: d1 = |v(w) - v(w - K)|. This difference measures the error of the LOWER pass. As a bound for v(w) it is
  safe but pessimistic, by about K digits when the error falls with w. It is also written as `abs_err_lower_pass`;
- with `--two-precision K --third-pass` (a third pass at w - 2K): d2 = |v(w - K) - v(w - 2K)|, q = d1/d2. The passes lower w,
  wf and the series orders together, so the error falls geometrically but not one digit per digit of w: measured 5 digits per
  step of K = 8 on the 2L fixture. If the measured rate log10(1/q) is at least 2 digits per step, v(w) is bounded by the
  geometric tail d1 q/(1 - q) with a safety factor 10. Otherwise, for example when an error floor that does not depend on w
  dominates (boundary digits), d1 is kept.

`abs_err_kind` names the rule used for each value. Contracts check `abs_err` against the need: the bound of the value they consume.
Production samples need at least two passes; a third pass makes the bound of the highest pass sharp instead of K digits too
pessimistic.

## Content seals and source digests (every producer, every reader)
- Canonical serialisation (`pflow::cli::canonical_json` = `tools/pflow/pflowio.py canonical`, the same bytes): object keys sorted,
  no spaces, strings JSON-escaped with UTF-8 kept, every number as its exact rational `p/q`. So the float formatting of the writer
  does not matter.
- `content_sha256` = sha256 of the canonical serialisation of the document without that field. Written by every producer
  (RustFlow samples, answer2prov / sample2prov tables, masters_from_sample, mkefree) and recomputed by every reader. Any edit of a
  value, an uncertainty, the frame or a header is refused. The frame digest is recomputed from the frame body.
- Recorded source digests (`{"file", "sha256"}` records anywhere in a header or in `inputs`) are recomputed whenever the file is
  reachable. A mismatch is refused (a stale or edited source). An unreachable source is an explicit downgrade, never silent.
- Python producers write no floats into sealed documents (strings instead), so that no reader depends on float parsing.
- Typed refusals: integers beyond 64 bits (RustFlow cannot represent them exactly: write them as strings), NaN / Infinity,
  strings with lone surrogates.

## Fit contracts (`tools/pflow/fitgrid.py`, schema `pflow-fit-spec.v1`)
- `nodes`: `{"p/q": sha256}` pins every node's ACCEPTED sample file. A different file for a node is refused, and so are missing,
  duplicate and extra nodes (never a smaller fit).
- `frame_sha256` is required; every v2 sample's recomputed frame digest must equal it.
- `production` must be `true` or `false` (null, 0 or a missing field are refused). A production contract refuses every downgraded
  sample; in particular a production sample needs a two-precision pass. The contract class is printed on every result line.
- `legacy_samples: {"p/q": sha256}`: the ONLY legacy path, for the 12 archived four-loop nodes of the blind eps^1 test, which are
  pinned in the code (`pflowio.LEGACY_4L_NODES`). Refused at loops >= 5 and in production. Every printed coefficient line of such a
  fit says `LEGACY (frame unverified)`.
- `masters_from_sample.py` needs `--frame-sha256 NAME=HEX` for every master, and refuses downgraded samples unless
  `--measurement`, which records the downgrade.

## Recognition of a coefficient (`tools/pflow/recognize.py`, GATE 3)
- The dictionary must be fixed before any value is looked at: a pflow-pslq-dict.v1 file, or generated. `--standard W` is the
  products of zeta(2), the odd zetas, zeta(5,3), zeta(7,3) and zeta(3,5,3): 20 / 27 / 36 elements through weight 9 / 10 / 11.
  `--hoffman W` is the Hoffman words.
  - A dictionary file is accepted only with `--dict-sha256`, its pinned digest: commit the dictionary and log its sha256 before
    the search; a changed file is refused. The sha256 is in the result.
- The constants come from one algorithm: Hoelder convolution of the iterated integral, with nested power sums.
  - Its two split points, 1/2 and 1/3, catch truncation errors at the full precision, but not a convention or algebra error.
  - Independently, every depth-1 factor is compared with mpmath's zeta and every depth-2 factor with a polylog quadrature, to 200
    digits.
  - Depth-3 factors have no independent evaluation; the tests check them through stuffle and duality identities and
    zeta(2,2,2) = pi^6/7!.
- Search: PSLQ on the target and the dictionary, all rounded to the search digits p1, with 20 guard digits.
  - The tolerance is n H 10^(1-p1), where H is the frozen height cap. A fixed 10^(3-p1), as in the Astra review's controls, rejects
    every true relation of height above about 10^3, whatever the precision.
  - In the blind regime p1 > (n+1) log10 H, no relation below the cap exists by chance.
- Confirmation on an INDEPENDENT value: the same coefficient from a fit of disjoint nodes (`--fit-a`, `--fit-b`), or from a
  higher-precision rerun.
  - The fit files must be sealed, have the same target, frame, scheme, loops and pole order, and share no node (d0) and no sample.
  - Fit A must carry at least p1 informative digits.
  - ACCEPTED only if the residual is below the propagated uncertainty (abs_err_B + sum |r| delta) and the independent value has
    >= 40 informative digits beyond p1.
  - A candidate that changes at p1 + 10 or p1 + 20 is AMBIGUOUS. Search digits are never used for the acceptance.
- Measured search digits (published controls, H = 10^12, s1 campaign 18-19 Nov):

  | control | true support | blind, standard products | blind, Hoffman words |
  |---|---|---|---|
  | weight 9 (GGPPSS int200, hatted zeta(3,5)) | 70 | 180 (blind bound 252) | not found <= 460 |
  | weight 10 (LSS M6,2 eps^3) | 150 | 320 (bound 336) | not found <= 460 |
  | weight 11 (LSS M6,2 eps^4) | 230 | not found <= 700 | not found <= 460 |

  - At weight 11 the true relation in the standard basis has height 10^12.4, above the cap. A blind weight-11 search needs
    H >= 10^13 and about 37 x 12.4 ~ 460-480 search digits, plus 40 for the confirmation.
  - A Hoffman-word dictionary is impractical: the coefficients of physical quantities in it are far larger.

## Status and limits (19 Nov 2026)
- Validated through FOUR loops:
  - the P4L_NO1 family end to end: boundary from p-integral's d-exact reductions and our own 3L samples, route B, fit;
  - the blind eps^1 test of no1 and no6, about 15 digits, with the exact c_-1 and c_0 reproduced;
  - every acceptance gate, the provenance chain and the uncertainty propagation.
- FIVE loops: NOT run end to end. No T266 master list, no five-loop DE block and no T266 boundary request exist yet. The T266
  reconstruction was paused on 19 Nov (it needs a dedicated machine and engine changes).
- Gates:
  - GATE 1 (input authentication): done, on main. Seals, pins, frame digests, confined legacy path, typed refusals, production
    contracts with two-precision samples.
  - GATE 2 (digit budget of the coefficient): partly done.
    - Done: the exact amplification and the propagated abs_err per coefficient (fitgrid), the third-pass bound of the highest pass.
    - Design only, not implemented: the truncation bound from the fitted orders, hold-out prediction nodes and the disjoint fit B
      in the driver.
  - GATE 3 (blind recognition): the tool and its measured requirements (above). Not yet applied to a five-loop coefficient.
    - The reviewer's items to fix before its first production use are done (19-20 Nov): sealed and matched fit files, disjoint
      nodes, digits of fit A, the dictionary pin, the independent depth <= 2 constants.
    - Left: no independent evaluation of depth >= 3 constants (identity tests only). The weight-11 test target uses the tool's own
      zeta(2,1,8).
- Known limits:
  - The threshold continuity is a conservative accuracy estimate. On the 4L default basis it is 15-50 digits below the
    closed-form masters.
  - At small |eps| the cost of the Frobenius families grows, strongly on the minimal (spurious-free) basis: no family finished in
    23 h at eps = 6.25e-6, w = 200.
  - Provided values the evaluator can reduce are cross-checked against the evaluator at TWICE the run precision.
    - The entry is re-read at four times the precision.
    - The tolerance is the entry's abs_err plus 2 x that evaluation's own error (measured against four times the precision).
    - Entries wrong far below the run precision are caught (1e-250 at 181 digits), while the evaluator's cancellations at small
      eps (11-15 digits) do not flag correct entries.
  - The evaluator's measured cancellation loss enters `digits_cap` (diag `provided.evaluator_lost_digits`).

## Checkpoints (`--workdir DIR`)
- Frobenius families write their exact state every 900 s. Every transport writes after every patch.
- Each precision pass uses its own subdirectory `p<w>-<wf>-<m>`. The family fingerprint is the DE content (d0, masters, every
  term).
- Rerunning the same command resumes.

## DE blocks
`--de` takes one RustRed harness solution, or several per-cut BLOCKS on one common master basis: a comma list, or a manifest
`{"blocks": [...]}`.
- The blocks are united.
- RustRed's identity rules (master = 1 x master) are ignored.
- The same target must have the same rule in every block.
- A master of one block must not be reduced by another block.

## Tests
`cargo test --release --lib pflow:: -- --include-ignored` (about 10 minutes at 2 threads without the four-loop regression; that
one takes 30-60 minutes).
Tools: `python3 -I -B -m unittest discover -s tools/pflow -p 'test_*.py'` (mpmath; seconds). They include the attacks of the
GATE 1 review: edited values, uncertainties and frames, stale sources, unpinned or foreign samples, legacy at five loops.
