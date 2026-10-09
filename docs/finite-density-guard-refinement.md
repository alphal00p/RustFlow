# Native finite-guard refinement

`WeightedClosureOptions.guard_refinement` submits narrower integer domains to
the existing RustRed guarded source engine. It does not introduce another
elimination algorithm, change integral ordering, or supply a master list. The
validated default profile enables three refinement passes per discovery program
and requests every provisional frontier label, including auxiliary-constant
lower sectors. Set `max_passes = 0` and `search_frontier_sectors = false` to
reproduce discovery without these coverage improvements.

The trigger is a concrete `NoApplicableRule` remainder reached while reducing a
requested target, frontier label, or derivative. With frontier-sector search
disabled, provisional refinement requires a nonzero auxiliary derivative;
final-audit remainders are checked even when constant. Only an encountered
native `UnprovedDescent` box containing that remainder is eligible. The frontend
chooses one finite nonfixed interval of bounded width and submits every integer
singleton face along that axis. All other axes retain their original bounds,
including unbounded symbolic axes and any fixed-zero storage tails.

The subdivision is atomic under the global added-domain budget. Learned faces
persist across closure rounds. The unsplit parents remain submitted, and their
native gaps remain inspectable. Successful face discovery is not recorded as a
proof that an entire parent domain is solved. All new rules must pass native
source replay; a closed result still requires one final program that reconstructs
every requested target and every derivative in its explicit stopping basis.
Conditions from all applications are retained. Newly requested frontier labels
force another bounded discovery round before this final audit, even when their
auxiliary derivatives vanish.

The default limits are:

| Option | Meaning | Default |
| --- | --- | --- |
| `guard_refinement.max_passes` | Additional discovery passes per provisional or final program | 3 |
| `guard_refinement.max_added_domains` | Distinct added faces across the entire preparation | 256 |
| `guard_refinement.max_interval_width` | Integer interval width, upper minus lower | 2 |
| `search_frontier_sectors` | Request all provisional labels for native discovery | true |
| `max_rounds` | Outer closure iterations | 12 |
| `discovery.max_depth` | Native source search depth | 3 |
| `discovery.max_domains` | Native discovery domain budget | 8192 |
| `max_frontier` | Provisional basis size | 256 |
| `max_requested` | Distinct requested integral labels | 4096 |

The low-level `GuardedDiscoveryOptions` constructor retains its smaller 2/256
limits; the table describes the closed-system service. The explicitly invoked
closure and full-flow tests expose the refinement values as
`RUSTFLOW_WEIGHTED_GUARD_PASSES`, `RUSTFLOW_WEIGHTED_GUARD_DOMAINS`, and
`RUSTFLOW_WEIGHTED_GUARD_WIDTH`, and frontier search as
`RUSTFLOW_WEIGHTED_FRONTIER_SECTORS`. Input records retain their values.
Checkpoints include `guard-refinements.json`, with exact parent bounds, native
failure, split axis, face values, and triggering index. The program files
retain the native source proofs and guards.

## Recorded validation

The opt-in `requested_index_rays` policy starts each positive index at its
requested value, each negative index at its requested upper bound, and fixes
zero indices to zero. Nonzero axes remain symbolic rays. Its boxes are subsets
of the usual sign sectors and are intersected with the same physical admitted
domain. The policy only directs native discovery; all subsequent remainders,
conditions, and final derivative/target replay obligations remain unchanged.
It defaults to false. The runtime graph harness accepts
`RUSTFLOW_WEIGHTED_REQUESTED_RAYS=true`; each round checkpoint records the
policy. Narrowed coverage is not a claim that the surrounding sign sector has
been solved.

The independent opt-in `prioritize_requested_indices` policy passes the current
requested labels to RustRed's native domain queue. Domains containing a hint
are visited before other domains; every split retains all children and all
unvisited pieces remain explicit `DomainBudget` gaps. Hints neither become
terminals nor alter identities, guards, ordering or exact replay. Empty hints
preserve the original FIFO traversal. The runtime graph harness accepts
`RUSTFLOW_WEIGHTED_PRIORITIZE_REQUESTED=true`, and each round checkpoint records
the policy. It defaults to false and can be combined with requested-index rays.

The initial 2026-10-09 standalone validation used `legacy-lorentz` sources,
polynomial completion powers, native depth 3/domain budget 8192, and refinement
limits 3/256/2, without requesting constant frontier sectors:

| Occupied cut | Spanning basis | Requested indices | Added faces | Closure time |
| --- | ---: | ---: | ---: | ---: |
| First charged line, N=7 | 7 | 24 | 22 | 2.884 s |
| Both charged lines, N=9 | 64 | 43 | 41 | 75.277 s |

Both closures passed their final exact audit. Native discovery gaps remained
outside the finite requested-target/derivative coverage claim. This historical
stage passed 91 selected native and ordinary regressions, recorded in
[`guard-refinement-regressions.json`](../reports/validation/2026-10-09-finite-density-native-assembly/guard-refinement-regressions.json).

Requesting every provisional frontier label lets the same native owner reduce
constant lower sectors as well. With this enabled, the first single cut closes
to six basis labels in 2.830 s. The double cut closes to eleven in 112.951 s,
with provisional sizes 66, 62 and 11, 143 requested labels, 118 additional
frontier requests, and 71 added faces. Its connection has 30 nonzero entries.
These are closed spanning sets; independence and minimality are not claimed.
Exact source, program and result evidence is retained in
[`frontier-sectors/sunset-N9/result.json`](../reports/validation/2026-10-09-finite-density-native-assembly/frontier-sectors/sunset-N9/result.json),
with 103 selected regressions recorded in
[`frontier-sector-regressions.json`](../reports/validation/2026-10-09-finite-density-native-assembly/frontier-sector-regressions.json).

The full massive-sunset fixed-dimension evaluation at epsilon=4/5 subsequently
passed four independent precision, series-order and starting-scale profiles,
including every occupied sector, the nonzero vacuum contribution and assembled
total. The separate reference comparison passed 40 component checks and 30
refinement checks; the maximum relative reference difference was 9.45e-36.
See
[`full-sunset-frontier-reference-comparison.json`](../reports/validation/2026-10-09-finite-density-native-assembly/full-sunset-frontier-reference-comparison.json).
The older 64-coordinate spanning system also passes at 40 guard digits; its
18-digit/20-guard run exposed a cancellation residual only 1.63 times the
strict endpoint zero threshold, not an exact physical divergence. No endpoint
tolerance was loosened. Laurent and four-loop evaluations are separate gates.
