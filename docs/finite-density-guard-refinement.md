# Native finite-guard refinement

`WeightedClosureOptions.guard_refinement` optionally submits narrower integer
domains to the existing RustRed guarded source engine. It does not introduce
another elimination algorithm, change integral ordering, or supply a master
list. The default remains disabled (`max_passes = 0`) pending full numerical
validation of the newly closed physical system.

The trigger is a concrete `NoApplicableRule` remainder reached while reducing a
requested target or derivative. For provisional discovery the remainder must
have a nonzero auxiliary derivative; final-audit remainders are checked even
when constant. Only an encountered native `UnprovedDescent` box containing that
remainder is eligible. The frontend chooses one finite nonfixed interval of
bounded width and submits every integer singleton face along that axis. All
other axes retain their original bounds, including unbounded symbolic axes.

The subdivision is atomic under the global added-domain budget. Learned faces
persist across closure rounds. The unsplit parents remain submitted, and their
native gaps remain inspectable. Successful face discovery is not recorded as a
proof that an entire parent domain is solved. All new rules must pass native
source replay; a closed result still requires one final program that reconstructs
every requested target and every derivative in its explicit stopping basis.
Conditions from all applications are retained.

The limits are:

| Option | Meaning | Default |
| --- | --- | --- |
| `max_passes` | Additional discovery passes per provisional or final program | 0 |
| `max_added_domains` | Distinct added faces across the entire preparation | 256 |
| `max_interval_width` | Inclusive integer interval width, upper minus lower | 2 |

The explicitly invoked closure and full-flow tests expose these as
`RUSTFLOW_WEIGHTED_GUARD_PASSES`, `RUSTFLOW_WEIGHTED_GUARD_DOMAINS`, and
`RUSTFLOW_WEIGHTED_GUARD_WIDTH`. Input records retain their values. Checkpoints
include `guard-refinements.json`, with the exact parent bounds, native failure,
split axis, face values, and triggering index. The program files retain the
native source proofs and guards.

The 2026-10-09 standalone validation used the original `legacy-lorentz` sources,
polynomial completion powers, native depth 3/domain budget 8192, and refinement
limits 3/256/2. For the massive sunset example:

| Occupied cut | Spanning basis | Requested indices | Added faces | Closure time |
| --- | ---: | ---: | ---: | ---: |
| First charged line, N=7 | 7 | 24 | 22 | 2.884 s |
| Both charged lines, N=9 | 64 | 43 | 41 | 75.277 s |

Both closures took three outer rounds and passed the final exact audit. No
refinement budget was exhausted. N=7 retains 21 native discovery gaps and N=9
retains 47, outside the claim of finite requested-target/derivative coverage.
The bases are closed spanning sets; independence and minimality are not claimed.
The N=9 matrix has 95 nonzero entries.

The selected native and ordinary regression gates passed 91 tests. The physical
pilot explicitly tests the parent-plus-faces production path and decode/replay,
including retention of the unsolved parent and lower-sector residuals. These
results establish source closure, not full numerical amplitude agreement.
Numerical transport and independent-reference comparisons are separate gates.
See the exact source hashes, commands, resources, and saved programs in
[`guard-refinement-regressions.json`](../reports/validation/2026-10-09-finite-density-native-assembly/guard-refinement-regressions.json).
