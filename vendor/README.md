# Vendored RustRed

`rustred/` is an editable copy of `../rustred_fermi` at commit
`78969aab524b7d6a2eec36f59b01e9af1e04cc04`, originally from
[`alphal00p/rustred`](https://github.com/alphal00p/rustred). The source checkout was
clean when copied. RustFlow's root manifest selects
`rustred/crates/rustred-core` as a path dependency; its `rustred-order` dependency
also resolves within this directory. Edit these files directly and commit the
changes in the RustFlow repository. No sibling checkout or Git submodule is
needed to build RustFlow.

The initial import contains all 3,071 tracked regular files and symlinks,
byte-checked against the source checkout, including the upstream workspace,
tests, examples, documentation, and [MIT license](rustred/LICENSE). Upstream
files are unmodified. Git metadata and untracked files, including build outputs,
were not copied. The only omitted tracked entry is the uninitialized
`vendor/symbolica` submodule gitlink. RustFlow already supplies Symbolica,
Numerica, and Graphica through its root manifest and lock.

Run builds and RustFlow tests from the **RustFlow repository root**, for example:

```sh
cargo check --locked
cargo test --locked --test source_fingerprint
cargo test --locked --release -p rustred --test guarded_source_api
```

Use `nix develop --command` before each command when using the repository's
Nix development environment. The root selects the shared Symbolica dependency
and disables RustRed's experimental `reconstruction` feature. RustRed's retained
standalone workspace patch table refers to its omitted Symbolica submodule, so
running Cargo with `vendor/rustred/Cargo.toml` as the root is a separate upstream
build setup, not the RustFlow build workflow.

RustFlow's existing source fingerprint hashes both resolved RustRed packages,
their manifests, source files, and enabled features. Moving from Git to this
path dependency changes the cache identity, and later edits here invalidate it
without a version bump. Historical numerical results retain their original
dependency provenance.

If an embedding host also uses RustRed, align all of its RustRed dependencies
with this copy as described in
[dependency embedding](../docs/dependency-embedding.md). To update the import,
review changes against the upstream revision above, preserve any local edits,
and record the new baseline here.

Import validation passed with the unchanged shared dependency pins:

```sh
cargo test --locked --offline --release -p rustred --test guarded_source_api
cargo check --locked --offline --release --features python --all-targets
```

Local changes after the initial import:

- Guarded recursive reduction combines exact coefficients of repeated
  `NoApplicableRule` residuals, including cancellations reached through
  different paths. It retains all nonzero conditions and preserves other
  failure records such as `WorkLimit`. Native replayed-diamond regressions
  cover cancellation, a nonzero sum, and exhaustion of the application budget.
- Guarded domain discovery exposes `solve_domains_with_priority_points` to
  visit pieces containing requested labels before unrelated pieces. Hints
  affect traversal only: source guards, native search, sealing and exact replay
  are unchanged, and all unvisited boxes remain explicit `DomainBudget`
  results. Calling the existing solver without hints retains its FIFO order.
- Guarded programs expose bounded `union_replayed` to keep earlier proved
  rules as an ordered fallback to fresh rules. The complete original source
  binding must agree; duplicate rules count against the supplied budget, and
  only the caller's current terminal set is retained. Native replay and the
  common coordinate order are checked again for the combined program. The
  operation does not certify unresolved discovery domains or family closure.
- Guarded programs also expose bounded `with_terminals_replayed` to replace
  a provisional stopping set before residual discovery. It replays the existing
  rules without reordering them, validates the supplied terminals, and leaves
  any newly exposed uncovered integrals explicitly unresolved.
- Guarded discovery also exposes `solve_domains_partitioned`: each requested
  box receives a bounded allocation from one total budget. It concatenates
  native rules and gaps without changing the search or replay contracts;
  unsubmitted boxes remain explicit `DomainBudget` gaps. The caller replays
  the combined rules when constructing its program.
- Guarded recursive application visits pending numeric labels in the common
  replay-validated integral order, hardest first. Descending paths therefore
  combine before a shared child is expanded, avoiding repeated applications
  and allowing exact cancellation before an unused rule acquires conditions.
  Conditions from every actually applied rule remain retained; unsupported
  labels and exhausted work budgets remain explicit. This changes scheduling,
  not rule selection, source identities, or exact replay.
- Guarded programs keep their source, rule and terminal state private and expose
  bounded `union_verified` and `with_terminals_verified` for immutable programs
  whose proofs were already replayed. Composition checks the complete source
  identity, common coordinate order, replacement terminals and duplicate-inclusive
  rule budget; it retains exact rule precedence and conditions. New construction,
  cold decoding, the existing replayed methods and RustFlow's final complete
  closure audit still replay proofs. RustFlow uses the verified methods only for
  internal scheduling, avoiding repeated replay of unchanged intermediate rules.
  Native and wrapper regressions compare exact persisted bytes and reductions,
  reject incompatible sources/storage/terminals, and exercise cold decoding.
