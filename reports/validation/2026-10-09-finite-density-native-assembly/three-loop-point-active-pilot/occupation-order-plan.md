# Isolated native occupation-order experiment

The larger canonical pilot does not close: original-output reductions alone grow54 →210 →647 →1188 leaves. Only0,1,4,8actual basis derivatives are processed. This supports investigating the native recurrence direction, but does not prove that ordering is the obstruction.

## Exact witness

`ordering-witness.json` decodes and natively replays the complete1035-rule program, examines921requested labels, and selects a minimal source-count witness among254occupation-raising RHS edges. Rule525 uses one original `lorentz/0/0/occupation-case-0` source and no nonzero condition. On its saved fixed-index domain,

```
A = B + (2 epsilon - 2) C
A = I[1,1,-1,2,1,0,-1,-1,0,0,0,0]
B = I[1,1, 0,1,1,0,-2,-1,0,1,0,0]
C = I[1,1, 0,1,1,0,-1,-1,0,0,0,0]
```

Current native comparison makes bulk A harder than surface B. An occupation-degree-first order would orient this row as `B = A + (2 - 2 epsilon) C`, with both RHS occupation degrees0. Its candidate discovery domain is the exact B point; the saved source seed has occupation0 and remains the same physical identity. That reversed rule has **not** yet been discovered or sealed under an alternate order. No rule is installed by hand.

## Native invariants and call sites

- `solver/index.rs::IntegralOrder::compare` is the common comparison for row sorting, pivot selection, candidate canonicalization, descent replay, concrete application, and the native pending queue. A guarded-only prefix there reaches each owner through the same existing order object.
- `solver/search.rs::SectorSolver::new_guarded` constructs that order with the exact source roles. `with_roles` forbids combining those roles with the legacy compiled-order program. A frontend coordinate permutation only changes physical ties and cannot move occupation ahead of them.
- `guarded/replay.rs::replay_candidate` validates exact roles and arity, re-instantiates original rows, sorts columns with the same order, and checks the candidate. `prove_descent` requires nonnegative occupation images, stable ordinary signs, and a strictly lower RHS under that same comparison. Those checks remain unchanged. Stable-sign checks may be conservative after reordering; they must not be weakened for this experiment.
- `guarded/lifecycle.rs::GuardedProgram::new` replays every rule and enforces a common coordinate permutation with no custom program. Its harder-first pending queue borrows the replayed order. Concrete rule application again checks descent. Source identities, zero domains, nonzero conditions, coefficient poles, and role admission are unchanged.
- A nonnegative total occupation grade followed by the existing well-founded physical order is well founded. In symbolic rules, the common occupation base cancels in comparisons; replay separately checks admissibility over the whole domain. With no occupation roles, the prefix is equal and ordinary ordering remains unchanged.
- No old rule object or proof container may enter the alternate library. Reconstruct only the exact source data, discover fresh candidates, replay and encode them under a unique experimental schema, then test reload under that library. The ordinary v2 loader must reject the experimental schema, and the experimental loader must reject v1/v2 proofs. A unique rustc metadata value prevents accidental linking as the production crate.

## Two distinct proposed controls

Parent is first building the narrower **total occupation degree before physical order** policy. Ties in total degree still use the entire existing physical order before occupation lex. Its tag is `rustred.guarded-source-program.experiment-occupation-degree-first.v1`; sources/builds live outside production under the dedicated occupation-order report.

`experimental-occupation-first.patch` is a separate, unapplied stronger control: total occupation degree **and occupation lex** before physical components. It uses another unique tag. It must not be confused with the first build. Both orient the displayed H0/H1 witness the same way, so the narrower control should be tested first.

## Standalone build and probe plan

`experimental-build-inputs.json` resolves the exact six extern rlibs by the native Cargo dependency fingerprints and records their hashes. The only generated include is `OUT_DIR/runtime_arities.rs`; the matching production directory is `target/release/build/rustred-790d7faed7c42925/out`. Copy that generated file, hash it, and use an isolated crate-source copy. Suggested environment is `CARGO_CRATE_NAME=rustred`, `CARGO_PKG_NAME=rustred`, `CARGO_PKG_VERSION=0.1.0`, isolated `CARGO_MANIFEST_DIR` and `OUT_DIR`. Use edition2024, native feature only, explicit externs and dependency/native-library search paths, a unique crate metadata suffix, and an output rlib outside Cargo’s target directory. No production Cargo graph or source is modified.

First rediscover/replay the fixed B-point witness with the original physical source corpus. Then compile the same corrected-history pilot against the experimental rlib and run matched bounded controls. Each run reconstructs sources only, retains strict non-NoRule failures and conditions, and uses exact weighted original targets and derivatives. Only a final common-program stopping-set reconstruction of all original targets and every final basis derivative can establish closure. Existing original-order unit expectations are not evidence for the alternate order; test its own invariant, replay, and schema-isolation obligations explicitly.

The pilot wrappers were compiled with `opt-level=0` against release dependencies. Their timings are diagnostic process costs and must not be compared as optimized production performance. Matched old/new wrappers use the same flags.
