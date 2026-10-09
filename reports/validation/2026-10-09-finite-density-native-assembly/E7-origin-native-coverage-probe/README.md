# Historical E7 native coverage probe

These diagnostics explain the first residual in each **historical origin-only** E7 run. They do not evaluate an amplitude and do not use the new optional free-virtual zero domains. The historical source snapshot is `../massless-blocks-source-hashes.json`; complete source-bound programs are archived by `../E7-origin-pilot-archives.json`.

The saved single and double programs contain 1650 and 874 rules respectively. Neither program has a rule domain containing its first recorded residual. Each target is admitted by 20 original Lorentz-source guards. This excludes coefficient-pole or descent rejection as the cause of those two saved application failures. The nearest matching target patterns miss one guard axis: ordinary powers requiring at least 3 or 4, or a separate compact-loop rule requiring a raised cut. Virtual-loop rules themselves do not uniformly require raised cuts. Full axis counts and original proof seeds are preserved in `structural-inspection.json.gz`.

`native-domain-probe.rs` reconstructs the **complete historical native source context** from the locally generated program's embedded original source records, roles, parameter map, conditions, measure identity and zero domains. It uses the pinned transport record schema only to inspect this owned artifact. `GuardedProgram::decode_generated` then independently replays every stored rule against that reconstructed context. No identity, guard, condition, order or source row is relaxed. New discovery and application use only the existing RustRed guarded APIs.

| Target | Requested box | Domains | Rules | Uncovered discovery boxes | Applies to target |
|---|---|---:|---:|---:|---|
| single | exact point | 256 | 1 | 0 | yes |
| single | one symbolic ray | 256 | 2 | 0 | yes |
| single | ordinary zero faces fixed | 256 | 129 | 143 | yes |
| single | broad sign sector | 8192 | 2464 | 6380 | yes |
| double | exact point | 256 | 1 | 0 | yes |
| double | one symbolic ray | 256 | 4 | 0 | yes |
| double | ordinary zero faces fixed | 8192 | 3652 | 0 | yes |
| double | broad sign sector | 8192 | 1865 | 6251 | no |

All discovery used depth 3 and sample seed 0. The lower-budget broad controls are also in `summary.json`. Fixing ordinary zero faces retains positive-index and negative-numerator rays; native discovery may refine those rays further. The one-ray probe keeps index 6 symbolic and fixes the remaining coordinates.

These results locate concrete presentation and queue-coverage effects. They do **not** establish full closure: the double zero-face program still leaves 59 explicit residual terms after recursively applying its rules to the target. Exact-point success also does not justify replacing symbolic discovery with a separate finite-row elimination backend. No native correctness defect follows from this experiment.

The single target has physical slot 2 absent and slot 6 squared. Its virtual `q3` direction appears only through `(q3-q4)^2-eta`, with the `g1_2^2` numerator independent of `q3`. A translated virtual Euler source therefore admits the familiar tadpole relation, but the native full-corpus search selected another descending recurrence. The saved source-proof records identify the actual native result; no manually chosen tadpole recurrence was installed.

## Reproduction

Restore any detailed JSON/native program with `gzip -dk FILE.gz`. `archives.json` records both raw and archive SHA256 values. Runtime/resource records are measured **prebuilt child processes** and exclude compilation and Nix startup. The standalone compile was not timed. `compile-provenance.json` records compiler, source and dependency artifact identities.

Compile the diagnostic against the same resolved native stack (paths below are the artifacts used in this workspace):

```sh
nix develop --command rustc --edition=2024 -C opt-level=1 -C debuginfo=0 \
  native-domain-probe.rs -o /tmp/e7_native_domain_probe \
  --extern rustred=target/release/deps/librustred-42398b4d6555f617.rlib \
  --extern bincode=target/release/deps/libbincode-30fae4fd0bc9d1b5.rlib \
  --extern serde_json=target/release/deps/libserde_json-84255895f7dd0f90.rlib \
  -L dependency=target/release/deps \
  -L native=target/release/build/gmp-mpfr-sys-2e6d668657be102d/out/lib \
  -L native=target/release/build/gmp-mpfr-sys-c789de89ed67d3db/out/lib \
  -L native=target/release/build/blake3-222d768860627915/out \
  -L native=target/release/build/blake3-9252326cebd21167/out \
  -L native=target/release/build/psm-2e54985d752c3d12/out
```

The source path should be adjusted to this report directory when compiling from the repository root. Decompress the historical round-002 program and run:

```sh
PROBE_DOMAINS=8192 /tmp/e7_native_domain_probe \
  /tmp/e7-origin-double.bin double zero-face /tmp/double-zero-face.json
```

Case arguments are `single` or `double`; domain arguments are `point`, `one-ray`, `sector`, or `zero-face`. The diagnostic writes the complete native program to `OUTPUT.json.bin` alongside its detailed report. It uses no reference values.
