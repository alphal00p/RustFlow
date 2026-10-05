# Native ggHg boundary runs

The `gg_hg_boundaries` example evaluates the complete planar or nonplanar
canonical vector using ordinary physical integrals, native AMF boundaries and
independently refined finite-epsilon fits. It requires neither Python nor
Mathematica. This is a long numerical acceptance runner; a successful exact
connection certificate alone does not establish successful boundary generation.

Build once, then run the exact Euclidean planar anchor:

```sh
cargo build --release --example gg_hg_boundaries
target/release/examples/gg_hg_boundaries \
  --family pl --digits 20 --guard-digits 40 --order 80 --workers 1 \
  --cache-directory /path/to/gg-hg-run \
  --cancel-file /path/to/gg-hg-run.cancel
```

Use `--family np` for the nonplanar anchor. The default points are respectively
`(-1/10,-1/25,-1/50)` and `(-1/10,-1/5,-1)` in unit vector-boson-mass coordinates,
with principal root germs, `D=4-2 epsilon`, and `+i0`. The default fit includes
powers zero through four. `--last-epsilon N` changes its final order; comparison
is restricted to the reference's available powers zero through four.

`--describe` prints the exact input, options, mathematical identity and build
fingerprints without loading numerical caches or starting a calculation.
`--list-configurations` lists the eight physical source labels for the selected
topology. A physical source run uses its exact published point and root germ:

```sh
target/release/examples/gg_hg_boundaries --family pl --list-configurations
target/release/examples/gg_hg_boundaries \
  --family pl --configuration W-Planar_EW1-1 \
  --cache-directory /path/to/gg-hg-run
```

Both topologies and all eight labels per topology are needed for the 16-source
physical seed gate. This runner generates boundaries; subsequent physical
transport and amplitude evaluation use the library or the community notebook.
The physical-source mode does not compare its values with an unrelated
Euclidean reference.

The directory layout separates types of state:

- `boundaries/`: the native binary `BoundaryCache`, saved atomically only after
  successful fitting and any applicable independent reference comparison.
- `samples/`: completed finite-epsilon samples, written by the library with
  strict source, point, precision, order and sample identities. A cancellation
  can leave reusable complete samples while producing no verified boundary.
- `exact/` and `ibp/`: exact reduction caches and native reducer checkpoints.
  `--checkpoint-directory PATH` overrides the latter.
- `runs/LABEL-TIMESTAMP-PID/`: an initial/final `result.json` and append-only
  `progress.jsonl`. Each invocation retains its own evidence.

Rerunning the same command reuses compatible boundaries and completed samples.
`--force` bypasses both numerical reuse paths while retaining exact reduction
work; the new numerical result and its evidence are retained in its run report.
The native bank keeps its strongest existing evidence when merging entries.
For a stability check, increase guard digits and order as well as using
`--force`, then compare the separate reports. Corrupt or incompatible caches
are reported as errors. They are never silently replaced with a fresh bank.
Run only one process per cache directory; use `--workers` for bounded sample
parallelism within that process.

Create the cancellation file to request cancellation, including before a run:

```sh
touch /path/to/gg-hg-run.cancel
```

The monitor polls every 250 ms; native algebra/reducer operations can defer
cancellation until their next cancellation check. A cancelled run exits with
status 130, saves a cancellation report, and does not save an unfinished
boundary. Remove that file before resuming. Other typed failures exit with
status 1. The runner prints the run-report path on stdout; progress uses stderr
and the JSONL log.

Every successful report contains canonical coefficients and absolute error
estimates as decimal strings, their basis and epsilon indices, working bits,
achieved mixed-scale digits, native provenance, build/cache identities and
elapsed times. Physical numbers never pass through binary64. Accuracy means
`10^(-digits) * max(1, |coefficient|)` and is supported by independent numerical
refinement, not rigorous interval bounds.

At a Euclidean anchor, reference loading occurs only after
`HiggsJetIntegralSystem::generate_boundary` has returned a complete verified
native result. The comparison-only loader reads the independently regenerated
40-digit records and requires each coefficient's complex difference to fit
inside the sum of native and recorded reference allowances. A failure is
retained in the report and prevents saving the verified bank. There is no
reference-boundary fallback. Planar reference evidence remains conditional on
its separately checked GPL constants, as explained in
[the anchor report](planar-anchor.md).

A degree-four anchor polynomial does not bound the unknown `O(epsilon^5)`
remainder of an individual sample. This runner compares fitted Laurent
coefficients, avoiding that truncation ambiguity. The
[planar Euclidean run](../reports/validation/2026-10-05-native-planar-boundary.json)
passes all 240 coefficients with 20 verified digits, including interrupted-sample
restart and an identical warm binary reload. Nonplanar and full physical-amplitude
acceptance remain separate gates. See the
[performance audit](native-boundary-performance-audit.md) for the scope of these
timings.

To compare a completed physical-source report with the independently recorded
plugin starting values, run:

```sh
python3 scripts/compare_gg_hg_native_source.py /path/to/result.json \
  --output /path/to/source-comparison.json
```

The comparison opens references only after checking native success, validates
the exact point, basis and root germs, and checks every coefficient using the
sum of recorded native and reference error allowances. Its default is 20 digits;
the reference's 24-digit evidence cap cannot be upgraded by its longer printed
mantissas. This script does not generate or inject numerical seeds.
