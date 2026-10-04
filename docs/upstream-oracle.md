# Live original AMFlow oracle

The pinned original AMFlow 2.0 package now runs on this host with Wolfram 15.0.1,
Kira 3.1 (`dad964cd79a39a9fca728bb2606f59537a95109c`) and Fermat 7.9b. Two automatic
sample calculations passed independent analytic checks at 20 digits on
2026-10-04. The [provenance report](../reports/validation/2026-10-04-live-upstream-oracles.json)
retains exact inputs, returned values and precision, source and output hashes,
commands, and log locations.

| Integral | Exact epsilon | Analytic reference | Original `SolveIntegrals` time |
|---|---:|---|---:|
| Massless bubble, `p² = -1`, powers `[1,1]` | `1/10` | `Gamma(eps) Gamma(1-eps)^2 / Gamma(2-2eps)` | 26.738822 s |
| Vacuum sunset, squared masses `[1,1,0]`, powers `[1,1,1]` | `3/4` | `Gamma(eps)^2 / ((1-eps)(1-2eps))` | 9.058479 s |

These are automatic correctness checks, not a matched performance comparison
with Rust. Both used one worker pinned to CPU 28, the original Mathematica DE
solver (`"MMA"`), no prior symbolic cache, 48 working decimal digits, expansion
order 96 plus extra order 50, and 24 requested output digits. Assertions required
20-digit scaled absolute agreement with analytic values evaluated at 80 digits.
Stored-number residuals do not establish accuracy beyond the returned 24-digit
precision. The normalization is `d^D l / (i pi^(D/2))` per loop, with
`D = 4 - 2 eps` and no `exp(EulerGamma eps)` factor.

Independent Rust evaluations with automatic recursive boundaries agree with
both live original results at 20 digits. Increasing Rust working precision from
60 to 80 digits and expansion order from 80 to 112 changed the bubble by
`8.82e-59` and the sunset by `1.71e-55` in absolute complex norm. Fresh boundary
providers were used for refinement. The sunset was prepared with symbolic
epsilon: preparing directly at `eps = 3/4` correctly rejected an ambiguous
dimensional indicial sector. The typed failure and successful symbolic
preparation, full numerical values, library hash and source provenance are
included in the report. These checks use different bases and arithmetic settings
from the original calculations, so their runtimes are not compared.

The bubble constructs a two-master auxiliary-mass differential system and its
terminal boundary automatically. The sunset is a genuine two-loop topology, but
Kira reduces this particular target to `-1/2` times the product of two tadpoles
at `D = 5/2`; its original evaluation does not require coupled sunset transport.
These checks do not constitute a live original calculation of all four paper
Laurent expansions. That Rust acceptance and its recorded upstream references
are documented separately in [coverage](coverage.md).

## Reproduction

Use an already licensed kernel launcher that also configures its child kernels.
On this host the launcher is `/home/ben/.local/bin/bern-wolfram`. AMFlow uses one
worker, but the parent and a sequential solver kernel may hold licenses at the
same time. Blade is a separate package; the bundled Kira adapter requires neither
Blade nor LiteRed. Its derivative helpers are included in the pinned source.

The runner verifies the hashes of `AMFlow.m`, `DESolver.m` and the Kira adapter,
copies the upstream tree into a new output directory, and changes only
`ibp_interface/Kira/install.m` to configure local executable paths. It bounds the
kernel process group to 600 seconds, followed by a 20-second termination grace
period, and retains the script, cache, logs, comparison and exit status.

```sh
source .dev-env
python scripts/upstream_oracle.py \
  --upstream target/upstream-performance/reproducible-build/amflow-26005517a288086c4cb4d1b26d829691bc088485 \
  --kernel /home/ben/.local/bin/bern-wolfram \
  --kira /nix/store/rlxy7f90ifbzr77a6hvnsc7a13vmzzay-kira-3.1-git-dad964c/bin/kira \
  --fermat /nix/store/0582y8gjxv1j01phy7i7l3z07ky72bry-fermat-7.9b/bin/fer64 \
  --output target/live-oracle-bubble --case bubble --mode sample --cpu 28
```

Use `--case sunset` and a different output directory for the second calculation.
Each output directory must be new. `--mode laurent` requests coefficients through
`eps^0` and checks their analytic values; that optional mode has not been run in
the recorded smoke checks. The original API's third `SolveIntegrals` argument is
the number of orders above `eps^(-2L)`, so the driver passes 2 for the one-loop
Laurent calculation and 4 for the two-loop one.

The two initial calculations used a target-local version of this driver and a
shell timeout wrapper. The reusable Python runner then passed the sunset sample
end to end (11.410213 seconds including startup). Its subsequent timeout cleanup
hardening passed a focused mock with a TERM-ignoring child; no additional
Wolfram calculation was needed. Failed harness attempts are retained in the report:
one path guard rejected Mathematica's trailing directory separator before any
calculation, and one assertion rejected an additional solved master returned by
Kira. Both harness issues were corrected without changing upstream algorithms.
