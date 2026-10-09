# Independent exact prism virtual coefficient identities

All five proposed generic-epsilon virtual coefficient identities now pass an
independent native exact proof. This is **not a complete finite-density prism
reference**, a thermal-contour proof, or an AMF prediction comparison.

The initial optional AMFlow/Kira symbolic reduction took 24.674871 seconds.
Its exact inputs, candidate coefficients, generic master expressions, runtime
metadata and upstream hashes remain under `candidate-origin/`. The original
candidate JSON still says that native verification is pending: that is the
unaltered historical proposal, not the final status recorded here. No external
solver is required to verify the resulting rational identities.

The independent verifier reconstructs the original five numerator/mass-jet
polynomials, lowers them exactly, and checks all 25 distinct original/master
integrals in a common RustRed exact reduction. It uses native source replay and
verified routing symmetries with analytic bubble shortcuts disabled. Every
coefficient of each of the five final differences is exactly zero.
The successful bounded process took **0.165157 seconds, RSS18492 KiB**;
compilation and library hashing are excluded. `native-replay/` contains its
exact reduction table, nonzero conditions, checkpoints, log, executable/source
hashes and result. The first run with symmetries disabled left only the
K↔L sunrise routing difference and is retained as an inconclusive diagnostic
under `native-replay-without-symmetry/`.

The ordinary native family is

```
D0=K^2, D1=L^2, D2=(K-p)^2, D3=(L-p)^2,
D4=(K-q)^2, D5=(K-L)^2, D6=(L-q)^2,
p^2=-1, q^2=0, p.q=-1/2, D=4-2*eps.
```

D6 is a numerator completion, with zero index in all physical targets. The
five original polynomials and their fixed-original-numerator off-shell jets
are derived in [the construction document](../../../../docs/finite-density-prism-reference-construction.md).
Their rational coefficients are in the
[validation-only fixture](../../../../fixtures/finite_density/prism_virtual_candidate_relations.json).
The [native test](../../../../tests/finite_density_prism_virtual_replay.rs) uses
no supplied oracle data, numerical reference values, or feature evaluator.
The persisted native proof is bound to the source283 snapshot; the integrated
test differs only in portable fixture/report paths and test harness behavior.
Its coordinated release build now passes both tests in 0.315822 seconds
(RSS18472 KiB); `native-cargo-replay/` records the fresh executable/source hashes
and exact proof table, bound to the free-virtual source snapshot.

The three remaining ordinary masters have independent sequential-integration
Gamma representations in the native `d^Dk/(i*pi^(D/2))` normalization. Set

```
B = Gamma(eps)*Gamma(1-eps)^2/Gamma(2-2*eps),
M1 = I[0,1,0,1,1,1,0]
   = B*Gamma(2*eps)*Gamma(1-2*eps)^2/Gamma(2-3*eps),
M2 = I[0,1,1,0,0,1,0]
   = -Gamma(1-eps)^3*Gamma(2*eps-1)/Gamma(3-3*eps),
M3 = I[1,1,1,1,0,0,0] = B^2.
```

For M1, first integrate the K bubble and then the remaining massless triangle
with two on-shell external legs. M2 is the sequential massless two-point
convolution; its odd denominator normalization supplies the displayed minus.
M3 factorizes into two ordinary virtual bubbles. The previous independent
finite-epsilon Gamma check remains in the
[partial numerical report](../upstream-prism-virtual-attempt-2/README.md).
The native proof establishes the rational reduction at generic epsilon,
retaining raw candidate denominators `eps`, `1+eps`, `1+2*eps` and every native
nonzero condition. Laurent use is by meromorphic continuation, not by replacing
these singular factors with values at epsilon=0.

`candidate-origin/formal-two-cut-combination.wl` combines the five verified
virtual coefficients with the separately derived compact Beta moments.
`formal-two-cut-msbar-series.wl` is a **formal two-cut-sector expression only**:
it starts at epsilon^-6 and lacks the three-cut amplitude and its required
cancellations. It must not be used as a complete prism reference. The common
regulated off-shell mass-jet continuation for the coupled two-virtual sectors,
the complete three-cut compact/simplex subtraction forest, complete cut
assembly and refined Laurent integration are still missing. The
[high-D triangle bound](three-cut-high-D-bound.md) resolves only a sufficient
compact domination domain after virtual UV continuation for the three-cut
one-virtual-loop sector.

To rerun the native proof after the coordinated build:

```sh
RUSTFLOW_PRISM_VIRTUAL_REPLAY_REPORT=/tmp/prism-native-proof \
  nix develop --command cargo test --locked --release \
  --test finite_density_prism_virtual_replay -- --nocapture
```

No Mathematica, AMFlow, Kira or Fermat installation is used by this command.

## Safe-pause handoff

The user requested a safe pause after the completed native proof. The restricted
coupled-virtual endpoint derivation and exact sector/audit scripts are preserved
in `joint-eta-certificate/`. All 3,008 maps, positive units, Jacobians, complete
branching/volume coverage and the conservative D>18 rank-six bound received an
independent exact audit. A subsequent independent review accepted the localized
UV-Taylor/Mellin Euler-derivative remainder argument; the saved draft predates
that final prose update. This applies only to the original full six-positive-
slot virtual integrands. Generated active supports require fresh certificates.
The compact/thermal completion write-up and complete refined three-cut Laurent
reference remain unfinished. No complete prism reference was generated.
