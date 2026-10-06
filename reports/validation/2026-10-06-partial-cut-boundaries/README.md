# Partial auxiliary-mass placements in cut families

The frozen native implementation passes 32 selected release checks and strict
all-target Clippy. `validation.json` records source/dependency fingerprints,
binary hashes, commands, numerical profiles and scope; `frozen-source-files.json`
was checked again after the tests. The generic library now permits explicit
nonempty subsets of uncut physical propagators. Python/CLI explicit-slot settings
are a separate interface milestone.

The nonfactorized three-loop example has massless cuts `r²`, `l²`, `(P-r-l)²`,
virtual denominators `k²-2`, `(k-r-l)²-3`, and `P²=1`. Only the first virtual line
is deformed in the partial calculation. The virtual bubble depends on the
varying Dalitz invariant `(r+l)²`. No numerical reference boundary is supplied.

Native All and partial placements agree with an independent double-beta series
to 20 relative digits at epsilon `1/13` and `1/17`. Increasing guard precision
from 50 to 70 digits and order from 80 to 112 also passes at `1/13`. The exact
140-term comparison series has a normalized tail below `3^(-140)/(1-1/3)`.
The two-loop massive virtual-bubble example independently exercises a nonzero
soft virtual child. Exact tests cover common hard virtual directions with a
soft difference, cut orientations, raised-cut normalization, and rejection of
dependent scalar-product poles with different affine constants.

| Three-body example, epsilon 1/13 | Native preparation | Native first evaluation | Original AMFlow SolveIntegrals |
|---|---:|---:|---:|
| One virtual line deformed | 0.431 s | 1.235 s | 146.248 s |
| Both virtual lines deformed | 0.592 s | 0.845 s | 125.876 s |

These are single runs with empty numerical/symbolic caches, excluding compilation
and runtime startup. Native calculations request 20 digits with 50 guard digits
and order 80. Original AMFlow requests and records 24 digits, using 48 working
digits, order 96 and 50 extra orders. Work and precision profiles differ; this is
a case comparison, not broad performance parity.

Both implementations use Lorentz-invariant positive-energy phase space times
`d^D k/(i*pi^(D/2))`, with no EulerGamma, flux or particle-symmetry factors. The
unchanged pinned original AMFlow `26005517a288086c4cb4d1b26d829691bc088485` agrees
with the comparison series at the asserted 20 scaled absolute digits. Its raw
residuals do not upgrade the recorded 24-digit reference precision. Metadata
retains every actual system configuration; root-system deformations are checked
as one-based positions `[4]` and `[4,5]`.

`initial-partial-oracle` preserves the initial successful numerical run and its
overbroad post-solve config collector. `partial-oracle` and `all-oracle` use the
corrected regular-file collector. Development compile/Clippy logs are retained
separately; final evidence is `unit-tests.log.gz`, `release-tests.log.gz` and
`clippy.log.gz`.

Remaining limits include dependent soft-pole decompositions, incomplete final
states, mixed virtual signs, unproved real phase-space pole bounds, massive
N-body leaves beyond two bodies, and linear cut-flow propagators. The native
region, tensor, reduction, ordinary recursion, terminal and Laurent-fit owners
are reused. No persistent numerical cache schema or dependency pin changed.
