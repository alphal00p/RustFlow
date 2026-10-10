# Optional PolynomialClosure policy draft

This is a separate, **unapplied and uncompiled draft**. The earlier raw-only
`implementation-draft/` remains byte-for-byte unchanged. Rustfmt parsed the
staged files, and `git apply --check` checks applicability only. No source or
numerical acceptance follows from those checks.

The requested enum variant is `WeightedSourcePolicy::PolynomialClosure`, with
the explicit serialized/string value `polynomial-closure-v1`. Existing options
already pass this enum through the source factory, sealed endpoint/origin
evidence, reports and source/cache identity; no parallel flow implementation is
introduced. Legacy remains the default. This policy rejects the inverse-energy
option, preserving its original polynomial contract.

## Order and exact owners

The factory assembles these groups, in order:

1. The entire unchanged `lorentz_ibps_legacy` corpus.
2. Temporal U rows for each matched compact energy completion, using the same
   full occupation partition as the original legacy temporal rows, then shifted
   by +1 on that completion.
3. Polynomial common-boost rows for each such completion, then the same +1
   shift. The actual action is zero on all Gram coordinates and
   `delta E_i = E_i E_a - g_ia`, with divergence `(D-1) E_a`. Every physical
   factor's action is checked to vanish exactly before returning sources.
4. The raw, unrecentered C1 Ward pair, **only** when the bound sealed class is
   `SingletonGerm`. Other bound classes deliberately omit this pair and retain
   the valid polynomial groups. Absence of bound origin evidence also omits
   this additional theorem; it does not infer a new prescription.

The singleton three-loop fixture has 52 + 4 + 4 + 2 = 62 rows. Those numbers
appear only in a focused regression, never in production selection. There is
no filtering by source ID, graph name, cut slot, or fixture coordinate slot.

The new private `PolynomialEnergyCompletion` means **only** the exact ordinary
factor identity `F_slot = c E_a` with nonzero rational c and occupied a. It is
not the existing positive-mass nonsingularity certificate. It grants no
inverse-energy integral or rational-vector-field source. Virtual energies,
mixed energy sums, constant offsets and nonordinary factors receive no such
map. When a coordinate basis has no direct cE completion, that compact loop
gets no shifted rows; metadata retains the exact matched list. The raw Ward
pair can still use the full inverse basis to convert its polynomial energy.

## Reindexing and conditions

Before shifting a generated identity, intersect its original domain with the
box where **all** ordinary completions are nonpositive. The private helper then
applies old `a = new a + shift` simultaneously to every coefficient, every
nonzero condition, every integral shift and the complete domain. In particular
`a_energy + 1 <= 0` gives the new base guard `a_energy <= -1`. The theorem is
evaluated at an original polynomial source seed; no forbidden reciprocal seed
is imported.

The helper preserves unsimplified coefficient and condition expressions so
the existing guarded adapter can retain raw denominator requirements. It uses
checked integral-shift arithmetic and native exact domain pullback. Every
backend storage tail remains fixed zero; support-zero boxes are never widened.
The existing final admission intersection remains in place.

Domain budgets apply per generated vector field. A checked aggregate row bound
is `2 * matched_completions * domain_budget`. Affine inverse factors and fields
of degree at most two give the checked monomial bound
`(physical_arity+1) * (coordinate_count+1)^2` per source. Exceeding a bound
returns an error without partial source data. These bounds supplement, not
replace, later native matrix/search budgets.

## Pending tests and review

The patch retains the raw-only draft's private proof tests, adapts its factory test to the actual policy path, and adds:

* A complete two-axis identity shift with nontrivial coefficient/condition
  substitutions, negative-domain pullback, retained raw pole, and integral and
  integer-domain overflow failures.
* Full 62-row order for the actual input, identical legacy prefix, unchanged
  zero domains, original polynomial bounds, native identity separation,
  source-option mismatch rejection and inverse-option rejection.
* A two-compact class that retains shifted temporal/common-boost rows while
  omitting the raw singleton theorem explicitly.
* Nonunit cE matching and rejection of a mixed-energy completion after an actual
  routing change.
* Explicit shifted-energy base -1 specialization (old energy power zero), with
  C2 and all other energy powers zero: every formal positive completion image
  must have zero coefficient before any required-cut or support-zero projection.
  Base zero is rejected; base -2 also retains polynomial images.

All tests remain unexecuted. The new policy must be compared against the frozen
62-source experiment by native exact rows, guards, conditions and context
replay, allowing only nonzero rational normalization when the adapter clears
denominators differently. The actual original target sums and their weighted
derivatives still need the final native closure audit, followed by the shared
boundary/transport/precision/Laurent gates. The current timed-out but improved
frontier history is evidence for continued investigation, not completed closure.

This revised coherent draft exposes only the existing source factory/options
path and the new enum policy. It removes the raw-only public factory and its
boolean branch; raw-only tests use the private generator, while source/cache
admission tests exercise the real policy factory. The previously reviewed
revision1 patch and source map are preserved under `revisions/revision1/`.
