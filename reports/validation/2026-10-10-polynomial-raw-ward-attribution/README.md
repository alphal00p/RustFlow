# Raw polynomial Ward attribution

The apparent reciprocal-domain improvement is completely reproduced within the original polynomial family. With the unrecentered raw Ward rows and four whole-boost polynomial shifts, all 19 physical RHS coefficient and condition maps exactly match the widened normalized-boost control: 143 terms versus 192 originally. Both positive-energy auxiliary probes remain explicitly `NoApplicableRule`. No reciprocal admission is needed for this finite result.

| Polynomial source presentation | Rows | Physical RHS terms |
|---|---:|---:|
| Original | 52 | 192 |
| Recentered raw Ward only | 54 | 177 |
| Shifted global boost only | 56 | 183 |
| Recentered raw Ward + shifted global boost | 58 | 168 |
| Unrecentered raw Ward only | 54 | 152 |
| Unrecentered raw Ward + shifted global boost | 58 | 143 |
| Same + shifted temporal U | 62 | 143 |

Both added-group orders agree in each combined comparison. The two 62-row orders provide the complete tested polynomial corpus; adding temporal U has no further effect on these first-hit reductions, although its independently valid rows remain available for other native searches. No active-closure pilot was run by this report; subsequent consumers must record their own bounds/results. The separate native source-selection experiment consumes the frozen corpus for matched point and closure probes.

The raw Ward source has the independent shell-cutoff proof and exact guard: one future occupied massless C1 shell, all energy-completion base indices zero, lower occupation zero, and original scalar/tail bounds. Before upper multiplication the rows are `(D−2) I0 − E0 I1 = 0` and `(D−2) Is + s E0 I(s+1) = 0`. The energy factor is a polynomial numerator image shift, not a widened base domain. Unrecentered and H-recentered rows are mathematically equivalent; the native first-hit search reaches the simpler bulk combination earlier with unrecentered rows.

The four shifted global-boost rows and four shifted temporal-U rows use the exact original-domain transformation: substitute `a6 → a6+1` in every coefficient and declared condition, shift every integral label by +1 in slot 6, and pull back the domain to `a6 <= −1`. All 52 original sources, zero domains, physical targets, phase coefficients and deformation remain. No unshifted raw boost is added to the 62-row corpus.

`source-restriction-comparison.json` checks the normalized reciprocal source directly against the polynomial raw Ward row. The raw guard is an exact subset. After coefficient specialization and the existing cut/lower-contact zero projection, coefficients and label shifts match up to the fixed nonzero scale 8. Three extra normalized-source conditions are the retained nonzero constants 2 from denominator clearing; there are no additional symbolic conditions. An initial checker assertion incorrectly assumed these condition lists were empty; its source and explicit failed assertion are retained before the corrected constant-condition verification. This algebraic comparison does not derive a physical source by using a forbidden reciprocal seed: the polynomial source has its own independent shell proof.

Every small arm constructs fresh native rules at depth 3 with one exact-point domain per label, encodes/decodes the program, and checks exact applied terms and conditions. There are 19 native rules for 21 requested labels: the two positive-E labels remain unresolved in all polynomial arms. Eight source-only guard tests verify the three surface recurrences and rejection of raised cut, excessive E0 numerator, other energy numerators and lower-contact cases. The source-only test deliberately omits the original zero boxes; complete controls retain them. Native replay establishes consequences of supplied sources; the physical derivations and independent review are separately hash-bound.

The frozen preferred corpus is `native-source-probe/strongest62-boost-ward/program.bin`, SHA256 `f482d60e6aee6f3187be2eac039b26fec8fb8cc7cef185d577b5d2f8ddad14a2`, 58,458 bytes. It remains raw and must not move while downstream probes consume it. Consumers should reconstruct its sources and perform fresh native search; its 19 small-probe rules are not a closed basis. `strongest62-summary.json` records both source orders and their exact hashes.

Raw result `additional_sources` includes the intermediate polynomial boost generator rows before index shifting; these are not all added directly. Actual additions are selected by the mode and recorded shifted-source groups. The unchanged production-order native library is SHA `ac85472324303a84e909023a9aac6fe77f4993934985519fd88275cae0e1aad0`. Standalone wrappers use opt-level 0 with release dependencies. Build and execution scopes are separated. `source-revisions.json` maps every historical source hash to its preserved version, and per-run provenance binds each binary. Production-owner drafts are excluded from the frozen experiment manifest.
