# Partial independent prism virtual reference

One bounded original AMFlow run completed a requested 18-digit profile for five
ordinary two-loop virtual coefficient functions at epsilon=1/5. Its independent
26-digit run did not finish before the 600-second limit (620.099887 seconds
including termination grace). No complete finite-density prism amplitude or
Laurent reference was generated. The earlier unreduced attempt also timed out
without a complete profile and remains in the adjacent attempt-1 directory.

The exact numerator/raised target combinations in `exact-targets.wl` were
reduced together with upstream Kira before numerical AMFlow evaluation. The
surviving masters, the full exact reduction and all five combinations are saved
in `reduced-combinations.wl` and `exact-combination-reduction.wl`. No supplied
oracle or RustFlow prediction was read. This optional external stack is used
only for reference validation and is not a feature dependency.

The common native family has denominators
`K²,L²,(K-p)²,(L-p)²,(K-q)²,(K-L)²`, with zero-power completion `(L-q)²`,
`p²=-1`, `q²=0`, `p.q=-1/2`. The first four coefficients are the simple and
off-shell derivative terms for occupied pairs `[1,5]` and `[1,7]`; `CC` is the
raised uncut original line for pair `[5,7]`. The exact maps and fixed-original-
numerator differentiation are documented in
[the subtraction construction](../../../../docs/finite-density-prism-reference-construction.md).
The initial `definition.json` prose describes only pair A; the exact saved
target list and coefficient labels define all five computations.

With `e=epsilon` and `B=Gamma(e) Gamma(1-e)^2/Gamma(2-2e)`, the independent
ordinary integral derivation gives

```text
M1 = B Gamma(2e) Gamma(1-2e)^2/Gamma(2-3e),
M2 = -Gamma(1-e)^3 Gamma(2e-1)/Gamma(3-3e),
M3 = B².
```

`M1` follows by integrating the K bubble and then the one-scale triangle of
indices `(1,1,e)`. `M2` is the three-denominator massless sunset with its native
Minkowski sign; `M3` is the product of two bubbles. These formulas were derived
independently of the upstream numerical answers. The separate `gamma-check.wl`
compares 60- and 90-digit Gamma evaluation and then the saved 18-digit upstream
profile. All three master comparisons and all five virtual coefficient
comparisons pass 1e-12 relative agreement. The largest observed relative
discrepancies are 8.314104e-26 for masters and 3.8960075e-25 for combinations.
Those discrepancies are observations, not rigorous uncertainty intervals or a
claim that AMFlow returned 25 reliable digits; its returned combination
precision labels range from approximately 16.7 to 17.3 digits.

The successful independent Gamma check ran in 7.750001 seconds. Its initial
run also computed matching numbers but emitted a precedence error in the extra
thread assertion. That log/script remain as `initial-gamma-*`; the corrected
run explicitly asserts one thread and finishes without that message. Both
parent and child Wolfram launches use the saved mandatory single-thread wrapper,
one CPU affinity, and a task-owned temporary license file removed on exit. No
other user's Wolfram configuration was accessed.

`profiles.json` preserves the actual virtual coefficients in native and raw
Euclidean normalization. The first four multiply `(4 pi)^(-D)`; `CC` has an
additional minus sign from its odd total virtual denominator power. These
virtual values cannot be inserted as production boundary conditions. Complete
compact assembly, regulated endpoint subtractions, common thermal continuation
and actual Laurent precision refinements remain missing.
