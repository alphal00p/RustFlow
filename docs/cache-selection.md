# Physical boundary cache selection

`RustFlowCache` validates identity, epsilon range, coordinate precision and verified accuracy before pricing candidates. `TransportCost::lower_bound` is optional and returns `None` by default, preserving exhaustive evaluation for existing policies. A supplied bound must be finite, nonnegative and no greater than the rounded cost returned by that policy at the requested working precision. Bounds may describe routes that the full policy later rejects. Evaluated violations produce typed errors; skipped routes necessarily rely on the policy contract.

`ScaledDistance` computes the same MPFR squared distance for its bound and accepted cost, subtracting exact coordinates before evaluation. Physical path guards are checked only when a candidate can still improve or tie the accepted cost. Unknown-bound candidates are all evaluated. Pruning uses strictly greater bounds, so equal rounded costs still receive rational-complex exact comparison or higher-precision algebraic comparison, followed by the existing accuracy and insertion-order preference. No binary64 coordinates, distance screening, spatial graph or fixed nearest-neighbour cutoff is used. Cancellation is checked throughout physical cache selection and again before an exact hit returns.

The benchmark uses 33,000 exact points on a 220 by 150 positive-real grid and the analytic solution `Y(s,t,epsilon)=1` of a zero connection, with the physical guard `s*t != 0`. Seed accuracy follows from the exact formula. It requests 20 digits, stores three epsilon coefficients at 315 working bits and prices routes at 216 bits. This measures cache infrastructure, not integral evaluation or Monte Carlo throughput.

```sh
cargo run --release --example cache_bank_benchmark -- 33000 target/cache-bank-benchmark
```

One standalone CPU31 run of the implemented selection algorithm measured:

| Operation | Elapsed |
| --- | ---: |
| Cold indexed insertion | 0.545 s |
| Binary save (20,289,497 bytes) | 0.896 s |
| Binary load | 0.722 s |
| Exhaustive guarded nearest selection | 6.048 s |
| Lower-bound nearest selection | 0.357 s |
| Exhaustive guarded exact-hit selection | 5.910 s |
| Lower-bound exact-hit selection | 0.346 s |

Each selected point and MPFR cost matched exactly. The optimized queries checked one path instead of 33,000; all 33,000 compatible candidates still received cheap distance pricing. Peak process RSS was 159,744 KiB, including storage and allocator retention. These are single-run observations, not a statistical speedup guarantee. Full settings, numerical scope, source/binary hashes and raw measurements are preserved in [the validation report](../reports/validation/2026-10-04-cache-selection.json).

No exact/MPFR spatial index was found in the consumed Symbolica/HEPKit APIs or dependencies. Future indexing should follow measurements and preserve rigorous distance bounds, identity/range filters, route admissibility and the same tie rules. A generic custom transport policy need not assign an exact coordinate hit the cheapest cost.

Registered-root systems use this same selection and binary bank. Their point
keys include the local root germ, so opposite sheets at identical coordinates
remain distinct and exact hits cannot cross sheets. The regular-real admission
rules, shared precision checks and schema-4 compatibility policy are documented
in [Registered roots in the physical boundary cache](algebraic-cache.md).

A candidate that passes the initial evidence filter can still fail after its
uncertainty is propagated. RustFlow now tries compatible alternatives within
`max_boundary_attempts`, preserving this selector's cost bounds and tie rules.
See [Accuracy-aware cached boundary selection](accuracy-fallback.md).
