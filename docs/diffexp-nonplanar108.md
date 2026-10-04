# Full nonplanar 108-master transport

This benchmark uses all 108 masters and epsilon coefficients 0–4 from
[1812.11160v2](https://arxiv.org/abs/1812.11160v2), with the unchanged DiffExp
notebook's X0→X1 path. The normalization is `epsilon^4 exp(2 EulerGamma epsilon)`.
The boundary is supplied scientific data; this benchmark does not construct it
with IBP reduction or auxiliary-mass recursion.

X0 lies on three vanishing letters, W10, W12 and W23. Their residue is
`R = M10 + M12 + M23`, of rational rank 7. The Gram radicand is −3 there, so its
chosen square-root sheet is locally analytic. All five encoded boundary rows
satisfy `R c[k,0] = 0` exactly over rational coefficients with Pi and log(3) as
formal constants. Treating finite decimal mantissas as rational encodings does
not increase their recorded physical accuracy.

`AlgebraicSystem::analytic_origin` constructs the analytic, log-free sector. For
`Y' = epsilon A(x) Y` and `G(x)=x A(x)`, it solves

```
n c[k,n] = sum(l=0..n) G[l] c[k-1,n-l]
```

in increasing epsilon order at each positive coordinate order. It reuses the
existing rational square-root Taylor solver, native Symbolica series products,
and root-chart acceptance checks. Tail and midpoint/endpoint differential
defects provide consistency estimates. It rejects higher-order poles, origin
root branch points, non-pure-epsilon connections, and incompatible encoded
residues. Formal residue compatibility can conservatively reject a cancellation
that holds only after choosing a particular root sheet.

`AnalyticOriginOptions` uses the existing shared `CancellationToken`. Pass a
clone of `RunContext::cancellation` to cancel initialization and later transport
together. Checks occur between exact operations, compiles and recurrence work;
a single underlying Symbolica operation is not preempted.

The offset boundary is derived from this convergent local series. It is not a
copy of the constant X0 vector at a small nonzero coordinate. Original source
domain guards remain active for the subsequent ordinary transport.

The supplied X0 data have minimum recorded absolute Accuracy about 64.17 digits;
the X1 reference has about 49.85. The fixture preserves these metadata and uses a
conservative `1e-62` input floor. No propagated uncertainty certificate is claimed
for an arbitrary error box at the singular origin: such a box may contain
components in forbidden logarithmic sectors. Numerical acceptance instead uses
all 540 ancillary coefficients, a fresh original DiffExp run, and independent
precision/order/offset recomputation.

The isolated native prototype passed all 540 comparisons at 20 absolute digits.
The independent runs use 40 working digits/order48 at `x=1/128` and 60/order64
at `x=1/256`; their maximum absolute difference is `1.37e-36`. The refined
run differs from the ancillary endpoint by `2.88e-52` and the accurate live
original by `2.53e-57`. These observed residuals do not promote the recorded
reference accuracy. The integrated production scientific test also passed all
540 comparisons in both profiles (1060.56 seconds overall, 158.9 MB peak polled
process RSS). Its transport phases took 248.62 and 686.51 seconds.
[The original report](../reports/diffexp/nonplanar108-validation.json) retains the
isolated checks; [the integrated report](../reports/validation/2026-10-04-prescribed-integration.json)
records the separate production binary and source provenance.

Measured phases on one CPU, with concurrent host workloads:

- Native 40/48: 18.12s exact preparation, 32.39s local initialization, 233.89s
  ordinary transport, 292.04s total including the remaining compilation work.
- Native 60/64: 17.98s exact preparation, 39.81s local initialization, 684.63s
  ordinary transport, 750.04s total.
- Original accurate: 200.73s in public `TransportTo`, which includes its own
  singular-origin processing; 206.07s for the process including kernel setup.
- Original fast: 44.08s in `TransportTo`, 50.02s for the process; checked to
  14 absolute digits.

The native scientific test is explicitly ignored in ordinary suites because it
is a substantial numerical run:

```sh
cargo test --release --locked --test nonplanar108 -- --ignored --nocapture
```

The original oracle script takes `DIFFEXP_CHECKOUT`,
`DIFFEXP_NP108_ANCILLARY`, `DIFFEXP_OUTPUT` and `DIFFEXP_PROFILE` (`accurate` or
`fast`). It checks pinned source hashes and executes the original notebook's
alphabet definition without modifying the upstream checkout. Run it through a
licensed Wolfram kernel with external time/memory limits and one worker.

The accurate original profile uses WP150, order80, goal30, Mobius and Pade.
The notebook's fast profile uses order25 and an unspecified accuracy goal; its
roughly `5e-18` reference difference is checked to 14 absolute digits, not 20.
Different arithmetic precision, continuation strategies and initialization work
make a broad timing ratio inappropriate. The report preserves per-phase timings,
process limits, source hashes and actual errors separately.
