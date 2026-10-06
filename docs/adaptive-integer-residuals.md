# Adaptive exact-source residual arithmetic

Ordinary rational and epsilon Taylor charts can opt into native Symbolica
integer polynomial products for their source defect. The default remains
`ResidualArithmetic::Ball`. Set `FlowOptions.residual_arithmetic` to
`ResidualArithmetic::AdaptiveInteger` in Rust, use
`EvaluationOptions(residual_arithmetic="adaptive_integer")` in Python, or set
`options.residual_arithmetic` to `"adaptive_integer"` in CLI JSON. The other
accepted string is `"ball"`; unknown values return a typed input error.

This changes arithmetic ownership
inside the existing whole-disk defect check, not the differential equation,
Taylor recurrence, conditioning tests, source domain, or error contract.
Registered-root charts retain their existing ball/MPFR owner. Optional Padé
candidates retain their separate rational-candidate defect calculation.

For the exact stored Taylor polynomial `Y`, each shifted source row is
`D_i Y'_i - sum_j A_ij Y_j`. Clear all real and imaginary Taylor coefficients
with one positive denominator `S` across every component and epsilon channel.
Clear `D_i` and every `A_ij` in that row with one positive denominator `T_i`.
The native paired integer polynomial products then represent
`S T_i (D_i Y'_i - sum_j A_ij Y_j)`. Division by `S T_i` occurs exactly before
native outward ball conversion. In an epsilon hierarchy, a source entry in
channel shift `r` multiplies solution channel `k-r`; omitted negative channels
remain zero. This includes the exact stored input precision, even if it exceeds
the current arithmetic precision.

The denominator enclosure remains that of the original shifted `D_i`, including
its native certified factors. Cancellation in a residual numerator cannot remove
a source hole. These are local differential-defect enclosures; boundary error,
parameter uncertainty, global amplification, and singular matching have separate
contracts. No new achieved-accuracy evidence is inferred from this arithmetic.

## Admission and fallback

The default internal selector admits at most 4096 scalar components, 1024 stored
Taylor rows, and a maximum input representation height of 4096 bits. It checks a
Float's stored precision and binary exponent **before** converting it to an exact
rational. Nonzero-center source shifts are limited to degree 24. Zero-center
shifts are identity copies and admit degree 1024, preserving sparse high-degree
source terms even for short solution polynomials.

A trial also has cumulative estimates of 100 million limb-weighted native
operations, 256 MiB of allocations/copies, and 131072 intermediate coefficient
bits. Preflights cover rational conversion, retained fallback inputs, denominator
LCMs, integer clearing, source shifts, derivatives, complex products, subtraction,
exact rational restoration, and final enclosures. Native operation outputs are
checked for observed coefficient growth. These are conservative admission
estimates, **not hard allocator quotas or interruptible native-call time limits**.
A resource rejection selects the established ball constructor. Limits affect
which arithmetic method is used, not which source degrees the library supports.

If an integer enclosure passes the supplied local error budget, the old ball
chart is never constructed. If it fails that budget or produces an inconclusive
numerical check, the same chart lazily constructs the old ball residual and
reuses it for subsequent trials. This preserves useful `InsufficientPrecision`
hints arising from arithmetic enclosure width; otherwise tighter exact arithmetic
could turn a precision problem into repeated step shrinking. There is no claim
that the two controllers make identical decisions: a genuinely tighter passing
enclosure may admit a step the old calculation rejected.

Only admitted integer charts retain cloned source rows and stored Taylor inputs
for this fallback; their estimated storage is charged before conversion. Clones
of the chart share one lazy owner. Its lifetime ends with the trial chart; saved
Taylor/Padé segments do not retain it. Cancellation is checked during bounded
operations, before waiting for the lazy-owner lock, after acquiring it, and
before/after the old construction and checks. A canceled construction is never
stored. Existing chart methods remain available; the shared controller uses
variants carrying options and context so cancellation reaches the optional
arithmetic path. Plain internal chart calls retain the default ball method.

Completed finite-epsilon sample checkpoint keys include this option through the
existing complete numerical-settings fingerprint. It does not change symbolic
reduction settings. Within a compatible build, physical boundary identities do
not depend on the arithmetic choice; reuse still requires their existing accuracy
contract. The existing `PORT_SOURCE_DIGEST` changes with the implementation and
rejects binary snapshots from an incompatible build. No cache schema is changed,
and transient residual charts are not serialized.

## Validation scope

The focused regressions compare exact integer results with an independent native
Gaussian-rational polynomial reference before outward conversion, and verify
that both old and new enclosures contain those exact coefficients. They cover
complex multichannel sources, unequal coefficient scales, higher stored precision,
source holes, selector limits, extreme exponents, cancellation, and same-chart
precision fallback. Test-only witnesses check whether fallback storage remains
empty or is reused; these are not public accuracy diagnostics.

An ignored benchmark runs the same complete transport controller with either
owner and records selected integer/ball chart counts and lazy fallback use, in
addition to constructor times. The paper-family benchmark uses a deterministic
supplied boundary and is not a measurement of automatic physical boundary
construction. Performance numbers must be read with the retained fixtures,
precision, step counts, outcomes, and independently checked transport agreement.

## Measured complete transports

The following medians use three alternating repetitions, order 80, and a
20-digit local target on the existing 12-master and paper-family 27-integral
rational systems. Times include Taylor construction and the complete continuation
controller; compiling the exact DE is excluded. The nearby runs start from their
recorded endpoint vectors: `0.1+0.2i -> 0.1001+0.2i` and
`-256i -> -255.9i`, respectively. These fixed boundary files do not provide an
independent input-uncertainty cap. Agreement and refinement statements are
therefore conditional on those stored vectors, not verified physical integral
accuracy. In particular, the paper-family vector is a deterministic supplied
boundary, not an automatically generated physical boundary.

| System | Bits | Long path, ball / adaptive | Nearby, ball / adaptive |
| --- | ---: | ---: | ---: |
| 12 masters | 201 | 3.827 / 3.545 s | 336.56 / 90.73 ms |
| 12 masters | 333 | 4.468 / 4.306 s | 392.35 / 128.60 ms |
| 27 integrals | 201 | 3.752 / 3.652 s | 928.66 / 150.36 ms |
| 27 integrals | 333 | 4.519 / 4.036 s | 1083.13 / 169.18 ms |

Long-path improvement is only 3–12% in these profiles: eight of the eleven
12-master charts construct their lazy ball fallback, and three of the four
27-integral charts use resource fallback. Nearby runs have one integer chart,
one accepted proposal, and no fallback, giving 3.05–6.40x improvement. Constructors
alone can improve much more; those microbenchmarks are not full-transport or
full-physics speedups. Separate fresh 201-bit/order-80 and 397-bit/order-112
runs compare both methods at 20 digits under the same fixed-boundary scope.

The shifted degree-64 and 512-bit degree-32 adverse sources select ball
arithmetic with essentially unchanged measured cost. The degree-64 case keeps
its typed minimum-244-bit precision request at 201 bits and succeeds at 397 bits.
Registered-root gg→hg transport uses a different owner and receives no speedup
claim from this experiment. Cache lookup, automatic boundary generation, and
full amplitude assembly are excluded from these timings.

Exact fixtures, raw timings, final-source hashes, validation logs, and reproduction
commands are retained in the [performance report](../reports/performance/2026-10-06-adaptive-integer-residual/README.md).
