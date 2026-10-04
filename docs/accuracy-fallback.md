# Accuracy-aware cached boundary selection

RustFlow first tries the cheapest admissible cached source that covers the
requested epsilon range and has sufficient recorded input evidence. That
initial evidence may still be insufficient after propagation along the path.
When an attempt returns an accuracy failure, RustFlow tries the next compatible
source under the same cost policy. Rational, registered-root and canonical
connections use this one selection and transport implementation.

`FlowOptions::max_boundary_attempts` limits attempted sources and defaults to 8.
It must be positive. Each attempt retains the existing independent precision
and order checks, weighted propagation of supplied errors and exact domain/germ
checks. Source uncertainty is never improved merely by increasing working
precision. A successful result includes `boundary_attempts`, recording each
actual source point, its MPFR cost, supplied verified digits, and accepted or
accuracy-rejected outcome. An exact cache hit records one accepted attempt and
performs no numerical transport.

Failed attempts do not insert intermediate points or change existing entries.
Only the final successful attempt commits its validated batch. Selection skips
already attempted entry indices before pricing; it preserves the existing cost
lower bounds, exact tie handling and source ordering. The cache is borrowed
mutably throughout, so those indices remain stable until the successful commit.
The full bank is neither cloned nor serialized during retry selection.

Cancellation is checked while selecting candidates, before transport and
between attempts. Invalid input, unsupported domains and numerical/backend
errors propagate immediately; only accuracy failures trigger another source.
No initial compatible source returns `IncompleteReduction`. If all compatible
sources fail accuracy, the result is `Accuracy`; reaching the configured attempt
budget returns `Limit`. Neither failure commits a partial attempt. On these
failures the existing typed error includes the number of attempted sources (or
the exhausted budget) and the last accuracy reason; the full attempt array is
returned on success.

The regression uses `Y'=i Y`: a nearby source with 20-digit evidence fails the
propagated uncertainty check, while a farther independently known source with
60-digit evidence reaches the 20-digit target and agrees with `exp(2i)`. Further
checks cover budget exhaustion, all-source failure, exact reuse after fallback
and cancellation between attempts without cache mutation.
