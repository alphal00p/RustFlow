# Symbolica root-scaling regression

Standalone check of the root-finding nonconvergence originally observed with
Symbolica 3.0.0. It now uses upstream `main` commit
[`75f8350094b90254ee71dc2a391fde0d14b0204a`](https://github.com/symbolica-dev/symbolica/commit/75f8350094b90254ee71dc2a391fde0d14b0204a)
(package version 3.0.1), which fixes the convergence check.
It depends only on Symbolica with GMP integers and MPFR floats, without RustRed,
AMFlow, Mathematica, or reference fixtures.

## Run

Install Rust 1.96 or newer and the native build prerequisites for Symbolica's
GMP/MPFR features. Configure any required Symbolica license through your normal
local environment; no license is included here. From this directory:

```sh
cargo run --release --locked
```

The program calls `UnivariatePolynomial::roots(1000, tolerance)` with
`tolerance = 1e-50`. All three cases must return `Ok`: `x^3 - 2` at 256 bits,
`10^40 * (x^3 - 2)` at 256 bits, and the scaled polynomial at 512 bits.
It also checks that three roots are returned and their unscaled residuals,
reevaluated at 768 bits, are below `1e-48`. All inputs and diagnostics use
arbitrary precision, without conversion through `f64`.

Current output is recorded in [output.txt](output.txt); historical output from
Symbolica 3.0.0 is preserved in [output-3.0.0.txt](output-3.0.0.txt).
Transitive dependencies are pinned in `Cargo.lock`.

## Original failure and upstream fix

With Symbolica 3.0.0, the polynomials `x^3 - 2` and `10^40 * (x^3 - 2)` had
identical roots, but at 256-bit precision the first returned `Ok` and the second
returned `Err`. The returned estimates had an unscaled residual near `1e-76`.
Increasing precision to 512 bits made the scaled call succeed.

The old stopping test used the absolute polynomial residual, so multiplying
coefficients by `1e40` made the requested tolerance unattainable at that working
precision. The
[updated implementation](https://github.com/symbolica-dev/symbolica/blob/75f8350094b90254ee71dc2a391fde0d14b0204a/src/poly/univariate.rs)
checks Newton and Aberth corrections relative to `max(1, |z|)`, together with
relative coefficient backward error and a roundoff allowance. It also restores
the working precision of each iterate after the update.

The diagnostic reevaluates returned estimates at higher precision to avoid
reporting only working-precision cancellation. Residuals are not certified
root-error bounds. The historical failure was a convergence report, not evidence
of incorrect root values.

## Relation to AMFlow

The original failure involved a degree-12 exact denominator factor with large
integer coefficients. AMFlow continues to make factors monic, scale their
variable, and validate the returned root set. Its degree-12 and extreme-scale
root regression tests remain in place. The separate tensor and boundary-budget
limits are in the AMFlow implementation and unrelated to this check.
