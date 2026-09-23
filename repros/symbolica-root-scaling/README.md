# Symbolica 3.0.0: coefficient scaling changes root-finder convergence

Standalone reproducer for the root-finding nonconvergence encountered while
implementing AMFlow. It depends only on Symbolica, with GMP integers and MPFR
floats. It does not need RustRed, AMFlow, Mathematica, or reference fixtures.

## Run

Install Rust with edition 2024 support and the native build prerequisites for
Symbolica's GMP/MPFR features. Configure any required Symbolica license through
your normal local environment; no license is included here.

From this directory:

```sh
cargo run --release --locked
```

The program calls `UnivariatePolynomial::roots(1000, tolerance)` with
`tolerance = 1e-50`. Coefficients and diagnostics use arbitrary precision
throughout, with no conversion through `f64`.

## Observed behavior

The polynomials `x^3 - 2` and `10^40 * (x^3 - 2)` have identical roots, but at
256-bit precision the first returns `Ok` and the second returns `Err`.
The returned estimates from both calls are accurate for the unscaled polynomial.
Increasing precision to 512 bits makes the scaled call succeed.

Recorded output is in [output.txt](output.txt). Tested on Linux x86_64,
rustc 1.97.1, with Symbolica pinned to 3.0.0 and transitive dependencies recorded
in `Cargo.lock`.

## Interpretation

This reproduces a scale-dependent convergence report, not evidence of incorrect
root values. Symbolica's `roots_hot_start` checks
`self.evaluate(x).norm_squared() < tolerance * tolerance` for every root.
Thus the tolerance is an absolute polynomial-residual tolerance; multiplying
the polynomial by a constant changes the stopping criterion. See the
[3.0.0 implementation](https://docs.rs/crate/symbolica/3.0.0/source/src/poly/univariate.rs)
and [API documentation](https://docs.rs/symbolica/3.0.0/symbolica/poly/univariate/struct.UnivariatePolynomial.html#method.roots).

At 256 bits, rounding leaves a residual near `1e-76` for `x^3 - 2`.
Multiplication by `1e40` pushes that residual above the requested `1e-50`.
The diagnostic reevaluates the returned estimates at 768 bits so the reported
residual is not just a cancellation artifact of the original working precision.
Residuals are not claimed as certified root-error bounds.

Whether the API should normalize its polynomial or document residual scaling
more explicitly is a separate question. This MRE does not establish a violation
of a specified root-location accuracy guarantee.

## Relation to AMFlow

The original failure involved a degree-12 exact denominator factor with large
integer coefficients. AMFlow now makes factors monic and scales their variable
before root finding, then validates the root set. Making the cubic monic is
enough for this smaller example; the original factor also needed variable
scaling at the chosen precision. The subsequent tensor-rank limit was in the
AMFlow implementation and is unrelated to this reproducer.
