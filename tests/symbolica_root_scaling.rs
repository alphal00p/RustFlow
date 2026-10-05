use std::sync::Arc;
use symbolica::domains::float::FloatField;
use symbolica::poly::univariate::UnivariatePolynomial;
use symbolica::prelude::*;

#[test]
fn root_convergence_is_independent_of_coefficient_scale() {
    // Multiplication by a nonzero constant leaves the roots unchanged.
    for (bits, scale_text) in [(256, "1"), (256, "1e40"), (512, "1e40")] {
        let f = |n| Float::with_val(bits, n);
        let scale = Float::parse(scale_text, Some(bits)).unwrap();
        let tolerance = Float::parse("1e-50", Some(bits)).unwrap();
        let field = FloatField::from_rep(Complex::new(f(0), f(0)));
        let polynomial = UnivariatePolynomial::from_coefficients(
            &field,
            [-2, 0, 0, 1]
                .into_iter()
                .map(|n| Complex::new(f(n) * &scale, f(0)))
                .collect(),
            Arc::new(PolyVariable::Symbol(symbol!("x"))),
        );
        let result = polynomial.roots(1000, &tolerance);
        println!(
            "{bits} bits, {scale_text} * (x^3 - 2), tolerance 1e-50: {}",
            if result.is_ok() { "Ok" } else { "Err" }
        );

        // The upstream fix must converge even with the scaled coefficients.
        let roots = result.expect("root scaling regression: expected convergence");
        assert_eq!(roots.len(), 3);
        // Check the UNscaled polynomial using higher-precision arithmetic.
        let mut max_residual = Float::with_val(768, 0);
        let two = Complex::new(Float::with_val(768, 2), Float::with_val(768, 0));
        for root in roots {
            let z = Complex::new(
                Float::with_val(768, root.re.as_raw()),
                Float::with_val(768, root.im.as_raw()),
            );
            let residual = (&z * &z * &z - &two).norm().re;
            if residual > max_residual {
                max_residual = residual;
            }
        }
        println!("  max |z^3 - 2|, reevaluated at 768 bits: {max_residual:.5e}");
        assert!(max_residual < Float::parse("1e-48", Some(768)).unwrap());
    }
}
