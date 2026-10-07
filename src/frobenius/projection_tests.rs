use super::*;

#[test]
fn weighted_endpoint_preserves_pole_cancellation_across_distinct_exponents() -> Result<()> {
    let epsilon = symbol!("projection_prefix_tests::epsilon");
    let eta = symbol!("projection_prefix_tests::eta");
    let parameter = symbol!("projection_prefix_tests::a");
    let z = Atom::var(eta);
    for digits in [30, 60, 100] {
        let p = Precision::decimal(digits + 12)?;
        let column = |exponent, rows: [[i64; 3]; 3]| FrobeniusColumn {
            exponent: Atom::num(exponent),
            coefficients: rows
                .into_iter()
                .map(|row| vec![row.into_iter().map(|v| p.i(v)).collect()])
                .collect(),
        };
        let basis = FrobeniusBasis {
            columns: vec![
                column(0, [[1, 1, 0], [0, -1, 0], [0, 0, 0]]),
                column(-1, [[0, 0, 0], [2, 2, 0], [0, -3, 0]]),
                column(0, [[3, 3, 0], [0, -5, 0], [0, 0, 0]]),
            ],
            precision: p,
        };
        let parameters = ahash::HashMap::from_iter([(
            Atom::var(parameter),
            p.rational(&Rational::from((3, 7))),
        )]);
        let weights = [
            1 / (&z * (1 - Atom::var(parameter) * &z)),
            -1 / &z,
            Atom::new(),
        ];
        // The matched functions are y0=14, y1=14-22*eta, y2=0.
        // Their weighted limit is 14*a+22=28, after cancelling 14/eta.
        let value = project_limit(
            &basis,
            &[p.i(1), p.i(2), p.i(3)],
            &weights,
            epsilon,
            eta,
            &parameters,
        )?;
        assert!(p.close(&value, &p.i(28), digits));
    }
    Ok(())
}
