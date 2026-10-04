use symbolica::prelude::*;
use symbolica_amflow::{
    BoundaryData, ComplexFloat, DifferentialSystem, FlowOptions, Precision, Result, RunContext,
};

#[test]
fn exact_complex_matrix_literals_compile_and_follow_the_analytic_solution() -> Result<()> {
    let x = symbol!("complex_literal_ode::x");
    let i = Atom::num(Complex::new(Rational::from(0), Rational::from(1)));
    let system = DifferentialSystem {
        variable: x,
        matrix: vec![vec![&i + Atom::one() / (Atom::var(x) + 10 - &i)]],
    };
    for (digits, order) in [(50, 48), (70, 80)] {
        let p = Precision::decimal(digits)?;
        let compiled = system.compile(p, &Default::default())?;
        assert_eq!(compiled.poles.len(), 1);
        assert!(p.close(
            &compiled.poles[0],
            &ComplexFloat::new(p.real(-10), p.real(1)),
            35
        ));
        let result = compiled.transport(
            &BoundaryData {
                point: p.zero(),
                values: vec![p.i(1)],
            },
            &[p.i(2)],
            &FlowOptions {
                digits: 30,
                guard_digits: digits - 30,
                series_order: order,
                ..Default::default()
            },
            &RunContext::default(),
        )?;
        let expected = p.mul(
            &p.exp(&ComplexFloat::new(p.real(0), p.real(2))),
            &p.div(
                &ComplexFloat::new(p.real(12), p.real(-1)),
                &ComplexFloat::new(p.real(10), p.real(-1)),
            ),
        );
        assert!(p.close(&result.values[0], &expected, 30));
    }
    Ok(())
}

#[test]
fn reserved_imaginary_symbol_cannot_alias_the_ode_variable() -> Result<()> {
    let system = DifferentialSystem {
        variable: symbol!("symbolica_amflow::imaginary_unit"),
        matrix: vec![vec![Atom::num(Complex::new(
            Rational::from(0),
            Rational::from(1),
        ))]],
    };
    assert!(matches!(
        system.validate(),
        Err(symbolica_amflow::Error::InvalidInput(_))
    ));
    assert!(matches!(
        system.compile(Precision::decimal(40)?, &Default::default()),
        Err(symbolica_amflow::Error::InvalidInput(_))
    ));
    Ok(())
}
