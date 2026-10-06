use std::sync::Arc;
use symbolica::prelude::*;
use symbolica_amflow::*;

#[test]
fn symbolic_epsilon_preparation_reuses_resonant_logarithms_at_increased_precision() {
    let x = symbol!("prepared_logs::x");
    let epsilon = parse!("prepared_logs::epsilon");
    let variable = Atom::var(x);
    let system = DifferentialSystem {
        variable: x,
        matrix: vec![
            vec![&epsilon / &variable, Atom::num(1) / &variable, Atom::new()],
            vec![Atom::new(), &epsilon / &variable, Atom::new()],
            vec![
                Atom::num(1),
                Atom::new(),
                (&epsilon + Atom::num(1)) / &variable,
            ],
        ],
    };
    let prepared = system.prepare_frobenius(&RunContext::default()).unwrap();
    assert_eq!(prepared.exponents().count(), 2);
    for (digits, sample, order) in [(55, 31, 12), (75, 43, 20)] {
        let p = Precision::decimal(digits).unwrap();
        let parameters = ahash::HashMap::from_iter([(
            epsilon.clone(),
            p.rational(&Rational::from((1, sample))),
        )]);
        let basis = prepared
            .evaluate(p, &parameters, order, &RunContext::default())
            .unwrap();
        assert!(basis.columns.iter().all(|column| {
            let difference = (&column.exponent - &epsilon).together().cancel();
            difference.is_zero() || difference.is_one()
        }));
        assert!(
            basis
                .columns
                .iter()
                .any(|column| column.coefficients.iter().any(|logs| logs.len() > 1))
        );
        let start = p.rational(&Rational::from((1, 10)));
        let end = p.rational(&Rational::from((1, 5)));
        let initial = basis.evaluate(&start, &parameters).unwrap();
        let expected = basis.evaluate(&end, &parameters).unwrap();
        let compiled = system.compile(p, &parameters).unwrap();
        let options = FlowOptions {
            digits: 30,
            series_order: 64,
            ..Default::default()
        };
        for column in 0..3 {
            let result = compiled
                .transport(
                    &BoundaryData {
                        point: start.clone(),
                        values: initial.iter().map(|row| row[column].clone()).collect(),
                    },
                    std::slice::from_ref(&end),
                    &options,
                    &RunContext::default(),
                )
                .unwrap();
            for (row, actual) in result.values.iter().enumerate() {
                assert!(p.close(actual, &expected[row][column], 28));
            }
        }
    }
}

#[test]
fn uncancelled_denominator_domain_survives_exact_preparation() {
    let system = DifferentialSystem {
        variable: symbol!("prepared_domain::x"),
        matrix: vec![vec![parse!(
            "(prepared_domain::epsilon^2-1)/((prepared_domain::epsilon-1)*prepared_domain::x)"
        )]],
    };
    let prepared = system.prepare_frobenius(&RunContext::default()).unwrap();
    assert!(!prepared.nonzero_conditions().is_empty());
    assert_eq!(prepared.system().matrix, system.matrix);
    let p = Precision::decimal(60).unwrap();
    let parameters = ahash::HashMap::from_iter([(parse!("prepared_domain::epsilon"), p.i(1))]);
    assert!(
        matches!(prepared.evaluate(p, &parameters, 12, &RunContext::default()), Err(Error::Numerical(message)) if message.contains("denominator condition"))
    );
    let parameters = ahash::HashMap::from_iter([(parse!("prepared_domain::epsilon"), p.i(2))]);
    // The excluded expansion center x=0 is still valid for a singular series.
    let basis = prepared
        .evaluate(p, &parameters, 12, &RunContext::default())
        .unwrap();
    let point = p.rational(&Rational::from((1, 2)));
    assert!(p.close(
        &basis.evaluate(&point, &parameters).unwrap()[0][0],
        &p.rational(&Rational::from((1, 8))),
        50
    ));
}

#[test]
fn cancellation_of_numerical_evaluation_does_not_change_preparation() {
    let system = DifferentialSystem {
        variable: symbol!("prepared_cancel::x"),
        matrix: vec![vec![parse!("1/prepared_cancel::x")]],
    };
    let prepared = system.prepare_frobenius(&RunContext::default()).unwrap();
    let cancellation = CancellationToken::default();
    let cancel = cancellation.clone();
    let context = RunContext {
        cancellation,
        progress: Some(Arc::new(move |_| cancel.cancel())),
    };
    let p = Precision::decimal(40).unwrap();
    assert!(matches!(
        prepared.evaluate(p, &Default::default(), 10, &context),
        Err(Error::Cancelled)
    ));
    assert!(
        prepared
            .evaluate(p, &Default::default(), 10, &RunContext::default())
            .is_ok()
    );
}

#[test]
fn nongeneric_parameter_collision_cannot_return_degenerate_columns() {
    let system = DifferentialSystem {
        variable: symbol!("prepared_collision::x"),
        matrix: vec![
            vec![
                parse!("prepared_collision::a/prepared_collision::x"),
                Atom::new(),
            ],
            vec![parse!("1/prepared_collision::x"), Atom::new()],
        ],
    };
    let prepared = system.prepare_frobenius(&RunContext::default()).unwrap();
    let p = Precision::decimal(60).unwrap();
    let parameters = ahash::HashMap::from_iter([(parse!("prepared_collision::a"), p.zero())]);
    assert!(
        matches!(prepared.evaluate(p, &parameters, 12, &RunContext::default()), Err(Error::Unsupported(message)) if message.contains("resonance"))
    );
    let parameters = ahash::HashMap::from_iter([(
        parse!("prepared_collision::a"),
        p.rational(&Rational::from((1, 7))),
    )]);
    assert!(
        prepared
            .evaluate(p, &parameters, 12, &RunContext::default())
            .is_ok()
    );
}

#[test]
fn near_integer_resonance_is_rejected_from_both_sides_and_signs() {
    let x = symbol!("portable_resonance::x");
    let a = symbol!("portable_resonance::a");
    let system = DifferentialSystem {
        variable: x,
        matrix: vec![
            vec![Atom::var(a) / Atom::var(x), Atom::new()],
            vec![Atom::num(1) / Atom::var(x), Atom::new()],
        ],
    };
    let prepared = system.prepare_frobenius(&RunContext::default()).unwrap();
    let p = Precision::decimal(60).unwrap();
    let displacement = Rational::from((Integer::one(), Integer::from(10).pow(55)));
    for sign in [-1, 1] {
        for side in [-1, 1] {
            let value = Rational::from(sign) + &(&displacement * &Rational::from(side));
            let parameters = ahash::HashMap::from_iter([(Atom::var(a), p.rational(&value))]);
            assert!(
                matches!(prepared.evaluate(p, &parameters, 12, &RunContext::default()), Err(Error::Unsupported(message)) if message.contains("resonance")),
                "near integer {sign}, side {side} was admitted"
            );
        }
    }
    // The integer candidate is zero; obtaining it must not materialize a
    // denominator proportional to this finite float's very negative exponent.
    let tiny = p.parse("1e-10000000", "0").unwrap();
    assert!(!tiny.re.is_zero());
    let parameters = ahash::HashMap::from_iter([(Atom::var(a), tiny)]);
    assert!(
        matches!(prepared.evaluate(p, &parameters, 12, &RunContext::default()), Err(Error::Unsupported(message)) if message.contains("resonance"))
    );
}

#[test]
fn tiny_nonzero_parameters_do_not_acquire_an_artificial_exclusion_radius() {
    let system = DifferentialSystem {
        variable: symbol!("prepared_small::x"),
        matrix: vec![vec![parse!(
            "(prepared_small::a^2+prepared_small::a)/(prepared_small::a*prepared_small::x)"
        )]],
    };
    let prepared = system.prepare_frobenius(&RunContext::default()).unwrap();
    let p = Precision::decimal(60).unwrap();
    let parameters =
        ahash::HashMap::from_iter([(parse!("prepared_small::a"), p.parse("1e-200", "0").unwrap())]);
    assert!(
        prepared
            .evaluate(p, &parameters, 8, &RunContext::default())
            .is_ok()
    );
}
