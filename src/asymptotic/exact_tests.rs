use super::*;
use crate::{DifferentialSystem, Precision};

fn atom(value: &str) -> Atom {
    Atom::parse(value, "exact_asymptotic_tests", Default::default()).unwrap()
}
fn prepared(matrix: &[&[&str]]) -> crate::frobenius::PreparedFrobenius {
    DifferentialSystem {
        variable: symbol!("exact_asymptotic_tests::z"),
        matrix: matrix
            .iter()
            .map(|row| row.iter().map(|a| atom(a)).collect())
            .collect(),
    }
    .prepare_frobenius(&RunContext::default())
    .unwrap()
}
fn selector(component: usize, power: &str, log_power: usize) -> AsymptoticSelector {
    AsymptoticSelector {
        component,
        power: atom(power),
        log_power,
    }
}
fn constraints(rows: Vec<ExactAsymptoticRelation>) -> ExactAsymptoticConstraints {
    ExactAsymptoticConstraints {
        relations: rows,
        provenance: "explicit analytic coefficient identities".into(),
    }
}
fn coefficient(component: usize, power: &str, log: usize, value: &str) -> ExactAsymptoticRelation {
    ExactAsymptoticRelation::coefficient(selector(component, power, log), atom(value))
}

#[test]
fn exact_log_constraint_preserves_free_amplitude_and_rejects_any_nonzero_log() {
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    let basis = prepared(&[&["0", "1/z"], &["0", "0"]])
        .exact_coefficients(8, &limits, &context)
        .unwrap();
    for value in ["0", "1/10^300"] {
        let space = basis
            .constrain(
                &constraints(vec![coefficient(0, "0", 1, value)]),
                &limits,
                &context,
            )
            .unwrap();
        assert_eq!(space.free_amplitudes(), 1);
        let finite = basis.finite_projection(&space, 2, &context);
        if value == "0" {
            let (offset, directions) = finite.unwrap();
            assert!(offset.iter().all(Gaussian::is_zero));
            assert_eq!(directions[0], vec![one()]);
            assert_eq!(directions[1], vec![zero()]);
        } else {
            assert!(matches!(finite, Err(Error::Accuracy(_))));
        }
    }
    let fixed = basis
        .constrain(
            &constraints(vec![
                coefficient(0, "0", 1, "0"),
                coefficient(0, "0", 0, "3"),
            ]),
            &limits,
            &context,
        )
        .unwrap();
    assert_eq!(fixed.free_amplitudes(), 0);
    let (offset, directions) = basis.finite_projection(&fixed, 2, &context).unwrap();
    assert_eq!(offset, vec![gaussian(&atom("3")).unwrap(), zero()]);
    assert!(directions.iter().all(Vec::is_empty));
}

#[test]
fn correlated_negative_power_identity_cancels_only_the_asserted_subspace() {
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    let basis = prepared(&[&["-1/(2*z)", "1/(2*z)"], &["1/(2*z)", "-1/(2*z)"]])
        .exact_coefficients(8, &limits, &context)
        .unwrap();
    let relation = ExactAsymptoticRelation {
        terms: vec![
            (selector(0, "-1", 0), atom("1")),
            (selector(1, "-1", 0), atom("-1")),
        ],
        value: atom("0"),
    };
    let space = basis
        .constrain(&constraints(vec![relation]), &limits, &context)
        .unwrap();
    assert_eq!(space.free_amplitudes(), 1);
    let (_, directions) = basis.finite_projection(&space, 2, &context).unwrap();
    assert_eq!(directions[0], directions[1]);
    assert!(directions[0].iter().any(|a| !a.is_zero()));
    // Fixing a finite coefficient alone cannot eliminate the independent pole.
    let unconstrained = basis
        .constrain(
            &constraints(vec![coefficient(0, "0", 0, "1")]),
            &limits,
            &context,
        )
        .unwrap();
    assert!(matches!(
        basis.finite_projection(&unconstrained, 2, &context),
        Err(Error::Accuracy(_))
    ));
}

#[test]
fn exact_later_resonance_and_normalization_prefix_are_not_rounded_away() {
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    let basis = prepared(&[&["0", "1/10^300"], &["0", "-1/z"]])
        .exact_coefficients(8, &limits, &context)
        .unwrap();
    assert!(
        basis
            .coefficient_row(&selector(0, "0", 1))
            .unwrap()
            .iter()
            .any(|a| !a.is_zero())
    );
    let space = basis
        .constrain(
            &constraints(vec![coefficient(0, "0", 1, "0")]),
            &limits,
            &context,
        )
        .unwrap();
    basis.finite_projection(&space, 2, &context).unwrap();
    let prepared = prepared(&[&["0", "1/z^21"], &["0", "0"]]);
    let short = prepared.exact_coefficients(8, &limits, &context).unwrap();
    assert!(matches!(
        short.coefficient_row(&selector(0, "0", 0)),
        Err(Error::Accuracy(_))
    ));
    let short_space = short
        .constrain(
            &constraints(vec![coefficient(0, "-20", 0, "0")]),
            &limits,
            &context,
        )
        .unwrap();
    assert!(matches!(
        short.finite_projection(&short_space, 2, &context),
        Err(Error::Accuracy(_))
    ));
    let long = prepared.exact_coefficients(40, &limits, &context).unwrap();
    let space = long
        .constrain(
            &constraints(vec![coefficient(0, "-20", 0, "0")]),
            &limits,
            &context,
        )
        .unwrap();
    long.finite_projection(&space, 2, &context).unwrap();
}

#[test]
fn exact_and_legacy_recurrences_agree_across_logs_resonances_and_shears() {
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    let p = Precision::decimal(90).unwrap();
    for matrix in [
        vec![vec!["0", "1/z"], vec!["0", "0"]],
        vec![vec!["0", "1"], vec!["0", "-1/z"]],
        vec![vec!["1/(2*z)", "1/z"], vec!["0", "1/(2*z)"]],
        vec![vec!["0", "1/z^3"], vec!["0", "0"]],
        vec![vec!["1/(1-z)", "1/(1+z)"], vec!["0", "0"]],
    ] {
        let rows = matrix.iter().map(Vec::as_slice).collect::<Vec<_>>();
        let prepared = prepared(&rows);
        let exact = prepared
            .exact_coefficients(24, &limits, &context)
            .unwrap()
            .numerical(p);
        let numerical = prepared
            .evaluate(p, &Default::default(), 24, &context)
            .unwrap();
        for winding in [-1, 0, 1] {
            let point = p.rational(&Rational::from((1, 16)));
            let exact = exact
                .evaluate_with_winding(&point, &Default::default(), winding)
                .unwrap();
            let numerical = numerical
                .evaluate_with_winding(&point, &Default::default(), winding)
                .unwrap();
            for (a, b) in exact.iter().flatten().zip(numerical.iter().flatten()) {
                assert!(p.close(a, b, 60), "{a} != {b}");
            }
        }
    }
}

#[test]
fn exact_relations_reject_contradictions_missing_orders_and_resource_excess() {
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    let prepared = prepared(&[&["0", "1/z"], &["0", "0"]]);
    let basis = prepared.exact_coefficients(8, &limits, &context).unwrap();
    assert!(matches!(
        basis.constrain(
            &constraints(vec![
                coefficient(0, "0", 1, "0"),
                coefficient(0, "0", 1, "1/10^300")
            ]),
            &limits,
            &context
        ),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        basis.constrain(
            &constraints(vec![coefficient(0, "9", 0, "0")]),
            &limits,
            &context
        ),
        Err(Error::Accuracy(_))
    ));
    assert!(matches!(
        prepared.exact_coefficients(
            8,
            &ExactFrobeniusLimits {
                max_dimension: 1,
                ..limits.clone()
            },
            &context
        ),
        Err(Error::Limit(_))
    ));
    context.cancellation.cancel();
    assert!(matches!(
        basis.constrain(
            &constraints(vec![coefficient(0, "0", 1, "0")]),
            &limits,
            &context
        ),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        prepared.exact_coefficients(8, &limits, &context),
        Err(Error::Cancelled)
    ));
}
#[test]
fn exact_complex_relations_share_native_gaussian_arithmetic_and_bound_input_height() {
    let z = symbol!("exact_asymptotic_tests::z");
    let imaginary = Atom::num(Complex::new(Rational::zero(), Rational::one()));
    let prepared = DifferentialSystem {
        variable: z,
        matrix: vec![
            vec![Atom::new(), &imaginary / Atom::var(z)],
            vec![Atom::new(), Atom::new()],
        ],
    }
    .prepare_frobenius(&RunContext::default())
    .unwrap();
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    let basis = prepared.exact_coefficients(8, &limits, &context).unwrap();
    let space = basis
        .constrain(
            &constraints(vec![ExactAsymptoticRelation {
                terms: vec![(selector(0, "0", 1), imaginary)],
                value: Atom::new(),
            }]),
            &limits,
            &context,
        )
        .unwrap();
    let (_, directions) = basis.finite_projection(&space, 2, &context).unwrap();
    assert_eq!(space.free_amplitudes(), 1);
    assert!(directions[0].iter().any(|value| !value.is_zero()));
    assert!(directions[1].iter().all(Gaussian::is_zero));
    let low = ExactFrobeniusLimits {
        max_coefficient_bits: 128,
        ..limits
    };
    assert!(matches!(
        basis.constrain(
            &constraints(vec![coefficient(0, "0", 1, "1/10^300")]),
            &low,
            &context
        ),
        Err(Error::Limit(_))
    ));
}
