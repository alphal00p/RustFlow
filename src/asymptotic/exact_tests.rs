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
fn redundant_analytic_origin_selectors_preserve_the_exact_affine_space() {
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    let basis = prepared(&[
        &["0", "1/z", "0", "0", "0", "0", "0"],
        &[
            "-(z+3)/((z+1)*(z+9))",
            "-2*(z+5)/((z+1)*(z+9))",
            "0",
            "0",
            "0",
            "0",
            "6/((z+1)*(z+9))",
        ],
        &["0", "0", "-(z+2)/(z*(z+4))", "0", "0", "0", "2/(z*(z+4))"],
        &["2/(z*(z+4))", "0", "0", "-(z+2)/(z*(z+4))", "1/z", "0", "0"],
        &[
            "0",
            "2/(z*(z+4))",
            "6/((z+1)*(z+9))",
            "-(z+3)/((z+1)*(z+9))",
            "-(z+2)/(z*(z+4))-2*(z+5)/((z+1)*(z+9))",
            "0",
            "0",
        ],
        &["0", "0", "0", "1/(1+z)", "0", "0", "0"],
        &["0", "0", "0", "0", "0", "0", "0"],
    ])
    .exact_coefficients(2, &limits, &context)
    .unwrap();
    let mut relations = Vec::new();
    for component in 0..7 {
        for log in 0..7 {
            relations.push(coefficient(component, "-1/2", log, "0"));
        }
        for log in 1..7 {
            relations.push(coefficient(component, "0", log, "0"));
        }
    }
    for (component, value) in [(1, "0"), (2, "1"), (4, "0"), (5, "0"), (6, "1")] {
        relations.push(coefficient(component, "0", 0, value));
    }
    relations.push(ExactAsymptoticRelation {
        terms: vec![
            (selector(3, "0", 0), atom("1")),
            (selector(0, "0", 0), atom("-1")),
        ],
        value: Atom::new(),
    });
    assert_eq!(relations.len(), 97);
    let constraints = constraints(relations);
    let space = basis.constrain(&constraints, &limits, &context).unwrap();
    assert_eq!(space.free_amplitudes(), 1);
    let (offset, directions) = basis.finite_projection(&space, 7, &context).unwrap();
    for component in [1, 4, 5] {
        assert!(offset[component].is_zero());
    }
    for component in [2, 6] {
        assert_eq!(offset[component], one());
    }
    assert_eq!(offset[0], offset[3]);
    assert_eq!(directions[0], directions[3]);
    assert!(!directions[0][0].is_zero());
    for component in [1, 2, 4, 5, 6] {
        assert!(directions[component].iter().all(Gaussian::is_zero));
    }
    // Compaction does not waive the allocation bound on the original input.
    let too_small = ExactFrobeniusLimits {
        max_scalar_cells: constraints.relations.len() * 8 - 1,
        ..limits
    };
    assert!(matches!(
        basis.constrain(&constraints, &too_small, &context),
        Err(Error::Limit(_))
    ));
}

#[test]
fn exact_duplicate_rows_do_not_hide_inconsistent_right_hand_sides() {
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    let basis = prepared(&[&["0", "1/z"], &["0", "0"]])
        .exact_coefficients(2, &limits, &context)
        .unwrap();
    for contradiction in [
        coefficient(0, "0", 1, "1/10^300"),
        coefficient(1, "0", 1, "1/10^300"),
    ] {
        let mut relations = vec![coefficient(0, "0", 1, "0"); 97];
        relations.push(contradiction);
        assert!(matches!(
            basis.constrain(&constraints(relations), &limits, &context),
            Err(Error::InvalidInput(_))
        ));
    }
}

#[test]
fn tautological_rows_leave_all_amplitudes_free_but_distinct_rows_keep_the_guard() {
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    let basis = prepared(&[&["0", "1/z"], &["0", "0"]])
        .exact_coefficients(2, &limits, &context)
        .unwrap();
    let space = basis
        .constrain(
            &constraints(vec![coefficient(1, "0", 1, "0"); 97]),
            &limits,
            &context,
        )
        .unwrap();
    assert_eq!(space.free_amplitudes(), 2);
    let distinct = (1..=97)
        .map(|weight| ExactAsymptoticRelation {
            terms: vec![(selector(0, "0", 1), Atom::num(weight))],
            value: Atom::new(),
        })
        .collect();
    assert!(matches!(
        basis.constrain(&constraints(distinct), &limits, &context),
        Err(Error::Limit(_))
    ));
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
