use super::*;
use crate::asymptotic::{AsymptoticSelector, ExactAsymptoticConstraints, ExactAsymptoticRelation};
use crate::frobenius::ExactFrobeniusLimits;

fn a(text: &str) -> Atom {
    Atom::parse(text, "sheet_proof_tests", Default::default()).unwrap()
}
fn z() -> Symbol {
    symbol!("sheet_proof_tests::z")
}
fn r() -> Symbol {
    symbol!("sheet_proof_tests::r")
}
fn s() -> Symbol {
    symbol!("sheet_proof_tests::s")
}
fn system(entries: &[&[&str]], roots: Vec<SquareRoot>) -> AlgebraicSystem {
    AlgebraicSystem::ordinary(
        DifferentialSystem {
            variable: z(),
            matrix: entries
                .iter()
                .map(|row| row.iter().map(|e| a(e)).collect())
                .collect(),
        },
        roots,
    )
}
fn prepared(entries: &[&[&str]], roots: Vec<SquareRoot>) -> PreparedAlgebraicFrobenius {
    system(entries, roots)
        .prepare_frobenius(32, &RunContext::default())
        .unwrap()
}
fn root(symbol: Symbol, radicand: &str) -> SquareRoot {
    SquareRoot {
        symbol,
        radicand: a(radicand),
    }
}
fn constraints(power: &str, log_power: usize) -> ExactAsymptoticConstraints {
    ExactAsymptoticConstraints {
        provenance: "analytic exact physical support".into(),
        relations: vec![ExactAsymptoticRelation::coefficient(
            AsymptoticSelector {
                component: 0,
                power: a(power),
                log_power,
            },
            Atom::new(),
        )],
    }
}

#[test]
fn root_sheet_action_retains_one_physical_direction_on_both_sheets() {
    let prepared = prepared(&[&["1/r"]], vec![root(r(), "z")]);
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    let p = Precision::decimal(70).unwrap();
    let exact = prepared
        .prepared
        .exact_coefficients(16, &limits, &context)
        .unwrap();
    for seed in [RootSeed::Principal, RootSeed::Opposite] {
        let rows = prepared
            .sheet_equations(
                &exact,
                &a("1/16"),
                &BTreeMap::from([(r(), seed)]),
                0,
                p,
                &limits,
                &context,
            )
            .unwrap()
            .unwrap();
        let space = exact
            .constrain_with_rows(&constraints("-1", 0), rows, &limits, &context)
            .unwrap();
        assert_eq!(space.free_amplitudes(), 1);
        assert!(exact.finite_projection(&space, 1, &context).is_ok());
    }
}

#[test]
fn complex_phase_and_winding_use_power_of_principal_coordinate_root() {
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    let p = Precision::decimal(70).unwrap();
    for (radicand, point, winding, expected) in [
        ("z", "1/16", 0, "1"),
        ("z", "1/16", 1, "-1"),
        ("z^2", "-1/16", 0, "-1"),
        ("z^3", "-1/16", 0, "-1"),
        ("z^3", "-1/16", 1, "1"),
        ("z^-3", "-1/16", 0, "1"),
        ("-z", "𝑖/16", 0, "-𝑖"),
    ] {
        let lifted = system(&[&["r"]], vec![root(r(), radicand)])
            .rational_lift(8, &context)
            .unwrap();
        let (_, leading, _) = super::sheet_germ::selected_germ(
            &lifted,
            1,
            &a(point),
            &BTreeMap::new(),
            winding,
            p,
            &limits,
            &context,
        )
        .unwrap();
        assert_eq!(
            leading,
            crate::frobenius::exact::gaussian(&a(expected)).unwrap(),
            "{radicand}, {point}, winding={winding}"
        );
    }
}

#[test]
fn product_only_generator_keeps_declared_individual_germ_signs() {
    let prepared = prepared(&[&["1/(r*s)"]], vec![root(r(), "2*z"), root(s(), "8*z")]);
    assert_eq!(prepared.lifted.monomials, vec![0, 3]);
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    let p = Precision::decimal(70).unwrap();
    for (point, seed, expected) in [
        ("1/16", RootSeed::Principal, "4"),
        ("1/16", RootSeed::Opposite, "-4"),
        ("-1/16", RootSeed::Principal, "4"),
    ] {
        let (_, leading, _) = super::sheet_germ::selected_germ(
            &prepared.lifted,
            3,
            &a(point),
            &BTreeMap::from([(s(), seed)]),
            0,
            p,
            &limits,
            &context,
        )
        .unwrap();
        assert_eq!(
            leading,
            crate::frobenius::exact::gaussian(&a(expected)).unwrap()
        );
    }
}

#[test]
fn all_positive_integer_resonances_must_be_crossed() {
    let prepared = prepared(&[&["0", "r"], &["0", "21/z"]], vec![root(r(), "z")]);
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    assert!(!prepared.complete_order(8, &limits, &context).unwrap());
    assert!(prepared.complete_order(32, &limits, &context).unwrap());
    let p = Precision::decimal(70).unwrap();
    let short = prepared
        .prepared
        .exact_coefficients(8, &limits, &context)
        .unwrap();
    assert!(
        prepared
            .sheet_equations(
                &short,
                &a("1/16"),
                &BTreeMap::new(),
                0,
                p,
                &limits,
                &context
            )
            .unwrap()
            .is_none()
    );
    let complete = prepared
        .prepared
        .exact_coefficients(32, &limits, &context)
        .unwrap();
    assert!(
        prepared
            .sheet_equations(
                &complete,
                &a("1/16"),
                &BTreeMap::new(),
                0,
                p,
                &limits,
                &context
            )
            .unwrap()
            .is_some()
    );
}

#[test]
fn product_generator_sheet_changes_which_physical_directions_are_finite() {
    let prepared = prepared(&[&["1/(r*s)"]], vec![root(r(), "2*z"), root(s(), "8*z")]);
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    let p = Precision::decimal(70).unwrap();
    let exact = prepared
        .prepared
        .exact_coefficients(16, &limits, &context)
        .unwrap();
    for (seed, finite) in [(RootSeed::Principal, true), (RootSeed::Opposite, false)] {
        let rows = prepared
            .sheet_equations(
                &exact,
                &a("1/16"),
                &BTreeMap::from([(s(), seed)]),
                0,
                p,
                &limits,
                &context,
            )
            .unwrap()
            .unwrap();
        let space = exact
            .constrain_with_rows(&constraints("-1", 0), rows.clone(), &limits, &context)
            .unwrap();
        assert_eq!(space.free_amplitudes(), 1);
        assert_eq!(exact.finite_projection(&space, 1, &context).is_ok(), finite);
        if !finite {
            let vanishing = exact
                .constrain_with_rows(&constraints("-1/4", 0), rows, &limits, &context)
                .unwrap();
            assert_eq!(vanishing.free_amplitudes(), 0);
            assert!(exact.finite_projection(&vanishing, 1, &context).is_ok());
        }
    }
}

#[test]
fn dependent_registered_roots_select_a_single_physical_space() {
    let prepared = prepared(&[&["1/r+1/s"]], vec![root(r(), "z"), root(s(), "4*z")]);
    assert_eq!(prepared.lifted.monomials.len(), 4);
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    let p = Precision::decimal(70).unwrap();
    let exact = prepared
        .prepared
        .exact_coefficients(16, &limits, &context)
        .unwrap();
    for seed in [RootSeed::Principal, RootSeed::Opposite] {
        let rows = prepared
            .sheet_equations(
                &exact,
                &a("1/16"),
                &BTreeMap::from([(s(), seed)]),
                0,
                p,
                &limits,
                &context,
            )
            .unwrap()
            .unwrap();
        let space = exact
            .constrain_with_rows(&constraints("-1", 0), rows, &limits, &context)
            .unwrap();
        assert_eq!(space.free_amplitudes(), 1);
        assert!(exact.finite_projection(&space, 1, &context).is_ok());
    }
}

#[test]
fn nonconstant_normalized_root_series_has_exact_constant_sheet_action() {
    let prepared = prepared(&[&["1/r"]], vec![root(r(), "z*(1-z)")]);
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    let p = Precision::decimal(70).unwrap();
    let exact = prepared
        .prepared
        .exact_coefficients(16, &limits, &context)
        .unwrap();
    let rows = prepared
        .sheet_equations(
            &exact,
            &a("1/16"),
            &BTreeMap::new(),
            0,
            p,
            &limits,
            &context,
        )
        .unwrap()
        .unwrap();
    let space = exact
        .constrain_with_rows(&constraints("-1", 0), rows, &limits, &context)
        .unwrap();
    assert_eq!(space.free_amplitudes(), 1);
    assert!(exact.finite_projection(&space, 1, &context).is_ok());
}

#[test]
fn aggregate_growth_is_refused_before_root_conversion_and_sign_selection() {
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    assert!(matches!(
        super::sheet_germ::valuation(&a("1+z^257"), z(), &limits, &context),
        Err(Error::Limit(_))
    ));
    let source = system(&[&["r"]], vec![root(r(), "z")])
        .rational_lift(8, &context)
        .unwrap();
    let p = Precision::decimal(70).unwrap();
    assert!(matches!(
        super::sheet_germ::selected_germ(
            &source,
            1,
            &a("2^-3000+𝑖*2^-3001"),
            &BTreeMap::new(),
            1,
            p,
            &limits,
            &context
        ),
        Err(Error::Limit(_))
    ));
    context.cancellation.cancel();
    assert!(matches!(
        super::sheet_germ::valuation(&a("1+z"), z(), &limits, &context),
        Err(Error::Cancelled)
    ));
}

#[test]
fn pinned_native_algebraic_field_supports_future_non_gaussian_matrix_work() {
    use symbolica::domains::algebraic::AlgebraicContext;
    let generator = a("2^(1/2)");
    let mut context = AlgebraicContext::from_generators(std::slice::from_ref(&generator)).unwrap();
    let coefficient = context.convert_atom(a("1+2^(1/2)").as_view()).unwrap();
    let field = context.field().clone();
    let matrix = Matrix::from_nested_vec(vec![vec![coefficient]], field.clone()).unwrap();
    assert!(field.is_one(&(&matrix * &matrix.inv().unwrap())[(0, 0)]));
    // Capability evidence only: endpoint admission still refuses that field.
}

#[test]
fn non_gaussian_leading_constant_and_unproved_numeric_seed_remain_typed() {
    let limits = ExactFrobeniusLimits::default();
    let context = RunContext::default();
    let p = Precision::decimal(70).unwrap();
    let irrational = prepared(&[&["1/r"]], vec![root(r(), "2*z")]);
    assert!(matches!(
        super::sheet_germ::selected_germ(
            &irrational.lifted,
            1,
            &a("1/16"),
            &BTreeMap::new(),
            0,
            p,
            &limits,
            &context
        ),
        Err(Error::Unsupported(_))
    ));
    let rational = prepared(&[&["1/r"]], vec![root(r(), "z")]);
    assert!(matches!(
        super::sheet_germ::selected_germ(
            &rational.lifted,
            1,
            &a("1/16"),
            &BTreeMap::from([(r(), RootSeed::Value(p.i(1)))]),
            0,
            p,
            &limits,
            &context
        ),
        Err(Error::Unsupported(_))
    ));
}
