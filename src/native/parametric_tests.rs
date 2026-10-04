use super::*;
use rustred::solver::{Case, ExceptionalConditions, Power, RuleCandidate, SectorConfig};
fn native_family() -> rustred::family::IntegralFamily {
    let (family, _) = crate::benchmarks::paper_two_loop().unwrap();
    let (family, _) = family
        .deform(
            symbol!("parametric_tests::eta"),
            &crate::MassMode::Propagators(vec![4]),
        )
        .unwrap();
    family
        .convert_at_epsilon(Some(&Rational::from((1, 2700))))
        .unwrap()
        .family
}
fn setup() -> SourceSystem<9> {
    SourceSystem::from_family_with_lorentz(&native_family(), false).unwrap()
}
fn sector() -> [bool; 9] {
    [false, true, false, false, true, false, true, false, false]
}
fn polynomial(
    sources: &SourceSystem<9>,
    axis: Option<usize>,
    constant: i64,
) -> CoefficientPolynomial {
    let mut p = sources.rows()[0][0].coefficient.zero();
    let mut exponents = vec![0; p.nvars()];
    if let Some(axis) = axis {
        exponents[sources.index_variables()[axis]] = 1;
        p.append_monomial(Integer::from(1), &exponents);
        exponents.fill(0);
    }
    if constant != 0 {
        p.append_monomial(Integer::from(constant), &exponents);
    }
    p
}
fn rule(sources: &SourceSystem<9>, delta: i16) -> SectorRule<9> {
    let case: Case<9> = CoordinateCase::new(std::array::from_fn(|i| (!sector()[i]).then_some(0)))
        .unwrap()
        .into();
    let target = case.integral();
    let mut powers = *target.powers();
    powers[1] = Power::new(true, delta).unwrap();
    let one = polynomial(sources, None, 1);
    SectorRule {
        candidate: RuleCandidate {
            case,
            target,
            rhs: vec![Term {
                integral: Integral::new(powers),
                coefficient: RationalPolynomial::from_num_den(one.clone(), one, &Z, false),
            }],
            sources: vec![],
            stats: Default::default(),
        },
        exceptions: ExceptionalConditions::default(),
        dispatch_policy: RuleDispatchPolicy::Partition,
    }
}
#[test]
fn fixed_indices_and_exception_conjunction_are_checked() {
    let sources = setup();
    let solver = SectorSolver::new(&sources, sector(), SectorConfig::default()).unwrap();
    let mut rule = rule(&sources, -1);
    rule.exceptions.branches = vec![vec![
        polynomial(&sources, Some(1), -3),
        polynomial(&sources, Some(4), -2),
    ]];
    assert!(
        apply(&rule, [0, 3, 0, 0, 2, 0, 2, 0, 0], &sources, &solver)
            .unwrap()
            .is_none()
    );
    assert!(
        apply(&rule, [0, 3, 0, 0, 3, 0, 2, 0, 0], &sources, &solver)
            .unwrap()
            .is_some()
    );
    assert!(
        apply(&rule, [-1, 3, 0, 0, 3, 0, 2, 0, 0], &sources, &solver)
            .unwrap()
            .is_none()
    );
}
#[test]
fn affine_equalities_and_original_denominators_are_checked() {
    let sources = setup();
    let solver = SectorSolver::new(&sources, sector(), SectorConfig::default()).unwrap();
    let mut rule = rule(&sources, -1);
    let equation = polynomial(&sources, Some(1), 0) - polynomial(&sources, Some(4), 0);
    rule.candidate.case = rule
        .candidate
        .case
        .intersect(&[equation], sources.index_variables(), &sector())
        .unwrap()
        .unwrap();
    rule.candidate.target = rule.candidate.case.integral();
    assert!(
        apply(&rule, [0, 3, 0, 0, 2, 0, 2, 0, 0], &sources, &solver)
            .unwrap()
            .is_none()
    );
    assert!(
        apply(&rule, [0, 3, 0, 0, 3, 0, 2, 0, 0], &sources, &solver)
            .unwrap()
            .is_some()
    );
    let zero = polynomial(&sources, Some(1), -3);
    rule.candidate.rhs[0].coefficient = RationalPolynomial {
        numerator: zero.clone(),
        denominator: zero,
    };
    assert!(
        apply(&rule, [0, 3, 0, 0, 3, 0, 2, 0, 0], &sources, &solver)
            .unwrap()
            .is_none()
    );
}
#[test]
fn parameter_denominators_survive_cancellation() {
    let sources = setup();
    let solver = SectorSolver::new(&sources, sector(), SectorConfig::default()).unwrap();
    let mut rule = rule(&sources, -1);
    let mut eta = polynomial(&sources, None, 0);
    let mut exponent = vec![0; eta.nvars()];
    exponent[0] = 1;
    eta.append_monomial(Integer::from(1), &exponent);
    rule.candidate.rhs[0].coefficient = RationalPolynomial {
        numerator: eta.clone(),
        denominator: eta.clone(),
    };
    let applied = apply(&rule, [0, 3, 0, 0, 3, 0, 2, 0, 0], &sources, &solver)
        .unwrap()
        .unwrap();
    assert!(applied.conditions.contains(&eta));
    assert_eq!(
        applied.terms[0].coefficient.numerator,
        applied.terms[0].coefficient.denominator
    );
}
#[test]
fn strict_descent_index_bounds_and_variable_maps_are_checked() {
    let sources = setup();
    let solver = SectorSolver::new(&sources, sector(), SectorConfig::default()).unwrap();
    let mut rule = rule(&sources, 1);
    assert!(apply(&rule, [0, 3, 0, 0, 3, 0, 2, 0, 0], &sources, &solver).is_err());
    assert!(apply(&rule, [0, 63, 0, 0, 3, 0, 2, 0, 0], &sources, &solver).is_err());
    let old = rule.candidate.rhs[0].coefficient.numerator.variables()[0].clone();
    rule.candidate.rhs[0].coefficient.numerator.rename_variable(
        &old,
        &symbolica::poly::PolyVariable::Symbol(symbol!("wrong_map")),
    );
    assert!(apply(&rule, [0, 3, 0, 0, 3, 0, 2, 0, 0], &sources, &solver).is_err());
}
#[test]
fn bank_is_bound_to_source_and_order() {
    let sources = setup();
    let other = setup();
    let solver = SectorSolver::new(&sources, sector(), SectorConfig::default()).unwrap();
    let mut bank = Bank::new(&sources, &solver);
    assert!(
        bank.try_reduce([0, 3, 0, 0, 3, 0, 2, 0, 0], &other, &solver, 3)
            .is_err()
    );
    let solver = SectorSolver::new(
        &sources,
        [true, true, false, false, true, false, true, false, false],
        SectorConfig::default(),
    )
    .unwrap();
    assert!(
        bank.try_reduce([1, 3, 0, 0, 3, 0, 2, 0, 0], &sources, &solver, 3)
            .is_err()
    );
}

#[test]
fn generated_tadpole_ray_matches_exact_gamma_recurrence() {
    let eta = symbol!("parametric_tadpole_eta");
    let family = crate::IntegralFamily {
        name: "ray_tadpole".into(),
        loops: vec!["k".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![crate::Propagator {
            constant: -Atom::num(1) - Atom::var(eta),
            scalar_products: vec![Atom::num(1)],
        }],
        physical_propagators: 1,
        epsilon: symbol!("parametric_tadpole_eps"),
        dimension: 4,
    };
    let epsilon = Rational::from((1, 2700));
    let converted = family.convert_at_epsilon(Some(&epsilon)).unwrap();
    let sources = SourceSystem::<1>::from_family_with_lorentz(&converted.family, false).unwrap();
    let solver = SectorSolver::new(&sources, [true], SectorConfig::default()).unwrap();
    let mut bank = Bank::new(&sources, &solver);
    // An exhausted scalar corner remains uncovered by the formula bank.
    assert!(
        bank.try_reduce([1], &sources, &solver, 3)
            .unwrap()
            .is_none()
    );
    for n in 2..=12 {
        let applied = bank.try_reduce([n], &sources, &solver, 3).unwrap().unwrap();
        assert_eq!(applied.terms.len(), 1);
        assert_eq!(
            applied.terms[0].integral,
            Integral::numeric([n - 1]).unwrap()
        );
        let actual = crate::family::substitute(
            &applied.terms[0].coefficient.to_expression(),
            &converted.reverse,
        );
        let expected = (Atom::num(Rational::from(3) - &epsilon) - Atom::num(n))
            / (Atom::num(n - 1) * (Atom::num(1) + Atom::var(eta)));
        assert!((actual - expected).together().cancel().is_zero());
    }
    assert_eq!(bank.generated, 1);
    assert_eq!(bank.hits, 11);
}
