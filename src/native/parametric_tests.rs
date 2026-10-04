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

#[test]
fn budget_exhausted_domain_preserves_exact_rules_and_fallbacks() {
    use crate::ReductionBackend;

    let eta = symbol!("partial_domain_bubble_eta");
    let family = crate::IntegralFamily {
        name: "partial_domain_bubble".into(),
        loops: vec!["l".into()],
        external: vec!["p".into()],
        external_gram: vec![vec![Atom::num(-2)]],
        propagators: vec![
            crate::Propagator {
                constant: -Atom::num(1) - Atom::var(eta),
                scalar_products: vec![Atom::num(1), Atom::new()],
            },
            crate::Propagator {
                constant: -Atom::num(5) - Atom::var(eta),
                scalar_products: vec![Atom::num(1), Atom::num(2)],
            },
        ],
        physical_propagators: 2,
        epsilon: symbol!("partial_domain_bubble_eps"),
        dimension: 4,
    };
    let epsilon = Rational::from((1, 10));
    let converted = family.convert_at_epsilon(Some(&epsilon)).unwrap();
    let sources = SourceSystem::<2>::from_family_with_lorentz(&converted.family, false).unwrap();
    let solver = SectorSolver::new(&sources, [true; 2], SectorConfig::default()).unwrap();
    let bounded = solver.solve_domains(
        vec![CoordinateCase::new([None; 2]).unwrap().into()],
        SectorSolveOptions {
            symbolic: SearchOptions {
                max_depth: Some(3),
                ..Default::default()
            },
            numerical_depth: 3,
            max_symbolic_cases: Some(1),
            ..Default::default()
        },
    );
    assert!(matches!(bounded, Err(SectorSolveError::CaseBudget { .. })));
    let partial = generic_rules(&solver, 3, 1);
    assert_eq!(partial.len(), 1);
    assert!(!partial[0].candidate.sources.is_empty());

    // A one-case budget deliberately leaves exceptional faces uncovered.
    // Find one covered target and one genuinely reducible excluded target.
    let mut covered = None;
    let mut excluded = None;
    let mut fixed = Bank::new(&sources, &solver);
    for a in 1..=4 {
        for b in 1..=4 {
            let target = [a, b];
            if apply(&partial[0], target, &sources, &solver)
                .unwrap()
                .is_some()
            {
                covered.get_or_insert(target);
            } else if fixed
                .try_reduce(target, &sources, &solver, 3)
                .unwrap()
                .is_some()
            {
                excluded.get_or_insert(target);
            }
        }
    }
    let covered = covered.expect("the admitted rule did not cover a test target");
    let excluded = excluded.expect("the test did not exercise exceptional-face fallback");
    let mut bank = Bank::new(&sources, &solver);
    // Inject the deliberately partial two-line domain to exercise exactly the
    // same dispatch/fallback path as production's four-line generic tier.
    bank.generic = Some(partial);
    let direct = bank
        .try_reduce(covered, &sources, &solver, 3)
        .unwrap()
        .unwrap();
    assert_eq!((bank.domain_hits, bank.generated), (1, 0));
    let fallback = bank
        .try_reduce(excluded, &sources, &solver, 3)
        .unwrap()
        .unwrap();
    assert_eq!((bank.domain_hits, bank.generated), (1, 1));
    assert!(
        bank.try_reduce([1, 1], &sources, &solver, 3)
            .unwrap()
            .is_none()
    );
    assert_eq!(bank.missed, 1);
    assert!(direct.conditions.iter().any(|p| !p.is_constant()));

    // Independently reduce both identities with the default concrete backend.
    // Its rules need not have the same RHS or choose the same local equation.
    let identities = [(covered, direct), (excluded, fallback)];
    let mut targets = std::collections::BTreeSet::new();
    for (lhs, applied) in &identities {
        targets.insert(crate::Integral(lhs.to_vec()));
        for term in &applied.terms {
            targets.insert(crate::Integral(
                term.integral.powers().iter().map(|p| p.value()).collect(),
            ));
        }
    }
    let reduction = crate::RustRedBackend::default()
        .reduce_at_epsilon(
            &family,
            &targets.into_iter().collect::<Vec<_>>(),
            &epsilon,
            &crate::RunContext::default(),
        )
        .unwrap();
    for (lhs, applied) in identities {
        let mut residual = reduction.expand(&crate::Integral(lhs.to_vec())).unwrap();
        for term in applied.terms {
            let coefficient =
                crate::family::substitute(&term.coefficient.to_expression(), &converted.reverse);
            let child = crate::Integral(term.integral.powers().iter().map(|p| p.value()).collect());
            for (master, reduced) in reduction.expand(&child).unwrap() {
                let old = residual.remove(&master).unwrap_or_default();
                residual.insert(master, old - &coefficient * reduced);
            }
        }
        assert!(residual.values().all(|c| c.together().cancel().is_zero()));
    }
}

#[test]
fn four_line_domains_are_reused_while_three_line_sectors_keep_rays() {
    let sources = setup();
    let four = [true, true, false, false, true, false, true, false, false];
    let solver = SectorSolver::new(&sources, four, SectorConfig::default()).unwrap();
    let mut bank = Bank::new(&sources, &solver);
    for rank in [3, 4] {
        assert!(
            bank.try_reduce([2, 2, -rank, 0, 2, 0, 2, 0, 0], &sources, &solver, 3)
                .unwrap()
                .is_some()
        );
    }
    assert_eq!((bank.domains, bank.domain_hits, bank.generated), (1, 2, 0));

    let solver = SectorSolver::new(&sources, sector(), SectorConfig::default()).unwrap();
    let mut bank = Bank::new(&sources, &solver);
    assert!(
        bank.try_reduce([-3, 2, 0, 0, 2, 0, 2, 0, 0], &sources, &solver, 3)
            .unwrap()
            .is_some()
    );
    assert_eq!((bank.domains, bank.domain_hits, bank.generated), (0, 0, 1));
}
