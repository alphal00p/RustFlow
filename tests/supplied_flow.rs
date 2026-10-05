use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::engine::SuppliedAuxiliarySystem;
use symbolica_amflow::reduction::ReducedSystem;
use symbolica_amflow::*;

fn problem() -> (IntegralFamily, Vec<Integral>, SuppliedAuxiliarySystem) {
    let epsilon = symbol!("supplied_test::epsilon");
    let eta = symbol!("supplied_test::eta");
    let gram = vec![vec![Atom::num(-1)]];
    let family = IntegralFamily {
        name: "two routed tadpoles".into(),
        loops: vec!["l".into()],
        external: vec!["p".into()],
        external_gram: gram.clone(),
        propagators: vec![
            Propagator::quadratic(&[1], &[0], Atom::num(2), &gram).unwrap(),
            Propagator::quadratic(&[1], &[1], Atom::num(3), &gram).unwrap(),
        ],
        physical_propagators: 2,
        epsilon,
        dimension: 4,
    };
    let targets = vec![
        Integral(vec![1, 0]),
        Integral(vec![0, 1]),
        Integral(vec![2, 0]),
    ];
    let derivative = (Atom::num(1) - Atom::var(epsilon)) / (Atom::num(2) + Atom::var(eta));
    let reduced = ReducedSystem {
        basis: targets[..2].to_vec(),
        matrix: vec![
            vec![derivative.clone(), Atom::new()],
            vec![Atom::new(), Atom::new()],
        ],
        targets: vec![
            BTreeMap::from([(targets[0].clone(), Atom::num(1))]),
            BTreeMap::from([(targets[1].clone(), Atom::num(1))]),
            BTreeMap::from([(targets[0].clone(), derivative)]),
        ],
        nonzero_conditions: vec![],
        candidates: BTreeMap::new(),
        transformations: vec![],
    };
    let supplied = SuppliedAuxiliarySystem {
        variable: eta,
        reduced,
        deformation_mask: vec![true, false],
        provenance: "analytic ordinary tadpoles; d_eta I(a) = a I(a+1)".into(),
    };
    (family, targets, supplied)
}

fn prepare(
    family: &IntegralFamily,
    targets: &[Integral],
    supplied: SuppliedAuxiliarySystem,
    options: &FlowOptions,
) -> Result<PreparedFlow> {
    PreparedFlow::from_supplied(
        family,
        targets,
        &KinematicPoint::default(),
        supplied,
        options,
        &RunContext::default(),
    )
}

#[test]
fn partial_deformation_stays_fixed_and_targets_precede_endpoint_extraction() {
    let (family, targets, supplied) = problem();
    // Both construction and provider options intentionally request all lines;
    // the supplied connection explicitly shifts just the first physical line.
    let options = FlowOptions {
        mass_mode: MassMode::All,
        ..Default::default()
    };
    let flow = prepare(&family, &targets, supplied, &options).unwrap();
    assert_eq!(flow.deformation_mask(), &[true, false]);
    assert!(
        flow.supplied_provenance()
            .unwrap()
            .contains("analytic ordinary tadpoles")
    );
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let provider = recursive::RecursiveBoundary::new(&backend, &options, &context);
    let epsilon = Rational::from((1, 10));
    let values = flow
        .evaluate(&epsilon, &options, &provider, &context)
        .unwrap();
    let p = Precision::decimal(60).unwrap();
    let dimension = Rational::from(4) - &Rational::from(2) * &epsilon;
    for (value, (power, mass)) in values.iter().zip([(1, 2), (1, 3), (2, 2)]) {
        let expected = vacuum::tadpole(power, &p.i(mass), &dimension, p).unwrap();
        assert!(p.close(value, &expected, 20), "{value} != {expected}");
    }
    assert!(matches!(
        flow.evaluate(&epsilon, &options, &boundary::OneLoopBoundary, &context),
        Err(Error::Unsupported(message)) if message.contains("all physical propagators"),
    ));
}

#[test]
fn malformed_declared_systems_fail_before_numerical_work() {
    let (family, targets, supplied) = problem();
    let options = FlowOptions::default();
    let mut malformed = supplied.clone();
    malformed.deformation_mask.pop();
    assert!(matches!(
        prepare(&family, &targets, malformed, &options),
        Err(Error::InvalidInput(_))
    ));
    let mut malformed = supplied.clone();
    malformed.deformation_mask.fill(false);
    assert!(matches!(
        prepare(&family, &targets, malformed, &options),
        Err(Error::InvalidInput(_))
    ));
    let mut malformed = supplied.clone();
    malformed.reduced.basis[1] = malformed.reduced.basis[0].clone();
    assert!(matches!(
        prepare(&family, &targets, malformed, &options),
        Err(Error::InvalidInput(_))
    ));
    let mut malformed = supplied.clone();
    malformed.reduced.matrix[0].pop();
    assert!(matches!(
        prepare(&family, &targets, malformed, &options),
        Err(Error::InvalidInput(_))
    ));
    let mut malformed = supplied.clone();
    malformed.reduced.targets.pop();
    assert!(matches!(
        prepare(&family, &targets, malformed, &options),
        Err(Error::InvalidInput(_))
    ));
    let mut malformed = supplied.clone();
    malformed.reduced.targets[0].insert(Integral(vec![2, 2]), Atom::num(1));
    assert!(matches!(
        prepare(&family, &targets, malformed, &options),
        Err(Error::IncompleteReduction(_))
    ));
    let mut malformed = supplied.clone();
    malformed.reduced.nonzero_conditions.push(Atom::new());
    assert!(matches!(
        prepare(&family, &targets, malformed, &options),
        Err(Error::Reduction(_))
    ));
    let mut malformed = supplied.clone();
    malformed.provenance.clear();
    assert!(matches!(
        prepare(&family, &targets, malformed, &options),
        Err(Error::InvalidInput(_))
    ));
    let mut malformed = supplied.clone();
    malformed.variable = family.epsilon;
    assert!(matches!(
        prepare(&family, &targets, malformed, &options),
        Err(Error::InvalidInput(_))
    ));
    let options = FlowOptions {
        dimension: 6,
        ..options
    };
    assert!(matches!(
        prepare(&family, &targets, supplied, &options),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn exact_specialization_is_simultaneous_and_applied_only_to_family() {
    let (mut family, targets, mut supplied) = problem();
    let mass = Atom::var(symbol!("supplied_test::mass"));
    let intermediary = Atom::var(symbol!("supplied_test::intermediary"));
    family.propagators[0].constant = -mass.clone();
    // The family needs a single a -> b, b -> 2 substitution. Supplying a -> b
    // leaves a nonnumeric coefficient, which the constructor must reject.
    let point = KinematicPoint(BTreeMap::from([
        (mass.clone(), intermediary.clone()),
        (intermediary, Atom::num(2)),
    ]));
    assert!(matches!(
        PreparedFlow::from_supplied(
            &family,
            &targets,
            &point,
            supplied.clone(),
            &FlowOptions::default(),
            &RunContext::default()
        ),
        Err(Error::InvalidInput(_))
    ));
    let point = KinematicPoint(BTreeMap::from([
        (mass, Atom::num(2)),
        (Atom::var(supplied.variable), Atom::num(42)),
    ]));
    let flow = PreparedFlow::from_supplied(
        &family,
        &targets,
        &point,
        supplied.clone(),
        &FlowOptions::default(),
        &RunContext::default(),
    )
    .unwrap();
    assert_eq!(flow.system.matrix, supplied.reduced.matrix);
    supplied.reduced.matrix[0][0] = Atom::var(symbol!("supplied_test::unresolved"));
    assert!(matches!(
        PreparedFlow::from_supplied(
            &family,
            &targets,
            &point,
            supplied,
            &FlowOptions::default(),
            &RunContext::default()
        ),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn public_mutation_and_vanishing_nonzero_conditions_are_rejected() {
    let (family, targets, supplied) = problem();
    let options = FlowOptions::default();
    let context = RunContext::default();
    let epsilon = Rational::from((1, 10));
    let mut flow = prepare(&family, &targets, supplied.clone(), &options).unwrap();
    flow.system.matrix[0][0] = Atom::num(1);
    assert!(
        matches!(flow.evaluate(&epsilon, &options, &boundary::OneLoopBoundary, &context), Err(Error::InvalidInput(message)) if message.contains("modified"))
    );
    let mut supplied = supplied;
    supplied
        .reduced
        .nonzero_conditions
        .push(Atom::var(family.epsilon) - Atom::num((1, 10)));
    let flow = prepare(&family, &targets, supplied, &options).unwrap();
    assert!(matches!(
        flow.evaluate(&epsilon, &options, &boundary::OneLoopBoundary, &context),
        Err(Error::Reduction(_))
    ));
}

#[test]
fn supplied_coefficients_require_rational_exact_arithmetic() {
    let (family, targets, supplied) = problem();
    let options = FlowOptions::default();
    for coefficient in [
        Atom::var(supplied.variable).pow(Atom::num((1, 2))),
        Atom::var(symbol!("supplied_test::unresolved")),
        Atom::num(Precision::decimal(40).unwrap().real(3)),
    ] {
        let mut malformed = supplied.clone();
        malformed.reduced.matrix[0][0] = coefficient;
        assert!(prepare(&family, &targets, malformed, &options).is_err());
    }
}

#[test]
fn eta_dependent_guards_are_specialized_exactly_even_for_zero_outputs() {
    let (family, targets, mut supplied) = problem();
    let eta = Atom::var(supplied.variable);
    supplied
        .reduced
        .nonzero_conditions
        .push((&eta + Atom::num(1)) * (Atom::var(family.epsilon) - Atom::num((1, 10))));
    let options = FlowOptions::default();
    for empty in [false, true] {
        let mut input = supplied.clone();
        if empty {
            input.reduced.basis.clear();
            input.reduced.matrix.clear();
            input.reduced.targets.iter_mut().for_each(BTreeMap::clear);
        }
        let flow = prepare(&family, &targets, input, &options).unwrap();
        assert!(matches!(
            flow.evaluate(&Rational::from((1, 10)), &options, &boundary::OneLoopBoundary, &RunContext::default()),
            Err(Error::Reduction(message)) if message.contains("vanishes identically"),
        ));
    }
}

#[test]
fn retained_guard_detours_preserve_constant_integrals_and_frobenius_endpoint() {
    let (family, targets, mut supplied) = problem();
    let integral = targets[1].clone(); // This tadpole uses the unshifted propagator.
    supplied.reduced.basis = vec![integral.clone()];
    supplied.reduced.matrix = vec![vec![Atom::new()]];
    supplied.reduced.targets = vec![BTreeMap::from([(integral.clone(), Atom::one())])];
    let eta = Atom::var(supplied.variable);
    let imaginary = Atom::num(Complex::new(Rational::from(0), Rational::from(1)));
    supplied.reduced.nonzero_conditions = vec![
        // The denominator's zero at -i survives cancellation of this ratio.
        (eta.clone().pow(2) + Atom::one()) / (&eta + imaginary),
        // Zero is excluded for ordinary continuation but remains a valid
        // Frobenius limiting endpoint, with an unambiguous constant solution.
        eta,
    ];
    let epsilon = Rational::from((1, 10));
    let p = Precision::decimal(60).unwrap();
    let expected = vacuum::tadpole(1, &p.i(3), &Rational::from((19, 5)), p).unwrap();
    for local_coordinate in [
        local_coordinates::LocalCoordinate::Identity,
        local_coordinates::LocalCoordinate::BalancedMobius,
    ] {
        let options = FlowOptions {
            local_coordinate,
            ..Default::default()
        };
        let flow = prepare(
            &family,
            std::slice::from_ref(&integral),
            supplied.clone(),
            &options,
        )
        .unwrap();
        let backend = RustRedBackend::default();
        let context = RunContext::default();
        let provider = recursive::RecursiveBoundary::new(&backend, &options, &context);
        let value = flow
            .evaluate(&epsilon, &options, &provider, &context)
            .unwrap();
        assert!(p.close(&value[0], &expected, 20));
    }
}
