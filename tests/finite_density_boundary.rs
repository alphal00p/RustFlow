use ahash::HashMap;
use symbolica::prelude::*;
use symbolica_amflow::finite_density::boundary::{
    IntegratedOccupiedBoundary, OccupiedBoundaryDistribution, OccupiedBoundaryLimits,
};
use symbolica_amflow::finite_density::compact::CompactShell;
use symbolica_amflow::*;

fn family() -> IntegralFamily {
    let gram = vec![vec![Atom::num(1)]];
    IntegralFamily {
        name: "occupied_boundary_sunset".into(),
        loops: vec!["q1".into(), "q2".into()],
        external: vec!["u".into()],
        external_gram: gram.clone(),
        propagators: vec![
            Propagator::quadratic(&[1, 0], &[0], Atom::num((1, 4)), &gram).unwrap(),
            Propagator::quadratic(&[0, 1], &[0], Atom::num((1, 4)), &gram).unwrap(),
            Propagator::quadratic(&[1, -1], &[0], Atom::num(1), &gram).unwrap(),
        ],
        physical_propagators: 3,
        epsilon: symbol!("occupied_boundary_eps"),
        dimension: 4,
    }
}

fn distribution(source_loop_index: usize) -> OccupiedBoundaryDistribution {
    OccupiedBoundaryDistribution {
        source_loop_index,
        shell: CompactShell {
            mass_squared: Rational::from((1, 4)),
            chemical_potential: Rational::from(1),
        },
        cut_index: 1,
        upper_index: 0,
        lower_index: 0,
    }
}

fn identity_region(hard: [bool; 2]) -> regions::LoopRegion {
    regions::LoopRegion {
        transformation: vec![
            vec![Atom::num(1), Atom::num(0)],
            vec![Atom::num(0), Atom::num(1)],
        ],
        hard: hard.to_vec(),
        hard_branches: vec![],
        jacobian_determinant: Atom::num(1),
    }
}

fn expected(expression: &str, p: Precision) -> ComplexFloat {
    p.eval(
        &Atom::parse(
            expression,
            "occupied_boundary_reference",
            Default::default(),
        )
        .unwrap(),
        &HashMap::default(),
    )
    .unwrap()
}

#[test]
fn inverse_energy_moments_require_explicit_positive_shell_admission() {
    let family = IntegralFamily {
        name: "compact_energy_laurent".into(),
        loops: vec!["q".into()],
        external: vec!["u".into()],
        external_gram: vec![vec![Atom::one()]],
        propagators: vec![],
        physical_propagators: 0,
        epsilon: symbol!("compact_energy_laurent_eps"),
        dimension: 4,
    };
    let coordinates = [parse!("laurent_g"), parse!("laurent_e")];
    let backend = RustRedBackend {
        bubble_subloops: false,
        ..Default::default()
    };
    let options = FlowOptions::default();
    let context = RunContext::default();
    let plain = IntegratedOccupiedBoundary::new(
        &backend,
        &options,
        &context,
        OccupiedBoundaryLimits::default(),
    )
    .unwrap();
    let enabled = IntegratedOccupiedBoundary::new(
        &backend,
        &options,
        &context,
        OccupiedBoundaryLimits::default(),
    )
    .unwrap()
    .with_positive_compact_energy_powers(true);
    let project = |expression: Atom| {
        integrand::projected_factor_region(&expression, &coordinates, &family, &[false], 100)
            .unwrap()
    };
    let inverse = project(coordinates[1].pow(-1));
    let epsilon = Rational::from((1, 2));
    let p = Precision::decimal(60).unwrap();
    assert!(matches!(
        plain.evaluate_projected(
            &inverse,
            &[distribution(0)],
            &epsilon,
            &HashMap::default(),
            p
        ),
        Err(Error::Unsupported(_))
    ));
    for (power, upper, expression, reference) in [
        (1, 0, coordinates[1].pow(-1), "log(2)/pi^(1/2)"),
        (1, 0, coordinates[1].pow(-2), "1/pi^(1/2)"),
        (2, 0, coordinates[1].pow(-1), "-2/pi^(1/2)"),
        (1, 1, coordinates[1].pow(-1), "1/pi^(1/2)"),
        // g*C2=C1+m²*C2 still holds with the original inverse-energy weight.
        (
            2,
            0,
            &coordinates[0] / &coordinates[1],
            "(log(2)-1/2)/pi^(1/2)",
        ),
    ] {
        let mut shell = distribution(0);
        shell.cut_index = power;
        shell.upper_index = upper;
        let value = enabled
            .evaluate_projected(
                &project(expression),
                &[shell],
                &epsilon,
                &HashMap::default(),
                p,
            )
            .unwrap();
        assert!(p.close(&value, &expected(reference, p), 40));
    }
    let mut massless = distribution(0);
    massless.shell.mass_squared = Rational::zero();
    assert!(matches!(
        enabled.evaluate_projected(&inverse, &[massless], &epsilon, &HashMap::default(), p),
        Err(Error::Unsupported(_))
    ));
    for expression in [
        coordinates[0].pow(-1),
        (&coordinates[1] + Atom::one()).pow(-1),
        (&coordinates[1] - Atom::one()).pow(-1),
    ] {
        assert!(matches!(
            enabled.evaluate_projected(
                &project(expression),
                &[distribution(0)],
                &epsilon,
                &HashMap::default(),
                p
            ),
            Err(Error::Unsupported(_))
        ));
    }
}

#[test]
fn actual_sunset_mixed_and_fully_occupied_region_coefficients_are_integrated() {
    let family = family();
    let backend = RustRedBackend {
        bubble_subloops: false,
        ..Default::default()
    };
    let options = FlowOptions::default();
    let context = RunContext::default();
    let boundary = IntegratedOccupiedBoundary::new(
        &backend,
        &options,
        &context,
        OccupiedBoundaryLimits::default(),
    )
    .unwrap();
    // D=3, d=2: raw native compact C_1 moment=(mu-m)/sqrt(pi).
    // The virtual hard factor is a uniform massive Gaussian/tadpole seed.
    let epsilon = Rational::from((1, 2));
    let p = Precision::decimal(60).unwrap();
    let mixed = regions::expand_region(
        &family,
        &Integral(vec![0, 1, 1]),
        &[false, true, true],
        &identity_region([false, true]),
        0,
    )
    .unwrap();
    assert!(
        (mixed.eta_power + Atom::var(family.epsilon))
            .together()
            .cancel()
            .is_zero()
    );
    let projected = integrand::projected_factor_region(
        &mixed.coefficients[0],
        &mixed.coordinates,
        &family,
        &[false, true],
        100,
    )
    .unwrap();
    let value = boundary
        .evaluate_projected(
            &projected,
            &[distribution(0)],
            &epsilon,
            &HashMap::default(),
            p,
        )
        .unwrap();
    assert!(p.close(&value, &p.scale(&p.i(1), 1, 2), 40), "{value}");

    // The all-soft expansion of the same one-cut integrand contains an
    // unrestricted polynomial virtual loop and is a measure-aware DR zero.
    let soft = regions::expand_region(
        &family,
        &Integral(vec![0, 1, 1]),
        &[false, true, true],
        &identity_region([false, false]),
        2,
    )
    .unwrap();
    for coefficient in &soft.coefficients {
        let projected = integrand::projected_factor_region(
            coefficient,
            &soft.coordinates,
            &family,
            &[false, false],
            100,
        )
        .unwrap();
        let report = boundary
            .evaluate_projected_with_provenance(
                &projected,
                &[distribution(0)],
                &epsilon,
                &HashMap::default(),
                p,
            )
            .unwrap();
        assert_eq!(report.value, p.zero());
        assert_eq!(report.integrated_products, 0);
        assert_eq!(report.scaleless_noncompact_products, projected.terms.len());
    }

    let full = regions::expand_region(
        &family,
        &Integral(vec![0, 0, 1]),
        &[false, false, true],
        &identity_region([false, false]),
        2,
    )
    .unwrap();
    assert_eq!(full.eta_power, Atom::num(-1));
    let shells = [distribution(0), distribution(1)];
    for (index, reference) in [(0, "-1/(4*pi)"), (1, "0"), (2, "13/(32*pi)")] {
        let projected = integrand::projected_factor_region(
            &full.coefficients[index],
            &full.coordinates,
            &family,
            &[false, false],
            100,
        )
        .unwrap();
        let value = boundary
            .evaluate_projected(&projected, &shells, &epsilon, &HashMap::default(), p)
            .unwrap();
        assert!(
            p.close(&value, &expected(reference, p), 40),
            "half-order {index}: {value}"
        );
    }
}

#[test]
fn compact_tensor_projection_and_raised_shell_keep_the_original_numerator() {
    let family = family();
    let coordinates = vec![
        parse!("occ_g11"),
        parse!("occ_g12"),
        parse!("occ_g22"),
        parse!("occ_e1"),
        parse!("occ_e2"),
    ];
    let backend = RustRedBackend {
        bubble_subloops: false,
        ..Default::default()
    };
    let options = FlowOptions::default();
    let context = RunContext::default();
    let boundary = IntegratedOccupiedBoundary::new(
        &backend,
        &options,
        &context,
        OccupiedBoundaryLimits::default(),
    )
    .unwrap();
    let epsilon = Rational::from((1, 2));
    for precision in [40, 60] {
        let p = Precision::decimal(precision).unwrap();
        let projected = integrand::projected_factor_region(
            &parse!("occ_g12^2"),
            &coordinates,
            &family,
            &[false, false],
            100,
        )
        .unwrap();
        let value = boundary
            .evaluate_projected(
                &projected,
                &[distribution(0), distribution(1)],
                &epsilon,
                &HashMap::default(),
                p,
            )
            .unwrap();
        // Independent angular average: (E1 E2)^2 + r1^2 r2^2/d.
        assert!(p.close(&value, &expected("19/(192*pi)", p), 28));
    }

    let mut one = family.clone();
    one.loops.truncate(1);
    let coordinates = vec![parse!("occ_g11"), parse!("occ_e1")];
    let p = Precision::decimal(60).unwrap();
    let mut shell = distribution(0);
    for (expression, cut, upper, lower, reference) in [
        ("1", 1, 0, 0, "1/(2*sqrt(pi))"),
        ("1", 2, 0, 0, "-1/sqrt(pi)"),
        // q^2 C_2=(q^2-m^2)C_2+m^2 C_2=C_1+m^2 C_2.
        // Substituting q^2=m^2 before raising would miss the C_1 term.
        ("occ_g11", 2, 0, 0, "1/(4*sqrt(pi))"),
        ("1", 1, 1, 0, "1/sqrt(pi)"),
        ("occ_g11", 1, 1, 0, "1/(4*sqrt(pi))"),
        ("1", 2, 1, 0, "0"),
        ("1", 1, 0, 1, "0"),
    ] {
        shell.cut_index = cut;
        shell.upper_index = upper;
        shell.lower_index = lower;
        let polynomial = parse!(expression);
        let projected =
            integrand::projected_factor_region(&polynomial, &coordinates, &one, &[false], 100)
                .unwrap();
        let value = boundary
            .evaluate_projected(
                &projected,
                &[shell.clone()],
                &epsilon,
                &HashMap::default(),
                p,
            )
            .unwrap();
        assert!(
            p.close(&value, &expected(reference, p), 38),
            "{expression}, C{cut}, H{upper}, H{lower}: {value}"
        );
    }
}

#[test]
fn retained_virtual_soft_factors_and_hard_occupied_loops_fail_explicitly() {
    let family = family();
    let coordinates = vec![
        parse!("occ_g11"),
        parse!("occ_g12"),
        parse!("occ_g22"),
        parse!("occ_e1"),
        parse!("occ_e2"),
    ];
    let backend = RustRedBackend {
        bubble_subloops: false,
        ..Default::default()
    };
    let options = FlowOptions::default();
    let insufficient_guard = FlowOptions {
        guard_digits: 7,
        ..options.clone()
    };
    let context = RunContext::default();
    assert!(matches!(
        IntegratedOccupiedBoundary::new(
            &backend,
            &insufficient_guard,
            &context,
            OccupiedBoundaryLimits::default()
        ),
        Err(Error::InvalidInput(_))
    ));
    let boundary = IntegratedOccupiedBoundary::new(
        &backend,
        &options,
        &context,
        OccupiedBoundaryLimits::default(),
    )
    .unwrap();
    let epsilon = Rational::from((1, 2));
    let p = Precision::decimal(40).unwrap();
    let projected = integrand::projected_factor_region(
        &parse!("1/(occ_g12-1)"),
        &coordinates,
        &family,
        &[false, false],
        100,
    )
    .unwrap();
    assert!(matches!(
        boundary.evaluate_projected(
            &projected,
            &[distribution(0), distribution(1)],
            &epsilon,
            &HashMap::default(),
            p
        ),
        Err(Error::Unsupported(_))
    ));
    assert!(matches!(
        boundary.evaluate_projected(
            &projected,
            &[distribution(0)],
            &epsilon,
            &HashMap::default(),
            p
        ),
        Err(Error::Unsupported(_))
    ));
    let projected = integrand::projected_factor_region(
        &Atom::num(1),
        &coordinates,
        &family,
        &[true, false],
        100,
    )
    .unwrap();
    assert!(matches!(
        boundary.evaluate_projected(
            &projected,
            &[distribution(0), distribution(1)],
            &epsilon,
            &HashMap::default(),
            p
        ),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn raised_compact_boundary_seeds_obey_native_multiplication_and_endpoint_constraints() {
    let family = family();
    let coordinates = vec![
        parse!("constraint_g11"),
        parse!("constraint_g12"),
        parse!("constraint_g22"),
        parse!("constraint_e1"),
        parse!("constraint_e2"),
    ];
    let backend = RustRedBackend {
        bubble_subloops: false,
        ..Default::default()
    };
    let options = FlowOptions::default();
    let context = RunContext::default();
    let boundary = IntegratedOccupiedBoundary::new(
        &backend,
        &options,
        &context,
        OccupiedBoundaryLimits::default(),
    )
    .unwrap();
    let p = Precision::decimal(70).unwrap();
    let epsilon = Rational::from((4, 5));
    let evaluate = |expression: Atom, shells: &[OccupiedBoundaryDistribution]| {
        let projected = integrand::projected_factor_region(
            &expression,
            &coordinates,
            &family,
            &[false, false],
            100,
        )
        .unwrap();
        boundary
            .evaluate_projected(&projected, shells, &epsilon, &HashMap::default(), p)
            .unwrap()
    };
    let mut shells = [distribution(0), distribution(1)];
    shells[1].upper_index = 1;
    // These are distribution identities, independently of the moment formula:
    // (q²-m²) C_n=C_(n-1), and E H_s=mu H_s-H_(s-1) for s>=2.
    // The H1 right-hand correction vanishes because h delta(h)=0.
    for cut in 1..=3 {
        for upper in 1..=3 {
            shells[0].cut_index = cut;
            shells[0].upper_index = upper;
            let value = evaluate(Atom::one(), &shells);
            let energy = evaluate(coordinates[3].clone(), &shells);
            let mut lower_upper = shells.clone();
            lower_upper[0].upper_index -= 1;
            let correction = if upper == 1 {
                p.zero()
            } else {
                evaluate(Atom::one(), &lower_upper)
            };
            assert!(p.close(&energy, &p.sub(&value, &correction), 40));
            let inverse_shell = evaluate(&coordinates[0] - Atom::num((1, 4)), &shells);
            let mut lower_cut = shells.clone();
            lower_cut[0].cut_index -= 1;
            let correction = if cut == 1 {
                p.zero()
            } else {
                evaluate(Atom::one(), &lower_cut)
            };
            assert!(p.close(&inverse_shell, &correction, 40));
        }
    }
    // Actual constant rows of the native 64-label fully occupied sunset
    // system. At eps=4/5 its eta^-2 forcing is
    // (I37-I45)/4 + 5*(I49-I58)/16; its eta^-4 forcing is 5/12 of this.
    shells[0].cut_index = 1;
    shells[0].upper_index = 1;
    shells[1].cut_index = 2;
    shells[1].upper_index = 1;
    let i37 = evaluate(coordinates[4].clone(), &shells);
    let i45 = evaluate(Atom::one(), &shells);
    shells[0].upper_index = 2;
    let i49 = evaluate(Atom::one(), &shells);
    shells.swap(0, 1);
    shells[0].source_loop_index = 0;
    shells[1].source_loop_index = 1;
    let i58 = evaluate(Atom::one(), &shells);
    assert!(p.close(&i37, &i45, 40));
    assert!(p.close(&i49, &i58, 40));
    let principal = p.add(
        &p.scale(&p.sub(&i37, &i45), 1, 4),
        &p.scale(&p.sub(&i49, &i58), 5, 16),
    );
    assert!(p.close(&principal, &p.zero(), 40));
}
