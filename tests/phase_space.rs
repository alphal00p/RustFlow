use symbolica::prelude::*;
use symbolica_amflow::cuts::*;
use symbolica_amflow::phase_space::*;
use symbolica_amflow::*;

fn family(s: Rational, reverse: bool, scale: i64) -> CutFamily {
    let gram = vec![vec![Atom::num(s)]];
    let sign = if reverse { -1 } else { 1 };
    let family = IntegralFamily {
        name: "physical_two_body_terminal".into(),
        loops: vec!["k".into()],
        external: vec!["p".into()],
        external_gram: gram.clone(),
        propagators: vec![
            Propagator::quadratic(&[scale], &[0], Atom::num(1), &gram).unwrap(),
            Propagator::quadratic(&[-scale], &[1], Atom::num(4), &gram).unwrap(),
        ],
        physical_propagators: 2,
        epsilon: symbol!("phase_space_test::eps"),
        dimension: 4,
    };
    CutFamily::new(
        family,
        CutMetadata::new(
            [
                CutLine {
                    propagator: 0,
                    definition: CutDefinition::PositiveEnergy {
                        momentum: MomentumRouting {
                            loops: vec![(sign * scale).into()],
                            external: vec![0.into()],
                        },
                    },
                },
                CutLine {
                    propagator: 1,
                    definition: CutDefinition::PositiveEnergy {
                        momentum: MomentumRouting {
                            loops: vec![(-sign * scale).into()],
                            external: vec![sign.into()],
                        },
                    },
                },
            ],
            vec![LoopPrescription::Insensitive],
        )
        .unwrap(),
    )
    .unwrap()
}
fn channel() -> FutureTimelikeChannel {
    FutureTimelikeChannel {
        external: vec![1.into()],
    }
}
fn targets() -> Vec<Integral> {
    vec![
        Integral(vec![1, 1]),
        Integral(vec![2, 1]),
        Integral(vec![1, 2]),
        Integral(vec![0, 1]),
    ]
}
fn relative_close(p: Precision, actual: &ComplexFloat, expected: &ComplexFloat, digits: u32) {
    let scale = p.norm(expected);
    let scale = if scale == p.real(0) { p.real(1) } else { scale };
    assert!(
        p.norm(&p.sub(actual, expected)) <= p.tolerance(digits) * scale,
        "{actual} != {expected}"
    );
}

#[test]
fn native_two_body_laurent_coefficients_match_analytic_volume_and_mass_derivatives() {
    let family = family(25.into(), false, 1);
    let result = solve_two_body_phase_space(
        &family,
        &channel(),
        &targets(),
        &KinematicPoint::default(),
        1,
        &FlowOptions::default(),
        &RustRedBackend::default(),
        &RunContext::default(),
    )
    .unwrap();
    let p = Precision::decimal(80).unwrap();
    let pi = ComplexFloat::new(p.real(1).pi(), p.real(0));
    let volume = p.div(
        &p.pow(&p.i(6), &p.rational(&Rational::from((1, 2)))),
        &p.scale(&pi, 25, 1),
    );
    let euler = p
        .parse(
            "0.5772156649015328606065120900824024310421593359399235988057672348848677",
            "0",
        )
        .unwrap();
    let logarithm = p.add(&p.sub(&p.i(2), &euler), &p.log(&p.scale(&pi, 100, 384)));
    for (index, factor) in [
        Rational::from(1),
        Rational::from((-7, 96)),
        Rational::from((-11, 192)),
    ]
    .iter()
    .enumerate()
    {
        let expansion = &result[index];
        assert_eq!(expansion.verified_digits, Some(20));
        assert!(expansion.validation_samples > 0 && expansion.refinements > 0);
        let finite = p.mul(&p.rational(factor), &volume);
        relative_close(p, &expansion.coefficients[&0], &finite, 20);
        let log = if index == 0 {
            logarithm.clone()
        } else {
            p.sub(&logarithm, &p.i(2))
        };
        relative_close(p, &expansion.coefficients[&1], &p.mul(&finite, &log), 20);
    }
    for coefficient in result[3].coefficients.values() {
        assert_eq!(*coefficient, p.zero());
    }
}

#[test]
fn phase_space_samples_match_independent_automatic_uncut_discontinuities() {
    let family = family(25.into(), false, 1);
    let targets = targets();
    let options = FlowOptions::default();
    let context = RunContext::default();
    let backend = RustRedBackend::default();
    let cut = PreparedTwoBodyPhaseSpace::new(
        &family,
        &channel(),
        &targets,
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    // This second computation uses the ordinary auxiliary-mass DE and its
    // recursively/analytically constructed infinity boundary, not this terminal.
    let uncut = PreparedFlow::new(
        family.family(),
        &targets[..3],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    for epsilon in [Rational::from((1, 13)), Rational::from((1, 17))] {
        let actual = cut.evaluate(&epsilon, &options, &context).unwrap();
        let expected = uncut
            .evaluate(&epsilon, &options, &boundary::OneLoopBoundary, &context)
            .unwrap();
        let p = Precision::decimal(options.digits + options.guard_digits).unwrap();
        let pi = ComplexFloat::new(p.real(1).pi(), p.real(0));
        let dimension = Rational::from(4) - Rational::from(2) * &epsilon;
        let factor = p.scale(
            &p.pow(
                &p.scale(&pi, 4, 1),
                &p.rational(&(-dimension / Rational::from(2))),
            ),
            2,
            1,
        );
        for (actual, uncut) in actual.iter().zip(expected) {
            let discontinuity = p.mul(&factor, &ComplexFloat::new(uncut.im, p.real(0)));
            relative_close(p, actual, &discontinuity, 20);
        }
        let refined = FlowOptions {
            guard_digits: 60,
            series_order: 112,
            ..options.clone()
        };
        let check = cut.evaluate(&epsilon, &refined, &context).unwrap();
        for (actual, check) in actual.iter().zip(check) {
            relative_close(p, actual, &check, 40);
        }
    }
}

#[test]
fn energy_orientation_threshold_support_and_routing_jacobian_are_explicit() {
    let options = FlowOptions::default();
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let epsilon = Rational::from((1, 13));
    let p = Precision::decimal(60).unwrap();
    for (s, reversed) in [
        (Rational::from(25), true),
        (Rational::from((1, 2)), false),
        (Rational::from(4), false),
    ] {
        let prepared = PreparedTwoBodyPhaseSpace::new(
            &family(s, reversed, 1),
            &channel(),
            &targets(),
            &KinematicPoint::default(),
            &backend,
            &options,
            &context,
        )
        .unwrap();
        assert!(
            prepared
                .evaluate(&epsilon, &options, &context)
                .unwrap()
                .iter()
                .all(|v| *v == p.zero())
        );
    }
    assert!(matches!(
        PreparedTwoBodyPhaseSpace::new(
            &family(9.into(), false, 1),
            &channel(),
            &targets(),
            &KinematicPoint::default(),
            &backend,
            &options,
            &context
        ),
        Err(Error::Unsupported(_))
    ));
    let regular = PreparedTwoBodyPhaseSpace::new(
        &family(25.into(), false, 1),
        &channel(),
        &targets(),
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    let scaled = PreparedTwoBodyPhaseSpace::new(
        &family(25.into(), false, 2),
        &channel(),
        &targets(),
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )
    .unwrap();
    let ordinary = regular.evaluate(&epsilon, &options, &context).unwrap();
    let scaled = scaled.evaluate(&epsilon, &options, &context).unwrap();
    let dimension = Rational::from(4) - Rational::from(2) * &epsilon;
    let jacobian = p.pow(&p.i(2), &p.rational(&(-dimension)));
    for (ordinary, scaled) in ordinary.iter().zip(scaled) {
        relative_close(p, &scaled, &p.mul(ordinary, &jacobian), 40);
    }
}
