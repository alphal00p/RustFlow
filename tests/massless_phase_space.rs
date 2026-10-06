use feynkit_kinematics::Kinematics;
use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::{cuts::*, phase_space::*, reduction::Reduction, *};

fn family(loops: usize, scaled: bool, reversed: bool, permuted: bool) -> CutFamily {
    let gram = vec![vec![Atom::num(7)]];
    let slots = loops * (loops + 1) / 2 + loops;
    let mut routings = Vec::new();
    for i in 0..loops {
        let mut q = vec![0_i64; loops];
        q[i] = if scaled { (i + 2) as i64 } else { 1 };
        // Affine translations do not change the phase-space measure.
        routings.push((q, if scaled && i == 0 { 2 } else { 0 }));
    }
    let total = (0..loops)
        .map(|j| -routings.iter().map(|q| q.0[j]).sum::<i64>())
        .collect();
    let external = 1 - routings.iter().map(|q| q.1).sum::<i64>();
    routings.push((total, external));
    if permuted {
        routings.rotate_left(1);
    }
    let mut propagators = routings
        .iter()
        .map(|(q, e)| Propagator::quadratic(q, &[*e], Atom::new(), &gram).unwrap())
        .collect::<Vec<_>>();
    let mut index = 0;
    for i in 0..loops {
        for j in i..loops {
            if i != j {
                let mut c = vec![Atom::new(); slots];
                c[index] = Atom::one();
                propagators.push(Propagator {
                    constant: Atom::new(),
                    scalar_products: c,
                });
            }
            index += 1;
        }
    }
    for i in 0..loops - 1 {
        let mut c = vec![Atom::new(); slots];
        c[index + i] = Atom::one();
        propagators.push(Propagator {
            constant: Atom::new(),
            scalar_products: c,
        });
    }
    let ordinary = IntegralFamily {
        name: format!("massless_{}_body", loops + 1),
        loops: (0..loops).map(|i| format!("k{i}")).collect(),
        external: vec!["p".into()],
        external_gram: gram,
        propagators,
        physical_propagators: loops + 1,
        epsilon: symbol!("massless_phase_test::eps"),
        dimension: 4,
    };
    let sign = if reversed { -1 } else { 1 };
    let cuts = CutMetadata::new(
        routings
            .into_iter()
            .enumerate()
            .map(|(propagator, (q, external))| CutLine {
                propagator,
                definition: CutDefinition::PositiveEnergy {
                    momentum: MomentumRouting {
                        loops: q.into_iter().map(|x| Rational::from(sign * x)).collect(),
                        external: vec![Rational::from(sign * external)],
                    },
                },
            }),
        vec![LoopPrescription::Insensitive; loops],
    )
    .unwrap();
    CutFamily::new(ordinary, cuts).unwrap()
}
fn channel() -> FutureTimelikeChannel {
    FutureTimelikeChannel {
        external: vec![Rational::one()],
    }
}
fn unit(f: &CutFamily) -> Integral {
    Integral(
        (0..f.family().propagators.len())
            .map(|i| i16::from(f.cuts().is_cut(i)))
            .collect(),
    )
}
fn options() -> FlowOptions {
    FlowOptions {
        digits: 20,
        guard_digits: 50,
        series_order: 80,
        ..Default::default()
    }
}
fn close(p: Precision, a: &ComplexFloat, b: &ComplexFloat, digits: u32) {
    let norm = p.norm(b);
    let scale = if norm == p.real(0) { p.real(1) } else { norm };
    assert!(
        p.norm(&p.sub(a, b)) <= p.tolerance(digits) * scale,
        "{a} != {b}"
    );
}
struct NoReduction;
impl ReductionBackend for NoReduction {
    fn identity(&self) -> String {
        "no-reduction-for-proved-unit-or-zero".into()
    }
    fn reduce(&self, _: &IntegralFamily, _: &[Integral], _: &RunContext) -> Result<Reduction> {
        panic!("uncut reducer must never be used")
    }
}

#[test]
fn massless_two_and_three_loop_samples_match_automatic_uncut_discontinuities() -> Result<()> {
    let options = options();
    let context = RunContext::default();
    let backend = RustRedBackend::default();
    let p = Precision::decimal(80)?;
    for loops in [2, 3] {
        let f = family(loops, false, false, false);
        let u = unit(&f);
        let phase = PreparedMasslessPhaseSpace::new(
            &f,
            &channel(),
            std::slice::from_ref(&u),
            &KinematicPoint::default(),
            &backend,
            &options,
            &context,
        )?;
        assert_eq!(phase.leading_power(), -2 * loops as i32);
        let uncut = PreparedFlow::new(
            f.family(),
            &[u],
            &KinematicPoint::default(),
            &backend,
            &options,
            &context,
        )?;
        for eps in [Rational::from((1, 13)), Rational::from((1, 17))] {
            let actual = phase.evaluate(&eps, &options, &context)?;
            let uncut = uncut.evaluate(
                &eps,
                &options,
                &recursive::RecursiveBoundary::new(&backend, &options, &context),
                &context,
            )?;
            let d = Rational::from(4) - Rational::from(2) * &eps;
            let pi = ComplexFloat::new(p.real(1).pi(), p.real(0));
            let factor = p.scale(
                &p.pow(
                    &p.scale(&pi, 4, 1),
                    &p.rational(&(-Rational::from(loops as i64) * d / Rational::from(2))),
                ),
                if (loops + 1) % 2 == 0 { 2 } else { -2 },
                1,
            );
            close(
                p,
                &actual[0],
                &p.mul(&factor, &ComplexFloat::new(uncut[0].im.clone(), p.real(0))),
                20,
            );
            let refined = phase.evaluate(
                &eps,
                &FlowOptions {
                    guard_digits: 70,
                    series_order: 112,
                    ..options.clone()
                },
                &context,
            )?;
            close(p, &actual[0], &refined[0], 50);
        }
    }
    Ok(())
}

#[test]
fn exact_routing_permutation_energy_and_pinched_cuts_use_proved_geometry() -> Result<()> {
    let options = options();
    let context = RunContext::default();
    let p = Precision::decimal(80)?;
    let eps = Rational::from((1, 13));
    for loops in [1, 2, 3, 4] {
        let regular = family(loops, false, false, false);
        let u = unit(&regular);
        let baseline = PreparedMasslessPhaseSpace::new(
            &regular,
            &channel(),
            std::slice::from_ref(&u),
            &KinematicPoint::default(),
            &NoReduction,
            &options,
            &context,
        )?
        .evaluate(&eps, &options, &context)?;
        let scaled = family(loops, true, false, true);
        let value = PreparedMasslessPhaseSpace::new(
            &scaled,
            &channel(),
            std::slice::from_ref(&u),
            &KinematicPoint::default(),
            &NoReduction,
            &options,
            &context,
        )?
        .evaluate(&eps, &options, &context)?;
        let determinant = (2..=loops + 1).product::<usize>() as i64;
        let jacobian = p.pow(
            &p.i(determinant),
            &p.rational(&(-Rational::from(4) + Rational::from(2) * &eps)),
        );
        close(p, &value[0], &p.mul(&baseline[0], &jacobian), 50);
        let reversed = family(loops, true, true, true);
        let zeros = PreparedMasslessPhaseSpace::new(
            &reversed,
            &channel(),
            std::slice::from_ref(&u),
            &KinematicPoint::default(),
            &NoReduction,
            &options,
            &context,
        )?
        .evaluate(&eps, &options, &context)?;
        assert_eq!(zeros[0], p.zero());
        let mut pinched = u.clone();
        pinched.0[0] = 0;
        let zeros = PreparedMasslessPhaseSpace::new(
            &regular,
            &channel(),
            &[pinched],
            &KinematicPoint::default(),
            &NoReduction,
            &options,
            &context,
        )?
        .evaluate(&eps, &options, &context)?;
        assert_eq!(zeros[0], p.zero());
    }
    Ok(())
}

#[test]
fn native_raised_cuts_and_numerators_close_and_respect_missing_cut_zero() -> Result<()> {
    let f = family(2, false, false, false);
    let u = unit(&f);
    let mut raised = u.clone();
    raised.0[0] = 2;
    let mut numerator = u.clone();
    numerator.0[3] = -1;
    let mut pinched = raised.clone();
    pinched.0[1] = -1;
    let targets = vec![u, raised, numerator, pinched];
    let options = options();
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let p = Precision::decimal(80)?;
    let phase = PreparedMasslessPhaseSpace::new(
        &f,
        &channel(),
        &targets,
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )?;
    let uncut = PreparedFlow::new(
        f.family(),
        &targets[..3],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )?;
    let eps = Rational::from((1, 13));
    let values = phase.evaluate(&eps, &options, &context)?;
    assert_eq!(values[3], p.zero());
    let actual = uncut.evaluate(
        &eps,
        &options,
        &recursive::RecursiveBoundary::new(&backend, &options, &context),
        &context,
    )?;
    let pi = ComplexFloat::new(p.real(1).pi(), p.real(0));
    let d = Rational::from(4) - Rational::from(2) * &eps;
    let factor = p.scale(&p.pow(&p.scale(&pi, 4, 1), &p.rational(&(-d))), -2, 1);
    for (v, i) in values[..3].iter().zip(actual) {
        close(
            p,
            v,
            &p.mul(&factor, &ComplexFloat::new(i.im, p.real(0))),
            20,
        );
    }
    Ok(())
}

#[test]
fn public_massless_laurent_fit_matches_native_three_body_dalitz_measure() -> Result<()> {
    let f = family(2, false, false, false);
    let context = RunContext::default();
    let p = Precision::decimal(80)?;
    let answer = solve_massless_phase_space(
        &f,
        &channel(),
        &[unit(&f)],
        &KinematicPoint::default(),
        0,
        &options(),
        &NoReduction,
        &context,
    )?;
    assert_eq!(answer[0].verified_digits, Some(20));
    assert!(answer[0].validation_samples > 0 && answer[0].refinements > 0);
    let [a, b, c] = [
        parse!("massless_phase_test::p1"),
        parse!("massless_phase_test::p2"),
        parse!("massless_phase_test::p3"),
    ];
    let mut native = Kinematics::new();
    for q in [&a, &b, &c] {
        native = native.with_mass_squared(q, Atom::new()).unwrap();
    }
    native = native
        .with_scalar_product(&a, &b, Atom::one())
        .unwrap()
        .with_scalar_product(&b, &c, Atom::one())
        .unwrap()
        .with_scalar_product(&a, &c, Atom::num((3, 2)))
        .unwrap();
    let volume = native.three_body_phase_space(&a, &b, &c).unwrap() * Atom::num((49, 2));
    close(
        p,
        &answer[0].coefficients[&0],
        &p.eval(&volume, &Default::default())?,
        20,
    );
    for power in -4..0 {
        assert!(p.norm(&answer[0].coefficients[&power]) < p.tolerance(20));
    }
    Ok(())
}

#[derive(Clone, Copy)]
enum BadReduction {
    Omitted,
    Residual,
    Pole,
    Guard,
}
impl ReductionBackend for BadReduction {
    fn identity(&self) -> String {
        "cut-table-domain-regression".into()
    }
    fn reduce(&self, _: &IntegralFamily, _: &[Integral], _: &RunContext) -> Result<Reduction> {
        panic!("uncut backend")
    }
    fn reduce_cut(&self, f: &CutFamily, targets: &[Integral], _: &RunContext) -> Result<Reduction> {
        let u = unit(f);
        let eps = Atom::var(f.family().epsilon);
        match self {
            Self::Omitted => Ok(Reduction::default()),
            Self::Residual => Ok(Reduction {
                residuals: targets.to_vec(),
                ..Default::default()
            }),
            Self::Pole | Self::Guard => Ok(Reduction {
                rules: targets
                    .iter()
                    .cloned()
                    .map(|t| {
                        (
                            t,
                            BTreeMap::from([(
                                u.clone(),
                                if matches!(self, Self::Pole) {
                                    eps.clone().pow(-3)
                                } else {
                                    Atom::one()
                                },
                            )]),
                        )
                    })
                    .collect(),
                residuals: vec![u],
                nonzero_conditions: if matches!(self, Self::Guard) {
                    vec![eps - Atom::num((1, 13))]
                } else {
                    vec![]
                },
            }),
        }
    }
}
#[test]
fn reduction_coverage_conditions_and_additional_epsilon_poles_are_not_dropped() -> Result<()> {
    let f = family(2, false, false, false);
    let u = unit(&f);
    let mut raised = u.clone();
    raised.0[0] = 2;
    let options = options();
    let context = RunContext::default();
    for backend in [BadReduction::Omitted, BadReduction::Residual] {
        assert!(matches!(
            PreparedMasslessPhaseSpace::new(
                &f,
                &channel(),
                std::slice::from_ref(&raised),
                &KinematicPoint::default(),
                &backend,
                &options,
                &context
            ),
            Err(Error::IncompleteReduction(_))
        ));
    }
    let guarded = PreparedMasslessPhaseSpace::new(
        &f,
        &channel(),
        std::slice::from_ref(&raised),
        &KinematicPoint::default(),
        &BadReduction::Guard,
        &options,
        &context,
    )?;
    assert!(matches!(
        guarded.evaluate(&Rational::from((1, 13)), &options, &context),
        Err(Error::Numerical(_))
    ));
    let pole = PreparedMasslessPhaseSpace::new(
        &f,
        &channel(),
        &[raised],
        &KinematicPoint::default(),
        &BadReduction::Pole,
        &options,
        &context,
    )?;
    assert_eq!(pole.leading_power(), -7);
    let p = Precision::decimal(80)?;
    let eps = Rational::from((1, 17));
    let base = PreparedMasslessPhaseSpace::new(
        &f,
        &channel(),
        &[u],
        &KinematicPoint::default(),
        &NoReduction,
        &options,
        &context,
    )?
    .evaluate(&eps, &options, &context)?;
    close(
        p,
        &pole.evaluate(&eps, &options, &context)?[0],
        &p.scale(&base[0], 17 * 17 * 17, 1),
        50,
    );
    Ok(())
}

#[test]
fn mixed_massive_or_extra_denominator_inputs_fail_explicitly() -> Result<()> {
    let f = family(2, false, false, false);
    let options = options();
    let context = RunContext::default();
    let u = unit(&f);
    let mut extra = u.clone();
    extra.0[3] = 1;
    assert!(matches!(
        PreparedMasslessPhaseSpace::new(
            &f,
            &channel(),
            &[extra],
            &KinematicPoint::default(),
            &NoReduction,
            &options,
            &context
        ),
        Err(Error::InvalidInput(_))
    ));
    let mut ordinary = f.family().clone();
    ordinary.propagators[0].constant -= Atom::one();
    let massive = CutFamily::new(ordinary, f.cuts().clone())?;
    assert!(matches!(
        PreparedMasslessPhaseSpace::new(
            &massive,
            &channel(),
            std::slice::from_ref(&u),
            &KinematicPoint::default(),
            &NoReduction,
            &options,
            &context
        ),
        Err(Error::Unsupported(_))
    ));
    let (_, cuts) = f.clone().into_parts();
    let lines = cuts.lines().iter().skip(1).map(|(&slot, def)| CutLine {
        propagator: slot,
        definition: def.clone(),
    });
    let mixed = CutFamily::new(
        f.family().clone(),
        CutMetadata::new(lines, vec![LoopPrescription::PlusI0; 2])?,
    )?;
    assert!(matches!(
        PreparedMasslessPhaseSpace::new(
            &mixed,
            &channel(),
            std::slice::from_ref(&u),
            &KinematicPoint::default(),
            &NoReduction,
            &options,
            &context
        ),
        Err(Error::Unsupported(_))
    ));
    let wrong = FutureTimelikeChannel {
        external: vec![Rational::from(2)],
    };
    assert!(matches!(
        PreparedMasslessPhaseSpace::new(
            &f,
            &wrong,
            &[u],
            &KinematicPoint::default(),
            &NoReduction,
            &options,
            &context
        ),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}

#[test]
fn several_external_coordinates_and_nondefault_dimension_match_the_same_channel() -> Result<()> {
    let gram = vec![
        vec![Atom::new(), Atom::num((7, 2))],
        vec![Atom::num((7, 2)), Atom::new()],
    ];
    let epsilon = symbol!("massless_phase_test::eps");
    let ordinary = IntegralFamily {
        name: "two_external_massless_channel".into(),
        loops: vec!["k".into()],
        external: vec!["p1".into(), "p2".into()],
        external_gram: gram.clone(),
        propagators: vec![
            Propagator::quadratic(&[1], &[1, 0], Atom::new(), &gram)?,
            Propagator::quadratic(&[-1], &[0, 1], Atom::new(), &gram)?,
            Propagator {
                constant: Atom::new(),
                scalar_products: vec![Atom::one(), Atom::new(), Atom::new()],
            },
        ],
        physical_propagators: 2,
        epsilon,
        dimension: 4,
    };
    let cuts = CutMetadata::new(
        [
            CutLine {
                propagator: 0,
                definition: CutDefinition::PositiveEnergy {
                    momentum: MomentumRouting {
                        loops: vec![1.into()],
                        external: vec![1.into(), 0.into()],
                    },
                },
            },
            CutLine {
                propagator: 1,
                definition: CutDefinition::PositiveEnergy {
                    momentum: MomentumRouting {
                        loops: vec![(-1).into()],
                        external: vec![0.into(), 1.into()],
                    },
                },
            },
        ],
        vec![LoopPrescription::Insensitive],
    )?;
    let f = CutFamily::new(ordinary, cuts)?;
    let channel = FutureTimelikeChannel {
        external: vec![1.into(), 1.into()],
    };
    let options = FlowOptions {
        dimension: 6,
        ..options()
    };
    let context = RunContext::default();
    let eps = Rational::from((1, 13));
    let p = Precision::decimal(80)?;
    let volume = PreparedMasslessPhaseSpace::new(
        &f,
        &channel,
        &[unit(&f)],
        &KinematicPoint::default(),
        &NoReduction,
        &options,
        &context,
    )?
    .evaluate(&eps, &options, &context)?;
    let reference = family(1, false, false, false);
    let backend = RustRedBackend::default();
    let flow = PreparedFlow::new(
        reference.family(),
        &[unit(&reference)],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )?;
    let uncut = flow.evaluate(&eps, &options, &boundary::OneLoopBoundary, &context)?;
    let pi = ComplexFloat::new(p.real(1).pi(), p.real(0));
    let dimension = Rational::from(6) - Rational::from(2) * eps;
    let factor = p.scale(
        &p.pow(
            &p.scale(&pi, 4, 1),
            &p.rational(&(-dimension / Rational::from(2))),
        ),
        2,
        1,
    );
    close(
        p,
        &volume[0],
        &p.mul(&factor, &ComplexFloat::new(uncut[0].im.clone(), p.real(0))),
        20,
    );
    Ok(())
}

#[test]
fn two_particle_gamma_cancellation_retains_the_finite_dimension_two_value() -> Result<()> {
    let f = family(1, false, false, false);
    let context = RunContext::default();
    let p = Precision::decimal(80)?;
    let prepared = PreparedMasslessPhaseSpace::new(
        &f,
        &channel(),
        &[unit(&f)],
        &KinematicPoint::default(),
        &NoReduction,
        &options(),
        &context,
    )?;
    // D=4-2epsilon=2: the native two-body angular measure gives Phi2(s)=1/s.
    // A separately evaluated Gamma(a)/Gamma(a) would spuriously fail at a=0.
    let result = prepared.evaluate(&Rational::one(), &options(), &context)?;
    close(p, &result[0], &p.rational(&Rational::from((1, 7))), 50);
    Ok(())
}

#[test]
fn generic_cut_dispatch_preserves_nbody_routing_and_weight_poles() -> Result<()> {
    let context = RunContext::default();
    let options = options();
    let p = Precision::decimal(70)?;
    let epsilon = Rational::from((1, 13));
    for loops in [2, 3] {
        let family = family(loops, true, false, true);
        let target = unit(&family);
        let row = BTreeMap::from([(
            target.clone(),
            Atom::one() / Atom::var(family.family().epsilon).pow(3),
        )]);
        let projections = PreparedCutProjections::new(
            &family,
            &channel(),
            &[row],
            &KinematicPoint::default(),
            &NoReduction,
            &options,
            &context,
        )?;
        let terminal = PreparedMasslessPhaseSpace::new(
            &family,
            &channel(),
            &[target],
            &KinematicPoint::default(),
            &NoReduction,
            &options,
            &context,
        )?;
        assert_eq!(projections.leading_power(), terminal.leading_power() - 3);
        let expected = p.div(
            &terminal.evaluate(&epsilon, &options, &context)?[0],
            &p.powi(&p.rational(&epsilon), 3),
        );
        close(
            p,
            &projections.evaluate(&epsilon, &options, &context)?[0],
            &expected,
            50,
        );
    }
    Ok(())
}
