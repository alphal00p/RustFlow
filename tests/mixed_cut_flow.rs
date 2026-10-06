use symbolica::prelude::*;
use symbolica_amflow::{cuts::*, *};

fn options() -> FlowOptions {
    FlowOptions {
        digits: 20,
        guard_digits: 50,
        series_order: 80,
        ..Default::default()
    }
}

fn channel() -> FutureTimelikeChannel {
    FutureTimelikeChannel {
        external: vec![1.into()],
    }
}

fn mixed_family() -> CutFamily {
    let gram = vec![vec![Atom::num(25)]];
    let mut propagators = vec![
        Propagator::quadratic(&[1, 0], &[0], Atom::num(1), &gram).unwrap(),
        Propagator::quadratic(&[-1, 0], &[1], Atom::num(4), &gram).unwrap(),
        Propagator::quadratic(&[0, 1], &[0], Atom::zero(), &gram).unwrap(),
        Propagator::quadratic(&[-1, 1], &[0], Atom::zero(), &gram).unwrap(),
    ];
    propagators.push(Propagator {
        constant: Atom::zero(),
        scalar_products: vec![
            Atom::zero(),
            Atom::zero(),
            Atom::zero(),
            Atom::zero(),
            Atom::one(),
        ],
    });
    let ordinary = IntegralFamily {
        name: "mixed_cut_massless_virtual_bubble".into(),
        loops: vec!["r".into(), "k".into()],
        external: vec!["P".into()],
        external_gram: gram,
        propagators,
        physical_propagators: 4,
        epsilon: symbol!("mixed_cut_test::eps"),
        dimension: 4,
    };
    CutFamily::new(
        ordinary,
        CutMetadata::new(
            [
                CutLine {
                    propagator: 0,
                    definition: CutDefinition::PositiveEnergy {
                        momentum: MomentumRouting {
                            loops: vec![1.into(), 0.into()],
                            external: vec![0.into()],
                        },
                    },
                },
                CutLine {
                    propagator: 1,
                    definition: CutDefinition::PositiveEnergy {
                        momentum: MomentumRouting {
                            loops: vec![(-1).into(), 0.into()],
                            external: vec![1.into()],
                        },
                    },
                },
            ],
            vec![LoopPrescription::Insensitive, LoopPrescription::PlusI0],
        )
        .unwrap(),
    )
    .unwrap()
}

fn mixed_analytic(epsilon: &Rational, p: Precision) -> ComplexFloat {
    let pi = ComplexFloat::new(p.real(1).pi(), p.real(0));
    let gamma = |r: Rational| p.gamma_real(&p.rational(&r).re).unwrap();
    let volume = p.mul(
        &p.div(
            &p.pow(&p.i(6), &p.rational(&Rational::from((1, 2)))),
            &p.scale(&pi, 25, 1),
        ),
        &p.mul(
            &p.pow(&p.scale(&pi, 25, 24), &p.rational(epsilon)),
            &p.div(
                &gamma(Rational::from((3, 2))),
                &gamma(Rational::from((3, 2)) - epsilon),
            ),
        ),
    );
    let bubble = p.div(
        &p.mul(
            &gamma(epsilon.clone()),
            &p.powi(&gamma(Rational::one() - epsilon), 2),
        ),
        &gamma(Rational::from(2) - Rational::from(2) * epsilon),
    );
    let phase = p.exp(&ComplexFloat::new(
        p.real(0),
        p.mul(&pi, &p.rational(epsilon)).re,
    ));
    p.mul(&volume, &p.mul(&bubble, &phase))
}

fn close(p: Precision, actual: &ComplexFloat, expected: &ComplexFloat, digits: u32) {
    let norm = p.norm(expected);
    let scale = if norm == p.real(0) { p.real(1) } else { norm };
    assert!(
        p.norm(&p.sub(actual, expected)) < p.tolerance(digits) * scale,
        "{actual} != {expected}"
    );
}

#[test]
fn connected_real_virtual_flow_matches_gamma_product() -> Result<()> {
    let family = mixed_family();
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let options = options();
    let targets = vec![
        Integral(vec![1, 1, 1, 1, 0]),
        Integral(vec![2, 1, 1, 1, 0]),
        Integral(vec![1, 2, 1, 1, 0]),
        Integral(vec![1, 1, 1, 1, -1]),
        Integral(vec![1, 1, 1, 1, -2]),
        Integral(vec![0, 1, 1, 1, 0]),
        Integral(vec![-1, 1, 1, 1, 0]),
        Integral(vec![1, 1, 0, 1, 0]),
    ];
    let flow = PreparedCutFlow::new(
        &family,
        &channel(),
        &targets,
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )?;
    let p = Precision::decimal(90)?;
    for epsilon in [Rational::from((1, 13)), Rational::from((1, 17))] {
        let actual = flow.evaluate(&epsilon, &options, &context)?;
        let dimension = Rational::from(4) - Rational::from(2) * &epsilon;
        let factors = [
            Rational::one(),
            Rational::from((-7, 96)) * (&dimension - &Rational::from(3)) - &epsilon,
            Rational::from((-11, 192)) * (&dimension - &Rational::from(3)),
            Rational::from((11, 2)),
            (Rational::from(121) * &dimension - Rational::from(25))
                / (Rational::from(4) * (&dimension - &Rational::one())),
            Rational::zero(),
            Rational::zero(),
            Rational::zero(),
        ];
        let base = mixed_analytic(&epsilon, p);
        for (actual, factor) in actual.iter().zip(factors) {
            close(p, actual, &p.mul(&base, &p.rational(&factor)), 20);
        }
        let mut refined = options.clone();
        refined.guard_digits += 20;
        refined.series_order += 32;
        for (actual, repeated) in actual
            .iter()
            .zip(flow.evaluate(&epsilon, &refined, &context)?)
        {
            close(p, actual, &repeated, 20);
        }
    }
    Ok(())
}

fn three_body(mass_squared: i64) -> CutFamily {
    let gram = vec![vec![Atom::one()]];
    let mut propagators = vec![
        Propagator::quadratic(&[1, 0], &[0], Atom::zero(), &gram).unwrap(),
        Propagator::quadratic(&[0, 1], &[0], Atom::zero(), &gram).unwrap(),
        Propagator::quadratic(&[-1, -1], &[1], Atom::zero(), &gram).unwrap(),
        Propagator::quadratic(&[1, 1], &[0], Atom::num(mass_squared), &gram).unwrap(),
    ];
    propagators.push(Propagator {
        constant: Atom::zero(),
        scalar_products: vec![
            Atom::zero(),
            Atom::zero(),
            Atom::zero(),
            Atom::one(),
            Atom::zero(),
        ],
    });
    let family = IntegralFamily {
        name: "weighted_three_body".into(),
        loops: vec!["q1".into(), "q2".into()],
        external: vec!["P".into()],
        external_gram: gram,
        propagators,
        physical_propagators: 4,
        epsilon: symbol!("weighted_phase_test::eps"),
        dimension: 4,
    };
    let cuts = [[1, 0], [0, 1], [-1, -1]]
        .into_iter()
        .enumerate()
        .map(|(propagator, routing)| CutLine {
            propagator,
            definition: CutDefinition::PositiveEnergy {
                momentum: MomentumRouting {
                    loops: routing.into_iter().map(Rational::from).collect(),
                    external: vec![Rational::from(i64::from(propagator == 2))],
                },
            },
        });
    CutFamily::new(
        family,
        CutMetadata::new(cuts, vec![LoopPrescription::Insensitive; 2]).unwrap(),
    )
    .unwrap()
}

// Independent convergent beta-integral series, only an analytic test oracle.
fn weighted_three_body_analytic(epsilon: &Rational, power: i64, p: Precision) -> ComplexFloat {
    let pi = ComplexFloat::new(p.real(1).pi(), p.real(0));
    let gamma = |r: Rational| p.gamma_real(&p.rational(&r).re).unwrap();
    let a = Rational::one() - epsilon;
    let base = p.div(
        &p.mul(
            &p.pow(&p.scale(&pi, 4, 1), &p.rational(epsilon)),
            &gamma(a.clone()),
        ),
        &p.mul(&p.scale(&pi, 8, 1), &gamma(Rational::from(2) * &a)),
    );
    let volume = p.mul(
        &p.div(&p.powi(&base, 2), &p.scale(&pi, 2, 1)),
        &p.div(
            &p.mul(&gamma(a.clone()), &gamma(Rational::from(2) * &a)),
            &gamma(Rational::from(3) * &a),
        ),
    );
    let mut term = p.i(1);
    let mut sum = term.clone();
    for n in 0..1000 {
        term = p.mul(
            &term,
            &p.rational(
                &((Rational::from(power + n) * (&a + &Rational::from(n)))
                    / (Rational::from(2 * (n + 1)) * (Rational::from(3) * &a + Rational::from(n)))),
            ),
        );
        sum = p.add(&sum, &term);
        if p.norm(&term) < p.tolerance(80) {
            return p.mul(
                &p.mul(&volume, &sum),
                &p.powi(&p.rational(&Rational::from((-1, 2))), power),
            );
        }
    }
    panic!("beta-integral oracle did not converge")
}

#[test]
fn weighted_three_body_flow_matches_beta_integral() -> Result<()> {
    let family = three_body(2);
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let options = options();
    let flow = PreparedCutFlow::new(
        &family,
        &channel(),
        &[Integral(vec![1, 1, 1, 1, 0]), Integral(vec![1, 1, 1, 2, 0])],
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )?;
    let p = Precision::decimal(90)?;
    let epsilon = Rational::from((1, 13));
    for (index, actual) in flow
        .evaluate(&epsilon, &options, &context)?
        .iter()
        .enumerate()
    {
        close(
            p,
            actual,
            &weighted_three_body_analytic(&epsilon, index as i64 + 1, p),
            20,
        );
    }
    let result = flow.solve(0, &options, &context)?;
    let pi = ComplexFloat::new(p.real(1).pi(), p.real(0));
    let finite = p.div(
        &p.sub(&p.log(&p.i(2)), &p.i(1)),
        &p.scale(&p.powi(&pi, 3), 128, 1),
    );
    assert_eq!(result[0].verified_digits, Some(20));
    close(p, &result[0].coefficients[&0], &finite, 20);
    Ok(())
}

#[test]
fn phase_space_pole_and_partial_cut_deformation_are_rejected() {
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let target = Integral(vec![1, 1, 1, 1, 0]);
    assert!(
        matches!(PreparedCutFlow::new(&three_body(1), &channel(), std::slice::from_ref(&target),
        &KinematicPoint::default(), &backend, &options(), &context), Err(Error::Unsupported(message)) if message.contains("away from zero"))
    );
    let mut options = options();
    options.mass_mode = MassMode::Propagators(vec![0]);
    assert!(
        matches!(PreparedCutFlow::new(&mixed_family(), &channel(), &[target],
        &KinematicPoint::default(), &backend, &options, &context), Err(Error::Unsupported(message)) if message.contains("every uncut"))
    );
}

fn routed_family(reversed: bool) -> CutFamily {
    let family = mixed_family();
    let transformation = vec![
        vec![Atom::num(2), Atom::zero()],
        vec![Atom::one(), Atom::num(3)],
    ];
    let mut ordinary = family.family().transform_loops(&transformation).unwrap();
    // Include a denominator permutation which interleaves cuts and virtual lines.
    let permutation = [2, 1, 3, 0, 4];
    ordinary.propagators = permutation
        .iter()
        .map(|&slot| ordinary.propagators[slot].clone())
        .collect();
    let lines = family.cuts().lines().iter().map(|(&old, definition)| {
        let CutDefinition::PositiveEnergy { momentum } = definition else {
            unreachable!()
        };
        let sign = if reversed { -1 } else { 1 };
        let loops = vec![
            Rational::from(sign) * (Rational::from(2) * &momentum.loops[0] + &momentum.loops[1]),
            Rational::from(3 * sign) * &momentum.loops[1],
        ];
        CutLine {
            propagator: permutation.iter().position(|&slot| slot == old).unwrap(),
            definition: CutDefinition::PositiveEnergy {
                momentum: MomentumRouting {
                    loops,
                    external: momentum
                        .external
                        .iter()
                        .map(|a| Rational::from(sign) * a)
                        .collect(),
                },
            },
        }
    });
    CutFamily::new(
        ordinary,
        CutMetadata::new(lines, family.cuts().loop_prescriptions().to_vec()).unwrap(),
    )
    .unwrap()
}

#[test]
fn triangular_routing_permutation_and_reversed_support_preserve_measure() -> Result<()> {
    let backend = RustRedBackend::default();
    let options = options();
    let context = RunContext::default();
    let epsilon = Rational::from((1, 17));
    let p = Precision::decimal(90)?;
    let targets = [Integral(vec![1, 1, 1, 1, 0]), Integral(vec![1, 1, 1, 2, 0])];
    for reversed in [false, true] {
        let family = routed_family(reversed);
        let flow = PreparedCutFlow::new(
            &family,
            &channel(),
            &targets,
            &KinematicPoint::default(),
            &backend,
            &options,
            &context,
        )?;
        let values = flow.evaluate(&epsilon, &options, &context)?;
        if reversed {
            assert!(values.iter().all(|value| *value == p.zero()));
            assert!(flow.differential_system().is_none());
        } else {
            let jacobian = p.pow(
                &p.i(6),
                &p.rational(&(Rational::from(4) - Rational::from(2) * &epsilon)),
            );
            let base = mixed_analytic(&epsilon, p);
            close(p, &p.mul(&values[0], &jacobian), &base, 20);
            let raised = Rational::from((-7, 96))
                * (Rational::one() - Rational::from(2) * &epsilon)
                - &epsilon;
            close(
                p,
                &p.mul(&values[1], &jacobian),
                &p.mul(&base, &p.rational(&raised)),
                20,
            );
        }
    }
    Ok(())
}

struct CheckedBackend {
    calls: std::sync::atomic::AtomicUsize,
    guard: bool,
}
impl ReductionBackend for CheckedBackend {
    fn identity(&self) -> String {
        format!(
            "mixed-cut-check:{}:{}",
            self.guard,
            RustRedBackend::default().identity()
        )
    }
    fn reduce(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        context: &RunContext,
    ) -> Result<reduction::Reduction> {
        self.calls.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        RustRedBackend::default().reduce(family, targets, context)
    }
    fn reduce_cut(
        &self,
        family: &CutFamily,
        targets: &[Integral],
        context: &RunContext,
    ) -> Result<reduction::Reduction> {
        self.calls.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        let mut result = RustRedBackend::default().reduce_cut(family, targets, context)?;
        if self.guard {
            result
                .nonzero_conditions
                .push(Atom::var(family.family().epsilon) - Atom::num((1, 13)));
        }
        Ok(result)
    }
}

#[test]
fn cut_system_restart_nonzero_conditions_and_cancellation_remain_explicit() -> Result<()> {
    use std::sync::atomic::{AtomicUsize, Ordering};
    let family = mixed_family();
    let target = Integral(vec![1, 1, 1, 1, 0]);
    let backend = CheckedBackend {
        calls: AtomicUsize::new(0),
        guard: false,
    };
    let context = RunContext::default();
    let directory =
        std::env::temp_dir().join(format!("mixed-cut-flow-cache-{}", std::process::id()));
    if directory.exists() {
        std::fs::remove_dir_all(&directory)?;
    }
    let mut options = options();
    options.cache_directory = Some(directory.clone());
    let first = PreparedCutFlow::new(
        &family,
        &channel(),
        std::slice::from_ref(&target),
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )?;
    assert!(backend.calls.load(Ordering::SeqCst) > 0);
    backend.calls.store(0, Ordering::SeqCst);
    let restarted = PreparedCutFlow::new(
        &family,
        &channel(),
        std::slice::from_ref(&target),
        &KinematicPoint::default(),
        &backend,
        &options,
        &context,
    )?;
    assert_eq!(backend.calls.load(Ordering::SeqCst), 0);
    assert_eq!(
        first.differential_system().unwrap().matrix,
        restarted.differential_system().unwrap().matrix
    );
    let epsilon = Rational::from((1, 17));
    close(
        Precision::decimal(90)?,
        &restarted.evaluate(&epsilon, &options, &context)?[0],
        &mixed_analytic(&epsilon, Precision::decimal(90)?),
        20,
    );
    std::fs::remove_dir_all(directory)?;
    options.cache_directory = None;
    let guarded = CheckedBackend {
        calls: AtomicUsize::new(0),
        guard: true,
    };
    let flow = PreparedCutFlow::new(
        &family,
        &channel(),
        std::slice::from_ref(&target),
        &KinematicPoint::default(),
        &guarded,
        &options,
        &context,
    )?;
    assert!(
        matches!(flow.evaluate(&Rational::from((1, 13)), &options, &context),
        Err(Error::Reduction(message)) if message.contains("vanishes"))
    );
    let cancelled = RunContext::default();
    cancelled.cancellation.cancel();
    assert!(matches!(
        flow.evaluate(&epsilon, &options, &cancelled),
        Err(Error::Cancelled)
    ));
    let (ordinary, cuts) = family.into_parts();
    let opposite = CutFamily::new(
        ordinary,
        CutMetadata::new(
            cuts.lines()
                .iter()
                .map(|(&propagator, definition)| CutLine {
                    propagator,
                    definition: definition.clone(),
                }),
            vec![LoopPrescription::Insensitive, LoopPrescription::MinusI0],
        )?,
    )?;
    assert!(
        matches!(PreparedCutFlow::new(&opposite, &channel(), &[target],
        &KinematicPoint::default(), &backend, &options, &context), Err(Error::Unsupported(message)) if message.contains("uniform +i0"))
    );
    Ok(())
}
