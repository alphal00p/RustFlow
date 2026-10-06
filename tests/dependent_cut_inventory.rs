use feynkit_graph::IntegralFamily as NativeFamily;
use feynkit_kinematics::Kinematics;
use symbolica::prelude::*;
use symbolica_amflow::{cuts::*, *};

fn eps() -> Symbol {
    symbol!("dependent_cut_test::eps")
}
fn channel() -> FutureTimelikeChannel {
    FutureTimelikeChannel {
        external: vec![1.into()],
    }
}
fn options() -> FlowOptions {
    FlowOptions {
        digits: 20,
        guard_digits: 50,
        series_order: 80,
        mass_mode: MassMode::All,
        ..Default::default()
    }
}
fn native(three_body: bool, mass: i64, reverse: bool) -> (NativeFamily, CutMetadata) {
    let dimension = Atom::var(symbol!("dependent_cut_test::D"));
    let r = Atom::var(symbol!("dependent_cut_test::r"));
    let l = Atom::var(symbol!("dependent_cut_test::l"));
    let k = Atom::var(symbol!("dependent_cut_test::k"));
    let p = Atom::var(symbol!("dependent_cut_test::p"));
    let kin = Kinematics::in_dimension(&dimension)
        .unwrap()
        .with_momenta([r.clone(), l.clone(), k.clone(), p.clone()])
        .unwrap()
        .with_scalar_product(&p, &p, Atom::one())
        .unwrap();
    let (loops, mut routing) = if three_body {
        (
            vec![r.clone(), l.clone(), k.clone()],
            vec![(vec![1, 0, 0], 0), (vec![0, 1, 0], 0), (vec![-1, -1, 0], 1)],
        )
    } else {
        (vec![r.clone()], vec![(vec![1], 0), (vec![-1], 1)])
    };
    let cut_count = routing.len();
    let masses = if three_body {
        routing.extend([(vec![0, 0, 1], 0), (vec![0, 0, 1], 0), (vec![-1, -1, 1], 0)]);
        vec![0, 0, 0, 2, mass, 5]
    } else {
        routing.push((vec![1], 0));
        vec![0, 0, mass]
    };
    let denominators = routing
        .iter()
        .zip(masses)
        .map(|((q, e), m)| {
            let momentum = loops
                .iter()
                .zip(q)
                .fold(Atom::num(*e) * &p, |sum, (q, c)| sum + Atom::num(*c) * q);
            kin.scalar_product(&momentum, &momentum).unwrap() - Atom::num(m)
        })
        .collect();
    let mut prescriptions = vec![LoopPrescription::Insensitive; loops.len()];
    if three_body {
        *prescriptions.last_mut().unwrap() = LoopPrescription::PlusI0;
    }
    let sign = if reverse { -1 } else { 1 };
    let cuts = CutMetadata::new(
        routing
            .into_iter()
            .take(cut_count)
            .enumerate()
            .map(|(propagator, (q, e))| CutLine {
                propagator,
                definition: CutDefinition::PositiveEnergy {
                    momentum: MomentumRouting {
                        loops: q.into_iter().map(|n| Rational::from(sign * n)).collect(),
                        external: vec![(sign * e).into()],
                    },
                },
            }),
        prescriptions,
    )
    .unwrap();
    (
        NativeFamily::new(loops, vec![p], denominators, &kin).unwrap(),
        cuts,
    )
}
fn volume(epsilon: &Rational, three_body: bool, precision: Precision) -> Result<ComplexFloat> {
    let p = precision;
    let a = Rational::one() - epsilon;
    let gamma = |x: Rational| p.gamma_real(&p.rational(&x).re);
    let pi = ComplexFloat::new(p.real(1).pi(), p.real(0));
    let phi = p.div(
        &p.mul(
            &p.pow(&p.scale(&pi, 4, 1), &p.rational(epsilon)),
            &gamma(a.clone())?,
        ),
        &p.mul(&p.scale(&pi, 8, 1), &gamma(Rational::from(2) * &a)?),
    );
    if three_body {
        Ok(p.mul(
            &p.div(&p.mul(&phi, &phi), &p.scale(&pi, 2, 1)),
            &p.div(
                &p.mul(&gamma(a.clone())?, &gamma(Rational::from(2) * &a)?),
                &gamma(Rational::from(3) * &a)?,
            ),
        ))
    } else {
        Ok(phi)
    }
}
fn close(p: Precision, actual: &ComplexFloat, expected: &ComplexFloat, digits: u32) {
    let scale = p.norm(expected);
    let scale = if scale == p.real(0) { p.real(1) } else { scale };
    assert!(
        p.norm(&p.sub(actual, expected)) < p.tolerance(digits) * scale,
        "{actual} != {expected}"
    );
}

#[test]
fn smooth_dependent_pole_raised_cuts_and_both_energy_orientations() -> Result<()> {
    let backend = RustRedBackend::default();
    let options = options();
    let context = RunContext::default();
    let p = Precision::decimal(90)?;
    for reversed in [false, true] {
        let (family, cuts) = native(false, 2, reversed);
        assert!(!family.is_independent());
        for (powers, factor) in [
            (vec![1, 1, 1], Atom::num((-1, 2))),
            (vec![2, 1, 1], Atom::num((1, 4)) - Atom::var(eps())),
            (vec![0, 1, 1], Atom::zero()),
        ] {
            let prepared = PreparedCutCombination::from_hepkit(
                &family,
                &powers,
                &Atom::one(),
                cuts.clone(),
                &channel(),
                eps(),
                &backend,
                &options,
                10000,
                &context,
            )?;
            for epsilon in [Rational::from((1, 13)), Rational::from((1, 17))] {
                let actual = prepared.evaluate(&epsilon, &options, &context)?.remove(0);
                let expected = if reversed {
                    p.zero()
                } else {
                    p.mul(
                        &volume(&epsilon, false, p)?,
                        &p.eval(
                            &factor,
                            &[(Atom::var(eps()), p.rational(&epsilon))]
                                .into_iter()
                                .collect(),
                        )?,
                    )
                };
                close(p, &actual, &expected, 35);
            }
        }
    }
    Ok(())
}

#[test]
fn original_domain_is_admitted_before_cut_zero_or_partial_fractions() -> Result<()> {
    let backend = RustRedBackend::default();
    let options = options();
    let context = RunContext::default();
    // This pole coincides with a cut and is not a well-defined smooth multiplier.
    let (family, cuts) = native(false, 0, false);
    assert!(matches!(
        PreparedCutCombination::from_hepkit(
            &family,
            &[0, 1, 1],
            &Atom::zero(),
            cuts,
            &channel(),
            eps(),
            &backend,
            &options,
            10000,
            &context
        ),
        Err(Error::Unsupported(_))
    ));
    let (family, cuts) = native(false, 2, false);
    assert!(matches!(
        PreparedCutCombination::from_hepkit(
            &family,
            &[1, 1, 1],
            &Atom::one(),
            cuts.clone(),
            &channel(),
            eps(),
            &backend,
            &options,
            0,
            &context
        ),
        Err(Error::Limit(_))
    ));
    context.cancellation.cancel();
    assert!(matches!(
        PreparedCutCombination::from_hepkit(
            &family,
            &[1, 1, 1],
            &Atom::one(),
            cuts,
            &channel(),
            eps(),
            &backend,
            &options,
            10000,
            &context
        ),
        Err(Error::Cancelled)
    ));
    Ok(())
}

#[test]
fn exact_epsilon_weights_are_combined_before_one_verified_fit() -> Result<()> {
    let backend = RustRedBackend::default();
    let options = options();
    let context = RunContext::default();
    let (family, cuts) = native(false, 2, false);
    let prepared = PreparedCutCombination::from_hepkit(
        &family,
        &[2, 1, 1],
        &(Atom::one() / Atom::var(eps())),
        cuts,
        &channel(),
        eps(),
        &backend,
        &options,
        10000,
        &context,
    )?;
    assert_eq!(prepared.leading_power(), -3);
    let result = prepared.solve(0, &options, &context)?.remove(0);
    assert_eq!(result.verified_digits, Some(20));
    let p = Precision::decimal(90)?;
    let pi = ComplexFloat::new(p.real(1).pi(), p.real(0));
    let phi0 = p.div(&p.i(1), &p.scale(&pi, 8, 1));
    let gamma = ComplexFloat::new(
        Float::from_raw(rug::Float::with_val(p.bits, rug::float::Constant::Euler)),
        p.real(0),
    );
    let phi_derivative = p.sub(&p.add(&p.i(2), &p.log(&p.scale(&pi, 4, 1))), &gamma);
    close(p, &result.coefficients[&-1], &p.scale(&phi0, 1, 4), 20);
    close(
        p,
        &result.coefficients[&0],
        &p.mul(&phi0, &p.sub(&p.scale(&phi_derivative, 1, 4), &p.i(1))),
        20,
    );
    for power in [-3, -2] {
        close(p, &result.coefficients[&power], &p.zero(), 20);
    }
    assert!(result.validation_samples > 0);
    Ok(())
}

// Independent Feynman-parameter oracle. With t=(r+l)^2 the three-body
// density is Beta(1-eps,2-2eps). Expanding the bubble about mass 5 gives
// delta=[(5-a)u+u(1-u)t]/5 <= 3/5 for a=2,3. For 0<eps<1 the
// omitted normalized tail after N terms is <= (3/5)^N/(1-3/5).
fn bubble_difference(epsilon: &Rational, precision: Precision) -> Result<ComplexFloat> {
    let a = Rational::one() - epsilon;
    let mut sum = Rational::zero();
    let mut binomial = Rational::one();
    for n in 0..220_i64 {
        let mut inner = Rational::zero();
        let mut choose = Rational::one();
        let mut beta = Rational::from((1, n + 1));
        let mut moment = Rational::one();
        let mut power_two = Rational::from(2).pow(n as u64);
        let mut power_three = Rational::from(3).pow(n as u64);
        for j in 0..=n {
            inner += &choose * &beta * &moment * (&power_two - &power_three);
            choose *= Rational::from((n - j, j + 1));
            beta *= Rational::from((j + 1, n + j + 2));
            moment *= (&a + &Rational::from(j)) / (Rational::from(3) * &a + Rational::from(j));
            power_two /= Rational::from(2);
            power_three /= Rational::from(3);
        }
        sum += &binomial * &inner;
        binomial *= (epsilon + &Rational::from(n)) / Rational::from(5 * (n + 1));
    }
    let p = precision;
    Ok(p.mul(
        &volume(epsilon, true, p)?,
        &p.mul(
            &p.gamma_real(&p.rational(epsilon).re)?,
            &p.mul(
                &p.pow(&p.i(5), &p.neg(&p.rational(epsilon))),
                &p.rational(&sum),
            ),
        ),
    ))
}

#[test]
fn nonfactorized_overcomplete_virtual_inventory_matches_independent_double_beta_series()
-> Result<()> {
    let (family, cuts) = native(true, 3, false);
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let p = Precision::decimal(90)?;
    let mut baseline = None;
    for mass_mode in [MassMode::All, MassMode::Propagators(vec![5])] {
        let mut settings = options();
        settings.mass_mode = mass_mode;
        let start = std::time::Instant::now();
        let prepared = PreparedCutCombination::from_hepkit(
            &family,
            &[1; 6],
            &Atom::one(),
            cuts.clone(),
            &channel(),
            eps(),
            &backend,
            &settings,
            10000,
            &context,
        )?;
        eprintln!(
            "dependent cut placement={:?} preparation={:?}",
            settings.mass_mode,
            start.elapsed()
        );
        assert_eq!(prepared.original_slots().len(), 2);
        assert!(
            prepared
                .original_slots()
                .iter()
                .all(|slots| slots.contains(&5))
        );
        for epsilon in [Rational::from((1, 13)), Rational::from((1, 17))] {
            let start = std::time::Instant::now();
            let value = prepared.evaluate(&epsilon, &settings, &context)?.remove(0);
            let elapsed = start.elapsed();
            let expected = bubble_difference(&epsilon, p)?;
            eprintln!(
                "dependent cut epsilon={epsilon} placement={:?} elapsed={elapsed:?} value={value} relative_error={}",
                settings.mass_mode,
                p.norm(&p.sub(&value, &expected)) / p.norm(&expected)
            );
            close(p, &value, &expected, 20);
            if epsilon == (1, 13) {
                if let Some(previous) = &baseline {
                    close(p, &value, previous, 20);
                }
                baseline = Some(value.clone());
                let mut refined = settings.clone();
                refined.guard_digits += 20;
                refined.series_order += 32;
                let refined = prepared.evaluate(&epsilon, &refined, &context)?.remove(0);
                close(p, &refined, &expected, 20);
                close(p, &refined, &value, 20);
            }
        }
    }
    Ok(())
}

#[test]
fn original_slot_selection_and_unresolved_regulators_are_not_reinterpreted() -> Result<()> {
    let (family, cuts) = native(true, 3, false);
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let mut options = options();
    options.mass_mode = MassMode::Propagators(vec![3]);
    assert!(matches!(
        PreparedCutCombination::from_hepkit(
            &family,
            &[1; 6],
            &Atom::one(),
            cuts.clone(),
            &channel(),
            eps(),
            &backend,
            &options,
            10000,
            &context
        ),
        Err(Error::Unsupported(_))
    ));
    options.mass_mode = MassMode::All;
    let cuts = CutMetadata::new(
        cuts.lines()
            .iter()
            .map(|(&propagator, definition)| CutLine {
                propagator,
                definition: definition.clone(),
            }),
        vec![
            LoopPrescription::Insensitive,
            LoopPrescription::Insensitive,
            LoopPrescription::MinusI0,
        ],
    )?;
    assert!(matches!(
        PreparedCutCombination::from_hepkit(
            &family,
            &[0, 1, 1, 1, 1, 1],
            &Atom::zero(),
            cuts,
            &channel(),
            eps(),
            &backend,
            &options,
            10000,
            &context
        ),
        Err(Error::Unsupported(_))
    ));
    Ok(())
}

#[test]
fn permuted_poles_and_complex_exact_numerator_preserve_distribution_normalization() -> Result<()> {
    let (family, cuts) = native(false, 2, false);
    let reversed = NativeFamily::new(
        family.loop_momenta().to_vec(),
        family.external_momenta().to_vec(),
        family.denominators().iter().rev().cloned().collect(),
        family.kinematics(),
    )
    .unwrap();
    let cuts = CutMetadata::new(
        cuts.lines()
            .iter()
            .map(|(&propagator, definition)| CutLine {
                propagator: 2 - propagator,
                definition: definition.clone(),
            }),
        vec![LoopPrescription::Insensitive],
    )?;
    let backend = RustRedBackend::default();
    let options = options();
    let context = RunContext::default();
    let p = Precision::decimal(90)?;
    let weight = Atom::num(Complex::new(Rational::one(), Rational::from((1, 3))));
    let prepared = PreparedCutCombination::from_hepkit(
        &reversed,
        &[1, 1, 2],
        &weight,
        cuts,
        &channel(),
        eps(),
        &backend,
        &options,
        10000,
        &context,
    )?;
    let epsilon = Rational::from((1, 13));
    let expected = p.mul(
        &volume(&epsilon, false, p)?,
        &p.mul(
            &p.rational(&(Rational::from((1, 4)) - &epsilon)),
            &p.eval(&weight, &Default::default())?,
        ),
    );
    close(
        p,
        &prepared.evaluate(&epsilon, &options, &context)?.remove(0),
        &expected,
        35,
    );
    Ok(())
}

#[test]
fn equal_mass_specialization_uses_raised_poles_without_generic_mass_difference() -> Result<()> {
    let (family, cuts) = native(true, 2, false);
    let parts = family.partial_fraction(&[1; 6], 10000).unwrap();
    assert_eq!(parts.len(), 1);
    assert!(parts[0].0.is_one());
    assert_eq!(parts[0].1.iter().sum::<i32>(), 6);
    let backend = RustRedBackend::default();
    let options = options();
    let context = RunContext::default();
    let prepared = PreparedCutCombination::from_hepkit(
        &family,
        &[1; 6],
        &Atom::one(),
        cuts,
        &channel(),
        eps(),
        &backend,
        &options,
        10000,
        &context,
    )?;
    assert_eq!(prepared.original_slots().len(), 1);
    let epsilon = Rational::from((1, 13));
    let a = Rational::one() - &epsilon;
    let mut sum = Rational::zero();
    let mut binomial = Rational::one();
    // Differentiate the convergent exact beta series with respect to mass a at a=2.
    for n in 0..220_i64 {
        let mut inner = Rational::zero();
        let mut choose = Rational::one();
        let mut beta = Rational::from((1, n + 1));
        let mut moment = Rational::one();
        let mut power = Rational::from(3).pow(n as u64);
        for j in 0..=n {
            inner -= &choose * &beta * &moment * &power * Rational::from((n - j, 3));
            choose *= Rational::from((n - j, j + 1));
            beta *= Rational::from((j + 1, n + j + 2));
            moment *= (&a + &Rational::from(j)) / (Rational::from(3) * &a + Rational::from(j));
            power /= Rational::from(3);
        }
        sum += &binomial * &inner;
        binomial *= (&epsilon + &Rational::from(n)) / Rational::from(5 * (n + 1));
    }
    let p = Precision::decimal(90)?;
    let expected = p.mul(
        &volume(&epsilon, true, p)?,
        &p.mul(
            &p.gamma_real(&p.rational(&epsilon).re)?,
            &p.mul(
                &p.pow(&p.i(5), &p.neg(&p.rational(&epsilon))),
                &p.rational(&sum),
            ),
        ),
    );
    close(
        p,
        &prepared.evaluate(&epsilon, &options, &context)?.remove(0),
        &expected,
        20,
    );
    Ok(())
}

#[test]
fn incompatible_original_quadratic_forms_fail_even_for_zero_targets() -> Result<()> {
    let (family, cuts) = native(true, 3, false);
    let backend = RustRedBackend::default();
    let options = options();
    let context = RunContext::default();
    let mut denominators = family.denominators().to_vec();
    denominators[3] = -denominators[3].clone();
    let negative = NativeFamily::new(
        family.loop_momenta().to_vec(),
        family.external_momenta().to_vec(),
        denominators,
        family.kinematics(),
    )
    .unwrap();
    assert!(matches!(
        PreparedCutCombination::from_hepkit(
            &negative,
            &[0, 1, 1, 1, 1, 1],
            &Atom::zero(),
            cuts.clone(),
            &channel(),
            eps(),
            &backend,
            &options,
            10000,
            &context
        ),
        Err(Error::Unsupported(_))
    ));
    let mut lines = cuts
        .lines()
        .iter()
        .filter(|(slot, _)| **slot != 0)
        .map(|(&propagator, definition)| CutLine {
            propagator,
            definition: definition.clone(),
        })
        .collect::<Vec<_>>();
    lines.push(CutLine {
        propagator: 3,
        definition: CutDefinition::PositiveEnergy {
            momentum: MomentumRouting {
                loops: vec![0.into(), 0.into(), 1.into()],
                external: vec![0.into()],
            },
        },
    });
    let virtual_cut = CutMetadata::new(lines, cuts.loop_prescriptions().to_vec())?;
    assert!(matches!(
        PreparedCutCombination::from_hepkit(
            &family,
            &[0, 1, 1, 1, 1, 1],
            &Atom::zero(),
            virtual_cut,
            &channel(),
            eps(),
            &backend,
            &options,
            10000,
            &context
        ),
        Err(Error::Unsupported(_))
    ));
    Ok(())
}

#[test]
fn original_coefficient_domains_survive_removed_cuts_and_exact_cancellations() -> Result<()> {
    let (family, cuts) = native(false, 2, false);
    let backend = RustRedBackend::default();
    let options = options();
    let context = RunContext::default();
    let e = Atom::var(eps());
    let forbidden = Rational::from((1, 13));
    let a = Atom::num(forbidden.clone());
    let guarded_zero = (e.clone().pow(2) - a.clone().pow(2)) / (&e - &a) - &e - &a;
    assert!(guarded_zero.together().cancel().is_zero());
    for (powers, numerator) in [
        (vec![0, 1, 1], Atom::one() / (&e - &a)),
        (vec![1, 1, 1], family.denominators()[0].clone() / (&e - &a)),
        (vec![1, 1, 1], guarded_zero),
    ] {
        let prepared = PreparedCutCombination::from_hepkit(
            &family,
            &powers,
            &numerator,
            cuts.clone(),
            &channel(),
            eps(),
            &backend,
            &options,
            10000,
            &context,
        )?;
        assert!(!prepared.nonzero_conditions().is_empty());
        assert!(matches!(
            prepared.evaluate(&forbidden, &options, &context),
            Err(Error::InvalidInput(_))
        ));
        let regular = prepared.evaluate(&Rational::from((1, 17)), &options, &context)?;
        assert!(Precision::decimal(90)?.norm(&regular[0]).is_zero());
    }
    assert!(matches!(
        PreparedCutCombination::from_hepkit(
            &family,
            &[0, 1, 1],
            &Atom::var(symbol!("dependent_cut_test::unbound")),
            cuts,
            &channel(),
            eps(),
            &backend,
            &options,
            10000,
            &context
        ),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}

#[test]
fn virtual_partial_fraction_sum_cancels_uv_poles_before_laurent_reconstruction() -> Result<()> {
    let (family, cuts) = native(true, 3, false);
    let backend = RustRedBackend::default();
    let options = options();
    let context = RunContext::default();
    let prepared = PreparedCutCombination::from_hepkit(
        &family,
        &[1; 6],
        &Atom::one(),
        cuts,
        &channel(),
        eps(),
        &backend,
        &options,
        10000,
        &context,
    )?;
    let started = std::time::Instant::now();
    let fitted = prepared.solve(0, &options, &context)?.remove(0);
    eprintln!(
        "dependent cut Laurent fit elapsed={:?} samples={} validation_samples={} refinements={}",
        started.elapsed(),
        fitted.samples,
        fitted.validation_samples,
        fitted.refinements
    );
    assert_eq!(fitted.verified_digits, Some(20));
    // At eps=0, Gamma(eps)*(eps)_n/n! tends to1/n for n>=1.
    // The n=0 UV residue cancels between the two virtual bubble children.
    let mut sum = Rational::zero();
    let mut five = Rational::one();
    for n in 1..220_i64 {
        five *= Rational::from(5);
        let mut inner = Rational::zero();
        let mut choose = Rational::one();
        let mut beta = Rational::from((1, n + 1));
        let mut two = Rational::from(2).pow(n as u64);
        let mut three = Rational::from(3).pow(n as u64);
        for j in 0..=n {
            inner += &choose * &beta * Rational::from((2, (j + 1) * (j + 2))) * (&two - &three);
            choose *= Rational::from((n - j, j + 1));
            beta *= Rational::from((j + 1, n + j + 2));
            two /= Rational::from(2);
            three /= Rational::from(3);
        }
        sum += inner / (Rational::from(n) * &five);
    }
    let p = Precision::decimal(90)?;
    let expected = p.mul(&volume(&Rational::zero(), true, p)?, &p.rational(&sum));
    close(p, &fitted.coefficients[&0], &expected, 20);
    for (power, coefficient) in &fitted.coefficients {
        if *power < 0 {
            close(p, coefficient, &p.zero(), 20);
        }
    }
    Ok(())
}
