use symbolica::prelude::*;
use symbolica_amflow::{cuts::*, *};

fn family(three_body: bool) -> Result<CutFamily> {
    let gram = vec![vec![Atom::one()]];
    let (loops, momenta, virtuals, completion) = if three_body {
        (
            vec!["r", "l", "k"],
            vec![(vec![1, 0, 0], 0), (vec![0, 1, 0], 0), (vec![-1, -1, 0], 1)],
            vec![vec![0, 0, 1], vec![-1, -1, 1]],
            vec![
                (vec![-1, 0, 0], 1),
                (vec![0, -1, 0], 1),
                (vec![-1, 0, 1], 0),
                (vec![0, 0, -1], 1),
            ],
        )
    } else {
        (
            vec!["r", "k"],
            vec![(vec![1, 0], 0), (vec![-1, 0], 1)],
            vec![vec![0, 1], vec![-1, 1]],
            vec![(vec![0, -1], 1)],
        )
    };
    let mut propagators = momenta
        .iter()
        .map(|(q, e)| Propagator::quadratic(q, &[*e], Atom::zero(), &gram))
        .collect::<Result<Vec<_>>>()?;
    for (q, mass) in virtuals.iter().zip([2, 3]) {
        propagators.push(Propagator::quadratic(q, &[0], Atom::num(mass), &gram)?);
    }
    let physical_propagators = propagators.len();
    for (q, e) in completion {
        propagators.push(Propagator::quadratic(&q, &[e], Atom::zero(), &gram)?);
    }
    let mut prescriptions = vec![LoopPrescription::Insensitive; loops.len()];
    *prescriptions.last_mut().unwrap() = LoopPrescription::PlusI0;
    let cuts = CutMetadata::new(
        momenta
            .into_iter()
            .enumerate()
            .map(|(propagator, (q, e))| CutLine {
                propagator,
                definition: CutDefinition::PositiveEnergy {
                    momentum: MomentumRouting {
                        loops: q.into_iter().map(Rational::from).collect(),
                        external: vec![e.into()],
                    },
                },
            }),
        prescriptions,
    )?;
    CutFamily::new(
        IntegralFamily {
            name: format!("partial_cut_{}", loops.len()),
            loops: loops.into_iter().map(str::to_owned).collect(),
            external: vec!["P".into()],
            external_gram: gram,
            propagators,
            physical_propagators,
            epsilon: symbol!("partial_cut::eps"),
            dimension: 4,
        },
        cuts,
    )
}

// Test-only independent Feynman-parameter average. For three bodies t=(r+l)^2
// has Beta(1-eps,2-2eps) density. The virtual bubble is
// Gamma(eps) integral_0^1 [3-u-u(1-u)t]^(-eps) du.
// Expanding about 3 has 0 <= delta <= 1/3, and for 0<eps<1 the
// normalized tail is bounded by 3^(-N-1)/(1-1/3).
fn analytic(eps: &Rational, three_body: bool, p: Precision) -> Result<ComplexFloat> {
    let a = Rational::one() - eps;
    let gamma = |x: Rational| p.gamma_real(&p.rational(&x).re);
    let pi = ComplexFloat::new(p.real(1).pi(), p.real(0));
    let phi2 = p.div(
        &p.mul(
            &p.pow(&p.scale(&pi, 4, 1), &p.rational(eps)),
            &gamma(a.clone())?,
        ),
        &p.mul(&p.scale(&pi, 8, 1), &gamma(Rational::from(2) * &a)?),
    );
    let volume = if three_body {
        p.mul(
            &p.div(&p.mul(&phi2, &phi2), &p.scale(&pi, 2, 1)),
            &p.div(
                &p.mul(&gamma(a.clone())?, &gamma(Rational::from(2) * &a)?),
                &gamma(Rational::from(3) * &a)?,
            ),
        )
    } else {
        phi2
    };
    let mut sum = Rational::zero();
    let mut binomial = Rational::one(); // (eps)_n / (n! 3^n)
    for n in 0..140_i64 {
        let mut inner = Rational::zero();
        let mut choose = Rational::one();
        let mut beta = Rational::from((1, n + 1)); // B(n+1,j+1)
        let mut moment = Rational::one();
        for j in 0..=n {
            if three_body || j == 0 {
                inner += &choose * &beta * &moment;
            }
            choose *= Rational::from((n - j, j + 1));
            beta *= Rational::from((j + 1, n + j + 2));
            moment *= (&a + &Rational::from(j)) / (Rational::from(3) * &a + Rational::from(j));
        }
        sum += &binomial * &inner;
        binomial *= (eps + &Rational::from(n)) / Rational::from(3 * (n + 1));
    }
    Ok(p.mul(
        &volume,
        &p.mul(
            &gamma(eps.clone())?,
            &p.mul(&p.pow(&p.i(3), &p.neg(&p.rational(eps))), &p.rational(&sum)),
        ),
    ))
}

fn compare(three_body: bool) -> Result<()> {
    let family = family(three_body)?;
    let mut indices = vec![0; family.family().propagators.len()];
    indices[..family.family().physical_propagators].fill(1);
    let target = Integral(indices);
    let channel = FutureTimelikeChannel {
        external: vec![1.into()],
    };
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let p = Precision::decimal(90)?;
    let mut baseline = None;
    for mass_mode in [
        MassMode::Propagators(vec![family.cuts().lines().len()]),
        MassMode::All,
    ] {
        let options = FlowOptions {
            digits: 20,
            guard_digits: 50,
            series_order: 80,
            mass_mode,
            ..Default::default()
        };
        let started = std::time::Instant::now();
        let flow = PreparedCutFlow::new(
            &family,
            &channel,
            std::slice::from_ref(&target),
            &KinematicPoint::default(),
            &backend,
            &options,
            &context,
        )?;
        eprintln!(
            "cut benchmark L={} placement={:?} preparation={:?}",
            family.family().loops.len(),
            options.mass_mode,
            started.elapsed()
        );
        for eps in [Rational::from((1, 13)), Rational::from((1, 17))] {
            let started = std::time::Instant::now();
            let actual = flow.evaluate(&eps, &options, &context)?.remove(0);
            eprintln!(
                "cut benchmark L={} placement={:?} epsilon={eps} evaluation={:?}",
                family.family().loops.len(),
                options.mass_mode,
                started.elapsed()
            );
            let expected = analytic(&eps, three_body, p)?;
            eprintln!(
                "cut comparison epsilon={eps} value={actual} relative_error={}",
                p.norm(&p.sub(&actual, &expected)) / p.norm(&expected)
            );
            assert!(
                p.norm(&p.sub(&actual, &expected)) < p.tolerance(20) * p.norm(&expected),
                "{actual} != {expected}"
            );
            if eps == (1, 13) && matches!(options.mass_mode, MassMode::Propagators(_)) {
                let mut refined = options.clone();
                refined.guard_digits += 20;
                refined.series_order += 32;
                let repeated = flow.evaluate(&eps, &refined, &context)?.remove(0);
                assert!(p.norm(&p.sub(&repeated, &expected)) < p.tolerance(20) * p.norm(&expected));
                assert!(p.norm(&p.sub(&repeated, &actual)) < p.tolerance(20) * p.norm(&expected));
            }
            if eps == (1, 13) {
                if let Some(value) = &baseline {
                    assert!(p.close(&actual, value, 20));
                }
                baseline = Some(actual);
            }
        }
    }
    Ok(())
}

#[test]
fn partial_mass_placement_retains_nonscaleless_virtual_soft_family() -> Result<()> {
    compare(false)
}

#[test]
fn nonfactorized_three_body_bubble_matches_independent_beta_series() -> Result<()> {
    compare(true)
}
