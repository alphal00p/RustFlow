//! Public compile regressions for a genuinely irreducible cubic and known-root
//! polynomials kept as one formal factor until exact parameter specialization.
use symbolica::prelude::*;
use symbolica_amflow::{ComplexFloat, DifferentialSystem, Error, Precision, Result};

fn variable() -> Symbol {
    symbol!("pole_scale::t")
}

fn ten_power(n: u64) -> Rational {
    Rational::from(Integer::from(10).pow(n))
}

fn compile(p: Precision, polynomial: Atom) -> Result<Vec<ComplexFloat>> {
    let auxiliary = Atom::var(symbol!("pole_scale::auxiliary"));
    let system = DifferentialSystem {
        variable: variable(),
        matrix: vec![vec![Atom::one() / (polynomial + &auxiliary)]],
    };
    Ok(system
        .compile(p, &[(auxiliary, p.zero())].into_iter().collect())?
        .poles)
}

#[test]
fn coprime_far_linear_poles_retain_the_exact_unit_gap() -> Result<()> {
    let p = Precision { bits: 256 };
    let t = Atom::var(variable());
    for n in [1, 20, 40, 80, 100] {
        let h = ten_power(n);
        let next = &h + &Rational::one();
        let compiled = DifferentialSystem {
            variable: variable(),
            matrix: vec![
                vec![Atom::one() / (&t - Atom::num(h.clone())), Atom::new()],
                vec![Atom::new(), Atom::one() / (&t - Atom::num(next.clone()))],
            ],
        }
        .compile(p, &Default::default())?;
        assert_eq!(compiled.poles.len(), 2);
        // Exact integer comparison, even when H+1 cannot be stored at p bits.
        let actual = compiled
            .poles
            .iter()
            .map(|root| {
                assert!(root.im.is_zero());
                root.re.to_rational()
            })
            .collect::<Vec<_>>();
        assert!(actual.contains(&h));
        assert!(actual.contains(&next));
    }
    Ok(())
}

#[test]
fn exact_coprime_gap_beyond_geometry_budget_still_fails_closed() {
    let p = Precision { bits: 256 };
    let t = Atom::var(variable());
    let delta = Rational::from((Integer::one(), Integer::from(2).pow(5000)));
    let system = DifferentialSystem {
        variable: variable(),
        matrix: vec![
            vec![Atom::one() / (&t - 2), Atom::new()],
            vec![
                Atom::new(),
                Atom::one() / (&t - Atom::num(Rational::from(2) + delta)),
            ],
        ],
    };
    assert!(
        matches!(system.compile(p, &Default::default()), Err(Error::Numerical(message))
        if message.contains("degree-1 geometry at 4352 bits")
            && message.contains("overlapping exact-coefficient isolation disks")
            && message.contains("bounded extra precision exhausted"))
    );
}

#[test]
fn extra_domains_share_exact_poles_without_merging_close_coprime_ones() -> Result<()> {
    let p = Precision { bits: 256 };
    let t = Atom::var(variable());
    let h = ten_power(100);
    let first = &t - Atom::num(h.clone());
    let second = &t - Atom::num(h + Rational::one());
    let mut compiled = DifferentialSystem {
        variable: variable(),
        matrix: vec![vec![Atom::one() / &first]],
    }
    .compile(p, &Default::default())?;
    compiled.exclude_polynomials(&[second.clone()])?;
    assert_eq!(compiled.poles.len(), 2);
    assert_ne!(compiled.poles[0], compiled.poles[1]);
    compiled.exclude_polynomials(&[first.pow(2), second.pow(3)])?;
    assert_eq!(compiled.poles.len(), 2);
    Ok(())
}

#[test]
fn scaled_far_pole_physical_transport_preserves_finite_source_floors() -> Result<()> {
    use std::collections::BTreeMap;
    use symbolica_amflow::kinematics::KinematicSystem;
    use symbolica_amflow::physical_transport::PhysicalTransport;
    use symbolica_amflow::transport_cache::*;
    use symbolica_amflow::{FlowOptions, Prescription, RunContext};

    let s = symbol!("pole_scale::s");
    let x = Atom::var(s);
    let h = ten_power(100);
    let next = &h + &Rational::one();
    // Y'_i = H_i/(H_i-s) Y_i, Y_i(0)=1. Both solutions change by O(1)
    // on [0,1], despite poles at H and H+1. Their exact endpoints are
    // exp(sum_{k>=1} 1/(k H_i^(k-1))); 0 < log(Y_i(1))-1 < 1/H_i.
    // Thus exp(1) is an independent oracle within 3/H_i, far below every
    // requested comparison here; no subtraction of nearly equal logs is used.
    let engine = PhysicalTransport::new(
        KinematicSystem {
            epsilon: symbol!("pole_scale::eps"),
            derivatives: BTreeMap::from([(
                s,
                vec![
                    vec![Atom::num(h.clone()) / (Atom::num(h) - &x), Atom::new()],
                    vec![
                        Atom::new(),
                        Atom::num(next.clone()) / (Atom::num(next) - &x),
                    ],
                ],
            )]),
        },
        &[parse!("pole_scale::I1"), parse!("pole_scale::I2")],
        &Atom::one(),
        Prescription::PlusI0,
        "real regular interval [0,1] below both far poles",
    )?;
    for (bits, verified, digits) in [(80, 16, 10), (201, 48, 32)] {
        let storage = Precision { bits };
        let mut cache = RustFlowCache::default();
        cache.insert(CachedBoundary {
            identity: engine.identity().clone(),
            point: CachedPoint::Exact(BTreeMap::from([(s, Atom::new())])),
            kind: PointKind::Physical,
            range: EpsilonRange::new(0, 0)?,
            coefficients: vec![vec![storage.i(1), storage.i(1)]],
            accuracy: BoundaryAccuracy::supplied(
                verified,
                bits,
                vec![vec![storage.real(0); 2]],
                "analytic unit seed with finite declared source evidence",
            )?,
        })?;
        let result = engine.evaluate_to(
            &mut cache,
            &BTreeMap::from([(s, Atom::one())]),
            EpsilonRange::new(0, 0)?,
            &FlowOptions {
                digits,
                guard_digits: 20,
                series_order: 64,
                ..Default::default()
            },
            &RunContext::default(),
            &ScaledDistance {
                scales: BTreeMap::new(),
                admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
            },
        )?;
        let accuracy = &result.boundary.accuracy;
        assert_eq!(accuracy.input_verified_digits(), verified);
        assert!(accuracy.verified_digits() <= verified);
        let p = Precision {
            bits: accuracy.working_bits(),
        };
        for value in &result.boundary.coefficients[0] {
            assert!(p.close(value, &p.exp(&p.i(1)), digits));
            assert!(value.re > p.real(2));
        }
        assert!(
            accuracy
                .comparison_errors()
                .iter()
                .flatten()
                .all(|error| *error >= p.tolerance(verified))
        );
    }
    Ok(())
}

/// Relative comparison even below unit scale; each expected physical pole
/// must occur exactly once. A fixed absolute comparison would miss this bug.
fn expect_roots(p: Precision, polynomial: Atom, expected: &[ComplexFloat]) -> Result<()> {
    let actual = compile(p, polynomial)?;
    assert_eq!(actual.len(), expected.len());
    for root in &actual {
        assert_eq!(root.re.prec(), p.bits);
        assert_eq!(root.im.prec(), p.bits);
    }
    for root in expected {
        assert_eq!(
            actual
                .iter()
                .filter(|candidate| {
                    p.norm(&p.sub(candidate, root)) <= p.tolerance(48) * p.norm(root)
                })
                .count(),
            1,
            "missing or duplicate physical pole {root:?}; found {actual:?}"
        );
    }
    Ok(())
}

#[test]
fn irreducible_cubic_reproduces_the_multiscale_failure() -> Result<()> {
    let p = Precision { bits: 256 };
    let t = Atom::var(variable());
    for n in [1, 20, 40, 80, 100] {
        let h = ten_power(n);
        // Frozen e732 passes N1/20/40 and fails N80/100 here. This cubic
        // itself is irreducible; no auxiliary parameter is needed.
        let polynomial =
            t.pow(3) / Atom::num(h.clone()) - Atom::num(3) * t.pow(2) + Atom::num(2) * &t - 1;
        let actual = DifferentialSystem {
            variable: variable(),
            matrix: vec![vec![Atom::one() / polynomial]],
        }
        .compile(p, &Default::default())?
        .poles;
        assert_eq!(actual.len(), 3);
        let h = p.rational(&h);
        let inverse_h = p.div(&p.i(1), &h);
        let upper = p.scale(&ComplexFloat::new(p.real(1), p.real(2).sqrt()), 1, 3);
        let lower = ComplexFloat::new(upper.re.clone(), -upper.im.clone());
        let large = p.sub(&p.scale(&h, 3, 1), &p.rational(&Rational::from((2, 3))));
        // For H>=10, Rouche on the disks of radius 1/H around
        // (1 +/- i*sqrt(2))/3 places one root in each: on their boundaries
        // |3t^2-2t+1| >= (2sqrt(2)-3/H)/H > |t^3/H|.
        // Vieta then bounds the third root's distance from 3H-2/3 by 2/H.
        // Add only the fixed comparison tolerance for the returned p-bit roots.
        for (expected, analytic_radius) in [
            (upper, inverse_h.clone()),
            (lower, inverse_h.clone()),
            (large, p.scale(&inverse_h, 2, 1)),
        ] {
            let allowed = analytic_radius.re + p.tolerance(48) * p.norm(&expected);
            assert_eq!(
                actual
                    .iter()
                    .filter(|root| p.norm(&p.sub(root, &expected)) <= allowed)
                    .count(),
                1,
                "wrong physical pole at N={n}: {actual:?}"
            );
        }
        for root in actual {
            assert_eq!(root.re.prec(), p.bits);
            assert_eq!(root.im.prec(), p.bits);
        }
    }
    Ok(())
}

#[test]
fn physical_poles_one_two_and_ten_to_n_remain_distinct() -> Result<()> {
    let p = Precision { bits: 256 };
    let t = Atom::var(variable());
    for n in [1, 20, 40, 80, 100] {
        let h = ten_power(n);
        let polynomial = (&t - 1) * (&t - 2) * (&t / Atom::num(h.clone()) - 1);
        expect_roots(p, polynomial, &[p.i(1), p.i(2), p.rational(&h)])?;
    }
    Ok(())
}

#[test]
fn reciprocal_and_uniformly_tiny_poles_use_relative_checks() -> Result<()> {
    let p = Precision { bits: 256 };
    let t = Atom::var(variable());
    for n in [20, 40, 80, 100] {
        let h = ten_power(n);
        let ht = Atom::num(h.clone()) * &t;
        let tiny = p.rational(&(Rational::one() / &h));
        expect_roots(
            p,
            (&t - 1) * (Atom::num(2) * &t - 1) * (&ht - 1),
            &[p.i(1), p.scale(&p.i(1), 1, 2), tiny.clone()],
        )?;
        expect_roots(
            p,
            (&ht - 1) * (&ht - 2) * (&ht - 3),
            &[tiny.clone(), p.scale(&tiny, 2, 1), p.scale(&tiny, 3, 1)],
        )?;
    }
    Ok(())
}

#[test]
fn zero_and_complex_poles_keep_their_physical_locations() -> Result<()> {
    let p = Precision { bits: 256 };
    let t = Atom::var(variable());
    let h = ten_power(100);
    let tiny = p.rational(&(Rational::one() / &h));
    expect_roots(
        p,
        &t * (Atom::num(h.clone()) * &t - 1) * (&t - Atom::num(h.clone())),
        &[p.zero(), tiny, p.rational(&h)],
    )?;
    let j = Atom::num(Complex::new(Rational::zero(), Rational::one()));
    expect_roots(
        p,
        (&t - &j) * (&t - Atom::num(2) * &j) * (&t / Atom::num(h.clone()) - &j),
        &[
            p.complex(0, 1),
            p.complex(0, 2),
            p.mul(&p.complex(0, 1), &p.rational(&h)),
        ],
    )?;
    Ok(())
}

#[test]
fn exact_multiplicity_is_removed_before_the_guarded_solve() -> Result<()> {
    let p = Precision { bits: 256 };
    let t = Atom::var(variable());
    let h = ten_power(100);
    expect_roots(
        p,
        (&t - 1).pow(2) * (&t - 2).pow(3) * (&t / Atom::num(h.clone()) - 1).pow(2),
        &[p.i(1), p.i(2), p.rational(&h)],
    )?;
    Ok(())
}

#[test]
fn a_resolvable_cluster_keeps_both_poles() -> Result<()> {
    let p = Precision { bits: 256 };
    let t = Atom::var(variable());
    let delta = Rational::from((Integer::one(), Integer::from(2).pow(60)));
    let nearby = Rational::one() + delta;
    expect_roots(
        p,
        (&t - 1) * (&t - Atom::num(nearby.clone())),
        &[p.i(1), p.rational(&nearby)],
    )?;
    Ok(())
}

#[test]
fn coefficient_rounding_cannot_manufacture_resolved_poles() {
    let p = Precision { bits: 256 };
    let t = Atom::var(variable());
    let tiny = Rational::from((Integer::one(), Integer::from(2).pow(700)));
    // Exact square-free roots are 1 +/- i*2^-350. At p bits the tiny
    // constant perturbation disappears. Extra root precision cannot fix it.
    let unresolved = (&t - 1).pow(2) + Atom::num(tiny);
    assert!(matches!(compile(p, unresolved), Err(Error::Numerical(_))));

    // This split survives coefficient rounding, but the coefficient-induced
    // forward uncertainty is too large for the requested root tolerance.
    let delta = Rational::from((Integer::one(), Integer::from(2).pow(100)));
    let nearby = Rational::one() + delta;
    assert!(matches!(
        compile(p, (&t - 1) * (&t - Atom::num(nearby))),
        Err(Error::Numerical(_))
    ));
}

#[test]
fn extreme_root_spread_hits_a_bounded_precision_error() {
    let p = Precision { bits: 256 };
    let t = Atom::var(variable());
    let huge = Atom::num(ten_power(2000));
    let result = compile(p, (&t - 1) * (&t - 2) * (&t / huge - 1));
    assert!(matches!(result, Err(Error::Numerical(message))
        if message.contains("bounded extra precision")));
}

#[test]
fn exact_zero_and_tiny_linear_poles_need_no_absolute_radius_floor() -> Result<()> {
    let p = Precision { bits: 256 };
    let t = Atom::var(variable());
    let tiny = Rational::from((Integer::one(), Integer::from(2).pow(5000)));
    let compiled = DifferentialSystem {
        variable: variable(),
        matrix: vec![
            vec![Atom::one() / &t, Atom::new()],
            vec![Atom::new(), Atom::one() / (&t - Atom::num(tiny.clone()))],
        ],
    }
    .compile(p, &Default::default())?;
    assert_eq!(compiled.poles.len(), 2);
    let roots = compiled
        .poles
        .iter()
        .map(|root| {
            assert!(root.im.is_zero());
            root.re.to_rational()
        })
        .collect::<Vec<_>>();
    assert!(roots.contains(&Rational::zero()));
    assert!(roots.contains(&tiny));
    Ok(())
}
