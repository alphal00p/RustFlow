use crate::algebraic::{AlgebraicSystem, SquareRoot};
use crate::diffexp::{EpsilonBoundary, EpsilonSystem};
use crate::kinematics::KinematicSystem;
use crate::transport_cache::*;
use crate::*;
use std::collections::BTreeMap;
use symbolica::prelude::*;

#[test]
fn noncommuting_physical_blocks_have_a_finite_dyson_bound() -> Result<()> {
    let p = Precision::decimal(70)?;
    let x = symbol!("finite_dyson::x");
    // B=[[0,1],[0,0]], C=[[0,0],[1,0]] do not commute.
    let system = EpsilonSystem {
        variable: x,
        matrices: vec![
            vec![vec![Atom::new(); 2]; 2],
            vec![
                vec![Atom::new(), Atom::one()],
                vec![Atom::var(x), Atom::new()],
            ],
            vec![vec![Atom::new(); 2]; 2],
        ],
    }
    .compile(p, &Default::default())?;
    let weights = [1, 2, 3, 4, 5, 6].map(|v| p.real(v));
    assert_eq!(system.epsilon_product_limit(), Some(3));
    let integral = system.error_norm_integral_weighted(&p.zero(), &p.i(1), &weights)?;
    assert!(p.close(
        &ComplexFloat::new(integral.clone(), p.real(0)),
        &p.rational(&Rational::from((4, 5))),
        60
    ));
    let bound = system.error_amplification_weighted(&p.zero(), &p.i(1), &weights)?;
    assert!(p.close(
        &ComplexFloat::new(bound.clone(), p.real(0)),
        &p.rational(&Rational::from((53, 25))),
        60
    ));
    assert!(bound < p.exp(&ComplexFloat::new(integral, p.real(0))).re);
    let result = system.transport(
        &EpsilonBoundary {
            point: p.zero(),
            leading: 0,
            coefficients: vec![
                vec![p.i(1), p.i(2)],
                vec![p.i(3), p.i(4)],
                vec![p.i(5), p.i(6)],
            ],
        },
        &[p.i(1)],
        &FlowOptions {
            digits: 40,
            guard_digits: 30,
            series_order: 16,
            ..Default::default()
        },
        &RunContext::default(),
        false,
    )?;
    // Exact iterated integrals contain BC/6 + CB/3, not (BC+CB)/4.
    let expected = [
        p.i(1),
        p.i(2),
        p.i(5),
        p.rational(&Rational::from((9, 2))),
        p.rational(&Rational::from((55, 6))),
        p.rational(&Rational::from((49, 6))),
    ];
    for ((actual, expected), weight) in result
        .coefficients
        .iter()
        .flatten()
        .zip(expected)
        .zip(weights)
    {
        assert!(p.close(actual, &expected, 60));
        assert!(p.norm(actual) / weight <= bound);
    }
    Ok(())
}

#[test]
fn registered_root_positive_epsilon_kernel_uses_same_finite_bound() -> Result<()> {
    let p = Precision::decimal(70)?;
    let r = symbol!("finite_dyson::r");
    let system = AlgebraicSystem {
        system: EpsilonSystem {
            variable: symbol!("finite_dyson::x"),
            matrices: vec![
                vec![vec![Atom::new()]],
                vec![vec![Atom::num(100) / Atom::var(r)]],
                vec![vec![Atom::new()]],
            ],
        },
        roots: vec![SquareRoot {
            symbol: r,
            radicand: Atom::num(4),
        }],
        nonzero_conditions: Vec::new(),
    }
    .compile(p)?;
    assert_eq!(system.epsilon_product_limit(), Some(3));
    let bound = system.error_amplification_weighted(&p.zero(), &p.i(1), &vec![p.real(1); 3])?;
    assert!(p.close(&ComplexFloat::new(bound, p.real(0)), &p.i(1301), 60));
    Ok(())
}

#[test]
fn whole_path_integral_is_subdivision_independent_for_constant_kernels() -> Result<()> {
    let p = Precision::decimal(70)?;
    let system = EpsilonSystem {
        variable: symbol!("finite_dyson::x"),
        matrices: vec![
            vec![vec![Atom::new()]],
            vec![vec![Atom::num(100)]],
            vec![vec![Atom::new()]],
            vec![vec![Atom::new()]],
            vec![vec![Atom::new()]],
        ],
    }
    .compile(p, &Default::default())?;
    let weights = vec![p.real(1); 5];
    let whole = system.error_amplification_weighted(&p.zero(), &p.i(1), &weights)?;
    let mut integral = p.real(0);
    for index in 0..16 {
        integral += system.error_norm_integral_weighted(
            &p.rational(&Rational::from((index, 16))),
            &p.rational(&Rational::from((index + 1, 16))),
            &weights,
        )?;
    }
    let subdivided =
        crate::diffexp::amplification_from_integral(p, &integral, system.epsilon_product_limit())?;
    assert_eq!(whole, subdivided);
    assert!(p.close(
        &ComplexFloat::new(whole.clone(), p.real(0)),
        &p.rational(&Rational::from((13_015_303, 3))),
        60
    ));
    // Multiplying a truncated scalar bound at every small segment would tend
    // back toward exp(100) and would lose the global nilpotent structure.
    assert!(whole < p.real(5_000_000));
    Ok(())
}

#[test]
fn arbitrarily_small_nonzero_epsilon_zero_block_keeps_exponential_bound() -> Result<()> {
    let p = Precision::decimal(70)?;
    let tiny = Atom::one() / Atom::num(10).pow(400);
    let system = EpsilonSystem {
        variable: symbol!("finite_dyson::x"),
        matrices: vec![vec![vec![tiny]], vec![vec![Atom::num(100)]]],
    }
    .compile(p, &Default::default())?;
    assert_eq!(system.epsilon_product_limit(), None);
    let bound = system.error_amplification_weighted(&p.zero(), &p.i(1), &vec![p.real(1); 2])?;
    assert!(p.close(&ComplexFloat::new(bound, p.real(0)), &p.exp(&p.i(100)), 60));
    Ok(())
}

#[test]
fn a_single_retained_epsilon_coefficient_has_no_positive_order_feedback() -> Result<()> {
    let p = Precision::decimal(70)?;
    let system = EpsilonSystem {
        variable: symbol!("finite_dyson::x"),
        matrices: vec![vec![vec![Atom::new()]]],
    }
    .compile(p, &Default::default())?;
    assert_eq!(system.epsilon_product_limit(), Some(1));
    assert_eq!(
        system.error_amplification(&p.zero(), &p.complex(100, 100))?,
        p.real(1)
    );
    Ok(())
}

#[test]
fn long_analytic_cache_transport_preserves_source_evidence_with_finite_growth() -> Result<()> {
    let p = Precision::decimal(80)?;
    let s = symbol!("finite_dyson_cache::s");
    let eps = symbol!("finite_dyson_cache::eps");
    let flow = RustFlow::new(
        KinematicSystem {
            epsilon: eps,
            derivatives: BTreeMap::from([(
                s,
                vec![vec![
                    Atom::num(100) * Atom::var(eps) / (Atom::var(s) + Atom::num(2)),
                ]],
            )]),
        },
        &[Atom::var(symbol!("finite_dyson_cache::Y"))],
        &Atom::one(),
        Prescription::PlusI0,
        "real positive s+2",
    )?;
    let range = EpsilonRange::new(0, 4)?;
    let source = CachedBoundary {
        identity: flow.identity().clone(),
        point: CachedPoint::Exact(BTreeMap::from([(s, Atom::new())])),
        kind: PointKind::Physical,
        range,
        coefficients: vec![
            vec![p.i(1)],
            vec![p.zero()],
            vec![p.zero()],
            vec![p.zero()],
            vec![p.zero()],
        ],
        accuracy: BoundaryAccuracy::supplied(
            36,
            p.bits,
            vec![vec![p.tolerance(44)]; 5],
            "independent analytic source, conservatively36 digits",
        )?,
    };
    let mut bank = RustFlowCache::default();
    bank.insert(source)?;
    let options = FlowOptions {
        digits: 20,
        guard_digits: 30,
        series_order: 64,
        ..Default::default()
    };
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    let result = flow.evaluate_to(
        &mut bank,
        &BTreeMap::from([(s, Atom::one())]),
        range,
        &options,
        &RunContext::default(),
        &policy,
    )?;
    assert!(result.transport.as_ref().unwrap().segments.len() > 1);
    assert!(result.boundary.accuracy.verified_digits() >= 20);
    assert!(result.boundary.accuracy.verified_digits() <= 36);
    let log = p.scale(&p.log(&p.rational(&Rational::from((3, 2)))), 100, 1);
    let mut expected = p.i(1);
    for (order, row) in result.boundary.coefficients.iter().enumerate() {
        if order > 0 {
            expected = p.scale(&p.mul(&expected, &log), 1, order as i64);
        }
        assert!(p.close(&row[0], &expected, 20));
    }
    assert!(result.boundary.accuracy.comparison_errors()[0][0] < p.tolerance(28));
    assert!(bank.len() > 1);
    // exp(integral||A||) is at least exp(100 ln(3/2)); that obsolete
    // all-orders estimate cannot support20 absolute digits for epsilon0.
    assert!(p.exp(&log).re * p.tolerance(36) > p.tolerance(20));
    Ok(())
}
