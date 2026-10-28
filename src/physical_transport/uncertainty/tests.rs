use super::*;
use crate::diffexp::EpsilonSystem;
use crate::kinematics::KinematicSystem;
use crate::local_coordinates::TaylorCoordinate;
use crate::transport_cache::{BoundaryAccuracy, CachedPoint, EpsilonRange, PointKind};
use crate::{ComplexFloat as C, FlowDiagnostics, Prescription, RustFlow};
use std::collections::BTreeMap;

fn fixture(
    a0: Atom,
    a1: Atom,
    initial: [i64; 5],
    digits: u32,
) -> Result<(CompiledConnection, CachedBoundary)> {
    let p = Precision::decimal(90)?;
    let x = symbol!("boundary_disks::x");
    let epsilon = symbol!("boundary_disks::epsilon");
    let flow = RustFlow::new(
        KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(x, vec![vec![&a0 + Atom::var(epsilon) * &a1]])]),
        },
        &[Atom::var(symbol!("boundary_disks::I"))],
        &Atom::one(),
        Prescription::PlusI0,
        "regular test interval",
    )?;
    let compiled = EpsilonSystem {
        variable: x,
        matrices: vec![
            vec![vec![a0]],
            vec![vec![a1]],
            vec![vec![Atom::new()]],
            vec![vec![Atom::new()]],
            vec![vec![Atom::new()]],
        ],
    }
    .compile(p, &Default::default())?;
    let source = CachedBoundary {
        identity: flow.identity().clone(),
        point: CachedPoint::Exact(BTreeMap::from([(x, Atom::new())])),
        kind: PointKind::Physical,
        range: EpsilonRange::new(-2, 2)?,
        coefficients: initial.map(|v| vec![p.i(v)]).to_vec(),
        accuracy: BoundaryAccuracy::supplied(
            digits,
            p.bits,
            vec![vec![p.tolerance(60)]; 5],
            "finite analytic source evidence, independent of storage precision",
        )?,
    };
    source.validate()?;
    Ok((CompiledConnection::Rational(compiled), source))
}

// Exact polynomial solutions let the admission test retain a whole regular
// chart without involving the central integrator's choice of step size.
fn linear_solution(source: &CachedBoundary, slope: &[C], split: bool) -> EpsilonSolution {
    let p = Precision {
        bits: source.accuracy.working_bits(),
    };
    let initial = source
        .coefficients
        .iter()
        .flatten()
        .cloned()
        .collect::<Vec<_>>();
    let values = |point: &C| {
        initial
            .iter()
            .zip(slope)
            .map(|(a, b)| p.add(a, &p.mul(b, point)))
            .map(|v| vec![v])
            .collect::<Vec<_>>()
    };
    let endpoints = if split {
        vec![p.rational(&Rational::from((1, 2))), p.i(1)]
    } else {
        vec![p.i(1)]
    };
    let mut start = p.zero();
    let mut segments = Vec::new();
    let mut checkpoints = Vec::new();
    for (segment, end) in endpoints.into_iter().enumerate() {
        segments.push(TaylorSegment {
            center: start.clone(),
            end: end.clone(),
            coefficients: vec![
                values(&start).into_iter().flatten().collect(),
                slope.to_vec(),
            ],
            pade: None,
            working_bits: p.bits,
            coordinate: TaylorCoordinate::Identity,
            conditioning_digits: 70,
        });
        checkpoints.push(VerifiedCheckpoint {
            segment,
            point: end.clone(),
            coefficients: values(&end),
            comparison_errors: vec![vec![p.tolerance(70)]; 5],
        });
        start = end;
    }
    EpsilonSolution {
        point: p.i(1),
        leading: -2,
        coefficients: values(&p.i(1)),
        diagnostics: FlowDiagnostics {
            working_bits: p.bits,
            expansion_order: 80,
            steps: segments.len(),
            ..Default::default()
        },
        segments,
        verified_digits: Some(12),
        comparison_errors: vec![vec![p.tolerance(70)]; 5],
        checkpoints,
    }
}

fn problem<'a>(
    prepared: &'a CompiledConnection,
    source: &'a CachedBoundary,
    options: &'a FlowOptions,
    context: &'a RunContext,
) -> Problem<'a> {
    let p = Precision {
        bits: source.accuracy.working_bits(),
    };
    let weights = source
        .coefficients
        .iter()
        .flatten()
        .map(|v| {
            let norm = p.norm(v);
            if norm > p.real(1) { norm } else { p.real(1) }
        })
        .collect::<Vec<_>>();
    let relative = source
        .accuracy
        .comparison_errors()
        .iter()
        .flatten()
        .zip(&weights)
        .map(|(e, w)| e.clone() / w)
        .reduce(|a, b| if a > b { a } else { b })
        .unwrap();
    Problem {
        prepared,
        source,
        options,
        context,
        p,
        count: 5,
        weights,
        relative_error: relative + p.tolerance(source.accuracy.verified_digits()),
    }
}

#[test]
fn positive_loose_disks_admit_the_same_polynomial_chart_after_refinement() -> Result<()> {
    let x = Atom::var(symbol!("boundary_disks::x"));
    // y_k = initial_k*(1-99*x/100). The source pole 100/99 is outside
    // the chart, but its positive whole-disk denominator bound is only 1.
    let (prepared, source) = fixture(
        Atom::num(-99) / (100 - 99 * x),
        Atom::new(),
        [1, 2, 3, 4, 5],
        32,
    )?;
    let p = Precision {
        bits: source.accuracy.working_bits(),
    };
    let slopes = (1..=5)
        .map(|i| p.rational(&Rational::from((-99 * i, 100))))
        .collect::<Vec<_>>();
    let mut solution = linear_solution(&source, &slopes, false);
    let before = solution.clone();
    let source_before = format!("{source:?}");
    let context = RunContext::default();
    for strategy in [
        BoundaryErrorStrategy::ScalarNorm,
        BoundaryErrorStrategy::Automatic,
    ] {
        let options = FlowOptions {
            digits: 12,
            boundary_error_strategy: strategy,
            ..Default::default()
        };
        let problem = problem(&prepared, &source, &options, &context);
        assert!(matches!(
            problem.scalar_disks(&before, false),
            Err(Error::Accuracy(_))
        ));
        solution = before.clone();
        admit(&prepared, &source, &mut solution, &options, &context)?;
        assert_eq!(solution.coefficients, before.coefficients);
        assert_eq!(solution.point, before.point);
        assert_eq!(solution.leading, before.leading);
        assert_eq!(solution.verified_digits, before.verified_digits);
        assert_eq!(
            format!("{:?}", solution.segments),
            format!("{:?}", before.segments)
        );
        assert_eq!(
            format!("{:?}", solution.diagnostics),
            format!("{:?}", before.diagnostics)
        );
        assert_eq!(solution.checkpoints.len(), 1);
        assert_eq!(
            solution.checkpoints[0].segment,
            before.checkpoints[0].segment
        );
        assert_eq!(solution.checkpoints[0].point, before.checkpoints[0].point);
        assert_eq!(
            solution.checkpoints[0].coefficients,
            before.checkpoints[0].coefficients
        );
        assert_eq!(
            solution.checkpoints[0].comparison_errors,
            solution.comparison_errors
        );
        assert!(errors_meet(
            p,
            &solution.coefficients,
            &solution.comparison_errors,
            12
        ));
        for (error, weight) in solution
            .comparison_errors
            .iter()
            .flatten()
            .zip(&problem.weights)
        {
            assert!(*error > weight.clone() * p.tolerance(32));
        }
        assert_eq!(format!("{source:?}"), source_before);
    }
    Ok(())
}

#[test]
fn constant_nilpotent_hierarchy_keeps_one_global_dyson_bound() -> Result<()> {
    let (prepared, source) = fixture(Atom::new(), Atom::num(100), [0; 5], 24)?;
    let p = Precision {
        bits: source.accuracy.working_bits(),
    };
    let solution = linear_solution(&source, &vec![p.zero(); 5], true);
    let options = FlowOptions {
        digits: 12,
        ..Default::default()
    };
    let context = RunContext::default();
    let problem = problem(&prepared, &source, &options, &context);
    let coarse = problem.scalar_disks(&solution, false)?;
    let refined = problem.scalar_disks(&solution, true)?;
    assert_eq!(coarse.errors, refined.errors);
    let expected =
        p.rational(&Rational::from((13_015_303, 3))).re * &problem.relative_error + p.tolerance(70);
    assert!(p.close(
        &C::new(refined.errors[4][0].clone(), p.real(0)),
        &C::new(expected, p.real(0)),
        70
    ));
    assert_eq!(refined.checkpoints.len(), 2);
    Ok(())
}

#[test]
fn constant_growth_retains_the_source_floor_failure_and_transaction() -> Result<()> {
    let (prepared, source) = fixture(Atom::num(20), Atom::new(), [0; 5], 18)?;
    let p = Precision {
        bits: source.accuracy.working_bits(),
    };
    // The central solution is exactly zero; only the finite source floor
    // makes exp(20)*1e-18 exceed the requested 1e-12 tolerance.
    let original = linear_solution(&source, &vec![p.zero(); 5], true);
    let context = RunContext::default();
    for strategy in [
        BoundaryErrorStrategy::ScalarNorm,
        BoundaryErrorStrategy::Automatic,
    ] {
        let options = FlowOptions {
            digits: 12,
            boundary_error_strategy: strategy,
            ..Default::default()
        };
        let problem = problem(&prepared, &source, &options, &context);
        let coarse = problem.scalar_disks(&original, false).err().unwrap();
        let refined = problem.scalar_disks(&original, true).err().unwrap();
        assert_eq!(coarse.to_string(), refined.to_string());
        let mut solution = original.clone();
        let error = admit(&prepared, &source, &mut solution, &options, &context).unwrap_err();
        assert!(matches!(error, Error::Accuracy(_)));
        if strategy == BoundaryErrorStrategy::Automatic {
            assert!(
                error
                    .to_string()
                    .contains("epsilon hierarchies retain scalar bounds")
            );
        } else {
            assert_eq!(error.to_string(), coarse.to_string());
        }
        assert_eq!(format!("{solution:?}"), format!("{original:?}"));
    }
    Ok(())
}

#[test]
fn cancellation_and_invalid_bound_inputs_do_not_commit() -> Result<()> {
    let (prepared, source) = fixture(Atom::num(20), Atom::new(), [0; 5], 18)?;
    let p = Precision {
        bits: source.accuracy.working_bits(),
    };
    let options = FlowOptions {
        digits: 12,
        ..Default::default()
    };
    let context = RunContext::default();
    context.cancellation.cancel();
    let mut solution = linear_solution(&source, &vec![p.zero(); 5], true);
    let before = format!("{solution:?}");
    assert!(matches!(
        admit(&prepared, &source, &mut solution, &options, &context),
        Err(Error::Cancelled)
    ));
    assert_eq!(format!("{solution:?}"), before);
    let context = RunContext::default();
    let mut problem = problem(&prepared, &source, &options, &context);
    problem.weights[0] = p.real(-1);
    assert!(matches!(
        problem.scalar(&solution),
        Err(Error::InvalidInput(_))
    ));
    assert_eq!(format!("{solution:?}"), before);
    Ok(())
}
