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
    let dimension = source.identity.dimension();
    let count = source.coefficients.len();
    let values = |point: &C| {
        initial
            .iter()
            .zip(slope)
            .map(|(a, b)| p.add(a, &p.mul(b, point)))
            .collect::<Vec<_>>()
            .chunks(dimension)
            .map(|row| row.to_vec())
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
            comparison_errors: vec![vec![p.tolerance(70); dimension]; count],
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
        comparison_errors: vec![vec![p.tolerance(70); dimension]; count],
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
        count: source.coefficients.len(),
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
            problem.scalar_disks(&before, 1),
            Err(Error::Accuracy(_))
        ));
        let four = problem.scalar_disks(&before, 4)?;
        solution = before.clone();
        admit(&prepared, &source, &mut solution, &options, &context)?;
        assert_eq!(solution.comparison_errors, four.errors);
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
fn sixteen_disks_admit_when_four_fail_without_changing_the_central_chart() -> Result<()> {
    let x = Atom::var(symbol!("boundary_disks::x"));
    // The same exact linear solutions, with less source evidence: the positive
    // four-disk bound is still too loose, while sixteen disks admit the source.
    let (prepared, source) = fixture(
        Atom::num(-99) / (100 - 99 * x),
        Atom::new(),
        [1, 2, 3, 4, 5],
        20,
    )?;
    let p = Precision {
        bits: source.accuracy.working_bits(),
    };
    let slopes = (1..=5)
        .map(|i| p.rational(&Rational::from((-99 * i, 100))))
        .collect::<Vec<_>>();
    let original = linear_solution(&source, &slopes, false);
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
        for subdivisions in [1, 4] {
            assert!(matches!(
                problem.scalar_disks(&original, subdivisions),
                Err(Error::Accuracy(_))
            ));
        }
        let sixteen = problem.scalar_disks(&original, 16)?;
        let mut expected = original.clone();
        expected.comparison_errors = sixteen.errors;
        expected.checkpoints = sixteen.checkpoints;
        let mut solution = original.clone();
        admit(&prepared, &source, &mut solution, &options, &context)?;
        // Every central field, segment and checkpoint value remains identical;
        // only the admitted source-error candidate replaces comparison errors.
        assert_eq!(format!("{solution:?}"), format!("{expected:?}"));
        assert_eq!(solution.checkpoints.len(), original.checkpoints.len());
        assert!(errors_meet(
            p,
            &solution.coefficients,
            &solution.comparison_errors,
            12
        ));
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
    let coarse = problem.scalar_disks(&solution, 1)?;
    let expected =
        p.rational(&Rational::from((13_015_303, 3))).re * &problem.relative_error + p.tolerance(70);
    for subdivisions in [4, 16] {
        let refined = problem.scalar_disks(&solution, subdivisions)?;
        assert_eq!(coarse.errors, refined.errors);
        assert!(p.close(
            &C::new(refined.errors[4][0].clone(), p.real(0)),
            &C::new(expected.clone(), p.real(0)),
            70
        ));
        assert_eq!(refined.checkpoints.len(), 2);
    }
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
        let coarse = problem.scalar_disks(&original, 1).err().unwrap();
        for subdivisions in [4, 16] {
            let refined = problem.scalar_disks(&original, subdivisions).err().unwrap();
            assert_eq!(coarse.to_string(), refined.to_string());
        }
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

// Ordinary constant 2x2 nilpotent coupling: a large off-diagonal entry
// defeats exp(||A||), although the independent exact transfer is I + x*A.
fn coupled_fixture(
    count: usize,
    digits: u32,
    complex: bool,
) -> Result<(CompiledConnection, CachedBoundary)> {
    let p = Precision::decimal(90)?;
    let x = symbol!("boundary_balanced::x");
    let epsilon = symbol!("boundary_balanced::epsilon");
    let coupling = if complex {
        Atom::num(100) * Atom::i()
    } else {
        Atom::num(100)
    };
    let first = vec![
        vec![Atom::new(), coupling.clone()],
        vec![Atom::new(), Atom::new()],
    ];
    let zero = vec![vec![Atom::new(); 2]; 2];
    let mut matrices = vec![zero.clone(); count];
    matrices[0] = first.clone();
    if complex && count > 1 {
        matrices[1][0][1] = Atom::num(30);
    }
    let mut derivative = first;
    if complex && count > 1 {
        derivative[0][1] += Atom::num(30) * Atom::var(epsilon);
    }
    let flow = RustFlow::new(
        KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(x, derivative)]),
        },
        &[
            Atom::var(symbol!("boundary_balanced::I")),
            Atom::var(symbol!("boundary_balanced::J")),
        ],
        &Atom::one(),
        Prescription::PlusI0,
        "regular constant nilpotent interval",
    )?;
    let compiled = EpsilonSystem {
        variable: x,
        matrices,
    }
    .compile(p, &Default::default())?;
    let source = CachedBoundary {
        identity: flow.identity().clone(),
        point: CachedPoint::Exact(BTreeMap::from([(x, Atom::new())])),
        kind: PointKind::Physical,
        range: EpsilonRange::new(-2, -2 + count as i32 - 1)?,
        coefficients: vec![vec![p.zero(), p.i(1)]; count],
        accuracy: BoundaryAccuracy::supplied(
            digits,
            p.bits,
            vec![vec![p.tolerance(49); 2]; count],
            "independent finite source comparison",
        )?,
    };
    source.validate()?;
    Ok((CompiledConnection::Rational(compiled), source))
}

#[test]
fn balanced_constant_coupling_admits_nonzero_uncertainty_at_unchanged_24_digits() -> Result<()> {
    let (prepared, source) = coupled_fixture(1, 48, false)?;
    let p = Precision {
        bits: source.accuracy.working_bits(),
    };
    // y(0)=(0,1), y(x)=(100*x,1); the source perturbation
    // (delta,delta) independently transfers to (101*delta,delta).
    let before = linear_solution(&source, &[p.i(100), p.zero()], true);
    let source_before = format!("{source:?}");
    let context = RunContext::default();
    for strategy in [
        BoundaryErrorStrategy::ScalarNorm,
        BoundaryErrorStrategy::Automatic,
    ] {
        let options = FlowOptions {
            digits: 24,
            boundary_error_strategy: strategy,
            ..Default::default()
        };
        let problem = problem(&prepared, &source, &options, &context);
        for disks in [1, 4, 16] {
            assert!(matches!(
                problem.scalar_disks(&before, disks),
                Err(Error::Accuracy(_))
            ));
        }
        let mut solution = before.clone();
        let provenance =
            admit_with_gauge(&prepared, &source, &mut solution, &options, &context, None)?.unwrap();
        assert!(
            provenance.contains("fixed balanced weights")
                && provenance.contains("1 disks per segment")
        );
        assert_eq!(solution.coefficients, vec![vec![p.i(100), p.i(1)]]);
        assert!(errors_meet(
            p,
            &solution.coefficients,
            &solution.comparison_errors,
            24
        ));
        assert!(solution.comparison_errors[0][0] > p.real(101) * p.tolerance(49));
        assert!(solution.comparison_errors[0][1] > p.tolerance(48));
        let mut expected = before.clone();
        expected.comparison_errors = solution.comparison_errors.clone();
        expected.checkpoints = solution.checkpoints.clone();
        assert_eq!(format!("{solution:?}"), format!("{expected:?}"));
        assert_eq!(solution.checkpoints.len(), 2);
        for (a, b) in solution.checkpoints.iter().zip(&before.checkpoints) {
            assert_eq!(a.coefficients, b.coefficients);
            assert_eq!(a.point, b.point);
            assert!(errors_meet(p, &a.coefficients, &a.comparison_errors, 24));
        }
        assert_eq!(format!("{source:?}"), source_before);
    }
    // FundamentalMatrix still tries its existing proof before scalar fallback.
    let options = FlowOptions {
        digits: 24,
        boundary_error_strategy: BoundaryErrorStrategy::FundamentalMatrix,
        ..Default::default()
    };
    let mut solution = before;
    assert!(
        admit_with_gauge(&prepared, &source, &mut solution, &options, &context, None)?.is_none()
    );
    assert!(solution.diagnostics.fundamental_boundary_charts > 0);
    Ok(())
}

#[test]
fn balanced_multi_epsilon_complex_weights_cover_source_evidence() -> Result<()> {
    let (prepared, mut source) = coupled_fixture(3, 48, true)?;
    let p = Precision {
        bits: source.accuracy.working_bits(),
    };
    source.coefficients = vec![
        vec![p.zero(), C::new(p.real(0), p.real(1))],
        vec![C::new(p.real(3), p.real(4)), p.zero()],
        vec![p.zero(), p.i(1)],
    ];
    // y0=(-100*x,i), y1=(3+i*(4+30*x),0), y2=(100*i*x,1).
    // These follow directly from J'_k=0, I'_k=100*i*J_k+30*J_(k-1).
    let slopes = vec![
        p.i(-100),
        p.zero(),
        C::new(p.real(0), p.real(30)),
        p.zero(),
        C::new(p.real(0), p.real(100)),
        p.zero(),
    ];
    let before = linear_solution(&source, &slopes, true);
    assert_eq!(
        before.coefficients,
        vec![
            vec![p.i(-100), C::new(p.real(0), p.real(1))],
            vec![C::new(p.real(3), p.real(34)), p.zero()],
            vec![C::new(p.real(0), p.real(100)), p.i(1)]
        ]
    );
    let options = FlowOptions {
        digits: 24,
        boundary_error_strategy: BoundaryErrorStrategy::ScalarNorm,
        ..Default::default()
    };
    let context = RunContext::default();
    let problem = problem(&prepared, &source, &options, &context);
    let balanced = problem.balanced(&before)?;
    // Independent epsilon-major row expectations, including zero entries.
    assert_eq!(balanced.weights[0], p.real(100));
    assert_eq!(balanced.weights[1], p.real(1));
    assert!(balanced.weights[2] > p.real(340) && balanced.weights[2] < p.real(350));
    assert_eq!(balanced.weights[3], p.real(10));
    assert_eq!(balanced.weights[4], p.real(10000));
    assert_eq!(balanced.weights[5], p.real(100));
    for ((weight, base), error) in balanced
        .weights
        .iter()
        .zip(&problem.weights)
        .zip(source.accuracy.comparison_errors().iter().flatten())
    {
        assert!(weight >= base);
        let envelope = weight.clone() * &balanced.relative_error;
        assert!(envelope >= error.clone() + base.clone() * p.tolerance(48));
    }
    for strategy in [
        BoundaryErrorStrategy::ScalarNorm,
        BoundaryErrorStrategy::Automatic,
        BoundaryErrorStrategy::FundamentalMatrix,
    ] {
        let options = FlowOptions {
            boundary_error_strategy: strategy,
            ..options.clone()
        };
        let mut solution = before.clone();
        assert!(
            admit_with_gauge(&prepared, &source, &mut solution, &options, &context, None)?
                .is_some()
        );
        assert_eq!(solution.coefficients, before.coefficients);
        assert_eq!(
            format!("{:?}", solution.segments),
            format!("{:?}", before.segments)
        );
        assert!(errors_meet(
            p,
            &solution.coefficients,
            &solution.comparison_errors,
            24
        ));
        assert_eq!(solution.diagnostics.fundamental_boundary_charts, 0);
        assert!(solution.diagnostics.fundamental_boundary_retry.is_none());
        assert_eq!(
            solution.diagnostics.fundamental_boundary_fallback.is_some(),
            strategy == BoundaryErrorStrategy::FundamentalMatrix
        );
        // Analytic independent perturbation propagation: |delta I0|=sqrt(10001)*delta,
        // |delta I1|=sqrt(31^2+100^2)*delta < 105*delta.
        assert!(solution.comparison_errors[1][0] > p.real(105) * p.tolerance(49));
    }
    Ok(())
}

#[test]
fn balanced_rejects_inaccurate_source_and_floor_without_committing() -> Result<()> {
    for (digits, comparison_digits) in [(24, 24), (24, 90)] {
        let (prepared, mut source) = coupled_fixture(5, digits, false)?;
        let p = Precision {
            bits: source.accuracy.working_bits(),
        };
        // The first case is finite actual source error too large after transfer;
        // the second has tiny comparisons but insufficient verified source digits.
        source.accuracy = BoundaryAccuracy::supplied(
            digits,
            p.bits,
            vec![vec![p.tolerance(comparison_digits); 2]; 5],
            "finite source error/floor rejection control",
        )?;
        source.validate()?;
        let slopes = (0..5)
            .flat_map(|_| [p.i(100), p.zero()])
            .collect::<Vec<_>>();
        let before = linear_solution(&source, &slopes, true);
        let context = RunContext::default();
        for strategy in [
            BoundaryErrorStrategy::ScalarNorm,
            BoundaryErrorStrategy::Automatic,
            BoundaryErrorStrategy::FundamentalMatrix,
        ] {
            let options = FlowOptions {
                digits: 24,
                boundary_error_strategy: strategy,
                ..Default::default()
            };
            let problem = problem(&prepared, &source, &options, &context);
            let original = problem.scalar_disks(&before, 1).err().unwrap();
            let balanced = problem.balanced(&before)?;
            for disks in [1, 4, 16] {
                assert!(matches!(
                    balanced.scalar_disks(&before, disks),
                    Err(Error::Accuracy(_))
                ));
            }
            let mut solution = before.clone();
            let error = admit(&prepared, &source, &mut solution, &options, &context).unwrap_err();
            assert!(matches!(error, Error::Accuracy(_)));
            assert!(error.to_string().starts_with(&original.to_string()));
            assert_eq!(format!("{solution:?}"), format!("{before:?}"));
        }
    }
    Ok(())
}

#[test]
fn balanced_zero_constant_and_malformed_inputs_keep_guards() -> Result<()> {
    let (prepared, mut source) = coupled_fixture(1, 48, false)?;
    let p = Precision {
        bits: source.accuracy.working_bits(),
    };
    source.coefficients[0] = vec![p.zero(); 2];
    let before = linear_solution(&source, &[p.zero(), p.zero()], false);
    let context = RunContext::default();
    let options = FlowOptions {
        digits: 24,
        boundary_error_strategy: BoundaryErrorStrategy::ScalarNorm,
        ..Default::default()
    };
    let problem = problem(&prepared, &source, &options, &context);
    let balanced = problem.balanced(&before)?;
    assert_eq!(balanced.weights, vec![p.real(1); 2]);
    assert_eq!(balanced.relative_error, problem.relative_error);
    let mut solution = before.clone();
    assert!(matches!(
        admit(&prepared, &source, &mut solution, &options, &context),
        Err(Error::Accuracy(_))
    ));
    assert_eq!(format!("{solution:?}"), format!("{before:?}"));
    for variant in 0..5 {
        let mut source = source.clone();
        let mut solution = before.clone();
        match variant {
            0 => source.coefficients.clear(),
            1 => solution.coefficients.clear(),
            2 => solution.coefficients[0].pop().map(|_| ()).unwrap(),
            3 => solution.comparison_errors[0].clear(),
            _ => solution.checkpoints[0].coefficients[0].clear(),
        }
        let saved = format!("{solution:?}");
        assert!(matches!(
            admit(&prepared, &source, &mut solution, &options, &context),
            Err(Error::InvalidInput(_))
        ));
        assert_eq!(format!("{solution:?}"), saved);
    }
    let mut bad = before.clone();
    bad.coefficients[0][0].re = p.real(0) / p.real(0);
    let saved = format!("{bad:?}");
    assert!(matches!(
        admit(&prepared, &source, &mut bad, &options, &context),
        Err(Error::InvalidInput(_))
    ));
    assert_eq!(format!("{bad:?}"), saved);
    context.cancellation.cancel();
    assert!(matches!(balanced.balanced(&before), Err(Error::Cancelled)));
    let mut solution = before.clone();
    assert!(matches!(
        admit(&prepared, &source, &mut solution, &options, &context),
        Err(Error::Cancelled)
    ));
    assert_eq!(format!("{solution:?}"), format!("{before:?}"));
    Ok(())
}
