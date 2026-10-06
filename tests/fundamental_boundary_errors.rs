use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::{
    kinematics::KinematicSystem, physical_transport::PhysicalTransport, transport_cache::*, *,
};

fn problem(
    matrix: Vec<Vec<Atom>>,
    values: Vec<ComplexFloat>,
    digits: u32,
) -> Result<(PhysicalTransport, BoundaryCache)> {
    let x = symbol!("fundamental_public::x");
    let epsilon = symbol!("fundamental_public::epsilon");
    let n = matrix.len();
    let basis = (0..n)
        .map(|i| symbol!("fundamental_public::I").call(i as i64))
        .collect::<Vec<_>>();
    let engine = PhysicalTransport::new(
        KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(x, matrix)]),
        },
        &basis,
        &Atom::one(),
        Prescription::PlusI0,
        "generic regular physical path",
    )?;
    let p = Precision::decimal(90)?;
    let mut cache = BoundaryCache::default();
    cache.insert(CachedBoundary {
        identity: engine.identity().clone(),
        point: CachedPoint::Exact(BTreeMap::from([(x, Atom::new())])),
        kind: PointKind::Physical,
        range: EpsilonRange::new(0, 0)?,
        coefficients: vec![values],
        accuracy: BoundaryAccuracy::supplied(
            digits,
            p.bits,
            vec![vec![p.tolerance(digits + 5); n]],
            "supplied analytic values with declared finite accuracy; storage is not evidence",
        )?,
    })?;
    Ok((engine, cache))
}
fn options() -> FlowOptions {
    FlowOptions {
        boundary_error_strategy: BoundaryErrorStrategy::FundamentalMatrix,
        ..Default::default()
    }
}
fn destination(a: Atom) -> BTreeMap<Symbol, Atom> {
    BTreeMap::from([(symbol!("fundamental_public::x"), a)])
}
fn policy() -> impl TransportCost {
    ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    }
}

#[test]
fn nilpotent_physical_samples_preserve_source_floor_and_restart() -> Result<()> {
    let p = Precision::decimal(90)?;
    for rate in [100_i64, 2500, 5800] {
        let matrix = vec![
            vec![Atom::num(rate), Atom::num(rate)],
            vec![Atom::num(-rate), Atom::num(-rate)],
        ];
        let (engine, mut cache) = problem(matrix, vec![p.i(1), p.zero()], 30)?;
        let mut scalar_cache = cache.clone();
        let rejected = engine.evaluate_to(
            &mut scalar_cache,
            &destination(Atom::one()),
            EpsilonRange::new(0, 0)?,
            &FlowOptions {
                boundary_error_strategy: BoundaryErrorStrategy::ScalarNorm,
                ..Default::default()
            },
            &RunContext::default(),
            &policy(),
        );
        assert!(matches!(rejected, Err(Error::Accuracy(_))));
        assert_eq!(scalar_cache.len(), 1);
        let result = engine.evaluate_to(
            &mut cache,
            &destination(Atom::one()),
            EpsilonRange::new(0, 0)?,
            &options(),
            &RunContext::default(),
            &policy(),
        )?;
        let diagnostics = &result.transport.as_ref().unwrap().diagnostics;
        assert!(diagnostics.fundamental_boundary_charts > 0);
        assert!(diagnostics.fundamental_boundary_fallback.is_none());
        assert_eq!(result.boundary.accuracy.input_verified_digits(), 30);
        assert!((20..=30).contains(&result.boundary.accuracy.verified_digits()));
        for (actual, expected) in result.boundary.coefficients[0]
            .iter()
            .zip([p.i(rate + 1), p.i(-rate)])
        {
            assert!(p.close(actual, &expected, 45));
        }
        assert!(
            result.boundary.accuracy.comparison_errors()[0]
                .iter()
                .all(|error| *error > p.tolerance(30))
        );
        let dir =
            std::env::temp_dir().join(format!("fundamental-restart-{}-{rate}", std::process::id()));
        cache.save(&dir)?;
        let mut loaded = BoundaryCache::load(&dir)?;
        std::fs::remove_dir_all(&dir)?;
        // Verification policy does not create a new mathematical cache identity.
        let hit = engine.evaluate_to(
            &mut loaded,
            &destination(Atom::one()),
            EpsilonRange::new(0, 0)?,
            &FlowOptions::default(),
            &RunContext::default(),
            &policy(),
        )?;
        assert!(hit.transport.is_none());
        assert_eq!(hit.boundary.coefficients, result.boundary.coefficients);
        let nearby = engine.evaluate_to(
            &mut loaded,
            &destination(Atom::num((1001, 1000))),
            EpsilonRange::new(0, 0)?,
            &options(),
            &RunContext::default(),
            &policy(),
        )?;
        assert_eq!(
            nearby.starting_point.rounded_coordinates_as_exact()?,
            destination(Atom::one())
        );
        assert!(
            nearby
                .transport
                .as_ref()
                .unwrap()
                .diagnostics
                .fundamental_boundary_charts
                > 0
        );
        assert!((20..=30).contains(&nearby.boundary.accuracy.input_verified_digits()));
    }
    Ok(())
}

#[test]
fn exponential_instability_rejects_weak_evidence_with_unchanged_cache() -> Result<()> {
    let p = Precision::decimal(90)?;
    let matrix = vec![
        vec![Atom::num(20), Atom::new()],
        vec![Atom::new(), Atom::new()],
    ];
    let (strong, mut strong_cache) = problem(matrix.clone(), vec![p.zero(), p.i(1)], 60)?;
    let accepted = strong.evaluate_to(
        &mut strong_cache,
        &destination(Atom::one()),
        EpsilonRange::new(0, 0)?,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert!(
        accepted
            .transport
            .as_ref()
            .unwrap()
            .diagnostics
            .fundamental_boundary_charts
            > 0
    );
    assert!(accepted.boundary.accuracy.comparison_errors()[0][0] > p.tolerance(53));
    let (weak, mut weak_cache) = problem(matrix, vec![p.zero(), p.i(1)], 25)?;
    let result = weak.evaluate_to(
        &mut weak_cache,
        &destination(Atom::one()),
        EpsilonRange::new(0, 0)?,
        &options(),
        &RunContext::default(),
        &policy(),
    );
    assert!(matches!(result, Err(Error::Accuracy(_))));
    assert_eq!(weak_cache.len(), 1);
    Ok(())
}

#[test]
fn multiple_rational_charts_keep_verified_intermediates_and_source_cap() -> Result<()> {
    let p = Precision::decimal(90)?;
    let x = Atom::var(symbol!("fundamental_public::x"));
    let matrix = vec![vec![Atom::one() / (&x + 1)]];
    let (engine, mut cache) = problem(matrix, vec![p.i(1)], 40)?;
    let result = engine.evaluate_to(
        &mut cache,
        &destination(Atom::num(8)),
        EpsilonRange::new(0, 0)?,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    let diagnostics = &result.transport.as_ref().unwrap().diagnostics;
    assert!(diagnostics.steps > 1);
    assert_eq!(diagnostics.fundamental_boundary_charts, diagnostics.steps);
    assert!(result.inserted_points > 1);
    assert!(p.close(&result.boundary.coefficients[0][0], &p.i(9), 30));
    assert_eq!(result.boundary.accuracy.input_verified_digits(), 40);
    for entry in cache.entries() {
        assert!(entry.accuracy.verified_digits() <= 40);
    }
    Ok(())
}

#[test]
fn unsupported_hierarchy_and_dense_limits_retain_scalar_fallback() -> Result<()> {
    let p = Precision::decimal(90)?;
    let n = 17;
    let (engine, mut cache) = problem(vec![vec![Atom::new(); n]; n], vec![p.i(1); n], 40)?;
    let result = engine.evaluate_to(
        &mut cache,
        &destination(Atom::one()),
        EpsilonRange::new(0, 0)?,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    let diagnostics = &result.transport.as_ref().unwrap().diagnostics;
    assert_eq!(diagnostics.fundamental_boundary_charts, 0);
    assert!(
        diagnostics
            .fundamental_boundary_fallback
            .as_ref()
            .unwrap()
            .contains("dimension")
    );
    let (engine, original) = problem(vec![vec![Atom::new()]], vec![p.i(1)], 40)?;
    let mut entry = original.entries()[0].clone();
    entry.range = EpsilonRange::new(0, 1)?;
    entry.coefficients.push(vec![p.zero()]);
    entry.accuracy = BoundaryAccuracy::supplied(
        40,
        p.bits,
        vec![vec![p.tolerance(45)]; 2],
        "two supplied coefficients",
    )?;
    let mut cache = BoundaryCache::default();
    cache.insert(entry)?;
    let result = engine.evaluate_to(
        &mut cache,
        &destination(Atom::one()),
        EpsilonRange::new(0, 1)?,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert!(
        result
            .transport
            .as_ref()
            .unwrap()
            .diagnostics
            .fundamental_boundary_fallback
            .as_ref()
            .unwrap()
            .contains("epsilon hierarchies")
    );
    Ok(())
}

#[test]
fn mapped_rational_central_candidates_use_the_same_physical_proof() -> Result<()> {
    let p = Precision::decimal(90)?;
    let x = Atom::var(symbol!("fundamental_public::x"));
    let (engine, mut cache) = problem(vec![vec![Atom::one() / (&x + 100)]], vec![p.i(1)], 40)?;
    let opts = FlowOptions {
        local_coordinate: LocalCoordinate::BalancedMobius,
        pade: Some(PadeOptions {
            degree: 8,
            ..Default::default()
        }),
        ..options()
    };
    let result = engine.evaluate_to(
        &mut cache,
        &destination(Atom::one()),
        EpsilonRange::new(0, 0)?,
        &opts,
        &RunContext::default(),
        &policy(),
    )?;
    let solution = result.transport.as_ref().unwrap();
    assert!(solution.diagnostics.fundamental_boundary_charts > 0);
    assert!(solution.diagnostics.pade_steps > 0);
    assert!(solution.segments.iter().any(|segment| !matches!(
        segment.coordinate,
        symbolica_amflow::local_coordinates::TaylorCoordinate::Identity
    )));
    assert!(p.close(
        &result.boundary.coefficients[0][0],
        &p.rational(&Rational::from((101, 100))),
        30
    ));
    Ok(())
}

#[test]
fn cancellation_during_optional_proof_is_typed_and_does_not_commit_cache() -> Result<()> {
    let p = Precision::decimal(90)?;
    let (engine, mut cache) = problem(vec![vec![Atom::one()]], vec![p.i(1)], 40)?;
    let cancellation = CancellationToken::default();
    let trigger = cancellation.clone();
    let context = RunContext {
        cancellation,
        progress: Some(std::sync::Arc::new(move |event| {
            if matches!(event,Progress::Stage { name } if name.contains("signed fundamental")) {
                trigger.cancel();
            }
        })),
    };
    let result = engine.evaluate_to(
        &mut cache,
        &destination(Atom::one()),
        EpsilonRange::new(0, 0)?,
        &options(),
        &context,
        &policy(),
    );
    assert!(matches!(result, Err(Error::Cancelled)));
    assert_eq!(cache.len(), 1);
    Ok(())
}

#[test]
fn automatic_scalar_success_and_explicit_policy_keep_identical_results() -> Result<()> {
    let p = Precision::decimal(90)?;
    let (engine, mut first) = problem(vec![vec![Atom::one()]], vec![p.i(1)], 40)?;
    let mut second = first.clone();
    let explicit = FlowOptions {
        boundary_error_strategy: BoundaryErrorStrategy::ScalarNorm,
        ..Default::default()
    };
    let a = engine.evaluate_to(
        &mut first,
        &destination(Atom::one()),
        EpsilonRange::new(0, 0)?,
        &FlowOptions::default(),
        &RunContext::default(),
        &policy(),
    )?;
    let b = engine.evaluate_to(
        &mut second,
        &destination(Atom::one()),
        EpsilonRange::new(0, 0)?,
        &explicit,
        &RunContext::default(),
        &policy(),
    )?;
    assert_eq!(a.boundary.coefficients, b.boundary.coefficients);
    assert_eq!(
        a.boundary.accuracy.comparison_errors(),
        b.boundary.accuracy.comparison_errors()
    );
    assert_eq!(
        a.transport
            .as_ref()
            .unwrap()
            .diagnostics
            .fundamental_boundary_charts,
        0
    );
    Ok(())
}

#[test]
fn automatic_retry_reuses_the_central_trajectory_and_verified_errors() -> Result<()> {
    use std::sync::{Arc, Mutex};
    let p = Precision::decimal(90)?;
    let x = Atom::var(symbol!("fundamental_public::x"));
    let entry = Atom::num(5800) / (&x + 1);
    let (engine, original) = problem(
        vec![
            vec![entry.clone(), entry.clone()],
            vec![-entry.clone(), -entry],
        ],
        vec![p.i(1), p.zero()],
        30,
    )?;
    let mut runs = Vec::new();
    for selected in [FlowOptions::default(), options()] {
        let events = Arc::new(Mutex::new(Vec::new()));
        let capture = events.clone();
        let context = RunContext {
            progress: Some(Arc::new(move |event| {
                if let Progress::Step { index } = event {
                    capture.lock().unwrap().push(index);
                }
            })),
            ..Default::default()
        };
        let mut cache = original.clone();
        let result = engine.evaluate_to(
            &mut cache,
            &destination(Atom::num(8)),
            EpsilonRange::new(0, 0)?,
            &selected,
            &context,
            &policy(),
        )?;
        let progress = events.lock().unwrap().clone();
        runs.push((result, progress, cache));
    }
    let (automatic, steps, cache) = &runs[0];
    let (forced, forced_steps, _) = &runs[1];
    assert!(!steps.is_empty());
    assert_eq!(
        steps, forced_steps,
        "an uncertainty retry must not solve again"
    );
    assert_eq!(
        automatic.boundary.coefficients,
        forced.boundary.coefficients
    );
    assert_eq!(
        automatic.boundary.accuracy.comparison_errors(),
        forced.boundary.accuracy.comparison_errors()
    );
    let automatic_solution = automatic.transport.as_ref().unwrap();
    let forced_solution = forced.transport.as_ref().unwrap();
    assert!(automatic_solution.diagnostics.steps > 1);
    assert!(automatic_solution.checkpoints.len() > 1);
    let shift = p.scale(&p.log(&p.i(9)), 5800, 1);
    assert!(p.close(
        &automatic.boundary.coefficients[0][0],
        &p.add(&p.i(1), &shift),
        40
    ));
    assert!(p.close(&automatic.boundary.coefficients[0][1], &p.neg(&shift), 40));
    assert!(
        automatic_solution
            .diagnostics
            .fundamental_boundary_retry
            .is_some()
    );
    assert!(
        forced_solution
            .diagnostics
            .fundamental_boundary_retry
            .is_none()
    );
    assert_eq!(
        automatic_solution.checkpoints.len(),
        forced_solution.checkpoints.len()
    );
    for (a, b) in automatic_solution
        .checkpoints
        .iter()
        .zip(&forced_solution.checkpoints)
    {
        assert_eq!(a.segment, b.segment);
        assert_eq!(a.coefficients, b.coefficients);
        assert_eq!(a.comparison_errors, b.comparison_errors);
    }
    assert!(
        cache
            .entries()
            .iter()
            .all(|entry| entry.accuracy.verified_digits() <= 30)
    );
    Ok(())
}

#[test]
fn automatic_handles_nonfinite_scalar_amplification() -> Result<()> {
    let p = Precision::decimal(90)?;
    // exp(2*rate) exceeds the native MPFR exponent range, while the signed
    // fundamental matrix and its source-error bounds remain finite.
    let rate = 1_000_000_000_i64;
    let (engine, original) = problem(
        vec![
            vec![Atom::num(rate), Atom::num(rate)],
            vec![Atom::num(-rate), Atom::num(-rate)],
        ],
        vec![p.i(1), p.zero()],
        30,
    )?;
    let mut scalar = original.clone();
    let rejected = engine.evaluate_to(
        &mut scalar,
        &destination(Atom::one()),
        EpsilonRange::new(0, 0)?,
        &FlowOptions {
            boundary_error_strategy: BoundaryErrorStrategy::ScalarNorm,
            ..Default::default()
        },
        &RunContext::default(),
        &policy(),
    );
    assert!(
        matches!(&rejected, Err(Error::Accuracy(message)) if message.contains("amplification exceeds numerical range")),
        "{:?}",
        rejected.as_ref().err()
    );
    assert_eq!(scalar.len(), 1);
    let mut cache = original;
    let result = engine.evaluate_to(
        &mut cache,
        &destination(Atom::one()),
        EpsilonRange::new(0, 0)?,
        &FlowOptions::default(),
        &RunContext::default(),
        &policy(),
    )?;
    assert!(
        result
            .transport
            .as_ref()
            .unwrap()
            .diagnostics
            .fundamental_boundary_retry
            .as_ref()
            .unwrap()
            .contains("amplification exceeds numerical range")
    );
    assert_eq!(result.boundary.accuracy.input_verified_digits(), 30);
    assert!((20..=30).contains(&result.boundary.accuracy.verified_digits()));
    for (actual, expected) in result.boundary.coefficients[0]
        .iter()
        .zip([p.i(rate + 1), p.i(-rate)])
    {
        assert!(p.close(actual, &expected, 45));
    }
    Ok(())
}

#[test]
fn automatic_optional_limit_and_real_instability_preserve_failure_and_cache() -> Result<()> {
    let p = Precision::decimal(90)?;
    for n in [2, 17] {
        let mut matrix = vec![vec![Atom::new(); n]; n];
        let mut values = vec![p.zero(); n];
        if n == 2 {
            matrix[0][0] = Atom::num(20);
            values[1] = p.i(1);
        } else {
            for row in &mut matrix[..2] {
                row[0] = Atom::num(5800);
                row[1] = Atom::num(5800);
            }
            matrix[1][0] = Atom::num(-5800);
            matrix[1][1] = Atom::num(-5800);
            values[0] = p.i(1);
        }
        let (engine, mut cache) = problem(matrix, values, 25)?;
        let result = engine.evaluate_to(
            &mut cache,
            &destination(Atom::one()),
            EpsilonRange::new(0, 0)?,
            &FlowOptions::default(),
            &RunContext::default(),
            &policy(),
        );
        assert!(
            matches!(&result, Err(Error::Accuracy(_))),
            "{:?}",
            result.as_ref().err()
        );
        if n == 17 {
            assert!(result.err().unwrap().to_string().contains("dimension"));
        }
        assert_eq!(cache.len(), 1);
    }
    Ok(())
}

#[test]
fn automatic_cancellation_during_retry_keeps_cache_unchanged() -> Result<()> {
    let p = Precision::decimal(90)?;
    let (engine, mut cache) = problem(
        vec![
            vec![Atom::num(5800), Atom::num(5800)],
            vec![Atom::num(-5800), Atom::num(-5800)],
        ],
        vec![p.i(1), p.zero()],
        30,
    )?;
    let cancellation = CancellationToken::default();
    let trigger = cancellation.clone();
    let context = RunContext {
        cancellation,
        progress: Some(std::sync::Arc::new(move |event| {
            if matches!(event,Progress::Stage { name } if name.contains("signed fundamental")) {
                trigger.cancel();
            }
        })),
    };
    let result = engine.evaluate_to(
        &mut cache,
        &destination(Atom::one()),
        EpsilonRange::new(0, 0)?,
        &FlowOptions::default(),
        &context,
        &policy(),
    );
    assert!(matches!(result, Err(Error::Cancelled)));
    assert_eq!(cache.len(), 1);
    Ok(())
}

#[test]
fn automatic_does_not_retry_central_domain_or_step_failures() -> Result<()> {
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };
    let p = Precision::decimal(90)?;
    let x = Atom::var(symbol!("fundamental_public::x"));
    for domain in [true, false] {
        let matrix = if domain {
            Atom::one() / (&x - Atom::num((1, 2)))
        } else {
            Atom::one() / (&x + 1)
        };
        let (engine, mut cache) = problem(vec![vec![matrix]], vec![p.i(1)], 40)?;
        let attempted = Arc::new(AtomicUsize::new(0));
        let counter = attempted.clone();
        let context = RunContext {
            progress: Some(Arc::new(move |event| {
                if matches!(event,Progress::Stage { name } if name.contains("signed fundamental")) {
                    counter.fetch_add(1, Ordering::SeqCst);
                }
            })),
            ..Default::default()
        };
        let selected = FlowOptions {
            max_steps: 1,
            ..Default::default()
        };
        let result = engine.evaluate_to(
            &mut cache,
            &destination(Atom::num(if domain { 1 } else { 8 })),
            EpsilonRange::new(0, 0)?,
            &selected,
            &context,
            &policy(),
        );
        if domain {
            assert!(
                matches!(
                    result,
                    Err(Error::Unsupported(_) | Error::IncompleteReduction(_))
                ),
                "{:?}",
                result.as_ref().err()
            );
        } else {
            assert!(
                matches!(&result, Err(Error::Limit(_))),
                "{:?}",
                result.as_ref().err()
            );
        }
        assert_eq!(attempted.load(Ordering::SeqCst), 0);
        assert_eq!(cache.len(), 1);
    }
    Ok(())
}
