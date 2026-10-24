use std::collections::BTreeMap;
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use std::time::Duration;
use symbolica::prelude::*;
use symbolica_amflow::{cuts::*, phase_space::PreparedMasslessPhaseSpace, reduction::Reduction, *};

struct UnitCutBackend;

impl ReductionBackend for UnitCutBackend {
    fn identity(&self) -> String {
        "resource-test-unit-cut-residual".into()
    }

    fn reduce(&self, _: &IntegralFamily, _: &[Integral], _: &RunContext) -> Result<Reduction> {
        panic!("the ordinary reduction measure must not be used")
    }

    fn reduce_cut(&self, _: &CutFamily, targets: &[Integral], _: &RunContext) -> Result<Reduction> {
        assert!(targets.iter().all(|target| target.0 == [1, 1]));
        Ok(Reduction {
            residuals: targets.to_vec(),
            ..Default::default()
        })
    }
}

fn family() -> Result<CutFamily> {
    let gram = vec![vec![Atom::num(25)]];
    let ordinary = IntegralFamily {
        name: "resource_scoped_unit_cut".into(),
        loops: vec!["k".into()],
        external: vec!["p".into()],
        external_gram: gram.clone(),
        propagators: vec![
            Propagator::quadratic(&[1], &[0], Atom::zero(), &gram)?,
            Propagator::quadratic(&[-1], &[1], Atom::zero(), &gram)?,
        ],
        physical_propagators: 2,
        epsilon: symbol!("cut_resource_test::eps"),
        dimension: 4,
    };
    CutFamily::new(
        ordinary,
        CutMetadata::new(
            [(1, 0), (-1, 1)].into_iter().enumerate().map(
                |(propagator, (loop_coefficient, external_coefficient))| CutLine {
                    propagator,
                    definition: CutDefinition::PositiveEnergy {
                        momentum: MomentumRouting {
                            loops: vec![Rational::from(loop_coefficient)],
                            external: vec![Rational::from(external_coefficient)],
                        },
                    },
                },
            ),
            vec![LoopPrescription::Insensitive],
        )?,
    )
}

fn assert_expired<T>(result: Result<T>) {
    assert!(matches!(
        result,
        Err(Error::Resource(ResourceLimitError::Elapsed { .. }))
    ));
}

fn worker_expiry_preserves_prepared_values(
    evaluate: impl Fn(&RunContext) -> Result<Vec<Vec<ComplexFloat>>>,
) -> Result<()> {
    let baseline = evaluate(&RunContext::default())?;
    let expired = Arc::new(AtomicBool::new(false));
    let observed = Arc::clone(&expired);
    let budget = ResourceBudget::with_memory_monitor(
        ResourceLimits {
            max_process_rss_bytes: Some(1_000_000),
            ..Default::default()
        },
        move || {
            Ok(if observed.load(Ordering::SeqCst) {
                2_000_000
            } else {
                0
            })
        },
    );
    let caller_thread = std::thread::current().id();
    let saw_inherited_error = Arc::new(AtomicBool::new(false));
    let saw_worker = Arc::new(AtomicBool::new(false));
    let expires_in_callback = Arc::clone(&expired);
    let worker_error = Arc::clone(&saw_inherited_error);
    let worker_thread = Arc::clone(&saw_worker);
    let context = RunContext {
        progress: Some(Arc::new(move |event| {
            if matches!(event, Progress::Sample { .. }) {
                worker_thread.fetch_or(
                    std::thread::current().id() != caller_thread,
                    Ordering::SeqCst,
                );
                expires_in_callback.store(true, Ordering::SeqCst);
                worker_error.fetch_or(
                    matches!(
                        symbolica::utils::resources::checkpoint(),
                        Err(ResourceLimitError::Memory {
                            rss_bytes: 2_000_000,
                            additional_bytes: 0,
                            limit_bytes: 1_000_000,
                        })
                    ),
                    Ordering::SeqCst,
                );
            }
        })),
        ..Default::default()
    };
    assert!(matches!(
        budget.run(|| evaluate(&context)),
        Err(Error::Resource(ResourceLimitError::Memory {
            rss_bytes: 2_000_000,
            additional_bytes: 0,
            limit_bytes: 1_000_000,
        }))
    ));
    assert!(saw_worker.load(Ordering::SeqCst));
    assert!(saw_inherited_error.load(Ordering::SeqCst));
    expired.store(false, Ordering::SeqCst);
    assert_eq!(evaluate(&RunContext::default())?, baseline);
    Ok(())
}

#[test]
fn prepared_phase_space_scopes_options_and_inherits_worker_budgets() -> Result<()> {
    let family = family()?;
    let channel = FutureTimelikeChannel {
        external: vec![Rational::one()],
    };
    let target = [Integral(vec![1, 1])];
    let samples = [Rational::from((1, 13)), Rational::from((1, 17))];
    let options = FlowOptions {
        workers: 2,
        ..Default::default()
    };
    let limited = FlowOptions {
        resource_limits: ResourceLimits {
            max_elapsed: Some(Duration::ZERO),
            ..Default::default()
        },
        ..options.clone()
    };
    let context = RunContext::default();
    let point = KinematicPoint::default();
    let prepared = PreparedMasslessPhaseSpace::new(
        &family,
        &channel,
        &target,
        &point,
        &UnitCutBackend,
        &options,
        &context,
    )?;
    assert_expired(PreparedMasslessPhaseSpace::new(
        &family,
        &channel,
        &target,
        &point,
        &UnitCutBackend,
        &limited,
        &context,
    ));
    assert_expired(prepared.evaluate(&samples[0], &limited, &context));
    assert_expired(prepared.evaluate_samples(&samples, &limited, &context));
    assert_expired(prepared.solve(0, &limited, &context));
    worker_expiry_preserves_prepared_values(|context| {
        prepared.evaluate_samples(&samples, &options, context)
    })
}

#[test]
fn prepared_cut_projections_scope_options_and_inherit_worker_budgets() -> Result<()> {
    let family = family()?;
    let channel = FutureTimelikeChannel {
        external: vec![Rational::one()],
    };
    let projections = [BTreeMap::from([(Integral(vec![1, 1]), Atom::one())])];
    let samples = [Rational::from((1, 13)), Rational::from((1, 17))];
    let options = FlowOptions {
        workers: 2,
        ..Default::default()
    };
    let limited = FlowOptions {
        resource_limits: ResourceLimits {
            max_elapsed: Some(Duration::ZERO),
            ..Default::default()
        },
        ..options.clone()
    };
    let context = RunContext::default();
    let point = KinematicPoint::default();
    let prepared = PreparedCutProjections::new(
        &family,
        &channel,
        &projections,
        &point,
        &UnitCutBackend,
        &options,
        &context,
    )?;
    assert_expired(PreparedCutProjections::new(
        &family,
        &channel,
        &projections,
        &point,
        &UnitCutBackend,
        &limited,
        &context,
    ));
    assert_expired(prepared.evaluate(&samples[0], &limited, &context));
    assert_expired(prepared.evaluate_samples(&samples, &limited, &context));
    assert_expired(prepared.solve(0, &limited, &context));
    worker_expiry_preserves_prepared_values(|context| {
        prepared.evaluate_samples(&samples, &options, context)
    })
}
