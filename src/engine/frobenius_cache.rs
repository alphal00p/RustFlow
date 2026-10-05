//! Lazy exact endpoint/infinity preparation shared by the flow's numerical samples.
use crate::frobenius::PreparedFrobenius;
use crate::{DifferentialSystem, Error, Result, RunContext};
use std::sync::{Arc, Mutex, TryLockError};
use std::time::Duration;
use symbolica::prelude::*;

#[derive(Default)]
pub(super) struct FrobeniusPreparations {
    zero: Mutex<Option<CachedPreparation>>,
    infinity: Mutex<Option<CachedPreparation>>,
}

struct CachedPreparation {
    // Deliberately retain the original expressions, before any cancellation.
    // PreparedFlow.system is public and may be replaced between evaluations.
    source: DifferentialSystem,
    prepared: Arc<PreparedFrobenius>,
}

impl FrobeniusPreparations {
    pub(super) fn get(
        &self,
        system: &DifferentialSystem,
        at_infinity: bool,
        context: &RunContext,
    ) -> Result<Arc<PreparedFrobenius>> {
        let slot = if at_infinity {
            &self.infinity
        } else {
            &self.zero
        };
        // Progress callbacks must not reenter the same calculation while this
        // preparation is in progress. Cancellation from callbacks is supported.
        // Only one sample constructs each exact preparation. Waiting callers
        // still observe their own cancellation token; cancellation of the
        // producer returns Err without installing a partial cache entry.
        let mut entry = loop {
            context.cancellation.check()?;
            match slot.try_lock() {
                Ok(entry) => break entry,
                Err(TryLockError::WouldBlock) => std::thread::sleep(Duration::from_millis(10)),
                Err(TryLockError::Poisoned(_)) => {
                    return Err(Error::Numerical(
                        "Frobenius preparation cache lock poisoned".into(),
                    ));
                }
            }
        };
        if let Some(cached) = entry.as_ref()
            && cached.source.variable == system.variable
            && cached.source.matrix == system.matrix
        {
            return Ok(cached.prepared.clone());
        }
        let prepared = if at_infinity {
            system
                .invert_variable(symbol!("symbolica_amflow::z"))
                .prepare_frobenius(context)?
        } else {
            system.prepare_frobenius(context)?
        };
        context.cancellation.check()?;
        let prepared = Arc::new(prepared);
        *entry = Some(CachedPreparation {
            source: system.clone(),
            prepared: prepared.clone(),
        });
        Ok(prepared)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CancellationToken, Precision, Progress};
    use std::sync::atomic::{AtomicUsize, Ordering};

    fn system() -> DifferentialSystem {
        DifferentialSystem {
            variable: symbol!("cached_frobenius::x"),
            matrix: vec![vec![parse!(
                "cached_frobenius::epsilon/cached_frobenius::x"
            )]],
        }
    }

    #[test]
    fn concurrent_samples_and_precision_refinements_share_exact_work() {
        let cache = FrobeniusPreparations::default();
        let system = system();
        let starts = Arc::new(AtomicUsize::new(0));
        let count = starts.clone();
        let context = RunContext {
            progress: Some(Arc::new(move |progress| {
                if matches!(progress, Progress::Stage { name } if name == "normalizing the exact Frobenius system")
                {
                    count.fetch_add(1, Ordering::Relaxed);
                }
            })),
            ..Default::default()
        };
        std::thread::scope(|scope| {
            let jobs = (0..4)
                .map(|sample| {
                    let cache = &cache;
                    let system = &system;
                    let context = &context;
                    scope.spawn(move || {
                        let prepared = cache.get(system, false, context).unwrap();
                        let p = Precision::decimal(40 + 10 * sample).unwrap();
                        let epsilon = Rational::from((1, 31 + sample as i64));
                        let params = ahash::HashMap::from_iter([(
                            parse!("cached_frobenius::epsilon"),
                            p.rational(&epsilon),
                        )]);
                        let basis = prepared
                            .evaluate(p, &params, 8 + sample as usize, context)
                            .unwrap();
                        assert_eq!(
                            basis.columns[0].exponent,
                            parse!("cached_frobenius::epsilon")
                        );
                        let x = p.rational(&Rational::from((1, 2)));
                        let actual = basis.evaluate(&x, &params).unwrap();
                        assert!(p.close(&actual[0][0], &p.pow(&x, &p.rational(&epsilon)), 30));
                        prepared
                    })
                })
                .collect::<Vec<_>>();
            let preparations = jobs
                .into_iter()
                .map(|job| job.join().unwrap())
                .collect::<Vec<_>>();
            assert!(
                preparations
                    .iter()
                    .all(|prepared| Arc::ptr_eq(prepared, &preparations[0]))
            );
        });
        assert_eq!(starts.load(Ordering::Relaxed), 1);
        let infinity = cache.get(&system, true, &context).unwrap();
        assert_eq!(
            infinity.exponents().next().unwrap(),
            &-parse!("cached_frobenius::epsilon")
        );
        assert!(Arc::ptr_eq(
            &infinity,
            &cache.get(&system, true, &context).unwrap()
        ));
        assert_eq!(starts.load(Ordering::Relaxed), 2);
    }

    #[test]
    fn changed_public_matrix_variable_and_uncancelled_domain_invalidate_preparation() {
        let cache = FrobeniusPreparations::default();
        let context = RunContext::default();
        let mut system = system();
        let original = cache.get(&system, false, &context).unwrap();
        system.matrix[0][0] = parse!("2*cached_frobenius::epsilon/cached_frobenius::x");
        let changed = cache.get(&system, false, &context).unwrap();
        assert!(!Arc::ptr_eq(&original, &changed));
        assert_eq!(
            changed.exponents().next().unwrap(),
            &parse!("2*cached_frobenius::epsilon")
        );
        system.variable = symbol!("cached_frobenius::other");
        let variable = cache.get(&system, false, &context).unwrap();
        assert!(!Arc::ptr_eq(&changed, &variable));
        assert_eq!(variable.exponents().next().unwrap(), &Atom::new());
        system.matrix[0][0] =
            parse!("(cached_frobenius::epsilon^2-1)/(cached_frobenius::epsilon-1)");
        let guarded = cache.get(&system, false, &context).unwrap();
        // Inversion must collect the original excluded epsilon=1 domain before
        // normalizing away (epsilon^2-1)/(epsilon-1).
        let p = Precision::decimal(40).unwrap();
        let parameters = ahash::HashMap::from_iter([(parse!("cached_frobenius::epsilon"), p.i(1))]);
        let logarithmic = DifferentialSystem {
            variable: symbol!("cached_frobenius::x"),
            matrix: vec![vec![parse!(
                "(cached_frobenius::epsilon^2-1)/((cached_frobenius::epsilon-1)*cached_frobenius::x)"
            )]],
        };
        let infinity = cache.get(&logarithmic, true, &context).unwrap();
        assert!(
            matches!(infinity.evaluate(p, &parameters, 8, &context), Err(Error::Numerical(message)) if message.contains("denominator condition"))
        );
        system.matrix[0][0] = parse!("cached_frobenius::epsilon+1");
        let cancelled = cache.get(&system, false, &context).unwrap();
        assert!(!Arc::ptr_eq(&guarded, &cancelled));
        assert!(!guarded.nonzero_conditions().is_empty());
        assert!(cancelled.nonzero_conditions().is_empty());
    }

    #[test]
    fn cancelled_preparation_is_retryable_and_waiting_caller_can_cancel() {
        let cache = FrobeniusPreparations::default();
        let system = system();
        let token = CancellationToken::default();
        let cancel = token.clone();
        let context = RunContext {
            cancellation: token,
            progress: Some(Arc::new(move |progress| {
                if matches!(progress, Progress::Stage { name } if name.starts_with("constructing exact Frobenius eigenspace"))
                {
                    cancel.cancel();
                }
            })),
        };
        assert!(matches!(
            cache.get(&system, false, &context),
            Err(Error::Cancelled)
        ));
        assert!(cache.zero.lock().unwrap().is_none());
        assert!(cache.get(&system, false, &RunContext::default()).is_ok());
        let held = cache.zero.lock().unwrap();
        std::thread::scope(|scope| {
            let context = RunContext::default();
            let token = context.cancellation.clone();
            let cache = &cache;
            let system = &system;
            let worker = scope.spawn(move || cache.get(system, false, &context));
            token.cancel();
            assert!(matches!(worker.join().unwrap(), Err(Error::Cancelled)));
        });
        drop(held);
        assert!(cache.get(&system, false, &RunContext::default()).is_ok());
    }
}
