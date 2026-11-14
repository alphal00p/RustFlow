//! Primary reduction backend with a fallback for small boundary children.
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

use super::{Reduction, ReductionBackend};
use crate::{Error, Integral, IntegralFamily, Result, RunContext};
use symbolica::prelude::*;

/// Which backend answers every request of one family.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Choice {
    Primary,
    Fallback,
}

/// Reduce with `primary`; when it fails (other than by cancellation) on the
/// FIRST request of a family with fewer than `fallback_below_loops` loops,
/// reduce that family with `fallback` instead.
///
/// The choice is made once per family and kept: a differential closure calls
/// the backend several times for one family, and rules from two backends are
/// in different master bases, so they must never be mixed. A later failure of
/// the chosen backend for that family is returned as an error instead of
/// switching backends.
///
/// Typical use: a LiteRed lane for every family, with RustRed's native search
/// as a fallback for the lower-loop boundary children the lane rejects, but
/// never for the top family itself (`fallback_below_loops` = its loop count).
pub struct FallbackBackend<P, F> {
    pub primary: P,
    pub fallback: Option<F>,
    pub fallback_below_loops: usize,
    choices: Mutex<BTreeMap<String, Arc<Mutex<Option<Choice>>>>>,
}

impl<P: ReductionBackend, F: ReductionBackend> FallbackBackend<P, F> {
    pub fn new(primary: P, fallback: Option<F>, fallback_below_loops: usize) -> Self {
        Self {
            primary,
            fallback,
            fallback_below_loops,
            choices: Mutex::new(BTreeMap::new()),
        }
    }

    fn slot(&self, family: &IntegralFamily) -> Result<Arc<Mutex<Option<Choice>>>> {
        // The complete family declaration identifies the family; epsilon
        // samples of one family share its choice (and its master basis).
        let key = blake3::hash(format!("{family:?}").as_bytes())
            .to_hex()
            .to_string();
        let mut choices = self
            .choices
            .lock()
            .map_err(|_| Error::Reduction("fallback backend choice map poisoned".into()))?;
        Ok(Arc::clone(choices.entry(key).or_default()))
    }

    fn run(
        &self,
        family: &IntegralFamily,
        reduce: impl Fn(&dyn ReductionBackend) -> Result<Reduction>,
    ) -> Result<Reduction> {
        let slot = self.slot(family)?;
        // Holding the family slot serializes only the first, deciding request.
        let mut choice = slot
            .lock()
            .map_err(|_| Error::Reduction("fallback backend choice poisoned".into()))?;
        match *choice {
            Some(Choice::Primary) => {
                drop(choice);
                return reduce(&self.primary).map_err(|error| match error {
                    Error::Cancelled => Error::Cancelled,
                    error => Error::Reduction(format!(
                        "{} failed for family {} after it was chosen for that family \
                         (backends are not switched within a family): {error}",
                        self.primary.identity(),
                        family.name
                    )),
                });
            }
            Some(Choice::Fallback) => {
                drop(choice);
                // Decided only when a fallback exists.
                let fallback = self
                    .fallback
                    .as_ref()
                    .ok_or_else(|| Error::Reduction("fallback backend chosen but absent".into()))?;
                return reduce(fallback);
            }
            None => {}
        }
        match reduce(&self.primary) {
            Ok(reduction) => {
                *choice = Some(Choice::Primary);
                Ok(reduction)
            }
            Err(Error::Cancelled) => Err(Error::Cancelled),
            Err(primary_error) => match &self.fallback {
                Some(fallback) if family.loops.len() < self.fallback_below_loops => {
                    match reduce(fallback) {
                        Ok(reduction) => {
                            *choice = Some(Choice::Fallback);
                            Ok(reduction)
                        }
                        Err(Error::Cancelled) => Err(Error::Cancelled),
                        Err(fallback_error) => Err(Error::Reduction(format!(
                            "{} failed ({primary_error}); fallback {} failed ({fallback_error})",
                            self.primary.identity(),
                            fallback.identity()
                        ))),
                    }
                }
                _ => Err(primary_error),
            },
        }
    }
}

impl<P: ReductionBackend, F: ReductionBackend> ReductionBackend for FallbackBackend<P, F> {
    fn identity(&self) -> String {
        format!(
            "fallback[{}|below {} loops: {}]",
            self.primary.identity(),
            self.fallback_below_loops,
            self.fallback
                .as_ref()
                .map_or_else(|| "none".into(), |f| f.identity())
        )
    }

    fn reduce_at_epsilon(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        epsilon: &Rational,
        context: &RunContext,
    ) -> Result<Reduction> {
        self.run(family, |backend| {
            backend.reduce_at_epsilon(family, targets, epsilon, context)
        })
    }

    fn reduce(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        context: &RunContext,
    ) -> Result<Reduction> {
        self.run(family, |backend| backend.reduce(family, targets, context))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Propagator;
    use std::sync::atomic::{AtomicUsize, Ordering};

    /// Fails on the calls listed in `fail_on` (0-based) and counts its calls.
    struct Scripted {
        name: &'static str,
        fail_on: Vec<usize>,
        calls: AtomicUsize,
    }

    impl Scripted {
        fn new(name: &'static str, fail_on: &[usize]) -> Self {
            Self {
                name,
                fail_on: fail_on.to_vec(),
                calls: AtomicUsize::new(0),
            }
        }
        fn calls(&self) -> usize {
            self.calls.load(Ordering::SeqCst)
        }
    }

    impl ReductionBackend for Scripted {
        fn identity(&self) -> String {
            self.name.into()
        }
        fn reduce(&self, _: &IntegralFamily, _: &[Integral], _: &RunContext) -> Result<Reduction> {
            let call = self.calls.fetch_add(1, Ordering::SeqCst);
            if self.fail_on.contains(&call) {
                return Err(Error::Reduction(format!("{} call {call} fails", self.name)));
            }
            Ok(Reduction {
                residuals: vec![Integral(vec![self.name.len() as i16])],
                ..Default::default()
            })
        }
    }

    fn family(name: &str, loops: usize) -> IntegralFamily {
        let momentum = |i: usize| (0..loops).map(|j| i64::from(i == j)).collect::<Vec<_>>();
        IntegralFamily {
            name: name.into(),
            loops: (0..loops).map(|i| format!("l{i}")).collect(),
            external: vec![],
            external_gram: vec![],
            propagators: (0..loops)
                .map(|i| Propagator::quadratic(&momentum(i), &[], Atom::num(1), &[]).unwrap())
                .collect(),
            physical_propagators: loops,
            epsilon: symbol!("fallback_test::eps"),
            dimension: 4,
        }
    }

    #[test]
    fn the_first_answer_decides_the_backend_of_a_family() {
        let context = RunContext::default();
        let targets = [Integral(vec![1])];
        // Primary answers the first call of family A and fails on its second:
        // the error is returned, the fallback is never consulted for A.
        let backend = FallbackBackend::new(
            Scripted::new("primary", &[1, 2]),
            Some(Scripted::new("fallback", &[])),
            3,
        );
        let a = family("fallback_test_a", 1);
        assert_eq!(
            backend.reduce(&a, &targets, &context).unwrap().residuals,
            vec![Integral(vec![7])]
        );
        let error = backend.reduce(&a, &targets, &context).unwrap_err();
        assert!(error.to_string().contains("not switched"), "{error}");
        assert_eq!(backend.fallback.as_ref().unwrap().calls(), 0);
        // Family B: the primary's first answer fails, so B is reduced by the
        // fallback from then on, even though the primary would answer later.
        let b = family("fallback_test_b", 2);
        assert_eq!(
            backend.reduce(&b, &targets, &context).unwrap().residuals,
            vec![Integral(vec![8])]
        );
        let primary_calls = backend.primary.calls();
        assert_eq!(
            backend
                .reduce_at_epsilon(&b, &targets, &Rational::from((1, 7)), &context)
                .unwrap()
                .residuals,
            vec![Integral(vec![8])]
        );
        assert_eq!(backend.primary.calls(), primary_calls);
        assert_eq!(backend.fallback.as_ref().unwrap().calls(), 2);
        // The top family (loops >= fallback_below_loops) never falls back and
        // a failed first request decides nothing.
        let top = family("fallback_test_top", 3);
        let backend = FallbackBackend::new(
            Scripted::new("primary", &[0]),
            Some(Scripted::new("fallback", &[])),
            3,
        );
        assert!(backend.reduce(&top, &targets, &context).is_err());
        assert!(backend.reduce(&top, &targets, &context).is_ok());
        assert_eq!(backend.fallback.as_ref().unwrap().calls(), 0);
        // Both failing on the first request: both reasons are reported.
        let backend = FallbackBackend::new(
            Scripted::new("primary", &[0]),
            Some(Scripted::new("fallback", &[0])),
            3,
        );
        let error = backend
            .reduce(&a, &targets, &context)
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("primary call 0") && error.contains("fallback call 0"),
            "{error}"
        );
    }
}
