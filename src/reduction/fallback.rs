//! Primary reduction backend with a fallback for small boundary children.
use super::{Reduction, ReductionBackend};
use crate::{Error, Integral, IntegralFamily, Result, RunContext};
use symbolica::prelude::*;

/// Reduce with `primary`; when it fails (other than by cancellation) and the
/// family has fewer than `fallback_below_loops` loops, reduce with `fallback`.
///
/// Typical use: a LiteRed lane for every family, with RustRed's native search
/// as a fallback for the lower-loop boundary children the lane rejects, but
/// never for the top family itself (`fallback_below_loops` = its loop count).
pub struct FallbackBackend<P, F> {
    pub primary: P,
    pub fallback: Option<F>,
    pub fallback_below_loops: usize,
}

impl<P: ReductionBackend, F: ReductionBackend> FallbackBackend<P, F> {
    fn run(
        &self,
        family: &IntegralFamily,
        reduce: impl Fn(&dyn ReductionBackend) -> Result<Reduction>,
    ) -> Result<Reduction> {
        match (reduce(&self.primary), &self.fallback) {
            (Err(error), Some(fallback))
                if !matches!(error, Error::Cancelled)
                    && family.loops.len() < self.fallback_below_loops =>
            {
                reduce(fallback)
            }
            (outcome, _) => outcome,
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
