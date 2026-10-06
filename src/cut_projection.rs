//! Exact weighted cut integrals through the existing terminal and flow owners.
use crate::cuts::{CutFamily, FutureTimelikeChannel};
use crate::phase_space::{PreparedMasslessPhaseSpace, PreparedTwoBodyPhaseSpace};
use crate::reduction::LinearCombination;
use crate::*;
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet};
use symbolica::prelude::*;

enum Preparation<'a> {
    TwoBody(PreparedTwoBodyPhaseSpace),
    Massless(PreparedMasslessPhaseSpace),
    Mixed(Box<PreparedCutFlow<'a>>),
}

/// Several exact numerator projections evaluated on shared finite-epsilon values.
///
/// Dispatch is determined by the retained physical denominators before numerical
/// evaluation. Domain errors are propagated without falling back to another
/// measure. The cut routing, Jacobian and all reduction conditions remain owned
/// by the selected terminal or mixed flow. No numerical seed or cache is needed.
pub struct PreparedCutProjections<'a> {
    preparation: Preparation<'a>,
    coefficients: Vec<Vec<Atom>>,
    leading: i32,
}
impl<'a> PreparedCutProjections<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        family: &CutFamily,
        channel: &FutureTimelikeChannel,
        projections: &[LinearCombination],
        point: &KinematicPoint,
        backend: &'a dyn ReductionBackend,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Self> {
        options.validate()?;
        context.cancellation.check()?;
        if projections.is_empty() {
            return Err(Error::InvalidInput("empty cut projection list".into()));
        }
        let (mut ordinary, cuts) = family.at(point)?.into_parts();
        ordinary.dimension = options.dimension;
        let family = CutFamily::new(ordinary, cuts)?;
        let ordinary = family.family();
        let allowed = BTreeSet::from([
            Atom::var(ordinary.epsilon),
            Atom::var(crate::family::imaginary_parameter()),
        ]);
        let mut targets = BTreeSet::new();
        let mut rows = Vec::with_capacity(projections.len());
        let mut valuation = 0;
        for terms in projections {
            context.cancellation.check()?;
            let mut row = BTreeMap::new();
            for (integral, coefficient) in terms {
                context.cancellation.check()?;
                ordinary.validate_integral(integral)?;
                let coefficient = point.apply(coefficient).together().cancel();
                if coefficient.is_zero() {
                    continue;
                }
                valuation = valuation.min(crate::engine::projection_weight_valuation(
                    &coefficient,
                    ordinary.epsilon,
                    &allowed,
                )?);
                targets.insert(integral.clone());
                row.insert(integral.clone(), coefficient);
            }
            rows.push(row);
        }
        // Even identically zero numerators must pass the full measure/domain
        // admission. A unit target supplies that validation without inventing a
        // nonzero output; every projection coefficient stays exactly zero.
        if targets.is_empty() {
            targets.insert(Integral(
                (0..ordinary.propagators.len())
                    .map(|slot| i16::from(slot < ordinary.physical_propagators))
                    .collect(),
            ));
        }
        let targets = targets.into_iter().collect::<Vec<_>>();
        let coefficients = rows
            .iter()
            .map(|row| {
                targets
                    .iter()
                    .map(|i| row.get(i).cloned().unwrap_or_default())
                    .collect()
            })
            .collect();
        let bound = i32::try_from(ordinary.loops.len())
            .ok()
            .and_then(|loops| loops.checked_mul(-2))
            .ok_or_else(|| Error::Limit("cut projection pole bound overflow".into()))?;
        let default = KinematicPoint::default();
        let (preparation, bound) = if family.cuts().lines().len() == ordinary.physical_propagators {
            if ordinary.loops.len() == 1
                && ordinary.external.len() == 1
                && ordinary.propagators.len() == 2
            {
                let prepared = PreparedTwoBodyPhaseSpace::new(
                    &family, channel, &targets, &default, backend, options, context,
                )?;
                (Preparation::TwoBody(prepared), bound)
            } else {
                let prepared = PreparedMasslessPhaseSpace::new(
                    &family, channel, &targets, &default, backend, options, context,
                )?;
                let bound = bound.min(prepared.leading_power());
                (Preparation::Massless(prepared), bound)
            }
        } else {
            let prepared = PreparedCutFlow::new(
                &family, channel, &targets, &default, backend, options, context,
            )?;
            (Preparation::Mixed(Box::new(prepared)), bound)
        };
        let leading = bound
            .checked_add(valuation)
            .ok_or_else(|| Error::Limit("weighted cut Laurent bound overflow".into()))?;
        context.cancellation.check()?;
        Ok(Self {
            preparation,
            coefficients,
            leading,
        })
    }

    pub fn family(&self) -> &CutFamily {
        match &self.preparation {
            Preparation::TwoBody(p) => p.family(),
            Preparation::Massless(p) => p.family(),
            Preparation::Mixed(p) => p.family(),
        }
    }

    pub fn leading_power(&self) -> i32 {
        self.leading
    }

    /// All exact preparation conditions; the owner also checks them at samples.
    pub fn nonzero_conditions(&self) -> &[Atom] {
        match &self.preparation {
            Preparation::TwoBody(p) => p.nonzero_conditions(),
            Preparation::Massless(p) => p.nonzero_conditions(),
            Preparation::Mixed(p) => p.nonzero_conditions(),
        }
    }

    pub fn evaluate(
        &self,
        epsilon: &Rational,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<ComplexFloat>> {
        options.validate()?;
        context.cancellation.check()?;
        if epsilon.is_zero() {
            return Err(Error::InvalidInput(
                "epsilon samples must be nonzero".into(),
            ));
        }
        let values = match &self.preparation {
            Preparation::TwoBody(p) => p.evaluate(epsilon, options, context)?,
            Preparation::Massless(p) => p.evaluate(epsilon, options, context)?,
            Preparation::Mixed(p) => p.evaluate(epsilon, options, context)?,
        };
        let p = Precision::decimal(options.digits + options.guard_digits)?;
        let parameters = [(
            Atom::var(self.family().family().epsilon),
            p.rational(epsilon),
        )]
        .into_iter()
        .collect();
        self.coefficients
            .iter()
            .map(|row| {
                context.cancellation.check()?;
                let mut sum = p.zero();
                for (coefficient, value) in row.iter().zip(&values) {
                    context.cancellation.check()?;
                    if !coefficient.is_zero() {
                        sum = p.add(&sum, &p.mul(&p.eval(coefficient, &parameters)?, value));
                    }
                }
                context.cancellation.check()?;
                if !p.finite(&sum) {
                    return Err(Error::Numerical("nonfinite weighted cut sample".into()));
                }
                Ok(sum)
            })
            .collect()
    }

    pub fn evaluate_samples(
        &self,
        samples: &[Rational],
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<Vec<ComplexFloat>>> {
        options.validate()?;
        context.cancellation.check()?;
        if samples.is_empty() || samples.iter().any(Rational::is_zero) {
            return Err(Error::InvalidInput(
                "at least one nonzero epsilon sample is required".into(),
            ));
        }
        let evaluate = |(index, epsilon): (usize, &Rational)| {
            context.emit(Progress::Sample {
                index,
                total: samples.len(),
            })?;
            self.evaluate(epsilon, options, context)
        };
        if options.workers == 1 {
            samples.iter().enumerate().map(evaluate).collect()
        } else {
            rayon::ThreadPoolBuilder::new()
                .num_threads(options.workers)
                .build()
                .map_err(|e| Error::InvalidInput(e.to_string()))?
                .install(|| samples.par_iter().enumerate().map(evaluate).collect())
        }
    }

    pub fn solve(
        &self,
        last: i32,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<LaurentExpansion>> {
        context.cancellation.check()?;
        let result = crate::engine::fit_samples_refined_leading(
            self.coefficients.len(),
            self.leading,
            last,
            options,
            |samples, refined| self.evaluate_samples(samples, refined, context),
        )?;
        context.cancellation.check()?;
        Ok(result)
    }
}
