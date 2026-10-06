//! Transactional admission of supplied uncertainty after central refinement.
//!
//! Retrying a proof never repeats transport or mutates its central solution.
use super::{CompiledConnection, errors_meet};
use crate::diffexp::{EpsilonSolution, VerifiedCheckpoint};
use crate::ode::TaylorSegment;
use crate::transport_cache::CachedBoundary;
use crate::{BoundaryErrorStrategy, Error, FlowOptions, Precision, Result, RunContext};
use symbolica::domains::float::{RealBall, RoundingDirection};
use symbolica::prelude::*;

struct Admission {
    errors: Vec<Vec<Float>>,
    checkpoints: Vec<VerifiedCheckpoint>,
}

struct Problem<'a> {
    prepared: &'a CompiledConnection,
    source: &'a CachedBoundary,
    options: &'a FlowOptions,
    context: &'a RunContext,
    p: Precision,
    count: usize,
    weights: Vec<Float>,
    relative_error: Float,
}

enum Increment<'a> {
    Scalar(Float),
    Fundamental(&'a [Float]),
}
impl Increment<'_> {
    fn add(&self, problem: &Problem<'_>, errors: &mut [Vec<Float>]) {
        for (component, (error, weight)) in errors
            .iter_mut()
            .flatten()
            .zip(&problem.weights)
            .enumerate()
        {
            match self {
                // Retain the original scalar operation order and rounding.
                Self::Scalar(amplification) => {
                    *error += weight.clone() * &problem.relative_error * amplification;
                }
                Self::Fundamental(increments) => {
                    *error = error.add_round(
                        &increments[component],
                        problem.p.bits,
                        RoundingDirection::Up,
                    );
                }
            }
        }
    }
}

fn accuracy_failure() -> Error {
    Error::Accuracy("cached boundary uncertainty grows beyond the requested target accuracy; provide a more accurate or closer boundary".into())
}

fn proof_unavailable(original: Error, proof: Error) -> Error {
    match original {
        Error::Accuracy(message) => Error::Accuracy(format!(
            "{message}; optional fundamental proof unavailable: {proof}"
        )),
        // Preserve typed scalar failures unrelated to accuracy admission.
        other => other,
    }
}

impl Problem<'_> {
    fn candidate<'a>(
        &self,
        solution: &EpsilonSolution,
        mut increment: impl FnMut(usize, &TaylorSegment) -> Result<Increment<'a>>,
    ) -> Result<Admission> {
        let mut last = Increment::Scalar(self.p.real(1));
        let mut checkpoints = Vec::new();
        for (index, segment) in solution.segments.iter().enumerate() {
            self.context.cancellation.check()?;
            last = increment(index, segment)?;
            if let Some(checkpoint) = solution.checkpoints.iter().find(|c| c.segment == index) {
                let mut checkpoint = checkpoint.clone();
                last.add(self, &mut checkpoint.comparison_errors);
                if errors_meet(
                    self.p,
                    &checkpoint.coefficients,
                    &checkpoint.comparison_errors,
                    self.options.digits,
                ) {
                    checkpoints.push(checkpoint);
                }
            }
        }
        self.context.cancellation.check()?;
        let mut errors = solution.comparison_errors.clone();
        last.add(self, &mut errors);
        if !errors_meet(self.p, &solution.coefficients, &errors, self.options.digits) {
            return Err(accuracy_failure());
        }
        Ok(Admission {
            errors,
            checkpoints,
        })
    }

    fn scalar(&self, solution: &EpsilonSolution) -> Result<Admission> {
        let mut integrated_norm = self.p.real(0);
        let result = self.candidate(solution, |_, segment| {
            integrated_norm += self.prepared.error_norm_integral_weighted(
                &segment.center,
                &segment.end,
                &self.weights,
            )?;
            Ok(Increment::Scalar(
                crate::diffexp::amplification_from_integral(
                    self.p,
                    &integrated_norm,
                    self.prepared.epsilon_product_limit(),
                )?,
            ))
        });
        self.context.cancellation.check()?;
        result
    }

    fn fundamental(&self, solution: &EpsilonSolution) -> Result<(Admission, usize)> {
        self.context.cancellation.check()?;
        // Mandatory source evidence floor, with outward arithmetic. Storage
        // precision never supplies missing boundary accuracy.
        let floor =
            RealBall::from_rational_ball(&Rational::from((1, 10)), &Rational::zero(), self.p.bits)
                .pow(u64::from(self.source.accuracy.verified_digits()))
                .upper_bound();
        let input_errors = self.source.accuracy.comparison_errors()[..self.count]
            .iter()
            .flatten()
            .zip(&self.weights)
            .map(|(error, weight)| {
                error.add_round(
                    &weight.mul_round(&floor, self.p.bits, RoundingDirection::Up),
                    self.p.bits,
                    RoundingDirection::Up,
                )
            })
            .collect::<Vec<_>>();
        let proof = match self.prepared {
            CompiledConnection::Rational(system) => system.fundamental_boundary_errors(
                &solution.segments,
                &input_errors,
                &self.weights,
                solution.diagnostics.expansion_order,
                self.context,
            )?,
            CompiledConnection::Algebraic(_) => {
                return Err(Error::Unsupported(
                    "registered-root boundary errors retain their scalar owner".into(),
                ));
            }
        };
        let candidate = self.candidate(solution, |index, _| {
            Ok(Increment::Fundamental(&proof.endpoints[index]))
        })?;
        Ok((candidate, proof.endpoints.len()))
    }
}

pub(super) fn admit(
    prepared: &CompiledConnection,
    source: &CachedBoundary,
    solution: &mut EpsilonSolution,
    options: &FlowOptions,
    context: &RunContext,
) -> Result<()> {
    let p = Precision {
        bits: solution.diagnostics.working_bits,
    };
    let count = solution.coefficients.len();
    let weights = source.coefficients[..count]
        .iter()
        .flatten()
        .map(|value| {
            let norm = p.norm(value);
            if norm > p.real(1) { norm } else { p.real(1) }
        })
        .collect::<Vec<_>>();
    let mut relative_error = p.real(0);
    for (error, weight) in source.accuracy.comparison_errors()[..count]
        .iter()
        .flatten()
        .zip(&weights)
    {
        let relative = error.clone() / weight;
        if relative > relative_error {
            relative_error = relative;
        }
    }
    relative_error += p.tolerance(source.accuracy.verified_digits());
    let problem = Problem {
        prepared,
        source,
        options,
        context,
        p,
        count,
        weights,
        relative_error,
    };
    let (candidate, charts, retry, fallback) = match options.boundary_error_strategy {
        BoundaryErrorStrategy::ScalarNorm => (problem.scalar(solution)?, 0, None, None),
        BoundaryErrorStrategy::Automatic => match problem.scalar(solution) {
            Ok(candidate) => (candidate, 0, None, None),
            Err(original @ Error::Accuracy(_)) => {
                let retry = original.to_string();
                match problem.fundamental(solution) {
                    Ok((candidate, charts)) => (candidate, charts, Some(retry), None),
                    Err(error @ (Error::Cancelled | Error::InvalidInput(_))) => return Err(error),
                    Err(proof) => return Err(proof_unavailable(original, proof)),
                }
            }
            Err(error) => return Err(error),
        },
        BoundaryErrorStrategy::FundamentalMatrix => match problem.fundamental(solution) {
            Ok((candidate, charts)) => (candidate, charts, None, None),
            Err(error @ (Error::Cancelled | Error::InvalidInput(_))) => return Err(error),
            Err(proof) => {
                let fallback = proof.to_string();
                match problem.scalar(solution) {
                    Ok(candidate) => (candidate, 0, None, Some(fallback)),
                    Err(original) => return Err(proof_unavailable(original, proof)),
                }
            }
        },
    };
    context.cancellation.check()?;
    // Commit only the admitted uncertainty candidate. Central coefficients,
    // segments and refinement evidence were borrowed throughout both attempts.
    solution.comparison_errors = candidate.errors;
    solution.checkpoints = candidate.checkpoints;
    solution.diagnostics.fundamental_boundary_charts = charts;
    solution.diagnostics.fundamental_boundary_retry = retry;
    solution.diagnostics.fundamental_boundary_fallback = fallback;
    Ok(())
}
