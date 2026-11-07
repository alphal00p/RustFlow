//! Transactional admission of supplied uncertainty after central refinement.
//!
//! Retrying an estimate never repeats transport or mutates its central solution.
use super::{CompiledConnection, errors_meet};
use crate::diffexp::{EpsilonSolution, VerifiedCheckpoint};
use crate::ode::TaylorSegment;
use crate::transport_cache::CachedBoundary;
use crate::{BoundaryErrorStrategy, Error, FlowOptions, Precision, Result, RunContext};
use symbolica::domains::float::{RealBall, RoundingDirection};
use symbolica::prelude::*;

#[cfg(test)]
mod tests;

struct Admission {
    errors: Vec<Vec<Float>>,
    checkpoints: Vec<VerifiedCheckpoint>,
    provenance: Option<String>,
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
    Components(&'a [Float]),
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
                Self::Components(increments) => {
                    *error += &increments[component];
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

impl<'a> Problem<'a> {
    fn candidate<'b>(
        &self,
        solution: &EpsilonSolution,
        mut increment: impl FnMut(usize, &TaylorSegment) -> Result<Increment<'b>>,
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
            provenance: None,
        })
    }

    fn scalar(&self, solution: &EpsilonSolution) -> Result<Admission> {
        // Preserve every successful source-only attempt, including operation
        // order and early failure returns, before trying one fixed norm.
        let original = match self.scalar_disks(solution, 1) {
            Err(error @ Error::Accuracy(_)) => error,
            result => return result,
        };
        for subdivisions in [4, 16] {
            match self.scalar_disks(solution, subdivisions) {
                Ok(candidate) => return Ok(candidate),
                Err(Error::Accuracy(_)) => {}
                Err(error @ (Error::Cancelled | Error::InvalidInput(_))) => return Err(error),
                // Only accuracy failures try another disk size. Retain the
                // original failure and automatic fundamental route otherwise.
                Err(_) => return Err(original),
            }
        }
        self.context.cancellation.check()?;
        let balanced = self.balanced(solution)?;
        for subdivisions in [1, 4, 16] {
            match balanced.scalar_disks(solution, subdivisions) {
                Ok(mut candidate) => {
                    candidate.provenance = Some(format!(
                        "finite-MPFR scalar input-error estimate in fixed balanced weights max(1,|source|,|refined endpoint|)*10^epsilon_offset; {subdivisions} disks per segment"
                    ));
                    return Ok(candidate);
                }
                Err(Error::Accuracy(_)) => {}
                Err(error @ (Error::Cancelled | Error::InvalidInput(_))) => return Err(error),
                Err(_) => return Err(original),
            }
        }
        Err(original)
    }

    fn balanced(&self, solution: &EpsilonSolution) -> Result<Problem<'a>> {
        if solution.coefficients.len() != self.count
            || self.source.coefficients.len() < self.count
            || self.source.accuracy.comparison_errors().len() < self.count
            || self.source.coefficients[..self.count]
                .iter()
                .zip(&solution.coefficients)
                .zip(&self.source.accuracy.comparison_errors()[..self.count])
                .any(|((source, endpoint), errors)| {
                    source.len() != endpoint.len() || source.len() != errors.len()
                })
        {
            return Err(Error::InvalidInput("balanced norm dimensions".into()));
        }
        let mut weights = Vec::with_capacity(self.weights.len());
        let mut epsilon_scale = self.p.real(1);
        for (source, endpoint) in self.source.coefficients[..self.count]
            .iter()
            .zip(&solution.coefficients)
        {
            self.context.cancellation.check()?;
            for (source, endpoint) in source.iter().zip(endpoint) {
                if !self.p.finite(source) || !self.p.finite(endpoint) {
                    return Err(Error::InvalidInput(
                        "nonfinite balanced norm coefficient".into(),
                    ));
                }
                let mut weight = self.p.norm(source);
                let endpoint_norm = self.p.norm(endpoint);
                if endpoint_norm > weight {
                    weight = endpoint_norm;
                }
                if weight < self.p.real(1) {
                    weight = self.p.real(1);
                }
                weight *= &epsilon_scale;
                if !weight.is_finite() || weight <= self.p.real(0) {
                    return Err(Error::InvalidInput("invalid balanced norm weight".into()));
                }
                weights.push(weight);
            }
            epsilon_scale *= self.p.real(10);
        }
        if weights.len() != self.weights.len() {
            return Err(Error::InvalidInput("balanced norm dimensions".into()));
        }
        let mut relative_error = self.p.real(0);
        for (error, weight) in self.source.accuracy.comparison_errors()[..self.count]
            .iter()
            .flatten()
            .zip(&weights)
        {
            let relative = error.clone() / weight;
            if relative > relative_error {
                relative_error = relative;
            }
        }
        // Exactly the source-only floor: endpoint working precision supplies
        // no accuracy credit. Every enlarged weight also enlarges its absolute
        // evidence floor; the maximum covers each recorded comparison error.
        relative_error += self.p.tolerance(self.source.accuracy.verified_digits());
        if !relative_error.is_finite() || relative_error <= self.p.real(0) {
            return Err(Error::InvalidInput(
                "invalid balanced norm source error".into(),
            ));
        }
        Ok(Problem {
            prepared: self.prepared,
            source: self.source,
            options: self.options,
            context: self.context,
            p: self.p,
            count: self.count,
            weights,
            relative_error,
        })
    }

    fn scalar_disks(&self, solution: &EpsilonSolution, subdivisions: i64) -> Result<Admission> {
        let mut integrated_norm = self.p.real(0);
        let result = self.candidate(solution, |_, segment| {
            if subdivisions > 1 {
                // Tighten only the supplied-error disks, keeping the central
                // trajectory, fixed source weights and original checkpoints.
                let delta = self.p.sub(&segment.end, &segment.center);
                let mut start = segment.center.clone();
                for part in 1..=subdivisions {
                    self.context.cancellation.check()?;
                    let end = if part == subdivisions {
                        segment.end.clone()
                    } else {
                        self.p
                            .add(&segment.center, &self.p.scale(&delta, part, subdivisions))
                    };
                    integrated_norm +=
                        self.prepared
                            .error_norm_integral_weighted(&start, &end, &self.weights)?;
                    start = end;
                }
            } else {
                // Preserve the original scalar attempt's arithmetic and rounding.
                integrated_norm += self.prepared.error_norm_integral_weighted(
                    &segment.center,
                    &segment.end,
                    &self.weights,
                )?;
            }
            // Accumulate the whole path before amplification, including the
            // finite Dyson bound when the epsilon hierarchy is nilpotent.
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
    admit_with_gauge(prepared, source, solution, options, context, None).map(|_| ())
}

pub(super) fn admit_with_gauge(
    prepared: &CompiledConnection,
    source: &CachedBoundary,
    solution: &mut EpsilonSolution,
    options: &FlowOptions,
    context: &RunContext,
    gauge: Option<&super::diagonal_gauge::PreparedGauge>,
) -> Result<Option<String>> {
    let p = Precision {
        bits: solution.diagnostics.working_bits,
    };
    let count = solution.coefficients.len();
    context.cancellation.check()?;
    let dimension = source.identity.dimension();
    let shape = |rows: &[Vec<_>]| rows.len() == count && rows.iter().all(|r| r.len() == dimension);
    if count == 0
        || dimension == 0
        || source.coefficients.len() < count
        || source.accuracy.comparison_errors().len() < count
        || !shape(&source.coefficients[..count])
        || !shape(&solution.coefficients)
        || source.accuracy.comparison_errors()[..count]
            .iter()
            .any(|r| r.len() != dimension)
        || solution.comparison_errors.len() != count
        || solution
            .comparison_errors
            .iter()
            .any(|r| r.len() != dimension)
        || solution.checkpoints.iter().any(|c| {
            !shape(&c.coefficients)
                || c.comparison_errors.len() != count
                || c.comparison_errors.iter().any(|r| r.len() != dimension)
        })
    {
        return Err(Error::InvalidInput("input uncertainty dimensions".into()));
    }
    if source.coefficients[..count]
        .iter()
        .flatten()
        .chain(solution.coefficients.iter().flatten())
        .any(|v| !p.finite(v))
        || source.accuracy.comparison_errors()[..count]
            .iter()
            .flatten()
            .any(|e| !e.is_finite() || *e < p.real(0))
    {
        return Err(Error::InvalidInput(
            "nonfinite input uncertainty or negative source error".into(),
        ));
    }
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
    if let Some(gauge) = gauge {
        // Only input uncertainty changes. Borrow this already-refined central
        // trajectory and retain the existing mixed-scale candidate checks.
        let mut first_accuracy = None;
        for subdivisions in [1, 4, 16] {
            context.cancellation.check()?;
            let attempt = gauge
                .increments(
                    solution,
                    &problem.weights,
                    &problem.relative_error,
                    source.accuracy.verified_digits(),
                    subdivisions,
                    context,
                )
                .and_then(|(increments, report)| {
                    let candidate = problem.candidate(solution, |index, _| {
                        Ok(Increment::Components(&increments[index]))
                    })?;
                    Ok((candidate, report))
                });
            match attempt {
                Ok((candidate, report)) => {
                    context.cancellation.check()?;
                    solution.comparison_errors = candidate.errors;
                    solution.checkpoints = candidate.checkpoints;
                    return Ok(Some(report));
                }
                Err(error @ Error::Accuracy(_)) => {
                    if first_accuracy.is_none() {
                        first_accuracy = Some(error);
                    }
                }
                Err(error) => return Err(error),
            }
        }
        return Err(first_accuracy.unwrap());
    }
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
    Ok(candidate.provenance)
}
