//! Preparation, auxiliary-mass flow, and independent epsilon-fit refinement.
use crate::boundary::BoundaryProvider;
use crate::reduction::{ReducedSystem, differential_system};
use crate::*;
use rayon::prelude::*;
use std::collections::BTreeMap;
use symbolica::prelude::*;

pub struct PreparedFlow {
    pub family: IntegralFamily,
    pub reduced: ReducedSystem,
    pub system: DifferentialSystem,
    epsilon_sample: Option<Rational>,
    pub basis_refinement: Option<crate::refine::RefinementReport>,
}

impl PreparedFlow {
    pub fn new(
        family: &IntegralFamily,
        targets: &[Integral],
        point: &KinematicPoint,
        backend: &dyn ReductionBackend,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Self> {
        context.cancellation.check()?;
        options.validate()?;
        if options.recursion != RecursionMode::Amf {
            return Err(Error::Unsupported(
                "PreparedFlow prepares AMF systems; use solve_integrals for FT recursion".into(),
            ));
        }
        if targets.is_empty() {
            return Err(Error::InvalidInput("empty target list".into()));
        }
        let mut family = family.at(point);
        family.dimension = options.dimension;
        family.validate()?;
        for target in targets {
            family.validate_integral(target)?;
        }
        let eta = symbol!("symbolica_amflow::eta");
        let cache = options
            .cache_directory
            .as_ref()
            .map(|directory| -> Result<_> {
                Ok((
                    directory,
                    crate::cache::system_key(&family, targets, backend, options)?,
                ))
            })
            .transpose()?;
        if let Some((directory, key)) = &cache
            && let Some((reduced, basis_refinement)) = crate::cache::read_system(directory, key)?
        {
            let system = DifferentialSystem {
                variable: eta,
                matrix: reduced.matrix.clone(),
            };
            return Ok(Self {
                family,
                reduced,
                system,
                basis_refinement,
                epsilon_sample: None,
            });
        }
        let (auxiliary, mask) = family.deform(eta, &options.mass_mode)?;
        let session = crate::reduction::ReductionSession::new(backend);
        let backend = &session as &dyn ReductionBackend;
        let mut reduced = if options.skip_reduction {
            crate::reduction::differential_system_skip_initial(
                backend, &auxiliary, targets, &mask, 12, context,
            )?
        } else {
            differential_system(backend, &auxiliary, targets, &mask, 12, context)?
        };
        let basis_refinement = if options.refine_basis {
            Some(crate::refine::refine_basis(
                &mut reduced,
                eta,
                family.epsilon,
                32,
            )?)
        } else {
            None
        };
        if let Some((directory, key)) = cache {
            crate::cache::write_system(directory, &key, &reduced, &basis_refinement)?;
        }
        let system = DifferentialSystem {
            variable: eta,
            matrix: reduced.matrix.clone(),
        };
        context.emit(Progress::Prepared {
            basis_size: reduced.basis.len(),
            blocks: if reduced.basis.is_empty() {
                vec![]
            } else {
                system.blocks()?.iter().map(Vec::len).collect()
            },
        })?;
        Ok(Self {
            family,
            reduced,
            system,
            basis_refinement,
            epsilon_sample: None,
        })
    }
    /// Prepare at an exact regulator sample. Dimension-factorizing refinement
    /// retains symbolic epsilon instead, producing a reusable symbolic flow.
    pub fn new_at_epsilon(
        family: &IntegralFamily,
        targets: &[Integral],
        point: &KinematicPoint,
        backend: &dyn ReductionBackend,
        options: &FlowOptions,
        context: &RunContext,
        epsilon: &Rational,
    ) -> Result<Self> {
        if options.refine_basis {
            return Self::new(family, targets, point, backend, options, context);
        }
        let sampled = crate::reduction::SampledBackend {
            backend,
            epsilon: epsilon.clone(),
        };
        let mut flow = Self::new(family, targets, point, &sampled, options, context)?;
        flow.epsilon_sample = Some(epsilon.clone());
        Ok(flow)
    }

    fn lift_exponents(&self, basis: &mut frobenius::FrobeniusBasis) -> Result<()> {
        let Some(epsilon) = &self.epsilon_sample else {
            return Ok(());
        };
        let limit = 2 * self.family.loops.len() as i64;
        for column in &mut basis.columns {
            let value = if let AtomView::Num(n) = column.exponent.as_view() {
                if let symbolica::coefficient::Coefficient::Complex(c) =
                    n.get_coeff_view().to_owned()
                {
                    if !c.im.is_zero() {
                        return Err(Error::Unsupported("complex sampled indicial root".into()));
                    }
                    c.re
                } else {
                    return Err(Error::Unsupported(
                        "nonrational sampled indicial root".into(),
                    ));
                }
            } else {
                return Err(Error::Unsupported(
                    "nonnumeric sampled indicial root".into(),
                ));
            };
            let candidates = (-limit..=limit)
                .filter_map(|slope| {
                    let base = &value - &(epsilon * &Rational::from(slope));
                    (&base * &Rational::from(2))
                        .is_integer()
                        .then_some((base, slope))
                })
                .collect::<Vec<_>>();
            if candidates.len() != 1 {
                return Err(Error::Unsupported(
                    "sample does not uniquely identify dimensional indicial sectors".into(),
                ));
            }
            let (base, slope) = &candidates[0];
            column.exponent =
                Atom::num(base.clone()) + Atom::num(*slope) * Atom::var(self.family.epsilon);
        }
        Ok(())
    }

    pub fn evaluate(
        &self,
        eps: &Rational,
        options: &FlowOptions,
        boundary: &dyn BoundaryProvider,
        context: &RunContext,
    ) -> Result<Vec<ComplexFloat>> {
        if self
            .epsilon_sample
            .as_ref()
            .is_some_and(|sample| sample != eps)
        {
            return Err(Error::InvalidInput(
                "prepared system belongs to a different epsilon sample".into(),
            ));
        }
        options.validate()?;
        context.cancellation.check()?;
        let digits = options
            .digits
            .checked_add(options.guard_digits)
            .ok_or_else(|| Error::InvalidInput("precision overflow".into()))?;
        let p = Precision::decimal(digits)?;
        if self.reduced.basis.is_empty() {
            return Ok(vec![p.zero(); self.reduced.targets.len()]);
        }
        let parameters =
            ahash::HashMap::from_iter([(Atom::var(self.family.epsilon), p.rational(eps))]);
        for condition in &self.reduced.nonzero_conditions {
            if condition.derivative(self.system.variable).is_zero()
                && p.eval(condition, &parameters)? == p.zero()
            {
                return Err(Error::Reduction(format!(
                    "a reduction nonzero condition vanishes at epsilon={eps}"
                )));
            }
        }
        let compiled = self.system.compile(p, &parameters)?;
        let mut infinity = self
            .system
            .invert_variable(symbol!("symbolica_amflow::z"))
            .frobenius(p, &parameters, options.series_order)?;
        self.lift_exponents(&mut infinity)?;
        ensure_generic_indicial(&infinity, &parameters)?;
        let constants = boundary.constants(&self.family, &self.reduced.basis, eps, p, &infinity)?;
        let maximum = compiled
            .poles
            .iter()
            .map(|v| p.norm(v))
            .fold(p.real(1), |a, b| if a > b { a } else { b });
        let direction = if options.prescription == Prescription::PlusI0 {
            -1
        } else {
            1
        };
        let start = ComplexFloat::new(p.real(0), maximum * (8 * direction));
        let at_start = infinity.evaluate(&p.div(&p.i(1), &start), &parameters)?;
        let initial = at_start
            .iter()
            .map(|row| {
                row.iter()
                    .zip(&constants)
                    .fold(p.zero(), |s, (v, c)| p.add(&s, &p.mul(v, c)))
            })
            .collect();
        let minimum = compiled
            .poles
            .iter()
            .map(|v| p.norm(v))
            .filter(|v| *v > p.real(0))
            .fold(p.real(1), |a, b| if a < b { a } else { b });
        let end = ComplexFloat::new(p.real(0), (minimum / 8) * direction);
        let path = compiled.plan_path(&start, &end, -direction)?;
        let transported = compiled.transport(
            &BoundaryData {
                point: start,
                values: initial,
            },
            &path,
            options,
            context,
        )?;
        let mut endpoint = self
            .system
            .frobenius(p, &parameters, options.series_order)?;
        self.lift_exponents(&mut endpoint)?;
        ensure_generic_indicial(&endpoint, &parameters)?;
        let constants = endpoint.match_values(&end, &transported.values, &parameters)?;
        // Reduction coefficients are multiplied into endpoint series before
        // selecting the physical constant; poles in these coefficients matter.
        self.reduced
            .targets
            .iter()
            .map(|target| {
                let weights = self
                    .reduced
                    .basis
                    .iter()
                    .map(|i| target.get(i).cloned().unwrap_or_default())
                    .collect::<Vec<_>>();
                project_limit(
                    &endpoint,
                    &constants,
                    &weights,
                    self.family.epsilon,
                    self.system.variable,
                    &parameters,
                )
            })
            .collect()
    }
}

fn ensure_generic_indicial(
    basis: &frobenius::FrobeniusBasis,
    parameters: &ahash::HashMap<Atom, ComplexFloat>,
) -> Result<()> {
    let p = basis.precision;
    for (i, a) in basis.columns.iter().enumerate() {
        for b in &basis.columns[..i] {
            let difference = (&a.exponent - &b.exponent).together().cancel();
            if difference.to_string().parse::<i64>().is_ok() {
                continue;
            }
            let value = p.eval(&difference, parameters)?;
            if let Some(integer) = value.re.as_raw().to_integer()
                && let Ok(integer) = integer.to_string().parse::<i64>()
                && p.close(&value, &p.i(integer), p.bits / 5)
            {
                return Err(Error::Unsupported("epsilon specialization introduces an additional indicial resonance; choose a generic nonzero sample".into()));
            }
        }
    }
    Ok(())
}

pub(crate) fn project_limit(
    basis: &frobenius::FrobeniusBasis,
    constants: &[ComplexFloat],
    weights: &[Atom],
    epsilon: Symbol,
    eta: Symbol,
    parameters: &ahash::HashMap<Atom, ComplexFloat>,
) -> Result<ComplexFloat> {
    basis.validate()?;
    let p = basis.precision;
    if constants.len() != basis.columns.len() || weights.len() != basis.columns.len() {
        return Err(Error::InvalidInput("endpoint projection dimensions".into()));
    }
    let mut terms: BTreeMap<(Rational, usize), ComplexFloat> = BTreeMap::new();
    for (column, constant) in basis.columns.iter().zip(constants) {
        if *constant == p.zero() || !column.exponent.derivative(epsilon).is_zero() {
            continue;
        }
        let lambda = if let AtomView::Num(n) = column.exponent.as_view() {
            if let symbolica::coefficient::Coefficient::Complex(c) = n.get_coeff_view().to_owned() {
                if !c.im.is_zero() {
                    return Err(Error::Unsupported(
                        "complex physical endpoint exponent".into(),
                    ));
                }
                c.re
            } else {
                return Err(Error::Unsupported("nonrational endpoint exponent".into()));
            }
        } else {
            return Err(Error::Unsupported(
                "parameter-dependent physical endpoint exponent".into(),
            ));
        };
        for (i, weight) in weights.iter().enumerate() {
            if weight.is_zero() {
                continue;
            }
            let valuation = crate::frobenius::valuation(weight, eta)?;
            if &lambda + &Rational::from(valuation + column.coefficients.len() as i64)
                <= Rational::zero()
            {
                return Err(Error::Accuracy(
                    "endpoint series does not reach every term required by the target reduction"
                        .into(),
                ));
            }
            let through = (-lambda.clone())
                .floor()
                .to_string()
                .parse::<i64>()
                .map_err(|_| Error::Limit("endpoint exponent exceeds index range".into()))?;
            if through < valuation {
                continue;
            }
            if through.saturating_sub(valuation) > 10000 {
                return Err(Error::Limit(
                    "endpoint target expansion exceeds 10000 terms".into(),
                ));
            }
            let expansion = weight
                .series(eta, 0, through + 1)
                .map_err(|e| Error::Unsupported(e.to_string()))?;
            for (power, coefficient) in expansion.terms() {
                let power = power
                    .to_string()
                    .parse::<i64>()
                    .map_err(|_| Error::Unsupported("noninteger reduction exponent".into()))?;
                let coefficient = p.eval(coefficient, parameters)?;
                for (k, logs) in column.coefficients.iter().enumerate() {
                    let exponent = &lambda + &Rational::from(power + k as i64);
                    if exponent > Rational::zero() {
                        break;
                    }
                    for (l, row) in logs.iter().enumerate() {
                        let term = p.mul(constant, &p.mul(&coefficient, &row[i]));
                        let entry = terms
                            .entry((exponent.clone(), l))
                            .or_insert_with(|| p.zero());
                        *entry = p.add(entry, &term);
                    }
                }
            }
        }
    }
    for ((power, log), value) in &terms {
        if (*power < Rational::zero() || *log > 0) && p.norm(value) > p.tolerance(p.bits / 5) {
            return Err(Error::Numerical(
                "uncancelled physical endpoint divergence".into(),
            ));
        }
    }
    Ok(terms
        .remove(&(Rational::zero(), 0))
        .unwrap_or_else(|| p.zero()))
}

pub fn evaluate_samples(
    prepared: &PreparedFlow,
    samples: &[Rational],
    options: &FlowOptions,
    boundary: &dyn BoundaryProvider,
    context: &RunContext,
) -> Result<Vec<Vec<ComplexFloat>>> {
    options.validate()?;
    if samples.iter().any(|v| v.is_zero()) {
        return Err(Error::InvalidInput(
            "epsilon samples must be nonzero".into(),
        ));
    }
    if options.workers == 1 {
        return samples
            .iter()
            .enumerate()
            .map(|(index, eps)| {
                context.emit(Progress::Sample {
                    index,
                    total: samples.len(),
                })?;
                prepared.evaluate(eps, options, boundary, context)
            })
            .collect();
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.workers)
        .build()
        .map_err(|e| Error::InvalidInput(e.to_string()))?;
    pool.install(|| {
        samples
            .par_iter()
            .enumerate()
            .map(|(index, eps)| {
                context.emit(Progress::Sample {
                    index,
                    total: samples.len(),
                })?;
                prepared.evaluate(eps, options, boundary, context)
            })
            .collect()
    })
}

pub fn solve_integrals(
    family: &IntegralFamily,
    targets: &[Integral],
    point: &KinematicPoint,
    last: i32,
    options: &FlowOptions,
    backend: &dyn ReductionBackend,
    context: &RunContext,
) -> Result<Vec<LaurentExpansion>> {
    options.validate()?;
    if targets.is_empty() {
        return Err(Error::InvalidInput("empty target list".into()));
    }
    if options.recursion == RecursionMode::Ft {
        let mut family = family.at(point);
        family.dimension = options.dimension;
        let evaluator = crate::ft::FtEvaluator::new(backend, context);
        return fit_samples_refined(
            targets.len(),
            family.loops.len(),
            last,
            options,
            |samples, refined| {
                let evaluate = |(index, eps): (usize, &Rational)| {
                    context.emit(Progress::Sample {
                        index,
                        total: samples.len(),
                    })?;
                    targets
                        .iter()
                        .map(|target| evaluator.evaluate(&family, target, eps, refined))
                        .collect::<Result<Vec<_>>>()
                };
                if refined.workers == 1 {
                    samples.iter().enumerate().map(evaluate).collect()
                } else {
                    let pool = rayon::ThreadPoolBuilder::new()
                        .num_threads(refined.workers)
                        .build()
                        .map_err(|e| Error::InvalidInput(e.to_string()))?;
                    pool.install(|| samples.par_iter().enumerate().map(evaluate).collect())
                }
            },
        );
    }
    if options.sampled_reduction && !options.refine_basis && family.loops.len() > 1 {
        let boundary = crate::recursive::RecursiveBoundary::new(backend, options, context);
        return fit_samples_refined(
            targets.len(),
            family.loops.len(),
            last,
            options,
            |samples, refined| {
                let evaluate = |(index, epsilon): (usize, &Rational)| {
                    context.emit(Progress::Sample {
                        index,
                        total: samples.len(),
                    })?;
                    let prepared = PreparedFlow::new_at_epsilon(
                        family, targets, point, backend, refined, context, epsilon,
                    )?;
                    prepared.evaluate(epsilon, refined, &boundary, context)
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
            },
        );
    }
    let prepared = PreparedFlow::new(family, targets, point, backend, options, context)?;
    let boundary = crate::recursive::RecursiveBoundary::new(backend, options, context);
    solve_prepared(&prepared, last, options, &boundary, context)
}

pub fn solve_prepared(
    prepared: &PreparedFlow,
    last: i32,
    options: &FlowOptions,
    boundary: &dyn BoundaryProvider,
    context: &RunContext,
) -> Result<Vec<LaurentExpansion>> {
    fit_samples_refined(
        prepared.reduced.targets.len(),
        prepared.family.loops.len(),
        last,
        options,
        |samples, refined| evaluate_samples(prepared, samples, refined, boundary, context),
    )
}

fn fit_samples_refined(
    targets: usize,
    loops: usize,
    last: i32,
    options: &FlowOptions,
    evaluate: impl Fn(&[Rational], &FlowOptions) -> Result<Vec<Vec<ComplexFloat>>>,
) -> Result<Vec<LaurentExpansion>> {
    options.validate()?;
    let leading = -2 * (loops as i32);
    if last < leading {
        return Err(Error::InvalidInput(
            "last epsilon power precedes leading pole bound".into(),
        ));
    }
    let count =
        (i64::from(last) - i64::from(leading) + 1) as usize + (options.digits as usize / 2) + 12;
    let count = count + (1 - count % 2);
    if count > 10000 {
        return Err(Error::Limit(
            "epsilon reconstruction exceeds 10000 samples".into(),
        ));
    }
    let mut previous: Option<Vec<LaurentExpansion>> = None;
    for attempt in 0..=options.max_precision_attempts {
        let mut refined = options.clone();
        refined.guard_digits = refined
            .guard_digits
            .checked_add(attempt as u32 * 20)
            .ok_or_else(|| Error::Limit("precision refinement overflow".into()))?;
        refined.series_order += attempt * 16;
        let samples = epsilon::epsilon_samples(
            count + attempt * 4,
            100 * (1_i64
                .checked_shl(attempt as u32)
                .ok_or_else(|| Error::Limit("too many precision attempts".into()))?),
        )?;
        let p = Precision::decimal(refined.digits + refined.guard_digits)?;
        let values = match evaluate(&samples, &refined) {
            Err(Error::Accuracy(_) | Error::Numerical(_))
                if attempt < options.max_precision_attempts =>
            {
                previous = None;
                continue;
            }
            result => result?,
        };
        let fits = (0..targets)
            .map(|i| {
                fit_epsilon(
                    &samples,
                    &values.iter().map(|row| row[i].clone()).collect::<Vec<_>>(),
                    leading,
                    last,
                    p,
                )
            })
            .collect::<Result<Vec<_>>>();
        let mut fits = match fits {
            Err(Error::Accuracy(_) | Error::Numerical(_))
                if attempt < options.max_precision_attempts =>
            {
                previous = None;
                continue;
            }
            result => result?,
        };
        if let Some(old) = &previous
            && old.iter().zip(&fits).all(|(a, b)| {
                a.coefficients.iter().all(|(k, v)| {
                    let next = &b.coefficients[k];
                    if !p.finite(v) || !p.finite(next) {
                        return false;
                    }
                    let old_magnitude = p.norm(v);
                    let new_magnitude = p.norm(next);
                    let magnitude = if old_magnitude > new_magnitude {
                        old_magnitude
                    } else {
                        new_magnitude
                    };
                    let scale = if magnitude <= p.tolerance(options.digits) {
                        p.real(1)
                    } else {
                        magnitude
                    };
                    p.norm(&p.sub(v, next)) <= p.tolerance(options.digits) * scale
                })
            })
        {
            for (f, previous) in fits.iter_mut().zip(old) {
                f.verified_digits = Some(options.digits);
                f.validation_samples = previous.samples;
                f.refinements = attempt;
                f.comparison_errors = f
                    .coefficients
                    .iter()
                    .map(|(k, c)| (*k, p.norm(&p.sub(c, &previous.coefficients[k]))))
                    .collect();
            }
            return Ok(fits);
        }
        previous = Some(fits);
    }
    Err(Error::Accuracy(
        "epsilon coefficients did not stabilize under independent sample/precision refinement"
            .into(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fit_refinement_checks_relative_accuracy_for_small_coefficients() {
        for guard_digits in [40, 0] {
            let options = FlowOptions {
                guard_digits,
                ..Default::default()
            };
            let fits = fit_samples_refined(1, 1, 0, &options, |samples, options| {
                let p = Precision::decimal(options.digits + options.guard_digits)?;
                Ok(samples
                    .iter()
                    .map(|eps| {
                        let x = p.rational(eps);
                        vec![p.scale(
                            &p.add(&p.div(&p.i(1), &x), &p.div(&p.i(1), &p.sub(&p.i(1), &x))),
                            1,
                            1_000_000,
                        )]
                    })
                    .collect())
            })
            .unwrap();
            let p = Precision::decimal(80).unwrap();
            let expected = p.scale(&p.i(1), 1, 1_000_000);
            assert_eq!(fits[0].verified_digits, Some(20));
            for power in [-1, 0] {
                assert!(
                    p.norm(&p.sub(&fits[0].coefficients[&power], &expected))
                        < p.tolerance(20) * p.norm(&expected)
                );
                assert!(fits[0].comparison_errors[&power] < p.tolerance(20) * p.norm(&expected));
            }
        }
    }
}
