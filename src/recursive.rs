//! Recursive AMF boundary evaluation with exact region and tensor algebra.
use crate::boundary::{BoundaryProvider, LeadingBoundary, OneLoopBoundary, RegionBoundary};
use crate::frobenius::FrobeniusBasis;
use crate::*;
use std::{
    collections::BTreeMap,
    sync::{Arc, Mutex},
};
use symbolica::prelude::*;

/// Optional user terminal integrals. Returning None delegates to the native
/// terminals and recursive AMF construction.
pub trait TerminalProvider: Send + Sync {
    fn evaluate(
        &self,
        family: &IntegralFamily,
        integral: &Integral,
        epsilon: &Rational,
        precision: Precision,
    ) -> Result<Option<ComplexFloat>>;
}

#[derive(Default)]
struct Memo {
    flows: BTreeMap<String, Arc<PreparedFlow>>,
    values: BTreeMap<String, ComplexFloat>,
}

/// Automatic boundary provider. All recursive children share symbolic and
/// numerical memo tables; an active ancestry detects recursive cycles.
pub struct RecursiveBoundary<'a> {
    backend: &'a dyn ReductionBackend,
    options: &'a FlowOptions,
    context: &'a RunContext,
    memo: Arc<Mutex<Memo>>,
    ancestry: Vec<String>,
    terminal: Option<&'a dyn TerminalProvider>,
    vacuum_ft: Arc<crate::ft::FtEvaluator<'a>>,
}
impl<'a> RecursiveBoundary<'a> {
    pub fn new(
        backend: &'a dyn ReductionBackend,
        options: &'a FlowOptions,
        context: &'a RunContext,
    ) -> Self {
        Self {
            backend,
            options,
            context,
            memo: Arc::default(),
            ancestry: vec![],
            terminal: None,
            vacuum_ft: Arc::new(crate::ft::FtEvaluator::new(backend, context)),
        }
    }
    pub fn with_terminal(mut self, terminal: &'a dyn TerminalProvider) -> Self {
        self.terminal = Some(terminal);
        self
    }
    fn numerical_options(&self, p: Precision) -> FlowOptions {
        let mut options = self.options.clone();
        let working_digits = ((u64::from(p.bits.saturating_sub(16)) * 1000) / 3322) as u32;
        let extra = working_digits.saturating_sub(options.digits + options.guard_digits);
        options.guard_digits = options
            .guard_digits
            .max(working_digits.saturating_sub(options.digits));
        options.series_order += extra as usize * 4 / 5;
        options
    }
    fn child(&self, key: String) -> Self {
        let mut ancestry = self.ancestry.clone();
        ancestry.push(key);
        Self {
            backend: self.backend,
            options: self.options,
            context: self.context,
            memo: self.memo.clone(),
            ancestry,
            terminal: self.terminal,
            vacuum_ft: self.vacuum_ft.clone(),
        }
    }
    pub fn evaluate(
        &self,
        family: &IntegralFamily,
        target: &Integral,
        epsilon: &Rational,
        p: Precision,
    ) -> Result<ComplexFloat> {
        self.context.cancellation.check()?;
        family.validate_integral(target)?;
        let key = format!(
            "{}:{family:?}:{:?}",
            family.convert()?.family.fingerprint(),
            target.0
        );
        let problem = family.integral_key(target)?;
        let value_key = format!("{problem}:{epsilon}:{}", p.bits);
        if let Some(v) = self
            .memo
            .lock()
            .map_err(|_| Error::Numerical("boundary memo poisoned".into()))?
            .values
            .get(&value_key)
        {
            return Ok(v.clone());
        }
        if self.ancestry.contains(&problem) {
            return Err(Error::Limit("recursive boundary cycle detected".into()));
        }
        if self.ancestry.len() >= 32 {
            return Err(Error::Limit("boundary recursion exceeds 32 levels".into()));
        }
        let custom = if let Some(provider) = self.terminal {
            provider.evaluate(family, target, epsilon, p)?
        } else {
            None
        };
        let value = if let Some(value) = custom {
            if !p.finite(&value) {
                return Err(Error::Numerical(
                    "nonfinite terminal provider result".into(),
                ));
            }
            value
        } else if scaleless(family, target)? {
            p.zero()
        } else if let Some(value) = vacuum_terminal(family, target, epsilon, p)? {
            value
        } else if family.external.is_empty()
            && family
                .propagators
                .iter()
                .zip(&target.0)
                .filter(|(d, n)| **n > 0 && !d.constant.is_zero())
                .count()
                == 1
        {
            // Single-mass vacuum systems are scale invariant under AMF on the
            // massive line. Finish unrecognized topologies by Euclidean FT,
            // whose recursion terminates at a Gaussian integral.
            self.vacuum_ft
                .evaluate(family, target, epsilon, &self.numerical_options(p))?
        } else {
            let mut options = self.numerical_options(p);
            // Keep the requested sampling policy in recursive multiloop
            // families too. At large exceptional samples distinct dimensional
            // exponent classes can coincide; retain symbolic epsilon there.
            let sampled = options.sampled_reduction
                && !options.refine_basis
                && family.loops.len() > 1
                && !(1..=4 * family.loops.len()).any(|difference| {
                    (epsilon * &Rational::from(2 * difference as i64)).is_integer()
                });
            let key = if sampled {
                format!("{key}:epsilon={epsilon}")
            } else {
                key
            };
            // Recursive vacuum families must lose mass scales. Shifting every
            // line of an equal-mass vacuum would reproduce the same problem.
            options.mass_mode = if family.external.is_empty() {
                let selected = family
                    .propagators
                    .iter()
                    .zip(&target.0)
                    .position(|(d, &n)| n > 0 && !d.constant.is_zero())
                    .or_else(|| target.0.iter().position(|&n| n > 0))
                    .ok_or_else(|| Error::InvalidInput("empty boundary sector".into()))?;
                MassMode::Propagators(vec![selected])
            } else {
                MassMode::All
            };
            let cached = {
                self.memo
                    .lock()
                    .map_err(|_| Error::Numerical("boundary memo poisoned".into()))?
                    .flows
                    .get(&key)
                    .cloned()
            };
            let prepared = if let Some(v) = cached {
                v
            } else {
                let v = Arc::new(if sampled {
                    PreparedFlow::new_at_epsilon(
                        family,
                        std::slice::from_ref(target),
                        &KinematicPoint::default(),
                        self.backend,
                        &options,
                        self.context,
                        epsilon,
                    )?
                } else {
                    PreparedFlow::new(
                        family,
                        std::slice::from_ref(target),
                        &KinematicPoint::default(),
                        self.backend,
                        &options,
                        self.context,
                    )?
                });
                self.memo
                    .lock()
                    .map_err(|_| Error::Numerical("boundary memo poisoned".into()))?
                    .flows
                    .insert(key.clone(), v.clone());
                v
            };
            let child = self.child(problem);
            let provider = RecursiveBoundary {
                options: &options,
                ..child
            };
            prepared.evaluate(epsilon, &options, &provider, self.context)?[0].clone()
        };
        self.memo
            .lock()
            .map_err(|_| Error::Numerical("boundary memo poisoned".into()))?
            .values
            .insert(value_key, value.clone());
        Ok(value)
    }
}

pub fn scaleless(family: &IntegralFamily, integral: &Integral) -> Result<bool> {
    use rustred::sector::{
        Mask,
        zero::{Analyzer, Decision},
    };
    let converted = family.convert()?;
    let analyzer = Analyzer::try_unrestricted(&converted.family)
        .map_err(|e| Error::Reduction(e.to_string()))?;
    let mask = Mask::try_new(integral.0.iter().map(|&n| n > 0))
        .map_err(|e| Error::InvalidInput(e.to_string()))?;
    Ok(matches!(
        analyzer
            .analyze(&mask)
            .map_err(|e| Error::Reduction(e.to_string()))?,
        Decision::ProvedZero(_)
    ))
}

fn vacuum_terminal(
    family: &IntegralFamily,
    target: &Integral,
    epsilon: &Rational,
    p: Precision,
) -> Result<Option<ComplexFloat>> {
    if !family.external.is_empty() {
        return Ok(None);
    }
    let dimension = Rational::from(family.dimension) - epsilon * &Rational::from(2);
    let params = ahash::HashMap::from_iter([(Atom::var(family.epsilon), p.rational(epsilon))]);
    if family.loops.len() == 1 && target.0.len() == 1 && target.0[0] > 0 {
        let d = &family.propagators[0];
        let scale = p.eval(&d.scalar_products[0], &params)?;
        let m = p.neg(&p.div(&p.eval(&d.constant, &params)?, &scale));
        return Ok(Some(p.mul(
            &p.powi(&scale, -i64::from(target.0[0])),
            &vacuum::tadpole(target.0[0] as u32, &m, &dimension, p)?,
        )));
    }
    if family.loops.len() != 2 || target.0.len() != 3 || target.0.iter().any(|&n| n <= 0) {
        return Ok(None);
    }
    let masses = family
        .propagators
        .iter()
        .enumerate()
        .filter_map(|(i, d)| (!d.constant.is_zero()).then_some(i))
        .collect::<Vec<_>>();
    if masses.len() != 1 {
        return Ok(None);
    }
    let a = masses[0];
    let b = (0..3).find(|&j| j != a).unwrap();
    let c = (0..3).find(|&j| j != a && j != b).unwrap();
    let branches = family
        .propagators
        .iter()
        .map(|d| crate::regions::branch(d, 2))
        .collect::<Result<Vec<_>>>()?;
    let mut routing = crate::algebra::inverse(&[branches[a].clone(), branches[b].clone()])?;
    let third = crate::algebra::matmul(&[branches[c].clone()], &routing);
    if third[0].iter().any(|v| v.is_zero()) {
        return Ok(None);
    }
    let ratio = (&third[0][0] / &third[0][1]).together().cancel();
    for row in &mut routing {
        row[1] = (&row[1] * &ratio).together().cancel();
    }
    let transformed = family.transform_loops(&routing)?;
    let scales = [
        transformed.propagators[a].scalar_products[0].clone(),
        transformed.propagators[b].scalar_products[2].clone(),
        transformed.propagators[c].scalar_products[0].clone(),
    ];
    let mass = (-&transformed.propagators[a].constant / &scales[0])
        .together()
        .cancel();
    let indices = [target.0[a] as u32, target.0[b] as u32, target.0[c] as u32];
    let mut value = vacuum::single_mass_sunset(indices, &p.eval(&mass, &params)?, &dimension, p)?;
    for (s, n) in scales.iter().zip(indices) {
        value = p.mul(&value, &p.powi(&p.eval(s, &params)?, -i64::from(n)));
    }
    let determinant = p.eval(&crate::algebra::determinant(routing), &params)?;
    let jacobian = ComplexFloat::new(p.norm(&determinant), p.real(0));
    Ok(Some(
        p.mul(&value, &p.pow(&jacobian, &p.rational(&dimension))),
    ))
}

impl BoundaryProvider for RecursiveBoundary<'_> {
    fn leading(
        &self,
        family: &IntegralFamily,
        basis: &[Integral],
        epsilon: &Rational,
        p: Precision,
    ) -> Result<Vec<LeadingBoundary>> {
        OneLoopBoundary.leading(family, basis, epsilon, p)
    }
    fn constants(
        &self,
        family: &IntegralFamily,
        basis: &[Integral],
        epsilon: &Rational,
        p: Precision,
        solutions: &FrobeniusBasis,
    ) -> Result<Vec<ComplexFloat>> {
        if family.loops.len() == 1
            && matches!(self.options.mass_mode, MassMode::All | MassMode::Auto)
            && family
                .propagators
                .iter()
                .all(|d| d.scalar_products.first() == Some(&Atom::num(1)))
        {
            return solutions.match_leading(&self.leading(family, basis, epsilon, p)?);
        }
        let (_, shifted) =
            family.deform(symbol!("symbolica_amflow::eta"), &self.options.mass_mode)?;
        let regions = crate::regions::enumerate_regions(family, 10000)?;
        let parameters =
            ahash::HashMap::from_iter([(Atom::var(family.epsilon), p.rational(epsilon))]);
        let powers = basis
            .iter()
            .map(|integral| {
                regions
                    .iter()
                    .map(|region| {
                        crate::regions::expand_region(family, integral, &shifted, region, 0)
                            .map(|leading| -leading.eta_power)
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let orders = solutions.region_orders(&powers, 32)?;
        self.context.emit(Progress::BoundaryPlan {
            basis_size: basis.len(),
            region_series: orders
                .iter()
                .flatten()
                .filter(|order| order.is_some())
                .count(),
            coefficients: orders
                .iter()
                .flatten()
                .flatten()
                .map(|order| order + 1)
                .sum(),
            max_half_order: orders
                .iter()
                .flatten()
                .flatten()
                .copied()
                .max()
                .unwrap_or(0),
        })?;
        let mut data = vec![Vec::new(); basis.len()];
        for (i, integral) in basis.iter().enumerate() {
            for (r, region) in regions.iter().enumerate() {
                self.context.cancellation.check()?;
                let mut coefficients = [Vec::new(), Vec::new()];
                if let Some(half_order) = orders[i][r] {
                    let expansion = crate::regions::expand_region(
                        family, integral, &shifted, region, half_order,
                    )?;
                    let determinant = p.eval(&expansion.jacobian_determinant, &parameters)?;
                    let jacobian = p.pow(
                        &ComplexFloat::new(p.norm(&determinant), p.real(0)),
                        &p.rational(
                            &(Rational::from(family.dimension) - epsilon * &Rational::from(2)),
                        ),
                    );
                    for (k, expression) in expansion.coefficients.iter().enumerate() {
                        let factors = crate::integrand::factor_region(
                            expression,
                            &expansion.coordinates,
                            family,
                            &region.hard,
                            10000,
                        )?;
                        let mut value = p.zero();
                        for term in factors {
                            // A scaleless factor annihilates the whole product.
                            if term
                                .factors
                                .iter()
                                .map(|f| scaleless(&f.family, &f.integral))
                                .collect::<Result<Vec<_>>>()?
                                .contains(&true)
                            {
                                continue;
                            }
                            let mut v = p.eval(&term.coefficient, &parameters)?;
                            for factor in term.factors {
                                v = p.mul(
                                    &v,
                                    &self.evaluate(&factor.family, &factor.integral, epsilon, p)?,
                                );
                            }
                            value = p.add(&value, &v);
                        }
                        coefficients[k % 2].push(p.mul(&value, &jacobian));
                    }
                }
                // Retain uncomputed regions (and parity classes) with empty
                // coefficient arrays. Matching must treat them as unknown,
                // rather than silently assuming their contributions vanish.
                for (parity, coefficients) in coefficients.into_iter().enumerate() {
                    data[i].push(RegionBoundary {
                        exponent: (&powers[i][r] + Atom::num((parity as i64, 2)))
                            .together()
                            .cancel(),
                        coefficients,
                    });
                }
            }
        }
        solutions.match_regions(&data)
    }
}
