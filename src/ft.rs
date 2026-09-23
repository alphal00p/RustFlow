//! Feynman-parameter propagator combination and Euclidean series integration.
use crate::*;
use symbolica::prelude::*;

/// x D_left + (1-x) D_right, kept exact including its quadratic form.
pub fn combine(left: &Propagator, right: &Propagator, x: Symbol) -> Result<Propagator> {
    if left.scalar_products.len() != right.scalar_products.len() {
        return Err(Error::InvalidInput("propagator dimensions differ".into()));
    }
    let x = Atom::var(x);
    let y = Atom::num(1) - &x;
    Ok(Propagator {
        constant: (&x * &left.constant + &y * &right.constant).expand(),
        scalar_products: left
            .scalar_products
            .iter()
            .zip(&right.scalar_products)
            .map(|(a, b)| (&x * a + &y * b).expand())
            .collect(),
    })
}

/// Exact FT terminal reduction for a one-loop two-denominator family. The
/// parameter integral is propagated and integrated as power/log series.
pub fn evaluate_bubble(
    family: &IntegralFamily,
    target: &Integral,
    epsilon: &Rational,
    options: &FlowOptions,
    context: &RunContext,
) -> Result<ComplexFloat> {
    options.validate()?;
    family.validate()?;
    family.validate_integral(target)?;
    if family.loops.len() != 1
        || family.external.len() != 1
        || family.propagators.len() != 2
        || target.0.iter().any(|&n| n <= 0)
    {
        return Err(Error::Unsupported(
            "automatic FT currently supports positive-power one-loop bubbles".into(),
        ));
    }
    if family
        .propagators
        .iter()
        .any(|d| d.scalar_products[0] != Atom::num(1))
    {
        return Err(Error::Unsupported(
            "FT bubble needs unit loop-square coefficients".into(),
        ));
    }
    let p = Precision::decimal(options.digits + options.guard_digits)?;
    let params = ahash::HashMap::from_iter([(Atom::var(family.epsilon), p.rational(epsilon))]);
    let s = &family.external_gram[0][0];
    let sn = p.eval(s, &params)?;
    if sn.im != p.real(0) || sn.re > p.real(0) {
        return Err(Error::Unsupported(
            "FT automatic evaluation requires Euclidean kinematics".into(),
        ));
    }
    for d in &family.propagators {
        let mass = (&d.scalar_products[1].clone().pow(2) * s / Atom::num(4) - &d.constant)
            .together()
            .cancel();
        let m = p.eval(&mass, &params)?;
        if m.im != p.real(0) || m.re < p.real(0) {
            return Err(Error::Unsupported(
                "FT needs nonnegative real squared masses".into(),
            ));
        }
    }
    let x = symbol!("symbolica_amflow::ft_x");
    let z = Atom::var(x);
    let combined = combine(&family.propagators[0], &family.propagators[1], x)?;
    let mass = (combined.scalar_products[1].clone().pow(2) * s / Atom::num(4) - combined.constant)
        .together()
        .cancel();
    if mass.is_zero() {
        return Ok(p.zero());
    }
    let a = i64::from(target.0[0]);
    let b = i64::from(target.0[1]);
    let exponent = Atom::num((family.dimension, 2)) - Atom::var(family.epsilon) - Atom::num(a + b);
    let g = (Atom::num(a - 1) / &z - Atom::num(b - 1) / (Atom::num(1) - &z)
        + &exponent * mass.derivative(x) / &mass)
        .together()
        .cancel();
    let system = DifferentialSystem {
        variable: x,
        matrix: vec![vec![g, Atom::new()], vec![Atom::num(1), Atom::new()]],
    };
    let half = p.scale(&p.i(1), 1, 2);
    let mut at_half = params.clone();
    at_half.insert(z.clone(), half.clone());
    let dhalf = Rational::from((family.dimension, 2)) - epsilon;
    let gamma = p.div(
        &p.gamma_real(&p.rational(&(Rational::from(a + b) - dhalf)).re)?,
        &p.mul(&p.gamma_real(&p.real(a))?, &p.gamma_real(&p.real(b))?),
    );
    let mut value = p.mul(
        &gamma,
        &p.mul(
            &p.powi(&half, a + b - 2),
            &p.pow(&p.eval(&mass, &at_half)?, &p.eval(&exponent, &params)?),
        ),
    );
    if (a + b) % 2 != 0 {
        value = p.neg(&value);
    }
    let initial = BoundaryData {
        point: half,
        values: vec![value, p.zero()],
    };
    let compiled = system.compile(p, &params)?;
    let mut endpoints = Vec::new();
    for right in [false, true] {
        let anchor = p.i(i64::from(right));
        let distance = ComplexFloat::new(compiled.endpoint_radius(&anchor), p.real(0));
        let point = if right {
            p.sub(&anchor, &distance)
        } else {
            distance
        };
        let data = compiled.transport(&initial, std::slice::from_ref(&point), options, context)?;
        let endpoint = if right {
            let rules = std::collections::BTreeMap::from([(z.clone(), Atom::num(1) - &z)]);
            DifferentialSystem {
                variable: x,
                matrix: system
                    .matrix
                    .iter()
                    .map(|r| {
                        r.iter()
                            .map(|v| -crate::family::substitute(v, &rules))
                            .collect()
                    })
                    .collect(),
            }
        } else {
            system.clone()
        };
        let basis = endpoint.frobenius(p, &params, options.series_order)?;
        let distance = if right { p.sub(&p.i(1), &point) } else { point };
        let constants = basis.match_values(&distance, &data.values, &params)?;
        endpoints.push(crate::engine::project_limit(
            &basis,
            &constants,
            &[Atom::new(), Atom::num(1)],
            family.epsilon,
            x,
            &params,
        )?);
    }
    Ok(p.sub(&endpoints[1], &endpoints[0]))
}

use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

struct FtSystem {
    child: IntegralFamily,
    reduced: crate::reduction::ReducedSystem,
    augmented: DifferentialSystem,
    powers: [i16; 2],
}
#[derive(Default)]
struct FtMemo {
    systems: BTreeMap<String, Arc<FtSystem>>,
    values: BTreeMap<String, ComplexFloat>,
}

/// Recursive FT evaluator. Each propagator combination reduces the number of
/// active denominators in its midpoint boundary problem. The final Gaussian
/// integral is analytic, including polynomial tensor numerators.
pub struct FtEvaluator<'a> {
    backend: &'a dyn ReductionBackend,
    context: &'a RunContext,
    memo: Arc<Mutex<FtMemo>>,
    ancestry: Vec<String>,
}
impl<'a> FtEvaluator<'a> {
    pub fn new(backend: &'a dyn ReductionBackend, context: &'a RunContext) -> Self {
        Self {
            backend,
            context,
            memo: Arc::default(),
            ancestry: vec![],
        }
    }
    pub fn evaluate(
        &self,
        family: &IntegralFamily,
        integral: &Integral,
        epsilon: &Rational,
        options: &FlowOptions,
    ) -> Result<ComplexFloat> {
        options.validate()?;
        self.context.cancellation.check()?;
        family.validate_integral(integral)?;
        let p = Precision::decimal(options.digits + options.guard_digits)?;
        let key = format!(
            "{}:{family:?}:{:?}:{}:{}",
            family.convert()?.family.fingerprint(),
            integral.0,
            options.skip_reduction,
            options.refine_basis
        );
        let problem = family.integral_key(integral)?;
        let value_key = format!("{problem}:{epsilon}:{}:{}", p.bits, options.series_order);
        if let Some(value) = self
            .memo
            .lock()
            .map_err(|_| Error::Numerical("FT memo poisoned".into()))?
            .values
            .get(&value_key)
        {
            return Ok(value.clone());
        }
        if self.ancestry.contains(&problem) {
            return Err(Error::Limit("FT recursion cycle".into()));
        }
        if self.ancestry.len() >= 32 {
            return Err(Error::Limit("FT recursion exceeds 32 levels".into()));
        }
        let parameters =
            ahash::HashMap::from_iter([(Atom::var(family.epsilon), p.rational(epsilon))]);
        check_euclidean(family, integral, p, &parameters)?;
        let active = integral
            .0
            .iter()
            .enumerate()
            .filter_map(|(i, &n)| (n > 0).then_some(i))
            .collect::<Vec<_>>();
        let value = if active.is_empty() || crate::recursive::scaleless(family, integral)? {
            p.zero()
        } else if active.len() == 1 {
            crate::gaussian::terminal(family, integral, epsilon, p)?
        } else {
            let cached = {
                self.memo
                    .lock()
                    .map_err(|_| Error::Numerical("FT memo poisoned".into()))?
                    .systems
                    .get(&key)
                    .cloned()
            };
            let system = if let Some(system) = cached {
                system
            } else {
                let x = symbol!("symbolica_amflow::ft_parameter");
                // Prefer a pair with the same loop quadratic form.
                let mut pair = [active[0], active[1]];
                'pairs: for (position, &i) in active.iter().enumerate() {
                    for &j in &active[..position] {
                        let n = family.loops.len() * (family.loops.len() + 1) / 2;
                        if family.propagators[i].scalar_products[..n]
                            == family.propagators[j].scalar_products[..n]
                        {
                            pair = [j, i];
                            break 'pairs;
                        }
                    }
                }
                let [left, right] = pair;
                let mut combined = family.clone();
                combined.propagators[left] =
                    combine(&family.propagators[left], &family.propagators[right], x)?;
                let mut target = integral.clone();
                target.0[left] = target.0[left]
                    .checked_add(target.0[right])
                    .ok_or_else(|| Error::Limit("combined FT power overflow".into()))?;
                target.0[right] = 0;
                let cache = options
                    .cache_directory
                    .as_ref()
                    .map(|directory| -> Result<_> {
                        Ok((
                            directory,
                            format!(
                                "ft-{}",
                                crate::cache::system_key(
                                    &combined,
                                    std::slice::from_ref(&target),
                                    self.backend,
                                    options
                                )?
                            ),
                        ))
                    })
                    .transpose()?;
                let stored = if let Some((directory, key)) = &cache {
                    crate::cache::read_system(directory, key)?
                } else {
                    None
                };
                let reduced = if let Some((reduced, _)) = stored {
                    reduced
                } else {
                    let session = crate::reduction::ReductionSession::new(self.backend);
                    let build = if options.skip_reduction {
                        crate::reduction::parameter_differential_system_skip_initial
                    } else {
                        crate::reduction::parameter_differential_system
                    };
                    let mut reduced = build(&session, &combined, &[target], x, 12, self.context)?;
                    let refinement = if options.refine_basis {
                        Some(crate::refine::refine_basis(
                            &mut reduced,
                            x,
                            family.epsilon,
                            32,
                        )?)
                    } else {
                        None
                    };
                    if let Some((directory, key)) = cache {
                        crate::cache::write_system(directory, &key, &reduced, &refinement)?;
                    }
                    reduced
                };
                let z = Atom::var(x);
                let weight = z.clone().pow(i64::from(integral.0[left]) - 1)
                    * (Atom::num(1) - &z).pow(i64::from(integral.0[right]) - 1);
                let mut matrix = reduced
                    .matrix
                    .iter()
                    .map(|r| r.iter().cloned().chain([Atom::new()]).collect::<Vec<_>>())
                    .collect::<Vec<_>>();
                let last = reduced
                    .basis
                    .iter()
                    .map(|i| {
                        (&weight * reduced.targets[0].get(i).cloned().unwrap_or_default())
                            .together()
                            .cancel()
                    })
                    .chain([Atom::new()])
                    .collect();
                matrix.push(last);
                let child = combined;
                let system = Arc::new(FtSystem {
                    child,
                    reduced,
                    augmented: DifferentialSystem {
                        variable: x,
                        matrix,
                    },
                    powers: [integral.0[left], integral.0[right]],
                });
                self.memo
                    .lock()
                    .map_err(|_| Error::Numerical("FT memo poisoned".into()))?
                    .systems
                    .insert(key.clone(), system.clone());
                system
            };
            let mut ancestry = self.ancestry.clone();
            ancestry.push(problem);
            let child = Self {
                backend: self.backend,
                context: self.context,
                memo: self.memo.clone(),
                ancestry,
            };
            let x = system.augmented.variable;
            let z = Atom::var(x);
            let compiled = system.augmented.compile(p, &parameters)?;
            // Apparent singularities can lie at 1/2 even for a Euclidean
            // integral. Choose an exact interior boundary point satisfying
            // every inherited reduction condition before evaluating children.
            let mut reference = None;
            'denominators: for denominator in 2..=32i64 {
                'points: for numerator in 1..denominator {
                    let exact = Rational::from((numerator, denominator));
                    let value = p.rational(&exact);
                    let mut at_point = parameters.clone();
                    at_point.insert(z.clone(), value.clone());
                    for condition in &system.reduced.nonzero_conditions {
                        match p.eval(condition, &at_point) {
                            Ok(v) if p.finite(&v) && v != p.zero() => {}
                            _ => continue 'points,
                        }
                    }
                    if compiled
                        .poles
                        .iter()
                        .any(|pole| p.close(pole, &value, options.digits))
                    {
                        continue;
                    }
                    reference = Some((exact, value));
                    break 'denominators;
                }
            }
            let (exact_reference, reference) = reference.ok_or_else(|| {
                Error::Reduction("FT found no regular rational interior boundary point".into())
            })?;
            let boundary_family = system.child.at(&KinematicPoint(BTreeMap::from([(
                z.clone(),
                Atom::num(exact_reference),
            )])));
            let mut initial = system
                .reduced
                .basis
                .iter()
                .map(|i| child.evaluate(&boundary_family, i, epsilon, options))
                .collect::<Result<Vec<_>>>()?;
            initial.push(p.zero());
            let initial = BoundaryData {
                point: reference,
                values: initial,
            };
            let mut ends = Vec::new();
            for right in [false, true] {
                let anchor = p.i(i64::from(right));
                let distance = ComplexFloat::new(compiled.endpoint_radius(&anchor), p.real(0));
                let point = if right {
                    p.sub(&anchor, &distance)
                } else {
                    distance
                };
                let path = compiled.plan_path(&initial.point, &point, 1)?;
                let values = compiled.transport(&initial, &path, options, self.context)?;
                let endpoint = if right {
                    let rules = BTreeMap::from([(z.clone(), Atom::num(1) - &z)]);
                    DifferentialSystem {
                        variable: x,
                        matrix: system
                            .augmented
                            .matrix
                            .iter()
                            .map(|r| {
                                r.iter()
                                    .map(|a| -crate::family::substitute(a, &rules))
                                    .collect()
                            })
                            .collect(),
                    }
                } else {
                    system.augmented.clone()
                };
                let basis = endpoint.frobenius(p, &parameters, options.series_order)?;
                let distance = if right { p.sub(&p.i(1), &point) } else { point };
                let constants = basis.match_values(&distance, &values.values, &parameters)?;
                let mut weights = vec![Atom::new(); initial.values.len()];
                *weights.last_mut().unwrap() = Atom::num(1);
                ends.push(crate::engine::project_limit(
                    &basis,
                    &constants,
                    &weights,
                    family.epsilon,
                    x,
                    &parameters,
                )?);
            }
            let [a, b] = system.powers.map(i64::from);
            let prefactor = p.div(
                &p.gamma_real(&p.real(a + b))?,
                &p.mul(&p.gamma_real(&p.real(a))?, &p.gamma_real(&p.real(b))?),
            );
            p.mul(&prefactor, &p.sub(&ends[1], &ends[0]))
        };
        self.memo
            .lock()
            .map_err(|_| Error::Numerical("FT memo poisoned".into()))?
            .values
            .insert(value_key, value.clone());
        Ok(value)
    }
}

fn check_euclidean(
    family: &IntegralFamily,
    integral: &Integral,
    p: Precision,
    parameters: &ahash::HashMap<Atom, ComplexFloat>,
) -> Result<()> {
    let converted = family.convert()?;
    let symanzik = rustred::family::symanzik::SymanzikPolynomials::try_from_family_with_limits(
        &converted.family,
        Default::default(),
    )
    .map_err(|e| Error::InvalidInput(e.to_string()))?;
    for (polynomial, sign) in [(symanzik.u(), 1), (symanzik.f(), -1)] {
        for (coefficient, powers) in polynomial.terms() {
            if powers
                .iter()
                .zip(&integral.0)
                .any(|(&k, &n)| k > 0 && n <= 0)
            {
                continue;
            }
            let value = p.eval(
                &crate::family::substitute(&coefficient.to_expression(), &converted.reverse),
                parameters,
            )?;
            if value.im != p.real(0) || value.re * sign < p.real(0) {
                return Err(Error::Unsupported(
                    "FT requires nonnegative real U and Euclidean F coefficients in the active sector".into(),
                ));
            }
        }
    }
    Ok(())
}
