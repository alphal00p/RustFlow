use crate::family::{eta_derivative, substitute};
use crate::{Error, Integral, IntegralFamily, Progress, Result, RunContext};
use std::collections::{BTreeMap, BTreeSet};
use symbolica::prelude::*;

pub type LinearCombination = BTreeMap<Integral, Atom>;

#[derive(Clone, Debug, Default)]
pub struct Reduction {
    pub rules: BTreeMap<Integral, LinearCombination>,
    pub residuals: Vec<Integral>,
    pub nonzero_conditions: Vec<Atom>,
}

pub trait ReductionBackend: Send + Sync {
    fn identity(&self) -> String;
    fn reduce_at_epsilon(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        epsilon: &Rational,
        context: &RunContext,
    ) -> Result<Reduction> {
        let mut result = self.reduce(family, targets, context)?;
        let rules = BTreeMap::from([(Atom::var(family.epsilon), Atom::num(epsilon.clone()))]);
        for terms in result.rules.values_mut() {
            for coefficient in terms.values_mut() {
                *coefficient = substitute(coefficient, &rules).together().cancel();
            }
        }
        for condition in &mut result.nonzero_conditions {
            *condition = substitute(condition, &rules).together().cancel();
        }
        Ok(result)
    }

    fn reduce(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        context: &RunContext,
    ) -> Result<Reduction>;
}

#[derive(Clone, Debug)]
pub struct RustRedBackend {
    pub max_depth: u32,
    pub max_targets: usize,
    pub include_lorentz: bool,
    /// Above this frontier size, search a conservative dependency superset
    /// instead of materializing large temporary coefficients. Zero disables
    /// intermediate coefficient pruning; final reductions are always exact.
    pub max_exact_frontier: usize,
    /// Use RustRed's exact factorized rational-polynomial coefficient field.
    /// This does not enable its experimental reconstruction feature.
    pub factorized: bool,
    /// Optional restart checkpoints between native RHS-closure rounds.
    pub checkpoints: Option<std::path::PathBuf>,
}
impl Default for RustRedBackend {
    fn default() -> Self {
        Self {
            max_depth: 2,
            max_targets: 4096,
            include_lorentz: false,
            max_exact_frontier: 512,
            factorized: true,
            checkpoints: None,
        }
    }
}

impl ReductionBackend for RustRedBackend {
    fn identity(&self) -> String {
        format!(
            "rustred-exact:{}:{}:{}:{}:{}:{}",
            env!("RUSTRED_SOURCE_DIGEST"),
            self.max_depth,
            self.max_targets,
            self.include_lorentz,
            self.factorized,
            self.max_exact_frontier
        )
    }
    fn reduce_at_epsilon(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        epsilon: &Rational,
        context: &RunContext,
    ) -> Result<Reduction> {
        self.reduce_native(family, targets, Some(epsilon), context)
    }
    fn reduce(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        context: &RunContext,
    ) -> Result<Reduction> {
        self.reduce_native(family, targets, None, context)
    }
}
impl RustRedBackend {
    fn reduce_native(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        epsilon: Option<&Rational>,
        context: &RunContext,
    ) -> Result<Reduction> {
        context.emit(Progress::Reduction {
            integrals: targets.len(),
        })?;
        for i in targets {
            family.validate_integral(i)?;
        }
        let converted = family.convert_at_epsilon(epsilon)?;
        let targets = targets.iter().map(|i| i.0.clone()).collect::<Vec<_>>();
        if self.factorized {
            return crate::native::reduce(&converted, &targets, self, context);
        }
        let result = rustred::solver::bridge::solve_laporta(
            &converted.family,
            &targets,
            rustred::solver::bridge::DynamicSolveOptions {
                max_depth: self.max_depth,
                max_targets: self.max_targets,
                include_lorentz: self.include_lorentz,
            },
        )
        .map_err(|e| match e {
            rustred::solver::SolverError::SearchExhausted { .. } => Error::Limit(e.to_string()),
            rustred::solver::SolverError::InvalidInput(ref message)
                if message.starts_with("Laporta RHS closure exceeds") =>
            {
                Error::Limit(e.to_string())
            }
            rustred::solver::SolverError::InvalidInput(ref message)
                if message.starts_with("the runtime bridge supports") =>
            {
                Error::Unsupported(e.to_string())
            }
            _ => Error::Reduction(e.to_string()),
        })?;
        context.cancellation.check()?;
        let mut out = Reduction {
            residuals: result.residuals.into_iter().map(Integral).collect(),
            ..Default::default()
        };
        for rule in result.rules {
            if rule
                .target
                .iter()
                .chain(rule.rhs.iter().flat_map(|t| t.powers.iter()))
                .any(|p| p.symbolic)
            {
                return Err(Error::Reduction(
                    "concrete reduction returned symbolic indices".into(),
                ));
            }
            let target = Integral(rule.target.iter().map(|p| p.value).collect());
            let mut combination = LinearCombination::new();
            for term in rule.rhs {
                let coefficient = substitute(&term.coefficient.to_expression(), &converted.reverse);
                let key = Integral(term.powers.iter().map(|p| p.value).collect());
                let old = combination.remove(&key).unwrap_or_default();
                combination.insert(key, old + coefficient);
            }
            out.rules.insert(target, combination);
            for p in rule.nonzero_conditions {
                out.nonzero_conditions
                    .push(substitute(&p.to_expression(), &converted.reverse));
            }
        }
        for i in targets {
            let i = Integral(i);
            if !out.rules.contains_key(&i) && !out.residuals.contains(&i) {
                return Err(Error::Reduction("target omitted from reduction".into()));
            }
        }
        out.nonzero_conditions.sort();
        out.nonzero_conditions.dedup();
        Ok(out)
    }
}

/// Explicit reduction tables. Their algebraic validity is the caller's contract;
/// recursive substitutions and differential closure are checked by this crate.
#[derive(Clone, Debug)]
pub struct TableBackend {
    pub name: String,
    pub reduction: Reduction,
}
impl ReductionBackend for TableBackend {
    fn identity(&self) -> String {
        let content = format!("{:?}", self.reduction);
        format!(
            "table:{}:{}",
            self.name,
            blake3::hash(content.as_bytes()).to_hex()
        )
    }
    fn reduce(
        &self,
        _: &IntegralFamily,
        targets: &[Integral],
        context: &RunContext,
    ) -> Result<Reduction> {
        context.cancellation.check()?;
        for t in targets {
            if !self.reduction.rules.contains_key(t) && !self.reduction.residuals.contains(t) {
                return Err(Error::IncompleteReduction(format!(
                    "no table rule for {:?}",
                    t.0
                )));
            }
        }
        Ok(self.reduction.clone())
    }
}

impl Reduction {
    pub fn expand(&self, target: &Integral) -> Result<LinearCombination> {
        fn visit(
            r: &Reduction,
            i: &Integral,
            active: &mut BTreeSet<Integral>,
        ) -> Result<LinearCombination> {
            if !active.insert(i.clone()) {
                return Err(Error::IncompleteReduction("cyclic reduction table".into()));
            }
            let answer = if let Some(terms) = r.rules.get(i) {
                let mut out = LinearCombination::new();
                for (j, c) in terms {
                    for (k, d) in visit(r, j, active)? {
                        let old = out.remove(&k).unwrap_or_default();
                        out.insert(k, (old + c * &d).together().cancel());
                    }
                }
                out.retain(|_, v| !v.is_zero());
                out
            } else if r.residuals.contains(i) {
                BTreeMap::from([(i.clone(), Atom::num(1))])
            } else {
                return Err(Error::IncompleteReduction(format!(
                    "uncovered integral {:?}",
                    i.0
                )));
            };
            active.remove(i);
            Ok(answer)
        }
        visit(self, target, &mut BTreeSet::new())
    }
}

/// An exact basis change, with `previous_values = matrix * current_values`.
#[derive(Clone, Debug)]
pub struct BasisTransformation {
    pub previous: Vec<Integral>,
    pub current: Vec<Integral>,
    pub matrix: Vec<Vec<Atom>>,
}

#[derive(Clone, Debug)]
pub struct ReducedSystem {
    pub basis: Vec<Integral>,
    pub matrix: Vec<Vec<Atom>>,
    pub targets: Vec<LinearCombination>,
    pub nonzero_conditions: Vec<Atom>,
    pub candidates: BTreeMap<Integral, LinearCombination>,
    pub transformations: Vec<BasisTransformation>,
}

/// Build a derivative-closed spanning basis. This checks closure, not minimality.
pub fn differential_system(
    backend: &dyn ReductionBackend,
    family: &IntegralFamily,
    targets: &[Integral],
    shifted: &[bool],
    max_rounds: usize,
    context: &RunContext,
) -> Result<ReducedSystem> {
    build_differential_system(backend, family, targets, max_rounds, context, |i| {
        eta_derivative(i, shifted)
    })
}

pub fn parameter_derivative(
    family: &IntegralFamily,
    integral: &Integral,
    variable: Symbol,
) -> Result<Vec<(Integral, Atom)>> {
    family.validate_integral(integral)?;
    let inverse = crate::algebra::inverse(
        &family
            .propagators
            .iter()
            .map(|p| p.scalar_products.clone())
            .collect::<Vec<_>>(),
    )?;
    let mut terms = LinearCombination::new();
    for (i, (&power, propagator)) in integral.0.iter().zip(&family.propagators).enumerate() {
        if power == 0 {
            continue;
        }
        let coefficients = (0..family.propagators.len())
            .map(|k| {
                (0..inverse.len())
                    .fold(Atom::new(), |s, j| {
                        s + propagator.scalar_products[j].derivative(variable) * &inverse[j][k]
                    })
                    .together()
                    .cancel()
            })
            .collect::<Vec<_>>();
        let constant = (propagator.constant.derivative(variable)
            - coefficients
                .iter()
                .zip(&family.propagators)
                .fold(Atom::new(), |s, (a, d)| s + a * &d.constant))
        .together()
        .cancel();
        for (j, c) in coefficients.into_iter().map(Some).chain([None]).enumerate() {
            let coefficient = c.unwrap_or_else(|| constant.clone());
            if coefficient.is_zero() {
                continue;
            }
            let mut shifted = integral.clone();
            shifted.0[i] = shifted.0[i]
                .checked_add(1)
                .ok_or_else(|| Error::Limit("derivative index overflow".into()))?;
            if j < family.propagators.len() {
                shifted.0[j] = shifted.0[j]
                    .checked_sub(1)
                    .ok_or_else(|| Error::Limit("derivative index overflow".into()))?;
            }
            let old = terms.remove(&shifted).unwrap_or_default();
            terms.insert(
                shifted,
                (old - Atom::num(i64::from(power)) * coefficient)
                    .together()
                    .cancel(),
            );
        }
    }
    Ok(terms.into_iter().filter(|(_, c)| !c.is_zero()).collect())
}

pub fn parameter_differential_system(
    backend: &dyn ReductionBackend,
    family: &IntegralFamily,
    targets: &[Integral],
    variable: Symbol,
    max_rounds: usize,
    context: &RunContext,
) -> Result<ReducedSystem> {
    build_differential_system(backend, family, targets, max_rounds, context, |i| {
        parameter_derivative(family, i, variable)
    })
}

fn build_differential_system(
    backend: &dyn ReductionBackend,
    family: &IntegralFamily,
    targets: &[Integral],
    max_rounds: usize,
    context: &RunContext,
    derivative: impl Fn(&Integral) -> Result<Vec<(Integral, Atom)>>,
) -> Result<ReducedSystem> {
    let mut requested = targets.iter().cloned().collect::<BTreeSet<_>>();
    let mut previous = Vec::new();
    for _ in 0..max_rounds {
        context.cancellation.check()?;
        let reduction = backend.reduce(
            family,
            &requested.iter().cloned().collect::<Vec<_>>(),
            context,
        )?;
        let mut basis = BTreeSet::new();
        for i in &requested {
            basis.extend(reduction.expand(i)?.into_keys());
        }
        let basis = basis.into_iter().collect::<Vec<_>>();
        let mut derivatives = Vec::new();
        for i in &basis {
            derivatives.extend(derivative(i)?.into_iter().map(|(i, _)| i));
        }
        let derivatives_covered = derivatives.iter().all(|i| requested.contains(i));
        if derivatives_covered && basis == previous {
            let mut matrix = vec![vec![Atom::new(); basis.len()]; basis.len()];
            for (row, i) in basis.iter().enumerate() {
                for (j, c) in derivative(i)? {
                    for (k, d) in reduction.expand(&j)? {
                        let col = basis.binary_search(&k).map_err(|_| {
                            Error::IncompleteReduction("derivative introduced new residual".into())
                        })?;
                        matrix[row][col] = (&matrix[row][col] + &c * d).together().cancel();
                    }
                }
            }
            let candidates = requested
                .iter()
                .map(|i| Ok((i.clone(), reduction.expand(i)?)))
                .collect::<Result<_>>()?;
            return Ok(ReducedSystem {
                transformations: vec![],
                targets: targets
                    .iter()
                    .map(|i| reduction.expand(i))
                    .collect::<Result<_>>()?,
                basis,
                matrix,
                nonzero_conditions: reduction.nonzero_conditions,
                candidates,
            });
        }
        previous = basis;
        requested.extend(derivatives);
    }
    Err(Error::IncompleteReduction(format!(
        "no stable derivative-closed basis after {max_rounds} rounds"
    )))
}

/// Retain requested integrals as differential variables and reduce only their
/// derivatives. Redundant target directions are fixed by boundary matching.
pub fn differential_system_skip_initial(
    backend: &dyn ReductionBackend,
    family: &IntegralFamily,
    targets: &[Integral],
    shifted: &[bool],
    max_rounds: usize,
    context: &RunContext,
) -> Result<ReducedSystem> {
    build_retained_system(backend, family, targets, max_rounds, context, |i| {
        eta_derivative(i, shifted)
    })
}

pub fn parameter_differential_system_skip_initial(
    backend: &dyn ReductionBackend,
    family: &IntegralFamily,
    targets: &[Integral],
    variable: Symbol,
    max_rounds: usize,
    context: &RunContext,
) -> Result<ReducedSystem> {
    build_retained_system(backend, family, targets, max_rounds, context, |i| {
        parameter_derivative(family, i, variable)
    })
}

fn build_retained_system(
    backend: &dyn ReductionBackend,
    family: &IntegralFamily,
    targets: &[Integral],
    max_rounds: usize,
    context: &RunContext,
    derivative: impl Fn(&Integral) -> Result<Vec<(Integral, Atom)>>,
) -> Result<ReducedSystem> {
    let mut basis = targets.iter().cloned().collect::<BTreeSet<_>>();
    for _ in 0..max_rounds {
        context.cancellation.check()?;
        let derivatives = basis
            .iter()
            .map(&derivative)
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .flatten()
            .map(|(i, _)| i)
            .collect::<BTreeSet<_>>();
        let reduction = backend.reduce(
            family,
            &derivatives.iter().cloned().collect::<Vec<_>>(),
            context,
        )?;
        let expansions = derivatives
            .iter()
            .map(|i| Ok((i.clone(), reduction.expand(i)?)))
            .collect::<Result<BTreeMap<_, _>>>()?;
        let old = basis.len();
        for terms in expansions.values() {
            basis.extend(terms.keys().cloned());
        }
        if basis.len() != old {
            continue;
        }
        let basis = basis.into_iter().collect::<Vec<_>>();
        let mut matrix = vec![vec![Atom::new(); basis.len()]; basis.len()];
        for (i, integral) in basis.iter().enumerate() {
            for (derivative, c) in derivative(integral)? {
                for (j, a) in &expansions[&derivative] {
                    let col = basis
                        .binary_search(j)
                        .map_err(|_| Error::IncompleteReduction("unclosed derivative".into()))?;
                    matrix[i][col] = (&matrix[i][col] + &c * a).together().cancel();
                }
            }
        }
        return Ok(ReducedSystem {
            transformations: vec![],
            basis,
            matrix,
            targets: targets
                .iter()
                .map(|i| BTreeMap::from([(i.clone(), Atom::num(1))]))
                .collect(),
            nonzero_conditions: reduction.nonzero_conditions,
            candidates: expansions,
        });
    }
    Err(Error::IncompleteReduction(
        "retained-target system failed to close".into(),
    ))
}

/// Exact rational epsilon specialization, before native IBP elimination.
pub struct SampledBackend<'a> {
    pub backend: &'a dyn ReductionBackend,
    pub epsilon: Rational,
}
impl ReductionBackend for SampledBackend<'_> {
    fn identity(&self) -> String {
        format!("{}:epsilon={}", self.backend.identity(), self.epsilon)
    }
    fn reduce(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        context: &RunContext,
    ) -> Result<Reduction> {
        self.backend
            .reduce_at_epsilon(family, targets, &self.epsilon, context)
    }
}

/// Reuse exact rules between differential-closure rounds. Unsolved residuals are
/// deliberately searched again when requested; they are not cached as masters.
pub struct ReductionSession<'a> {
    backend: &'a dyn ReductionBackend,
    families: std::sync::Mutex<BTreeMap<String, Reduction>>,
}
impl<'a> ReductionSession<'a> {
    pub fn new(backend: &'a dyn ReductionBackend) -> Self {
        Self {
            backend,
            families: Default::default(),
        }
    }
}
impl ReductionBackend for ReductionSession<'_> {
    fn reduce_at_epsilon(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        epsilon: &Rational,
        context: &RunContext,
    ) -> Result<Reduction> {
        self.backend
            .reduce_at_epsilon(family, targets, epsilon, context)
    }
    fn identity(&self) -> String {
        self.backend.identity()
    }
    fn reduce(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        context: &RunContext,
    ) -> Result<Reduction> {
        context.cancellation.check()?;
        let key = format!("{}:{family:?}", family.convert()?.family.fingerprint());
        let mut cached = self
            .families
            .lock()
            .map_err(|_| Error::Reduction("reduction session poisoned".into()))?
            .get(&key)
            .cloned()
            .unwrap_or_default();
        let missing = targets
            .iter()
            .filter(|i| !cached.rules.contains_key(*i))
            .cloned()
            .collect::<Vec<_>>();
        if !missing.is_empty() {
            let new = self.backend.reduce(family, &missing, context)?;
            cached.rules.extend(new.rules);
            cached.residuals.extend(new.residuals);
            cached.nonzero_conditions.extend(new.nonzero_conditions);
            cached.residuals.retain(|i| !cached.rules.contains_key(i));
            cached.residuals.sort();
            cached.residuals.dedup();
            cached.nonzero_conditions.sort();
            cached.nonzero_conditions.dedup();
            // Merged tables must remain acyclic and cover every requested RHS.
            for target in targets {
                cached.expand(target)?;
            }
            self.families
                .lock()
                .map_err(|_| Error::Reduction("reduction session poisoned".into()))?
                .insert(key, cached.clone());
        }
        Ok(cached)
    }
}
