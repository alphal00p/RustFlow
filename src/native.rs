//! Runtime orchestration of RustRed's native exact sector solver with a
//! factorized coefficient field. Source derivation, zero certificates, seed
//! search and exact replay remain entirely in RustRed.
use crate::family::{ConvertedFamily, substitute};
use crate::reduction::{LinearCombination, Reduction, RustRedBackend};
use crate::{Error, Integral, Progress, Result, RunContext};
use rustred::sector::{
    Mask,
    zero::{Analyzer, Decision},
};
use rustred::solver::{
    CoordinateCase, NumericalExactBackend, SearchOptions, SectorConfig, SectorSolver, SourceSystem,
};
use std::collections::{BTreeMap, BTreeSet};
use symbolica::prelude::*;

pub(crate) fn reduce(
    family: &ConvertedFamily,
    targets: &[Vec<i16>],
    options: &RustRedBackend,
    context: &RunContext,
) -> Result<Reduction> {
    macro_rules! dispatch { ($($n:literal),*)=>{match family.family.denominator_count() { $($n=>solve::<$n>(family,targets,options,context),)* n=>Err(Error::Unsupported(format!("native runtime supports 1..=12 denominators; received {n}"))) }}; }
    dispatch!(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
}
fn native_error(error: impl std::fmt::Display) -> Error {
    Error::Reduction(error.to_string())
}
fn solve<const N: usize>(
    family: &ConvertedFamily,
    targets: &[Vec<i16>],
    options: &RustRedBackend,
    context: &RunContext,
) -> Result<Reduction> {
    let mut pending = targets
        .iter()
        .map(|t| {
            t.as_slice()
                .try_into()
                .map_err(|_| Error::InvalidInput("native integral arity".into()))
        })
        .collect::<Result<BTreeSet<[i16; N]>>>()?;
    let initial = pending.iter().map(|p| p.map(|n| n > 0)).collect::<Vec<_>>();
    let sources =
        SourceSystem::<N>::from_family_with_lorentz(&family.family, options.include_lorentz)
            .map_err(native_error)?;
    let analyzer = Analyzer::try_unrestricted(&family.family).map_err(native_error)?;
    let mut zero = Vec::new();
    for bits in 0..(1usize << N) {
        context.cancellation.check()?;
        let sector = std::array::from_fn(|i| bits & (1 << i) != 0);
        if !initial
            .iter()
            .any(|parent| sector.iter().zip(parent).all(|(a, b)| !a || *b))
        {
            continue;
        }
        let mask = Mask::try_new(sector).map_err(native_error)?;
        if matches!(
            analyzer.analyze(&mask).map_err(native_error)?,
            Decision::ProvedZero(_)
        ) {
            zero.push(sector);
        }
    }
    let zero: std::sync::Arc<[[bool; N]]> = zero.into();
    let mut result = Reduction::default();
    result.nonzero_conditions.extend(
        sources
            .conditions()
            .iter()
            .map(|c| substitute(&c.to_expression(), &family.reverse)),
    );
    result.nonzero_conditions.extend(
        analyzer
            .domain()
            .conditions()
            .iter()
            .map(|c| substitute(&c.polynomial().to_expression(), &family.reverse)),
    );
    let mut visited = BTreeSet::new();
    let stage_key = blake3::hash(
        format!(
            "native-factorized-v3:{}:{:?}:{targets:?}:{}:{}:{}",
            family.family.fingerprint(),
            family.reverse,
            options.max_depth,
            options.include_lorentz,
            env!("RUSTRED_SOURCE_DIGEST")
        )
        .as_bytes(),
    )
    .to_hex()
    .to_string();
    if let Some(directory) = &options.checkpoints
        && let Some((restored, done, remaining)) =
            crate::cache::read_native_stage(directory, &stage_key)?
    {
        let array = |i: Integral| {
            i.0.as_slice()
                .try_into()
                .map_err(|_| Error::Cache("checkpoint index arity".into()))
        };
        visited = done
            .into_iter()
            .map(array)
            .collect::<Result<BTreeSet<[i16; N]>>>()?;
        pending = remaining
            .into_iter()
            .map(array)
            .collect::<Result<BTreeSet<[i16; N]>>>()?;
        result = restored;
    }
    let mut solvers = BTreeMap::new();
    while !pending.is_empty() {
        if visited.len() + pending.len() > options.max_targets {
            return Err(Error::Limit(format!(
                "Laporta RHS closure exceeds max_targets={}",
                options.max_targets
            )));
        }
        // Finish higher sectors before introducing their descendants into
        // lower-sector searches. Otherwise direct rules in a higher sector
        // can keep generating increasingly complicated lower-sector targets
        // after those lower sectors have already been expanded.
        let active_lines = pending
            .iter()
            .map(|i| i.iter().filter(|&&n| n > 0).count())
            .max()
            .unwrap();
        let selected = pending
            .iter()
            .filter(|i| i.iter().filter(|&&n| n > 0).count() == active_lines)
            .copied()
            .collect::<Vec<_>>();
        let mut sectors = BTreeMap::<[bool; N], Vec<CoordinateCase<N>>>::new();
        for target in selected {
            pending.remove(&target);
            visited.insert(target);
            sectors
                .entry(target.map(|n| n > 0))
                .or_default()
                .push(CoordinateCase::new(target.map(Some)).map_err(native_error)?);
        }
        for (sector, cases) in sectors {
            context.emit(Progress::SectorReduction {
                active_lines: sector.iter().filter(|&&v| v).count(),
                integrals: cases.len(),
                visited: visited.len(),
            })?;
            if zero.contains(&sector) {
                for case in cases {
                    result.rules.insert(
                        Integral(case.integral().powers().iter().map(|p| p.value()).collect()),
                        BTreeMap::new(),
                    );
                }
                continue;
            }
            if let std::collections::btree_map::Entry::Vacant(entry) = solvers.entry(sector) {
                entry.insert(
                    SectorSolver::new(
                        &sources,
                        sector,
                        SectorConfig {
                            zero_sectors: zero.clone(),
                            numerical_exact_backend: NumericalExactBackend::SparseFactorized,
                            ..Default::default()
                        },
                    )
                    .map_err(native_error)?,
                );
            }
            let solved = solvers[&sector]
                .solve_numeric_cases(
                    cases,
                    SearchOptions {
                        max_depth: Some(options.max_depth),
                        ..Default::default()
                    },
                )
                .map_err(native_error)?;
            context.cancellation.check()?;
            result.residuals.extend(
                solved
                    .residuals
                    .into_iter()
                    .map(|c| Integral(c.integral().powers().iter().map(|p| p.value()).collect())),
            );
            for rule in solved.rules {
                if rule
                    .target
                    .powers()
                    .iter()
                    .chain(rule.rhs.iter().flat_map(|t| t.integral.powers()))
                    .any(|p| p.is_symbolic())
                {
                    return Err(Error::Reduction(
                        "concrete search returned symbolic indices".into(),
                    ));
                }
                let target = Integral(rule.target.powers().iter().map(|p| p.value()).collect());
                let mut combination = LinearCombination::new();
                for term in rule.rhs {
                    let powers: [i16; N] =
                        std::array::from_fn(|i| term.integral.powers()[i].value());
                    if !visited.contains(&powers) {
                        pending.insert(powers);
                    }
                    let integral = Integral(powers.to_vec());
                    let coefficient =
                        substitute(&term.coefficient.to_expression(), &family.reverse);
                    let previous = combination.remove(&integral).unwrap_or_default();
                    combination.insert(integral, (previous + coefficient).together().cancel());
                    if !term.coefficient.denominator.is_constant() {
                        result.nonzero_conditions.push(substitute(
                            &term.coefficient.denominator.to_expression(),
                            &family.reverse,
                        ));
                    }
                }
                combination.retain(|_, c| !c.is_zero());
                result.rules.insert(target, combination);
            }
        }
        // Back-substitute into the requested targets before searching the
        // next frontier. Exact cancellations can remove intermediate RHS
        // integrals; searching those cancelled terms is unnecessary.
        context.emit(Progress::Substitution {
            rules: result.rules.len(),
        })?;
        let mut expander = Expander::new(&result, context)?;
        let mut frontier = BTreeSet::new();
        for target in targets {
            frontier.extend(expander.expand(&Integral(target.clone()))?.keys().cloned());
        }
        pending = frontier
            .iter()
            .map(|i| {
                i.0.as_slice()
                    .try_into()
                    .map_err(|_| Error::Reduction("native frontier arity".into()))
            })
            .collect::<Result<BTreeSet<[i16; N]>>>()?;
        pending.retain(|i| !visited.contains(i));
        context.emit(Progress::ReductionFrontier {
            searched: visited.len(),
            remaining: pending.len(),
        })?;
        result.residuals = frontier.into_iter().collect();
        if let Some(directory) = &options.checkpoints {
            result.residuals.sort();
            result.residuals.dedup();
            result.nonzero_conditions.sort();
            result.nonzero_conditions.dedup();
            crate::cache::write_native_stage(
                directory,
                &stage_key,
                &result,
                &visited
                    .iter()
                    .map(|i| Integral(i.to_vec()))
                    .collect::<Vec<_>>(),
                &pending
                    .iter()
                    .map(|i| Integral(i.to_vec()))
                    .collect::<Vec<_>>(),
            )?;
        }
    }
    result.residuals.sort();
    result.residuals.dedup();
    result.residuals.retain(|i| !result.rules.contains_key(i));
    result.nonzero_conditions.sort();
    result.nonzero_conditions.dedup();
    let mut expander = Expander::new(&result, context)?;
    let mut flattened = BTreeMap::new();
    let mut residuals = BTreeSet::new();
    for target in targets {
        let integral = Integral(target.clone());
        let terms = expander.expand(&integral)?;
        residuals.extend(terms.keys().cloned());
        if result.rules.contains_key(&integral) {
            flattened.insert(
                integral,
                terms
                    .iter()
                    .map(|(i, c)| (i.clone(), coefficient_atom(c)))
                    .collect(),
            );
        }
    }
    for i in &residuals {
        let powers: [i16; N] =
            i.0.as_slice()
                .try_into()
                .map_err(|_| Error::Reduction("native residual arity".into()))?;
        if !visited.contains(&powers) {
            return Err(Error::IncompleteReduction(
                "native reduction left an unsearched target contribution".into(),
            ));
        }
    }
    result.rules = flattened;
    result.residuals = residuals.into_iter().collect();
    for target in targets {
        result.expand(&Integral(target.clone()))?;
    }
    Ok(result)
}
use std::sync::Arc;
use symbolica::domains::factorized_rational_polynomial::{
    FactorizedRationalPolynomial, FactorizedRationalPolynomialField,
    FromNumeratorAndFactorizedDenominator,
};
use symbolica::poly::PolyVariable;
type Coefficient = FactorizedRationalPolynomial<IntegerRing, u16>;
type Terms = BTreeMap<Integral, Coefficient>;
struct Expander<'a> {
    reduction: &'a Reduction,
    context: &'a RunContext,
    variables: Arc<Vec<PolyVariable>>,
    field: FactorizedRationalPolynomialField<IntegerRing, u16>,
    coefficients: ahash::HashMap<Atom, Arc<Coefficient>>,
}
impl<'a> Expander<'a> {
    fn new(reduction: &'a Reduction, context: &'a RunContext) -> Result<Self> {
        let mut symbols = BTreeSet::new();
        for coefficient in reduction.rules.values().flat_map(|r| r.values()) {
            crate::family::scalar_symbols(coefficient.as_view(), &mut symbols)?;
        }
        let variables: Arc<Vec<PolyVariable>> = Arc::new(
            symbols
                .into_iter()
                .map(|a| match a.as_view() {
                    AtomView::Var(v) => PolyVariable::Symbol(v.get_symbol()),
                    _ => unreachable!(),
                })
                .collect(),
        );
        Ok(Self {
            reduction,
            context,
            field: FactorizedRationalPolynomialField::new(Z, variables.clone()),
            variables,
            coefficients: Default::default(),
        })
    }
    fn coefficient(&mut self, a: &Atom) -> Result<Arc<Coefficient>> {
        if let Some(c) = self.coefficients.get(a) {
            return Ok(c.clone());
        }
        let polynomial: RationalPolynomial<IntegerRing, u16> = a
            .try_to_rational_polynomial(&Q, &Z, Some(self.variables.clone()))
            .map_err(|e| Error::Reduction(e.to_string()))?;
        let c = Arc::new(Coefficient::from_num_den(
            polynomial.numerator,
            vec![(polynomial.denominator, 1)],
            &Z,
            true,
        ));
        self.coefficients.insert(a.clone(), c.clone());
        Ok(c)
    }
    fn expand(&mut self, integral: &Integral) -> Result<Terms> {
        // Accumulate coefficients in dependency order. Expanding every child
        // into its own complete residual map stores a quadratic number of
        // large rational functions on realistic Laporta DAGs. Forward
        // substitution needs only one coefficient per live integral.
        let mut state = BTreeMap::<Integral, u8>::new();
        let mut stack = vec![(integral.clone(), false)];
        let mut order = Vec::new();
        while let Some((node, finish)) = stack.pop() {
            self.context.cancellation.check()?;
            if finish {
                state.insert(node.clone(), 2);
                order.push(node);
                continue;
            }
            match state.get(&node) {
                Some(2) => continue,
                Some(1) => {
                    return Err(Error::IncompleteReduction("cyclic native reduction".into()));
                }
                _ => {}
            }
            state.insert(node.clone(), 1);
            stack.push((node.clone(), true));
            if let Some(terms) = self.reduction.rules.get(&node) {
                stack.extend(terms.keys().cloned().map(|child| (child, false)));
            }
        }
        let mut out = Terms::from([(integral.clone(), self.field.one())]);
        for node in order.into_iter().rev() {
            self.context.cancellation.check()?;
            let Some(terms) = self.reduction.rules.get(&node) else {
                continue;
            };
            let Some(value) = out.remove(&node) else {
                continue;
            };
            if value.is_zero() {
                continue;
            }
            for (child, coefficient) in terms {
                let coefficient = self.coefficient(coefficient)?;
                let product = self.field.mul(&value, coefficient.as_ref());
                match out.entry(child.clone()) {
                    std::collections::btree_map::Entry::Vacant(entry) => {
                        if !product.is_zero() {
                            entry.insert(product);
                        }
                    }
                    std::collections::btree_map::Entry::Occupied(mut entry) => {
                        self.field.add_assign(entry.get_mut(), &product);
                        if entry.get().is_zero() {
                            entry.remove();
                        }
                    }
                }
            }
        }
        Ok(out)
    }
}
fn coefficient_atom(c: &Coefficient) -> Atom {
    let mut a = c.numerator.to_expression() * Atom::num(c.numer_coeff.clone())
        / Atom::num(c.denom_coeff.clone());
    for (denominator, power) in &c.denominators {
        a *= denominator.to_expression().pow(-(*power as i64));
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn forward_substitution_cancels_shared_dependencies_and_rejects_cycles() {
        let i = |n| Integral(vec![n]);
        let mut reduction = Reduction {
            rules: BTreeMap::from([
                (
                    i(1),
                    BTreeMap::from([(i(2), parse!("x")), (i(3), Atom::num(1))]),
                ),
                (i(2), BTreeMap::from([(i(4), parse!("x+1"))])),
                (
                    i(3),
                    BTreeMap::from([(i(4), parse!("-x*(x+1)")), (i(5), Atom::num(1))]),
                ),
                (i(4), BTreeMap::from([(i(6), parse!("1/(x+1)"))])),
            ]),
            residuals: vec![i(5), i(6)],
            ..Default::default()
        };
        let context = RunContext::default();
        let actual = Expander::new(&reduction, &context)
            .unwrap()
            .expand(&i(1))
            .unwrap();
        let actual: LinearCombination = actual
            .iter()
            .map(|(i, c)| (i.clone(), coefficient_atom(c)))
            .collect();
        assert_eq!(actual, BTreeMap::from([(i(5), Atom::num(1))]));
        assert_eq!(actual, reduction.expand(&i(1)).unwrap());
        reduction
            .rules
            .insert(i(5), BTreeMap::from([(i(1), Atom::num(1))]));
        assert!(matches!(
            Expander::new(&reduction, &context).unwrap().expand(&i(1)),
            Err(Error::IncompleteReduction(_))
        ));
    }
}
