//! Runtime orchestration of RustRed's native exact sector solver with a
//! factorized coefficient field. Source derivation, zero certificates, seed
//! search and exact replay remain entirely in RustRed.
use crate::family::{ConvertedFamily, substitute};
use crate::reduction::{LinearCombination, Reduction, RustRedBackend};
use crate::{Error, Integral, IntegralFamily, Progress, Result, RunContext};
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
    original: &IntegralFamily,
    dimension: &Atom,
    targets: &[Vec<i16>],
    options: &RustRedBackend,
    context: &RunContext,
) -> Result<Reduction> {
    macro_rules! dispatch { ($($n:literal),*)=>{match family.family.denominator_count() { $($n=>solve::<$n>(family,original,dimension,targets,options,context),)* n=>Err(Error::Unsupported(format!("native runtime supports 1..=12 denominators; received {n}"))) }}; }
    dispatch!(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
}
fn native_error(error: impl std::fmt::Display) -> Error {
    Error::Reduction(error.to_string())
}
fn save_stage<const N: usize>(
    options: &RustRedBackend,
    key: &str,
    reduction: &mut Reduction,
    visited: &BTreeSet<[i16; N]>,
    pending: &BTreeSet<[i16; N]>,
) -> Result<()> {
    if let Some(directory) = &options.checkpoints {
        reduction.residuals.sort();
        reduction.residuals.dedup();
        reduction.nonzero_conditions.sort();
        reduction.nonzero_conditions.dedup();
        crate::cache::write_native_stage(
            directory,
            key,
            reduction,
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
    Ok(())
}
fn save_if_due<const N: usize>(
    options: &RustRedBackend,
    key: &str,
    reduction: &mut Reduction,
    visited: &BTreeSet<[i16; N]>,
    pending: &BTreeSet<[i16; N]>,
    last_saved: &mut std::time::Instant,
) -> Result<()> {
    if options.checkpoints.is_some() && last_saved.elapsed() >= options.checkpoint_interval {
        save_stage(options, key, reduction, visited, pending)?;
        *last_saved = std::time::Instant::now();
    }
    Ok(())
}
fn solve<const N: usize>(
    family: &ConvertedFamily,
    original: &IntegralFamily,
    dimension: &Atom,
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
    let legacy_key_for_depth = |depth| {
        blake3::hash(
            format!(
                "native-factorized-v3:{}:{:?}:{targets:?}:{}:{}:{}",
                family.family.fingerprint(),
                family.reverse,
                depth,
                options.include_lorentz,
                env!("RUSTRED_SOURCE_DIGEST")
            )
            .as_bytes(),
        )
        .to_hex()
        .to_string()
    };
    let stage_key_for_depth = |depth| {
        let legacy_key = legacy_key_for_depth(depth);
        if options.bubble_subloops {
            blake3::hash(format!("bubble-subloops-v2:{legacy_key}").as_bytes())
                .to_hex()
                .to_string()
        } else {
            legacy_key
        }
    };
    let legacy_key = legacy_key_for_depth(options.max_depth);
    let stage_key = stage_key_for_depth(options.max_depth);
    let mut patterns = BTreeMap::<[bool; N], Option<crate::bubble::Bubble>>::new();
    let mut restart = if let Some(directory) = &options.checkpoints {
        crate::cache::read_native_stage(directory, &stage_key)?
    } else {
        None
    };
    // Exact identities survive a larger search radius, but a previous search
    // of a residual is not a search at the new depth. Reuse the raw DAG and
    // search all surviving leaves again, including formerly visited leaves.
    if restart.is_none()
        && let Some(depth) = options.max_depth.checked_sub(1)
        && let Some(directory) = &options.checkpoints
        && let Some((old, done, _)) =
            crate::cache::read_native_stage(directory, &stage_key_for_depth(depth))?
    {
        let frontier = structural_frontier(&old, targets, context)?;
        let done = done
            .into_iter()
            .filter(|i| old.rules.contains_key(i))
            .collect();
        restart = Some((old, done, frontier.into_iter().collect()));
    }
    // A legacy partial search can be reused only when none of its searched
    // sectors changes ordering under bubble elimination. In particular, the
    // expensive higher-sector searches need not be repeated on migration.
    if restart.is_none()
        && options.bubble_subloops
        && let Some(directory) = &options.checkpoints
        && let Some(old) = crate::cache::read_native_stage(directory, &legacy_key)?
    {
        let mut compatible = true;
        for integral in old.0.rules.keys().chain(&old.1) {
            let sector: [bool; N] = integral
                .0
                .iter()
                .map(|&n| n > 0)
                .collect::<Vec<_>>()
                .try_into()
                .map_err(|_| Error::Cache("checkpoint index arity".into()))?;
            if let std::collections::btree_map::Entry::Vacant(entry) = patterns.entry(sector) {
                entry.insert(crate::bubble::Bubble::find(
                    original,
                    &sector,
                    dimension.clone(),
                )?);
            }
            if patterns[&sector].is_some() {
                compatible = false;
                break;
            }
        }
        if compatible {
            restart = Some(old);
        }
    }
    if let Some((restored, done, remaining)) = restart {
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
    let mut last_saved = std::time::Instant::now();
    let search: Result<()> = (|| {
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
                context.cancellation.check()?;
                let sector = target.map(|n| n > 0);
                if zero.contains(&sector) {
                    result
                        .rules
                        .insert(Integral(target.to_vec()), BTreeMap::new());
                    pending.remove(&target);
                    visited.insert(target);
                    save_if_due(
                        options,
                        &stage_key,
                        &mut result,
                        &visited,
                        &pending,
                        &mut last_saved,
                    )?;
                    continue;
                }
                if options.bubble_subloops {
                    if let std::collections::btree_map::Entry::Vacant(entry) =
                        patterns.entry(sector)
                    {
                        entry.insert(crate::bubble::Bubble::find(
                            original,
                            &sector,
                            dimension.clone(),
                        )?);
                    }
                    if let Some(pattern) = &patterns[&sector]
                        && let Some(terms) = pattern.reduce(&Integral(target.to_vec()))?
                    {
                        let order = rustred::solver::IntegralOrder::new(sector, [false; N])
                            .with_permutation(
                                pattern
                                    .permutation()
                                    .try_into()
                                    .expect("bubble permutation arity"),
                            )
                            .map_err(native_error)?;
                        let lhs =
                            rustred::solver::Integral::numeric(target).map_err(native_error)?;
                        let decreasing = terms.keys().all(|i| {
                            let Ok(powers) = <[i16; N]>::try_from(i.0.as_slice()) else {
                                return false;
                            };
                            let Ok(rhs) = rustred::solver::Integral::numeric(powers) else {
                                return false;
                            };
                            order.compare(&lhs, &rhs) == std::cmp::Ordering::Less
                        });
                        if decreasing {
                            for (i, c) in &terms {
                                let powers: [i16; N] =
                                    i.0.as_slice().try_into().expect("bubble target arity");
                                if !visited.contains(&powers) {
                                    pending.insert(powers);
                                }
                                let rational: RationalPolynomial<IntegerRing, u16> = c
                                    .try_to_rational_polynomial(&Q, &Z, None)
                                    .map_err(native_error)?;
                                if !rational.denominator.is_constant() {
                                    result
                                        .nonzero_conditions
                                        .push(rational.denominator.to_expression());
                                }
                            }
                            result.rules.insert(Integral(target.to_vec()), terms);
                            pending.remove(&target);
                            visited.insert(target);
                            save_if_due(
                                options,
                                &stage_key,
                                &mut result,
                                &visited,
                                &pending,
                                &mut last_saved,
                            )?;
                            continue;
                        }
                    }
                }
                sectors
                    .entry(sector)
                    .or_default()
                    .push(CoordinateCase::new(target.map(Some)).map_err(native_error)?);
            }
            for (sector, mut cases) in sectors {
                if let std::collections::btree_map::Entry::Vacant(entry) = solvers.entry(sector) {
                    entry.insert(
                        SectorSolver::new(
                            &sources,
                            sector,
                            SectorConfig {
                                zero_sectors: zero.clone(),
                                permutation: patterns.get(&sector).and_then(|p| p.as_ref()).map(
                                    |p| {
                                        p.permutation()
                                            .try_into()
                                            .expect("bubble permutation arity")
                                    },
                                ),
                                numerical_exact_backend: NumericalExactBackend::SparseFactorized,
                                ..Default::default()
                            },
                        )
                        .map_err(native_error)?,
                    );
                }
                // Group nearby complexities so a very easy target does not
                // enlarge the exact elimination needed by every difficult target.
                let order = rustred::solver::IntegralOrder::new(sector, [false; N]);
                let order = if let Some(pattern) = patterns.get(&sector).and_then(|p| p.as_ref()) {
                    order
                        .with_permutation(
                            pattern
                                .permutation()
                                .try_into()
                                .expect("bubble permutation arity"),
                        )
                        .map_err(native_error)?
                } else {
                    order
                };
                cases.sort_by(|a, b| order.compare(&a.integral(), &b.integral()).reverse());
                for cases in cases.chunks(options.max_sector_batch) {
                    context.emit(Progress::SectorReduction {
                        active_lines: sector.iter().filter(|&&v| v).count(),
                        integrals: cases.len(),
                        visited: visited.len(),
                    })?;
                    // A progress callback may request cancellation itself.
                    context.cancellation.check()?;
                    let solved = solvers[&sector]
                        .solve_numeric_cases(
                            cases.to_vec(),
                            SearchOptions {
                                max_depth: Some(options.max_depth),
                                ..Default::default()
                            },
                        )
                        .map_err(native_error)?;
                    let completed = Progress::SectorReduced {
                        integrals: solved.stats.cases,
                        seeds: solved.stats.seeds,
                        rows: solved.stats.rows,
                        exact_trace_rows: solved.stats.exact_trace_rows,
                        elapsed_ms: solved.stats.elapsed.as_millis(),
                        exact_ms: solved.stats.exact_materialization.as_millis(),
                    };
                    result
                        .residuals
                        .extend(solved.residuals.into_iter().map(|c| {
                            Integral(c.integral().powers().iter().map(|p| p.value()).collect())
                        }));
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
                        let target =
                            Integral(rule.target.powers().iter().map(|p| p.value()).collect());
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
                            // Native coefficients are already reduced rational
                            // polynomials. Replacing their parameter names is
                            // bijective and needs no second polynomial gcd.
                            if let Some(previous) = combination.get_mut(&integral) {
                                *previous = (&*previous + coefficient).together().cancel();
                            } else {
                                combination.insert(integral, coefficient);
                            }
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
                    // Only completed searches enter the visited set. A checkpoint
                    // taken mid-round must leave every other case pending.
                    for case in cases {
                        let powers = std::array::from_fn(|i| case.integral().powers()[i].value());
                        pending.remove(&powers);
                        visited.insert(powers);
                    }
                    context.emit(completed)?;
                    save_if_due(
                        options,
                        &stage_key,
                        &mut result,
                        &visited,
                        &pending,
                        &mut last_saved,
                    )?;
                }
            }
            // Back-substitute into the requested targets before searching the
            // next frontier. Exact cancellations can remove intermediate RHS
            // integrals; searching those cancelled terms is unnecessary.
            context.emit(Progress::Substitution {
                rules: result.rules.len(),
            })?;
            let structural = structural_frontier(&result, targets, context)?;
            let unsearched = |i: &Integral| {
                <[i16; N]>::try_from(i.0.as_slice()).is_ok_and(|v| !visited.contains(&v))
            };
            let current_count = structural
                .iter()
                .filter(|i| line_count(i) == active_lines && unsearched(i))
                .count();
            let lower_count = structural
                .iter()
                .filter(|i| line_count(i) < active_lines && unsearched(i))
                .count();
            let mut frontier;
            if current_count > options.max_exact_frontier {
                // Extra searches are safe; unsearched leaves are never returned as
                // masters. Avoid constructing thousands of coefficients of lower
                // integrals that will subsequently be eliminated anyway.
                frontier = structural
                    .into_iter()
                    .filter(|i| line_count(i) <= active_lines)
                    .collect::<BTreeSet<_>>();
            } else if current_count == 0 && lower_count > options.max_exact_frontier {
                frontier = structural
                    .into_iter()
                    .filter(|i| line_count(i) < active_lines)
                    .collect();
            } else {
                let mut expander = Expander::new(&result, context)?;
                frontier = BTreeSet::new();
                for target in targets {
                    frontier.extend(
                        expander
                            .expand_filtered(&Integral(target.clone()), active_lines)?
                            .keys()
                            .cloned(),
                    );
                }
                let unfinished = frontier
                    .iter()
                    .any(|i| line_count(i) == active_lines && unsearched(i));
                if !unfinished {
                    if lower_count > options.max_exact_frontier {
                        frontier = structural
                            .into_iter()
                            .filter(|i| line_count(i) < active_lines)
                            .collect();
                    } else {
                        frontier.clear();
                        for target in targets {
                            frontier.extend(
                                expander.expand(&Integral(target.clone()))?.keys().cloned(),
                            );
                        }
                    }
                }
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
            save_stage(options, &stage_key, &mut result, &visited, &pending)?;
            last_saved = std::time::Instant::now();
        }
        Ok(())
    })();
    if let Err(error) = search {
        if !matches!(&error, Error::Io(_) | Error::Cache(_)) {
            save_stage(options, &stage_key, &mut result, &visited, &pending)?;
        }
        return Err(error);
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
fn line_count(integral: &Integral) -> usize {
    integral.0.iter().filter(|&&n| n > 0).count()
}
fn structural_frontier(
    reduction: &Reduction,
    targets: &[Vec<i16>],
    context: &RunContext,
) -> Result<BTreeSet<Integral>> {
    let mut seen = BTreeSet::new();
    let mut leaves = BTreeSet::new();
    let mut stack = targets.iter().cloned().map(Integral).collect::<Vec<_>>();
    while let Some(node) = stack.pop() {
        context.cancellation.check()?;
        if !seen.insert(node.clone()) {
            continue;
        }
        if let Some(terms) = reduction.rules.get(&node) {
            for child in terms.keys() {
                if line_count(child) > line_count(&node) {
                    return Err(Error::Reduction(
                        "native rewrite increases its sector; descending-sector search is invalid"
                            .into(),
                    ));
                }
                stack.push(child.clone());
            }
        } else {
            leaves.insert(node);
        }
    }
    Ok(leaves)
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
        self.expand_filtered(integral, 0)
    }
    fn expand_filtered(&mut self, integral: &Integral, minimum_lines: usize) -> Result<Terms> {
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
        // Retain exactly the dependency paths reaching requested terminal
        // sectors. This graph test is valid even for a caller's nontriangular
        // DAG; it does not assume that every rewrite lowers the sector.
        let mut retained = BTreeSet::new();
        for node in &order {
            let keep = if let Some(terms) = self.reduction.rules.get(node) {
                terms.keys().any(|child| retained.contains(child))
            } else {
                node.0.iter().filter(|&&n| n > 0).count() >= minimum_lines
            };
            if keep {
                retained.insert(node.clone());
            }
        }
        let mut out = Terms::new();
        if retained.contains(integral) {
            out.insert(integral.clone(), self.field.one());
        }
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
                if !retained.contains(child) {
                    continue;
                }
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
    #[test]
    fn filtered_substitution_preserves_all_paths_to_selected_sectors() {
        let target = Integral(vec![1, 1]);
        let low = Integral(vec![1, 0]);
        let high = Integral(vec![2, 1]);
        let child = Integral(vec![1, 2]);
        let reduction = Reduction {
            rules: BTreeMap::from([
                (
                    target.clone(),
                    BTreeMap::from([
                        (child.clone(), parse!("x")),
                        (low.clone(), parse!("1/(1+x^2)")),
                    ]),
                ),
                (
                    child,
                    BTreeMap::from([(high.clone(), parse!("1/x")), (low.clone(), parse!("x"))]),
                ),
            ]),
            residuals: vec![low.clone(), high.clone()],
            ..Default::default()
        };
        let context = RunContext::default();
        let mut expander = Expander::new(&reduction, &context).unwrap();
        let full = expander.expand(&target).unwrap();
        let filtered = expander.expand_filtered(&target, 2).unwrap();
        assert_eq!(filtered.len(), 1);
        assert_eq!(coefficient_atom(&filtered[&high]), Atom::num(1));
        assert_eq!(
            coefficient_atom(&filtered[&high]),
            coefficient_atom(&full[&high])
        );
        assert!(
            (coefficient_atom(&full[&low]) - parse!("x^2+1/(1+x^2)"))
                .together()
                .cancel()
                .is_zero()
        );
    }
}
