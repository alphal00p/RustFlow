//! Exact specialization of bounded index-domain and denominator-ray searches.
use rustred::algebra::CoefficientPolynomial;
use rustred::solver::{
    CoordinateCase, Integral, IntegralOrder, RuleCandidate, RuleDispatchPolicy, SearchOptions,
    SectorEvent, SectorRule, SectorSolveError, SectorSolveOptions, SectorSolver, SolverError,
    SourceSystem, SourceVisitOrder, Term,
};
use std::collections::BTreeMap;
use symbolica::domains::rational_polynomial::FromNumeratorAndDenominator;
use symbolica::prelude::*;

pub(super) struct Applied<const N: usize> {
    pub terms: Vec<Term<N, RationalPolynomial<IntegerRing, u16>>>,
    pub conditions: Vec<CoefficientPolynomial>,
}

pub(super) struct Bank<'a, const N: usize> {
    sources: &'a SourceSystem<N>,
    order: IntegralOrder<N>,
    visit_order: Option<SourceVisitOrder>,
    generic: Option<Vec<SectorRule<N>>>,
    rays: BTreeMap<[Option<i16>; N], Vec<SectorRule<N>>>,
    pub domains: usize,
    pub domain_hits: usize,
    pub generated: usize,
    pub hits: usize,
    pub missed: usize,
}

fn generic_rules<const N: usize>(
    solver: &SectorSolver<'_, N>,
    sources: &SourceSystem<N>,
    depth: u32,
    max_cases: usize,
) -> Vec<SectorRule<N>> {
    let mut admitted = Vec::new();
    let result = solver.solve_domains_with_observer(
        vec![
            CoordinateCase::new(*sources.fixed())
                .expect("symbolic coordinates are valid")
                .into(),
        ],
        SectorSolveOptions {
            symbolic: SearchOptions {
                max_depth: Some(depth),
                ..Default::default()
            },
            numerical_depth: depth,
            max_symbolic_cases: Some(max_cases),
            ..Default::default()
        },
        |event| {
            if let SectorEvent::RuleFound { rule, .. } = event {
                // RustRed emits this only after exact guards and the complete
                // exceptional geometry have been admitted. A later budget
                // stop does not invalidate this equation on its guarded domain.
                admitted.push(SectorRule {
                    candidate: RuleCandidate {
                        case: rule.candidate.case.clone(),
                        target: rule.candidate.target,
                        rhs: rule.candidate.rhs.clone(),
                        sources: rule.candidate.sources.clone(),
                        stats: rule.candidate.stats,
                    },
                    exceptions: rule.exceptions.clone(),
                    dispatch_policy: rule.dispatch_policy,
                });
            }
        },
    );
    match result {
        Ok(_) | Err(SectorSolveError::CaseBudget { .. }) => admitted,
        // Other failures do not install a partially searched domain. Neither
        // successful nor exhausted searches promote finite residuals to rules.
        Err(_) => Vec::new(),
    }
}

fn specialize<const N: usize>(
    p: &CoefficientPolynomial,
    indices: &[usize; N],
    values: &[i16; N],
) -> CoefficientPolynomial {
    let mut p = p.clone();
    // Last-variable evaluation is a cheaper polynomial path in Symbolica.
    let mut replacements = indices
        .iter()
        .copied()
        .zip(values.iter().copied())
        .collect::<Vec<_>>();
    replacements.sort_unstable_by_key(|v| std::cmp::Reverse(v.0));
    for (index, value) in replacements {
        if p.degree(index) != 0 {
            p = p.replace(index, &Integer::from(value));
        }
    }
    p
}

fn apply<const N: usize>(
    rule: &SectorRule<N>,
    target: [i16; N],
    sources: &SourceSystem<N>,
    solver: &SectorSolver<'_, N>,
) -> Result<Option<Applied<N>>, SolverError> {
    if rule.dispatch_policy != RuleDispatchPolicy::Partition
        || rule.candidate.target != rule.candidate.case.integral()
        || !rule.candidate.case.is_in_sector(solver.ordering().sector())
    {
        return Err(SolverError::InvalidInput(
            "incompatible parametric rule".into(),
        ));
    }
    if target
        .iter()
        .zip(solver.ordering().sector())
        .any(|(&n, &active)| (n > 0) != active)
        || rule
            .candidate
            .case
            .fixed()
            .iter()
            .zip(target)
            .any(|(fixed, n)| fixed.is_some_and(|v| v != n))
    {
        return Ok(None);
    }
    let indices = sources.index_variables();
    let mut polynomials = sources.conditions().iter().collect::<Vec<_>>();
    if let Some(affine) = rule.candidate.case.affine() {
        polynomials.extend(affine.equations());
    }
    polynomials.extend(rule.exceptions.branches.iter().flatten());
    for term in &rule.candidate.rhs {
        polynomials.extend([&term.coefficient.numerator, &term.coefficient.denominator]);
    }
    if polynomials
        .iter()
        .any(|p| p.variables().as_slice() != sources.coefficient_variables())
    {
        return Err(SolverError::InvalidInput(
            "parametric polynomial variable map mismatch".into(),
        ));
    }
    let at = |p: &CoefficientPolynomial| specialize(p, indices, &target);
    if let Some(affine) = rule.candidate.case.affine() {
        if affine.index_variables() != indices {
            return Err(SolverError::InvalidInput(
                "parametric affine index map mismatch".into(),
            ));
        }
        if affine.equations().iter().any(|p| !at(p).is_zero()) {
            return Ok(None);
        }
    }
    let mut conditions = Vec::new();
    for p in sources.conditions() {
        let p = at(p);
        if p.is_zero() {
            return Ok(None);
        }
        if !p.is_constant() {
            conditions.push(p);
        }
    }
    // Each branch is an AND of zero equalities. Preserve a nonzero witness
    // for every excluded branch; splitting an AND into individual vetoes
    // would change the rule's domain.
    for branch in &rule.exceptions.branches {
        let Some(witness) = branch.iter().map(at).find(|p| !p.is_zero()) else {
            return Ok(None);
        };
        if !witness.is_constant() {
            conditions.push(witness);
        }
    }
    // Check the original denominators before cancellation or zero terms.
    let mut denominators = Vec::new();
    for term in &rule.candidate.rhs {
        if term
            .integral
            .powers()
            .iter()
            .zip(rule.candidate.target.powers())
            .any(|(a, b)| a.is_symbolic() != b.is_symbolic())
        {
            return Err(SolverError::InvalidInput(
                "parametric RHS coordinate layout mismatch".into(),
            ));
        }
        let denominator = at(&term.coefficient.denominator);
        if denominator.is_zero() {
            return Ok(None);
        }
        if !denominator.is_constant() {
            conditions.push(denominator.clone());
        }
        denominators.push(denominator);
    }
    let lhs = Integral::numeric(target)?;
    let mut terms = Vec::new();
    for (term, denominator) in rule.candidate.rhs.iter().zip(denominators) {
        let numerator = at(&term.coefficient.numerator);
        if numerator.is_zero() {
            continue;
        }
        let mut child = [0; N];
        for axis in 0..N {
            let shift = term.integral[axis].value() - rule.candidate.target[axis].value();
            child[axis] = target[axis]
                .checked_add(shift)
                .ok_or_else(|| SolverError::InvalidInput("parametric index overflow".into()))?;
        }
        let integral = Integral::numeric(child)?;
        if solver.ordering().compare(&lhs, &integral) != std::cmp::Ordering::Less {
            return Err(SolverError::InvalidInput(
                "nondecreasing specialized parametric rule".into(),
            ));
        }
        let coefficient = RationalPolynomial::from_num_den(numerator, denominator, &Z, true);
        terms.push(Term {
            integral,
            coefficient,
        });
    }
    Ok(Some(Applied { terms, conditions }))
}

impl<'a, const N: usize> Bank<'a, N> {
    pub fn new(sources: &'a SourceSystem<N>, solver: &SectorSolver<'_, N>) -> Self {
        Self {
            sources,
            order: solver.ordering().clone(),
            visit_order: solver.source_visit_order().cloned(),
            generic: None,
            rays: BTreeMap::new(),
            domains: 0,
            domain_hits: 0,
            generated: 0,
            hits: 0,
            missed: 0,
        }
    }
    pub fn try_reduce(
        &mut self,
        target: [i16; N],
        sources: &SourceSystem<N>,
        solver: &SectorSolver<'_, N>,
        depth: u32,
    ) -> Result<Option<Applied<N>>, SolverError> {
        if !std::ptr::eq(self.sources, sources)
            || &self.order != solver.ordering()
            || self.visit_order.as_ref() != solver.source_visit_order()
        {
            return Err(SolverError::InvalidInput(
                "parametric bank source/order mismatch".into(),
            ));
        }
        // Four-line paper sectors benefit from shared numerator domains as
        // well as shorter RHSs. Three-line measurements instead favored fixed
        // rays, so retain that strategy outside this measured tier.
        if self.generic.is_none() && solver.ordering().sector().iter().filter(|&&v| v).count() == 4
        {
            self.generic = Some(generic_rules(solver, sources, depth, 32));
            self.domains += 1;
        }
        if let Some(rules) = &self.generic {
            for rule in rules {
                if let Some(result) = apply(rule, target, sources, solver)? {
                    self.domain_hits += 1;
                    self.hits += 1;
                    return Ok(Some(result));
                }
            }
        }
        let ray = std::array::from_fn(|i| (!solver.ordering().sector()[i]).then_some(target[i]));
        if let std::collections::btree_map::Entry::Vacant(entry) = self.rays.entry(ray) {
            let result = solver.solve_domains(
                vec![CoordinateCase::new(ray)?.into()],
                SectorSolveOptions {
                    symbolic: SearchOptions {
                        max_depth: Some(depth),
                        ..Default::default()
                    },
                    numerical_depth: depth,
                    max_symbolic_cases: Some(32),
                    ..Default::default()
                },
            );
            self.generated += 1;
            // No finite search residual is installed as a master or rule.
            // A failed ray leaves the original concrete solver responsible.
            entry.insert(result.map_or_else(|_| Vec::new(), |r| r.rules));
        }
        for rule in &self.rays[&ray] {
            if let Some(result) = apply(rule, target, sources, solver)? {
                self.hits += 1;
                return Ok(Some(result));
            }
        }
        self.missed += 1;
        Ok(None)
    }
}
#[cfg(test)]
#[path = "parametric_tests.rs"]
mod tests;
