//! Report-only prototype: complete unit-call cache, never descendant substitution.
use rustred::algebra::{Coefficient, CoefficientContext, CoefficientPolynomial};
use rustred::solver::guarded::{
    GuardedApplicationFailure, GuardedApplicationStatus, GuardedProgram, GuardedReduction,
    GuardedReductionLimits, GuardedSource, GuardedSourceSystem, IndexBounds, IndexDomain,
    IndexRole,
};
use rustred::solver::{Integral, SearchOptions, SolverError, Term};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::sync::Arc;
use symbolica::prelude::*;

#[derive(Clone, Copy)]
struct CacheCaps {
    entries: usize,
    polynomial_terms: usize,
    coefficient_bits: u64,
}
#[derive(Default, Debug, Clone, Copy)]
struct Stats {
    native_calls: usize,
    hits: usize,
    promotions: usize,
    declined_results: usize,
    evictions: usize,
}
#[derive(Debug)]
enum Error {
    Native(SolverError),
    Transition(&'static str),
}
impl From<SolverError> for Error {
    fn from(e: SolverError) -> Self {
        Self::Native(e)
    }
}
#[derive(Clone)]
struct Entry<const N: usize> {
    target: [i64; N],
    limits: GuardedReductionLimits,
    result: GuardedReduction<N>,
    terms: usize,
    bits: u64,
}
// Program ownership binds source, ordered rules, maps, native order and old terminals.
// There is no API for inserting a result computed under a foreign program.
struct CachedProgram<const N: usize> {
    program: GuardedProgram<N>,
    entries: VecDeque<Entry<N>>,
    caps: CacheCaps,
    stats: Stats,
}
fn add_cost(p: &CoefficientPolynomial, terms: &mut usize, bits: &mut u64) -> Option<()> {
    *terms = terms.checked_add(p.nterms())?;
    for c in &p.coefficients {
        *bits = bits.checked_add(c.significant_bits())?;
    }
    Some(())
}
fn cost<const N: usize>(r: &GuardedReduction<N>) -> Option<(usize, u64)> {
    let (mut terms, mut bits) = (0usize, 0u64);
    for c in r
        .terms
        .values()
        .chain(r.unresolved.iter().map(|r| &r.coefficient))
    {
        add_cost(&c.numerator, &mut terms, &mut bits)?;
        add_cost(&c.denominator, &mut terms, &mut bits)?;
    }
    for p in &r.nonzero_conditions {
        add_cost(p, &mut terms, &mut bits)?;
    }
    Some((terms, bits))
}
fn same_limits(a: GuardedReductionLimits, b: GuardedReductionLimits) -> bool {
    a.max_pending_integrals == b.max_pending_integrals
        && a.max_rule_applications == b.max_rule_applications
}
impl<const N: usize> CachedProgram<N> {
    fn new(program: GuardedProgram<N>, caps: CacheCaps) -> Self {
        Self {
            program,
            entries: VecDeque::new(),
            caps,
            stats: Stats::default(),
        }
    }
    fn reset(self, program: GuardedProgram<N>) -> Self {
        Self::new(program, self.caps)
    }
    fn reduce(
        &mut self,
        target: [i64; N],
        limits: GuardedReductionLimits,
    ) -> Result<GuardedReduction<N>, Error> {
        if let Some(i) = self
            .entries
            .iter()
            .position(|e| e.target == target && same_limits(e.limits, limits))
        {
            let e = self.entries.remove(i).unwrap();
            let result = e.result.clone();
            self.entries.push_back(e);
            self.stats.hits += 1;
            return Ok(result);
        }
        self.stats.native_calls += 1;
        let result = self.program.reduce(target, limits)?;
        // Incomplete failures are never transformed, and no arithmetic error is cached.
        let eligible = result
            .unresolved
            .iter()
            .all(|r| r.reason == GuardedApplicationFailure::NoApplicableRule);
        let c = if eligible { cost(&result) } else { None };
        if let Some((terms, bits)) = c.filter(|(t, b)| {
            self.caps.entries > 0
                && *t <= self.caps.polynomial_terms
                && *b <= self.caps.coefficient_bits
        }) {
            while self.entries.len() >= self.caps.entries
                || self.entries.iter().map(|e| e.terms).sum::<usize>()
                    > self.caps.polynomial_terms - terms
                || self.entries.iter().map(|e| e.bits).sum::<u64>()
                    > self.caps.coefficient_bits - bits
            {
                assert!(self.entries.pop_front().is_some());
                self.stats.evictions += 1;
            }
            self.entries.push_back(Entry {
                target,
                limits,
                result: result.clone(),
                terms,
                bits,
            });
        } else {
            self.stats.declined_results += 1;
        }
        Ok(result)
    }
    fn promote(
        mut self,
        terminals: impl IntoIterator<Item = [i64; N]>,
        max_rules: usize,
    ) -> Result<Self, Error> {
        let terminals: BTreeSet<_> = terminals.into_iter().collect();
        if !self.program.terminals().is_subset(&terminals) {
            return Err(Error::Transition(
                "removing a prior terminal invalidates cached results",
            ));
        }
        let added = terminals
            .difference(self.program.terminals())
            .copied()
            .collect::<BTreeSet<_>>();
        for label in &added {
            let applied = self.program.apply(label)?;
            if !matches!(
                applied.status,
                GuardedApplicationStatus::Unresolved(GuardedApplicationFailure::NoApplicableRule)
            ) || !applied.terms.is_empty()
                || !applied.nonzero_conditions.is_empty()
            {
                return Err(Error::Transition(
                    "added terminal was not exactly NoApplicableRule",
                ));
            }
        }
        // This is deliberately the replayed API, not verified rebind: independent proof replay stays.
        self.program = self.program.with_terminals_replayed(terminals, max_rules)?;
        for entry in &mut self.entries {
            let mut remaining = Vec::new();
            for r in std::mem::take(&mut entry.result.unresolved) {
                assert_eq!(r.reason, GuardedApplicationFailure::NoApplicableRule);
                if added.contains(&r.integral) {
                    // Such labels were not old terminals; native uncovered output already has one exact sum/label.
                    assert!(
                        entry
                            .result
                            .terms
                            .insert(r.integral, r.coefficient)
                            .is_none()
                    );
                } else {
                    remaining.push(r);
                }
            }
            entry.result.unresolved = remaining;
        }
        self.stats.promotions += 1;
        Ok(self)
    }
}
fn assert_equal<const N: usize>(a: &GuardedReduction<N>, b: &GuardedReduction<N>) {
    assert_eq!(a.rule_applications, b.rule_applications);
    assert_eq!(a.terms, b.terms);
    assert_eq!(a.nonzero_conditions, b.nonzero_conditions);
    assert_eq!(a.unresolved.len(), b.unresolved.len());
    for (x, y) in a.unresolved.iter().zip(&b.unresolved) {
        assert_eq!(x.integral, y.integral);
        assert_eq!(x.coefficient, y.coefficient);
        assert_eq!(x.reason, y.reason);
    }
}
fn caps() -> CacheCaps {
    CacheCaps {
        entries: 32,
        polynomial_terms: 10000,
        coefficient_bits: 1_000_000,
    }
}
fn discover(
    sources: Arc<GuardedSourceSystem<1>>,
    points: impl IntoIterator<Item = i64>,
    terminals: impl IntoIterator<Item = [i64; 1]>,
) -> GuardedProgram<1> {
    let mut rules = Vec::new();
    for point in points {
        let found = sources
            .solve_domains(
                vec![IndexDomain::new([IndexBounds::fixed(point)]).unwrap()],
                SearchOptions {
                    max_depth: Some(0),
                    sample_seed: 0,
                    ..Default::default()
                },
                1,
            )
            .unwrap();
        assert!(
            !found.rules.is_empty(),
            "missing direct row at {point}: {:?}",
            found.unresolved
        );
        rules.extend(found.rules);
    }
    GuardedProgram::new(sources, rules, terminals).unwrap()
}

fn diamond_program(cancel: bool) -> GuardedProgram<1> {
    let context = CoefficientContext::try_new(["guarded_diamond_n", "guarded_diamond_x"]).unwrap();
    let x = context.parameter("guarded_diamond_x").unwrap();
    let target = Integral::symbolic([0]).unwrap();
    let branches = [
        (4, vec![(-1, context.one()), (-2, context.one())]),
        (2, vec![(-1, x.clone())]),
        (3, vec![(-2, if cancel { -x.clone() } else { x.clone() })]),
    ];
    let sources = Arc::new(
        GuardedSourceSystem::new(
            "replayed-uncovered-diamond",
            [IndexRole::Occupation],
            [0],
            branches
                .iter()
                .enumerate()
                .map(|(ordinal, (point, rhs))| {
                    let mut row = vec![Term {
                        integral: target,
                        coefficient: context.one().numerator,
                    }];
                    row.extend(rhs.iter().map(|(shift, coefficient)| Term {
                        integral: Integral::symbolic([*shift]).unwrap(),
                        coefficient: (-coefficient.clone()).numerator,
                    }));
                    GuardedSource::new(
                        format!("diamond-{ordinal}"),
                        row,
                        IndexDomain::new([IndexBounds::fixed(*point)]).unwrap(),
                    )
                    .with_nonzero_conditions(vec![x.numerator.clone()])
                })
                .collect(),
        )
        .unwrap(),
    );
    discover(sources, branches.iter().map(|(point, ..)| *point), [])
}

fn guarded_pending_diamond(cancel: bool) -> GuardedProgram<1> {
    let context = CoefficientContext::try_new(["pending_n", "pending_x", "pending_y"]).unwrap();
    let x = context.parameter("pending_x").unwrap();
    let y = context.parameter("pending_y").unwrap();
    let target = Integral::symbolic([0]).unwrap();
    let branches = [
        (
            6,
            context.one(),
            vec![(-1, context.one()), (-4, context.one())],
        ),
        (5, context.one(), vec![(-1, context.one())]),
        (
            4,
            context.one(),
            vec![(
                -2,
                if cancel {
                    -context.one()
                } else {
                    context.one()
                },
            )],
        ),
        (2, y.clone(), vec![(-1, &context.one() / &y)]),
    ];
    let sources = Arc::new(
        GuardedSourceSystem::new(
            "pending-native-order-guarded-diamond",
            [IndexRole::Occupation],
            [0],
            branches
                .iter()
                .enumerate()
                .map(|(ordinal, (point, pivot, rhs))| {
                    let mut row = vec![Term {
                        integral: target,
                        coefficient: pivot.numerator.clone(),
                    }];
                    row.extend(rhs.iter().map(|(shift, coefficient)| Term {
                        integral: Integral::symbolic([*shift]).unwrap(),
                        coefficient: (-(coefficient * pivot)).numerator,
                    }));
                    let source = GuardedSource::new(
                        format!("pending-{ordinal}"),
                        row,
                        IndexDomain::new([IndexBounds::fixed(*point)]).unwrap(),
                    );
                    if *point == 5 {
                        source.with_nonzero_conditions(vec![x.numerator.clone()])
                    } else {
                        source
                    }
                })
                .collect(),
        )
        .unwrap(),
    );
    discover(sources, branches.iter().map(|(point, ..)| *point), [])
}

fn capacity_fixture() -> GuardedProgram<1> {
    let context =
        CoefficientContext::try_new(["balanced_capacity_n", "balanced_capacity_x"]).unwrap();
    let target = Integral::symbolic([0]).unwrap();
    let branches = [
        (8, vec![(-1, 1), (-2, 1)]),
        (7, vec![(-6, 1)]),
        (6, vec![(-5, -1), (-4, 1), (-3, 1)]),
    ];
    let sources = Arc::new(
        GuardedSourceSystem::new(
            "balanced-budget-cancellation",
            [IndexRole::Occupation],
            [0],
            branches
                .iter()
                .enumerate()
                .map(|(ordinal, (point, rhs))| {
                    let mut row = vec![Term {
                        integral: target,
                        coefficient: context.one().numerator,
                    }];
                    row.extend(rhs.iter().map(|(shift, c)| Term {
                        integral: Integral::symbolic([*shift]).unwrap(),
                        coefficient: context.integer(-i64::from(*c)).numerator,
                    }));
                    GuardedSource::new(
                        format!("budget-{ordinal}"),
                        row,
                        IndexDomain::new([IndexBounds::fixed(*point)]).unwrap(),
                    )
                })
                .collect(),
        )
        .unwrap(),
    );
    discover(sources, branches.iter().map(|(point, _)| *point), [])
}
#[test]
fn completed_whole_calls_promote_only_no_rule_residuals() {
    let p = capacity_fixture();
    let mut cache = CachedProgram::new(p, caps());
    let limits = GuardedReductionLimits::default();
    for t in [[8], [7], [6]] {
        cache.reduce(t, limits).unwrap();
    }
    let old_rules = cache.program.rules().len();
    let mut cache = cache.promote([[1], [2]], old_rules).unwrap();
    for t in [[8], [7], [6]] {
        let expected = cache.program.reduce(t, limits).unwrap();
        let actual = cache.reduce(t, limits).unwrap();
        assert_equal(&actual, &expected);
    }
    assert_eq!(cache.stats.hits, 3);
    assert_eq!(cache.stats.native_calls, 3);
    assert!(
        cache
            .reduce([8], limits)
            .unwrap()
            .unresolved
            .iter()
            .any(|r| r.integral == [3])
    );
}
#[test]
fn cancellation_and_ordered_guards_match_after_replay_and_decode() {
    for cancel in [true, false] {
        let p = guarded_pending_diamond(cancel);
        let mut cache = CachedProgram::new(p, caps());
        let before = cache.reduce([6], Default::default()).unwrap();
        let mut cache = cache.promote([[1]], usize::MAX).unwrap();
        let after = cache.reduce([6], Default::default()).unwrap();
        assert_equal(
            &after,
            &cache.program.reduce([6], Default::default()).unwrap(),
        );
        assert_eq!(before.nonzero_conditions, after.nonzero_conditions);
        let ctx = CoefficientContext::try_new(["pending_n", "pending_x", "pending_y"]).unwrap();
        assert!(
            after
                .nonzero_conditions
                .contains(&ctx.parameter("pending_x").unwrap().numerator)
        );
        assert_eq!(
            after
                .nonzero_conditions
                .contains(&ctx.parameter("pending_y").unwrap().numerator),
            !cancel
        );
        let decoded = GuardedProgram::decode_generated(
            &cache.program.encode_native(Default::default()).unwrap(),
            cache.program.sources().clone(),
            Default::default(),
        )
        .unwrap();
        assert_equal(&after, &decoded.reduce([6], Default::default()).unwrap());
    }
}
#[test]
fn cancelling_no_rule_branches_do_not_create_terms_or_drop_used_guards() {
    let mut cache = CachedProgram::new(diamond_program(true), caps());
    let before = cache.reduce([4], Default::default()).unwrap();
    assert!(before.unresolved.is_empty());
    assert!(before.terms.is_empty());
    assert!(!before.nonzero_conditions.is_empty());
    let mut cache = cache.promote([[1]], usize::MAX).unwrap();
    let after = cache.reduce([4], Default::default()).unwrap();
    assert_equal(
        &after,
        &cache.program.reduce([4], Default::default()).unwrap(),
    );
    assert_equal(&after, &before);
}
#[test]
fn work_limit_and_failed_conditions_are_never_cached_or_promoted() {
    // At the exact application cap a newly added terminal can succeed where old NoRule failed.
    let p = diamond_program(false);
    let mut cache = CachedProgram::new(p, caps());
    let lim = GuardedReductionLimits {
        max_rule_applications: 3,
        max_pending_integrals: 100,
    };
    let old = cache.reduce([4], lim).unwrap();
    assert!(
        old.unresolved
            .iter()
            .any(|r| r.reason == GuardedApplicationFailure::WorkLimit)
    );
    assert!(cache.entries.is_empty());
    let mut cache = cache.promote([[1]], usize::MAX).unwrap();
    let new = cache.reduce([4], lim).unwrap();
    assert!(new.unresolved.is_empty());
    assert_eq!(cache.stats.hits, 0);
    assert_equal(&new, &cache.program.reduce([4], lim).unwrap());
    let c = CoefficientContext::try_new(["cache_guard_n", "cache_guard_m"]).unwrap();
    let n = c.parameter("cache_guard_n").unwrap();
    let m = c.parameter("cache_guard_m").unwrap();
    let d = IndexDomain::new([
        IndexBounds::new(Some(1), None).unwrap(),
        IndexBounds::new(Some(1), None).unwrap(),
    ])
    .unwrap();
    let sources = Arc::new(
        GuardedSourceSystem::new(
            "cache-coupled-guard",
            [IndexRole::Occupation; 2],
            [0, 1],
            vec![GuardedSource::new(
                "coupled",
                vec![
                    Term {
                        integral: Integral::symbolic([0, 0]).unwrap(),
                        coefficient: (&n - &m).numerator,
                    },
                    Term {
                        integral: Integral::symbolic([-1, 0]).unwrap(),
                        coefficient: (-c.one()).numerator,
                    },
                ],
                d.clone(),
            )],
        )
        .unwrap(),
    );
    let found = sources
        .solve_domains(
            vec![d],
            SearchOptions {
                max_depth: Some(1),
                sample_seed: 0,
                ..Default::default()
            },
            1,
        )
        .unwrap();
    let p = GuardedProgram::new(sources, found.rules, []).unwrap();
    let mut cache = CachedProgram::new(p, caps());
    for _ in 0..2 {
        let r = cache.reduce([2, 2], Default::default()).unwrap();
        assert!(matches!(
            r.unresolved[0].reason,
            GuardedApplicationFailure::ConditionVanished { .. }
        ));
    }
    assert_eq!(cache.stats.hits, 0);
    assert!(cache.promote([[2, 2]], usize::MAX).is_err());
}
#[test]
fn application_pending_limits_and_cache_cap_matrix_match_native() {
    for applications in 0..=5 {
        for pending in 0..=4 {
            for maximum_entries in [0, 1, 4] {
                let p = capacity_fixture();
                let lim = GuardedReductionLimits {
                    max_rule_applications: applications,
                    max_pending_integrals: pending,
                };
                let expected_old = p.reduce([8], lim).unwrap();
                let mut c = caps();
                c.entries = maximum_entries;
                let mut cache = CachedProgram::new(p, c);
                assert_equal(&cache.reduce([8], lim).unwrap(), &expected_old);
                let mut cache = cache.promote([[1], [2], [3]], usize::MAX).unwrap();
                let expected = cache.program.reduce([8], lim).unwrap();
                assert_equal(&cache.reduce([8], lim).unwrap(), &expected);
                assert_equal(&cache.reduce([8], lim).unwrap(), &expected);
            }
        }
    }
}
#[test]
fn terminal_removal_reducible_addition_and_rule_budget_are_rejected() {
    let p = capacity_fixture()
        .with_terminals_replayed([[1]], usize::MAX)
        .unwrap();
    assert!(
        CachedProgram::new(p, caps())
            .promote([], usize::MAX)
            .is_err()
    );
    assert!(
        CachedProgram::new(capacity_fixture(), caps())
            .promote([[6]], usize::MAX)
            .is_err()
    );
    assert!(
        CachedProgram::new(capacity_fixture(), caps())
            .promote([[-1]], usize::MAX)
            .is_err()
    );
    assert!(
        CachedProgram::new(capacity_fixture(), caps())
            .promote([[1]], 0)
            .is_err()
    );
    let p = capacity_fixture()
        .with_terminals_replayed([[1]], usize::MAX)
        .unwrap();
    let mut cache = CachedProgram::new(p, caps());
    cache.reduce([8], Default::default()).unwrap();
    let mut cache = cache.promote([[1], [2]], usize::MAX).unwrap();
    assert_equal(
        &cache.reduce([8], Default::default()).unwrap(),
        &cache.program.reduce([8], Default::default()).unwrap(),
    );
}
#[test]
fn immutable_owner_reset_limits_and_retained_shape_caps_invalidate_safely() {
    let mut c = caps();
    c.entries = 1;
    let mut cache = CachedProgram::new(diamond_program(false), c);
    cache.reduce([4], Default::default()).unwrap();
    cache.reduce([2], Default::default()).unwrap();
    assert_eq!(cache.entries.len(), 1);
    assert_eq!(cache.stats.evictions, 1);
    let mut altered = cache.reset(diamond_program(true));
    assert!(altered.entries.is_empty());
    let r = altered.reduce([4], Default::default()).unwrap();
    assert!(r.unresolved.is_empty());
    assert_eq!(altered.stats.hits, 0);
    let mut other_limits = GuardedReductionLimits::default();
    other_limits.max_pending_integrals -= 1;
    altered.reduce([4], other_limits).unwrap();
    assert_eq!(altered.stats.hits, 0);
    for cap in [
        CacheCaps {
            polynomial_terms: 0,
            ..caps()
        },
        CacheCaps {
            coefficient_bits: 0,
            ..caps()
        },
    ] {
        let mut cache = CachedProgram::new(diamond_program(false), cap);
        for _ in 0..2 {
            cache.reduce([4], Default::default()).unwrap();
        }
        assert_eq!(cache.stats.hits, 0);
        assert!(cache.entries.is_empty());
    }
}
#[test]
fn required_cut_zero_cannot_be_promoted_and_complete_zero_is_cacheable() {
    let c = CoefficientContext::try_new(["cache_cut_n", "cache_cut_h"]).unwrap();
    let sources = Arc::new(
        GuardedSourceSystem::new(
            "cache-cut-zero",
            [IndexRole::RequiredCut, IndexRole::Occupation],
            [0, 1],
            vec![GuardedSource::new(
                "zero",
                vec![Term {
                    integral: Integral::symbolic([0, 0]).unwrap(),
                    coefficient: c.one().numerator,
                }],
                IndexDomain::new([IndexBounds::fixed(1), IndexBounds::fixed(1)]).unwrap(),
            )],
        )
        .unwrap(),
    );
    let p = GuardedProgram::new(sources, vec![], []).unwrap();
    let mut cache = CachedProgram::new(p, caps());
    let z = cache.reduce([0, 1], Default::default()).unwrap();
    assert!(z.terms.is_empty() && z.unresolved.is_empty());
    assert_equal(&z, &cache.reduce([0, 1], Default::default()).unwrap());
    assert_eq!(cache.stats.hits, 1);
    assert!(cache.promote([[0, 1]], usize::MAX).is_err());
}
