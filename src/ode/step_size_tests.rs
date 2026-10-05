use super::*;
use std::cell::{Cell, RefCell};

#[test]
fn analytic_pole_two_precisions_and_complex_waypoint_reversal() -> Result<()> {
    let x = symbol!("bracket_probe::x");
    for (working, guard, order) in [(45, 15, 24), (60, 30, 32)] {
        let p = Precision::decimal(working)?;
        let system = DifferentialSystem {
            variable: x,
            matrix: vec![vec![Atom::num(10) / (Atom::num(2) - Atom::var(x))]],
        }
        .compile(p, &Default::default())?;
        let end = p.scale(&p.i(19), 1, 10);
        let options = FlowOptions {
            step_size_strategy: StepSizeStrategy::Bracketed,
            digits: 20,
            guard_digits: guard,
            series_order: order,
            max_steps: 3000,
            ..Default::default()
        };
        let answer = system.transport(
            &BoundaryData {
                point: p.zero(),
                values: vec![p.i(1)],
            },
            &[p.complex(0, 1), p.complex(-1, 1), p.zero(), end.clone()],
            &options,
            &RunContext::default(),
        )?;
        assert_eq!(answer.point, end);
        assert!(p.close(&answer.values[0], &p.powi(&p.i(20), 10), 25));
        assert!(answer.diagnostics.superseded_successes > 0);
        assert_eq!(
            answer.diagnostics.predicate_evaluations,
            answer.diagnostics.steps
                + answer.diagnostics.rejected_steps
                + answer.diagnostics.superseded_successes
        );
    }
    Ok(())
}

#[test]
fn sparse_entire_defect_guard_remains_active() -> Result<()> {
    let p = Precision::decimal(45)?;
    let x = symbol!("bracket_sparse::x");
    let system = DifferentialSystem {
        variable: x,
        matrix: vec![vec![Atom::var(x).pow(20)]],
    }
    .compile(p, &Default::default())?;
    let answer = system.transport(
        &BoundaryData {
            point: p.zero(),
            values: vec![p.i(1)],
        },
        &[p.i(1)],
        &FlowOptions {
            step_size_strategy: StepSizeStrategy::Bracketed,
            digits: 15,
            guard_digits: 10,
            series_order: 16,
            max_steps: 3000,
            ..Default::default()
        },
        &RunContext::default(),
    )?;
    assert!(answer.diagnostics.rejected_steps > 0);
    let expected = p.scale(&p.i(1), 1, 21).exp();
    assert!(p.close(&answer.values[0], &expected, 15));
    Ok(())
}

struct Mock {
    p: Precision,
    calls: RefCell<Vec<C>>,
    reject_all: bool,
    cancel_at: Option<(usize, crate::CancellationToken)>,
    fatal_at: Option<usize>,
    charts: Cell<usize>,
}
impl Mock {
    fn new() -> Self {
        Self {
            p: Precision::decimal(40).unwrap(),
            calls: Default::default(),
            reject_all: false,
            cancel_at: None,
            fatal_at: None,
            charts: Cell::new(0),
        }
    }
}
impl SeriesSystem for Mock {
    fn whole_segment_residual(&self, _: &Self::Chart, _: &C) -> Result<Option<Vec<Float>>> {
        Ok(Some(vec![self.precision().real(0); self.dimension()]))
    }

    type State = C;
    type Chart = C;
    fn initial_state(&self, b: &BoundaryData) -> Result<C> {
        Ok(b.point.clone())
    }
    fn precision(&self) -> Precision {
        self.p
    }
    fn dimension(&self) -> usize {
        1
    }
    fn poles(&self) -> &[C] {
        &[]
    }
    fn local_chart(
        &self,
        center: &C,
        values: &[C],
        order: usize,
        state: &C,
    ) -> Result<(Vec<Vec<C>>, C)> {
        assert_eq!(state, center, "only committed state may seed next chart");
        self.charts.set(self.charts.get() + 1);
        let mut v = vec![vec![self.p.zero()]; order + 1];
        v[0] = values.to_vec();
        Ok((v, center.clone()))
    }
    fn rhs(&self, _: &C, _: &[C], _: &C) -> Result<Vec<C>> {
        Ok(vec![self.p.zero()])
    }
    fn accepted_state(&self, center: &C, point: &C, _: &Float) -> Result<Option<C>> {
        self.calls.borrow_mut().push(point.clone());
        let count = self.calls.borrow().len();
        if self.fatal_at == Some(count) {
            return Err(Error::Numerical("mock fatal predicate".into()));
        }
        if let Some((n, token)) = &self.cancel_at
            && count == *n
        {
            token.cancel();
        }
        let length = self.p.norm(&self.p.sub(point, center));
        // Explicitly nonmonotone: an additional distant island is acceptable.
        let accepted = length <= self.p.real(5) / 8
            || (length >= self.p.real(7) / 8 && length <= self.p.real(9) / 10);
        Ok((!self.reject_all && accepted).then(|| point.clone()))
    }
}
fn mock_run(
    system: &Mock,
    budget: usize,
    context: &RunContext,
    observed: &RefCell<Vec<C>>,
) -> Result<FlowResult> {
    let p = system.p;
    let (answer, _) = transport_series_observed(
        system,
        &BoundaryData {
            point: p.zero(),
            values: vec![p.i(1)],
        },
        &[p.i(1)],
        &FlowOptions {
            step_size_strategy: StepSizeStrategy::Bracketed,
            digits: 20,
            guard_digits: 10,
            series_order: 8,
            max_steps: budget,
            ..Default::default()
        },
        context,
        None,
        |_, state| observed.borrow_mut().push(state.clone()),
    )?;
    Ok(answer)
}
#[test]
fn nonmonotone_branch_predicate_commits_only_individually_valid_states() -> Result<()> {
    let system = Mock::new();
    let observed = RefCell::new(vec![]);
    let answer = mock_run(&system, 100, &RunContext::default(), &observed)?;
    assert_eq!(observed.borrow()[0], system.p.scale(&system.p.i(5), 1, 8));
    assert_eq!(observed.borrow().len(), answer.diagnostics.steps);
    assert_eq!(answer.point, system.p.i(1));
    assert_eq!(
        answer.diagnostics.predicate_evaluations,
        system.calls.borrow().len()
    );
    assert_eq!(
        answer.diagnostics.predicate_evaluations,
        answer.diagnostics.steps
            + answer.diagnostics.rejected_steps
            + answer.diagnostics.superseded_successes
    );
    assert_eq!(answer.diagnostics.superseded_successes, 1);
    Ok(())
}
#[test]
fn tight_budget_charges_every_trial_and_skips_optional_work() {
    for budget in [1, 2, 3, 4] {
        let system = Mock::new();
        let observed = RefCell::new(vec![]);
        assert!(matches!(
            mock_run(&system, budget, &RunContext::default(), &observed),
            Err(Error::Limit(_))
        ));
        assert_eq!(system.calls.borrow().len(), budget);
        assert_eq!(observed.borrow().len(), usize::from(budget > 1));
    }
}
#[test]
fn mandatory_local_failure_limit_is_unchanged() {
    let mut system = Mock::new();
    system.reject_all = true;
    assert!(matches!(
        mock_run(&system, 100, &RunContext::default(), &RefCell::new(vec![])),
        Err(Error::Accuracy(_))
    ));
    assert_eq!(system.calls.borrow().len(), 32);
}
#[test]
fn cancellation_before_or_during_optional_search_never_commits() {
    for cancel_at in [0, 3] {
        let context = RunContext::default();
        let mut system = Mock::new();
        if cancel_at == 0 {
            context.cancellation.cancel();
        } else {
            system.cancel_at = Some((cancel_at, context.cancellation.clone()));
        }
        let observed = RefCell::new(vec![]);
        assert!(matches!(
            mock_run(&system, 100, &context, &observed),
            Err(Error::Cancelled)
        ));
        assert_eq!(system.calls.borrow().len(), cancel_at);
        assert!(observed.borrow().is_empty());
    }
}
#[test]
fn fatal_optional_error_is_not_replaced_by_earlier_success() {
    let mut system = Mock::new();
    system.fatal_at = Some(3);
    assert!(matches!(
        mock_run(&system, 100, &RunContext::default(), &RefCell::new(vec![])),
        Err(Error::Numerical(_))
    ));
}

#[test]
fn default_halving_and_multiple_waypoint_budget_are_explicit() -> Result<()> {
    assert_eq!(
        FlowOptions::default().step_size_strategy,
        StepSizeStrategy::Halving
    );
    let system = Mock::new();
    let p = system.p;
    let options = FlowOptions {
        digits: 20,
        guard_digits: 10,
        series_order: 8,
        max_steps: 6,
        ..Default::default()
    };
    let boundary = BoundaryData {
        point: p.zero(),
        values: vec![p.i(1)],
    };
    let waypoints = [p.i(1), p.i(2)];
    let answer = transport_series(
        &system,
        &boundary,
        &waypoints,
        &options,
        &RunContext::default(),
        None,
    )?;
    assert_eq!(answer.point, p.i(2));
    assert_eq!(answer.diagnostics.predicate_evaluations, 6);
    assert_eq!(answer.diagnostics.superseded_successes, 0);
    let system = Mock::new();
    let observed = RefCell::new(vec![]);
    let options = FlowOptions {
        step_size_strategy: StepSizeStrategy::Bracketed,
        ..options
    };
    let error = transport_series_observed(
        &system,
        &boundary,
        &waypoints,
        &options,
        &RunContext::default(),
        None,
        |_, state| observed.borrow_mut().push(state.clone()),
    )
    .unwrap_err();
    assert!(matches!(error, Error::Limit(_)));
    assert_eq!(system.calls.borrow().len(), 6);
    assert_eq!(observed.borrow().last(), Some(&p.i(1)));
    Ok(())
}

#[test]
fn successful_initial_proposal_uses_no_upward_trials() -> Result<()> {
    let system = Mock::new();
    let p = system.p;
    let options = FlowOptions {
        digits: 20,
        guard_digits: 10,
        series_order: 8,
        max_steps: 1,
        step_size_strategy: StepSizeStrategy::Bracketed,
        ..Default::default()
    };
    let answer = transport_series(
        &system,
        &BoundaryData {
            point: p.zero(),
            values: vec![p.i(1)],
        },
        &[p.scale(&p.i(1), 1, 2)],
        &options,
        &RunContext::default(),
        None,
    )?;
    assert_eq!(answer.diagnostics.steps, 1);
    assert_eq!(answer.diagnostics.predicate_evaluations, 1);
    assert_eq!(answer.diagnostics.superseded_successes, 0);
    Ok(())
}
