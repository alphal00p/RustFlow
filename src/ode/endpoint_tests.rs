//! Geometry invariants at the shared acceptance predicate. This analytic
//! fixture has Y=[C*h+C^2*h^2/2,1], h=x-1, C=2^400. Its point-dependent RHS
//! makes evaluating a different physical point an observable error.
use super::*;
use std::cell::RefCell;

struct EndpointProbe {
    p: Precision,
    center: C,
    large: C,
    rhs_points: RefCell<Vec<C>>,
    state_points: RefCell<Vec<C>>,
}

impl EndpointProbe {
    fn coefficients(&self, values: &[C], order: usize, jacobian: &C) -> Vec<Vec<C>> {
        let p = self.p;
        let mut coefficients = vec![vec![p.zero(); 2]; order + 1];
        coefficients[0] = values.to_vec();
        let slope = p.mul(&self.large, jacobian);
        coefficients[1][0] = slope.clone();
        coefficients[2][0] = p.scale(&p.mul(&slope, &slope), 1, 2);
        coefficients
    }
}

impl SeriesSystem for EndpointProbe {
    type State = C;
    type Chart = TaylorCoordinate;

    fn initial_state(&self, boundary: &BoundaryData) -> Result<C> {
        Ok(boundary.point.clone())
    }
    fn precision(&self) -> Precision {
        self.p
    }
    fn dimension(&self) -> usize {
        2
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
    ) -> Result<(Vec<Vec<C>>, Self::Chart)> {
        assert_eq!(center, &self.center);
        assert_eq!(state, center);
        Ok((
            self.coefficients(values, order, &self.p.i(1)),
            TaylorCoordinate::Identity,
        ))
    }
    fn mapped_chart(
        &self,
        coordinate: &TaylorCoordinate,
        center: &C,
        values: &[C],
        order: usize,
        state: &C,
    ) -> Result<MappedTaylorChart<Self::Chart>> {
        assert_eq!(center, &self.center);
        assert_eq!(state, center);
        // This fixture selects the linear member of the Mobius family, so its
        // transformed analytic polynomial remains exactly quadratic.
        let variable = symbol!("declared_endpoint_probe::local");
        assert!(
            coordinate
                .expression(variable)?
                .derivative(variable)
                .derivative(variable)
                .together()
                .cancel()
                .is_zero()
        );
        let jacobian = coordinate.jacobian_at(self.p, center, center)?;
        Ok((
            self.coefficients(values, order, &jacobian),
            coordinate.clone(),
            Vec::new(),
        ))
    }
    fn whole_segment_residual(&self, _: &Self::Chart, _: &C) -> Result<Option<Vec<Float>>> {
        // Native exact algebra below proves this fixture's quadratic solves
        // the stated ODE. The test isolates endpoint geometry, not residuals.
        Ok(Some(vec![self.p.real(0); 2]))
    }
    fn rhs(&self, point: &C, values: &[C], _: &Self::Chart) -> Result<Vec<C>> {
        self.rhs_points.borrow_mut().push(point.clone());
        let p = self.p;
        let offset = p.sub(point, &self.center);
        let coefficient = p.mul(&self.large, &p.add(&p.i(1), &p.mul(&self.large, &offset)));
        Ok(vec![p.mul(&coefficient, &values[1]), p.zero()])
    }
    fn accepted_state(&self, _: &Self::Chart, point: &C, _: &Float) -> Result<Option<C>> {
        self.state_points.borrow_mut().push(point.clone());
        Ok(Some(point.clone()))
    }
}

fn fixture() -> Result<(EndpointProbe, BoundaryData, C, C)> {
    let p = Precision::decimal(60)?;
    let high = Precision::decimal(160)?;
    let exact_large = Rational::from(Integer::from(2).pow(400));
    let step = high.rational(&(Rational::one() / &exact_large));
    let center = high.i(1);
    let target = high.add(&center, &step);
    let midpoint = high.add(&center, &high.scale(&step, 1, 2));
    assert_ne!(center, target);
    assert_ne!(center, midpoint);
    assert_eq!(p.add(&center, &p.sub(&target, &center)), p.i(1));
    let x = symbol!("declared_endpoint_probe::x");
    let offset = Atom::var(x) - Atom::one();
    let scaled = Atom::num(exact_large.clone()) * &offset;
    let exact_solution = &scaled + scaled.clone().pow(2) / Atom::num(2);
    let exact_rhs = Atom::num(exact_large.clone()) * (Atom::one() + scaled);
    assert!(
        (exact_solution.derivative(x) - exact_rhs)
            .expand()
            .is_zero()
    );
    let probe = EndpointProbe {
        p,
        center: center.clone(),
        large: p.rational(&exact_large),
        rhs_points: RefCell::new(Vec::new()),
        state_points: RefCell::new(Vec::new()),
    };
    let boundary = BoundaryData {
        point: center,
        values: vec![p.zero(), p.i(1)],
    };
    Ok((probe, boundary, target, midpoint))
}

#[test]
fn high_storage_target_matches_rhs_state_and_saved_polynomial() -> Result<()> {
    let (probe, boundary, target, midpoint) = fixture()?;
    let p = probe.p;
    let mut saved = Vec::new();
    let observed = RefCell::new(Vec::new());
    let (result, state) = transport_series_observed(
        &probe,
        &boundary,
        std::slice::from_ref(&target),
        &FlowOptions {
            digits: 20,
            guard_digits: 40,
            series_order: 16,
            max_steps: 8,
            ..Default::default()
        },
        &RunContext::default(),
        Some(&mut saved),
        |_, state| observed.borrow_mut().push(state.clone()),
    )?;
    assert_eq!(result.point, target);
    assert_eq!(state, target);
    assert_eq!(
        result.values,
        vec![p.rational(&Rational::from((3, 2))), p.i(1)]
    );
    assert_eq!(saved.len(), 1);
    assert_eq!(saved[0].center, boundary.point);
    assert_eq!(saved[0].end, target);
    assert_eq!(*probe.rhs_points.borrow(), vec![target.clone(), midpoint]);
    assert_eq!(*probe.state_points.borrow(), vec![target.clone()]);
    assert_eq!(*observed.borrow(), vec![target]);
    let local = saved[0]
        .coordinate
        .local_point(p, &saved[0].center, &saved[0].end)?;
    assert_eq!(
        evaluate_taylor(p, &saved[0].coefficients, &local).0,
        result.values
    );
    Ok(())
}

#[test]
fn mapped_trial_uses_declared_point_for_local_argument_jacobian_and_state() -> Result<()> {
    let (probe, boundary, target, midpoint) = fixture()?;
    let p = probe.p;
    let high = Precision::decimal(160)?;
    let delta = p.sub(&target, &boundary.point);
    // Symmetric geometry chooses b=0 and a non-unit physical Jacobian. These
    // points select the chart directly; the analytic fixture itself is entire.
    let twice = high.scale(&delta, 2, 1);
    let coordinate = TaylorCoordinate::balanced(
        p,
        &boundary.point,
        &target,
        &[
            high.sub(&boundary.point, &twice),
            high.add(&boundary.point, &twice),
        ],
    )?;
    assert!(matches!(coordinate, TaylorCoordinate::Mobius(_)));
    let (coefficients, chart, _) = probe.mapped_chart(
        &coordinate,
        &boundary.point,
        &boundary.values,
        16,
        &boundary.point,
    )?;
    let tolerance = p.tolerance(50);
    let conditioning = ConditioningChart::new(p, &coefficients)?;
    let candidate = StepTrial {
        pade: None,
        system: &probe,
        chart: &chart,
        coordinate: &coordinate,
        center: &boundary.point,
        delta: &delta,
        target: &target,
        coefficients: &coefficients,
        tolerance: &tolerance,
        conditioning: &conditioning,
        conditioning_digits: 20,
    }
    .evaluate(&delta, &mut String::new(), &mut FlowDiagnostics::default())?
    .expect("the exact mapped quadratic must pass the shared predicate");
    assert_eq!(candidate.0, target);
    assert_eq!(
        candidate.1,
        vec![p.rational(&Rational::from((3, 2))), p.i(1)]
    );
    assert_eq!(candidate.2, target);
    assert_eq!(*probe.rhs_points.borrow(), vec![target.clone(), midpoint]);
    assert_eq!(*probe.state_points.borrow(), vec![target]);
    Ok(())
}

#[test]
#[cfg(feature = "native")]
fn saved_segment_rejects_nonfinite_span_or_interval_coordinate() -> Result<()> {
    let p = Precision::decimal(60)?;
    let make = |end: C| crate::diffexp::EpsilonSolution {
        point: end.clone(),
        leading: 0,
        coefficients: vec![vec![p.i(7)]],
        diagnostics: FlowDiagnostics {
            working_bits: p.bits,
            ..Default::default()
        },
        segments: vec![TaylorSegment {
            pade: None,
            center: p.zero(),
            end,
            coefficients: vec![vec![p.i(7)]],
            working_bits: p.bits,
            coordinate: TaylorCoordinate::Identity,
            conditioning_digits: 20,
        }],
        verified_digits: None,
        comparison_errors: Vec::new(),
        checkpoints: Vec::new(),
    };
    for special in [rug::float::Special::Nan, rug::float::Special::Infinity] {
        let end = C::new(Float::with_val(p.bits, special), p.real(0));
        assert!(matches!(
            make(end).evaluate_segment(0, &p.zero()),
            Err(Error::InvalidInput(message)) if message.contains("nonfinite")
        ));
    }
    // MPFR stores this enormous finite exponent compactly. The current
    // complex division forms its squared norm, which overflows; NaN must not
    // bypass all interval comparisons. No large exact integer is constructed.
    let mut raw = rug::Float::with_val(p.bits, 1);
    raw <<= rug::float::exp_max() / 2 + 32;
    let end = C::new(Float::from_raw(raw), p.real(0));
    assert!(p.finite(&end));
    assert!(!p.finite(&p.div(&end, &end)));
    assert!(matches!(
        make(end.clone()).evaluate_segment(0, &end),
        Err(Error::InvalidInput(message)) if message.contains("nonfinite")
    ));
    Ok(())
}
