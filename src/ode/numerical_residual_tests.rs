use super::*;
use crate::ode::source::{ExactPolynomialRow, ExactSpecialization};
use crate::ode::{SeriesSystem, TaylorCoordinate};
use crate::{DifferentialSystem, FlowOptions, ResidualArithmetic, RunContext};

fn atom(s: &str) -> Atom {
    Atom::parse(s, "numerical_residual", Default::default()).unwrap()
}
fn variable(s: &str) -> Symbol {
    match atom(s).as_view() {
        AtomView::Var(v) => v.get_symbol(),
        _ => unreachable!(),
    }
}
fn row(
    p: Precision,
    denominator: &str,
    factors: &[(&str, usize)],
    numerator: &str,
) -> ExactPolynomialRow {
    let exact = ExactSpecialization::new(p, &Default::default()).unwrap();
    let polynomial = |s: &str| exact.polynomial(&atom(s), variable("x")).unwrap();
    ExactPolynomialRow {
        denominator: polynomial(denominator),
        denominator_factors: factors.iter().map(|(a, k)| (polynomial(a), *k)).collect(),
        entries: vec![(0, polynomial(numerator))],
    }
}
fn chart(
    p: Precision,
    row: ExactPolynomialRow,
    center: &C,
    coefficients: &[Vec<C>],
    channels: usize,
) -> RationalResidualChart {
    RationalResidualChart::new_with_options(
        p,
        &[row],
        center,
        coefficients,
        channels,
        &FlowOptions {
            residual_arithmetic: ResidualArithmetic::NumericalEstimate,
            ..Default::default()
        },
        &RunContext::default(),
    )
    .unwrap()
}

#[test]
fn repeated_high_degree_factors_rescue_inconclusive_expansion_without_roots() {
    let p = Precision { bits: 256 };
    let chart = chart(
        p,
        row(p, "(1+x)^32", &[("1+x", 32)], "1"),
        &p.zero(),
        &[vec![p.i(1)]],
        1,
    );
    let radius = p.real(1) / p.real(4);
    let denominators = chart.numerical_denominators.as_ref().unwrap();
    assert!(numerical_lower(p, &denominators[0].expanded, &radius).is_none());
    let expected = (p.real(3) / p.real(4)).pow(32);
    assert_eq!(denominators[0].lower(p, &radius), Some(expected.clone()));
    let bounds = chart.defect_bounds(p, &p.scale(&p.i(1), 1, 4)).unwrap();
    assert_eq!(bounds, [radius / expected]);
    assert!(bounds[0] > p.real(0));
}

#[test]
fn complex_and_nonzero_centers_keep_content_sign_and_multiplicity() {
    let p = Precision { bits: 256 };
    for center in [p.scale(&p.i(1), 1, 2), p.scale(&p.complex(1, 1), 1, 2)] {
        let c = chart(
            p,
            row(p, "-6*(1+x)^12", &[("-6", 1), ("1+x", 12)], "1"),
            &center,
            &[vec![p.i(1)]],
            1,
        );
        let radius = p.real(1) / p.real(4);
        let d = &c.numerical_denominators.as_ref().unwrap()[0];
        assert!(numerical_lower(p, &d.expanded, &radius).is_none());
        let mut expected = p.real(6);
        let factor = p.norm(&p.add(&p.i(1), &center)) - &radius;
        for _ in 0..12 {
            expected *= &factor;
        }
        assert_eq!(d.lower(p, &radius), Some(expected));
        assert!(c.defect_bounds(p, &p.scale(&p.i(1), 1, 4)).unwrap()[0] > p.real(0));
    }
    let c = chart(
        p,
        row(p, "-3", &[("-3", 1)], "1"),
        &p.complex(2, 1),
        &[vec![p.i(1)]],
        1,
    );
    assert_eq!(
        c.numerical_denominators.as_ref().unwrap()[0].lower(p, &p.real(10)),
        Some(p.real(3))
    );
    let mut complex_content = row(p, "3", &[("3", 1)], "1");
    let content = crate::ode::source::Gaussian::new(Rational::zero(), Rational::from(3));
    complex_content.denominator[0] = content.clone();
    complex_content.denominator_factors[0].0[0] = content;
    let c = chart(p, complex_content, &p.complex(2, 1), &[vec![p.i(1)]], 1);
    assert_eq!(
        c.numerical_denominators.as_ref().unwrap()[0].lower(p, &p.real(10)),
        Some(p.real(3))
    );
}

#[test]
fn expanded_lower_is_kept_when_factors_are_inconclusive_or_weaker() {
    let p = Precision { bits: 160 };
    // The expanded triangle bound is stronger than the product estimate.
    let c = chart(
        p,
        row(p, "1-x^2", &[("1-x", 1), ("1+x", 1)], "1"),
        &p.zero(),
        &[vec![p.i(1)]],
        1,
    );
    let radius = p.real(1) / p.real(2);
    let d = &c.numerical_denominators.as_ref().unwrap()[0];
    assert_eq!(d.lower(p, &radius), Some(p.real(3) / p.real(4)));
    // Both grouped factors are inconclusive, but their expansion is positive.
    let c = chart(
        p,
        row(p, "(1-x^2)^2", &[("(1-x)^2", 1), ("(1+x)^2", 1)], "1"),
        &p.zero(),
        &[vec![p.i(1)]],
        1,
    );
    let d = &c.numerical_denominators.as_ref().unwrap()[0];
    assert!(
        d.factors
            .iter()
            .all(|(a, _)| numerical_lower(p, a, &radius).is_none())
    );
    assert_eq!(d.lower(p, &radius), Some(p.real(7) / p.real(16)));
}

#[test]
fn high_source_degrees_are_not_truncated_to_candidate_order() {
    let p = Precision { bits: 160 };
    let c = chart(p, row(p, "1", &[], "x^32"), &p.zero(), &[vec![p.i(1)]], 1);
    assert_eq!(c.residuals[0].numerator.len(), 33);
    assert_eq!(
        c.defect_bounds(p, &p.scale(&p.i(1), 1, 2)).unwrap(),
        [(p.real(1) / p.real(2)).pow(33)]
    );
}

#[test]
fn known_low_order_solution_and_genuine_defect_use_original_budget() {
    let p = Precision { bits: 160 };
    let candidate = [vec![p.i(1)], vec![p.i(1)]];
    let exact = chart(
        p,
        row(p, "1+x", &[("1+x", 1)], "1"),
        &p.zero(),
        &candidate,
        1,
    );
    assert_eq!(
        exact.defect_bounds(p, &p.scale(&p.i(1), 1, 4)).unwrap(),
        [p.real(0)]
    );
    let bad = chart(p, row(p, "1", &[], "1"), &p.zero(), &candidate, 1);
    let bounds = bad
        .defect_bounds_with_budget(p, &p.scale(&p.i(1), 1, 4), &[p.i(1)], &p.tolerance(30))
        .unwrap();
    assert_eq!(bounds, [p.real(1) / p.real(16)]);
    assert!(bounds[0] > p.tolerance(30));
}

#[test]
fn factors_are_shared_by_physical_row_across_epsilon_channels() {
    let p = Precision { bits: 160 };
    let c = chart(
        p,
        row(p, "(1+x)^12", &[("1+x", 12)], "1"),
        &p.zero(),
        &[vec![p.i(1), p.i(2), p.i(3)]],
        3,
    );
    assert_eq!(c.numerical_denominators.as_ref().unwrap().len(), 1);
    assert!(c.residuals.iter().all(|r| r.denominator.is_empty()));
    let bounds = c.defect_bounds(p, &p.scale(&p.i(1), 1, 4)).unwrap();
    assert_eq!(bounds[1], bounds[0].clone() * p.real(2));
    assert_eq!(bounds[2], bounds[0].clone() * p.real(3));
}

#[test]
fn inconclusive_factors_zero_denominators_and_nonfinite_values_reject() {
    let p = Precision { bits: 160 };
    for r in [
        row(p, "1-x^2", &[("1-x", 1), ("1+x", 1)], "0"),
        row(p, "(1-x)^2", &[("1-x", 2)], "0"),
        row(p, "0", &[("0", 1)], "0"),
    ] {
        let c = chart(p, r, &p.zero(), &[vec![p.i(1)]], 1);
        assert!(matches!(
            c.defect_bounds(p, &p.i(2)),
            Err(Error::Accuracy(_))
        ));
    }
    let c = chart(
        p,
        row(p, "0", &[("0", 1)], "0"),
        &p.zero(),
        &[vec![p.i(1)]],
        1,
    );
    assert!(matches!(
        c.defect_bounds(p, &p.zero()),
        Err(Error::Accuracy(_))
    ));
    let nan = C::new(p.real(0) / p.real(0), p.real(0));
    let inf = C::new(p.real(1) / p.real(0), p.real(0));
    let c = chart(p, row(p, "1", &[], "0"), &p.zero(), &[vec![p.i(1)]], 1);
    for value in [nan, inf] {
        assert!(c.defect_bounds(p, &value).is_err());
        let o = FlowOptions {
            residual_arithmetic: ResidualArithmetic::NumericalEstimate,
            ..Default::default()
        };
        assert!(
            RationalResidualChart::new_with_options(
                p,
                &[row(p, "1", &[], "0")],
                &p.zero(),
                &[vec![value]],
                1,
                &o,
                &RunContext::default()
            )
            .is_err()
        );
    }
    let mut c = c;
    c.residuals[0].numerator[0] = C::new(p.real(1) / p.real(0), p.real(0));
    assert!(matches!(
        c.defect_bounds(p, &p.i(1)),
        Err(Error::Accuracy(_))
    ));
}

#[test]
fn numerical_selection_preserves_recurrence_default_and_mapped_epsilon_owners() {
    let p = Precision { bits: 192 };
    assert_eq!(ResidualArithmetic::default(), ResidualArithmetic::Ball);
    assert_eq!(
        "numerical_estimate".parse::<ResidualArithmetic>().unwrap(),
        ResidualArithmetic::NumericalEstimate
    );
    assert_eq!(
        ResidualArithmetic::NumericalEstimate.as_str(),
        "numerical_estimate"
    );
    let ordinary = DifferentialSystem {
        variable: variable("x"),
        matrix: vec![vec![atom("1/(1-x)")]],
    }
    .compile(p, &Default::default())
    .unwrap();
    let center = p.scale(&p.complex(1, 1), 1, 8);
    let (baseline, default) = ordinary.local_chart(&center, &[p.i(1)], 12, &()).unwrap();
    assert_eq!(
        default.selected_arithmetic(),
        Some(ResidualArithmetic::Ball)
    );
    let o = FlowOptions {
        residual_arithmetic: ResidualArithmetic::NumericalEstimate,
        ..Default::default()
    };
    let context = RunContext::default();
    let (candidate, numerical) = ordinary
        .local_chart_with_options(&center, &[p.i(1)], 12, &(), &o, &context)
        .unwrap();
    assert_eq!(candidate, baseline);
    assert_eq!(
        numerical.selected_arithmetic(),
        Some(ResidualArithmetic::NumericalEstimate)
    );
    let coordinate =
        TaylorCoordinate::balanced(p, &center, &p.scale(&p.i(1), 1, 4), &ordinary.poles).unwrap();
    let (_, mapped, _) = ordinary
        .mapped_chart_with_options(&coordinate, &center, &[p.i(1)], 12, &(), &o, &context)
        .unwrap();
    assert_eq!(
        mapped.selected_arithmetic(),
        Some(ResidualArithmetic::NumericalEstimate)
    );
    let eps = crate::diffexp::EpsilonSystem::from_differential_system(
        &DifferentialSystem {
            variable: variable("x"),
            matrix: vec![vec![atom("(1+eps)/(1-x)")]],
        },
        variable("eps"),
        2,
    )
    .unwrap()
    .compile(p, &Default::default())
    .unwrap();
    let values = [p.i(1), p.zero(), p.zero()];
    let (_, c) = eps
        .local_chart_with_options(&center, &values, 12, &(), &o, &context)
        .unwrap();
    assert_eq!(
        c.selected_arithmetic(),
        Some(ResidualArithmetic::NumericalEstimate)
    );
    let (_, c, _) = eps
        .mapped_chart_with_options(&coordinate, &center, &values, 12, &(), &o, &context)
        .unwrap();
    assert_eq!(
        c.selected_arithmetic(),
        Some(ResidualArithmetic::NumericalEstimate)
    );
    context.cancellation.cancel();
    assert!(matches!(
        ordinary.local_chart_with_options(&center, &[p.i(1)], 12, &(), &o, &context),
        Err(Error::Cancelled)
    ));
}

#[test]
fn numerical_physical_transport_keeps_refinement_source_floor_and_honest_provenance() {
    use crate::Prescription;
    use crate::kinematics::KinematicSystem;
    use crate::physical_transport::PhysicalTransport;
    use crate::transport_cache::*;
    use std::collections::BTreeMap;
    let s = variable("s");
    let engine = PhysicalTransport::new(
        KinematicSystem {
            epsilon: variable("eps"),
            derivatives: BTreeMap::from([(s, vec![vec![atom("eps/s")]])]),
        },
        &[atom("I1")],
        &Atom::num(1),
        Prescription::PlusI0,
        "positive s",
    )
    .unwrap();
    let p = Precision::decimal(90).unwrap();
    let range = EpsilonRange::new(0, 2).unwrap();
    let coefficients = vec![vec![p.i(1)], vec![p.zero()], vec![p.zero()]];
    let source = CachedBoundary {
        identity: engine.identity().clone(),
        point: CachedPoint::Exact(BTreeMap::from([(s, Atom::num(1))])),
        kind: PointKind::Physical,
        range,
        coefficients,
        accuracy: BoundaryAccuracy::supplied(
            70,
            p.bits,
            vec![vec![p.real(0)]; 3],
            "exact boundary",
        )
        .unwrap(),
    };
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    let o = FlowOptions {
        digits: 20,
        guard_digits: 20,
        series_order: 32,
        residual_arithmetic: ResidualArithmetic::NumericalEstimate,
        ..Default::default()
    };
    let destination = BTreeMap::from([(s, Atom::num(2))]);
    let mut cache = RustFlowCache::default();
    cache.insert(source.clone()).unwrap();
    let result = engine
        .evaluate_to(
            &mut cache,
            &destination,
            range,
            &o,
            &RunContext::default(),
            &policy,
        )
        .unwrap();
    let value = &result.boundary;
    assert!(result.transport.as_ref().unwrap().verified_digits.is_some());
    assert!(value.accuracy.verified_digits() >= o.digits);
    let log = p.log(&p.i(2));
    assert!(p.close(&value.coefficients[1][0], &log, 20));
    assert!(p.close(
        &value.coefficients[2][0],
        &p.scale(&p.mul(&log, &log), 1, 2),
        20
    ));
    assert!(
        value
            .accuracy
            .provenance()
            .contains("independently refined")
    );
    assert!(
        value
            .accuracy
            .provenance()
            .contains("propagated cached input uncertainty")
    );
    assert!(
        value
            .accuracy
            .provenance()
            .contains("finite-MPFR numerical local-defect estimate")
    );
    assert!(
        value
            .accuracy
            .provenance()
            .contains("source conversion and residual rounding unenclosed")
    );
    assert!(
        cache
            .entries()
            .iter()
            .filter(|b| {
                b.kind == PointKind::Physical
                    && b.point.rounded_coordinates_as_exact().unwrap()
                        != source.point.rounded_coordinates_as_exact().unwrap()
            })
            .all(|b| b
                .accuracy
                .provenance()
                .contains("finite-MPFR numerical local-defect estimate"))
    );
    // Source accuracy cannot be promoted by choosing numerical residuals.
    let mut weak = source;
    weak.accuracy = BoundaryAccuracy::supplied(
        20,
        p.bits,
        vec![vec![p.real(0)]; 3],
        "only twenty source digits",
    )
    .unwrap();
    let mut cache = RustFlowCache::default();
    cache.insert(weak).unwrap();
    assert!(matches!(
        engine.evaluate_to(
            &mut cache,
            &destination,
            range,
            &o,
            &RunContext::default(),
            &policy
        ),
        Err(Error::Accuracy(_))
    ));
    assert_eq!(cache.len(), 1);
}
