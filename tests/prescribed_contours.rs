use symbolica::prelude::*;
use symbolica_amflow::contour::*;
use symbolica_amflow::*;

fn polynomial(text: &str) -> Atom {
    Atom::parse(text, "contour_geometry_audit", Default::default()).unwrap()
}
fn planner(polynomials: &[&str]) -> PrescribedContour {
    PrescribedContour {
        variable: symbol!("contour_geometry_audit::x"),
        prescriptions: polynomials
            .iter()
            .map(|text| PolynomialPrescription {
                polynomial: polynomial(text),
                prescription: Prescription::PlusI0,
            })
            .collect(),
        unprescribed_side: Prescription::PlusI0,
    }
}
fn assert_upper_apex_below(plan: &PlannedContour, limit: &Float, p: Precision) {
    assert_eq!(plan.crossings.len(), 1);
    assert_eq!(plan.crossings[0].side, Prescription::PlusI0);
    assert_eq!(plan.waypoints.len(), 4);
    let apex = &plan.waypoints[1];
    assert!(apex.im > p.real(0));
    assert!(
        apex.im < *limit,
        "detour must not enclose a separate complex singularity"
    );
}

#[test]
fn a_close_complex_pole_is_not_merged_with_a_real_prescribed_zero() {
    let planner = planner(&["x-1/2"]);
    let singularities = [polynomial("(x-1/2)^2+1/10^100")];
    let low = Precision::decimal(60).unwrap();
    let limit = low.parse("2.5e-51", "0").unwrap().re;
    match planner.plan(low, &singularities, &low.zero(), &low.i(1)) {
        Ok(plan) => assert_upper_apex_below(&plan, &limit, low),
        Err(Error::Accuracy(_)) => {}
        Err(error) => panic!("expected a resolved safe contour or a precision error: {error}"),
    }
    let high = Precision::decimal(120).unwrap();
    let result = planner
        .plan(high, &singularities, &high.zero(), &high.i(1))
        .unwrap();
    assert_upper_apex_below(&result, &high.parse("2.5e-51", "0").unwrap().re, high);
}

#[test]
fn every_complex_root_of_a_declared_polynomial_constrains_its_detour() {
    let p = Precision::decimal(60).unwrap();
    let result = planner(&["(x-1/2)*((x-1/2)^2+1/10000)"])
        .plan(p, &[], &p.zero(), &p.i(1))
        .unwrap();
    assert_upper_apex_below(&result, &p.parse("0.0025", "0").unwrap().re, p);
}

#[test]
fn shared_exact_factors_merge_across_different_polynomials_without_near_tests() {
    let p = Precision::decimal(60).unwrap();
    let result = planner(&["x^2-2", "3*(x^2-2)*(x+5)"])
        .plan(p, &[], &p.i(1), &p.i(2))
        .unwrap();
    assert_eq!(result.crossings.len(), 1);
    assert_eq!(result.crossings[0].prescription_indices, vec![0, 1]);
    assert_eq!(result.crossings[0].side, Prescription::PlusI0);
    assert!(matches!(
        planner(&["x^2-2", "-3*(x^2-2)*(x+5)"]).plan(p, &[], &p.i(1), &p.i(2)),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn input_exactness_and_cancellation_are_checked_before_root_search() {
    let p = Precision::decimal(60).unwrap();
    let original = planner(&["x-1/2"]);
    let mut approximate = original.clone();
    approximate.prescriptions[0].polynomial += Atom::num(p.parse("0.125", "0").unwrap().re);
    assert!(matches!(
        approximate.plan(p, &[], &p.zero(), &p.i(1)),
        Err(Error::InvalidInput(_))
    ));
    let mut complex = original.clone();
    complex.prescriptions[0].polynomial +=
        Atom::num(Complex::new(Rational::zero(), Rational::one()));
    assert!(matches!(
        complex.plan(p, &[], &p.zero(), &p.i(1)),
        Err(Error::Unsupported(_))
    ));
    let context = RunContext::default();
    context.cancellation.cancel();
    assert!(matches!(
        original.plan_with_context(p, &[], &p.zero(), &p.i(1), &context),
        Err(Error::Cancelled)
    ));
}

#[test]
fn polynomial_prescription_transports_scalar_system_and_reverses_with_refinement() -> Result<()> {
    let x = symbol!("contour_scientific::x");
    let a = Atom::var(x) - Atom::num((1, 5));
    let b = Atom::var(x) - Atom::num((4, 5));
    let declaration = &a * &b;
    let planner = PrescribedContour {
        variable: x,
        prescriptions: vec![PolynomialPrescription {
            polynomial: declaration.clone(),
            prescription: Prescription::PlusI0,
        }],
        unprescribed_side: Prescription::PlusI0,
    };
    let system = DifferentialSystem {
        variable: x,
        matrix: vec![vec![Atom::num((1, 2)) / a + Atom::num((1, 4)) / b]],
    };
    let mut previous = None;
    for (working, order, digits) in [(60, 64, 24), (80, 96, 36)] {
        let p = Precision::decimal(working)?;
        let options = FlowOptions {
            digits,
            guard_digits: working - digits,
            series_order: order,
            ..Default::default()
        };
        let compiled = system.compile(p, &Default::default())?;
        let forward = planner.plan_with_context(
            p,
            std::slice::from_ref(&declaration),
            &p.zero(),
            &p.i(1),
            &RunContext::default(),
        )?;
        assert_eq!(forward.crossings[0].side, Prescription::MinusI0);
        assert_eq!(forward.crossings[1].side, Prescription::PlusI0);
        let result = compiled.transport(
            &BoundaryData {
                point: p.zero(),
                values: vec![p.i(1)],
            },
            &forward.waypoints,
            &options,
            &RunContext::default(),
        )?;
        let error = p.norm(&p.sub(&result.values[0], &p.complex(1, 1)));
        assert!(error < p.tolerance(digits));
        let reverse = planner.plan(p, std::slice::from_ref(&declaration), &p.i(1), &p.zero())?;
        let restored = compiled.transport(
            &BoundaryData {
                point: p.i(1),
                values: result.values.clone(),
            },
            &reverse.waypoints,
            &options,
            &RunContext::default(),
        )?;
        let reverse_error = p.norm(&p.sub(&restored.values[0], &p.i(1)));
        assert!(reverse_error < p.tolerance(digits));
        eprintln!(
            "scalar prescription working={working} order={order} requested={digits} value={} error={} reverse_error={} forward_steps={} reverse_steps={}",
            result.values[0],
            error,
            reverse_error,
            result.diagnostics.steps,
            restored.diagnostics.steps
        );
        if let Some(previous) = previous {
            let difference = p.norm(&p.sub(&p.round(&previous), &result.values[0]));
            assert!(difference < p.tolerance(24));
            eprintln!("scalar prescription refinement_difference={difference}");
        }
        previous = Some(result.values[0].clone());
    }
    Ok(())
}

#[test]
#[ignore = "Symbolica75f8350 isolate_roots stalls on this exact quadratic; separate bounded MRE and direct exact-disk geometry test retained"]
fn a_margin_sized_obstacle_cannot_collapse_triangle_waypoints_onto_the_zero() {
    let p = Precision::decimal(60).unwrap();
    assert_eq!(p.bits / 5, 43);
    let x = symbol!("contour_rounding_audit::x");
    let center = Atom::num((1, 2));
    // Strictly above the proposed distance margin, but the remaining radius
    // would be below the MPFR spacing around x=1/2 after subtraction.
    let height = Atom::num(16) / Atom::num(10).pow(43) + Atom::one() / Atom::num(10).pow(68);
    let singularities = [(Atom::var(x) - &center).pow(2) + height.pow(2)];
    let planner = PrescribedContour {
        variable: x,
        prescriptions: vec![PolynomialPrescription {
            polynomial: Atom::var(x) - center,
            prescription: Prescription::PlusI0,
        }],
        unprescribed_side: Prescription::PlusI0,
    };
    match planner.plan(p, &singularities, &p.zero(), &p.i(1)) {
        Err(Error::Accuracy(_)) => {}
        Ok(result) => {
            assert_eq!(result.crossings.len(), 1);
            let center = &result.crossings[0].point.re;
            assert!(result.waypoints[0].re < *center);
            assert!(result.waypoints[2].re > *center);
            assert!(result.waypoints[1].im > p.real(0));
        }
        Err(error) => panic!("expected a resolved triangle or a precision error: {error}"),
    }
}

mod prescription_validation {
    use super::*;
    fn plan(polynomials: &[&str]) -> PrescribedContour {
        PrescribedContour {
            variable: symbol!("contour_test::x"),
            prescriptions: polynomials
                .iter()
                .map(|s| PolynomialPrescription {
                    polynomial: Atom::parse(s, "contour_test", Default::default()).unwrap(),
                    prescription: Prescription::PlusI0,
                })
                .collect(),
            unprescribed_side: Prescription::PlusI0,
        }
    }
    #[test]
    fn varying_polynomial_sides_and_reverse_path() {
        let p = Precision::decimal(60).unwrap();
        let planner = plan(&["(x-1/5)*(x-4/5)", "2*(x-1/5)*(x-4/5)"]);
        let f = planner.plan(p, &[], &p.zero(), &p.i(1)).unwrap();
        assert_eq!(f.crossings.len(), 2);
        assert_eq!(f.crossings[0].side, Prescription::MinusI0);
        assert_eq!(f.crossings[1].side, Prescription::PlusI0);
        assert_eq!(f.crossings[0].prescription_indices, vec![0, 1]);
        let r = planner.plan(p, &[], &p.i(1), &p.zero()).unwrap();
        let mut forward = vec![p.zero()];
        forward.extend(f.waypoints);
        let mut reverse = vec![p.i(1)];
        reverse.extend(r.waypoints);
        reverse.reverse();
        assert_eq!(forward, reverse);
    }
    #[test]
    fn conflicts_multiplicity_and_endpoint_are_errors() {
        let p = Precision::decimal(60).unwrap();
        assert!(matches!(
            plan(&["x-1/2", "1/2-x"]).plan(p, &[], &p.zero(), &p.i(1)),
            Err(Error::Unsupported(_))
        ));
        assert!(matches!(
            plan(&["(x-1/2)^2"]).plan(p, &[], &p.zero(), &p.i(1)),
            Err(Error::Unsupported(_))
        ));
        assert!(matches!(
            plan(&["x"]).plan(p, &[], &p.zero(), &p.i(1)),
            Err(Error::InvalidInput(_))
        ));
        assert!(matches!(
            plan(&["x+y"]).plan(p, &[], &p.zero(), &p.i(1)),
            Err(Error::Unsupported(_))
        ));
        assert!(matches!(
            plan(&["1/x"]).plan(p, &[], &p.i(-1), &p.i(1)),
            Err(Error::Unsupported(_))
        ));
    }
}

#[test]
fn a_coarse_real_obstacle_from_another_polynomial_is_refined_before_detouring() {
    // Extracted from the PH1-to-PH6 source-domain obstacles: native isolation
    // gives 357/577 an initial radius 1/4, covering the distinct crossing 5/6.
    // The exact roots are separated by more than 1/5, so precision is ample.
    let p = Precision::decimal(50).unwrap();
    let plan = planner(&["577*x-357"])
        .plan(p, &[polynomial("x-5/6")], &p.zero(), &p.i(1))
        .unwrap();
    assert_eq!(plan.crossings.len(), 2);
    assert!(p.close(
        &plan.crossings[0].point,
        &p.rational(&Rational::from((357, 577))),
        14
    ));
    assert!(p.close(
        &plan.crossings[1].point,
        &p.rational(&Rational::from((5, 6))),
        14
    ));
    assert!(plan.waypoints[2].re < plan.waypoints[3].re);
    assert!(plan.waypoints[1].im > p.real(0));
    assert!(plan.waypoints[4].im > p.real(0));
}
