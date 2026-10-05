use symbolica::prelude::*;
use symbolica_amflow::contour::*;
use symbolica_amflow::*;
fn atom(s: &str) -> Atom {
    Atom::parse(s, "endpoint_separation_tests", Default::default()).unwrap()
}
fn planner(s: &str) -> PrescribedContour {
    PrescribedContour {
        variable: symbol!("endpoint_separation_tests::x"),
        prescriptions: vec![PolynomialPrescription {
            polynomial: atom(s),
            prescription: Prescription::PlusI0,
        }],
        unprescribed_side: Prescription::PlusI0,
    }
}
#[test]
fn near_endpoint_roots_inside_and_outside_are_exactly_classified() {
    let p = Precision::decimal(80).unwrap();
    for expression in ["x-1/10^20", "x-1+1/10^20", "(x-1/10^12)^2-2/10^26"] {
        let plan = planner(expression)
            .plan(p, &[], &p.zero(), &p.i(1))
            .unwrap();
        let count = if expression.starts_with('(') { 2 } else { 1 };
        assert_eq!(plan.crossings.len(), count);
        for c in &plan.crossings {
            assert!(c.point.re > p.real(0) && c.point.re < p.real(1));
        }
    }
    for expression in ["x+1/10^20", "x-1-1/10^20"] {
        assert!(
            planner(expression)
                .plan(p, &[], &p.zero(), &p.i(1))
                .unwrap()
                .crossings
                .is_empty()
        );
    }
}
#[test]
fn endpoint_roots_remain_rejected_and_nearby_reverse_path_matches() {
    let p = Precision::decimal(80).unwrap();
    for expression in ["x", "x-1"] {
        assert!(matches!(
            planner(expression).plan(p, &[], &p.zero(), &p.i(1)),
            Err(Error::InvalidInput(_))
        ));
    }
    let instance = planner("(x-1/10^15)*(x-1+1/10^15)");
    let forward = instance.plan(p, &[], &p.zero(), &p.i(1)).unwrap();
    let reverse = instance.plan(p, &[], &p.i(1), &p.zero()).unwrap();
    assert_eq!(forward.crossings.len(), 2);
    assert_eq!(reverse.crossings.len(), 2);
    for (a, b) in forward.crossings.iter().zip(reverse.crossings.iter().rev()) {
        assert_eq!(a.side, b.side);
        assert!(p.norm(&p.sub(&a.point, &b.point)) < p.tolerance(40));
    }
}
