use symbolica::prelude::*;
use symbolica_amflow::contour::{PolynomialPrescription, PrescribedContour};
use symbolica_amflow::{Precision, Prescription, Result};

// RootCache keys only exact polynomial coefficients, not namespaces. These
// centers differ from every other contour regression; each is exercised once.
fn check(center: (i64, i64), coefficient: (i64, i64), sign: i64) -> Result<()> {
    let x = symbol!("contour_fresh_sign::x");
    let z = Atom::var(x) - Atom::num(center);
    let planner = PrescribedContour {
        variable: x,
        prescriptions: vec![PolynomialPrescription {
            polynomial: Atom::num(sign) * &z * (&z * &z + Atom::num(coefficient)),
            prescription: Prescription::PlusI0,
        }],
        unprescribed_side: Prescription::PlusI0,
    };
    let p = Precision::decimal(60)?;
    let plan = planner.plan(p, &[], &p.zero(), &p.i(1))?;
    assert_eq!(plan.crossings.len(), 1);
    assert_eq!(plan.waypoints.len(), 4);
    assert_eq!(
        plan.crossings[0].side,
        if sign > 0 {
            Prescription::PlusI0
        } else {
            Prescription::MinusI0
        }
    );
    assert_eq!(plan.crossings[0].prescription_indices, vec![0]);
    assert!(p.norm(&plan.waypoints[1]) > p.real(0));
    assert!((plan.waypoints[1].im > p.real(0)) == (sign > 0));
    assert!(plan.waypoints[1].im.to_rational().abs() < p.parse("0.0025", "0")?.re.to_rational());
    Ok(())
}
#[test]
fn fresh_real_certificate_retains_small_positive_full_derivative() -> Result<()> {
    check((137, 311), (1, 10000), 1)
}
#[test]
fn fresh_real_certificate_retains_small_negative_full_derivative() -> Result<()> {
    check((139, 313), (1, 20000), -1)
}

#[test]
fn fresh_real_factor_certifies_tiny_companion_derivative() -> Result<()> {
    let x = symbol!("contour_fresh_tiny_sign::x");
    let z = Atom::var(x) - Atom::num((73, 163));
    let coefficient =
        Atom::parse("1/10^30", "contour_fresh_tiny_sign", Default::default()).unwrap();
    let planner = PrescribedContour {
        variable: x,
        prescriptions: vec![PolynomialPrescription {
            polynomial: &z * (&z * &z + coefficient),
            prescription: Prescription::PlusI0,
        }],
        unprescribed_side: Prescription::PlusI0,
    };
    let p = Precision::decimal(70)?;
    let plan = planner.plan(p, &[], &p.zero(), &p.i(1))?;
    assert_eq!(plan.crossings.len(), 1);
    assert_eq!(plan.crossings[0].side, Prescription::PlusI0);
    assert!(plan.waypoints[1].im > p.real(0));
    assert!(plan.waypoints[1].im < p.parse("2.5e-16", "0")?.re);
    Ok(())
}

#[test]
fn cached_root_certificates_rebind_polynomial_variables() -> Result<()> {
    let p = Precision::decimal(60)?;
    for variable in [
        symbol!("contour_cache_rebind::first"),
        symbol!("contour_cache_rebind::second"),
    ] {
        let z = Atom::var(variable) - Atom::num((149, 337));
        let planner = PrescribedContour {
            variable,
            prescriptions: vec![PolynomialPrescription {
                polynomial: -&z * (&z * &z + Atom::num((1, 10000))),
                prescription: Prescription::PlusI0,
            }],
            unprescribed_side: Prescription::PlusI0,
        };
        let plan = planner.plan(p, &[], &p.zero(), &p.i(1))?;
        assert_eq!(plan.crossings.len(), 1);
        assert_eq!(plan.crossings[0].side, Prescription::MinusI0);
        assert!(plan.waypoints[1].im < p.real(0));
    }
    Ok(())
}
