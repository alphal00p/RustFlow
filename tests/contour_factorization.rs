use symbolica::domains::rational::RationalField;
use symbolica::poly::univariate::UnivariatePolynomial;
use symbolica::prelude::*;
use symbolica_amflow::contour::{PolynomialPrescription, PrescribedContour};
use symbolica_amflow::{Precision, Prescription};

fn atom(text: &str) -> Atom {
    Atom::parse(text, "np_w_long", Default::default()).unwrap()
}
fn polynomial(expression: &Atom) -> UnivariatePolynomial<RationalField> {
    let rational: RationalPolynomial<IntegerRing, u16> =
        expression.try_to_rational_polynomial(&Q, &Z, None).unwrap();
    assert!(rational.denominator.is_constant());
    let denominator = Rational::from(rational.denominator.get_constant());
    rational
        .numerator
        .to_univariate_from_univariate(0)
        .map_coeff(|c| Rational::from(c.clone()) / &denominator, Q)
}

#[test]
fn native_factors_preserve_the_degree_eight_crossings_and_full_polynomial_signs() {
    // This exact domain polynomial exposed the native cross-factor isolation
    // bottleneck. Here P+i0 is deliberately declared to check the full P' sign;
    // this is a polynomial-contour test, not an NP physical-sheet assertion.
    let input = atom(include_str!(
        "../repros/symbolica-cross-factor-isolation/polynomial.txt"
    ));
    let factors = [
        atom("3988330731937*path-2553573700"),
        atom(
            "257792982677968681370575937105*path^2+836584409914389788218059108*path-10045668600261177304974800",
        ),
        atom(
            "53348179342863732718927370685625*path^2+2101068496438385809561323778424*path+9658194937855178207169320800",
        ),
        atom(
            "40851405991964612066156792461895*path^3-2868377422595018253030284632700060*path^2+11972652032967197358219580100071*path+186698944121788600620212522900",
        ),
    ];
    let rational: RationalPolynomial<IntegerRing, u16> =
        input.try_to_rational_polynomial(&Q, &Z, None).unwrap();
    assert!(rational.denominator.get_constant() > Integer::zero());
    let product = factors.iter().fold(Atom::one(), |a, b| a * b);
    assert!(
        (rational.numerator.to_expression() - product)
            .expand()
            .is_zero()
    );
    let exact = factors.iter().map(polynomial).collect::<Vec<_>>();
    let intervals = [
        (Rational::from((1, 2000)), Rational::from((1, 1000))),
        (Rational::from((1, 250)), Rational::from((1, 200))),
        (Rational::from((1, 100)), Rational::from((11, 1000))),
    ];
    let brackets = [
        (0, vec![intervals[0].clone()]),
        (
            1,
            vec![
                (Rational::from((-1, 10)), Rational::zero()),
                intervals[1].clone(),
            ],
        ),
        (
            3,
            vec![
                (Rational::from((-1, 10)), Rational::zero()),
                intervals[2].clone(),
                (Rational::one(), Rational::from(100)),
            ],
        ),
    ];
    // Each of these factors changes sign in as many disjoint intervals as its
    // degree, proving all of its roots are simple and accounted for. The third
    // factor has strictly positive coefficients, hence no nonnegative root.
    for (i, intervals) in brackets {
        assert_eq!(exact[i].degree(), intervals.len());
        for (left, right) in intervals {
            assert!(exact[i].evaluate(&left) * exact[i].evaluate(&right) < Rational::zero());
        }
    }
    assert!(
        exact[2]
            .coefficients()
            .iter()
            .all(|c| *c > Rational::zero())
    );
    let path = atom("path");
    let AtomView::Var(variable) = path.as_view() else {
        panic!("expected one path variable")
    };
    let planner = PrescribedContour {
        variable: variable.get_symbol(),
        prescriptions: vec![PolynomialPrescription {
            polynomial: input.clone(),
            prescription: Prescription::PlusI0,
        }],
        unprescribed_side: Prescription::PlusI0,
    };
    let p = Precision::decimal(78).unwrap();
    let plan = planner.plan(p, &[], &p.zero(), &p.i(1)).unwrap();
    assert_eq!(plan.crossings.len(), 3);
    assert_eq!(plan.waypoints.len(), 10);
    let expected = [
        Prescription::MinusI0,
        Prescription::PlusI0,
        Prescription::MinusI0,
    ];
    let whole = polynomial(&input);
    for (i, ((crossing, interval), side)) in plan
        .crossings
        .iter()
        .zip(intervals)
        .zip(expected)
        .enumerate()
    {
        assert!(crossing.point.re > p.rational(&interval.0).re);
        assert!(crossing.point.re < p.rational(&interval.1).re);
        assert_eq!(crossing.side, side);
        assert_eq!(crossing.prescription_indices, vec![0]);
        let left = whole.evaluate(&plan.waypoints[3 * i].re.to_rational());
        let right = whole.evaluate(&plan.waypoints[3 * i + 2].re.to_rational());
        assert!(left.clone() * &right < Rational::zero());
        assert_eq!(left > Rational::zero(), side == Prescription::MinusI0);
        assert_eq!(
            plan.waypoints[3 * i + 1].im > p.real(0),
            side == Prescription::PlusI0
        );
    }
}
