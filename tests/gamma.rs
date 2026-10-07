use symbolica::prelude::*;
use symbolica_amflow::{ComplexFloat, Error, Precision};

#[test]
fn gamma_half_integer_values_preserve_arbitrary_precision_and_sign() {
    for digits in [30, 80, 140] {
        let p = Precision::decimal(digits + 10).unwrap();
        let sqrt_pi = ComplexFloat::new(p.real(0).pi().sqrt(), p.real(0));
        for (numerator, scale_numerator, scale_denominator) in
            [(1, 1, 1), (3, 1, 2), (-1, -2, 1), (-3, 4, 3)]
        {
            let argument = p.rational(&Rational::from((numerator, 2)));
            let actual = p.gamma_real(&argument.re).unwrap();
            let expected = p.scale(&sqrt_pi, scale_numerator, scale_denominator);
            assert!(p.close(&actual, &expected, digits), "Gamma({argument})");
            assert_eq!(actual.re.prec(), p.bits);
            assert_eq!(actual.im, p.real(0));
        }
    }
}

#[test]
fn gamma_recurrence_and_refinement_hold_away_from_half_integers() {
    let low = Precision::decimal(70).unwrap();
    let high = Precision::decimal(110).unwrap();
    for value in [Rational::from((1, 10)), Rational::from((-13, 10))] {
        let argument = high.rational(&value);
        let gamma = high.gamma_real(&argument.re).unwrap();
        let successor = high
            .gamma_real(&high.add(&argument, &high.i(1)).re)
            .unwrap();
        assert!(high.close(&successor, &high.mul(&argument, &gamma), 100));
        let lower = low.gamma_real(&low.rational(&value).re).unwrap();
        assert!(high.close(&lower, &gamma, 65));
    }
}

#[test]
fn gamma_poles_are_typed_numerical_failures() {
    let p = Precision::decimal(50).unwrap();
    for pole in [0, -1, -2, -6] {
        assert!(matches!(
            p.gamma_real(&p.real(pole)),
            Err(Error::Numerical(_))
        ));
    }
}

#[test]
fn gamma_residues_remain_finite_for_tiny_and_near_pole_arguments() {
    let p = Precision::decimal(100).unwrap();
    // Conversion of the argument to a machine float would underflow to zero.
    let tiny = p.parse("1e-400", "0").unwrap();
    let residue = p.mul(&tiny, &p.gamma_real(&tiny.re).unwrap());
    assert!(p.close(&residue, &p.i(1), 90));

    let displacement = p.parse("1e-60", "0").unwrap();
    // Astro arithmetic can retain stored guard bits beyond the requested
    // precision. Match gamma_real's input rounding before measuring the pole
    // distance; otherwise these are two differently represented arguments.
    let near_minus_one = p.round(&p.add(&p.i(-1), &displacement));
    // Use the represented distance to the pole rather than assume decimal
    // conversion and subtraction were exact.
    let distance = p.add(&near_minus_one, &p.i(1));
    let residue = p.mul(&distance, &p.gamma_real(&near_minus_one.re).unwrap());
    assert!(p.close(&residue, &p.i(-1), 55));
}

#[cfg(feature = "native")]
#[test]
fn native_gamma_matches_mpfr_at_the_requested_working_precision() {
    for digits in [30, 80, 140] {
        let p = Precision::decimal(digits).unwrap();
        for text in ["0.1", "-1.3", "0.5", "1e-400", "9.75"] {
            let argument = p.parse(text, "0").unwrap();
            let expected = rug::Float::with_val(p.bits, argument.re.as_raw()).gamma();
            let actual = p.gamma_real(&argument.re).unwrap();
            assert_eq!(actual.re.as_raw(), &expected);
            assert_eq!(actual.re.prec(), p.bits);
        }
    }
}
