use symbolica::domains::float::{Complex, Float};
use symbolica_amflow::{ComplexFloat, Precision};

// The pre-migration MPFR implementation is an independent regression oracle.
// Keep its complex operation order: each real product is rounded before the
// components are combined, and division multiplies by the conjugate first.
mod previous {
    use super::*;

    pub fn add(p: Precision, a: &ComplexFloat, b: &ComplexFloat) -> ComplexFloat {
        Complex::new(
            Float::with_val(p.bits, a.re.as_raw() + b.re.as_raw()),
            Float::with_val(p.bits, a.im.as_raw() + b.im.as_raw()),
        )
    }

    pub fn sub(p: Precision, a: &ComplexFloat, b: &ComplexFloat) -> ComplexFloat {
        add(p, a, &Complex::new(-b.re.clone(), -b.im.clone()))
    }

    pub fn mul(p: Precision, a: &ComplexFloat, b: &ComplexFloat) -> ComplexFloat {
        let ac = rug::Float::with_val(p.bits, a.re.as_raw() * b.re.as_raw());
        let bd = rug::Float::with_val(p.bits, a.im.as_raw() * b.im.as_raw());
        let ad = rug::Float::with_val(p.bits, a.re.as_raw() * b.im.as_raw());
        let bc = rug::Float::with_val(p.bits, a.im.as_raw() * b.re.as_raw());
        Complex::new(
            Float::with_val(p.bits, ac - bd),
            Float::with_val(p.bits, ad + bc),
        )
    }

    pub fn div(p: Precision, a: &ComplexFloat, b: &ComplexFloat) -> ComplexFloat {
        let conjugate = Complex::new(b.re.clone(), -b.im.clone());
        let numerator = mul(p, a, &conjugate);
        let denominator = mul(p, b, &conjugate).re;
        Complex::new(
            Float::with_val(p.bits, numerator.re.as_raw() / denominator.as_raw()),
            Float::with_val(p.bits, numerator.im.as_raw() / denominator.as_raw()),
        )
    }
}

fn assert_same_bits(actual: &ComplexFloat, expected: &ComplexFloat, label: &str) {
    for (part, actual, expected) in [
        ("real", &actual.re, &expected.re),
        ("imaginary", &actual.im, &expected.im),
    ] {
        assert!(actual.is_finite() && expected.is_finite(), "{label}/{part}");
        assert_eq!(actual.prec(), expected.prec(), "{label}/{part} precision");
        assert_eq!(
            actual.as_raw().to_integer_exp(),
            expected.as_raw().to_integer_exp(),
            "{label}/{part} significand and exponent",
        );
        // Equality and integer conversion alone do not distinguish signed zero.
        assert_eq!(
            actual.is_sign_negative(),
            expected.is_sign_negative(),
            "{label}/{part} sign",
        );
    }
}

#[test]
fn native_rounding_matches_previous_mpfr_bits_for_mixed_complex_inputs() {
    let real_inputs = [
        (2, "0"),
        (53, "-0"),
        (2, "1"),
        (113, "-1"),
        (
            256,
            "1.000000000000000000000000000000000000000000000000000000000001",
        ),
        (
            256,
            "-1.000000000000000000000000000000000000000000000000000000000001",
        ),
        (24, "0.1"),
        (257, "-0.1"),
        (64, "1.125"),
        (127, "1.375"),
        (53, "1e100000"),
        (256, "-1e100000"),
        (53, "1e-100000"),
        (256, "-1e-100000"),
        (
            193,
            "3.14159265358979323846264338327950288419716939937510582097494",
        ),
        (
            17,
            "-2.71828182845904523536028747135266249775724709369995957496697",
        ),
    ]
    .map(|(bits, text)| Float::parse(text, Some(bits)).unwrap());
    let inputs = real_inputs
        .iter()
        .enumerate()
        .flat_map(|(i, real)| {
            [
                Complex::new(real.clone(), Float::with_val(53, 0)),
                Complex::new(
                    real.clone(),
                    real_inputs[(i + 5) % real_inputs.len()].clone(),
                ),
            ]
        })
        .collect::<Vec<_>>();

    for bits in [2, 3, 24, 53, 113, 256] {
        let p = Precision { bits };
        for (i, a) in inputs.iter().enumerate() {
            for (j, b) in inputs.iter().enumerate() {
                let label = format!("{bits} bits, inputs {i}/{j}");
                assert_same_bits(&p.add(a, b), &previous::add(p, a, b), &label);
                assert_same_bits(&p.sub(a, b), &previous::sub(p, a, b), &label);
                assert_same_bits(&p.mul(a, b), &previous::mul(p, a, b), &label);
                if !b.re.as_raw().is_zero() || !b.im.as_raw().is_zero() {
                    assert_same_bits(&p.div(a, b), &previous::div(p, a, b), &label);
                }
            }
        }
    }
}

#[test]
fn component_products_are_rounded_before_combining() {
    let p = Precision { bits: 3 };
    let a = Precision { bits: 128 }.parse("1.25", "1.5").unwrap();
    let b = Precision { bits: 128 }.parse("1.25", "1").unwrap();
    // At three bits 1.25^2 rounds to 1.5 before subtracting 1.5.
    // Rounding the exact complex result only at the end would leave real 1/16.
    assert_same_bits(&p.mul(&a, &b), &p.complex(0, 3), "component order");
}

#[test]
fn cancellation_preserves_input_bits_until_the_operation() {
    let p = Precision { bits: 24 };
    let high = Precision { bits: 512 };
    let large = high.parse("1e100", "-1e100").unwrap();
    let incremented = high.add(&large, &high.complex(1, -1));
    let opposite = high.neg(&large);
    assert_same_bits(
        &p.add(&incremented, &opposite),
        &p.complex(1, -1),
        "addition cancellation",
    );
    assert_same_bits(
        &p.sub(&incremented, &large),
        &p.complex(1, -1),
        "subtraction cancellation",
    );
    assert_eq!(p.round(&incremented), p.round(&large));
}
