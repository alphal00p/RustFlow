//! Exact rational rounding oracles and unchanged general complex formulas.
use super::*;

fn rounded(p: Precision, value: Rational) -> Float {
    // The portable arithmetic backend stores full words, whereas with_val
    // rounds to the requested bit count. Compare at that requested precision.
    let value = Float::from_rational_round(&value, p.bits, Nearest);
    Float::with_val(p.bits, value.as_raw())
}

fn legacy_mul(p: Precision, a: &ComplexFloat, b: &ComplexFloat) -> ComplexFloat {
    let ac = a.re.mul_round(&b.re, p.bits, Nearest);
    let bd = a.im.mul_round(&b.im, p.bits, Nearest);
    let ad = a.re.mul_round(&b.im, p.bits, Nearest);
    let bc = a.im.mul_round(&b.re, p.bits, Nearest);
    Complex::new(
        ac.sub_round(&bd, p.bits, Nearest),
        ad.add_round(&bc, p.bits, Nearest),
    )
}

fn legacy_div(p: Precision, a: &ComplexFloat, b: &ComplexFloat) -> ComplexFloat {
    let conjugate = Complex::new(b.re.clone(), -b.im.clone());
    let numerator = legacy_mul(p, a, &conjugate);
    let denominator = legacy_mul(p, b, &conjugate).re;
    Complex::new(
        numerator.re.div_round(&denominator, p.bits, Nearest),
        numerator.im.div_round(&denominator, p.bits, Nearest),
    )
}

#[test]
fn real_operand_shortcuts_use_correct_rounding_at_mixed_precisions() -> Result<()> {
    for bits in [53, 80, 201, 400] {
        let p = Precision { bits };
        for source_bits in [40, 201, 700] {
            let storage = Precision { bits: source_bits };
            let tiny = storage.powi(&storage.i(2), -1000);
            let large = storage.powi(&storage.i(2), 1000);
            let third = storage.rational(&Rational::from((1, 3)));
            let values = [
                storage.zero(),
                storage.i(1),
                storage.i(-3),
                third,
                storage.complex(0, -2),
                storage.complex(2, -3),
                tiny,
                large,
            ];
            for a in &values {
                let expected = if a.im.is_zero() {
                    rounded(p, a.re.to_rational().abs())
                } else {
                    legacy_mul(p, a, &Complex::new(a.re.clone(), -a.im.clone()))
                        .re
                        .sqrt()
                };
                assert_eq!(p.norm(a), expected);
                for b in &values {
                    assert_eq!(p.mul(a, b), legacy_mul(p, a, b));
                    if *b != storage.zero() {
                        let expected = if b.im.is_zero() {
                            Complex::new(
                                rounded(p, a.re.to_rational() / b.re.to_rational()),
                                rounded(p, a.im.to_rational() / b.re.to_rational()),
                            )
                        } else {
                            legacy_div(p, a, b)
                        };
                        let actual = p.div(a, b);
                        if b.im.is_zero() {
                            assert_eq!(p.round(&actual), expected);
                        } else {
                            assert_eq!(actual, expected);
                        }
                    }
                }
                for (n, d) in [(1, 3), (-7, 11), (0, 5)] {
                    let numerator = if n == 1 {
                        a.clone()
                    } else {
                        legacy_mul(p, a, &p.i(n))
                    };
                    let expected = Complex::new(
                        rounded(p, numerator.re.to_rational() / Rational::from(d)),
                        rounded(p, numerator.im.to_rational() / Rational::from(d)),
                    );
                    assert_eq!(p.round(&p.scale(a, n, d)), expected);
                }
            }
        }
    }
    Ok(())
}

#[test]
fn nonfinite_values_keep_general_complex_arithmetic_behavior() -> Result<()> {
    let p = Precision::decimal(40)?;
    for real in ["inf", "-inf", "NaN"] {
        let a = p.parse(real, "0")?;
        let b = p.i(3);
        let actual = p.mul(&a, &b);
        let expected = legacy_mul(p, &a, &b);
        assert_eq!(actual.re.to_string(), expected.re.to_string());
        assert_eq!(actual.im.to_string(), expected.im.to_string());
        let actual = p.div(&a, &b);
        let expected = legacy_div(p, &a, &b);
        assert_eq!(actual.re.to_string(), expected.re.to_string());
        assert_eq!(actual.im.to_string(), expected.im.to_string());
    }
    Ok(())
}
