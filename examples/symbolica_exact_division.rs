//! Symbolica-only reproducer for exact AtomField division ignoring its zero-test setting.
//!
//! Run `cargo run --release --example symbolica_exact_division -- scalar` or `matrix`.
//! Both fail at Symbolica 942bd2c. `normalized` and `polynomial` are passing controls.
use symbolica::prelude::*;
use symbolica::tensors::matrix::Matrix;

fn exact_field() -> AtomField {
    AtomField {
        statistical_zero_test: false,
        cancel_check_on_division: true,
        custom_normalization: Some(Box::new(|a: AtomView<'_>, out: &mut Atom| {
            *out = a.together().cancel();
            true
        })),
    }
}

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_else(|| "scalar".into());
    // Exactly x + 1 - y, represented with a large unexpanded integer factor.
    let denominator = parse!("x+(10^40-10^40*y)/10^40");
    let field = exact_field();
    assert!(!field.is_zero(&denominator));
    assert!(
        (denominator.expand() - parse!("x+1-y"))
            .together()
            .cancel()
            .is_zero()
    );

    match mode.as_str() {
        "scalar" => {
            let quotient = field.try_div(&Atom::num(1), &denominator);
            let inverse = field.try_inv(&denominator);
            println!("exact nonzero denominator: {denominator}");
            println!("try_div succeeded: {}", quotient.is_some());
            println!("try_inv succeeded: {}", inverse.is_some());
            for value in [quotient, inverse] {
                let value = value.expect("exact nonzero denominator must remain invertible");
                assert!(
                    (value * &denominator - Atom::num(1))
                        .together()
                        .cancel()
                        .is_zero()
                );
            }
        }
        "matrix" | "normalized" | "polynomial" => {
            // Four is the smallest size that uses Matrix::det's Bareiss path.
            let mut matrix = Matrix::identity(4, field);
            matrix[(0, 0)] = if mode == "normalized" {
                denominator.expand()
            } else {
                denominator.clone()
            };
            let determinant = if mode == "polynomial" {
                use symbolica::domains::rational_polynomial::RationalPolynomialField;
                let matrix = matrix.map(
                    |entry| entry.try_to_rational_polynomial(&Q, &Z, None).unwrap(),
                    RationalPolynomialField::<_, u16>::new(Z),
                );
                matrix.det().unwrap().to_expression()
            } else {
                matrix.det().unwrap()
            };
            assert!((&determinant - &denominator).together().cancel().is_zero());
            println!("determinant = {determinant}");
        }
        _ => panic!("expected scalar, matrix, normalized, or polynomial"),
    }
}
