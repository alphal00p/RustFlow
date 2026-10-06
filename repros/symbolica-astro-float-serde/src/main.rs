use symbolica::domains::float::RoundingDirection;
use symbolica::prelude::*;

fn main() {
    let value =
        Float::from_rational_round(&Rational::from((1, 3)), 415, RoundingDirection::Nearest);
    let config = bincode::config::standard();
    let bytes = bincode::encode_to_vec(&value, config).unwrap();
    let (binary, used): (Float, usize) = bincode::decode_from_slice(&bytes, config).unwrap();
    assert_eq!(used, bytes.len());
    assert_eq!(binary.prec(), value.prec());
    assert_eq!(binary.to_rational(), value.to_rational());

    let bytes = bincode::serde::encode_to_vec(&value, config).unwrap();
    let (serde, used): (Float, usize) = bincode::serde::decode_from_slice(&bytes, config).unwrap();
    assert_eq!(used, bytes.len());
    assert_eq!(serde.prec(), value.prec());
    println!(
        "precision={}, native_binary_exact=true, serde_exact={}",
        value.prec(),
        serde.to_rational() == value.to_rational()
    );
    assert_eq!(
        serde.to_rational(),
        value.to_rational(),
        "Astro serde changed the exact numerical value"
    );
}
