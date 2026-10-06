#![cfg(feature = "wasm")]
use symbolica::prelude::*;
use symbolica_amflow::{Error, FlowOptions, Precision};

#[test]
fn portable_constants_and_arithmetic_use_the_existing_numeric_owner() {
    for digits in [30, 80, 140] {
        let p = Precision::decimal(digits).unwrap();
        let zero = p.real(0);
        let pi = zero.pi();
        let gamma = zero.euler();
        assert_eq!(pi.prec(), p.bits);
        assert_eq!(gamma.prec(), p.bits);
        let through_exponential = p.exp(&Complex::new(zero, pi));
        assert!(p.close(&through_exponential, &p.i(-1), digits));
        let reference = p
            .parse("0.57721566490153286060651209008240243104215933593992", "0")
            .unwrap();
        assert!(p.close(&Complex::new(gamma, p.real(0)), &reference, digits.min(45)));
        assert!(matches!(
            p.gamma_real(&p.real(3)),
            Err(Error::Unsupported(_))
        ));
    }
}

#[test]
fn portable_transport_rejects_multiple_workers_explicitly() {
    FlowOptions::default().validate().unwrap();
    let options = FlowOptions {
        workers: 2,
        ..FlowOptions::default()
    };
    assert!(matches!(options.validate(), Err(Error::Unsupported(_))));
}
