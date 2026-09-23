use symbolica::prelude::*;
use symbolica_amflow::*;
fn sunset() -> IntegralFamily {
    IntegralFamily {
        name: "sunset".into(),
        loops: vec!["l1".into(), "l2".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![
            Propagator::quadratic(&[1, 0], &[], Atom::num(1), &[]).unwrap(),
            Propagator::quadratic(&[0, 1], &[], Atom::new(), &[]).unwrap(),
            Propagator::quadratic(&[1, 1], &[], Atom::new(), &[]).unwrap(),
        ],
        physical_propagators: 3,
        epsilon: symbol!("eps"),
        dimension: 4,
    }
}
#[test]
fn sunset_has_five_regions_including_alternative_routings() {
    let family = sunset();
    let regions = regions::enumerate_regions(&family, 1000).unwrap();
    assert_eq!(regions.len(), 5);
    for signature in [
        vec![false, false, false],
        vec![true, true, true],
        vec![true, true, false],
        vec![true, false, true],
        vec![false, true, true],
    ] {
        assert!(regions.iter().any(|r| r.hard_branches == signature));
    }
    let hard = regions.iter().find(|r| r.hard.iter().all(|&v| v)).unwrap();
    let expansion = regions::expand_region(
        &family,
        &Integral(vec![1, 1, 1]),
        &[true, false, false],
        hard,
        2,
    )
    .unwrap();
    assert!(
        (&expansion.eta_power - parse!("1-2*eps"))
            .together()
            .cancel()
            .is_zero()
    );
    assert!(expansion.coefficients[1].is_zero());
    assert!(!expansion.coefficients[0].is_zero());
    let soft = regions.iter().find(|r| r.hard.iter().all(|&v| !v)).unwrap();
    let expansion = regions::expand_region(
        &family,
        &Integral(vec![1, 1, 1]),
        &[true, false, false],
        soft,
        2,
    )
    .unwrap();
    assert_eq!(expansion.eta_power, Atom::num(-1));
}
