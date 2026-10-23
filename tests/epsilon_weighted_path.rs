//! Downstream regressions for weighted error bounds along explicit paths.
use symbolica_amflow::{diffexp::EpsilonSystem, symbolica::prelude::*, *};

#[test]
fn nonzero_epsilon_zero_block_path_composition_matches_constant_kernel() -> Result<()> {
    let p = Precision::decimal(70)?;
    let system = EpsilonSystem {
        variable: symbol!("weighted_path_api::x"),
        matrices: vec![vec![vec![Atom::one()]], vec![vec![Atom::num(2)]]],
    }
    .compile(p, &Default::default())?;
    let weights = vec![p.real(1); 2];
    let points = vec![
        p.zero(),
        p.rational(&Rational::from((1, 8))),
        p.rational(&Rational::from((1, 2))),
        p.i(1),
    ];
    // The weighted row norm is exactly 3. Nonzero A0 requires exp(3),
    // rather than the degree-one polynomial 1+3 of a two-coefficient hierarchy.
    let expected = p.exp(&p.i(3));
    let actual = system.error_amplification_weighted_path(&points, &weights)?;
    let whole = system.error_amplification_weighted(&points[0], &points[3], &weights)?;
    let mut product = p.i(1);
    for leg in points.windows(2) {
        let factor = system.error_amplification_weighted(&leg[0], &leg[1], &weights)?;
        product = p.mul(&product, &ComplexFloat::new(factor, p.real(0)));
    }
    assert!(p.close(&ComplexFloat::new(actual.clone(), p.real(0)), &expected, 60));
    assert!(p.close(&ComplexFloat::new(whole, p.real(0)), &expected, 60));
    assert!(p.close(&product, &expected, 60));
    assert!(actual > p.real(4));
    Ok(())
}

#[test]
fn weighted_path_rejects_too_short_paths_and_bad_weights() -> Result<()> {
    let p = Precision::decimal(60)?;
    let system = EpsilonSystem {
        variable: symbol!("weighted_path_api::invalid_x"),
        matrices: vec![vec![vec![Atom::one()]]],
    }
    .compile(p, &Default::default())?;
    let weights = vec![p.real(1)];
    for points in [Vec::new(), vec![p.zero()]] {
        assert!(matches!(
            system.error_amplification_weighted_path(&points, &weights),
            Err(Error::InvalidInput(_))
        ));
    }
    let points = vec![p.zero(), p.i(1)];
    for bad in [
        Vec::new(),
        vec![p.real(0)],
        vec![p.real(-1)],
        vec![p.real(1), p.real(1)],
    ] {
        assert!(matches!(
            system.error_amplification_weighted_path(&points, &bad),
            Err(Error::InvalidInput(_))
        ));
    }
    assert_eq!(
        system.error_amplification_weighted_path(&[p.zero(), p.zero()], &weights)?,
        p.real(1)
    );
    Ok(())
}
