use super::*;

#[test]
fn retained_guard_roots_detour_even_when_constant_matrix_has_no_poles() -> Result<()> {
    let variable = symbol!("retained_guard::eta");
    let eta = Atom::var(variable);
    let p = Precision::decimal(60)?;
    let mut compiled = DifferentialSystem {
        variable,
        matrix: vec![vec![Atom::new()]],
    }
    .compile(p, &Default::default())?;
    assert!(compiled.poles.is_empty());
    let imaginary = Atom::num(Complex::new(Rational::from(0), Rational::from(1)));
    let guards = crate::physical_conditions::canonical_conditions(
        &[(eta.clone().pow(2) + Atom::one()) / (&eta + imaginary)],
        &std::collections::BTreeSet::from([variable]),
    )?;
    compiled.exclude_polynomials(&guards)?;
    assert!(
        compiled
            .poles
            .iter()
            .any(|pole| p.close(pole, &p.complex(0, -1), 50))
    );
    assert!(
        compiled
            .poles
            .iter()
            .any(|pole| p.close(pole, &p.complex(0, 1), 50))
    );
    let path = compiled.plan_path(&p.complex(0, -8), &p.parse("0", "-0.125")?, 1)?;
    assert!(path.len() > 1, "retained -i must force a detour");
    assert!(path.iter().any(|point| point.re != p.real(0)));
    Ok(())
}

#[test]
fn retained_guard_survives_zero_jacobian_product_in_mobius_chart() -> Result<()> {
    let variable = symbol!("retained_map_guard::eta");
    let p = Precision::decimal(60)?;
    let mut compiled = DifferentialSystem {
        variable,
        matrix: vec![vec![Atom::new()]],
    }
    .compile(p, &Default::default())?;
    compiled.exclude_polynomials(&[Atom::var(variable) - Atom::num(2)])?;
    let coordinate = TaylorCoordinate::balanced(p, &p.zero(), &p.i(1), &compiled.poles)?;
    assert!(matches!(coordinate, TaylorCoordinate::Mobius(_)));
    // eta = 4 y/(1+y): retained eta=2 becomes y=1, while y=-1 is
    // the map's own excluded point. Multiplying the zero matrix by d eta/dy
    // must not erase either of them.
    let mapped = compiled.in_coordinate(&coordinate)?;
    assert_eq!(mapped.source.rows, vec![vec![Atom::new()]]);
    for pole in [p.i(1), p.i(-1)] {
        assert!(mapped.poles.iter().any(|found| p.close(found, &pole, 50)));
    }
    assert!(mapped.pole_polynomials.len() >= 2);
    Ok(())
}
