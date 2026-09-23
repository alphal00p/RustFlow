use symbolica::prelude::*;
use symbolica_amflow::*;
fn template() -> IntegralFamily {
    IntegralFamily {
        name: "boundary".into(),
        loops: vec!["l".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![],
        physical_propagators: 0,
        epsilon: symbol!("eps"),
        dimension: 4,
    }
}
#[test]
fn partial_fractions_and_numerators_reconstruct_exactly() {
    let x = Atom::var(symbol!("x"));
    let expr = parse!("x^3/((x-1)^2*(x-2))");
    let terms =
        integrand::to_integrals(&expr, std::slice::from_ref(&x), &template(), 1000).unwrap();
    let mut rebuilt = Atom::new();
    for term in terms {
        term.family.validate().unwrap();
        let mut v = term.coefficient;
        for (d, &n) in term.family.propagators.iter().zip(&term.integral.0) {
            v *= (&d.constant + &d.scalar_products[0] * &x).pow(-(n as i64));
        }
        rebuilt += v;
    }
    assert!((rebuilt - expr).together().cancel().is_zero());
}
#[test]
fn mixed_tensor_boundary_factorizes() {
    let coordinates = vec![parse!("hh"), parse!("hs"), parse!("ss")];
    let mut family = template();
    family.loops.push("q".into());
    let terms = integrand::factor_region(
        &parse!("hs^2/((hh-1)^3*(ss-2)^2)"),
        &coordinates,
        &family,
        &[true, false],
        1000,
    )
    .unwrap();
    assert!(!terms.is_empty());
    assert!(
        terms
            .iter()
            .all(|t| t.factors.len() == 2 && t.factors.iter().all(|f| f.family.loops.len() == 1))
    );
    let d = parse!("4-2*eps");
    let mut rebuilt = Atom::new();
    for term in terms {
        let mut v = term.coefficient;
        for (i, factor) in term.factors.into_iter().enumerate() {
            let x = &coordinates[i * 2];
            for (p, n) in factor.family.propagators.iter().zip(factor.integral.0) {
                v *= (&p.constant + &p.scalar_products[0] * x).pow(-(n as i64));
            }
        }
        rebuilt += v;
    }
    assert!(
        (rebuilt - parse!("hh*ss/((hh-1)^3*(ss-2)^2)") / d)
            .together()
            .cancel()
            .is_zero()
    );
}
