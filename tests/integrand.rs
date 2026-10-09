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
fn large_mixed_rank_is_bounded_before_gram_allocation() {
    let coordinates = vec![
        parse!("large_rank_hh"),
        parse!("large_rank_hs"),
        parse!("large_rank_ss"),
    ];
    let mut family = template();
    family.loops.push("q".into());
    // Without the preflight bound this tiny monomial requests two matrices
    // with more than a billion Atom entries each before the native rejection.
    let even = coordinates[1].clone().pow(32766);
    assert!(matches!(
        integrand::projected_factor_region(&even, &coordinates, &family, &[true, false], 10),
        Err(Error::Limit(message)) if message.contains("before Gram construction")
    ));
    // Odd vacuum tensors still vanish exactly without allocating those arrays.
    let odd = coordinates[1].clone().pow(32767);
    let projected =
        integrand::projected_factor_region(&odd, &coordinates, &family, &[true, false], 10)
            .unwrap();
    assert!(
        projected
            .terms
            .iter()
            .fold(Atom::zero(), |value, term| {
                value + &term.hard * &term.soft
            })
            .is_zero()
    );
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

#[test]
fn projected_region_preserves_soft_denominators_and_coordinate_maps() {
    let mut family = template();
    family.loops.extend(["q".into(), "r".into()]);
    family.external = vec!["u".into()];
    family.external_gram = vec![vec![Atom::num(1)]];
    let coordinates = vec![
        parse!("q00"),
        parse!("q01"),
        parse!("q02"),
        parse!("q11"),
        parse!("q12"),
        parse!("q22"),
        parse!("q0u"),
        parse!("q1u"),
        parse!("q2u"),
    ];
    // The soft factor is nonpolynomial and its first denominator is not an
    // ordinary affine propagator. Projection must preserve it without invoking
    // ordinary partial fractions or treating it as a compact polynomial moment.
    let expression = parse!("q01*q0u/((q00-1)^3*(q11^2+1)*(q1u+3))");
    let projected = integrand::projected_factor_region(
        &expression,
        &coordinates,
        &family,
        &[true, false, true],
        10,
    )
    .unwrap();
    assert_eq!(projected.hard.source_loop_indices, [0, 2]);
    assert_eq!(projected.hard.source_coordinate_indices, [0, 2, 5]);
    assert_eq!(
        projected.hard.coordinates,
        [parse!("q00"), parse!("q02"), parse!("q22")]
    );
    assert_eq!(projected.hard.template.loops, ["l", "r"]);
    assert!(projected.hard.template.external.is_empty());
    assert_eq!(projected.soft.source_loop_indices, [1]);
    assert_eq!(projected.soft.source_coordinate_indices, [3, 7]);
    assert_eq!(projected.soft.coordinates, [parse!("q11"), parse!("q1u")]);
    assert_eq!(projected.soft.template.loops, ["q"]);
    assert_eq!(projected.soft.template.external, ["u"]);
    assert!(projected.soft.template.propagators.is_empty());
    let rebuilt = projected
        .terms
        .iter()
        .fold(Atom::new(), |sum, term| sum + &term.hard * &term.soft);
    let expected = parse!("q00*q1u/((4-2*eps)*(q00-1)^3*(q11^2+1)*(q1u+3))");
    assert!((rebuilt - expected).together().cancel().is_zero());
    assert!(matches!(
        integrand::factor_region(&expression, &coordinates, &family, &[true, false, true], 10),
        Err(Error::Unsupported(message)) if message.contains("nonaffine boundary denominator")
    ));
}

#[test]
fn projected_region_bounds_products_and_keeps_unmixed_sides() {
    let mut family = template();
    family.loops.push("q".into());
    let coordinates = vec![parse!("hh"), parse!("hs"), parse!("ss")];
    let expression = parse!("(hh+ss)/((hh-1)*(ss-2))");
    assert!(matches!(
        integrand::projected_factor_region(&expression, &coordinates, &family, &[true, false], 1),
        Err(Error::Limit(_))
    ));
    for hard in [[true, false], [true, true], [false, false]] {
        let projected =
            integrand::projected_factor_region(&expression, &coordinates, &family, &hard, 2)
                .unwrap();
        let rebuilt = projected
            .terms
            .iter()
            .fold(Atom::new(), |sum, term| sum + &term.hard * &term.soft);
        assert!((rebuilt - &expression).together().cancel().is_zero());
        if hard == [true, true] {
            assert!(projected.soft.source_loop_indices.is_empty());
            assert!(projected.soft.coordinates.is_empty());
        }
        if hard == [false, false] {
            assert!(projected.hard.source_loop_indices.is_empty());
            assert!(projected.hard.coordinates.is_empty());
        }
    }
}

#[test]
fn factorized_budget_counts_distinct_integral_products() {
    let coordinates = vec![parse!("hh"), parse!("hs"), parse!("ss")];
    let mut family = template();
    family.loops.push("q".into());
    let expression = parse!("(hh+ss)^12/((hh-1)^15*(ss-2)^15)");
    // Monomial-by-monomial conversion produces 455 contributions; only 91
    // different integral products remain after exact collection.
    for hard in [[true, false], [true, true], [false, false]] {
        let terms =
            integrand::factor_region(&expression, &coordinates, &family, &hard, 100).unwrap();
        assert_eq!(terms.len(), 91);
        let mut reconstructed = Atom::new();
        for term in terms {
            let mut value = term.coefficient;
            for (side, factor) in term.factors.into_iter().enumerate() {
                let factor_coordinates = if hard == [true, false] {
                    vec![coordinates[side * 2].clone()]
                } else {
                    coordinates.clone()
                };
                for (d, n) in factor.family.propagators.iter().zip(factor.integral.0) {
                    let denominator = d
                        .scalar_products
                        .iter()
                        .zip(&factor_coordinates)
                        .fold(d.constant.clone(), |a, (c, x)| a + c * x);
                    value *= denominator.pow(-i64::from(n));
                }
            }
            reconstructed += value;
        }
        assert!((reconstructed - &expression).together().cancel().is_zero());
        assert!(matches!(
            integrand::factor_region(&expression, &coordinates, &family, &hard, 80),
            Err(Error::Limit(message)) if message.contains("91 distinct terms exceed 80")
        ));
    }
}

#[test]
fn paper_boundary_all_hard_order_eight_fits_the_term_budget() {
    let (paper, targets) = benchmarks::paper_two_loop().unwrap();
    let (_, shifted) = paper
        .deform(symbol!("budget_eta"), &MassMode::Auto)
        .unwrap();
    let regions = regions::enumerate_regions(&paper, 10000).unwrap();
    let all_soft = regions.iter().find(|r| r.hard.iter().all(|v| !v)).unwrap();
    let leading = regions::expand_region(&paper, &targets[0], &shifted, all_soft, 0).unwrap();
    let factors = integrand::factor_region(
        &leading.coefficients[0],
        &leading.coordinates,
        &paper,
        &all_soft.hard,
        10000,
    )
    .unwrap();
    let family = &factors[0].factors[0].family;
    // The three-line subtopology from the actual six-line all-soft boundary:
    // (l1+l2)^2, (l1+p1+p2)^2, and (l2-p1-p2)^2.
    let integral = Integral(
        family
            .propagators
            .iter()
            .map(|d| i16::from(d.constant == Atom::num(30) || d.scalar_products[1] == Atom::num(2)))
            .collect(),
    );
    assert_eq!(integral.0.iter().sum::<i16>(), 3);
    let (_, shifted) = family
        .deform(symbol!("budget_eta"), &MassMode::All)
        .unwrap();
    let regions = regions::enumerate_regions(family, 10000).unwrap();
    let all_hard = regions.iter().find(|r| r.hard.iter().all(|v| *v)).unwrap();
    let expansion = regions::expand_region(family, &integral, &shifted, all_hard, 8).unwrap();
    // The previous algorithm generated over 10,000 duplicate contributions here.
    let terms = integrand::factor_region(
        &expansion.coefficients[8],
        &expansion.coordinates,
        family,
        &all_hard.hard,
        10000,
    )
    .unwrap();
    assert_eq!(terms.len(), 75);
    assert!(terms.iter().all(|t| t.factors.len() == 1));
}
