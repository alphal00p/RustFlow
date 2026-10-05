use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{AlgebraicKinematicSystem, AlgebraicSystem, SquareRoot};
use symbolica_amflow::diffexp::EpsilonSystem;
use symbolica_amflow::transport_cache::*;
use symbolica_amflow::*;

fn symbols() -> (Symbol, Symbol, Symbol) {
    (
        symbol!("algebraic_bank::s"),
        symbol!("algebraic_bank::eps"),
        symbol!("algebraic_bank::r"),
    )
}
fn germ(sheet: RootSheet) -> RootGerm {
    RootGerm {
        sheets: BTreeMap::from([(symbols().2, sheet)]),
    }
}
fn flow(radicand: Atom, conditions: &[Atom]) -> Result<RustFlow<AlgebraicKinematicSystem>> {
    let (s, epsilon, r) = symbols();
    RustFlow::with_algebraic_conditions(
        AlgebraicKinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(s, vec![vec![Atom::var(epsilon) / Atom::var(r)]])]),
            roots: vec![SquareRoot {
                symbol: r,
                radicand,
            }],
        },
        &[parse!("algebraic_bank::I")],
        &Atom::one(),
        Prescription::PlusI0,
        "specified local root sheet in a simply connected regular real component",
        conditions,
    )
}
fn point(value: Atom, sheet: RootSheet) -> Result<CachedPoint> {
    CachedPoint::Exact(BTreeMap::from([(symbols().0, value)])).with_root_germ(germ(sheet))
}
fn boundary(
    flow: &RustFlow<AlgebraicKinematicSystem>,
    value: Atom,
    sheet: RootSheet,
    coefficients: &[i64],
    digits: u32,
) -> Result<CachedBoundary> {
    let p = Precision::decimal(80)?;
    Ok(CachedBoundary {
        identity: flow.identity().clone(),
        point: point(value, sheet)?,
        kind: PointKind::Physical,
        range: EpsilonRange::new(0, coefficients.len() as i32 - 1)?,
        coefficients: coefficients.iter().map(|&a| vec![p.i(a)]).collect(),
        accuracy: BoundaryAccuracy::supplied(
            digits,
            p.bits,
            vec![vec![p.tolerance(60)]; coefficients.len()],
            "independent analytic boundary",
        )?,
    })
}
fn options() -> FlowOptions {
    FlowOptions {
        digits: 20,
        guard_digits: 30,
        series_order: 64,
        ..Default::default()
    }
}
fn policy() -> impl TransportCost {
    ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    }
}

#[test]
fn progressive_algebraic_cache_keeps_germs_and_restarts_binary() -> Result<()> {
    let s = symbols().0;
    let flow = flow(Atom::var(s), &[])?;
    let mut cache = RustFlowCache::default();
    cache.insert(boundary(
        &flow,
        Atom::one(),
        RootSheet::Principal,
        &[1, 0, 0],
        35,
    )?)?;
    let range = EpsilonRange::new(0, 2)?;
    let first = flow.evaluate_to(
        &mut cache,
        &BTreeMap::from([(s, Atom::num(4))]),
        &germ(RootSheet::Principal),
        range,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert!(first.inserted_points > 1);
    assert!(cache.entries().iter().any(|b|matches!(&b.point,CachedPoint::Algebraic {point,..} if matches!(point.as_ref(),CachedPoint::Derived{..}))));
    assert!(
        cache
            .entries()
            .iter()
            .all(|b| b.point.root_germ() == Some(&germ(RootSheet::Principal)))
    );
    let path =
        std::env::temp_dir().join(format!("rustflow-algebraic-cache-{}", std::process::id()));
    cache.save(&path)?;
    let mut restored = RustFlowCache::load(&path)?;
    std::fs::remove_dir_all(&path)?;
    assert_eq!(cache.len(), restored.len());
    let second = flow.evaluate_to(
        &mut restored,
        &BTreeMap::from([(s, Atom::num(9))]),
        &germ(RootSheet::Principal),
        range,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert!((second.starting_point.restart_coordinates()?[&s].clone() - 4).is_zero());
    let p = Precision::decimal(80)?;
    // exp(2 epsilon (sqrt(s)-1)) gives [1,4,8] at s=9.
    for (actual, expected) in second.boundary.coefficients.iter().flatten().zip([1, 4, 8]) {
        assert!(p.close(actual, &p.i(expected), 30));
    }
    let hit = flow.evaluate_to(
        &mut restored,
        &BTreeMap::from([(s, Atom::num(9))]),
        &germ(RootSheet::Principal),
        range,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert!(hit.transport.is_none());
    assert_eq!(hit.inserted_points, 0);
    Ok(())
}

#[test]
fn opposite_root_sheets_coexist_and_exact_hits_do_not_cross() -> Result<()> {
    let s = symbols().0;
    let flow = flow(Atom::var(s), &[])?;
    let mut cache = RustFlowCache::default();
    cache.insert(boundary(
        &flow,
        Atom::num(9),
        RootSheet::Principal,
        &[1, 4, 8],
        35,
    )?)?;
    cache.insert(boundary(
        &flow,
        Atom::num(9),
        RootSheet::Opposite,
        &[1, -4, 8],
        40,
    )?)?;
    assert_eq!(cache.len(), 2);
    for sheet in [RootSheet::Principal, RootSheet::Opposite] {
        let hit = flow.evaluate_to(
            &mut cache,
            &BTreeMap::from([(s, Atom::num(9))]),
            &germ(sheet),
            EpsilonRange::new(0, 2)?,
            &options(),
            &RunContext::default(),
            &policy(),
        )?;
        assert!(hit.transport.is_none());
        let expected = if sheet == RootSheet::Principal { 4 } else { -4 };
        assert_eq!(
            hit.boundary.coefficients[1][0],
            Precision::decimal(80)?.i(expected)
        );
    }
    Ok(())
}

#[test]
fn nearest_source_across_even_root_zero_is_rejected_before_ranking() -> Result<()> {
    let s = symbols().0;
    let flow = flow((Atom::var(s) - 3).pow(2), &[])?;
    let mut cache = RustFlowCache::default();
    cache.insert(boundary(
        &flow,
        Atom::num(Rational::from((29, 10))),
        RootSheet::Principal,
        &[1, 0],
        35,
    )?)?;
    cache.insert(boundary(
        &flow,
        Atom::num(4),
        RootSheet::Principal,
        &[1, 0],
        35,
    )?)?;
    let result = flow.evaluate_to(
        &mut cache,
        &BTreeMap::from([(s, Atom::num(Rational::from((31, 10))))]),
        &germ(RootSheet::Principal),
        EpsilonRange::new(0, 1)?,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert_eq!(
        result.starting_point.restart_coordinates()?[&s],
        Atom::num(4)
    );
    Ok(())
}

#[test]
fn algebraic_cache_preserves_original_reduction_guards_and_input_accuracy() -> Result<()> {
    let s = symbols().0;
    let flow = flow(Atom::var(s), &[Atom::var(s) - 2])?;
    let mut cache = RustFlowCache::default();
    cache.insert(boundary(
        &flow,
        Atom::one(),
        RootSheet::Principal,
        &[1, 0],
        35,
    )?)?;
    let before = cache.len();
    assert!(matches!(
        flow.evaluate_to(
            &mut cache,
            &BTreeMap::from([(s, Atom::num(3))]),
            &germ(RootSheet::Principal),
            EpsilonRange::new(0, 1)?,
            &options(),
            &RunContext::default(),
            &policy()
        ),
        Err(Error::IncompleteReduction(_))
    ));
    assert_eq!(before, cache.len());
    let flow = flow_without_guard()?;
    let mut cache = RustFlowCache::default();
    cache.insert(boundary(
        &flow,
        Atom::one(),
        RootSheet::Principal,
        &[1, 0],
        20,
    )?)?;
    assert!(matches!(
        flow.evaluate_to(
            &mut cache,
            &BTreeMap::from([(s, Atom::num(4))]),
            &germ(RootSheet::Principal),
            EpsilonRange::new(0, 1)?,
            &options(),
            &RunContext::default(),
            &policy()
        ),
        Err(Error::Accuracy(_))
    ));
    assert_eq!(cache.len(), 1);
    Ok(())
}
fn flow_without_guard() -> Result<RustFlow<AlgebraicKinematicSystem>> {
    flow(Atom::var(symbols().0), &[])
}

#[test]
fn algebraic_cache_rejects_foreign_germs_and_branch_points_but_accepts_complex_points() -> Result<()>
{
    let s = symbols().0;
    let flow = flow(Atom::var(s), &[])?;
    let mut wrong = boundary(&flow, Atom::one(), RootSheet::Principal, &[1], 35)?;
    wrong.point = CachedPoint::Exact(BTreeMap::from([(s, Atom::one())]));
    assert!(matches!(wrong.validate(), Err(Error::InvalidInput(_))));
    wrong.point =
        CachedPoint::Exact(BTreeMap::from([(s, Atom::one())])).with_root_germ(RootGerm {
            sheets: BTreeMap::from([(symbol!("foreign_root"), RootSheet::Principal)]),
        })?;
    assert!(matches!(wrong.validate(), Err(Error::InvalidInput(_))));
    wrong.point = point(Atom::new(), RootSheet::Principal)?;
    assert!(matches!(wrong.validate(), Err(Error::InvalidInput(_))));
    wrong.point = point(
        Atom::num(Complex::new(Rational::from(1), Rational::from(1))),
        RootSheet::Principal,
    )?;
    wrong.validate()?;
    assert_ne!(
        flow.identity().key(),
        flow_without_guard_squared()?.identity().key()
    );
    Ok(())
}
fn flow_without_guard_squared() -> Result<RustFlow<AlgebraicKinematicSystem>> {
    flow(Atom::var(symbols().0).pow(2), &[])
}

#[test]
fn registered_root_weighted_amplification_matches_constant_analytic_norm() -> Result<()> {
    let p = Precision::decimal(70)?;
    let x = symbol!("algebraic_bound::x");
    let r = symbol!("algebraic_bound::r");
    let compiled = AlgebraicSystem {
        system: EpsilonSystem {
            variable: x,
            matrices: vec![
                vec![
                    vec![Atom::one() / Atom::var(r), Atom::new()],
                    vec![Atom::new(), Atom::num(2) / Atom::var(r)],
                ],
                vec![
                    vec![Atom::new(), Atom::new()],
                    vec![Atom::num(3) / Atom::var(r), Atom::new()],
                ],
            ],
        },
        roots: vec![SquareRoot {
            symbol: r,
            radicand: Atom::num(4),
        }],
        nonzero_conditions: Vec::new(),
    }
    .compile(p)?;
    let amplification = compiled.error_amplification_weighted(
        &p.zero(),
        &p.rational(&Rational::from((1, 100))),
        &[p.real(100), p.real(1), p.real(1), p.real(2)],
    )?;
    // epsilon1,row1 has weighted norm 1 + (3/2)*100/2 = 76.
    assert!(p.close(
        &ComplexFloat::new(amplification, p.real(0)),
        &p.exp(&p.rational(&Rational::from((76, 100)))),
        60
    ));
    Ok(())
}

#[test]
fn negative_real_radicand_retains_its_imaginary_physical_sheet() -> Result<()> {
    let s = symbols().0;
    let flow = flow(-Atom::var(s), &[])?;
    let mut cache = RustFlowCache::default();
    cache.insert(boundary(
        &flow,
        Atom::one(),
        RootSheet::Principal,
        &[1, 0, 0],
        35,
    )?)?;
    let result = flow.evaluate_to(
        &mut cache,
        &BTreeMap::from([(s, Atom::num(9))]),
        &germ(RootSheet::Principal),
        EpsilonRange::new(0, 2)?,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    let p = Precision::decimal(80)?;
    for (actual, expected) in
        result
            .boundary
            .coefficients
            .iter()
            .flatten()
            .zip([p.i(1), p.complex(0, -4), p.i(-8)])
    {
        assert!(p.close(actual, &expected, 30));
    }
    Ok(())
}

#[test]
fn real_endpoint_radicands_admit_a_regular_complex_radicand_inside_the_path() -> Result<()> {
    let s = symbols().0;
    let imaginary = Atom::num(Complex::new(Rational::from(0), Rational::from(1)));
    let flow = flow(
        Atom::one() + imaginary * Atom::var(s) * (Atom::var(s) - 1),
        &[],
    )?;
    assert!(flow.identity().conditions_admit_straight_path(
        &point(Atom::new(), RootSheet::Principal)?,
        &point(Atom::one(), RootSheet::Principal)?,
        Precision::decimal(60)?,
        20
    )?);
    Ok(())
}

#[test]
fn cache_uses_registered_root_sum_norms_and_generic_epsilon_domains() -> Result<()> {
    let (s, epsilon, r) = symbols();
    let system = AlgebraicKinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([(s, vec![vec![Atom::var(epsilon) / (Atom::var(r) + 1)]])]),
        roots: vec![SquareRoot {
            symbol: r,
            radicand: Atom::var(s),
        }],
    };
    let flow = RustFlow::new_algebraic(
        system,
        &[parse!("algebraic_bank::I")],
        &Atom::one(),
        Prescription::PlusI0,
        "regular real root-sum domain",
    )?;
    assert!(matches!(
        boundary(&flow, Atom::one(), RootSheet::Principal, &[1, 0, 0], 35)?.validate(),
        Err(Error::InvalidInput(_))
    ));
    let mut cache = RustFlowCache::default();
    cache.insert(boundary(
        &flow,
        Atom::num(4),
        RootSheet::Principal,
        &[1, 0, 0],
        35,
    )?)?;
    let result = flow.evaluate_to(
        &mut cache,
        &BTreeMap::from([(s, Atom::num(9))]),
        &germ(RootSheet::Principal),
        EpsilonRange::new(0, 2)?,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    let p = Precision::decimal(80)?;
    let value = p.sub(
        &p.i(2),
        &p.scale(&p.log(&p.rational(&Rational::from((4, 3)))), 2, 1),
    );
    assert!(p.close(&result.boundary.coefficients[1][0], &value, 30));
    assert!(p.close(
        &result.boundary.coefficients[2][0],
        &p.scale(&p.mul(&value, &value), 1, 2),
        30
    ));
    let valuation = RustFlow::new_algebraic(
        AlgebraicKinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(
                s,
                vec![vec![Atom::one() / (Atom::var(r) + 1 + Atom::var(epsilon))]],
            )]),
            roots: vec![SquareRoot {
                symbol: r,
                radicand: Atom::var(s) + 1,
            }],
        },
        &[parse!("algebraic_bank::I")],
        &Atom::one(),
        Prescription::PlusI0,
        "fixed generic epsilon valuation",
    )?;
    // Norm epsilon^2 + 2 epsilon - s has generic leading coefficient -s.
    // Its specialization at s=0 is nonzero as an epsilon polynomial, but the
    // Laurent structure has changed and must not become an exact cache hit.
    assert!(matches!(
        boundary(&valuation, Atom::new(), RootSheet::Principal, &[1], 35)?.validate(),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}

#[test]
fn internal_path_parameters_do_not_collide_with_registered_roots() -> Result<()> {
    let (s, epsilon, _) = symbols();
    let a = symbol!("symbolica_amflow::physical_path_parameter_0");
    let b = symbol!("symbolica_amflow::guard_path_0");
    let system = AlgebraicKinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([(
            s,
            vec![vec![
                Atom::var(epsilon) * (Atom::one() / Atom::var(a) + Atom::one() / Atom::var(b)),
            ]],
        )]),
        roots: vec![
            SquareRoot {
                symbol: a,
                radicand: Atom::var(s),
            },
            SquareRoot {
                symbol: b,
                radicand: Atom::var(s),
            },
        ],
    };
    let flow = RustFlow::new_algebraic(
        system,
        &[parse!("algebraic_bank::I")],
        &Atom::one(),
        Prescription::PlusI0,
        "two principal positive roots",
    )?;
    let germ = RootGerm {
        sheets: BTreeMap::from([(a, RootSheet::Principal), (b, RootSheet::Principal)]),
    };
    let p = Precision::decimal(80)?;
    let range = EpsilonRange::new(0, 1)?;
    let mut cache = RustFlowCache::default();
    cache.insert(CachedBoundary {
        identity: flow.identity().clone(),
        point: CachedPoint::Exact(BTreeMap::from([(s, Atom::one())]))
            .with_root_germ(germ.clone())?,
        kind: PointKind::Physical,
        range,
        coefficients: vec![vec![p.i(1)], vec![p.zero()]],
        accuracy: BoundaryAccuracy::supplied(
            35,
            p.bits,
            vec![vec![p.tolerance(60)]; 2],
            "analytic independent source",
        )?,
    })?;
    let result = flow.evaluate_to(
        &mut cache,
        &BTreeMap::from([(s, Atom::num(4))]),
        &germ,
        range,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert!(p.close(&result.boundary.coefficients[1][0], &p.i(4), 30));
    Ok(())
}
