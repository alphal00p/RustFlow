use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::{algebraic::*, transport_cache::*, *};

fn a(text: &str) -> Atom {
    Atom::parse(text, "epsilon_algebraic", Default::default()).unwrap()
}
fn sym(text: &str) -> Symbol {
    let value = a(text);
    let AtomView::Var(v) = value.as_view() else {
        panic!()
    };
    v.get_symbol()
}
fn roots() -> Vec<SquareRoot> {
    vec![SquareRoot {
        symbol: sym("r"),
        radicand: a("s"),
    }]
}
fn system(entry: Atom, lower: Atom) -> AlgebraicKinematicSystem {
    AlgebraicKinematicSystem {
        epsilon: sym("eps"),
        derivatives: BTreeMap::from([(
            sym("s"),
            vec![vec![Atom::new(), entry], vec![lower, Atom::new()]],
        )]),
        roots: roots(),
    }
}
fn flow(system: AlgebraicKinematicSystem) -> Result<RustFlow<AlgebraicKinematicSystem>> {
    RustFlow::new_algebraic(
        system,
        &[a("I(1)"), a("I(2)")],
        &Atom::one(),
        Prescription::PlusI0,
        "specified analytic continuation with named roots",
    )
}
fn germ(sheet: RootSheet) -> RootGerm {
    RootGerm {
        sheets: BTreeMap::from([(sym("r"), sheet)]),
    }
}
fn point(value: Atom) -> BTreeMap<Symbol, Atom> {
    BTreeMap::from([(sym("s"), value)])
}
fn seed(
    flow: &RustFlow<AlgebraicKinematicSystem>,
    at: Atom,
    sheet: RootSheet,
    last: i32,
) -> Result<CachedBoundary> {
    let p = Precision::decimal(90)?;
    Ok(CachedBoundary {
        identity: flow.identity().clone(),
        point: CachedPoint::Exact(point(at)).with_root_germ(germ(sheet))?,
        kind: PointKind::Physical,
        range: EpsilonRange::new(-1, last)?,
        coefficients: (-1..=last)
            .map(|k| vec![p.zero(), p.i(i64::from(k == 0 || k == 1))])
            .collect(),
        accuracy: BoundaryAccuracy::supplied(
            45,
            p.bits,
            vec![vec![p.tolerance(65); 2]; (last + 2) as usize],
            "Analytic source (0,1+epsilon), justified global original pole bound -1",
        )?,
    })
}
fn options() -> FlowOptions {
    FlowOptions {
        digits: 25,
        guard_digits: 45,
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
fn quotient_relations_precede_epsilon_valuation_and_expansion() -> Result<()> {
    let input = system(a("1/(r^2-s+eps)"), a("eps"));
    let (normal, gauge) = EpsilonShearing::regularize_algebraic(&input, &RunContext::default())?;
    assert_eq!(gauge.weights(), &[-1, 0]);
    assert_eq!(
        normal.derivatives[&sym("s")],
        vec![
            vec![Atom::new(), Atom::one()],
            vec![Atom::one(), Atom::new()]
        ]
    );
    let path = kinematics::KinematicPath {
        parameter: sym("x"),
        coordinates: BTreeMap::from([(sym("s"), a("1+x"))]),
    };
    normal
        .pullback(&path, 2)?
        .compile(Precision::decimal(60)?)?;
    let input = system(a("(r^2-s)/eps+1"), a("eps"));
    let (normal, gauge) = EpsilonShearing::regularize_algebraic(&input, &RunContext::default())?;
    assert_eq!(gauge.weights(), &[0, 0]);
    assert_eq!(normal.derivatives[&sym("s")][0][1], Atom::one());
    Ok(())
}

#[test]
fn native_gaussian_coefficients_and_named_root_relations_stay_exact() -> Result<()> {
    let i = Atom::num(Complex::new(Rational::from(0), Rational::from(1)));
    let input = system(
        (Atom::one() + &i) * a("r^-3/eps") + a("(r^2-s)/eps^7"),
        a("eps*r"),
    );
    let (normal, gauge) = EpsilonShearing::regularize_algebraic(&input, &RunContext::default())?;
    assert_eq!(gauge.weights(), &[-1, 0]);
    assert!(
        (normal.derivatives[&sym("s")][0][1].clone() - (Atom::one() + i) * a("r/s^2"))
            .together()
            .cancel()
            .is_zero()
    );
    let mut nonunit = system(a("1/(eps*(r-q))"), a("eps"));
    nonunit.roots.push(SquareRoot {
        symbol: sym("q"),
        radicand: a("s"),
    });
    let nonunit_result = EpsilonShearing::regularize_algebraic(&nonunit, &RunContext::default());
    assert!(
        matches!(
            &nonunit_result,
            Err(Error::Unsupported(_) | Error::InvalidInput(_))
        ),
        "{nonunit_result:?}"
    );
    let mut independent = system(a("(r-q)/eps"), a("eps"));
    independent.roots.push(SquareRoot {
        symbol: sym("q"),
        radicand: a("s+1"),
    });
    let (_, gauge) = EpsilonShearing::regularize_algebraic(&independent, &RunContext::default())?;
    assert_eq!(gauge.weights(), &[-1, 0]);
    Ok(())
}

#[test]
fn unsupported_root_gauges_and_common_partial_constraints_are_explicit() -> Result<()> {
    let mut input = system(a("1/(eps*r)"), a("eps*r"));
    input.derivatives.insert(
        sym("t"),
        vec![
            vec![Atom::new(), a("-1/(eps*r)")],
            vec![Atom::new(), Atom::new()],
        ],
    );
    let (_, gauge) = EpsilonShearing::regularize_algebraic(&input, &RunContext::default())?;
    assert_eq!(gauge.weights(), &[-1, 0]);
    let mut cycle = system(a("1/(eps*r)"), a("1/r"));
    assert!(matches!(
        EpsilonShearing::regularize_algebraic(&cycle, &RunContext::default()),
        Err(Error::Unsupported(_))
    ));
    cycle.roots[0].radicand = a("s+eps");
    assert!(EpsilonShearing::regularize_algebraic(&cycle, &RunContext::default()).is_err());
    assert!(
        EpsilonShearing::regularize_algebraic(
            &system(a("1/(r^2-s)"), a("eps")),
            &RunContext::default()
        )
        .is_err()
    );
    let context = RunContext::default();
    context.cancellation.cancel();
    assert!(EpsilonShearing::regularize_algebraic(&input, &context).is_err());
    Ok(())
}

#[test]
fn original_root_holes_and_generic_norm_domains_survive_cache_conversion() -> Result<()> {
    let original = flow(system(a("(r^2-1)/(eps*(r-1))"), a("eps")))?;
    let gauge = original.regularize_epsilon(&RunContext::default())?;
    assert!(
        gauge
            .to_sheared_boundary(&seed(&original, a("1"), RootSheet::Principal, 1)?, 1)
            .is_err()
    );
    let valid = gauge.to_sheared_boundary(&seed(&original, a("2"), RootSheet::Principal, 1)?, 1)?;
    let mut bank = RustFlowCache::default();
    bank.insert(valid)?;
    assert!(
        gauge
            .flow()
            .evaluate_to(
                &mut bank,
                &point(a("1")),
                &germ(RootSheet::Principal),
                EpsilonRange::new(-1, 1)?,
                &options(),
                &RunContext::default(),
                &policy()
            )
            .is_err()
    );
    let mut input = system(a("1/(eps*(eps+1+r))"), a("eps"));
    input.roots[0].radicand = a("1+s");
    let original = flow(input)?;
    let gauge = original.regularize_epsilon(&RunContext::default())?;
    assert!(
        gauge
            .to_sheared_boundary(&seed(&original, a("0"), RootSheet::Principal, 1)?, 1)
            .is_err()
    );
    Ok(())
}

fn check_analytic(boundary: &CachedBoundary, delta: &ComplexFloat) -> Result<()> {
    let p = Precision::decimal(100)?;
    let ep = p.exp(delta);
    let em = p.exp(&p.neg(delta));
    let sinh = p.div(&p.sub(&ep, &em), &p.i(2));
    let cosh = p.div(&p.add(&ep, &em), &p.i(2));
    for (row, expected) in boundary
        .coefficients
        .iter()
        .zip([vec![sinh.clone(), p.zero()], vec![sinh, cosh]])
    {
        for (found, expected) in row.iter().zip(expected) {
            assert!(p.close(found, &expected, 25), "{found} != {expected}");
        }
    }
    Ok(())
}

#[test]
fn complex_root_cache_maps_extra_orders_errors_sheets_and_binary_restart() -> Result<()> {
    let original = flow(system(a("1/(2*eps*r)"), a("eps/(2*r)")))?;
    let gauge = original.regularize_epsilon(&RunContext::default())?;
    let original_range = EpsilonRange::new(-1, 0)?;
    let range = gauge.required_range(original_range)?;
    assert_eq!(range, EpsilonRange::new(-1, 1)?);
    assert!(
        gauge
            .to_sheared_boundary(&seed(&original, a("1"), RootSheet::Principal, 0)?, 1)
            .is_err()
    );
    let p = Precision::decimal(100)?;
    let i = Atom::num(Complex::new(Rational::from(0), Rational::from(1)));
    let destination = a("2") + &i;
    for sheet in [RootSheet::Principal, RootSheet::Opposite] {
        let input = seed(&original, a("1"), sheet, 1)?;
        let shifted = gauge.to_sheared_boundary(&input, 1)?;
        assert_eq!(
            shifted.point.restart_coordinates()?,
            input.point.restart_coordinates()?
        );
        assert_eq!(shifted.point.root_germ(), input.point.root_germ());
        assert_eq!(
            shifted.accuracy.verified_digits(),
            input.accuracy.verified_digits()
        );
        assert_eq!(
            shifted.accuracy.comparison_errors()[1][0],
            input.accuracy.comparison_errors()[0][0]
        );
        let mut bank = RustFlowCache::default();
        bank.insert(shifted)?;
        let first = gauge.flow().evaluate_to(
            &mut bank,
            &point(a("1") + &i),
            &germ(sheet),
            range,
            &options(),
            &RunContext::default(),
            &policy(),
        )?;
        assert!(first.transport.is_some());
        let result = gauge.flow().evaluate_to(
            &mut bank,
            &point(destination.clone()),
            &germ(sheet),
            range,
            &options(),
            &RunContext::default(),
            &policy(),
        )?;
        let restored = gauge.to_original_boundary(&result.boundary, original_range)?;
        let root = p.pow(
            &p.eval(&destination, &Default::default())?,
            &p.rational(&Rational::from((1, 2))),
        );
        let delta = p.sub(&root, &p.i(1));
        let delta = if sheet == RootSheet::Principal {
            delta
        } else {
            p.neg(&delta)
        };
        check_analytic(&restored, &delta)?;
        assert_eq!(restored.identity.key(), original.identity().key());
        assert_ne!(restored.identity.key(), result.boundary.identity.key());
        let directory = std::env::temp_dir().join(format!(
            "epsilon-algebraic-{}-{sheet:?}",
            std::process::id()
        ));
        bank.save(&directory)?;
        let mut reloaded = RustFlowCache::load(&directory)?;
        std::fs::remove_dir_all(directory)?;
        let hit = gauge.flow().evaluate_to(
            &mut reloaded,
            &point(destination.clone()),
            &germ(sheet),
            range,
            &options(),
            &RunContext::default(),
            &policy(),
        )?;
        assert!(hit.transport.is_none());
        assert_eq!(hit.boundary.coefficients, result.boundary.coefficients);
        let other = if sheet == RootSheet::Principal {
            RootSheet::Opposite
        } else {
            RootSheet::Principal
        };
        let before = reloaded.len();
        assert!(
            gauge
                .flow()
                .evaluate_to(
                    &mut reloaded,
                    &point(destination.clone()),
                    &germ(other),
                    range,
                    &options(),
                    &RunContext::default(),
                    &policy()
                )
                .is_err()
        );
        assert_eq!(reloaded.len(), before);
    }
    Ok(())
}

#[test]
fn two_coordinate_root_connection_shares_one_gauge_and_rejects_conflicting_partials() -> Result<()>
{
    let mut input = system(a("1/(2*eps*r)"), a("eps/(2*r)"));
    input.roots.push(SquareRoot {
        symbol: sym("q"),
        radicand: a("t"),
    });
    input.derivatives.insert(
        sym("t"),
        vec![
            vec![Atom::new(), a("1/(2*eps*q)")],
            vec![a("eps/(2*q)"), Atom::new()],
        ],
    );
    let original = flow(input.clone())?;
    let regular = original.regularize_epsilon(&RunContext::default())?;
    assert_eq!(regular.shearing().weights(), &[-1, 0]);
    let sheets = RootGerm {
        sheets: BTreeMap::from([
            (sym("r"), RootSheet::Principal),
            (sym("q"), RootSheet::Principal),
        ]),
    };
    let mut source = seed(&original, a("1"), RootSheet::Principal, 1)?;
    source.point = CachedPoint::Exact(BTreeMap::from([(sym("s"), a("1")), (sym("t"), a("1"))]))
        .with_root_germ(sheets.clone())?;
    let mut bank = RustFlowCache::default();
    bank.insert(regular.to_sheared_boundary(&source, 1)?)?;
    let result = regular.flow().evaluate_to(
        &mut bank,
        &BTreeMap::from([(sym("s"), a("4")), (sym("t"), a("9"))]),
        &sheets,
        EpsilonRange::new(-1, 1)?,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    let restored = regular.to_original_boundary(&result.boundary, EpsilonRange::new(-1, 0)?)?;
    check_analytic(&restored, &Precision::decimal(100)?.i(3))?;
    input.derivatives.insert(
        sym("s"),
        vec![
            vec![Atom::new(), a("1/(eps*r)")],
            vec![Atom::new(), Atom::new()],
        ],
    );
    input.derivatives.insert(
        sym("t"),
        vec![
            vec![Atom::new(), Atom::new()],
            vec![a("1/(eps*q)"), Atom::new()],
        ],
    );
    for variable in [sym("s"), sym("t")] {
        let mut single = input.clone();
        // Retain the other derivative as zero so both physical root parameters
        // remain declared, without contributing a second shearing constraint.
        for (&v, matrix) in &mut single.derivatives {
            if v != variable {
                *matrix = vec![vec![Atom::new(); 2]; 2];
            }
        }
        EpsilonShearing::regularize_algebraic(&single, &RunContext::default())?;
    }
    assert!(matches!(
        EpsilonShearing::regularize_algebraic(&input, &RunContext::default()),
        Err(Error::Unsupported(_))
    ));
    Ok(())
}

#[test]
fn prescribed_root_transport_keeps_continuation_identity_and_restores_original_values() -> Result<()>
{
    for side in [Prescription::PlusI0, Prescription::MinusI0] {
        let original = flow(system(a("1/(2*eps*r)"), a("eps/(2*r)")))?
            .with_prescribed_continuation(PhysicalContinuation {
                prescriptions: vec![contour::PolynomialPrescription {
                    polynomial: a("s"),
                    prescription: side,
                }],
                unprescribed_side: side,
                domain: "one specified detour with no additional winding".into(),
            })?;
        let gauge = original.regularize_epsilon(&RunContext::default())?;
        let range = EpsilonRange::new(-1, 1)?;
        let mut bank = RustFlowCache::default();
        bank.insert(
            gauge.to_sheared_boundary(&seed(&original, a("1"), RootSheet::Principal, 1)?, 1)?,
        )?;
        let sheet = if side == Prescription::PlusI0 {
            RootSheet::Principal
        } else {
            RootSheet::Opposite
        };
        assert!(
            gauge
                .flow()
                .evaluate_prescribed_to(
                    &mut bank,
                    &point(a("-1")),
                    &germ(sheet),
                    range,
                    &options(),
                    &RunContext::default(),
                    &policy(),
                    &|_: &CachedBoundary,
                      _: &CachedPoint,
                      _: &physical_transport::PhysicalRoute| Ok(false)
                )
                .is_err()
        );
        let result = gauge.flow().evaluate_prescribed_to(
            &mut bank,
            &point(a("-1")),
            &germ(sheet),
            range,
            &options(),
            &RunContext::default(),
            &policy(),
            &|_: &CachedBoundary, _: &CachedPoint, _: &physical_transport::PhysicalRoute| Ok(true),
        )?;
        let restored = gauge.to_original_boundary(&result.boundary, EpsilonRange::new(-1, 0)?)?;
        let p = Precision::decimal(100)?;
        let delta = p.parse(
            "-1",
            if side == Prescription::PlusI0 {
                "1"
            } else {
                "-1"
            },
        )?;
        check_analytic(&restored, &delta)?;
        assert_eq!(restored.identity.key(), original.identity().key());
        let directory = std::env::temp_dir().join(format!(
            "epsilon-prescribed-root-{}-{side:?}",
            std::process::id()
        ));
        bank.save(&directory)?;
        let mut loaded = RustFlowCache::load(&directory)?;
        std::fs::remove_dir_all(directory)?;
        let hit = gauge.flow().evaluate_prescribed_to(
            &mut loaded,
            &point(a("-1")),
            &germ(sheet),
            range,
            &options(),
            &RunContext::default(),
            &policy(),
            &|_: &CachedBoundary, _: &CachedPoint, _: &physical_transport::PhysicalRoute| Ok(true),
        )?;
        assert!(hit.transport.is_none());
    }
    Ok(())
}
