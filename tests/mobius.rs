use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{AlgebraicSystem, RootSeed, SquareRoot};
use symbolica_amflow::diffexp::{EpsilonBoundary, EpsilonSystem};
use symbolica_amflow::local_coordinates::{LocalCoordinate, TaylorCoordinate};
use symbolica_amflow::*;

fn options(working: u32, order: usize, mapped: bool) -> FlowOptions {
    FlowOptions {
        digits: 20,
        guard_digits: working - 20,
        series_order: order,
        local_coordinate: if mapped {
            LocalCoordinate::BalancedMobius
        } else {
            LocalCoordinate::Identity
        },
        max_steps: 10000,
        ..Default::default()
    }
}
fn scalar() -> DifferentialSystem {
    let x = symbol!("mobius_test::x");
    DifferentialSystem {
        variable: x,
        matrix: vec![vec![
            Atom::num(1) / (Atom::var(x) + Atom::num((1, 10)))
                - Atom::num(1) / (Atom::num(1) - Atom::var(x)),
        ]],
    }
}

#[test]
fn asymmetric_rational_map_extends_chart_and_refines_in_physical_coordinates() -> Result<()> {
    let system = EpsilonSystem {
        variable: scalar().variable,
        matrices: vec![scalar().matrix],
    };
    let mut reference = None;
    for (working, order) in [(50, 96), (70, 128)] {
        let p = Precision::decimal(working)?;
        let boundary = EpsilonBoundary {
            point: p.zero(),
            leading: 0,
            coefficients: vec![vec![p.i(1)]],
        };
        let compiled = system.compile(p, &Default::default())?;
        let target = p.rational(&Rational::from((1, 4)));
        let ordinary = compiled.transport(
            &boundary,
            std::slice::from_ref(&target),
            &options(working, order, false),
            &RunContext::default(),
            true,
        )?;
        let mapped = compiled.transport(
            &boundary,
            std::slice::from_ref(&target),
            &options(working, order, true),
            &RunContext::default(),
            true,
        )?;
        let expected = p.rational(&Rational::from((21, 8)));
        assert!(p.close(&mapped.coefficients[0][0], &expected, 35));
        assert!(p.close(&mapped.coefficients[0][0], &ordinary.coefficients[0][0], 35));
        assert_eq!(mapped.point, target);
        assert!(matches!(
            mapped.segments[0].coordinate,
            TaylorCoordinate::Mobius(_)
        ));
        assert!(p.norm(&mapped.segments[0].end) > p.norm(&ordinary.segments[0].end));
        for (i, segment) in mapped.segments.iter().enumerate() {
            let mid = p.scale(&p.add(&segment.center, &segment.end), 1, 2);
            let exact = p.mul(
                &p.add(&p.i(1), &p.scale(&mid, 10, 1)),
                &p.sub(&p.i(1), &mid),
            );
            assert!(p.close(&mapped.evaluate_segment(i, &mid)?[0][0], &exact, 35));
            assert!(
                mapped
                    .evaluate_segment(i, &p.add(&segment.end, &p.i(1)))
                    .is_err()
            );
        }
        if let Some(old) = reference {
            assert!(p.close(&mapped.coefficients[0][0], &old, 35));
        }
        reference = Some(mapped.coefficients[0][0].clone());
    }
    Ok(())
}

#[test]
fn coupled_epsilon_hierarchy_uses_the_same_mapped_recurrence() -> Result<()> {
    let x = scalar().variable;
    let f = Atom::num(1) / (Atom::var(x) + Atom::num((1, 10)));
    let g = Atom::num(1) / (Atom::num(1) - Atom::var(x));
    let zero = vec![vec![Atom::zero(); 2]; 2];
    let system = EpsilonSystem {
        variable: x,
        matrices: vec![
            zero.clone(),
            vec![vec![Atom::zero(), f], vec![g, Atom::zero()]],
            zero.clone(),
            zero.clone(),
            zero,
        ],
    };
    let mut baseline = None;
    for (working, order, mapped) in [(50, 56, false), (50, 56, true), (70, 80, true)] {
        let p = Precision::decimal(working)?;
        let mut coefficients = vec![vec![p.zero(); 2]; 5];
        coefficients[0] = vec![p.i(1), p.i(2)];
        let result = system.compile(p, &Default::default())?.transport(
            &EpsilonBoundary {
                point: p.zero(),
                leading: -2,
                coefficients,
            },
            &[p.rational(&Rational::from((1, 4)))],
            &options(working, order, mapped),
            &RunContext::default(),
            true,
        )?;
        if let Some(ref old) = baseline {
            for (a, b) in result.coefficients.iter().flatten().zip(old) {
                assert!(p.close(a, b, 35));
            }
        }
        baseline = Some(
            result
                .coefficients
                .into_iter()
                .flatten()
                .collect::<Vec<_>>(),
        );
    }
    Ok(())
}

#[test]
fn complex_legs_and_root_windings_keep_the_accepted_physical_sheet() -> Result<()> {
    let x = symbol!("mobius_winding::x");
    let r = symbol!("mobius_winding::r");
    let system = AlgebraicSystem::ordinary(
        DifferentialSystem {
            variable: x,
            matrix: vec![vec![Atom::one() / Atom::var(r)]],
        },
        vec![SquareRoot {
            symbol: r,
            radicand: Atom::var(x),
        }],
    );
    for (working, order) in [(45, 56), (60, 80)] {
        let p = Precision::decimal(working)?;
        let compiled = system.compile(p)?;
        let route = vec![
            p.complex(1, 1),
            p.complex(-1, 1),
            p.complex(-1, -1),
            p.complex(1, -1),
            p.i(1),
        ];
        let result = compiled.transport_ordinary(
            &BoundaryData {
                point: p.i(1),
                values: vec![p.i(1)],
            },
            &route,
            &BTreeMap::from([(r, RootSeed::Principal)]),
            &options(working, order, true),
            &RunContext::default(),
            true,
        )?;
        assert_eq!(result.solution.point, p.i(1));
        assert!(p.close(&result.solution.coefficients[0][0], &p.exp(&p.i(-4)), 30));
        assert!(p.close(&result.branches.roots[0].1, &p.i(-1), 30));
        assert_eq!(result.branch_segments.len(), result.solution.segments.len());
        assert!(
            result
                .solution
                .segments
                .iter()
                .any(|s| matches!(s.coordinate, TaylorCoordinate::Mobius(_)))
        );
        for (i, segment) in result.solution.segments.iter().enumerate() {
            let value = result.solution.evaluate_segment(i, &segment.end)?;
            assert!(value[0][0].re.is_finite());
        }
        let twice = route
            .iter()
            .cloned()
            .chain(route.iter().cloned())
            .collect::<Vec<_>>();
        let result = compiled.transport_ordinary(
            &BoundaryData {
                point: p.i(1),
                values: vec![p.i(1)],
            },
            &twice,
            &BTreeMap::from([(r, RootSeed::Principal)]),
            &options(working, order, true),
            &RunContext::default(),
            false,
        )?;
        assert!(p.close(&result.solution.coefficients[0][0], &p.i(1), 30));
        assert!(p.close(&result.branches.roots[0].1, &p.i(1), 30));
    }
    Ok(())
}

#[test]
fn original_source_holes_and_budget_cancellation_survive_mapping() -> Result<()> {
    let x = symbol!("mobius_hole::x");
    let a = Atom::var(x);
    // The original negative power is retained even though exact rational
    // simplification of the first summand erases the hole at x=1.
    let hole = (Atom::one() + (a.clone().pow(2) - Atom::num(2))) / (a.clone() - Atom::one());
    let system = DifferentialSystem {
        variable: x,
        matrix: vec![vec![hole + Atom::one() / (a + Atom::one())]],
    };
    let p = Precision::decimal(50)?;
    let compiled = system.compile(p, &Default::default())?;
    assert!(
        compiled
            .transport(
                &BoundaryData {
                    point: p.i(1),
                    values: vec![p.i(1)]
                },
                &[p.i(2)],
                &options(50, 56, true),
                &RunContext::default()
            )
            .is_err()
    );
    let compiled = scalar().compile(p, &Default::default())?;
    let mut limited = options(50, 8, true);
    limited.max_steps = 1;
    assert!(matches!(
        compiled.transport(
            &BoundaryData {
                point: p.zero(),
                values: vec![p.i(1)]
            },
            &[p.rational(&Rational::from((1, 4)))],
            &limited,
            &RunContext::default()
        ),
        Err(Error::Limit(_))
    ));
    let context = RunContext::default();
    context.cancellation.cancel();
    assert!(matches!(
        compiled.transport(
            &BoundaryData {
                point: p.zero(),
                values: vec![p.i(1)]
            },
            &[p.i(1)],
            &options(50, 56, true),
            &context
        ),
        Err(Error::Cancelled)
    ));
    Ok(())
}

#[test]
fn entire_systems_keep_identity_fallback_and_sparse_defect_checks() -> Result<()> {
    let x = symbol!("mobius_entire::x");
    let p = Precision::decimal(40)?;
    let system = DifferentialSystem {
        variable: x,
        matrix: vec![vec![Atom::var(x).pow(20)]],
    }
    .compile(p, &Default::default())?;
    let boundary = BoundaryData {
        point: p.zero(),
        values: vec![p.i(1)],
    };
    let ordinary = system.transport(
        &boundary,
        &[p.i(1)],
        &options(40, 16, false),
        &RunContext::default(),
    )?;
    let mapped = system.transport(
        &boundary,
        &[p.i(1)],
        &options(40, 16, true),
        &RunContext::default(),
    )?;
    assert_eq!(ordinary.diagnostics.steps, mapped.diagnostics.steps);
    assert_eq!(
        ordinary.diagnostics.rejected_steps,
        mapped.diagnostics.rejected_steps
    );
    assert!(mapped.diagnostics.rejected_steps > 0);
    assert!(mapped.diagnostics.steps < 10000);
    assert_eq!(mapped.point, p.i(1));
    assert!(p.close(
        &mapped.values[0],
        &p.exp(&p.rational(&Rational::from((1, 21)))),
        25
    ));
    Ok(())
}

#[test]
fn reverse_and_complex_legs_keep_the_same_rational_solution() -> Result<()> {
    let system = EpsilonSystem {
        variable: scalar().variable,
        matrices: vec![scalar().matrix],
    };
    let p = Precision::decimal(50)?;
    let compiled = system.compile(p, &Default::default())?;
    let quarter = p.rational(&Rational::from((1, 4)));
    let reverse = compiled.transport(
        &EpsilonBoundary {
            point: quarter.clone(),
            leading: 0,
            coefficients: vec![vec![p.rational(&Rational::from((21, 8)))]],
        },
        &[p.zero()],
        &options(50, 96, true),
        &RunContext::default(),
        true,
    )?;
    assert_eq!(reverse.point, p.zero());
    assert!(p.close(&reverse.coefficients[0][0], &p.i(1), 35));
    let complex = p.add(&quarter, &p.mul(&quarter, &p.complex(0, 1)));
    let result = compiled.transport(
        &EpsilonBoundary {
            point: p.zero(),
            leading: 0,
            coefficients: vec![vec![p.i(1)]],
        },
        &[complex, quarter.clone()],
        &options(50, 96, true),
        &RunContext::default(),
        true,
    )?;
    assert_eq!(result.point, quarter);
    assert!(p.close(
        &result.coefficients[0][0],
        &p.rational(&Rational::from((21, 8))),
        35
    ));
    // The coordinate metadata is inseparable from its actual physical center.
    let mut malformed = result.clone();
    let bad = p.add(&malformed.segments[0].center, &p.i(1));
    malformed.segments[0].center = bad.clone();
    malformed.segments[0].end = p.add(&bad, &p.i(1));
    assert!(matches!(
        malformed.evaluate_segment(0, &bad),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}

#[test]
fn bounded_bracket_controller_uses_the_same_mapped_physical_checks() -> Result<()> {
    let p = Precision::decimal(40)?;
    let compiled = scalar().compile(p, &Default::default())?;
    let mut previous = None;
    for strategy in [StepSizeStrategy::Halving, StepSizeStrategy::Bracketed] {
        let mut opts = options(40, 24, true);
        opts.step_size_strategy = strategy;
        let result = compiled.transport(
            &BoundaryData {
                point: p.zero(),
                values: vec![p.i(1)],
            },
            &[p.rational(&Rational::from((1, 4)))],
            &opts,
            &RunContext::default(),
        )?;
        assert!(p.close(&result.values[0], &p.rational(&Rational::from((21, 8))), 25));
        let d = &result.diagnostics;
        assert_eq!(
            d.predicate_evaluations,
            d.steps + d.rejected_steps + d.superseded_successes
        );
        assert!(d.rejected_steps > 0);
        if let Some(old) = previous {
            assert!(p.close(&old, &result.values[0], 25));
        }
        previous = Some(result.values[0].clone());
    }
    Ok(())
}

#[test]
fn mapped_full_degree_residual_rejects_a_sparse_collocation_alias() -> Result<()> {
    let x = symbol!("mobius_hidden_degree::x");
    let a = Atom::var(x);
    // The geometric proposal ends at 2/3. The exact forcing vanishes at that
    // point and its physical midpoint 1/3, while every retained order is zero.
    // Its full-degree differential defect is not zero on the intervening path.
    let g = Atom::num((3, 2)).pow(1000)
        * Atom::num(1002001)
        * a.clone().pow(1000)
        * (a.clone() - Atom::num((1, 3)))
        * (a.clone() - Atom::num((2, 3)))
        / (a.clone().pow(2) - Atom::num(4));
    let system = DifferentialSystem {
        variable: x,
        matrix: vec![
            vec![Atom::zero(), g.clone()],
            vec![Atom::zero(), Atom::zero()],
        ],
    };
    let r = symbol!("mobius_hidden_degree::r");
    let algebraic = AlgebraicSystem::ordinary(
        DifferentialSystem {
            variable: x,
            matrix: vec![
                vec![Atom::zero(), g * Atom::var(r)],
                vec![Atom::zero(), Atom::zero()],
            ],
        },
        vec![SquareRoot {
            symbol: r,
            radicand: Atom::one(),
        }],
    );
    let p = Precision::decimal(50)?;
    let compiled = system.compile(p, &Default::default())?;
    for strategy in [StepSizeStrategy::Halving, StepSizeStrategy::Bracketed] {
        let mut opts = options(50, 16, true);
        opts.max_steps = 1;
        opts.step_size_strategy = strategy;
        let answer = compiled.transport(
            &BoundaryData {
                point: p.zero(),
                values: vec![p.zero(), p.i(1)],
            },
            &[p.rational(&Rational::from((2, 3)))],
            &opts,
            &RunContext::default(),
        );
        assert!(matches!(answer, Err(Error::Limit(_))), "{answer:?}");
        let answer = algebraic.compile(p)?.transport_ordinary(
            &BoundaryData {
                point: p.zero(),
                values: vec![p.zero(), p.i(1)],
            },
            &[p.rational(&Rational::from((2, 3)))],
            &BTreeMap::from([(r, RootSeed::Principal)]),
            &opts,
            &RunContext::default(),
            false,
        );
        assert!(matches!(answer, Err(Error::Limit(_))), "{answer:?}");
    }
    Ok(())
}
