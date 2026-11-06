use super::*;

// A real, exactly soluble transport keeps the supplied Precision on compiled
// arithmetic, diagnostics and retained segments, independently of FlowOptions.
fn constant_profile(
    p: Precision,
    options: &FlowOptions,
    value: C,
    leading: i32,
    save: bool,
) -> Result<EpsilonSolution> {
    EpsilonSystem {
        variable: symbol!("precision_recovery::x"),
        matrices: vec![vec![vec![Atom::zero()]]],
    }
    .compile(p, &Default::default())?
    .transport(
        &EpsilonBoundary {
            point: p.zero(),
            leading,
            coefficients: vec![vec![value]],
        },
        &[p.i(1)],
        options,
        &RunContext::default(),
        save,
    )
}

fn options() -> FlowOptions {
    FlowOptions {
        digits: 12,
        guard_digits: 48,
        series_order: 8,
        max_precision_attempts: 2,
        ..Default::default()
    }
}

#[test]
fn hinted_retry_keeps_local_budget_and_then_refines_after_success() -> Result<()> {
    let mut profiles = Vec::new();
    let minimum_bits = 451;
    let result = refine_epsilon_transport(
        &options(),
        &RunContext::default(),
        true,
        |p, local| {
            profiles.push((p.bits, local.digits, local.guard_digits, local.series_order));
            if p.bits < minimum_bits {
                return Err(Error::InsufficientPrecision {
                    minimum_bits,
                    context: "fixed local budget needs more arithmetic bits".into(),
                });
            }
            let delta = if profiles.len() == 2 { "1e-14" } else { "2e-14" };
            let value = p.add(&p.i(1), &p.parse(delta, "0")?);
            let result = constant_profile(p, local, value, 0, true)?;
            assert_eq!(result.diagnostics.working_bits, p.bits);
            assert!(result.segments.iter().all(|s| s.working_bits == p.bits));
            Ok(result)
        },
    )?;
    assert_eq!(profiles.len(), 3);
    assert_eq!(profiles.iter().map(|p| p.1).collect::<Vec<_>>(), [12, 12, 12]);
    assert_eq!(profiles.iter().map(|p| p.2).collect::<Vec<_>>(), [48, 48, 68]);
    assert_eq!(profiles.iter().map(|p| p.3).collect::<Vec<_>>(), [8, 40, 72]);
    assert_eq!(profiles[0].0, Precision::decimal(60)?.bits);
    assert!(profiles[1].0 >= minimum_bits && profiles[2].0 > profiles[1].0);
    // The retry adds arithmetic margin over the unchanged local profile.
    assert!(profiles[1].0 > Precision::decimal(profiles[1].1 + profiles[1].2)?.bits);
    assert!(profiles[2].0 > Precision::decimal(profiles[2].1 + profiles[2].2)?.bits);
    assert_eq!(result.verified_digits, Some(12));
    let p = Precision { bits: result.diagnostics.working_bits };
    assert!(result.comparison_errors[0][0] > p.real(0));
    assert!(result.comparison_errors[0][0] < p.tolerance(12));
    assert_eq!(result.checkpoints.len(), 1);
    assert_eq!(result.checkpoints[0].comparison_errors, result.comparison_errors);
    Ok(())
}

#[test]
fn failure_discards_previous_and_one_recovered_profile_is_not_verified() -> Result<()> {
    for max_precision_attempts in [2, 3] {
        let mut calls = 0;
        let mut guards = Vec::new();
        let opts = FlowOptions { max_precision_attempts, ..options() };
        let result = refine_epsilon_transport(&opts, &RunContext::default(), false, |p, local| {
            calls += 1;
            guards.push(local.guard_digits);
            if calls == 2 {
                return Err(Error::InsufficientPrecision {
                    minimum_bits: 451,
                    context: "intervening arithmetic failure".into(),
                });
            }
            // A stale first profile would also cause a spurious leading-power
            // mismatch. Only the two fresh leading=1 profiles may be compared.
            let leading = if calls == 1 { 0 } else { 1 };
            constant_profile(p, local, p.i(1), leading, false)
        });
        assert_eq!(calls, max_precision_attempts + 1);
        assert_eq!(guards, [48, 68, 68, 88][..calls]);
        if max_precision_attempts == 2 {
            assert!(matches!(result, Err(Error::Accuracy(_))));
        } else {
            let result = result?;
            assert_eq!(result.leading, 1);
            assert_eq!(result.verified_digits, Some(12));
            assert!(result.checkpoints.is_empty());
        }
    }
    Ok(())
}

#[test]
fn ordinary_success_still_requires_two_agreeing_independent_profiles() -> Result<()> {
    for max_precision_attempts in [1, 2] {
        let mut calls = 0;
        let opts = FlowOptions { max_precision_attempts, ..options() };
        let result = refine_epsilon_transport(&opts, &RunContext::default(), false, |p, local| {
            calls += 1;
            // With no arithmetic failures, retain the original +20/+32 profile.
            assert_eq!(p.bits, Precision::decimal(60 + 20 * (calls as u32 - 1))?.bits);
            assert_eq!(local.guard_digits, 48 + 20 * (calls as u32 - 1));
            assert_eq!(local.series_order, 8 + 32 * (calls - 1));
            constant_profile(p, local, p.i(if calls == 1 { 1 } else { 2 }), 0, false)
        });
        assert_eq!(calls, max_precision_attempts + 1);
        if max_precision_attempts == 1 {
            assert!(matches!(result, Err(Error::Accuracy(_))));
        } else {
            let result = result?;
            assert_eq!(result.verified_digits, Some(12));
            let p = Precision { bits: result.diagnostics.working_bits };
            assert_eq!(result.coefficients[0][0], p.i(2));
            assert_eq!(result.comparison_errors[0][0], p.real(0));
        }
    }
    Ok(())
}

#[test]
fn precision_exhaustion_retains_last_hint_and_never_verifies() {
    let mut calls = 0;
    let error = refine_epsilon_transport(&options(), &RunContext::default(), false, |p, local| {
        calls += 1;
        assert_eq!(local.guard_digits, 48);
        Err(Error::InsufficientPrecision {
            minimum_bits: p.bits + 100,
            context: format!("unresolved profile {calls}"),
        })
    })
    .unwrap_err();
    assert_eq!(calls, 3);
    match error {
        Error::InsufficientPrecision { context, .. } => {
            assert!(context.contains("without two successful profiles"));
            assert!(context.ends_with("unresolved profile 3"));
        }
        other => panic!("unexpected error: {other}"),
    }
}

#[test]
fn ordinary_errors_and_changed_leading_power_are_still_rejected() {
    for error in [Error::Accuracy("local rejection".into()), Error::InvalidInput("bad source".into()), Error::Cancelled] {
        let expected = error.to_string();
        let mut failure = Some(error);
        let mut calls = 0;
        let error = refine_epsilon_transport(&options(), &RunContext::default(), false, |_, _| {
            calls += 1;
            Err(failure.take().unwrap())
        })
        .unwrap_err();
        assert_eq!(error.to_string(), expected);
        assert_eq!(calls, 1);
    }
    let mut calls = 0;
    let error = refine_epsilon_transport(&options(), &RunContext::default(), false, |p, local| {
        let leading = calls;
        calls += 1;
        constant_profile(p, local, p.i(1), leading, false)
    })
    .unwrap_err();
    assert_eq!(calls, 2);
    assert!(matches!(error, Error::InvalidInput(message) if message.contains("leading epsilon power")));
}

// The pre-preparation epsilon adapter, retained as an independent compatibility
// oracle for flattening, exact denominator clearing and numerical compilation.
fn legacy_compile(
    system: &EpsilonSystem,
    p: Precision,
    values: &ahash::HashMap<Atom, C>,
) -> Result<CompiledEpsilonSystem> {
    system.validate()?;
    let rows = (0..system.matrices[0].len())
        .map(|i| {
            system
                .matrices
                .iter()
                .flat_map(|m| m[i].iter().cloned())
                .collect()
        })
        .collect::<Vec<_>>();
    Ok(CompiledEpsilonSystem {
        rows: crate::ode::compile_rows(system.variable, &rows, p, values)?,
        count: system.matrices.len(),
        strictly_positive_epsilon: system.matrices[0].iter().flatten().all(Atom::is_zero),
    })
}

fn same_compilation(a: &CompiledEpsilonSystem, b: &CompiledEpsilonSystem) {
    assert_eq!(a.count, b.count);
    assert_eq!(a.strictly_positive_epsilon, b.strictly_positive_epsilon);
    assert_eq!(a.rows.p.bits, b.rows.p.bits);
    assert_eq!(a.poles().len(), b.poles().len());
    assert!(a.poles().iter().all(|pole| b.poles().contains(pole)));
    for (a, b) in a
        .rows
        .matrix
        .iter()
        .flatten()
        .zip(b.rows.matrix.iter().flatten())
    {
        assert_eq!(a.numerator, b.numerator);
        assert_eq!(a.denominator, b.denominator);
    }
    for (a, b) in a.rows.polynomial_rows.iter().zip(&b.rows.polynomial_rows) {
        assert_eq!(a.denominator, b.denominator);
        assert_eq!(a.entries, b.entries);
    }
    for (a, b) in a
        .rows
        .exact_source_rows
        .iter()
        .zip(&b.rows.exact_source_rows)
    {
        assert_eq!(a.denominator, b.denominator);
        assert_eq!(a.denominator_factors, b.denominator_factors);
        assert_eq!(a.entries, b.entries);
    }
}

#[test]
fn prepared_epsilon_snapshot_and_fresh_parameters_match_legacy_rows() -> Result<()> {
    fn send_sync<T: Send + Sync>() {}
    send_sync::<PreparedEpsilonSystem>();
    send_sync::<crate::RustFlow>();
    send_sync::<crate::kinematics::PreparedKinematicEpsilon>();
    let x = symbol!("prepared_epsilon::x");
    let t = Atom::var(x);
    let a = Atom::var(symbol!("prepared_epsilon::a"));
    let i = Atom::num(Complex::new(Rational::zero(), Rational::one()));
    let mut system = EpsilonSystem {
        variable: x,
        matrices: vec![
            vec![
                vec![(&a + &i * &t) / (&t - 2).pow(2), Atom::one() / (&t + 3)],
                vec![Atom::new(), (&t - &a) / (&t - 2)],
            ],
            vec![
                vec![Atom::one() / (&t - 2), Atom::new()],
                vec![Atom::one() / (&t + &i), Atom::new()],
            ],
        ],
    };
    let original = system.clone();
    let context = RunContext::default();
    let prepared = system.prepare_with_context(&context)?;
    system.matrices = vec![vec![vec![Atom::num(99)]]];
    system.variable = symbol!("prepared_epsilon::changed_x");
    let storage = Precision { bits: 700 };
    for value in [
        storage.rational(&Rational::from((1, 3))),
        storage.complex(2, -3),
        storage.zero(),
    ] {
        let values = [(a.clone(), value)].into_iter().collect();
        for bits in [80, 201] {
            let p = Precision { bits };
            let actual = prepared.compile_with_context(p, &values, &context)?;
            same_compilation(&actual, &legacy_compile(&original, p, &values)?);
            same_compilation(&actual, &original.compile(p, &values)?);
            for pole in [p.i(2), p.i(-3), p.complex(0, -1)] {
                assert!(actual.poles().contains(&pole));
            }
        }
    }
    Ok(())
}

#[test]
fn prepared_epsilon_preserves_specialized_removable_poles_and_errors() -> Result<()> {
    let x = symbol!("prepared_epsilon_domains::x");
    let t = Atom::var(x);
    let a = Atom::var(symbol!("prepared_epsilon_domains::a"));
    let context = RunContext::default();
    // Specializing a=0 makes both entries identically one or zero, but their
    // original denominators still restrict the numerical domain.
    let system = EpsilonSystem {
        variable: x,
        matrices: vec![vec![vec![(&t + &a) / &t]], vec![vec![&a / (&t - 2).pow(3)]]],
    };
    let prepared = system.prepare_with_context(&context)?;
    for bits in [80, 201] {
        let p = Precision { bits };
        let values = [(a.clone(), p.zero())].into_iter().collect();
        let actual = prepared.compile(p, &values)?;
        same_compilation(&actual, &legacy_compile(&system, p, &values)?);
        assert_eq!(actual.poles().len(), 2);
        assert!(actual.poles().contains(&p.zero()));
        assert!(actual.poles().contains(&p.i(2)));
        assert_eq!(
            prepared
                .compile(p, &Default::default())
                .unwrap_err()
                .to_string(),
            legacy_compile(&system, p, &Default::default())
                .unwrap_err()
                .to_string(),
        );
        let singular = EpsilonSystem {
            variable: x,
            matrices: vec![vec![vec![Atom::one() / &a]]],
        };
        assert_eq!(
            singular
                .prepare_with_context(&context)?
                .compile(p, &values)
                .unwrap_err()
                .to_string(),
            legacy_compile(&singular, p, &values)
                .unwrap_err()
                .to_string(),
        );
    }
    Ok(())
}

#[test]
fn prepared_epsilon_validates_dimensions_precision_and_cancellation() -> Result<()> {
    let x = symbol!("prepared_epsilon_validation::x");
    let context = RunContext::default();
    for matrices in [vec![], vec![vec![]], vec![vec![vec![Atom::one()]], vec![]]] {
        let system = EpsilonSystem {
            variable: x,
            matrices,
        };
        assert_eq!(
            system
                .prepare_with_context(&context)
                .unwrap_err()
                .to_string(),
            system.validate().unwrap_err().to_string(),
        );
    }
    let system = EpsilonSystem {
        variable: x,
        matrices: vec![vec![vec![Atom::zero()]]],
    };
    let prepared = system.prepare_with_context(&context)?;
    assert!(matches!(
        prepared.compile(Precision { bits: 1 }, &Default::default()),
        Err(Error::InvalidInput(_))
    ));
    context.cancellation.cancel();
    assert!(matches!(
        system.prepare_with_context(&context),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        prepared.compile_with_context(Precision { bits: 80 }, &Default::default(), &context),
        Err(Error::Cancelled)
    ));
    Ok(())
}

#[test]
fn prepared_epsilon_complex_positive_orders_transport_at_80_and_201_bits() -> Result<()> {
    let x = symbol!("prepared_epsilon_transport::x");
    let t = Atom::var(x);
    let context = RunContext::default();
    for ordinary in [false, true] {
        let i = Atom::num(Complex::new(Rational::zero(), Rational::one()));
        let system = EpsilonSystem {
            variable: x,
            matrices: vec![
                vec![vec![if ordinary {
                    &i / (&t + 2)
                } else {
                    Atom::zero()
                }]],
                vec![vec![Atom::one() / (&t + 2)]],
                vec![vec![Atom::num(2) * &t]],
            ],
        };
        let prepared = system.prepare_with_context(&context)?;
        for (bits, digits) in [(80, 10), (201, 32)] {
            let p = Precision { bits };
            let compiled = prepared.compile(p, &Default::default())?;
            let boundary = EpsilonBoundary {
                point: p.zero(),
                leading: -2,
                coefficients: vec![
                    vec![p.complex(1, 1)],
                    vec![p.complex(2, -1)],
                    vec![p.complex(-1, 2)],
                ],
            };
            let options = FlowOptions {
                digits,
                guard_digits: 10,
                series_order: 64,
                ..Default::default()
            };
            let actual = compiled.transport(&boundary, &[p.i(1)], &options, &context, true)?;
            let legacy = legacy_compile(&system, p, &Default::default())?.transport(
                &boundary,
                &[p.i(1)],
                &options,
                &context,
                true,
            )?;
            assert_eq!(actual.coefficients, legacy.coefficients);
            assert_eq!(actual.leading, -2);
            assert_eq!(actual.diagnostics.working_bits, bits);
            assert!(!actual.segments.is_empty());
            let l = p.log(&p.rational(&Rational::from((3, 2))));
            let phase = if ordinary {
                p.exp(&p.mul(&p.complex(0, 1), &l))
            } else {
                p.i(1)
            };
            let c = &boundary.coefficients;
            let expected = [
                c[0][0].clone(),
                p.add(&c[1][0], &p.mul(&l, &c[0][0])),
                p.add(
                    &p.add(&c[2][0], &p.mul(&l, &c[1][0])),
                    &p.mul(&p.add(&p.i(1), &p.scale(&p.mul(&l, &l), 1, 2)), &c[0][0]),
                ),
            ];
            for (row, expected) in actual.coefficients.iter().zip(expected) {
                assert!(p.close(&row[0], &p.mul(&phase, &expected), digits));
            }
        }
    }
    Ok(())
}

#[test]
fn prepared_epsilon_refinement_returns_accepted_profile_and_checks_checkpoints() -> Result<()> {
    let system = EpsilonSystem {
        variable: symbol!("prepared_epsilon_refined::x"),
        matrices: vec![
            vec![vec![Atom::zero()]],
            vec![vec![Atom::one()]],
            vec![vec![Atom::zero()]],
        ],
    };
    let context = RunContext::default();
    let options = FlowOptions {
        digits: 12,
        guard_digits: 20,
        series_order: 16,
        ..Default::default()
    };
    let prepared = system.prepare_with_context(&context)?;
    // Both matching and mismatched initial precision must give the same accepted
    // trajectory; the latter must not lend its arithmetic to the first profile.
    let initial_p = Precision::decimal(options.digits + options.guard_digits)?;
    let mut reference: Option<EpsilonSolution> = None;
    for p in [initial_p, Precision { bits: 80 }] {
        let initial = prepared.compile(p, &Default::default())?;
        let (result, compiled) = prepared.transport_with_compiled(
            Some(initial),
            |p| {
                Ok(EpsilonBoundary {
                    point: p.zero(),
                    leading: 0,
                    coefficients: vec![vec![p.complex(1, 1)], vec![p.zero()], vec![p.zero()]],
                })
            },
            &[Atom::one()],
            &options,
            &context,
            true,
        )?;
        assert_eq!(result.verified_digits, Some(12));
        assert_eq!(
            result.diagnostics.working_bits,
            Precision::decimal(52)?.bits
        );
        assert_eq!(compiled.rows.p.bits, result.diagnostics.working_bits);
        assert!(!result.checkpoints.is_empty());
        let p = compiled.rows.p;
        for (row, denominator) in result.coefficients.iter().zip([1, 1, 2]) {
            assert!(p.close(&row[0], &p.scale(&p.complex(1, 1), 1, denominator), 12));
        }
        if let Some(old) = &reference {
            assert_eq!(result.coefficients, old.coefficients);
            assert_eq!(result.comparison_errors, old.comparison_errors);
        }
        reference = Some(result);
    }
    Ok(())
}

#[test]
fn prepared_rational_physical_transport_keeps_complex_source_and_finite_floor() -> Result<()> {
    use crate::physical_transport::PhysicalTransport;
    use crate::transport_cache::*;
    use std::collections::BTreeMap;
    let s = symbol!("prepared_epsilon_physical::s");
    let eps = symbol!("prepared_epsilon_physical::eps");
    let x = Atom::var(s);
    let e = Atom::var(eps);
    let i = Atom::num(Complex::new(Rational::zero(), Rational::one()));
    let engine = PhysicalTransport::new(
        crate::kinematics::KinematicSystem {
            epsilon: eps,
            derivatives: BTreeMap::from([(
                s,
                vec![vec![(&i + &e) / (&x + 2) + &e * &e * Atom::num(2) * &x]],
            )]),
        },
        &[Atom::var(symbol!("prepared_epsilon_physical::I"))],
        &Atom::one(),
        crate::Prescription::PlusI0,
        "regular interval [0,1]",
    )?;
    for (bits, verified, digits) in [(80, 16, 10), (201, 48, 32)] {
        let p = Precision { bits };
        let seed = p.complex(1, 1);
        let mut cache = RustFlowCache::default();
        cache.insert(CachedBoundary {
            identity: engine.identity().clone(),
            point: CachedPoint::Exact(BTreeMap::from([(s, Atom::zero())])),
            kind: PointKind::Physical,
            range: EpsilonRange::new(-2, 0)?,
            coefficients: vec![vec![seed], vec![p.zero()], vec![p.zero()]],
            accuracy: BoundaryAccuracy::supplied(
                verified,
                bits,
                vec![vec![p.real(0)]; 3],
                "exact complex analytic seed with finite source evidence",
            )?,
        })?;
        let result = engine.evaluate_to(
            &mut cache,
            &BTreeMap::from([(s, Atom::one())]),
            EpsilonRange::new(-2, 0)?,
            &FlowOptions {
                digits,
                guard_digits: 20,
                series_order: 64,
                ..Default::default()
            },
            &RunContext::default(),
            &ScaledDistance {
                scales: BTreeMap::new(),
                admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
            },
        )?;
        let accuracy = &result.boundary.accuracy;
        assert_eq!(accuracy.input_verified_digits(), verified);
        assert!(accuracy.verified_digits() <= verified);
        let transport = result.transport.as_ref().unwrap();
        assert!(!transport.checkpoints.is_empty());
        let p = Precision {
            bits: accuracy.working_bits(),
        };
        let l = p.log(&p.rational(&Rational::from((3, 2))));
        let base = p.mul(&p.complex(1, 1), &p.exp(&p.mul(&p.complex(0, 1), &l)));
        let expected = [
            base.clone(),
            p.mul(&base, &l),
            p.mul(&base, &p.add(&p.i(1), &p.scale(&p.mul(&l, &l), 1, 2))),
        ];
        for (row, expected) in result.boundary.coefficients.iter().zip(expected) {
            assert!(p.close(&row[0], &expected, digits));
        }
        assert!(
            accuracy
                .comparison_errors()
                .iter()
                .flatten()
                .all(|error| *error >= p.tolerance(verified))
        );
    }
    Ok(())
}
