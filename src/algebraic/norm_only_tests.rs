//! Small public-source comparisons; no trajectory, boundary or private fixtures.
use super::*;

// Native pole enumeration does not promise a vector order. Preserve exact
// values and multiplicities while comparing independent compilations.
fn assert_same_multiset<T: std::fmt::Debug + PartialEq>(a: &[T], b: &[T]) {
    assert_eq!(a.len(), b.len());
    let mut remaining = b.iter().collect::<Vec<_>>();
    for value in a {
        let index = remaining
            .iter()
            .position(|entry| *entry == value)
            .unwrap_or_else(|| panic!("missing exact value {value:?}"));
        remaining.swap_remove(index);
    }
    assert!(remaining.is_empty());
}

fn assert_same_error<T, U>(full: Result<T>, norm: Result<U>) {
    let full = match full {
        Err(error) => error,
        Ok(_) => panic!("full path unexpectedly accepted"),
    };
    let norm = match norm {
        Err(error) => error,
        Ok(_) => panic!("norm path unexpectedly accepted"),
    };
    assert_eq!(std::mem::discriminant(&full), std::mem::discriminant(&norm));
    assert_eq!(full.to_string(), norm.to_string());
}

#[test]
fn norm_only_matches_full_at_two_precisions() -> Result<()> {
    let x = symbol!("algebraic_norm_only_equal::x");
    let r = symbol!("algebraic_norm_only_equal::r");
    let t = Atom::var(x);
    let root = Atom::var(r);
    let i = Atom::num(Complex::new(Rational::zero(), Rational::one()));
    let source = AlgebraicSystem {
        system: EpsilonSystem {
            variable: x,
            matrices: vec![
                vec![
                    vec![Atom::zero(), (Atom::one() + &i * &t) * &root / (&t + 2)],
                    vec![Atom::zero(), Atom::zero()],
                ],
                vec![
                    vec![Atom::one() / (&t + &root + 3), Atom::zero()],
                    vec![&root / (&t + 4).pow(2), Atom::zero()],
                ],
            ],
        },
        roots: vec![SquareRoot {
            symbol: r,
            radicand: (&t + 5) / (&t + 6),
        }],
        nonzero_conditions: vec![(&t - 7) / (&t - 8)],
    };
    let context = RunContext::default();
    let prepared = source.prepare_with_context(&context)?;
    let norm_prepared = source.prepare_norm_with_context(&context)?;
    assert!(prepared.residual_entries.is_some());
    assert!(norm_prepared.inner.residual_entries.is_none());
    assert!(norm_prepared.inner.residual_terms.is_empty());
    // Even an accidental internal call cannot expose incomplete transport.
    assert!(matches!(
        norm_prepared.inner.compile(Precision { bits: 80 }),
        Err(Error::Unsupported(_))
    ));
    for bits in [80, 201] {
        let p = Precision { bits };
        let full = prepared.compile_with_context(p, &context)?;
        let norm = norm_prepared.compile_with_context(p, &context)?;
        let kernel = &norm.inner;
        assert!(kernel.residual_rows.is_empty());
        assert!(!full.residual_rows.is_empty());
        assert_eq!(kernel.p.bits, full.p.bits);
        assert_eq!(
            (kernel.variable, kernel.size, kernel.count),
            (full.variable, full.size, full.count)
        );
        assert_eq!(kernel.source.system.matrices, full.source.system.matrices);
        assert_eq!(
            kernel.source.nonzero_conditions,
            full.source.nonzero_conditions
        );
        assert_same_multiset(&kernel.poles, &full.poles);
        assert_same_multiset(&kernel.pole_polynomials, &full.pole_polynomials);
        assert_eq!(kernel.domain_guards, full.domain_guards);
        assert_eq!(kernel.terms.len(), full.terms.len());
        for (a, b) in kernel.terms.iter().zip(&full.terms) {
            assert_eq!(
                (a.shift, a.row, a.column, &a.roots),
                (b.shift, b.row, b.column, &b.roots)
            );
            assert_eq!(a.coefficient.numerator, b.coefficient.numerator);
            assert_eq!(a.coefficient.denominator, b.coefficient.denominator);
        }
        assert_eq!(kernel.roots.len(), full.roots.len());
        for (a, b) in kernel.roots.iter().zip(&full.roots) {
            assert_eq!(a.definition.symbol, b.definition.symbol);
            assert_eq!(a.definition.radicand, b.definition.radicand);
            assert_eq!(a.value.numerator, b.value.numerator);
            assert_eq!(a.value.denominator, b.value.denominator);
            assert_eq!(
                a.logarithmic_derivative.numerator,
                b.logarithmic_derivative.numerator
            );
            assert_eq!(
                a.logarithmic_derivative.denominator,
                b.logarithmic_derivative.denominator
            );
        }
        let a = kernel.root_system.as_ref().unwrap();
        let b = full.root_system.as_ref().unwrap();
        assert_same_multiset(&a.poles, &b.poles);
        assert_same_multiset(&a.pole_polynomials, &b.pole_polynomials);
        assert_eq!(a.polynomial_rows.len(), b.polynomial_rows.len());
        for (a, b) in a.polynomial_rows.iter().zip(&b.polynomial_rows) {
            assert_eq!(a.denominator, b.denominator);
            assert_eq!(a.entries, b.entries);
        }
        let weights = vec![p.real(1), p.real(2), p.real(3), p.real(5)];
        let start = p.parse("0.125", "0.03125")?;
        let end = p.parse("0.25", "0.0625")?;
        assert_eq!(
            full.error_norm_integral_weighted(&start, &end, &weights)?,
            norm.error_norm_integral_weighted(&start, &end, &weights)?
        );
        // Named branches and their signs are retained internally, although the
        // norm wrapper intentionally exposes no branch/transport operations.
        for seed in [RootSeed::Principal, RootSeed::Opposite] {
            let seeds = BTreeMap::from([(r, seed)]);
            let a = kernel.branch_state_at(&start, &seeds)?;
            let b = full.branch_state_at(&start, &seeds)?;
            assert_eq!(a.point, b.point);
            assert_eq!(a.roots, b.roots);
        }
        let seeds = BTreeMap::from([(r, RootSeed::Value(p.zero()))]);
        assert_same_error(
            full.branch_state_at(&start, &seeds),
            kernel.branch_state_at(&start, &seeds),
        );
        assert_same_error(
            full.error_norm_integral_weighted(&start, &end, &[]),
            norm.error_norm_integral_weighted(&start, &end, &[]),
        );
    }
    Ok(())
}

#[test]
fn norm_only_keeps_cancelled_domain_and_root_rejections() -> Result<()> {
    let x = symbol!("algebraic_norm_only_domain::x");
    let r = symbol!("algebraic_norm_only_domain::r");
    let t = Atom::var(x);
    let root = Atom::var(r);
    let source = AlgebraicSystem {
        system: EpsilonSystem {
            variable: x,
            // Zero kernel, with the uncancelled source inverse norm x-1.
            matrices: vec![vec![vec![(root.pow(2) - &t) / (Atom::var(r) - 1)]]],
        },
        roots: vec![SquareRoot {
            symbol: r,
            radicand: t.clone(),
        }],
        nonzero_conditions: vec![(&t - 3) / (&t - 4)],
    };
    let context = RunContext::default();
    let p = Precision { bits: 120 };
    let full = source.prepare()?.compile(p)?;
    let norm = source
        .prepare_norm_with_context(&context)?
        .compile_with_context(p, &context)?;
    assert!(full.terms.is_empty());
    assert!(norm.inner.terms.is_empty());
    assert_same_multiset(&full.poles, &norm.inner.poles);
    assert_same_multiset(&full.pole_polynomials, &norm.inner.pole_polynomials);
    assert_eq!(full.domain_guards, norm.inner.domain_guards);
    let seeds = BTreeMap::from([(r, RootSeed::Principal)]);
    for pole in [0, 1, 3, 4] {
        assert!(full.singularities().contains(&p.i(pole)));
        assert_same_error(
            full.branch_state_at(&p.i(pole), &seeds),
            norm.inner.branch_state_at(&p.i(pole), &seeds),
        );
    }
    for excluded in [1, 3, 4] {
        assert_same_error(
            full.check_domain(&p.i(excluded)),
            norm.inner.check_domain(&p.i(excluded)),
        );
    }
    full.check_domain(&p.i(2))?;
    norm.inner.check_domain(&p.i(2))?;
    // A root branch zero must still fail the shared disk norm, even for a
    // cancelled kernel. Source-domain exclusions are retained separately above.
    assert_same_error(
        full.error_norm_integral_weighted(&p.zero(), &p.parse("0.125", "0")?, &[p.real(1)]),
        norm.error_norm_integral_weighted(&p.zero(), &p.parse("0.125", "0")?, &[p.real(1)]),
    );
    let invalid_domain = AlgebraicSystem {
        nonzero_conditions: vec![Atom::zero()],
        ..source.clone()
    };
    assert_same_error(
        invalid_domain.prepare(),
        invalid_domain.prepare_norm_with_context(&context),
    );
    let invalid_root = AlgebraicSystem {
        roots: vec![SquareRoot {
            symbol: r,
            radicand: Atom::zero(),
        }],
        ..source
    };
    assert_same_error(
        invalid_root.prepare(),
        invalid_root.prepare_norm_with_context(&context),
    );
    Ok(())
}

#[test]
fn norm_only_keeps_unresolved_pole_precision_and_cancellation_rejections() -> Result<()> {
    let x = symbol!("algebraic_norm_only_poles::x");
    let r = symbol!("algebraic_norm_only_poles::r");
    let t = Atom::var(x);
    let first = &t - 1;
    let second = &first - Atom::num(Rational::from((1, 1_i64 << 60)));
    let source = AlgebraicSystem::ordinary(
        DifferentialSystem {
            variable: x,
            matrix: vec![
                vec![Atom::one() / &first, Atom::one() / &second],
                vec![Atom::zero(), Atom::zero()],
            ],
        },
        vec![SquareRoot {
            symbol: r,
            radicand: &first * &second,
        }],
    );
    let context = RunContext::default();
    let full = source.prepare_with_context(&context)?;
    let norm = source.prepare_norm_with_context(&context)?;
    let p = Precision { bits: 80 };
    assert!(matches!(
        full.compile_with_context(p, &context),
        Err(Error::Numerical(message)) if message.contains("ambiguous nearby poles")
    ));
    assert_same_error(
        full.compile_with_context(p, &context),
        norm.compile_with_context(p, &context),
    );
    assert_same_error(
        full.compile_with_context(Precision { bits: 1 }, &context),
        norm.compile_with_context(Precision { bits: 1 }, &context),
    );
    context.cancellation.cancel();
    assert_same_error(
        full.compile_with_context(p, &context),
        norm.compile_with_context(p, &context),
    );
    assert_same_error(
        source.prepare_with_context(&context),
        source.prepare_norm_with_context(&context),
    );
    Ok(())
}
