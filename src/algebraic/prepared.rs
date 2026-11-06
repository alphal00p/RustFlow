//! Immutable exact preparation for repeated numerical registered-root compilation.
use super::*;
use crate::ode::prepared::PreparedRows;

#[derive(Clone, Debug)]
struct PreparedKernelTerm {
    shift: usize,
    row: usize,
    column: usize,
    roots: Vec<usize>,
}

/// Validated exact source and native root normalization for one fixed path.
///
/// The source snapshot is immutable. Each [`Self::compile`] freshly specializes
/// coefficients, source residuals and poles at its requested working precision;
/// no rounded values or boundary data are cached here. Source domain exclusions
/// are retained before cancellations, as in direct compilation.
#[derive(Clone, Debug)]
pub struct PreparedAlgebraicSystem {
    source: std::sync::Arc<AlgebraicSystem>,
    entries: PreparedRows,
    normalized: Vec<PreparedKernelTerm>,
    residual_entries: Option<PreparedRows>,
    residual_terms: Vec<Vec<usize>>,
    term_count: usize,
    root_matrix: Option<PreparedRows>,
    domain_guards: Vec<Atom>,
}

// Only the full public preparation may expose transport compilation. These
// wrappers never expose their inner preparations/compiled systems: omitted
// defect residuals are valid only for weighted norm evaluation.
pub(crate) struct PreparedAlgebraicNorm {
    inner: PreparedAlgebraicSystem,
}

pub(crate) struct CompiledAlgebraicNorm {
    inner: CompiledAlgebraicSystem,
}

impl PreparedAlgebraicNorm {
    pub(crate) fn compile_with_context(
        &self,
        p: Precision,
        context: &RunContext,
    ) -> Result<CompiledAlgebraicNorm> {
        Ok(CompiledAlgebraicNorm {
            inner: self
                .inner
                .compile_kernel_with_context(p, context, Vec::new())?,
        })
    }
}

impl CompiledAlgebraicNorm {
    pub(crate) fn error_norm_integral_weighted(
        &self,
        start: &C,
        end: &C,
        weights: &[Float],
    ) -> Result<Float> {
        self.inner.error_norm_integral_weighted(start, end, weights)
    }
}

impl AlgebraicSystem {
    /// Snapshot and normalize the exact system for repeated precision profiles.
    pub fn prepare(&self) -> Result<PreparedAlgebraicSystem> {
        self.prepare_with_context(&RunContext::default())
    }

    /// Prepare exact data with cancellation checks between native operations.
    /// Individual native algebra calls retain their existing interruption granularity.
    pub fn prepare_with_context(&self, context: &RunContext) -> Result<PreparedAlgebraicSystem> {
        self.prepare_internal(context, true)
    }

    /// Keep all normalized entries, roots and domain/pole checks, but omit
    /// central-defect row LCMs. The result permits only norm compilation.
    pub(crate) fn prepare_norm_with_context(
        &self,
        context: &RunContext,
    ) -> Result<PreparedAlgebraicNorm> {
        Ok(PreparedAlgebraicNorm {
            inner: self.prepare_internal(context, false)?,
        })
    }

    // A preparation without residuals must stay inside PreparedAlgebraicNorm.
    fn prepare_internal(
        &self,
        context: &RunContext,
        with_residuals: bool,
    ) -> Result<PreparedAlgebraicSystem> {
        context.cancellation.check()?;
        self.validate()?;
        context.cancellation.check()?;
        let root_atoms = self
            .roots
            .iter()
            .map(|r| Atom::var(r.symbol))
            .collect::<Vec<_>>();
        let mut allowed = BTreeSet::from([
            Atom::var(self.system.variable),
            Atom::var(imaginary_parameter()),
        ]);
        allowed.extend(root_atoms.iter().cloned());
        let mut entries = Vec::new();
        let mut normalized = Vec::new();
        let mut domain_entries = Vec::new();
        for (shift, matrix) in self.system.matrices.iter().enumerate() {
            context.cancellation.check()?;
            for (row, values) in matrix.iter().enumerate() {
                context.cancellation.check()?;
                for (column, value) in values.iter().enumerate() {
                    context.cancellation.check()?;
                    let normalized_entry = normalized_root_entry(value, &self.roots, &allowed)?;
                    domain_entries.extend(
                        normalized_entry
                            .inverse_coefficients
                            .into_iter()
                            .map(|coefficient| vec![coefficient]),
                    );
                    for (roots, coefficient) in normalized_entry.terms {
                        entries.push(vec![coefficient]);
                        normalized.push(PreparedKernelTerm {
                            shift,
                            row,
                            column,
                            roots,
                        });
                    }
                }
            }
        }
        // Retain exact row-wise denominator LCMs across all epsilon shifts
        // and root monomials. Ragged rows avoid a dense augmented matrix.
        let (residual_terms, residual_entries) = if with_residuals {
            let mut terms = vec![Vec::<usize>::new(); self.system.matrices[0].len()];
            for (index, term) in normalized.iter().enumerate() {
                terms[term.row].push(index);
            }
            let rows = terms
                .iter()
                .map(|indices| {
                    if indices.is_empty() {
                        vec![Atom::new()]
                    } else {
                        indices.iter().map(|&i| entries[i][0].clone()).collect()
                    }
                })
                .collect::<Vec<_>>();
            (terms, rows)
        } else {
            (Vec::new(), Vec::new())
        };
        let term_count = entries.len();
        let mut root_matrix = vec![vec![Atom::new(); self.roots.len()]; self.roots.len()];
        for (index, root) in self.roots.iter().enumerate() {
            context.cancellation.check()?;
            let r = encode_complex(&root.radicand).together().cancel();
            entries.push(vec![r.clone()]);
            entries.push(vec![Atom::num(1) / &r]); // zeros are branch singularities
            let logarithmic_derivative = (r.derivative(self.system.variable) / (&r * 2)).cancel();
            entries.push(vec![logarithmic_derivative.clone()]);
            root_matrix[index][index] = logarithmic_derivative;
        }
        // Retain norm poles even when multiplication by a numerator cancels them.
        let mut domain_guards = domain_entries
            .iter()
            .map(|row| rational(&row[0], &allowed).map(|r| r.denominator.to_expression()))
            .collect::<Result<Vec<_>>>()?;
        let mut source_conditions = registered_domain_conditions(
            &self
                .system
                .matrices
                .iter()
                .flatten()
                .flatten()
                .cloned()
                .collect::<Vec<_>>(),
            &self.roots,
            &BTreeSet::from([self.system.variable]),
        )?;
        source_conditions.extend(self.nonzero_conditions.iter().cloned());
        for condition in &source_conditions {
            context.cancellation.check()?;
            let condition = rational(condition, &allowed)?;
            for part in [
                condition.numerator.to_expression(),
                condition.denominator.to_expression(),
            ] {
                domain_entries.push(vec![Atom::one() / &part]);
                domain_guards.push(part);
            }
        }
        domain_guards.sort();
        domain_guards.dedup();
        domain_guards.retain(|a| !matches!(a.as_view(), AtomView::Num(_)));
        entries.extend(domain_entries);
        // A zero differential system with no roots still needs a valid native
        // compilation row, while it retains no artificial coupling terms.
        if entries.is_empty() {
            entries.push(vec![Atom::new()]);
        }
        context.cancellation.check()?;
        let residual_entries = if with_residuals {
            Some(PreparedRows::new(
                self.system.variable,
                &residual_entries,
                context,
            )?)
        } else {
            None
        };
        let entries = PreparedRows::new(self.system.variable, &entries, context)?;
        let root_matrix = if root_matrix.is_empty() {
            None
        } else {
            Some(PreparedRows::new(
                self.system.variable,
                &root_matrix,
                context,
            )?)
        };
        context.cancellation.check()?;
        Ok(PreparedAlgebraicSystem {
            source: std::sync::Arc::new(self.clone()),
            entries,
            normalized,
            residual_entries,
            residual_terms,
            term_count,
            root_matrix,
            domain_guards,
        })
    }
}

impl PreparedAlgebraicSystem {
    /// Read the immutable exact snapshot from which all profiles are compiled.
    pub fn source(&self) -> &AlgebraicSystem {
        &self.source
    }

    /// Compile every numerical object freshly at the requested precision.
    pub fn compile(&self, p: Precision) -> Result<CompiledAlgebraicSystem> {
        self.compile_with_context(p, &RunContext::default())
    }

    /// Fresh numerical compilation with checks before and after native compilers.
    pub fn compile_with_context(
        &self,
        p: Precision,
        context: &RunContext,
    ) -> Result<CompiledAlgebraicSystem> {
        if p.bits < 2 {
            return Err(Error::InvalidInput(
                "algebraic precision must be at least two bits".into(),
            ));
        }
        context.cancellation.check()?;
        // Residuals need only cleared polynomials. The entry and root-matrix
        // compilers below retain the full domain and pole checks.
        let mut residual_rows = self
            .residual_entries
            .as_ref()
            .ok_or_else(|| {
                Error::Unsupported("norm-only preparation cannot compile transport".into())
            })?
            .compile_polynomial_rows(p, &Default::default(), context)?;
        for (row, indices) in residual_rows.iter_mut().zip(&self.residual_terms) {
            for (column, _) in &mut row.entries {
                *column = indices[*column];
            }
        }
        self.compile_kernel_with_context(p, context, residual_rows)
    }

    // Shared entry/root specialization, including all original pole scans.
    // Only the norm wrapper may call this with absent central-defect residuals.
    fn compile_kernel_with_context(
        &self,
        p: Precision,
        context: &RunContext,
        residual_rows: Vec<crate::ode::PolynomialRow>,
    ) -> Result<CompiledAlgebraicSystem> {
        if p.bits < 2 {
            return Err(Error::InvalidInput(
                "algebraic precision must be at least two bits".into(),
            ));
        }
        let source = &self.source;
        context.cancellation.check()?;
        let compiled = self.entries.compile(p, &Default::default(), context)?;
        context.cancellation.check()?;
        let terms = self
            .normalized
            .iter()
            .zip(&compiled.matrix)
            .map(|(term, values)| KernelTerm {
                shift: term.shift,
                row: term.row,
                column: term.column,
                roots: term.roots.clone(),
                coefficient: values[0].clone(),
            })
            .collect();
        let roots = source
            .roots
            .iter()
            .enumerate()
            .map(|(i, root)| RootKernel {
                definition: SquareRoot {
                    symbol: root.symbol,
                    radicand: decoded(&root.radicand),
                },
                value: compiled.matrix[self.term_count + 3 * i][0].clone(),
                logarithmic_derivative: compiled.matrix[self.term_count + 3 * i + 2][0].clone(),
            })
            .collect();
        context.cancellation.check()?;
        let root_system = self
            .root_matrix
            .as_ref()
            .map(|rows| rows.compile(p, &Default::default(), context))
            .transpose()?;
        context.cancellation.check()?;
        Ok(CompiledAlgebraicSystem {
            source: self.source.clone(),
            variable: source.system.variable,
            p,
            size: source.system.matrices[0].len(),
            count: source.system.matrices.len(),
            roots,
            root_system,
            terms,
            residual_rows,
            poles: compiled.poles,
            pole_polynomials: compiled.pole_polynomials,
            domain_guards: self.domain_guards.clone(),
        })
    }
}

#[cfg(test)]
mod residual_compile_tests {
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

    #[test]
    fn residual_rows_match_old_full_compiler_across_precision_profiles() -> Result<()> {
        let x = symbol!("algebraic_residual_compile::x");
        let r = symbol!("algebraic_residual_compile::r");
        let t = Atom::var(x);
        let root = Atom::var(r);
        let i = Atom::num(Complex::new(Rational::zero(), Rational::one()));
        let zero = vec![Atom::new(); 3];
        let prepared = AlgebraicSystem {
            system: EpsilonSystem {
                variable: x,
                matrices: vec![
                    vec![
                        vec![
                            (Atom::one() + &i * &t) * &root / (Atom::num(3) * (&t - 1).pow(2)),
                            Atom::one() / (&t + &i),
                            Atom::new(),
                        ],
                        zero.clone(),
                        vec![
                            &root / (&t + 2),
                            Atom::new(),
                            (&t + &i) / (Atom::num(5) * (&t - 1)),
                        ],
                    ],
                    vec![
                        vec![
                            Atom::new(),
                            &root / (Atom::num(7) * (&t + 2)),
                            (Atom::one() + &i) / (&t - 1),
                        ],
                        zero.clone(),
                        vec![
                            Atom::one() / (Atom::num(11) * (&t + 2)),
                            Atom::new(),
                            &root / (&t - 1).pow(2),
                        ],
                    ],
                ],
            },
            roots: vec![SquareRoot {
                symbol: r,
                radicand: (&t + 3) / (&t + 4),
            }],
            nonzero_conditions: vec![&t - 6],
        }
        .prepare()?;
        assert!(prepared.residual_terms[1].is_empty());
        assert!(
            prepared.residual_terms[2]
                .iter()
                .enumerate()
                .any(|(local, &global)| local != global)
        );
        let context = RunContext::default();
        // Revisit a lower profile to exercise fresh specialization of one snapshot.
        for bits in [80, 201, 80] {
            let p = Precision { bits };
            // Freeze the OLD residual compilation path as the reference.
            let old = prepared
                .residual_entries
                .as_ref()
                .unwrap()
                .compile(p, &Default::default(), &context)?
                .polynomial_rows;
            let actual = prepared.compile_with_context(p, &context)?;
            assert_eq!(actual.residual_rows.len(), old.len());
            assert_eq!(old.len(), 3);
            for ((actual, expected), indices) in actual
                .residual_rows
                .iter()
                .zip(&old)
                .zip(&prepared.residual_terms)
            {
                assert_eq!(actual.denominator, expected.denominator);
                let mapped = expected
                    .entries
                    .iter()
                    .map(|(column, coefficients)| (indices[*column], coefficients.clone()))
                    .collect::<Vec<_>>();
                assert_eq!(actual.entries, mapped);
            }
            assert!(actual.residual_rows[1].entries.is_empty());
            let entries = prepared.entries.compile(p, &Default::default(), &context)?;
            assert_same_multiset(&actual.poles, &entries.poles);
            assert_same_multiset(&actual.pole_polynomials, &entries.pole_polynomials);
            assert_eq!(actual.domain_guards, prepared.domain_guards);
            assert_eq!(actual.terms.len(), prepared.normalized.len());
            for ((term, normalized), values) in actual
                .terms
                .iter()
                .zip(&prepared.normalized)
                .zip(&entries.matrix)
            {
                assert_eq!(
                    (term.shift, term.row, term.column, &term.roots),
                    (
                        normalized.shift,
                        normalized.row,
                        normalized.column,
                        &normalized.roots
                    )
                );
                assert_eq!(term.coefficient.numerator, values[0].numerator);
                assert_eq!(term.coefficient.denominator, values[0].denominator);
            }
            let expected_roots =
                prepared
                    .root_matrix
                    .as_ref()
                    .unwrap()
                    .compile(p, &Default::default(), &context)?;
            let actual_roots = actual.root_system.as_ref().unwrap();
            assert_same_multiset(&actual_roots.poles, &expected_roots.poles);
            assert_same_multiset(
                &actual_roots.pole_polynomials,
                &expected_roots.pole_polynomials,
            );
            assert_eq!(
                actual_roots.polynomial_rows.len(),
                expected_roots.polynomial_rows.len()
            );
            for (actual, expected) in actual_roots
                .polynomial_rows
                .iter()
                .zip(&expected_roots.polynomial_rows)
            {
                assert_eq!(actual.denominator, expected.denominator);
                assert_eq!(actual.entries, expected.entries);
            }
        }
        Ok(())
    }

    #[test]
    fn cancelled_kernel_keeps_norm_domain_and_registered_root_poles() -> Result<()> {
        let x = symbol!("algebraic_residual_domain::x");
        let r = symbol!("algebraic_residual_domain::r");
        let t = Atom::var(x);
        let root = Atom::var(r);
        let prepared = AlgebraicSystem {
            system: EpsilonSystem {
                variable: x,
                matrices: vec![vec![vec![(root.clone().pow(2) - &t) / (&root - 1)]]],
            },
            roots: vec![SquareRoot {
                symbol: r,
                radicand: t.clone(),
            }],
            nonzero_conditions: vec![&t - 3],
        }
        .prepare()?;
        // The kernel vanishes under r^2=x; its inverse norm at x=1 must survive.
        assert_eq!(prepared.term_count, 0);
        let p = Precision { bits: 120 };
        let compiled = prepared.compile(p)?;
        assert!(compiled.residual_rows[0].entries.is_empty());
        assert!(compiled.terms.is_empty());
        for pole in [0, 1, 3] {
            assert!(compiled.singularities().contains(&p.i(pole)));
        }
        for excluded in [1, 3] {
            assert!(matches!(
                compiled.check_domain(&p.i(excluded)),
                Err(Error::InvalidInput(_))
            ));
        }
        compiled.check_domain(&p.i(2))?;
        let context = RunContext::default();
        let entries = prepared.entries.compile(p, &Default::default(), &context)?;
        assert_same_multiset(&compiled.poles, &entries.poles);
        assert_same_multiset(&compiled.pole_polynomials, &entries.pole_polynomials);
        let roots =
            prepared
                .root_matrix
                .as_ref()
                .unwrap()
                .compile(p, &Default::default(), &context)?;
        assert!(roots.poles.contains(&p.i(0)));
        Ok(())
    }

    fn prepared_pole_source(gap_bits: u64) -> Result<PreparedAlgebraicSystem> {
        let x = symbol!("algebraic_residual_failures::x");
        let r = symbol!("algebraic_residual_failures::r");
        let t = Atom::var(x);
        let first = &t - 1;
        let gap = Rational::from((Integer::one(), Integer::from(2).pow(gap_bits)));
        let second = &first - Atom::num(gap);
        AlgebraicSystem::ordinary(
            DifferentialSystem {
                variable: x,
                matrix: vec![
                    vec![Atom::one() / &first, Atom::one() / &second],
                    vec![Atom::new(), Atom::new()],
                ],
            },
            vec![SquareRoot {
                symbol: r,
                radicand: &first * &second,
            }],
        )
        .prepare()
    }

    #[test]
    fn entry_and_root_compilers_keep_resolvable_nearby_poles() -> Result<()> {
        // The retained Git79 public probe also accepts the original 2^-60 gap.
        let prepared = prepared_pole_source(60)?;
        let context = RunContext::default();
        let p = Precision { bits: 80 };
        let gap = Rational::from((Integer::one(), Integer::from(2).pow(60)));
        let expected = [p.i(1), p.rational(&(Rational::one() + gap))];
        assert_ne!(expected[0], expected[1]);
        let entries = prepared.entries.compile(p, &Default::default(), &context)?;
        let roots =
            prepared
                .root_matrix
                .as_ref()
                .unwrap()
                .compile(p, &Default::default(), &context)?;
        let full = prepared.compile_with_context(p, &context)?;
        assert_same_multiset(&entries.poles, &expected);
        assert_same_multiset(&roots.poles, &expected);
        assert_same_multiset(full.singularities(), &expected);
        Ok(())
    }

    #[test]
    fn row_factor_reuse_lcm_then_whole_entry_keeps_specialized_shared_poles() -> Result<()> {
        let x = symbol!("algebraic_row_factor_cache::x");
        let a = Atom::var(symbol!("algebraic_row_factor_cache::a"));
        let t = Atom::var(x);
        let f = &t - &a;
        let g = &t + &a;
        let h = &t + 3;
        let denominator = Atom::num(6) * f.clone().pow(2) * g.clone().pow(3) * h.clone().pow(2);
        // Row zero certifies a multi-entry LCM. Later entries encounter that
        // entire denominator in the same factor cache, including its content.
        let rows = PreparedRows::new(
            x,
            &[
                vec![
                    Atom::one() / (Atom::num(2) * f.clone().pow(2) * &h),
                    Atom::one() / (Atom::num(3) * g.clone().pow(3) * h.clone().pow(2)),
                ],
                vec![Atom::one() / &denominator, Atom::num(5) / &denominator],
            ],
            &RunContext::default(),
        )?;
        let context = RunContext::default();
        let p = Precision { bits: 80 };
        for value in [0, 1, 0] {
            let values = [(a.clone(), p.i(value))].into_iter().collect();
            let compiled = rows.compile(p, &values, &context)?;
            let expected = if value == 0 {
                vec![p.i(0), p.i(-3)]
            } else {
                vec![p.i(-1), p.i(1), p.i(-3)]
            };
            assert_same_multiset(&compiled.poles, &expected);
            // Specialization merges f and g at zero; all five powers survive.
            let expected_denominator = if value == 0 {
                Atom::num(6) * t.clone().pow(5) * h.clone().pow(2)
            } else {
                Atom::num(6) * (&t - 1).pow(2) * (&t + 1).pow(3) * h.clone().pow(2)
            };
            let reference = PreparedRows::new(
                x,
                &[vec![
                    Atom::one() / &expected_denominator,
                    Atom::num(5) / &expected_denominator,
                ]],
                &context,
            )?
            .compile(p, &Default::default(), &context)?;
            assert_eq!(
                compiled.polynomial_rows[1].denominator,
                reference.polynomial_rows[0].denominator
            );
            assert_eq!(
                compiled.polynomial_rows[1].entries,
                reference.polynomial_rows[0].entries
            );
            assert_eq!(
                compiled.polynomial_rows[0].denominator,
                compiled.polynomial_rows[1].denominator
            );
        }
        Ok(())
    }

    #[test]
    fn entry_and_root_compilers_keep_pole_failures_and_cancellation() -> Result<()> {
        let prepared = prepared_pole_source(5000)?;
        let context = RunContext::default();
        let p = Precision { bits: 80 };
        // Finite residual polynomials do not waive unresolved-pole rejection.
        prepared
            .residual_entries
            .as_ref()
            .unwrap()
            .compile_polynomial_rows(p, &Default::default(), &context)?;
        assert!(matches!(
            prepared.entries.compile(p, &Default::default(), &context),
            Err(Error::Numerical(message)) if message.contains("ambiguous nearby poles")
        ));
        let root_rows = prepared.root_matrix.as_ref().unwrap();
        assert!(matches!(
            root_rows.compile(p, &Default::default(), &context),
            Err(Error::Numerical(message)) if message.contains("ambiguous nearby poles")
        ));
        assert!(matches!(
            prepared.compile_with_context(p, &context),
            Err(Error::Numerical(message)) if message.contains("ambiguous nearby poles")
        ));
        assert!(matches!(
            prepared.compile_with_context(Precision { bits: 1 }, &context),
            Err(Error::InvalidInput(_))
        ));
        context.cancellation.cancel();
        assert!(matches!(
            prepared.entries.compile(p, &Default::default(), &context),
            Err(Error::Cancelled)
        ));
        assert!(matches!(
            root_rows.compile(p, &Default::default(), &context),
            Err(Error::Cancelled)
        ));
        assert!(matches!(
            prepared.compile_with_context(p, &context),
            Err(Error::Cancelled)
        ));
        Ok(())
    }
}

#[cfg(test)]
#[path = "norm_only_tests.rs"]
mod norm_only_tests;
