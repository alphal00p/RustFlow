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
    residual_entries: PreparedRows,
    residual_terms: Vec<Vec<usize>>,
    term_count: usize,
    root_matrix: Option<PreparedRows>,
    domain_guards: Vec<Atom>,
}

impl AlgebraicSystem {
    /// Snapshot and normalize the exact system for repeated precision profiles.
    pub fn prepare(&self) -> Result<PreparedAlgebraicSystem> {
        self.prepare_with_context(&RunContext::default())
    }

    /// Prepare exact data with cancellation checks between native operations.
    /// Individual native algebra calls retain their existing interruption granularity.
    pub fn prepare_with_context(&self, context: &RunContext) -> Result<PreparedAlgebraicSystem> {
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
        let mut residual_terms = vec![Vec::<usize>::new(); self.system.matrices[0].len()];
        for (index, term) in normalized.iter().enumerate() {
            residual_terms[term.row].push(index);
        }
        let residual_entries = residual_terms
            .iter()
            .map(|indices| {
                if indices.is_empty() {
                    vec![Atom::new()]
                } else {
                    indices.iter().map(|&i| entries[i][0].clone()).collect()
                }
            })
            .collect::<Vec<_>>();
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
        let residual_entries = PreparedRows::new(self.system.variable, &residual_entries, context)?;
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
        let source = &self.source;
        let mut residual_rows = self
            .residual_entries
            .compile(p, &Default::default(), context)?
            .polynomial_rows;
        for (row, indices) in residual_rows.iter_mut().zip(&self.residual_terms) {
            for (column, _) in &mut row.entries {
                *column = indices[*column];
            }
        }
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
