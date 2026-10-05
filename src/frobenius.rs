#![allow(clippy::needless_range_loop)] // Matrix row/column indexing mirrors the recurrence equations.
//! Generalized power/logarithm series at regular singular points.
use crate::algebra;
use crate::family::substitute;
use crate::numeric::solve;
use crate::ode::{polynomial_coefficients, quotient_series};
use crate::{
    ComplexFloat as C, DifferentialSystem, Error, Precision, Progress, Result, RunContext,
};
use std::collections::BTreeMap;
use symbolica::prelude::*;

#[derive(Clone, Debug)]
pub struct FrobeniusColumn {
    /// Symbolic exponent, retained to distinguish dimensional-regulator regions.
    pub exponent: Atom,
    /// coefficients[power][log_power][component]
    pub coefficients: Vec<Vec<Vec<C>>>,
}

#[derive(Clone, Debug)]
pub struct FrobeniusBasis {
    pub columns: Vec<FrobeniusColumn>,
    pub precision: Precision,
}

/// Parameter-independent, exact preparation of a regular-singular system.
///
/// Normalization, indicial factorization and generalized eigenspaces are computed
/// once with Symbolica's exact algebra. Numerical evaluations can then use any
/// precision, admissible generic parameter sample and truncation order without specializing the
/// symbolic exponents or sharing mutable numerical workspaces.
#[derive(Clone, Debug)]
pub struct PreparedFrobenius {
    source: DifferentialSystem,
    normal: DifferentialSystem,
    transformation: Vec<Vec<Atom>>,
    blocks: Vec<Vec<usize>>,
    eigenspaces: Vec<PreparedEigenspace>,
    nonzero_conditions: Vec<Atom>,
}

#[derive(Clone, Debug)]
struct PreparedEigenspace {
    exponent: Atom,
    /// columns[generalized eigenvector][log power][component]
    leading: Vec<Vec<Vec<Atom>>>,
}

fn phase(context: &RunContext, name: impl Into<String>) -> Result<()> {
    context.emit(Progress::Stage { name: name.into() })?;
    // A progress callback can request cancellation before the next exact owner call.
    context.cancellation.check()
}

/// Reject parameter specialization that changes the symbolic resonance pattern.
pub(crate) fn ensure_generic_exponents(
    exponents: &[&Atom],
    p: Precision,
    parameters: &ahash::HashMap<Atom, C>,
) -> Result<()> {
    for (i, a) in exponents.iter().enumerate() {
        for b in &exponents[..i] {
            let difference = (*a - *b).together().cancel();
            if difference.to_string().parse::<i64>().is_ok() {
                continue;
            }
            let value = p.eval(&difference, parameters)?;
            if let Some(integer) = value.re.as_raw().to_integer()
                && let Ok(integer) = integer.to_string().parse::<i64>()
                && p.close(&value, &p.i(integer), p.bits / 5)
            {
                return Err(Error::Unsupported("parameter specialization introduces an additional indicial resonance; choose a generic sample or prepare the exactly specialized system".into()));
            }
        }
    }
    Ok(())
}

pub(crate) fn valuation(a: &Atom, x: Symbol) -> Result<i64> {
    if a.is_zero() {
        return Ok(i64::MAX / 4);
    }
    let s = a
        .series(x, 0, 1)
        .map_err(|e| Error::Unsupported(e.to_string()))?;
    s.get_trailing_exponent()
        .to_string()
        .parse()
        .map_err(|_| Error::Unsupported("noninteger matrix valuation".into()))
}

impl DifferentialSystem {
    /// Remove higher poles when a diagonal integer-power shearing suffices.
    /// Returns shifts s such that original_y[i] = x^s[i] * normalized_y[i].
    pub fn diagonal_fuchsian_form(&self) -> Result<(Self, Vec<i64>)> {
        self.validate()?;
        let n = self.matrix.len();
        let orders = self
            .matrix
            .iter()
            .map(|row| {
                row.iter()
                    .map(|a| {
                        if a.is_zero() {
                            Ok(None)
                        } else {
                            valuation(a, self.variable).map(Some)
                        }
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let mut shifts = crate::integer_shearing::integer_potentials(&orders, 1, 10000, None)?;
        let minimum = *shifts.iter().min().unwrap();
        for s in &mut shifts {
            *s -= minimum;
        }
        if shifts.iter().any(|&s| s > 10000) {
            return Err(Error::Limit(
                "Fuchsian shearing exceeds 10000 powers".into(),
            ));
        }
        let x = Atom::var(self.variable);
        let transformation = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        if i == j {
                            x.clone().pow(shifts[i])
                        } else {
                            Atom::new()
                        }
                    })
                    .collect()
            })
            .collect::<Vec<Vec<_>>>();
        Ok((self.change_basis(&transformation)?, shifts))
    }
    /// Constant similarity built from the kernel filtration of the leading
    /// nilpotent coefficient. This exposes hidden triangular higher poles
    /// before integer-power shearing.
    pub fn nilpotent_leading_rotation(&self) -> Result<Vec<Vec<Atom>>> {
        self.validate()?;
        let n = self.matrix.len();
        let order = self
            .matrix
            .iter()
            .flatten()
            .map(|a| valuation(a, self.variable))
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .min()
            .unwrap();
        if order >= -1 {
            return Err(Error::InvalidInput("rotation needs a higher pole".into()));
        }
        let leading = self
            .matrix
            .iter()
            .map(|r| {
                r.iter()
                    .map(|a| {
                        Ok(a.series(self.variable, 0, order + 1)
                            .map_err(|e| Error::Unsupported(e.to_string()))?
                            .coefficient(Rational::from(order))
                            .unwrap_or_default())
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let mut power = leading.clone();
        let mut vectors = Vec::new();
        for _ in 0..n {
            for v in algebra::nullspace(power.clone()) {
                let mut trial = vectors.clone();
                trial.push(v.clone());
                if algebra::rref(trial).1.len() > vectors.len() {
                    vectors.push(v);
                }
            }
            if vectors.len() == n {
                break;
            }
            power = algebra::matmul(&power, &leading);
        }
        if vectors.len() != n {
            return Err(Error::Unsupported(
                "nonnilpotent irregular singularity".into(),
            ));
        }
        Ok((0..n)
            .map(|i| vectors.iter().map(|v| v[i].clone()).collect())
            .collect())
    }

    pub fn invert_variable(&self, new_variable: Symbol) -> Self {
        let z = Atom::var(new_variable);
        let rules = BTreeMap::from([(Atom::var(self.variable), z.clone().pow(-1))]);
        Self {
            variable: new_variable,
            matrix: self
                .matrix
                .iter()
                .map(|row| {
                    row.iter()
                        // Preserve original denominator domains until preparation
                        // collects them; normalization can cancel afterward.
                        .map(|a| -substitute(a, &rules) / z.clone().pow(2))
                        .collect()
                })
                .collect(),
        }
    }
    /// Exact residue for a system already in regular-singular form.
    pub fn residue(&self) -> Result<Vec<Vec<Atom>>> {
        self.residue_with_context(&RunContext::default())
    }

    fn residue_with_context(&self, context: &RunContext) -> Result<Vec<Vec<Atom>>> {
        self.validate()?;
        self.matrix
            .iter()
            .map(|row| {
                context.cancellation.check()?;
                row.iter()
                    .map(|a| {
                        if valuation(a, self.variable)? < -1 {
                            return Err(Error::Unsupported(
                                "matrix needs a nontrivial Fuchsian basis transformation".into(),
                            ));
                        }
                        let expansion = a
                            .series(self.variable, 0, 1)
                            .map_err(|e| Error::Unsupported(e.to_string()))?;
                        Ok(expansion
                            .coefficient(Rational::from(-1))
                            .unwrap_or_default())
                    })
                    .collect()
            })
            .collect()
    }

    /// Perform the exact work shared by numerical epsilon samples and precision
    /// refinements. Failed or cancelled preparations return no partial object.
    pub fn prepare_frobenius(&self, context: &RunContext) -> Result<PreparedFrobenius> {
        phase(context, "normalizing the exact Frobenius system")?;
        self.validate()?;
        let expressions = self.matrix.iter().flatten().cloned().collect::<Vec<_>>();
        let variables = expressions
            .iter()
            .flat_map(|a| a.get_all_symbols(true))
            .collect();
        // Collect before normalization so cancellable factors retain their domain.
        let nonzero_conditions =
            crate::physical_conditions::rational_denominator_conditions(&expressions, &variables)?;
        let (normal, transformation) = self.fuchsian_form(32)?;
        phase(context, "extracting the exact Frobenius residue")?;
        let residue = normal.residue_with_context(context)?;
        let n = residue.len();
        let blocks = normal.blocks()?;
        let mut roots = BTreeMap::new();
        for (index, block) in blocks.iter().enumerate() {
            phase(
                context,
                format!(
                    "factoring exact Frobenius indicial block {}/{}",
                    index + 1,
                    blocks.len()
                ),
            )?;
            let diagonal = block
                .iter()
                .map(|&i| block.iter().map(|&j| residue[i][j].clone()).collect())
                .collect::<Vec<_>>();
            for (root, multiplicity) in algebra::eigenvalues(&diagonal)? {
                *roots.entry(root).or_insert(0) += multiplicity;
            }
        }
        let mut eigenspaces = Vec::new();
        for (index, (lambda, multiplicity)) in roots.iter().enumerate() {
            phase(
                context,
                format!(
                    "constructing exact Frobenius eigenspace {}/{}",
                    index + 1,
                    roots.len()
                ),
            )?;
            let mut nilpotent = residue.clone();
            for (i, row) in nilpotent.iter_mut().enumerate() {
                context.cancellation.check()?;
                row[i] = (&row[i] - lambda).together().cancel();
            }
            let mut power = nilpotent.clone();
            let mut generalized = Vec::new();
            // Once the kernel reaches the algebraic multiplicity it is the full
            // generalized eigenspace. Large repeated eigenvalues need not have
            // equally long Jordan chains.
            for depth in 1..=*multiplicity {
                context.cancellation.check()?;
                generalized = algebra::nullspace(power.clone());
                if generalized.len() == *multiplicity {
                    break;
                }
                if depth < *multiplicity {
                    context.cancellation.check()?;
                    power = algebra::matmul(&power, &nilpotent);
                }
            }
            if generalized.len() != *multiplicity {
                return Err(Error::Numerical(
                    "generalized eigenspace dimension mismatch".into(),
                ));
            }
            let mut leading = Vec::new();
            for vector in generalized {
                let mut column = Vec::new();
                let mut v = vector;
                for k in 0..*multiplicity {
                    context.cancellation.check()?;
                    column.push(v.clone());
                    v = (0..n)
                        .map(|i| {
                            context.cancellation.check()?;
                            Ok(
                                ((0..n).fold(Atom::new(), |s, j| s + &nilpotent[i][j] * &v[j])
                                    / Atom::num((k + 1) as i64))
                                .together()
                                .cancel(),
                            )
                        })
                        .collect::<Result<Vec<_>>>()?;
                    if v.iter().all(|entry| entry.is_zero()) {
                        break;
                    }
                }
                leading.push(column);
            }
            eigenspaces.push(PreparedEigenspace {
                exponent: lambda.clone(),
                leading,
            });
        }
        context.cancellation.check()?;
        Ok(PreparedFrobenius {
            source: self.clone(),
            normal,
            transformation,
            blocks,
            eigenspaces,
            nonzero_conditions,
        })
    }

    /// Construct all power/logarithm solutions. Rational indicial roots are kept
    /// exact; integer resonances use a coupled logarithmic recurrence.
    /// For repeated samples, use [`Self::prepare_frobenius`] once instead.
    pub fn frobenius(
        &self,
        p: Precision,
        values: &ahash::HashMap<Atom, C>,
        order: usize,
    ) -> Result<FrobeniusBasis> {
        let context = RunContext::default();
        self.prepare_frobenius(&context)?
            .evaluate(p, values, order, &context)
    }
}

impl PreparedFrobenius {
    /// Original exact connection, including expressions before rational cancellation.
    pub fn system(&self) -> &DifferentialSystem {
        &self.source
    }

    /// Fuchsian connection used by the recurrence.
    pub fn normalized_system(&self) -> &DifferentialSystem {
        &self.normal
    }

    /// Exact transformation with original_y = T * normalized_y.
    pub fn basis_transformation(&self) -> &[Vec<Atom>] {
        &self.transformation
    }

    /// Original rational denominator restrictions. A singular expansion center
    /// is allowed; a parameter sample must not annihilate a whole restriction.
    pub fn nonzero_conditions(&self) -> &[Atom] {
        &self.nonzero_conditions
    }

    /// Distinct exact indicial exponents before the final basis transformation.
    pub fn exponents(&self) -> impl ExactSizeIterator<Item = &Atom> {
        self.eigenspaces.iter().map(|space| &space.exponent)
    }

    /// Evaluate with a fresh numerical workspace. The exact preparation can be
    /// shared concurrently across epsilon samples, precisions and series orders.
    pub fn evaluate(
        &self,
        p: Precision,
        values: &ahash::HashMap<Atom, C>,
        order: usize,
        context: &RunContext,
    ) -> Result<FrobeniusBasis> {
        phase(context, "evaluating the prepared Frobenius recurrence")?;
        if order > 10_000 {
            return Err(Error::Limit(
                "Frobenius expansion exceeds 10000 terms".into(),
            ));
        }
        for condition in &self.nonzero_conditions {
            context.cancellation.check()?;
            let coefficients = polynomial_coefficients(condition, self.source.variable, p, values)?;
            if coefficients
                .iter()
                .all(|coefficient| coefficient == &p.zero())
            {
                return Err(Error::Numerical(
                    "parameter sample annihilates an original Frobenius denominator condition"
                        .into(),
                ));
            }
        }
        ensure_generic_exponents(&self.exponents().collect::<Vec<_>>(), p, values)?;
        let mut basis = self.evaluate_normal(p, values, order, context)?;
        let mut terms = Vec::new();
        let mut low = 0_i64;
        let mut high = 0_i64;
        for (i, row) in self.transformation.iter().enumerate() {
            context.cancellation.check()?;
            for (j, entry) in row.iter().enumerate() {
                if entry.is_zero() {
                    continue;
                }
                let series = entry
                    .series(self.source.variable, 0, order as i64 + 1)
                    .map_err(|e| Error::Unsupported(e.to_string()))?;
                for (power, coefficient) in series.terms() {
                    let shift = power.to_string().parse::<i64>().map_err(|_| {
                        Error::Unsupported("noninteger basis transformation power".into())
                    })?;
                    low = low.min(shift);
                    high = high.max(shift);
                    terms.push((i, j, shift, p.eval(coefficient, values)?));
                }
            }
        }
        if high - low > 256 {
            return Err(Error::Limit("normalization power span exceeds 256".into()));
        }
        for column in &mut basis.columns {
            context.cancellation.check()?;
            let mut transformed = vec![
                vec![vec![p.zero(); self.source.matrix.len()]];
                column.coefficients.len() + (high - low) as usize
            ];
            for (k, logs) in column.coefficients.iter().enumerate() {
                for (l, row) in logs.iter().enumerate() {
                    for (i, j, shift, coefficient) in &terms {
                        let destination = &mut transformed[(k as i64 + shift - low) as usize];
                        destination.resize(
                            destination.len().max(l + 1),
                            vec![p.zero(); self.source.matrix.len()],
                        );
                        destination[l][*i] =
                            p.add(&destination[l][*i], &p.mul(coefficient, &row[*j]));
                    }
                }
            }
            column.exponent = (&column.exponent + Atom::num(low)).together().cancel();
            column.coefficients = transformed;
        }
        Ok(basis)
    }

    fn evaluate_normal(
        &self,
        p: Precision,
        values: &ahash::HashMap<Atom, C>,
        order: usize,
        context: &RunContext,
    ) -> Result<FrobeniusBasis> {
        let n = self.normal.matrix.len();
        let mut a = vec![vec![vec![p.zero(); n]; n]; order + 1];
        for (i, row) in self.normal.matrix.iter().enumerate() {
            context.cancellation.check()?;
            for (j, expr) in row.iter().enumerate() {
                let expr = (expr * Atom::var(self.source.variable)).together().cancel();
                let rat: RationalPolynomial<IntegerRing, u16> = expr
                    .try_to_rational_polynomial(&Q, &Z, None)
                    .map_err(|e| Error::InvalidInput(e.to_string()))?;
                let num = polynomial_coefficients(
                    &rat.numerator.to_expression(),
                    self.source.variable,
                    p,
                    values,
                )?;
                let den = polynomial_coefficients(
                    &rat.denominator.to_expression(),
                    self.source.variable,
                    p,
                    values,
                )?;
                let series = quotient_series(p, &num, &den, order)?;
                for k in 0..=order {
                    a[k][i][j] = series[k].clone();
                }
            }
        }
        let zero = p.zero();
        let sparse = a
            .iter()
            .map(|matrix| {
                matrix
                    .iter()
                    .enumerate()
                    .flat_map(|(i, row)| {
                        row.iter()
                            .enumerate()
                            .filter_map(|(j, value)| (value != &zero).then_some((i, j, value)))
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut columns = Vec::new();
        for space in &self.eigenspaces {
            context.cancellation.check()?;
            let lambda = &space.exponent;
            let exponent = p.eval(lambda, values)?;
            let mut solvers = (0..=order)
                .map(|_| None)
                .collect::<Vec<Option<crate::numeric::BlockSolve>>>();
            for vector in &space.leading {
                context.cancellation.check()?;
                let mut leading = vector
                    .iter()
                    .map(|log| {
                        log.iter()
                            .map(|x| p.eval(x, values))
                            .collect::<Result<Vec<_>>>()
                    })
                    .collect::<Result<Vec<_>>>()?;
                while leading.len() > 1 && leading.last().unwrap().iter().all(|v| v == &zero) {
                    leading.pop();
                }
                let mut coefficients = vec![leading];
                for k in 1..=order {
                    context.cancellation.check()?;
                    let mut m = a[0]
                        .iter()
                        .map(|r| r.iter().map(|v| p.neg(v)).collect::<Vec<_>>())
                        .collect::<Vec<_>>();
                    for (i, row) in m.iter_mut().enumerate() {
                        row[i] = p.add(&row[i], &p.add(&exponent, &p.i(k as i64)));
                    }
                    let mut resonance = false;
                    for other in self.eigenspaces.iter().map(|space| &space.exponent) {
                        let difference = (lambda + Atom::num(k as i64) - other).together().cancel();
                        if difference.is_zero()
                            || p.norm(&p.eval(&difference, values)?) < p.tolerance(p.bits / 4)
                        {
                            resonance = true;
                            break;
                        }
                    }
                    let logs = coefficients.iter().map(|c| c.len()).max().unwrap();
                    let mut rhs = vec![vec![p.zero(); n]; logs];
                    for j in 1..=k {
                        for (l, previous) in coefficients[k - j].iter().enumerate() {
                            for &(i, column, value) in &sparse[j] {
                                if previous[column] != zero {
                                    rhs[l][i] = p.add(&rhs[l][i], &p.mul(value, &previous[column]));
                                }
                            }
                        }
                    }
                    let mut next = vec![vec![p.zero(); n]; logs];
                    if resonance {
                        let mut result = None;
                        for extra in 0..=n {
                            context.cancellation.check()?;
                            let count = logs + extra;
                            let size = n * count;
                            let mut system = vec![vec![p.zero(); size]; size];
                            for l in 0..count {
                                for i in 0..n {
                                    for j in 0..n {
                                        system[l * n + i][l * n + j] = m[i][j].clone();
                                    }
                                    if l + 1 < count {
                                        system[l * n + i][(l + 1) * n + i] = p.i((l + 1) as i64);
                                    }
                                }
                            }
                            let mut extended = rhs.clone();
                            extended.resize(count, vec![p.zero(); n]);
                            match solve_any(
                                p,
                                system,
                                extended.into_iter().flatten().collect(),
                                context,
                            ) {
                                Ok(flat) => {
                                    result = Some(flat.chunks(n).map(|r| r.to_vec()).collect());
                                    break;
                                }
                                Err(Error::Numerical(_)) => continue,
                                Err(e) => return Err(e),
                            }
                        }
                        next = result.ok_or_else(|| {
                            Error::Numerical(
                                "resonant recurrence did not close with sufficient logarithms"
                                    .into(),
                            )
                        })?;
                    } else {
                        if solvers[k].is_none() {
                            solvers[k] =
                                Some(crate::numeric::BlockSolve::new(p, m.clone(), &self.blocks)?);
                        }
                        for l in (0..logs).rev() {
                            if l + 1 < logs {
                                for i in 0..n {
                                    rhs[l][i] = p.sub(
                                        &rhs[l][i],
                                        &p.scale(&next[l + 1][i], (l + 1) as i64, 1),
                                    );
                                }
                            }
                            next[l] = solvers[k].as_ref().unwrap().solve(&rhs[l]);
                        }
                    }
                    while next.len() > 1
                        && next
                            .last()
                            .unwrap()
                            .iter()
                            .all(|v| p.norm(v) < p.tolerance(p.bits / 4))
                    {
                        next.pop();
                    }
                    coefficients.push(next);
                }
                columns.push(FrobeniusColumn {
                    exponent: lambda.clone(),
                    coefficients,
                });
            }
        }
        Ok(FrobeniusBasis {
            columns,
            precision: p,
        })
    }
}

/// RREF with free parameters set to zero, used only at exact detected resonances.
fn solve_any(
    p: Precision,
    mut a: Vec<Vec<C>>,
    mut b: Vec<C>,
    context: &RunContext,
) -> Result<Vec<C>> {
    let n = b.len();
    let mut row = 0;
    let mut pivots = Vec::new();
    let tol = p.tolerance(p.bits / 4);
    for col in 0..n {
        context.cancellation.check()?;
        let Some(pivot) =
            (row..n).max_by(|&i, &j| p.norm(&a[i][col]).partial_cmp(&p.norm(&a[j][col])).unwrap())
        else {
            break;
        };
        if p.norm(&a[pivot][col]) < tol {
            continue;
        }
        a.swap(row, pivot);
        b.swap(row, pivot);
        let v = a[row][col].clone();
        for j in col..n {
            a[row][j] = p.div(&a[row][j], &v);
        }
        b[row] = p.div(&b[row], &v);
        for i in 0..n {
            if i != row && a[i][col] != p.zero() {
                let f = a[i][col].clone();
                for j in col..n {
                    a[i][j] = p.sub(&a[i][j], &p.mul(&f, &a[row][j]));
                }
                b[i] = p.sub(&b[i], &p.mul(&f, &b[row]));
            }
        }
        pivots.push(col);
        row += 1;
        if row == n {
            break;
        }
    }
    if b.iter().skip(row).any(|v| p.norm(v) > tol) {
        return Err(Error::Numerical(
            "inconsistent resonant Frobenius recurrence".into(),
        ));
    }
    let mut result = vec![p.zero(); n];
    for (i, &j) in pivots.iter().enumerate() {
        result[j] = b[i].clone();
    }
    Ok(result)
}

impl FrobeniusBasis {
    pub(crate) fn validate(&self) -> Result<()> {
        let n = self.columns.len();
        if n == 0
            || self.columns.iter().any(|c| {
                c.coefficients.is_empty()
                    || c.coefficients
                        .iter()
                        .any(|logs| logs.is_empty() || logs.iter().any(|row| row.len() != n))
            })
        {
            return Err(Error::InvalidInput(
                "Frobenius coefficient dimensions".into(),
            ));
        }
        if self
            .columns
            .iter()
            .flat_map(|c| &c.coefficients)
            .flatten()
            .flatten()
            .any(|value| !self.precision.finite(value))
        {
            return Err(Error::Numerical("nonfinite Frobenius coefficient".into()));
        }
        Ok(())
    }
    pub fn evaluate(&self, z: &C, values: &ahash::HashMap<Atom, C>) -> Result<Vec<Vec<C>>> {
        self.validate()?;
        let p = self.precision;
        if !p.finite(z) || *z == p.zero() {
            return Err(Error::InvalidInput(
                "evaluate Frobenius series at a nonzero matching point".into(),
            ));
        }
        let log = p.log(z);
        let n = self.columns.len();
        let mut matrix = vec![vec![p.zero(); n]; n];
        for (j, column) in self.columns.iter().enumerate() {
            let exponent = p.eval(&column.exponent, values)?;
            let mut power = p.pow(z, &exponent);
            for coeff in &column.coefficients {
                let mut logpower = p.i(1);
                for row in coeff {
                    for i in 0..n {
                        matrix[i][j] =
                            p.add(&matrix[i][j], &p.mul(&p.mul(&power, &logpower), &row[i]));
                    }
                    logpower = p.mul(&logpower, &log);
                }
                power = p.mul(&power, z);
            }
        }
        Ok(matrix)
    }
    pub fn match_values(
        &self,
        z: &C,
        values: &[C],
        parameters: &ahash::HashMap<Atom, C>,
    ) -> Result<Vec<C>> {
        solve(
            self.precision,
            self.evaluate(z, parameters)?,
            values.to_vec(),
        )
    }
    /// Select the constant Taylor coefficient in dimensional regularization.
    /// Nonconstant epsilon-dependent powers are discarded before epsilon fitting.
    pub fn physical_limit(&self, constants: &[C], epsilon: Symbol) -> Result<Vec<C>> {
        self.validate()?;
        let n = self.columns.len();
        (0..n)
            .map(|i| {
                let weights = (0..n)
                    .map(|j| Atom::num(i64::from(i == j)))
                    .collect::<Vec<_>>();
                crate::engine::project_limit(
                    self,
                    constants,
                    &weights,
                    epsilon,
                    symbol!("symbolica_amflow::endpoint_projection"),
                    &ahash::HashMap::default(),
                )
            })
            .collect()
    }
}
