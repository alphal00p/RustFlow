//! Exact coefficient operations for the shared Frobenius recurrence.
pub(crate) use super::recurrence::RecurrenceColumn;
use super::recurrence::*;
use super::*;
use symbolica::domains::float::FloatField;
use symbolica::tensors::matrix::Matrix;
use symbolica::utils::resources as process_resources;
mod resources;
mod sparse_policy;
#[cfg(test)]
mod sparse_policy_tests;
pub(crate) use resources::preflight_operation;
#[cfg(test)]
mod series_tests;

pub(crate) type Gaussian = Complex<Rational>;
type ExactMatrix = Matrix<FloatField<Gaussian>>;

/// Native elimination policy for exact Frobenius coefficients and their constraints.
///
/// The default retains conservative dense preflight estimates. The opt-in sparse
/// policy uses native sparse inversion with input/output height checks. For
/// resonant solves and constraints it also observes native reducer heights after
/// every inserted row and after back substitution. Inverse internals and arithmetic
/// inside a native row operation are unobserved; transient allocations are not bounded.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ExactEliminationMode {
    #[default]
    ConservativeDense,
    ObservedSparse,
}

/// Optional diagnostic caps for exact coefficient work, not allocator quotas.
#[derive(
    Clone, Debug, Default, serde::Serialize, serde::Deserialize, bincode::Encode, bincode::Decode,
)]
pub struct ExactFrobeniusLimits {
    pub max_dimension: Option<usize>,
    pub max_order: Option<usize>,
    pub max_coefficient_bits: Option<u64>,
    pub max_scalar_cells: Option<usize>,
}
impl ExactFrobeniusLimits {
    pub(crate) fn validate(&self) -> Result<()> {
        if self.max_dimension == Some(0)
            || self.max_order == Some(0)
            || self.max_coefficient_bits == Some(0)
            || self.max_scalar_cells == Some(0)
        {
            return Err(Error::InvalidInput(
                "invalid exact Frobenius resource limits".into(),
            ));
        }
        Ok(())
    }
}

fn scalar_count(count: Option<usize>) -> Result<usize> {
    count.ok_or_else(|| Error::Limit("exact Frobenius scalar count overflow".into()))
}

fn coefficient_bits(bits: Option<u64>, limits: &ExactFrobeniusLimits) -> Result<u64> {
    match bits {
        Some(bits) => Ok(bits),
        // This is a conservative admission estimate, not an index or a size
        // consumed by the exact algorithm. Disabled diagnostics cannot reject
        // compact truncated output because the estimate exceeds its counter.
        None if limits.max_coefficient_bits.is_none() => Ok(u64::MAX),
        None => Err(Error::Limit(
            "exact Frobenius coefficient-height estimate overflow".into(),
        )),
    }
}

fn allocation_bytes(count: usize, item_bytes: usize) -> Result<u64> {
    let bytes = count
        .checked_mul(item_bytes)
        .ok_or_else(|| Error::Limit("exact Frobenius allocation size overflow".into()))?;
    u64::try_from(bytes)
        .map_err(|_| Error::Limit("exact Frobenius allocation size exceeds u64".into()))
}

fn series_terms(order: usize) -> Result<usize> {
    scalar_count(order.checked_add(1))
}

fn native_series_order(order: usize) -> Result<i64> {
    i64::try_from(series_terms(order)?)
        .map_err(|_| Error::Limit("exact Frobenius native series order exceeds i64".into()))
}

/// A source-derived exact prefix. Padding from normalization does not extend
/// `known_order`; every column's uniformly known prefix has that same length.
pub struct ExactFrobeniusBasis {
    pub(crate) columns: Vec<RecurrenceColumn<Gaussian>>,
    pub(crate) known_order: usize,
    pub(crate) limits: ExactFrobeniusLimits,
    pub(crate) elimination_mode: ExactEliminationMode,
}
impl ExactFrobeniusBasis {
    /// Elimination policy retained by this exact basis for subsequent constraints.
    pub fn elimination_mode(&self) -> ExactEliminationMode {
        self.elimination_mode
    }
    pub fn dimension(&self) -> usize {
        self.columns.len()
    }
    pub fn known_order(&self) -> usize {
        self.known_order
    }
    /// Round exact coefficients for the existing branch-aware basis evaluator.
    /// This conversion makes no assertion about truncation or input accuracy.
    pub fn numerical(&self, p: Precision) -> FrobeniusBasis {
        FrobeniusBasis {
            columns: self
                .columns
                .iter()
                .map(|column| FrobeniusColumn {
                    exponent: column.exponent.clone(),
                    coefficients: column
                        .coefficients
                        .iter()
                        .map(|logs| {
                            logs.iter()
                                .map(|row| {
                                    row.iter()
                                        .map(|a| {
                                            Complex::new(
                                                a.re.to_multi_prec_float(p.bits),
                                                a.im.to_multi_prec_float(p.bits),
                                            )
                                        })
                                        .collect()
                                })
                                .collect()
                        })
                        .collect(),
                })
                .collect(),
            precision: p,
        }
    }
}

pub(crate) fn gaussian(expression: &Atom) -> Result<Gaussian> {
    let expression = crate::family::substitute(
        expression,
        &BTreeMap::from([(
            Atom::var(crate::family::imaginary_parameter()),
            Atom::num(Complex::new(Rational::zero(), Rational::one())),
        )]),
    )
    .together()
    .cancel();
    if let AtomView::Num(number) = expression.as_view()
        && let symbolica::coefficient::Coefficient::Complex(value) =
            number.get_coeff_view().to_owned()
    {
        return Ok(value);
    }
    Err(Error::Unsupported("exact Frobenius coefficients require specialized Gaussian-rational constants; rounded numbers and unresolved algebraic constants cannot prove exact relations".into()))
}
pub(crate) fn gaussian_bounded(
    expression: &Atom,
    limits: &ExactFrobeniusLimits,
    context: &RunContext,
) -> Result<Gaussian> {
    let domain = ExactDomain::dense(limits, context);
    domain.preflight_series(
        expression,
        symbol!("symbolica_amflow::exact_constant_preflight"),
        0,
    )?;
    let value = gaussian(expression)?;
    domain.check_value(&value)?;
    Ok(value)
}
pub(crate) fn height(a: &Gaussian) -> u64 {
    [&a.re, &a.im]
        .into_iter()
        .map(|v| {
            v.numerator_ref()
                .significant_bits()
                .max(v.denominator_ref().significant_bits())
        })
        .max()
        .unwrap_or(0)
}
pub(crate) fn field() -> FloatField<Gaussian> {
    FloatField::from_rep(Gaussian::new(Rational::zero(), Rational::zero()))
}

pub(crate) fn preflight_expression(
    expression: &Atom,
    variable: Symbol,
    order: usize,
    limits: &ExactFrobeniusLimits,
    context: &RunContext,
) -> Result<()> {
    ExactDomain::dense(limits, context).preflight_series(expression, variable, order)
}

/// The caller proves H(0)=1 exactly. Native series arithmetic then selects the
/// unique inverse square-root germ with constant one, over the existing field.
pub(crate) fn normalized_root_series(
    h: &Atom,
    factor: &Atom,
    variable: Symbol,
    order: usize,
    limits: &ExactFrobeniusLimits,
    context: &RunContext,
) -> Result<Vec<Gaussian>> {
    let domain = ExactDomain::dense(limits, context);
    domain.preflight_series(h, variable, order)?;
    domain.preflight_series(factor, variable, order)?;
    let h0 = crate::family::substitute(h, &BTreeMap::from([(Atom::var(variable), Atom::new())]))
        .together()
        .cancel();
    if gaussian_bounded(&h0, limits, context)? != Gaussian::new(Rational::one(), Rational::zero()) {
        return Err(Error::InvalidInput(
            "normalized exact root series must have constant one".into(),
        ));
    }
    // Binomial denominators and products can grow beyond the rational input.
    let mut nodes = 0;
    let bits = coefficient_bits(
        domain
            .series_bound(h.as_view(), variable, order, &mut nodes)?
            .checked_add(domain.series_bound(factor.as_view(), variable, order, &mut nodes)?)
            .and_then(|v| v.checked_add(4))
            .and_then(|v| v.checked_mul(series_terms(order).ok()?.try_into().ok()?))
            .and_then(|v| v.checked_mul(2)),
        limits,
    )?;
    if limits
        .max_coefficient_bits
        .is_some_and(|limit| bits > limit)
    {
        return Err(Error::Limit(
            "exact root series height estimate exceeds limit".into(),
        ));
    }
    let expression = factor * h.pow(Atom::num(Rational::from((-1, 2))));
    context.cancellation.check()?;
    process_resources::checkpoint()?;
    let series = expression
        .series(variable, 0, native_series_order(order)?)
        .map_err(|e| Error::Unsupported(format!("native exact root series: {e}")))?;
    context.cancellation.check()?;
    process_resources::checkpoint()?;
    regular_root_coefficients(series.terms(), order, &domain)
}

fn regular_root_coefficients<'a>(
    terms: impl Iterator<Item = (Rational, &'a Atom)>,
    order: usize,
    domain: &ExactDomain<'_>,
) -> Result<Vec<Gaussian>> {
    let count = series_terms(order)?;
    domain.check_allocation(count)?;
    let mut result = vec![domain.zero(); count];
    for (power, coefficient) in terms {
        domain.context.cancellation.check()?;
        // Native series retain their ramification grid, including exact zero
        // slots at half-integer powers of this regular unit-root series.
        if coefficient.is_zero() {
            continue;
        }
        let k = power.to_string().parse::<usize>().map_err(|_| {
            Error::Unsupported(format!(
                "normalized root product has nonzero nonregular power {power}"
            ))
        })?;
        if k <= order {
            result[k] = domain.evaluate(coefficient)?;
        }
    }
    Ok(result)
}

pub(crate) fn rank(
    rows: &[Vec<Gaussian>],
    limits: &ExactFrobeniusLimits,
    context: &RunContext,
) -> Result<usize> {
    let domain = ExactDomain::dense(limits, context);
    domain.matrix_bound(rows)?;
    let matrix = ExactMatrix::from_nested_vec(rows.to_vec(), field()).map_err(Error::Numerical)?;
    let result = matrix.rank();
    context.cancellation.check()?;
    Ok(result)
}

pub(crate) fn inverse(
    rows: &[Vec<Gaussian>],
    limits: &ExactFrobeniusLimits,
    context: &RunContext,
) -> Result<Vec<Vec<Gaussian>>> {
    let domain = ExactDomain::dense(limits, context);
    domain.matrix_bound(rows)?;
    let n = rows.len();
    let inverse = ExactMatrix::from_nested_vec(rows.to_vec(), field())
        .map_err(Error::Numerical)?
        .inv()
        .map_err(|e| Error::Numerical(format!("exact prefix inverse: {e}")))?;
    let result = inverse
        .into_vec()
        .chunks(n)
        .map(<[Gaussian]>::to_vec)
        .collect::<Vec<_>>();
    domain.check_rows(&result)?;
    Ok(result)
}

pub(crate) fn polynomial_product(
    mut a: Vec<Gaussian>,
    mut b: Vec<Gaussian>,
    variable: Symbol,
    limits: &ExactFrobeniusLimits,
    context: &RunContext,
) -> Result<Vec<Gaussian>> {
    let domain = ExactDomain::dense(limits, context);
    domain.check_rows(std::slice::from_ref(&a))?;
    domain.check_rows(std::slice::from_ref(&b))?;
    while a.len() > 1 && a.last().is_some_and(Gaussian::is_zero) {
        a.pop();
    }
    while b.len() > 1 && b.last().is_some_and(Gaussian::is_zero) {
        b.pop();
    }
    let terms = a.len().min(b.len());
    let bits = coefficient_bits(
        a.iter()
            .map(height)
            .max()
            .unwrap_or(0)
            .checked_add(b.iter().map(height).max().unwrap_or(0))
            .and_then(|v| v.checked_add(1))
            .and_then(|v| v.checked_mul(terms.try_into().ok()?))
            .and_then(|v| v.checked_mul(2)),
        limits,
    )?;
    let cells = scalar_count(a.len().checked_mul(b.len()))?;
    if limits.max_scalar_cells.is_some_and(|limit| cells > limit)
        || limits
            .max_coefficient_bits
            .is_some_and(|limit| bits > limit)
    {
        return Err(Error::Limit(
            "exact root polynomial product estimate exceeds limit".into(),
        ));
    }
    domain.check_allocation(scalar_count(a.len().checked_add(b.len()))?)?;
    let variable = std::sync::Arc::new(PolyVariable::Symbol(variable));
    let a = symbolica::poly::univariate::UnivariatePolynomial::from_coefficients(
        &field(),
        a,
        variable.clone(),
    );
    let b =
        symbolica::poly::univariate::UnivariatePolynomial::from_coefficients(&field(), b, variable);
    let product = a * &b;
    let result = product.coefficients().to_vec();
    domain.check_rows(std::slice::from_ref(&result))?;
    Ok(result)
}

pub(crate) fn row_reduce_with_mode(
    rows: Vec<Vec<Gaussian>>,
    columns: usize,
    limits: &ExactFrobeniusLimits,
    context: &RunContext,
    elimination_mode: ExactEliminationMode,
) -> Result<Vec<Vec<Gaussian>>> {
    let augmented_columns = series_terms(columns)?;
    if rows.is_empty() || rows.iter().any(|row| row.len() != augmented_columns) {
        return Err(Error::InvalidInput(
            "exact asymptotic relation matrix dimensions".into(),
        ));
    }
    let domain = ExactDomain {
        limits,
        context,
        elimination_mode,
    };
    limits.validate()?;
    if elimination_mode == ExactEliminationMode::ObservedSparse {
        return domain.sparse_reduced_rows(rows);
    }
    domain.matrix_bound(&rows)?;
    let mut matrix = ExactMatrix::from_nested_vec(rows, field()).map_err(Error::Numerical)?;
    let columns_u32 = u32::try_from(columns)
        .map_err(|_| Error::Limit("exact relation column count exceeds u32".into()))?;
    matrix.row_reduce(columns_u32);
    context.cancellation.check()?;
    let rows = matrix
        .into_vec()
        .chunks(augmented_columns)
        .map(<[Gaussian]>::to_vec)
        .collect::<Vec<_>>();
    domain.check_rows(&rows)?;
    Ok(rows)
}

pub(crate) fn multiply(
    left: Vec<Vec<Gaussian>>,
    right: Vec<Vec<Gaussian>>,
    limits: &ExactFrobeniusLimits,
    context: &RunContext,
) -> Result<Vec<Vec<Gaussian>>> {
    context.cancellation.check()?;
    let rows = left.len();
    let inner = right.len();
    let columns = right.first().map_or(0, Vec::len);
    if rows == 0
        || inner == 0
        || columns == 0
        || left.iter().any(|row| row.len() != inner)
        || right.iter().any(|row| row.len() != columns)
    {
        return Err(Error::InvalidInput(
            "exact asymptotic product dimensions".into(),
        ));
    }
    let domain = ExactDomain::dense(limits, context);
    domain.check_allocation(scalar_count(rows.checked_mul(columns))?)?;
    let bits = coefficient_bits(
        left.iter()
            .flatten()
            .map(height)
            .max()
            .unwrap_or(0)
            .checked_add(right.iter().flatten().map(height).max().unwrap_or(0))
            .and_then(|v| v.checked_add(1))
            .and_then(|v| v.checked_mul(inner.try_into().ok()?))
            .and_then(|v| v.checked_mul(2)),
        limits,
    )?;
    if limits
        .max_coefficient_bits
        .is_some_and(|limit| bits > limit)
    {
        return Err(Error::Limit(
            "exact asymptotic product height estimate exceeds limit".into(),
        ));
    }
    let left = ExactMatrix::from_nested_vec(left, field()).map_err(Error::Numerical)?;
    let right = ExactMatrix::from_nested_vec(right, field()).map_err(Error::Numerical)?;
    let result = (&left * &right)
        .into_vec()
        .chunks(columns)
        .map(<[Gaussian]>::to_vec)
        .collect::<Vec<_>>();
    ExactDomain::dense(limits, context).check_rows(&result)?;
    Ok(result)
}

struct ExactDomain<'a> {
    limits: &'a ExactFrobeniusLimits,
    context: &'a RunContext,
    elimination_mode: ExactEliminationMode,
}
impl<'a> ExactDomain<'a> {
    fn dense(limits: &'a ExactFrobeniusLimits, context: &'a RunContext) -> Self {
        Self {
            limits,
            context,
            elimination_mode: ExactEliminationMode::ConservativeDense,
        }
    }

    fn check_allocation(&self, cells: usize) -> Result<()> {
        self.context.cancellation.check()?;
        if self
            .limits
            .max_scalar_cells
            .is_some_and(|limit| cells > limit)
        {
            return Err(Error::Limit(
                "exact Frobenius scalar allocation exceeds diagnostic cap".into(),
            ));
        }
        process_resources::check_allocation(allocation_bytes(
            cells,
            std::mem::size_of::<Gaussian>(),
        )?)?;
        Ok(())
    }

    fn series_bound(
        &self,
        view: AtomView<'_>,
        variable: Symbol,
        order: usize,
        nodes: &mut usize,
    ) -> Result<u64> {
        crate::stack_resources::run(|| self.expression(view, variable, order, nodes))
    }

    fn expression(
        &self,
        view: AtomView<'_>,
        variable: Symbol,
        order: usize,
        nodes: &mut usize,
    ) -> Result<u64> {
        self.context.cancellation.check()?;
        *nodes = scalar_count(nodes.checked_add(1))?;
        if self
            .limits
            .max_scalar_cells
            .is_some_and(|limit| *nodes > limit)
        {
            return Err(Error::Limit(
                "exact Frobenius expression size exceeds diagnostic cap".into(),
            ));
        }
        let bound = match view {
            AtomView::Num(number) => match number.get_coeff_view().to_owned() {
                symbolica::coefficient::Coefficient::Complex(value) => height(&value),
                _ => {
                    return Err(Error::Unsupported(
                        "exact Frobenius source contains a nonexact coefficient".into(),
                    ));
                }
            },
            AtomView::Var(v)
                if v.get_symbol() == variable
                    || v.get_symbol() == crate::family::imaginary_parameter() =>
            {
                1
            }
            AtomView::Add(v) => v.iter().try_fold(0_u64, |sum, a| {
                let bound = self.series_bound(a, variable, order, nodes)?;
                coefficient_bits(
                    sum.checked_add(bound).and_then(|v| v.checked_add(1)),
                    self.limits,
                )
            })?,
            AtomView::Mul(v) => v.iter().try_fold(0_u64, |sum, a| {
                let bound = self.series_bound(a, variable, order, nodes)?;
                coefficient_bits(
                    sum.checked_add(bound).and_then(|v| v.checked_add(2)),
                    self.limits,
                )
            })?,
            AtomView::Pow(v) => {
                let (base, power) = v.get_base_exp();
                let AtomView::Num(number) = power else {
                    return Err(Error::Unsupported(
                        "exact Frobenius source has a noninteger power".into(),
                    ));
                };
                let symbolica::coefficient::Coefficient::Complex(power) =
                    number.get_coeff_view().to_owned()
                else {
                    return Err(Error::Unsupported(
                        "exact Frobenius source has a nonexact power".into(),
                    ));
                };
                if !power.im.is_zero() || !power.re.is_integer() {
                    return Err(Error::Unsupported(
                        "exact Frobenius source has a noninteger power".into(),
                    ));
                }
                let negative = power.re < Rational::zero();
                let magnitude = power.re.abs().to_string().parse::<u64>().ok();
                // Only native nonnegative integer series powers convert to i64;
                // negative rational powers retain arbitrary exact numerators.
                if !negative && magnitude.is_none_or(|n| i64::try_from(n).is_err()) {
                    return Err(Error::Limit(
                        "exact Frobenius native integer power exceeds i64 representation".into(),
                    ));
                }
                coefficient_bits(
                    self.series_bound(base, variable, order, nodes)?
                        .checked_add(1)
                        .and_then(|v| v.checked_mul(magnitude?))
                        .and_then(|v| {
                            v.checked_mul(if negative {
                                series_terms(order).ok()?.try_into().ok()?
                            } else {
                                1
                            })
                        }),
                    self.limits,
                )?
            }
            _ => {
                return Err(Error::Unsupported(
                    "exact Frobenius source requires a specialized rational coefficient field"
                        .into(),
                ));
            }
        };
        if self
            .limits
            .max_coefficient_bits
            .is_some_and(|limit| bound > limit)
        {
            return Err(Error::Limit(
                "exact Frobenius series coefficient-height estimate exceeds limit".into(),
            ));
        }
        Ok(bound)
    }
    fn matrix_bound(&self, a: &[Vec<Gaussian>]) -> Result<()> {
        self.context.cancellation.check()?;
        self.preflight_matrix(a.len().max(a.first().map_or(0, Vec::len)))?;
        self.check_rows(a)?;
        let n = u64::try_from(a.len().max(1))
            .map_err(|_| Error::Limit("exact Frobenius matrix dimension exceeds u64".into()))?;
        let bits = a.iter().flatten().map(height).max().unwrap_or(0);
        // Common denominator clearing plus determinant/minor growth. This is
        // deliberately conservative before a native elimination/inverse call.
        let predicted = coefficient_bits(
            bits.checked_add(1)
                .and_then(|v| v.checked_mul(n))
                .and_then(|v| v.checked_mul(n))
                .and_then(|v| v.checked_mul(4)),
            self.limits,
        )?;
        if self
            .limits
            .max_coefficient_bits
            .is_some_and(|limit| predicted > limit)
        {
            return Err(Error::Limit(
                "exact Frobenius elimination height estimate exceeds limit".into(),
            ));
        }
        Ok(())
    }
}
impl CoefficientDomain for ExactDomain<'_> {
    type Scalar = Gaussian;
    type Solver = ExactMatrix;
    fn zero(&self) -> Gaussian {
        Gaussian::new(Rational::zero(), Rational::zero())
    }
    fn integer(&self, value: i64) -> Gaussian {
        Gaussian::new(Rational::from(value), Rational::zero())
    }
    fn add(&self, a: &Gaussian, b: &Gaussian) -> Gaussian {
        a + b
    }
    fn sub(&self, a: &Gaussian, b: &Gaussian) -> Gaussian {
        a - b
    }
    fn mul(&self, a: &Gaussian, b: &Gaussian) -> Gaussian {
        a * b
    }
    fn neg(&self, a: &Gaussian) -> Gaussian {
        -a.clone()
    }
    fn scale(&self, a: &Gaussian, n: i64, d: i64) -> Gaussian {
        let scale = Rational::from((n, d));
        Gaussian::new(&a.re * &scale, &a.im * &scale)
    }
    fn evaluate(&self, expression: &Atom) -> Result<Gaussian> {
        self.context.cancellation.check()?;
        let value = gaussian_bounded(expression, self.limits, self.context)?;
        if self
            .limits
            .max_coefficient_bits
            .is_some_and(|limit| height(&value) > limit)
        {
            return Err(Error::Limit(
                "exact Frobenius coefficient height exceeds limit".into(),
            ));
        }
        Ok(value)
    }
    fn series(&self, expression: &Atom, variable: Symbol, order: usize) -> Result<Vec<Gaussian>> {
        self.preflight_series(expression, variable, order)?;
        let series = expression
            .series(variable, 0, native_series_order(order)?)
            .map_err(|e| Error::Unsupported(e.to_string()))?;
        process_resources::checkpoint()?;
        let mut coefficients = vec![self.zero(); series_terms(order)?];
        for (power, coefficient) in series.terms() {
            let power = power.to_string().parse::<usize>().map_err(|_| {
                Error::Unsupported("exact normalized Frobenius coefficient is not regular".into())
            })?;
            if power <= order {
                coefficients[power] = self.evaluate(coefficient)?;
            }
        }
        Ok(coefficients)
    }
    fn near_resonance(&self, _difference: &Atom) -> Result<bool> {
        Ok(false)
    }
    fn insignificant(&self, value: &Gaussian) -> bool {
        value.is_zero()
    }
    fn solver(
        &self,
        matrix: Vec<Vec<Gaussian>>,
        blocks: &[Vec<usize>],
        context: &RunContext,
    ) -> Result<ExactMatrix> {
        if self.elimination_mode == ExactEliminationMode::ObservedSparse {
            self.validate_sparse_blocks(&matrix, blocks)?;
            return self.sparse_inverse(matrix, context);
        }
        self.matrix_bound(&matrix)?;
        let matrix = ExactMatrix::from_nested_vec(matrix, field()).map_err(Error::Numerical)?;
        let inverse = matrix
            .inv()
            .map_err(|e| Error::Numerical(format!("exact Frobenius recurrence inverse: {e}")))?;
        context.cancellation.check()?;
        let rows = inverse
            .clone()
            .into_vec()
            .chunks(inverse.ncols())
            .map(<[Gaussian]>::to_vec)
            .collect::<Vec<_>>();
        self.check_rows(&rows)?;
        Ok(inverse)
    }
    fn solve(&self, solver: &ExactMatrix, rhs: &[Gaussian]) -> Result<Vec<Gaussian>> {
        let rows = solver
            .clone()
            .into_vec()
            .chunks(solver.ncols())
            .map(<[Gaussian]>::to_vec)
            .collect();
        Ok(multiply(
            rows,
            rhs.iter().cloned().map(|a| vec![a]).collect(),
            self.limits,
            self.context,
        )?
        .into_iter()
        .flatten()
        .collect())
    }
    fn solve_singular(
        &self,
        matrix: Vec<Vec<Gaussian>>,
        rhs: Vec<Gaussian>,
        context: &RunContext,
    ) -> Result<Vec<Gaussian>> {
        if self.elimination_mode == ExactEliminationMode::ObservedSparse {
            return self.sparse_solve_any(matrix, rhs, context);
        }
        let augmented = matrix
            .iter()
            .zip(&rhs)
            .map(|(row, b)| {
                let mut row = row.clone();
                row.push(b.clone());
                row
            })
            .collect::<Vec<_>>();
        self.matrix_bound(&augmented)?;
        let matrix = ExactMatrix::from_nested_vec(matrix, field()).map_err(Error::Numerical)?;
        let rhs = ExactMatrix::from_nested_vec(rhs.into_iter().map(|a| vec![a]).collect(), field())
            .map_err(Error::Numerical)?;
        let result = matrix
            .solve_any(&rhs)
            .map_err(|e| Error::Numerical(format!("exact resonant Frobenius recurrence: {e}")))?
            .into_vec();
        context.cancellation.check()?;
        self.check_rows(std::slice::from_ref(&result))?;
        Ok(result)
    }
    fn preflight(&self, n: usize, order: usize) -> Result<()> {
        self.context.cancellation.check()?;
        process_resources::checkpoint()?;
        self.limits.validate()?;
        if n == 0
            || self.limits.max_dimension.is_some_and(|limit| n > limit)
            || self.limits.max_order.is_some_and(|limit| order > limit)
        {
            return Err(Error::Limit(
                "exact Frobenius dimension or order exceeds limit".into(),
            ));
        }
        let cells = scalar_count(
            n.checked_mul(n)
                .and_then(|v| v.checked_mul(n.checked_add(1)?))
                .and_then(|v| v.checked_mul(order.checked_add(1)?)),
        )?;
        self.check_allocation(cells)
    }
    fn preflight_matrix(&self, n: usize) -> Result<()> {
        self.context.cancellation.check()?;
        self.check_allocation(scalar_count(n.checked_mul(n))?)
    }
    fn preflight_series(&self, expression: &Atom, variable: Symbol, order: usize) -> Result<()> {
        self.context.cancellation.check()?;
        self.limits.validate()?;
        if self.limits.max_order.is_some_and(|limit| order > limit) {
            return Err(Error::Limit(
                "exact Frobenius series order exceeds diagnostic cap".into(),
            ));
        }
        native_series_order(order)?;
        self.check_allocation(series_terms(order)?)?;
        process_resources::check_allocation(allocation_bytes(
            expression.as_view().get_byte_size(),
            1,
        )?)?;
        let mut nodes = 0;
        self.series_bound(expression.as_view(), variable, order, &mut nodes)?;
        process_resources::checkpoint()?;
        Ok(())
    }
    fn check_rows(&self, rows: &[Vec<Gaussian>]) -> Result<()> {
        self.context.cancellation.check()?;
        process_resources::checkpoint()?;
        if self
            .limits
            .max_coefficient_bits
            .is_some_and(|limit| rows.iter().flatten().any(|a| height(a) > limit))
        {
            return Err(Error::Limit(
                "exact Frobenius observed coefficient height exceeds limit".into(),
            ));
        }
        Ok(())
    }
    fn check_value(&self, value: &Gaussian) -> Result<()> {
        self.context.cancellation.check()?;
        process_resources::checkpoint()?;
        if self
            .limits
            .max_coefficient_bits
            .is_some_and(|limit| height(value) > limit)
        {
            return Err(Error::Limit(
                "exact Frobenius accumulation height exceeds limit".into(),
            ));
        }
        Ok(())
    }
}

impl PreparedFrobenius {
    /// Compute exact coefficients using the same recurrence and normalization
    /// as numerical Frobenius evaluation. Unresolved coefficient fields are
    /// rejected; numerical boundary values never enter this construction.
    pub fn exact_coefficients(
        &self,
        order: usize,
        limits: &ExactFrobeniusLimits,
        context: &RunContext,
    ) -> Result<ExactFrobeniusBasis> {
        self.exact_coefficients_with_elimination(
            order,
            limits,
            ExactEliminationMode::ConservativeDense,
            context,
        )
    }

    /// Compute an exact prefix with an explicitly selected native elimination policy.
    /// All existing dimension, series, allocation and observed coefficient limits
    /// remain active. Sparse mode substitutes operation-specific observed height
    /// checks for the conservative dense elimination estimate only: inverse inputs
    /// and outputs, plus native row boundaries for resonances and constraints.
    /// Internal inverse/elimination arithmetic remains unobserved. No numerical
    /// amplitudes or precision thresholds enter this operation.
    pub fn exact_coefficients_with_elimination(
        &self,
        order: usize,
        limits: &ExactFrobeniusLimits,
        elimination_mode: ExactEliminationMode,
        context: &RunContext,
    ) -> Result<ExactFrobeniusBasis> {
        let domain = ExactDomain {
            limits,
            context,
            elimination_mode,
        };
        domain.preflight(self.source.matrix.len(), order)?;
        for exponent in self.exponents() {
            let value = gaussian_bounded(exponent, limits, context)?;
            if !value.im.is_zero() {
                return Err(Error::Unsupported(
                    "exact endpoint coefficient support requires rational real exponents".into(),
                ));
            }
        }
        let normal = self.recurrence_normal(&domain, order, context)?;
        let columns = self.transform_recurrence(&domain, normal, order, context)?;
        Ok(ExactFrobeniusBasis {
            columns,
            known_order: order,
            limits: limits.clone(),
            elimination_mode,
        })
    }
}
