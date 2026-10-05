//! Exact diagonal epsilon gauges and lossless Laurent coefficient reindexing.
use crate::diffexp::EpsilonBoundary;
use crate::kinematics::KinematicSystem;
use crate::{ComplexFloat, DifferentialSystem, Error, Result, RunContext};
use symbolica::domains::algebraic::AlgebraicExtension;
use symbolica::domains::rational::RationalField;
use symbolica::prelude::*;

type CoefficientErrors = Vec<Vec<Float>>;
type ShiftedBoundary = (EpsilonBoundary, CoefficientErrors);

/// Original I = diag(epsilon^weights) J. Weights are normalized to max=0.
/// The source's declared common Laurent pole bound is preserved; individual
/// components may acquire additional proved-zero leading coefficients.
/// That bound must hold throughout the intended transport, not just at a
/// starting point where some coefficients happen to vanish.
#[derive(Clone, Debug)]
pub struct EpsilonShearing {
    epsilon: Symbol,
    weights: Vec<i32>,
    nonzero_conditions: Vec<Atom>,
}

fn valuation(entry: &Atom, epsilon: Symbol) -> Result<Option<i64>> {
    let gaussian = AlgebraicExtension::complex(Q);
    let rational: RationalPolynomial<AlgebraicExtension<RationalField>, u16> = entry
        .try_to_rational_polynomial(&gaussian, &gaussian, None)
        .map_err(|e| Error::Unsupported(format!("rational epsilon shearing required: {e}")))?;
    if rational.numerator.is_zero() {
        return Ok(None);
    }
    if rational
        .numerator
        .variables()
        .iter()
        .any(|v| !matches!(v, PolyVariable::Symbol(_)))
    {
        return Err(Error::Unsupported(
            "epsilon shearing requires rational symbol variables".into(),
        ));
    }
    let order = |poly: &MultivariatePolynomial<AlgebraicExtension<RationalField>, u16>| {
        poly.variables()
            .iter()
            .position(|v| *v == PolyVariable::Symbol(epsilon))
            .map_or(0, |index| i64::from(poly.degree_bounds(index).0))
    };
    Ok(Some(
        order(&rational.numerator) - order(&rational.denominator),
    ))
}

fn range_count(leading: i32, last: i32) -> Result<usize> {
    let count = i64::from(last) - i64::from(leading) + 1;
    if !(1..=1024).contains(&count) {
        return Err(Error::Limit(
            "epsilon shearing requires 1 through 1024 Laurent coefficients".into(),
        ));
    }
    Ok(count as usize)
}

impl EpsilonShearing {
    /// Regularize all physical partials with one exact epsilon-only diagonal
    /// transformation. A negative cycle rejects this restricted gauge class;
    /// it does not assert that every nondiagonal gauge is impossible.
    pub fn regularize(
        system: &KinematicSystem,
        context: &RunContext,
    ) -> Result<(KinematicSystem, Self)> {
        system.validate()?;
        context.cancellation.check()?;
        let n = system.derivatives.first_key_value().unwrap().1.len();
        let variables = system
            .derivatives
            .keys()
            .copied()
            .chain([system.epsilon])
            .collect();
        let entries = system
            .derivatives
            .values()
            .flatten()
            .flatten()
            .cloned()
            .collect::<Vec<_>>();
        // Capture raw poles before cancellation. The existing owner also keeps
        // leading epsilon denominators before physical path specialization.
        let mut conditions =
            crate::physical_conditions::rational_denominator_conditions(&entries, &variables)?;
        conditions.extend(crate::physical_conditions::matrix_domain_conditions(
            system,
        )?);
        let conditions = crate::physical_conditions::canonical_conditions(&conditions, &variables)?;
        let mut orders = vec![vec![None; n]; n];
        for matrix in system.derivatives.values() {
            for (i, row) in matrix.iter().enumerate() {
                context.cancellation.check()?;
                for (j, entry) in row.iter().enumerate() {
                    if let Some(order) = valuation(entry, system.epsilon)? {
                        orders[i][j] =
                            Some(orders[i][j].map_or(order, |previous: i64| previous.min(order)));
                    }
                }
            }
        }
        let result = Self::from_orders(system.epsilon, &orders, conditions, context)?;
        let transform = result.transformation();
        let mut derivatives = std::collections::BTreeMap::new();
        for (&variable, matrix) in &system.derivatives {
            context.cancellation.check()?;
            let transformed = DifferentialSystem {
                variable,
                matrix: matrix.clone(),
            }
            .change_basis(&transform)?;
            for entry in transformed.matrix.iter().flatten() {
                if valuation(entry, system.epsilon)?.is_some_and(|order| order < 0) {
                    return Err(Error::Numerical(
                        "epsilon shearing leaves a negative matrix valuation".into(),
                    ));
                }
            }
            derivatives.insert(variable, transformed.matrix);
        }
        Ok((
            KinematicSystem {
                epsilon: system.epsilon,
                derivatives,
            },
            result,
        ))
    }

    fn from_orders(
        epsilon: Symbol,
        orders: &[Vec<Option<i64>>],
        conditions: Vec<Atom>,
        context: &RunContext,
    ) -> Result<Self> {
        let weights = crate::integer_shearing::integer_potentials(orders, 0, 10000, Some(context))?
            .into_iter()
            .map(|w| {
                i32::try_from(w)
                    .map_err(|_| Error::Limit("epsilon shearing weight overflow".into()))
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            epsilon,
            weights,
            nonzero_conditions: conditions,
        })
    }

    /// Regularize a formal registered-square-root connection. The root relations
    /// are normalized before epsilon valuations or series expansion; their
    /// epsilon-independent radicands and independent named sheets are preserved.
    /// The returned descriptor retains source domains erased by normalization.
    pub fn regularize_algebraic(
        system: &crate::algebraic::AlgebraicKinematicSystem,
        context: &RunContext,
    ) -> Result<(crate::algebraic::AlgebraicKinematicSystem, Self)> {
        system.validate()?;
        context.cancellation.check()?;
        let n = system.derivatives.first_key_value().unwrap().1.len();
        let variables = system
            .derivatives
            .keys()
            .copied()
            .chain([system.epsilon])
            .collect::<std::collections::BTreeSet<_>>();
        let allowed = variables
            .iter()
            .copied()
            .chain([crate::family::imaginary_parameter()])
            .chain(system.roots.iter().map(|r| r.symbol))
            .map(Atom::var)
            .collect();
        // Capture original holes and formal norm domains before multiplying
        // numerators or applying the epsilon-only gauge.
        let mut conditions = system.nonzero_conditions()?;
        let mut normalized = system.clone();
        let mut orders = vec![vec![None; n]; n];
        let mut inverse_coefficients = Vec::new();
        for matrix in normalized.derivatives.values_mut() {
            for (i, row) in matrix.iter_mut().enumerate() {
                context.cancellation.check()?;
                for (j, entry) in row.iter_mut().enumerate() {
                    let native =
                        crate::algebraic::normalized_root_entry(entry, &system.roots, &allowed)?;
                    for coefficient in native.coefficients() {
                        if let Some(order) = valuation(&coefficient, system.epsilon)? {
                            orders[i][j] =
                                Some(orders[i][j].map_or(order, |old: i64| old.min(order)));
                        }
                    }
                    *entry = native.expression(&system.roots);
                    inverse_coefficients.extend(native.decoded_inverse_coefficients());
                }
            }
        }
        conditions.extend(crate::physical_conditions::rational_denominator_conditions(
            &inverse_coefficients,
            &variables,
        )?);
        conditions = crate::physical_conditions::canonical_conditions(&conditions, &variables)?;
        let gaussian = AlgebraicExtension::complex(Q);
        let mut generic = Vec::new();
        for condition in &conditions {
            let rational: RationalPolynomial<AlgebraicExtension<RationalField>, u16> = condition
                .try_to_rational_polynomial(&gaussian, &gaussian, None)
                .map_err(|e| {
                    Error::Unsupported(format!("rational root-norm domain required: {e}"))
                })?;
            for part in [
                rational.numerator.to_expression(),
                rational.denominator.to_expression(),
            ] {
                generic.push(crate::physical_conditions::epsilon_leading_coefficient(
                    &part,
                    system.epsilon,
                )?);
            }
        }
        conditions.extend(generic);
        conditions = crate::physical_conditions::canonical_conditions(&conditions, &variables)?;
        let result = Self::from_orders(system.epsilon, &orders, conditions, context)?;
        let transform = result.transformation();
        for (&variable, matrix) in &mut normalized.derivatives {
            context.cancellation.check()?;
            *matrix = DifferentialSystem {
                variable,
                matrix: std::mem::take(matrix),
            }
            .change_basis(&transform)?
            .matrix;
            for entry in matrix.iter_mut().flatten() {
                let native =
                    crate::algebraic::normalized_root_entry(entry, &system.roots, &allowed)?;
                for coefficient in native.coefficients() {
                    if valuation(&coefficient, system.epsilon)?.is_some_and(|order| order < 0) {
                        return Err(Error::Numerical(
                            "algebraic epsilon shearing leaves a negative coefficient valuation"
                                .into(),
                        ));
                    }
                }
                *entry = native.expression(&system.roots);
            }
        }
        normalized.validate()?;
        Ok((normalized, result))
    }

    pub fn epsilon(&self) -> Symbol {
        self.epsilon
    }
    pub fn weights(&self) -> &[i32] {
        &self.weights
    }
    /// Original raw denominator and generic epsilon-leading domain conditions.
    pub fn nonzero_conditions(&self) -> &[Atom] {
        &self.nonzero_conditions
    }
    /// Exact T in original I=T J. Physical derivatives of T vanish.
    pub fn transformation(&self) -> Vec<Vec<Atom>> {
        self.weights
            .iter()
            .enumerate()
            .map(|(i, &weight)| {
                (0..self.weights.len())
                    .map(|j| {
                        if i == j {
                            Atom::var(self.epsilon).pow(i64::from(weight))
                        } else {
                            Atom::new()
                        }
                    })
                    .collect()
            })
            .collect()
    }
    /// Cache labels explicitly describe J_i=epsilon^-w_i I_i, even for a zero
    /// connection. These must not replace labels on an unconverted old entry.
    pub fn scaled_basis_labels(&self, labels: &[Atom]) -> Result<Vec<Atom>> {
        if labels.len() != self.weights.len() {
            return Err(Error::InvalidInput(
                "epsilon shearing basis dimensions".into(),
            ));
        }
        Ok(labels
            .iter()
            .zip(&self.weights)
            .map(|(a, &w)| a * Atom::var(self.epsilon).pow(-i64::from(w)))
            .collect())
    }
    /// Required rectangular J range to restore every original component through
    /// `last`. The original common pole bound must already be justified.
    pub fn required_sheared_range(&self, leading: i32, last: i32) -> Result<(i32, i32)> {
        range_count(leading, last)?;
        let last = self.weights.iter().try_fold(last, |high, &w| {
            let next = last
                .checked_sub(w)
                .ok_or_else(|| Error::Limit("epsilon shearing Laurent order overflow".into()))?;
            Ok::<_, Error>(high.max(next))
        })?;
        range_count(leading, last)?;
        Ok((leading, last))
    }

    fn validate_boundary(&self, boundary: &EpsilonBoundary, errors: &[Vec<Float>]) -> Result<i32> {
        let n = self.weights.len();
        let source_count = i64::try_from(boundary.coefficients.len())
            .map_err(|_| Error::Limit("epsilon boundary coefficient count overflow".into()))?;
        let last = i64::from(boundary.leading) + source_count - 1;
        let last = i32::try_from(last)
            .map_err(|_| Error::Limit("epsilon boundary order overflow".into()))?;
        range_count(boundary.leading, last)?;
        if !boundary.point.re.is_finite()
            || !boundary.point.im.is_finite()
            || errors.len() != boundary.coefficients.len()
            || boundary.coefficients.iter().any(|row| row.len() != n)
            || errors.iter().any(|row| row.len() != n)
            || boundary
                .coefficients
                .iter()
                .flatten()
                .any(|v| !v.re.is_finite() || !v.im.is_finite())
            || errors
                .iter()
                .flatten()
                .any(|e| !e.is_finite() || *e < Float::with_val(32, 0))
        {
            return Err(Error::InvalidInput(
                "epsilon shearing boundary/error dimensions or values".into(),
            ));
        }
        Ok(last)
    }

    /// Reindex I coefficients into J, keeping the source's justified common
    /// leading bound. Only coefficients below that original bound are zero.
    /// Unknown positive orders are rejected; errors shift without rescaling.
    pub fn to_sheared(
        &self,
        boundary: &EpsilonBoundary,
        errors: &[Vec<Float>],
        last: i32,
    ) -> Result<ShiftedBoundary> {
        let source_last = self.validate_boundary(boundary, errors)?;
        let count = range_count(boundary.leading, last)?;
        let mut coefficients = Vec::with_capacity(count);
        let mut shifted_errors = Vec::with_capacity(count);
        for offset in 0..count {
            let power = i64::from(boundary.leading) + offset as i64;
            let mut row = Vec::with_capacity(self.weights.len());
            let mut error_row = Vec::with_capacity(self.weights.len());
            for (column, &weight) in self.weights.iter().enumerate() {
                let original = power + i64::from(weight);
                if original > i64::from(source_last) {
                    return Err(Error::InvalidInput(format!(
                        "sheared boundary requires original component {column} through epsilon^{original}, but source ends at {source_last}"
                    )));
                }
                if original < i64::from(boundary.leading) {
                    let value = &boundary.coefficients[0][column];
                    row.push(ComplexFloat::new(
                        Float::with_val(value.re.as_raw().prec(), 0),
                        Float::with_val(value.im.as_raw().prec(), 0),
                    ));
                    error_row.push(Float::with_val(errors[0][column].as_raw().prec(), 0));
                } else {
                    let index = (original - i64::from(boundary.leading)) as usize;
                    row.push(boundary.coefficients[index][column].clone());
                    error_row.push(errors[index][column].clone());
                }
            }
            coefficients.push(row);
            shifted_errors.push(error_row);
        }
        Ok((
            EpsilonBoundary {
                point: boundary.point.clone(),
                leading: boundary.leading,
                coefficients,
            },
            shifted_errors,
        ))
    }

    /// Restore I from J on an explicitly requested range. Unlike known leading
    /// zeros in to_sheared, absent source coefficients are always an error.
    pub fn to_original(
        &self,
        boundary: &EpsilonBoundary,
        errors: &[Vec<Float>],
        leading: i32,
        last: i32,
    ) -> Result<ShiftedBoundary> {
        let source_last = self.validate_boundary(boundary, errors)?;
        let count = range_count(leading, last)?;
        let mut coefficients = Vec::with_capacity(count);
        let mut shifted_errors = Vec::with_capacity(count);
        for offset in 0..count {
            let power = i64::from(leading) + offset as i64;
            let mut row = Vec::with_capacity(self.weights.len());
            let mut error_row = Vec::with_capacity(self.weights.len());
            for (column, &weight) in self.weights.iter().enumerate() {
                let source = power - i64::from(weight);
                if source < i64::from(boundary.leading) || source > i64::from(source_last) {
                    return Err(Error::InvalidInput(format!(
                        "original boundary requires sheared component {column} at epsilon^{source}, outside available {}..={source_last}",
                        boundary.leading
                    )));
                }
                let index = (source - i64::from(boundary.leading)) as usize;
                row.push(boundary.coefficients[index][column].clone());
                error_row.push(errors[index][column].clone());
            }
            coefficients.push(row);
            shifted_errors.push(error_row);
        }
        Ok((
            EpsilonBoundary {
                point: boundary.point.clone(),
                leading,
                coefficients,
            },
            shifted_errors,
        ))
    }
}
