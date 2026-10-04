//! Physical invariant derivatives with loop integration variables held fixed.
//!
//! An external Gram matrix `G(x)` determines the covariant vector field
//! `dp_i/dx = sum_j A_ij p_j`, where `A = (dG/dx) G^-1 / 2`.
//! Thus `A G + G A^T = dG/dx`. Explicit propagator coefficients are
//! differentiated too. This is distinct from the fixed-scalar-product
//! derivative used for Feynman parameters in [`crate::reduction`].
//!
//! Construct the derivative before substituting numerical kinematics. The
//! loop momenta and integration measure stay fixed, so there is no Jacobian.

use crate::family::{imaginary_parameter, substitute};
use crate::reduction::{LinearCombination, ReducedSystem, ReductionBackend};
use crate::{Error, Integral, IntegralFamily, Propagator, Result, RunContext};
use std::collections::BTreeMap;
use symbolica::domains::float::Complex;
use symbolica::prelude::*;

/// Reusable exact physical derivative for one family and one invariant.
///
/// A changing external Gram matrix must be invertible. A constant Gram matrix
/// needs no inverse, including for mass differentiation at lightlike momenta.
/// Every returned identity is subject to [`Self::nonzero_conditions`]; the
/// differential-system builder retains these conditions in its result.
#[derive(Clone, Debug)]
pub struct KinematicDerivative {
    family: IntegralFamily,
    variable: Symbol,
    vector_field: Vec<Vec<Atom>>,
    denominator_derivatives: Vec<Propagator>,
    coefficients: Vec<Vec<Atom>>,
    constants: Vec<Atom>,
    nonzero_conditions: Vec<Atom>,
}

impl KinematicDerivative {
    pub fn new(family: &IntegralFamily, variable: Symbol) -> Result<Self> {
        if variable == family.epsilon {
            return Err(Error::Unsupported(
                "dimensional-regulator differentiation is not a physical invariant derivative"
                    .into(),
            ));
        }
        let converted = family.convert()?;
        // Preserve original coefficient-denominator and denominator-map guards,
        // including conditions that disappear from the final derivative.
        let imaginary = BTreeMap::from([(
            Atom::var(imaginary_parameter()),
            Atom::num(Complex::new(Rational::from(0), Rational::from(1))),
        )]);
        let mut conditions = converted
            .family
            .domain()
            .conditions()
            .map(|condition| {
                substitute(
                    &substitute(&condition.polynomial().to_expression(), &converted.reverse),
                    &imaginary,
                )
                .together()
                .cancel()
            })
            .collect::<Vec<_>>();
        let external = family.external.len();
        let gram_derivative = family
            .external_gram
            .iter()
            .map(|row| {
                row.iter()
                    .map(|entry| entry.derivative(variable).together().cancel())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut vector_field = vec![vec![Atom::new(); external]; external];
        if gram_derivative
            .iter()
            .flatten()
            .any(|entry| !entry.is_zero())
        {
            let determinant = crate::algebra::determinant(family.external_gram.clone())
                .together()
                .cancel();
            if determinant.is_zero() {
                return Err(Error::Unsupported(
                    "a varying singular external Gram matrix has no inverse-Gram derivative chart; choose independent external momenta or differentiate before taking the singular limit".into(),
                ));
            }
            conditions.push(determinant);
            let inverse = crate::algebra::inverse(&family.external_gram)?;
            vector_field = crate::algebra::matmul(&gram_derivative, &inverse);
            for entry in vector_field.iter_mut().flatten() {
                *entry = (&*entry / Atom::num(2)).together().cancel();
            }
            for (i, row) in gram_derivative.iter().enumerate() {
                for (j, expected) in row.iter().enumerate() {
                    let actual = (0..external).fold(Atom::new(), |sum, k| {
                        sum + &vector_field[i][k] * &family.external_gram[k][j]
                            + &family.external_gram[i][k] * &vector_field[j][k]
                    });
                    if !(actual - expected).together().cancel().is_zero() {
                        return Err(Error::Numerical(
                            "exact external-momentum vector field failed Gram covariance".into(),
                        ));
                    }
                }
            }
        }
        let offset = family.loops.len() * (family.loops.len() + 1) / 2;
        let denominator_derivatives = family
            .propagators
            .iter()
            .map(|propagator| {
                let mut derivative = Propagator {
                    constant: propagator.constant.derivative(variable).together().cancel(),
                    scalar_products: propagator
                        .scalar_products
                        .iter()
                        .map(|coefficient| coefficient.derivative(variable))
                        .collect(),
                };
                for loop_index in 0..family.loops.len() {
                    for j in 0..external {
                        let coefficient =
                            &mut derivative.scalar_products[offset + loop_index * external + j];
                        for (k, row) in vector_field.iter().enumerate() {
                            *coefficient += &propagator.scalar_products
                                [offset + loop_index * external + k]
                                * &row[j];
                        }
                    }
                }
                for coefficient in &mut derivative.scalar_products {
                    *coefficient = coefficient.together().cancel();
                }
                derivative
            })
            .collect::<Vec<_>>();
        let inverse = crate::algebra::inverse(
            &family
                .propagators
                .iter()
                .map(|propagator| propagator.scalar_products.clone())
                .collect::<Vec<_>>(),
        )?;
        let mut coefficients = Vec::with_capacity(family.propagators.len());
        let mut constants = Vec::with_capacity(family.propagators.len());
        for derivative in &denominator_derivatives {
            let row = (0..family.propagators.len())
                .map(|k| {
                    (0..inverse.len())
                        .fold(Atom::new(), |sum, j| {
                            sum + &derivative.scalar_products[j] * &inverse[j][k]
                        })
                        .together()
                        .cancel()
                })
                .collect::<Vec<_>>();
            constants.push(
                (&derivative.constant
                    - row
                        .iter()
                        .zip(&family.propagators)
                        .fold(Atom::new(), |sum, (coefficient, propagator)| {
                            sum + coefficient * &propagator.constant
                        }))
                .together()
                .cancel(),
            );
            coefficients.push(row);
        }
        if conditions.iter().any(|condition| condition.is_zero()) {
            return Err(Error::Unsupported(
                "physical derivative specializes a required family condition to zero".into(),
            ));
        }
        conditions.retain(|condition| !matches!(condition.as_view(), AtomView::Num(_)));
        conditions.sort();
        conditions.dedup();
        Ok(Self {
            family: family.clone(),
            variable,
            vector_field,
            denominator_derivatives,
            coefficients,
            constants,
            nonzero_conditions: conditions,
        })
    }

    pub fn variable(&self) -> Symbol {
        self.variable
    }

    /// Coefficients in `dp_i/dx = sum_j A[i][j] p_j`.
    pub fn external_vector_field(&self) -> &[Vec<Atom>] {
        &self.vector_field
    }

    /// Total physical derivatives of the affine propagators in scalar-product
    /// coordinates, including their explicit coefficient dependence.
    pub fn denominator_derivatives(&self) -> &[Propagator] {
        &self.denominator_derivatives
    }

    pub fn nonzero_conditions(&self) -> &[Atom] {
        &self.nonzero_conditions
    }

    /// Differentiate an integral with the exact external-momentum vector field.
    /// The identity holds on the domain returned by `nonzero_conditions()`.
    pub fn integral(&self, integral: &Integral) -> Result<LinearCombination> {
        self.family.validate_integral(integral)?;
        let mut terms = LinearCombination::new();
        for (i, &power) in integral.0.iter().enumerate() {
            if power == 0 {
                continue;
            }
            for (j, coefficient) in self.coefficients[i]
                .iter()
                .chain(std::iter::once(&self.constants[i]))
                .enumerate()
            {
                if coefficient.is_zero() {
                    continue;
                }
                let mut shifted = integral.clone();
                // Combine the shifts before checking overflow when the raised
                // and lowered denominator coincide.
                if i != j {
                    shifted.0[i] = shifted.0[i]
                        .checked_add(1)
                        .ok_or_else(|| Error::Limit("physical derivative index overflow".into()))?;
                    if j < integral.0.len() {
                        shifted.0[j] = shifted.0[j].checked_sub(1).ok_or_else(|| {
                            Error::Limit("physical derivative index overflow".into())
                        })?;
                    }
                }
                let old = terms.remove(&shifted).unwrap_or_default();
                terms.insert(
                    shifted,
                    (old - Atom::num(i64::from(power)) * coefficient)
                        .together()
                        .cancel(),
                );
            }
        }
        terms.retain(|_, coefficient| !coefficient.is_zero());
        Ok(terms)
    }

    /// Reduce physical derivatives until the same exact closure criterion used
    /// for auxiliary-mass systems is met. This does not certify minimality.
    pub fn differential_system(
        &self,
        backend: &dyn ReductionBackend,
        targets: &[Integral],
        max_rounds: usize,
        context: &RunContext,
    ) -> Result<ReducedSystem> {
        let mut system = crate::reduction::build_differential_system(
            backend,
            &self.family,
            targets,
            max_rounds,
            context,
            |integral| Ok(self.integral(integral)?.into_iter().collect()),
        )?;
        system
            .nonzero_conditions
            .extend(self.nonzero_conditions.iter().cloned());
        system.nonzero_conditions.sort();
        system.nonzero_conditions.dedup();
        Ok(system)
    }
}
