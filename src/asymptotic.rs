//! Match selected coefficients of generalized power/logarithm solutions.
use crate::frobenius::FrobeniusBasis;
use crate::numeric::solve_constraints;
use crate::{ComplexFloat, Error, Result};
use symbolica::coefficient::Coefficient;
use symbolica::prelude::*;
mod exact;
pub use exact::{
    AsymptoticSelector, ExactAsymptoticConstraints, ExactAsymptoticRelation, ExactAsymptoticSpace,
};

/// A known coefficient of `x^power * log(x)^log_power` in one component.
///
/// Powers are compared symbolically, before numerical specialization. A power
/// belongs to a Frobenius column when their difference is an exact nonnegative
/// integer. Omitted coefficients remain unknown; they are not set to zero.
#[derive(Clone, Debug)]
pub struct AsymptoticConstraint {
    pub component: usize,
    pub power: Atom,
    pub log_power: usize,
    pub coefficient: ComplexFloat,
}

impl FrobeniusBasis {
    /// Determine integration constants from partial asymptotic information.
    ///
    /// Logarithmic and resonant columns contribute to the same equation when
    /// they have the specified power and logarithm. The constraints must fix
    /// every integration constant and must be mutually consistent at working
    /// precision. A coefficient beyond the computed series returns an accuracy
    /// error: missing series terms cannot be interpreted as exact zeros.
    ///
    /// Numerical agreement after increasing precision and series order is still
    /// required; a successful match is not an accuracy certificate.
    pub fn match_constraints(
        &self,
        constraints: &[AsymptoticConstraint],
    ) -> Result<Vec<ComplexFloat>> {
        self.validate()?;
        let n = self.columns.len();
        let p = self.precision;
        if constraints.is_empty() {
            return Err(Error::InvalidInput(
                "empty partial asymptotic matching problem".into(),
            ));
        }
        let exact_power = |power: &Atom| {
            crate::family::scalar_symbols(power.as_view(), &mut Default::default()).map_err(|_| {
                Error::InvalidInput("asymptotic powers must be exact scalar expressions".into())
            })
        };
        for column in &self.columns {
            exact_power(&column.exponent)?;
        }
        let mut matrix = Vec::with_capacity(constraints.len());
        let mut rhs = Vec::with_capacity(constraints.len());
        for constraint in constraints {
            if constraint.component >= n || !p.finite(&constraint.coefficient) {
                return Err(Error::InvalidInput(
                    "invalid asymptotic coefficient constraint".into(),
                ));
            }
            exact_power(&constraint.power)?;
            let mut row = vec![p.zero(); n];
            for (j, column) in self.columns.iter().enumerate() {
                let difference = (&constraint.power - &column.exponent).together().cancel();
                let AtomView::Num(number) = difference.as_view() else {
                    continue;
                };
                let Coefficient::Complex(value) = number.get_coeff_view().to_owned() else {
                    continue;
                };
                if !value.im.is_zero() || !value.re.is_integer() || value.re < 0 {
                    continue;
                }
                let offset = value.re.to_string().parse::<usize>().map_err(|_| {
                    Error::Accuracy(
                        "asymptotic coefficient order exceeds addressable series".into(),
                    )
                })?;
                let coefficient = column.coefficients.get(offset).ok_or_else(|| {
                    Error::Accuracy(
                        "Frobenius basis is too short for the requested asymptotic constraint"
                            .into(),
                    )
                })?;
                if let Some(log) = coefficient.get(constraint.log_power) {
                    row[j] = log[constraint.component].clone();
                }
            }
            matrix.push(row);
            rhs.push(constraint.coefficient.clone());
        }
        solve_constraints(p, matrix, rhs, n)
    }
}
