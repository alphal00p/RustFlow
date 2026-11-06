//! Whole-disk differential-defect majorants for rational Taylor charts.
//!
//! All source polynomial degrees are retained: a high-degree sparse connection
//! must not disappear merely because the requested solution order is smaller.
//! Symbolica owns polynomial differentiation, shifting and multiplication.
//! By default rational charts retain exact source coefficients and use native
//! directed ball arithmetic. Explicit NumericalEstimate charts instead use
//! finite MPFR arithmetic, with source conversion and residual rounding unenclosed.
//! Registered-root assembly retains its separate MPFR estimate scope.
//! These charts assess local differential defects, not
//! accumulated global solution error.
use crate::fixed_series::FixedComplexRing;
use crate::ode::{NumericRational, PolynomialRow};
use crate::{ComplexFloat as C, Error, Precision, Result};
use std::sync::Arc;
use symbolica::poly::univariate::UnivariatePolynomial;
use symbolica::prelude::*;

#[derive(Clone, Debug)]
pub(crate) struct RationalResidualChart {
    residuals: Vec<NumericRational>,
    exact_source: Option<super::integer_residual::AdaptiveResidual>,
    /// Only the explicit numerical mode; one denominator per physical row.
    numerical_denominators: Option<Vec<NumericalDenominator>>,
}

/// Finite-MPFR disk estimates, with certified source factor structure retained.
/// These triangle/product estimates do not enclose MPFR rounding.
#[derive(Clone, Debug)]
struct NumericalDenominator {
    expanded: Vec<C>,
    factors: Vec<(Vec<C>, usize)>,
}

fn numerical_lower(p: Precision, a: &[C], radius: &Float) -> Option<Float> {
    let mut lower = p.norm(a.first()?);
    let mut power = radius.clone();
    for coefficient in a.iter().skip(1) {
        lower -= p.norm(coefficient) * &power;
        power *= radius;
    }
    (lower.is_finite() && lower > p.real(0)).then_some(lower)
}

impl NumericalDenominator {
    fn lower(&self, p: Precision, radius: &Float) -> Option<Float> {
        let expanded = numerical_lower(p, &self.expanded, radius);
        // Empty factors are the certified unit denominator only (validated at
        // construction). Constant content is already one factor, used once.
        let mut factored = p.real(1);
        for (factor, multiplicity) in &self.factors {
            let Some(lower) = numerical_lower(p, factor, radius) else {
                // In particular, never multiply two negative triangle bounds.
                return expanded;
            };
            for _ in 0..*multiplicity {
                factored *= &lower;
            }
            if !factored.is_finite() || factored <= p.real(0) {
                return expanded;
            }
        }
        Some(
            expanded
                .filter(|value| *value > factored)
                .unwrap_or(factored),
        )
    }
}

impl RationalResidualChart {
    pub(crate) fn from_residuals(residuals: Vec<NumericRational>) -> Result<Self> {
        if residuals
            .iter()
            .any(|r| r.numerator.is_empty() || r.denominator.is_empty())
        {
            return Err(Error::InvalidInput("empty residual polynomial".into()));
        }
        Ok(Self {
            residuals,
            exact_source: None,
            numerical_denominators: None,
        })
    }
    pub(crate) fn new_with_options(
        p: Precision,
        rows: &[super::source::ExactPolynomialRow],
        center: &C,
        coefficients: &[Vec<C>],
        channels: usize,
        options: &crate::FlowOptions,
        context: &crate::RunContext,
    ) -> Result<Self> {
        let constructor = match options.residual_arithmetic {
            crate::ResidualArithmetic::Ball => super::integer_residual::AdaptiveResidual::ball,
            crate::ResidualArithmetic::AdaptiveInteger => {
                super::integer_residual::AdaptiveResidual::new
            }
            crate::ResidualArithmetic::NumericalEstimate => {
                return Self::numerical_estimate(p, rows, center, coefficients, channels, context);
            }
        };
        Ok(Self {
            residuals: Vec::new(),
            numerical_denominators: None,
            exact_source: Some(constructor(
                p,
                rows,
                center,
                coefficients,
                channels,
                &context.cancellation,
            )?),
        })
    }

    #[cfg(test)]
    pub(crate) fn selected_arithmetic(&self) -> Option<crate::ResidualArithmetic> {
        if self.numerical_denominators.is_some() {
            return Some(crate::ResidualArithmetic::NumericalEstimate);
        }
        self.exact_source.as_ref().map(|a| match a {
            super::integer_residual::AdaptiveResidual::Ball { .. } => {
                crate::ResidualArithmetic::Ball
            }
            super::integer_residual::AdaptiveResidual::Integer { .. } => {
                crate::ResidualArithmetic::AdaptiveInteger
            }
        })
    }

    pub(crate) fn defect_bounds_with_budget(
        &self,
        p: Precision,
        step: &C,
        values: &[C],
        tolerance: &Float,
    ) -> Result<Vec<Float>> {
        if let Some(source) = &self.exact_source {
            source.defect_bounds(p, step, Some((values, tolerance)))
        } else {
            self.defect_bounds(p, step)
        }
    }

    pub(crate) fn new(
        p: Precision,
        rows: &[PolynomialRow],
        center: &C,
        coefficients: &[Vec<C>],
        channels: usize,
    ) -> Result<Self> {
        let n = rows.len();
        let size = n
            .checked_mul(channels)
            .ok_or_else(|| Error::Limit("residual channel dimension overflow".into()))?;
        if n == 0
            || channels == 0
            || coefficients.is_empty()
            || coefficients.iter().any(|r| r.len() != size)
            || rows.iter().any(|r| {
                r.denominator.is_empty()
                    || r.entries.iter().any(|(j, a)| *j >= size || a.is_empty())
            })
        {
            return Err(Error::InvalidInput("residual polynomial dimensions".into()));
        }
        let ring = FixedComplexRing::new(p);
        let variable = Arc::new(PolyVariable::Temporary(0));
        let polynomial =
            |a: Vec<C>| UnivariatePolynomial::from_coefficients(&ring, a, variable.clone());
        let solutions = (0..size)
            .map(|i| polynomial(coefficients.iter().map(|r| r[i].clone()).collect()))
            .collect::<Vec<_>>();
        let mut residuals = vec![
            NumericRational {
                numerator: vec![p.zero()],
                denominator: vec![p.i(1)]
            };
            size
        ];
        for (i, row) in rows.iter().enumerate() {
            let denominator = polynomial(row.denominator.clone()).shift_var(center);
            let entries = row
                .entries
                .iter()
                .map(|(j, a)| (*j, polynomial(a.clone()).shift_var(center)))
                .collect::<Vec<_>>();
            for order in 0..channels {
                let index = order * n + i;
                let mut numerator = &denominator * &solutions[index].derivative();
                for (column, a) in &entries {
                    let shift = column / n;
                    if shift <= order {
                        let source = (order - shift) * n + column % n;
                        numerator = &numerator - &(a * &solutions[source]);
                    }
                }
                residuals[index] = NumericRational {
                    numerator: if numerator.is_zero() {
                        vec![p.zero()]
                    } else {
                        numerator.coefficients().to_vec()
                    },
                    denominator: if denominator.is_zero() {
                        vec![p.zero()]
                    } else {
                        denominator.coefficients().to_vec()
                    },
                };
            }
        }
        Ok(Self {
            residuals,
            exact_source: None,
            numerical_denominators: None,
        })
    }

    fn numerical_estimate(
        p: Precision,
        rows: &[super::source::ExactPolynomialRow],
        center: &C,
        coefficients: &[Vec<C>],
        channels: usize,
        context: &crate::RunContext,
    ) -> Result<Self> {
        context.cancellation.check()?;
        super::source::validate_inputs(p, rows, center, coefficients, channels)?;
        let convert = |a: &super::source::Gaussian| {
            C::new(
                a.re.to_multi_prec_float(p.bits),
                a.im.to_multi_prec_float(p.bits),
            )
        };
        let numeric_rows = rows
            .iter()
            .map(|row| PolynomialRow {
                denominator: row.denominator.iter().map(&convert).collect(),
                entries: row
                    .entries
                    .iter()
                    .map(|(j, a)| (*j, a.iter().map(&convert).collect()))
                    .collect(),
            })
            .collect::<Vec<_>>();
        if numeric_rows.iter().any(|row| {
            row.denominator
                .iter()
                .chain(row.entries.iter().flat_map(|(_, a)| a))
                .any(|a| !p.finite(a))
        }) {
            return Err(Error::Accuracy(
                "nonfinite numerical source conversion".into(),
            ));
        }
        // Reuse the existing full-degree native shift/derivative/product path.
        // The candidate and the rational Taylor recurrence are unchanged.
        let mut chart = Self::new(p, &numeric_rows, center, coefficients, channels)?;
        let ring = FixedComplexRing::new(p);
        let mut denominators = Vec::with_capacity(rows.len());
        for (i, row) in rows.iter().enumerate() {
            context.cancellation.check()?;
            let mut factors = Vec::with_capacity(row.denominator_factors.len());
            for (factor, multiplicity) in &row.denominator_factors {
                let polynomial = UnivariatePolynomial::from_coefficients(
                    &ring,
                    factor.iter().map(&convert).collect(),
                    Arc::new(PolyVariable::Temporary(0)),
                )
                .shift_var(center);
                if polynomial.coefficients().iter().any(|a| !p.finite(a)) {
                    return Err(Error::Accuracy(
                        "nonfinite numerical denominator factor".into(),
                    ));
                }
                factors.push((
                    if polynomial.is_zero() {
                        vec![p.zero()]
                    } else {
                        polynomial.coefficients().to_vec()
                    },
                    *multiplicity,
                ));
            }
            denominators.push(NumericalDenominator {
                expanded: chart.residuals[i].denominator.clone(),
                factors,
            });
        }
        // Denominators (including factors/content) are shared across channels.
        for residual in &mut chart.residuals {
            residual.denominator.clear();
        }
        chart.numerical_denominators = Some(denominators);
        context.cancellation.check()?;
        Ok(chart)
    }

    /// For every component, estimate |h| sup_{|z|<=|h|} |P'(z)-A(c+z)P(z)|.
    /// An inconclusive denominator lower bound requests a shorter step; it
    /// must never be converted into acceptance or a zero bound.
    pub(crate) fn defect_bounds(&self, p: Precision, step: &C) -> Result<Vec<Float>> {
        if let Some(source) = &self.exact_source {
            return source.defect_bounds(p, step, None);
        }
        let radius = p.norm(step);
        if !radius.is_finite() {
            return Err(Error::Accuracy("nonfinite residual disk radius".into()));
        }
        if let Some(denominators) = &self.numerical_denominators {
            let lowers = denominators
                .iter()
                .map(|d| {
                    d.lower(p, &radius).ok_or_else(|| {
                        Error::Accuracy(
                            "numerical denominator estimate is inconclusive on the trial disk"
                                .into(),
                        )
                    })
                })
                .collect::<Result<Vec<_>>>()?;
            return self
                .residuals
                .iter()
                .enumerate()
                .map(|(i, residual)| {
                    let mut upper = p.real(0);
                    let mut power = p.real(1);
                    for coefficient in &residual.numerator {
                        upper += p.norm(coefficient) * &power;
                        power *= &radius;
                    }
                    let bound = upper / &lowers[i % denominators.len()] * &radius;
                    if bound.is_finite() && bound >= p.real(0) {
                        Ok(bound)
                    } else {
                        Err(Error::Accuracy(
                            "nonfinite numerical differential-defect estimate".into(),
                        ))
                    }
                })
                .collect();
        }
        self.residuals
            .iter()
            .map(|residual| {
                let bound = crate::diffexp::rational_disk_bound(p, residual, &p.zero(), &radius)?
                    .ok_or_else(|| {
                    Error::Accuracy(
                        "cannot bound Taylor differential defect on the trial disk".into(),
                    )
                })? * &radius;
                if bound.is_finite() {
                    Ok(bound)
                } else {
                    Err(Error::Accuracy(
                        "nonfinite whole-segment differential defect".into(),
                    ))
                }
            })
            .collect()
    }
}

#[cfg(test)]
#[path = "numerical_residual_tests.rs"]
mod tests;
