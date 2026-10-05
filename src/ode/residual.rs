//! Whole-disk differential-defect majorants for rational Taylor charts.
//!
//! All source polynomial degrees are retained: a high-degree sparse connection
//! must not disappear merely because the requested solution order is smaller.
//! Symbolica owns polynomial differentiation, shifting and multiplication.
//! Ordinary rational charts retain exact source coefficients and use native
//! directed ball arithmetic. Registered-root assembly currently retains its
//! separate MPFR estimate scope. Both bound local differential defects, not
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
    exact_source: Option<super::source::ExactSourceResidual>,
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
        })
    }
    pub(crate) fn new_exact(
        p: Precision,
        rows: &[super::source::ExactPolynomialRow],
        center: &C,
        coefficients: &[Vec<C>],
        channels: usize,
    ) -> Result<Self> {
        Ok(Self {
            residuals: Vec::new(),
            exact_source: Some(super::source::ExactSourceResidual::new(
                p,
                rows,
                center,
                coefficients,
                channels,
            )?),
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
        })
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
