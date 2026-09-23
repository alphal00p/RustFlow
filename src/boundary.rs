//! Boundary providers and matching of asymptotic data to a Frobenius basis.
use crate::frobenius::FrobeniusBasis;
use crate::numeric::solve_constraints;
use crate::{ComplexFloat as C, Error, Integral, IntegralFamily, Precision, Result};
use symbolica::prelude::*;

#[derive(Clone, Debug)]
pub struct LeadingBoundary {
    pub exponent: Atom,
    pub coefficient: C,
}

pub trait BoundaryProvider: Send + Sync {
    fn constants(
        &self,
        family: &IntegralFamily,
        basis: &[Integral],
        epsilon: &Rational,
        p: Precision,
        solutions: &FrobeniusBasis,
    ) -> Result<Vec<C>> {
        solutions.match_leading(&self.leading(family, basis, epsilon, p)?)
    }
    /// Data at z=1/eta -> 0; one entry per integral in the supplied basis.
    fn leading(
        &self,
        family: &IntegralFamily,
        basis: &[Integral],
        epsilon: &Rational,
        p: Precision,
    ) -> Result<Vec<LeadingBoundary>>;
}

/// All-hard boundary for one-loop families with unit quadratic coefficients.
/// Its coefficients are analytic massive-vacuum integrals, not fitted values.
#[derive(Clone, Copy, Debug, Default)]
pub struct OneLoopBoundary;
impl BoundaryProvider for OneLoopBoundary {
    fn leading(
        &self,
        family: &IntegralFamily,
        basis: &[Integral],
        epsilon: &Rational,
        p: Precision,
    ) -> Result<Vec<LeadingBoundary>> {
        if family.loops.len() != 1 {
            return Err(Error::Unsupported(
                "OneLoopBoundary supports one-loop families; use RecursiveBoundary for multiloop families".into(),
            ));
        }
        for d in &family.propagators {
            if d.scalar_products.first() != Some(&Atom::num(1)) {
                return Err(Error::Unsupported(
                    "one-loop boundary requires unit coefficient of l^2".into(),
                ));
            }
        }
        let half_dimension = Rational::from((family.dimension, 2)) - epsilon;
        basis
            .iter()
            .map(|integral| {
                let power = integral.0[..family.physical_propagators]
                    .iter()
                    .map(|&v| v as i64)
                    .sum::<i64>();
                let rank = -integral.0[family.physical_propagators..]
                    .iter()
                    .map(|&v| v as i64)
                    .sum::<i64>();
                if power <= 0 {
                    return Ok(LeadingBoundary {
                        exponent: Atom::num(power - rank) - Atom::num((family.dimension, 2))
                            + Atom::var(family.epsilon),
                        coefficient: p.zero(),
                    });
                }
                let gamma = |q: Rational| p.gamma_real(&p.rational(&q).re);
                let coefficient = p.div(
                    &p.mul(
                        &gamma(&half_dimension + &Rational::from(rank))?,
                        &gamma(Rational::from(power - rank) - &half_dimension)?,
                    ),
                    &p.mul(
                        &gamma(half_dimension.clone())?,
                        &gamma(Rational::from(power))?,
                    ),
                );
                let coefficient = if (power + rank) % 2 == 0 {
                    coefficient
                } else {
                    p.neg(&coefficient)
                };
                Ok(LeadingBoundary {
                    exponent: Atom::num(power - rank) - Atom::num((family.dimension, 2))
                        + Atom::var(family.epsilon),
                    coefficient,
                })
            })
            .collect()
    }
}

impl FrobeniusBasis {
    /// Fix constants from asymptotic leading coefficients and absence of other
    /// power sectors. Only valid when the provider supplies all contributing regions.
    pub fn match_leading(&self, leading: &[LeadingBoundary]) -> Result<Vec<C>> {
        let p = self.precision;
        let n = self.columns.len();
        if leading.len() != n {
            return Err(Error::InvalidInput("boundary basis dimension".into()));
        }
        let mut equations = Vec::new();
        let mut rhs = Vec::new();
        for (j, column) in self.columns.iter().enumerate() {
            if !leading.iter().any(|b| {
                (&b.exponent - &column.exponent)
                    .together()
                    .cancel()
                    .to_string()
                    .parse::<i64>()
                    .is_ok()
            }) {
                let mut row = vec![p.zero(); n];
                row[j] = p.i(1);
                equations.push(row);
                rhs.push(p.zero());
            }
        }
        for (i, boundary) in leading.iter().enumerate() {
            let mut max_preceding = 0;
            for column in &self.columns {
                if let Ok(k) = (&boundary.exponent - &column.exponent)
                    .together()
                    .cancel()
                    .to_string()
                    .parse::<i64>()
                {
                    max_preceding = max_preceding.max(k);
                }
            }
            if max_preceding > self.columns[0].coefficients.len() as i64 {
                return Err(Error::Accuracy("insufficient boundary series order".into()));
            }
            for before in 0..=max_preceding.max(0) {
                let mut row = vec![p.zero(); n];
                for (j, column) in self.columns.iter().enumerate() {
                    if let Ok(k) = (&boundary.exponent - &column.exponent - Atom::num(before))
                        .together()
                        .cancel()
                        .to_string()
                        .parse::<i64>()
                        && k >= 0
                    {
                        let coeff = column.coefficients.get(k as usize).ok_or_else(|| {
                            Error::Accuracy("boundary expansion too short".into())
                        })?;
                        row[j] = coeff[0][i].clone();
                    }
                }
                equations.push(row);
                rhs.push(if before == 0 {
                    boundary.coefficient.clone()
                } else {
                    p.zero()
                });
            }
        }
        solve_constraints(p, equations, rhs, n)
    }
}

/// One Taylor series within a dimensional region at z=1/eta=0.
#[derive(Clone, Debug)]
pub struct RegionBoundary {
    pub exponent: Atom,
    pub coefficients: Vec<C>,
}

impl FrobeniusBasis {
    /// Match every supplied power and require logarithmic coefficients to vanish.
    /// Distinct epsilon-dependent exponent classes are kept separate.
    pub fn match_regions(&self, data: &[Vec<RegionBoundary>]) -> Result<Vec<C>> {
        let p = self.precision;
        let n = self.columns.len();
        if data.len() != n {
            return Err(Error::InvalidInput("region boundary dimension".into()));
        }
        let mut rows = Vec::new();
        let mut rhs = Vec::new();
        for (i, regions) in data.iter().enumerate() {
            let mut powers = std::collections::BTreeSet::new();
            for region in regions {
                for k in 0..region.coefficients.len() {
                    powers.insert((&region.exponent + Atom::num(k as i64)).together().cancel());
                }
            }
            // A region starts at its stated leading power. Lower powers in
            // the same indicial class are known zeros, not missing data.
            for column in &self.columns {
                let end = regions
                    .iter()
                    .filter_map(|r| {
                        (&r.exponent + Atom::num(r.coefficients.len() as i64 - 1)
                            - &column.exponent)
                            .together()
                            .cancel()
                            .to_string()
                            .parse::<usize>()
                            .ok()
                    })
                    .max();
                if let Some(end) = end {
                    if end >= column.coefficients.len() {
                        return Err(Error::Accuracy("insufficient boundary series".into()));
                    }
                    for k in 0..=end {
                        powers.insert((&column.exponent + Atom::num(k as i64)).together().cancel());
                    }
                }
            }
            for exponent in powers {
                // Do not mistake an uncomputed subleading term of another
                // region in the same class for a known zero.
                if regions.iter().any(|r| {
                    (&exponent - &r.exponent)
                        .together()
                        .cancel()
                        .to_string()
                        .parse::<usize>()
                        .is_ok_and(|k| k >= r.coefficients.len())
                }) {
                    continue;
                }
                let mut desired = p.zero();
                for region in regions {
                    if let Ok(k) = (&exponent - &region.exponent)
                        .together()
                        .cancel()
                        .to_string()
                        .parse::<usize>()
                        && let Some(value) = region.coefficients.get(k)
                    {
                        desired = p.add(&desired, value);
                    }
                }
                let mut logarithms = vec![vec![p.zero(); n]];
                for (j, column) in self.columns.iter().enumerate() {
                    if let Ok(k) = (&exponent - &column.exponent)
                        .together()
                        .cancel()
                        .to_string()
                        .parse::<usize>()
                    {
                        let terms = column.coefficients.get(k).ok_or_else(|| {
                            Error::Accuracy("insufficient Frobenius boundary order".into())
                        })?;
                        logarithms.resize(logarithms.len().max(terms.len()), vec![p.zero(); n]);
                        for (l, row) in terms.iter().enumerate() {
                            logarithms[l][j] = row[i].clone();
                        }
                    }
                }
                for (l, row) in logarithms.into_iter().enumerate() {
                    rows.push(row);
                    rhs.push(if l == 0 { desired.clone() } else { p.zero() });
                }
            }
        }
        for (j, column) in self.columns.iter().enumerate() {
            if !data.iter().flatten().any(|r| {
                (&r.exponent - &column.exponent)
                    .together()
                    .cancel()
                    .to_string()
                    .parse::<i64>()
                    .is_ok()
            }) {
                let mut row = vec![p.zero(); n];
                row[j] = p.i(1);
                rows.push(row);
                rhs.push(p.zero());
            }
        }
        solve_constraints(p, rows, rhs, n)
    }
}
