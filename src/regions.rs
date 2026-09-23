#![allow(clippy::needless_range_loop)] // Exact quadratic forms use explicit row/column coordinates.
//! Hard/soft loop regions and exact large-mass integrand expansions.
//!
//! Region enumeration includes independent branch momentum routings, rather than
//! just hard/soft assignments in the caller's original loop coordinates.
use crate::algebra::{determinant, inverse, matmul};
use crate::{Error, Integral, IntegralFamily, Propagator, Result};
use std::collections::BTreeSet;
use symbolica::prelude::*;

#[derive(Clone, Debug)]
pub struct LoopRegion {
    /// old_loop[i] = sum_j transformation[i][j] new_loop[j].
    pub transformation: Vec<Vec<Atom>>,
    pub hard: Vec<bool>,
    pub hard_branches: Vec<bool>,
    pub jacobian_determinant: Atom,
}

#[derive(Clone, Debug)]
pub struct RegionExpansion {
    /// Overall eta power, before the series in t = eta^(-1/2).
    pub eta_power: Atom,
    pub coordinates: Vec<Atom>,
    /// Coefficients of t^0, t^1, ... . The integration measure Jacobian is
    /// |det(transformation)|^(D0 - 2 epsilon) and is kept separately.
    pub coefficients: Vec<Atom>,
    pub jacobian_determinant: Atom,
}

fn quadratic_matrix(p: &Propagator, loops: usize) -> Vec<Vec<Atom>> {
    let mut matrix = vec![vec![Atom::new(); loops]; loops];
    let mut k = 0;
    for i in 0..loops {
        for j in i..loops {
            let a = if i == j {
                p.scalar_products[k].clone()
            } else {
                &p.scalar_products[k] / Atom::num(2)
            };
            matrix[i][j] = a.clone();
            matrix[j][i] = a;
            k += 1;
        }
    }
    matrix
}

fn transpose(a: &[Vec<Atom>]) -> Vec<Vec<Atom>> {
    (0..a[0].len())
        .map(|j| a.iter().map(|r| r[j].clone()).collect())
        .collect()
}

pub(crate) fn branch(p: &Propagator, loops: usize) -> Result<Vec<Atom>> {
    let q = quadratic_matrix(p, loops);
    let pivot = (0..loops)
        .find(|&i| !q[i][i].is_zero())
        .ok_or_else(|| Error::Unsupported("propagator has no squared loop momentum".into()))?;
    let vector = q[pivot]
        .iter()
        .map(|a| (a / &q[pivot][pivot]).together().cancel())
        .collect::<Vec<_>>();
    for i in 0..loops {
        for j in 0..loops {
            if !(&q[i][j] - &q[pivot][pivot] * &vector[i] * &vector[j])
                .together()
                .cancel()
                .is_zero()
            {
                return Err(Error::Unsupported(
                    "region analysis requires ordinary rank-one quadratic propagators".into(),
                ));
            }
        }
    }
    Ok(vector)
}

pub fn enumerate_regions(
    family: &IntegralFamily,
    max_candidates: usize,
) -> Result<Vec<LoopRegion>> {
    family.validate()?;
    let loops = family.loops.len();
    let assignments = 1_usize
        .checked_shl(loops as u32)
        .ok_or_else(|| Error::Limit("too many hard/soft assignments".into()))?;
    let mut branches = Vec::new();
    for p in &family.propagators[..family.physical_propagators] {
        let b = branch(p, loops)?;
        if !branches.contains(&b) {
            branches.push(b);
        }
    }
    let mut choices = Vec::new();
    fn choose(
        n: usize,
        k: usize,
        start: usize,
        current: &mut Vec<usize>,
        out: &mut Vec<Vec<usize>>,
        limit: usize,
    ) -> Result<()> {
        if current.len() == k {
            if out.len() >= limit {
                return Err(Error::Limit("region routing budget exhausted".into()));
            }
            out.push(current.clone());
            return Ok(());
        }
        for i in start..n {
            current.push(i);
            choose(n, k, i + 1, current, out, limit)?;
            current.pop();
        }
        Ok(())
    }
    if assignments > max_candidates {
        return Err(Error::Limit("hard/soft assignment budget exhausted".into()));
    }
    choose(
        branches.len(),
        loops,
        0,
        &mut Vec::new(),
        &mut choices,
        max_candidates / assignments,
    )?;
    let mut seen = BTreeSet::new();
    let mut regions = Vec::new();
    for choice in choices {
        let routing = choice
            .iter()
            .map(|&i| branches[i].clone())
            .collect::<Vec<_>>();
        if determinant(routing.clone()).is_zero() {
            continue;
        }
        let transformation = inverse(&routing)?;
        let transformed = matmul(&branches, &transformation);
        for mask in 0..assignments {
            let hard = (0..loops).map(|i| mask & (1 << i) != 0).collect::<Vec<_>>();
            let signature = transformed
                .iter()
                .map(|r| r.iter().zip(&hard).any(|(a, h)| *h && !a.is_zero()))
                .collect::<Vec<_>>();
            if seen.insert(signature.clone()) {
                regions.push(LoopRegion {
                    jacobian_determinant: determinant(transformation.clone()),
                    transformation: transformation.clone(),
                    hard,
                    hard_branches: signature,
                });
            }
        }
    }
    if regions.is_empty() {
        return Err(Error::InvalidInput(
            "loop branch momenta do not span loop space".into(),
        ));
    }
    Ok(regions)
}

impl IntegralFamily {
    pub fn transform_loops(&self, transformation: &[Vec<Atom>]) -> Result<Self> {
        let loops = self.loops.len();
        if transformation.len() != loops
            || transformation.iter().any(|r| r.len() != loops)
            || determinant(transformation.to_vec()).is_zero()
        {
            return Err(Error::InvalidInput("singular loop transformation".into()));
        }
        let transposed = transpose(transformation);
        let mut out = self.clone();
        let offset = loops * (loops + 1) / 2;
        for (old, new) in self.propagators.iter().zip(&mut out.propagators) {
            let q = matmul(
                &matmul(&transposed, &quadratic_matrix(old, loops)),
                transformation,
            );
            let mut coefficients = Vec::new();
            for (i, row) in q.iter().enumerate() {
                for (j, a) in row.iter().enumerate().skip(i) {
                    coefficients.push(if i == j { a.clone() } else { Atom::num(2) * a });
                }
            }
            for j in 0..loops {
                for e in 0..self.external.len() {
                    coefficients.push(
                        (0..loops)
                            .fold(Atom::new(), |sum, i| {
                                sum + &transformation[i][j]
                                    * &old.scalar_products[offset + i * self.external.len() + e]
                            })
                            .together()
                            .cancel(),
                    );
                }
            }
            new.scalar_products = coefficients;
        }
        Ok(out)
    }
}

pub fn expand_region(
    family: &IntegralFamily,
    integral: &Integral,
    shifted: &[bool],
    region: &LoopRegion,
    order: usize,
) -> Result<RegionExpansion> {
    family.validate_integral(integral)?;
    if shifted.len() != family.propagators.len() || region.hard.len() != family.loops.len() {
        return Err(Error::InvalidInput("region mask dimensions".into()));
    }
    if order > 100 {
        return Err(Error::Limit("boundary integrand order exceeds 100".into()));
    }
    let transformed = family.transform_loops(&region.transformation)?;
    let loops = family.loops.len();
    let mut coordinate_powers = Vec::new();
    for i in 0..loops {
        for j in i..loops {
            coordinate_powers.push(u8::from(region.hard[i]) + u8::from(region.hard[j]));
        }
    }
    for i in 0..loops {
        for _ in &family.external {
            coordinate_powers.push(u8::from(region.hard[i]));
        }
    }
    let coordinates = (0..coordinate_powers.len())
        .map(|i| Atom::var(symbol!(format!("symbolica_amflow::sp_{i}"))))
        .collect::<Vec<_>>();
    let t = symbol!("symbolica_amflow::region_t");
    let mut integrand = Atom::num(1);
    let mut scaling = 0_i64;
    for ((p, &nu), &q) in transformed.propagators.iter().zip(&integral.0).zip(shifted) {
        if nu == 0 {
            continue;
        }
        let mut terms = [
            p.constant.clone(),
            Atom::new(),
            if q { Atom::num(-1) } else { Atom::new() },
        ];
        for ((c, sp), &power) in p
            .scalar_products
            .iter()
            .zip(&coordinates)
            .zip(&coordinate_powers)
        {
            terms[power as usize] = &terms[power as usize] + c * sp;
        }
        let degree = (0..=2)
            .rev()
            .find(|&k| !terms[k].is_zero())
            .ok_or_else(|| Error::InvalidInput("zero propagator in region".into()))?;
        let polynomial = (0..=degree).fold(Atom::new(), |s, k| {
            s + &terms[k] * Atom::var(t).pow((degree - k) as i64)
        });
        integrand *= polynomial.pow(-(nu as i64));
        scaling += degree as i64 * nu as i64;
    }
    let expanded = integrand
        .series(t, 0, (order + 1) as i64)
        .map_err(|e| Error::Unsupported(e.to_string()))?;
    let hard = region.hard.iter().filter(|&&h| h).count() as i64;
    let eta_power = (Atom::num(hard)
        * (Atom::num((family.dimension, 2)) - Atom::var(family.epsilon))
        - Atom::num((scaling, 2)))
    .together()
    .cancel();
    Ok(RegionExpansion {
        eta_power,
        coordinates,
        coefficients: (0..=order)
            .map(|i| {
                expanded
                    .coefficient(Rational::from(i as i64))
                    .unwrap_or_default()
            })
            .collect(),
        jacobian_determinant: region.jacobian_determinant.clone(),
    })
}
