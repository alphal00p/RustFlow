#![allow(clippy::needless_range_loop)] // Exact quadratic forms use explicit row/column coordinates.
//! Hard/soft loop regions and exact large-mass integrand expansions.
//!
//! Region enumeration includes independent branch momentum routings, rather than
//! just hard/soft assignments in the caller's original loop coordinates.
use crate::algebra::{determinant, inverse, matmul};
use crate::{Error, Integral, IntegralFamily, Propagator, Result};
use std::collections::BTreeSet;
use symbolica::prelude::*;
use symbolica::utils::resources;

fn check_region_allocation<T>(count: usize) -> Result<()> {
    resources::check_allocation(
        u64::try_from(count.saturating_mul(std::mem::size_of::<T>())).unwrap_or(u64::MAX),
    )?;
    Ok(())
}

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
    check_region_allocation::<Atom>(loops.saturating_mul(loops))?;
    let q = quadratic_matrix(p, loops);
    let pivot = (0..loops)
        .find(|&i| !q[i][i].is_zero())
        .ok_or_else(|| Error::Unsupported("propagator has no squared loop momentum".into()))?;
    let vector = q[pivot]
        .iter()
        .map(|a| (a / &q[pivot][pivot]).together().cancel())
        .collect::<Vec<_>>();
    for i in 0..loops {
        resources::checkpoint()?;
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

/// Enumerate loop regions, with an optional diagnostic cap on routing/assignment
/// candidates. `None` leaves the count unrestricted; active resource budgets
/// still apply to the whole enumeration.
pub fn enumerate_regions(
    family: &IntegralFamily,
    max_candidates: Option<usize>,
) -> Result<Vec<LoopRegion>> {
    resources::checkpoint()?;
    family.validate()?;
    let loops = family.loops.len();
    let assignments = u32::try_from(loops)
        .ok()
        .and_then(|loops| 1_usize.checked_shl(loops))
        .ok_or_else(|| Error::Limit("too many hard/soft assignments".into()))?;
    let mut branches = Vec::new();
    for p in &family.propagators[..family.physical_propagators] {
        resources::checkpoint()?;
        let b = branch(p, loops)?;
        if !branches.contains(&b) {
            check_region_allocation::<Vec<Atom>>(1)?;
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
        limit: Option<usize>,
    ) -> Result<()> {
        crate::stack_resources::run(|| {
            resources::checkpoint()?;
            if current.len() == k {
                if limit.is_some_and(|limit| out.len() >= limit) {
                    return Err(Error::Limit("region routing budget exhausted".into()));
                }
                check_region_allocation::<Vec<usize>>(1)?;
                check_region_allocation::<usize>(current.len())?;
                out.push(current.clone());
                return Ok(());
            }
            for i in start..n {
                check_region_allocation::<usize>(1)?;
                current.push(i);
                choose(n, k, i + 1, current, out, limit)?;
                current.pop();
            }
            Ok(())
        })
    }
    if max_candidates.is_some_and(|max_candidates| assignments > max_candidates) {
        return Err(Error::Limit("hard/soft assignment budget exhausted".into()));
    }
    choose(
        branches.len(),
        loops,
        0,
        &mut Vec::new(),
        &mut choices,
        max_candidates.map(|max_candidates| max_candidates / assignments),
    )?;
    let mut seen = BTreeSet::new();
    let mut regions = Vec::new();
    for choice in choices {
        resources::checkpoint()?;
        check_region_allocation::<Atom>(loops.saturating_mul(loops))?;
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
            resources::checkpoint()?;
            check_region_allocation::<bool>(loops.saturating_add(transformed.len()))?;
            let hard = (0..loops).map(|i| mask & (1 << i) != 0).collect::<Vec<_>>();
            let signature = transformed
                .iter()
                .map(|r| r.iter().zip(&hard).any(|(a, h)| *h && !a.is_zero()))
                .collect::<Vec<_>>();
            if seen.insert(signature.clone()) {
                check_region_allocation::<LoopRegion>(1)?;
                check_region_allocation::<Atom>(loops.saturating_mul(loops))?;
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
    resources::checkpoint()?;
    Ok(regions)
}

impl IntegralFamily {
    pub fn transform_loops(&self, transformation: &[Vec<Atom>]) -> Result<Self> {
        resources::checkpoint()?;
        let loops = self.loops.len();
        if transformation.len() != loops
            || transformation.iter().any(|r| r.len() != loops)
            || determinant(transformation.to_vec()).is_zero()
        {
            return Err(Error::InvalidInput("singular loop transformation".into()));
        }
        let transposed = transpose(transformation);
        check_region_allocation::<Propagator>(self.propagators.len())?;
        let mut out = self.clone();
        let offset = loops * (loops + 1) / 2;
        for (old, new) in self.propagators.iter().zip(&mut out.propagators) {
            resources::checkpoint()?;
            check_region_allocation::<Atom>(loops.saturating_mul(loops).saturating_mul(3))?;
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
    resources::checkpoint()?;
    family.validate_integral(integral)?;
    if shifted.len() != family.propagators.len() || region.hard.len() != family.loops.len() {
        return Err(Error::InvalidInput("region mask dimensions".into()));
    }
    let coefficient_count = order
        .checked_add(1)
        .ok_or_else(|| Error::Limit("boundary integrand order overflow".into()))?;
    let series_order = i64::try_from(coefficient_count)
        .map_err(|_| Error::Limit("boundary integrand order exceeds i64".into()))?;
    check_region_allocation::<Atom>(coefficient_count)?;
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
        resources::checkpoint()?;
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
        .series(t, 0, series_order)
        .map_err(|e| Error::Unsupported(e.to_string()))?;
    resources::checkpoint()?;
    let hard = region.hard.iter().filter(|&&h| h).count() as i64;
    let eta_power = (Atom::num(hard)
        * (Atom::num((family.dimension, 2)) - Atom::var(family.epsilon))
        - Atom::num((scaling, 2)))
    .together()
    .cancel();
    Ok(RegionExpansion {
        eta_power,
        coordinates,
        coefficients: (0..coefficient_count)
            .map(|i| {
                resources::checkpoint()?;
                Ok(expanded
                    .coefficient(Rational::from(i as i64))
                    .unwrap_or_default())
            })
            .collect::<Result<_>>()?,
        jacobian_determinant: region.jacobian_determinant.clone(),
    })
}
