#![allow(clippy::needless_range_loop)] // Momentum and Gram matrix coordinates.
//! Exact massless two-point subloop integration, normalized to its scalar
//! (1,1) integral. Integer shifts of the beta-function factors are evaluated
//! as rational Pochhammer products, so no numerical gamma functions enter.
use crate::algebra::{inverse, rref};
use crate::family::substitute;
use crate::integrand::powers;
use crate::reduction::LinearCombination;
use crate::{Error, Integral, IntegralFamily, Result};
use std::collections::BTreeMap;
use symbolica::prelude::*;

#[derive(Clone)]
pub(crate) struct Bubble {
    lines: [usize; 2],
    dependent: Vec<bool>,
    shifted: Vec<Atom>,
    z: Vec<Atom>,
    denominators: Vec<Atom>,
    gram: Vec<Vec<Atom>>,
    qdot: Vec<Atom>,
    q_squared: Atom,
    q_line: Option<usize>,
    dimension: Atom,
}

fn clean(a: Atom) -> Atom {
    a.together().cancel()
}

fn pochhammer(z: &Atom, n: i64) -> Option<Atom> {
    if n >= 0 {
        Some((0..n).fold(Atom::num(1), |a, k| a * (z + Atom::num(k))))
    } else {
        let denominator = clean((1..=-n).fold(Atom::num(1), |a, k| a * (z - Atom::num(k))));
        (!denominator.is_zero()).then(|| Atom::num(1) / denominator)
    }
}

fn tensor_coefficient(d: &Atom, a: i64, b: i64, r: i64, j: i64) -> Option<Atom> {
    let half = d / Atom::num(2);
    let factorial = |n| (1..=n).fold(Atom::num(1), |a, k| a * Atom::num(k));
    let numerator = pochhammer(&(Atom::num(2) - &half), a + b - j - 2)?
        * pochhammer(&(&half - Atom::num(1)), 1 - a + r - j)?
        * pochhammer(&(&half - Atom::num(1)), 1 - b + j)?;
    let denominator = clean(
        Atom::num(2).pow(j)
            * factorial(a - 1)
            * factorial(b - 1)
            * pochhammer(&(d - Atom::num(2)), 2 - a - b + r)?,
    );
    (!denominator.is_zero())
        .then(|| clean(Atom::num(if r % 2 == 0 { 1 } else { -1 }) * numerator / denominator))
}

// Sum all partial pairings. Unpaired vectors contract with q; each paired
// pair contracts with the metric. Index j counts pairs, not tensor rank.
fn contractions(indices: &[usize], gram: &[Vec<Atom>], qdot: &[Atom]) -> Vec<Atom> {
    let Some((&first, rest)) = indices.split_first() else {
        return vec![Atom::num(1)];
    };
    let mut out = vec![Atom::new(); indices.len() / 2 + 1];
    for (j, a) in contractions(rest, gram, qdot).into_iter().enumerate() {
        out[j] += &qdot[first] * a;
    }
    for (position, &second) in rest.iter().enumerate() {
        let remaining = rest[..position]
            .iter()
            .chain(&rest[position + 1..])
            .copied()
            .collect::<Vec<_>>();
        for (j, a) in contractions(&remaining, gram, qdot).into_iter().enumerate() {
            out[j + 1] += &gram[first][second] * a;
        }
    }
    out
}

impl Bubble {
    /// Only accepts a massless bubble whose transfer square is external or
    /// an already active denominator. Thus it introduces no new propagators.
    pub(crate) fn find(
        family: &IntegralFamily,
        sector: &[bool],
        dimension: Atom,
    ) -> Result<Option<Self>> {
        if family
            .propagators
            .iter()
            .flat_map(|p| std::iter::once(&p.constant).chain(&p.scalar_products))
            .chain(family.external_gram.iter().flatten())
            .any(|a| crate::family::encode_complex(a) != *a)
        {
            return Ok(None);
        }
        let l = family.loops.len();
        let e = family.external.len();
        let coordinates = (0..l)
            .flat_map(|i| (i..l).map(move |j| (i, j)))
            .chain((0..l).flat_map(|i| (0..e).map(move |j| (i, l + j))))
            .collect::<Vec<_>>();
        let x = (0..coordinates.len())
            .map(|i| Atom::var(symbol!(&format!("symbolica_amflow::bubble_x{i}"))))
            .collect::<Vec<_>>();
        let denominators = (0..family.propagators.len())
            .map(|i| Atom::var(symbol!(&format!("symbolica_amflow::bubble_d{i}"))))
            .collect::<Vec<_>>();
        let expressions = family
            .propagators
            .iter()
            .map(|p| {
                p.scalar_products
                    .iter()
                    .zip(&x)
                    .fold(p.constant.clone(), |a, (c, x)| a + c * x)
            })
            .collect::<Vec<_>>();
        let mut gram = vec![vec![Atom::new(); l + e]; l + e];
        for (&(i, j), value) in coordinates.iter().zip(&x) {
            gram[i][j] = value.clone();
            gram[j][i] = value.clone();
        }
        for i in 0..e {
            for j in 0..e {
                gram[l + i][l + j] = family.external_gram[i][j].clone();
            }
        }
        let dot = |a: &[Atom], b: &[Atom]| {
            a.iter().enumerate().fold(Atom::new(), |sum, (i, a)| {
                b.iter()
                    .enumerate()
                    .fold(sum, |sum, (j, b)| sum + a * b * &gram[i][j])
            })
        };
        for loop_index in 0..l {
            let dependent = family
                .propagators
                .iter()
                .map(|p| {
                    coordinates
                        .iter()
                        .zip(&p.scalar_products)
                        .any(|(&(i, j), c)| (i == loop_index || j == loop_index) && !c.is_zero())
                })
                .collect::<Vec<_>>();
            let active = dependent
                .iter()
                .zip(sector)
                .enumerate()
                .filter_map(|(i, (&d, &s))| (d && s).then_some(i))
                .collect::<Vec<_>>();
            let [first, second] = active.as_slice() else {
                continue;
            };
            let route = |line: usize| -> Option<Vec<Atom>> {
                let p = &family.propagators[line];
                let diagonal = coordinates
                    .iter()
                    .position(|&ij| ij == (loop_index, loop_index))?;
                if !p.scalar_products[diagonal].is_one() {
                    return None;
                }
                let mut v = vec![Atom::new(); l + e];
                v[loop_index] = Atom::num(1);
                for (&(i, j), c) in coordinates.iter().zip(&p.scalar_products) {
                    if i == loop_index && j != loop_index {
                        v[j] = c / Atom::num(2);
                    } else if j == loop_index && i != loop_index {
                        v[i] = c / Atom::num(2);
                    }
                }
                clean(dot(&v, &v) - &expressions[line])
                    .is_zero()
                    .then_some(v)
            };
            let (Some(mut shift), Some(other)) = (route(*first), route(*second)) else {
                continue;
            };
            let q = other
                .iter()
                .zip(&shift)
                .map(|(a, b)| a - b)
                .collect::<Vec<_>>();
            shift[loop_index] = Atom::new();
            let q_squared = clean(dot(&q, &q));
            if q_squared.is_zero() {
                continue;
            }
            let q_line = expressions.iter().enumerate().find_map(|(i, a)| {
                (!dependent[i] && sector[i] && clean(a - &q_squared).is_zero()).then_some(i)
            });
            let remaining = coordinates
                .iter()
                .enumerate()
                .filter_map(|(k, &(i, j))| (i != loop_index && j != loop_index).then_some(k))
                .collect::<Vec<_>>();
            if q_line.is_none()
                && remaining
                    .iter()
                    .any(|&k| !q_squared.coefficient(&x[k]).is_zero())
            {
                continue;
            }
            let pure = dependent
                .iter()
                .enumerate()
                .filter_map(|(i, &d)| (!d).then_some(i))
                .collect::<Vec<_>>();
            let matrix = remaining
                .iter()
                .map(|&k| {
                    pure.iter()
                        .map(|&i| family.propagators[i].scalar_products[k].clone())
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let (_, pivots) = rref(matrix);
            if pivots.len() != remaining.len() {
                continue;
            }
            let chosen = pivots.iter().map(|&i| pure[i]).collect::<Vec<_>>();
            let mut replacement = BTreeMap::new();
            if !remaining.is_empty() {
                let inv = inverse(
                    &chosen
                        .iter()
                        .map(|&i| {
                            remaining
                                .iter()
                                .map(|&k| family.propagators[i].scalar_products[k].clone())
                                .collect()
                        })
                        .collect::<Vec<_>>(),
                )?;
                for (&k, row) in remaining.iter().zip(inv) {
                    let value = row.iter().zip(&chosen).fold(Atom::new(), |a, (c, &i)| {
                        a + c * (&denominators[i] - &family.propagators[i].constant)
                    });
                    replacement.insert(x[k].clone(), clean(value));
                }
            }
            let vectors = (0..l + e).filter(|&i| i != loop_index).collect::<Vec<_>>();
            let z = (0..=vectors.len())
                .map(|i| Atom::var(symbol!(&format!("symbolica_amflow::bubble_z{i}"))))
                .collect::<Vec<_>>();
            let kdot = |v: &[Atom]| {
                vectors
                    .iter()
                    .enumerate()
                    .fold(Atom::new(), |a, (k, &i)| a + &v[i] * &z[k + 1])
            };
            let mut shifted_coordinates = replacement.clone();
            for (k, &(i, j)) in coordinates.iter().enumerate() {
                let value = if i == loop_index && j == loop_index {
                    &z[0] - Atom::num(2) * kdot(&shift) + dot(&shift, &shift)
                } else if i == loop_index || j == loop_index {
                    let other = if i == loop_index { j } else { i };
                    let mut unit = vec![Atom::new(); l + e];
                    unit[other] = Atom::num(1);
                    kdot(&unit) - dot(&shift, &unit)
                } else {
                    continue;
                };
                shifted_coordinates.insert(x[k].clone(), substitute(&value, &replacement));
            }
            let shifted = expressions
                .iter()
                .map(|a| substitute(a, &shifted_coordinates).expand())
                .collect();
            let tensor_gram = vectors
                .iter()
                .map(|&i| {
                    vectors
                        .iter()
                        .map(|&j| substitute(&gram[i][j], &replacement))
                        .collect()
                })
                .collect();
            let qdot = vectors
                .iter()
                .map(|&i| {
                    let mut unit = vec![Atom::new(); l + e];
                    unit[i] = Atom::num(1);
                    substitute(&dot(&q, &unit), &replacement)
                })
                .collect();
            return Ok(Some(Self {
                lines: [*first, *second],
                dependent,
                shifted,
                z,
                denominators,
                gram: tensor_gram,
                qdot,
                q_squared: substitute(&q_squared, &replacement),
                q_line,
                dimension,
            }));
        }
        Ok(None)
    }

    pub(crate) fn permutation(&self) -> Vec<usize> {
        (0..self.dependent.len())
            .filter(|&i| !self.dependent[i])
            .chain((0..self.dependent.len()).filter(|&i| self.dependent[i]))
            .collect()
    }

    pub(crate) fn reduce(&self, integral: &Integral) -> Result<Option<LinearCombination>> {
        let [first, second] = self.lines;
        if integral.0[first] == 1
            && integral.0[second] == 1
            && integral
                .0
                .iter()
                .zip(&self.dependent)
                .all(|(&n, &d)| !d || n >= 0)
        {
            return Ok(None);
        }
        let rank: i64 = integral
            .0
            .iter()
            .zip(&self.dependent)
            .filter(|(_, d)| **d)
            .map(|(&n, _)| i64::from(n).min(0).abs())
            .sum();
        if rank > 8 || i64::from(integral.0[first]) + i64::from(integral.0[second]) > 32 {
            return Ok(None);
        }
        let mut base = integral.clone();
        base.0[first] = 1;
        base.0[second] = 1;
        let mut numerator = Atom::num(1);
        for (i, (&n, &dependent)) in integral.0.iter().zip(&self.dependent).enumerate() {
            if dependent && n < 0 {
                numerator *= self.shifted[i].clone().pow(-i64::from(n));
                base.0[i] = 0;
            }
        }
        let mut result = LinearCombination::new();
        for (monomial, c) in numerator.expand().coefficient_list::<i32>(&self.z) {
            let exponents = powers(&monomial, &self.z)?;
            let a = i64::from(integral.0[first]) - i64::from(exponents[0]);
            let b = i64::from(integral.0[second]);
            if a <= 0 || b <= 0 {
                continue;
            }
            let indices = exponents[1..]
                .iter()
                .enumerate()
                .flat_map(|(i, &n)| std::iter::repeat_n(i, n as usize))
                .collect::<Vec<_>>();
            for (j, contraction) in contractions(&indices, &self.gram, &self.qdot)
                .into_iter()
                .enumerate()
            {
                let Some(weight) =
                    tensor_coefficient(&self.dimension, a, b, indices.len() as i64, j as i64)
                else {
                    return Ok(None);
                };
                let q_power = j as i64 + 2 - a - b;
                let mut target = base.clone();
                let mut coefficient = &c * weight * contraction;
                if let Some(line) = self.q_line {
                    target.0[line] = i16::try_from(i64::from(target.0[line]) - q_power)
                        .map_err(|_| Error::Limit("bubble index overflow".into()))?;
                } else {
                    coefficient *= self.q_squared.clone().pow(q_power);
                }
                for (monomial, c) in coefficient
                    .expand()
                    .coefficient_list::<i32>(&self.denominators)
                {
                    let mut target = target.clone();
                    for (n, p) in target
                        .0
                        .iter_mut()
                        .zip(powers(&monomial, &self.denominators)?)
                    {
                        *n = n.checked_sub(p).ok_or_else(|| {
                            Error::Limit("bubble numerator index overflow".into())
                        })?;
                    }
                    let previous = result.remove(&target).unwrap_or_default();
                    result.insert(target, clean(previous + c));
                }
            }
        }
        result.retain(|_, c| !c.is_zero());
        Ok(Some(result))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Propagator, ReductionBackend, RunContext, RustRedBackend};

    fn assert_zero(a: Atom) {
        assert!(clean(a.clone()).is_zero(), "nonzero residual: {a}");
    }

    #[test]
    fn tensor_beta_factors_obey_low_rank_identities() {
        let d = parse!("bubble_test_d");
        assert_zero(tensor_coefficient(&d, 2, 1, 0, 0).unwrap() + &d - Atom::num(3));
        assert_zero(tensor_coefficient(&d, 1, 1, 1, 0).unwrap() + Atom::num((1, 2)));
        assert_zero(
            tensor_coefficient(&d, 1, 1, 2, 0).unwrap() - &d / (Atom::num(4) * (&d - Atom::num(1))),
        );
        assert_zero(
            tensor_coefficient(&d, 1, 1, 2, 1).unwrap()
                + Atom::num(1) / (Atom::num(4) * (&d - Atom::num(1))),
        );
        assert_zero(
            tensor_coefficient(&d, 1, 1, 3, 0).unwrap()
                + (&d + Atom::num(2)) / (Atom::num(8) * (&d - Atom::num(1))),
        );
        assert_zero(
            tensor_coefficient(&d, 1, 1, 3, 1).unwrap()
                - Atom::num(1) / (Atom::num(8) * (&d - Atom::num(1))),
        );
    }

    fn compare_native(family: &IntegralFamily, targets: Vec<Integral>) {
        let d = Atom::num(family.dimension) - Atom::num(2) * Atom::var(family.epsilon);
        let pattern = Bubble::find(
            family,
            &targets[0].0.iter().map(|&n| n > 0).collect::<Vec<_>>(),
            d,
        )
        .unwrap()
        .unwrap();
        let mut all = targets.clone();
        let rules = targets
            .iter()
            .map(|i| {
                let terms = pattern.reduce(i).unwrap().unwrap();
                all.extend(terms.keys().cloned());
                terms
            })
            .collect::<Vec<_>>();
        all.sort();
        all.dedup();
        let optimized = RustRedBackend {
            max_sector_batch: 1,
            ..Default::default()
        }
        .reduce(family, &targets, &RunContext::default())
        .unwrap();
        all.extend(optimized.residuals.iter().cloned());
        let reduced = RustRedBackend {
            factorized: false,
            max_targets: 4096,
            ..Default::default()
        }
        .reduce(family, &all, &RunContext::default())
        .unwrap();
        for (target, rhs) in targets.iter().zip(rules) {
            let mut residual = reduced.expand(target).unwrap();
            for (i, coefficient) in rhs {
                for (master, c) in reduced.expand(&i).unwrap() {
                    let old = residual.remove(&master).unwrap_or_default();
                    residual.insert(master, clean(old - &coefficient * c));
                }
            }
            for c in residual.into_values() {
                assert_zero(c);
            }
            let mut residual = reduced.expand(target).unwrap();
            for (i, coefficient) in optimized.expand(target).unwrap() {
                for (master, c) in reduced.expand(&i).unwrap() {
                    let old = residual.remove(&master).unwrap_or_default();
                    residual.insert(master, clean(old - &coefficient * c));
                }
            }
            for c in residual.into_values() {
                assert_zero(c);
            }
        }
    }

    #[test]
    fn external_bubble_tensor_identities_agree_with_native_ibp() {
        let gram = vec![
            vec![Atom::num(-2), Atom::num(1)],
            vec![Atom::num(1), Atom::num(-3)],
        ];
        let family = IntegralFamily {
            name: "tensor_bubble_check".into(),
            loops: vec!["l".into()],
            external: vec!["p".into(), "v".into()],
            propagators: vec![
                Propagator::quadratic(&[1], &[0, 0], Atom::new(), &gram).unwrap(),
                Propagator::quadratic(&[1], &[1, 0], Atom::new(), &gram).unwrap(),
                Propagator {
                    constant: Atom::new(),
                    scalar_products: vec![Atom::new(), Atom::new(), Atom::num(1)],
                },
            ],
            external_gram: gram,
            physical_propagators: 2,
            epsilon: symbol!("bubble_check_eps"),
            dimension: 4,
        };
        compare_native(
            &family,
            vec![
                Integral(vec![2, 1, 0]),
                Integral(vec![1, 2, -1]),
                Integral(vec![1, 1, -2]),
                Integral(vec![1, 1, -3]),
            ],
        );
    }

    #[test]
    fn bubble_insertion_tensor_identities_agree_with_native_ibp() {
        let gram = vec![vec![Atom::num(-1)]];
        let family = IntegralFamily {
            name: "insertion_bubble_check".into(),
            loops: vec!["l1".into(), "l2".into()],
            external: vec!["p".into()],
            propagators: vec![
                Propagator::quadratic(&[1, 0], &[0], Atom::new(), &gram).unwrap(),
                Propagator::quadratic(&[1, 1], &[0], Atom::new(), &gram).unwrap(),
                Propagator::quadratic(&[0, 1], &[0], Atom::new(), &gram).unwrap(),
                Propagator::quadratic(&[0, 1], &[1], Atom::num(1), &gram).unwrap(),
                Propagator::quadratic(&[1, 0], &[1], Atom::new(), &gram).unwrap(),
            ],
            external_gram: gram,
            physical_propagators: 4,
            epsilon: symbol!("insertion_check_eps"),
            dimension: 4,
        };
        compare_native(
            &family,
            vec![
                Integral(vec![2, 1, 1, 1, 0]),
                Integral(vec![1, 2, 1, 1, -1]),
                Integral(vec![1, 1, 1, 1, -2]),
            ],
        );
    }

    #[test]
    fn paper_bubble_subloop_rewrites_follow_native_order() {
        let (family, _) = crate::benchmarks::paper_two_loop().unwrap();
        let (family, _) = family
            .deform(symbol!("bubble_paper_eta"), &crate::MassMode::Auto)
            .unwrap();
        let sector = [false, false, true, true, true, true, true, false, false];
        let pattern = Bubble::find(
            &family,
            &sector,
            Atom::num(4) - Atom::num(2) * Atom::var(family.epsilon),
        )
        .unwrap()
        .unwrap();
        let order = rustred::solver::IntegralOrder::new(sector, [false; 9])
            .with_permutation(pattern.permutation().try_into().unwrap())
            .unwrap();
        for a in 1..=3 {
            for b in 1..=3 {
                for rank in 0..=3 {
                    let indices = [0, 0, a, 1, 1, 1, b, -rank, 0];
                    if a == 1 && b == 1 && rank == 0 {
                        continue;
                    }
                    let lhs = rustred::solver::Integral::numeric(indices).unwrap();
                    let terms = pattern
                        .reduce(&Integral(indices.to_vec()))
                        .unwrap()
                        .unwrap();
                    for rhs in terms.keys() {
                        let rhs = rustred::solver::Integral::numeric(
                            rhs.0.as_slice().try_into().unwrap(),
                        )
                        .unwrap();
                        assert_eq!(
                            order.compare(&lhs, &rhs),
                            std::cmp::Ordering::Less,
                            "{lhs} -> {rhs}"
                        );
                    }
                }
            }
        }
    }
}
