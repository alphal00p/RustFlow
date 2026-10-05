#![allow(clippy::needless_range_loop)] // Symmetric quadratic forms and contractions use explicit paired indices.
//! Gaussian FT terminals, including polynomial tensor numerators.
use crate::algebra::{determinant, inverse};
use crate::*;
use symbolica::prelude::*;

/// Integrate a family whose target has exactly one positive denominator power.
/// The positive denominator may be a sum of ordinary propagators and hence
/// have a full-rank quadratic form in several loop momenta.
pub fn terminal(
    family: &IntegralFamily,
    integral: &Integral,
    epsilon: &Rational,
    p: Precision,
) -> Result<ComplexFloat> {
    let positive = integral
        .0
        .iter()
        .enumerate()
        .filter_map(|(i, &n)| (n > 0).then_some(i))
        .collect::<Vec<_>>();
    if positive.len() != 1 {
        return Err(Error::InvalidInput(
            "Gaussian terminal needs exactly one positive denominator".into(),
        ));
    }
    let index = positive[0];
    let power = i64::from(integral.0[index]);
    let denominator = &family.propagators[index];
    let l = family.loops.len();
    let e = family.external.len();
    let mut q = vec![vec![Atom::new(); l]; l];
    let mut k = 0;
    for i in 0..l {
        for j in i..l {
            let a = &denominator.scalar_products[k] / Atom::num(if i == j { 1 } else { 2 });
            q[i][j] = a.clone();
            q[j][i] = a;
            k += 1;
        }
    }
    let det = determinant(q.clone());
    if det.is_zero() {
        return Ok(p.zero());
    }
    let inverse = inverse(&q)?;
    let b = (0..l)
        .map(|i| {
            (0..e)
                .map(|a| &denominator.scalar_products[k + i * e + a] / Atom::num(2))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let shifts = (0..l)
        .map(|i| {
            (0..e)
                .map(|a| -(0..l).fold(Atom::new(), |s, j| s + &inverse[i][j] * &b[j][a]))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut mass = -denominator.constant.clone();
    for i in 0..l {
        for j in 0..l {
            for a in 0..e {
                for c in 0..e {
                    mass += &b[i][a] * &inverse[i][j] * &b[j][c] * &family.external_gram[a][c];
                }
            }
        }
    }
    mass = mass.together().cancel();
    if mass.is_zero() {
        return Ok(p.zero());
    }
    let mut labels = Vec::new();
    for i in 0..l {
        for j in i..l {
            labels.push((i, j));
        }
    }
    for i in 0..l {
        for a in 0..e {
            labels.push((i, l + a));
        }
    }
    let variables = (0..labels.len())
        .map(|i| Atom::var(symbol!(format!("symbolica_amflow::gaussian_{i}"))))
        .collect::<Vec<_>>();
    let scalar = |i, j| {
        variables[labels
            .iter()
            .position(|&(a, b)| (a == i && b == j) || (a == j && b == i))
            .unwrap()]
        .clone()
    };
    let transformed = labels
        .iter()
        .enumerate()
        .map(|(k, &(i, j))| {
            let mut v = variables[k].clone();
            if j < l {
                for a in 0..e {
                    v += &shifts[i][a] * scalar(j, l + a) + &shifts[j][a] * scalar(i, l + a);
                    for b in 0..e {
                        v += &shifts[i][a] * &shifts[j][b] * &family.external_gram[a][b];
                    }
                }
            } else {
                for a in 0..e {
                    v += &shifts[i][a] * &family.external_gram[a][j - l];
                }
            }
            v
        })
        .collect::<Vec<_>>();
    let mut numerator = Atom::num(1);
    for (d, &n) in family.propagators.iter().zip(&integral.0) {
        if n < 0 {
            let polynomial = d
                .scalar_products
                .iter()
                .zip(&transformed)
                .fold(d.constant.clone(), |a, (b, c)| a + b * c);
            numerator *= polynomial.pow(-i64::from(n));
        }
    }
    let dimension = Atom::num(family.dimension) - Atom::num(2) * Atom::var(family.epsilon);
    let values = ahash::HashMap::from_iter([(Atom::var(family.epsilon), p.rational(epsilon))]);
    let half = Rational::from((family.dimension, 2)) - epsilon;
    let mut answer = p.zero();
    for (monomial, coefficient) in
        crate::coefficient::exact_coefficient_list(&numerator.expand(), &variables)?
    {
        let mut loops = Vec::new();
        let mut links = Vec::new();
        let mut external = Vec::new();
        for ((i, j), power) in labels
            .iter()
            .copied()
            .zip(crate::integrand::powers(&monomial, &variables)?)
        {
            for _ in 0..power {
                let a = loops.len();
                loops.push(i);
                external.push(if j < l { None } else { Some(j - l) });
                if j < l {
                    loops.push(j);
                    external.push(None);
                    links.push((a, a + 1));
                }
            }
        }
        if loops.len() % 2 != 0 {
            continue;
        }
        if loops.len() > 8 {
            return Err(Error::Limit(
                "Gaussian numerator tensor rank exceeds eight".into(),
            ));
        }
        let r = loops.len() / 2;
        let moment = wick(
            &loops,
            &links,
            &external,
            &inverse,
            &family.external_gram,
            &dimension,
        )?;
        let gamma_argument = Rational::from(power - r as i64) - &half * &Rational::from(l as i64);
        let gamma = p.gamma_real(&p.rational(&gamma_argument).re)?;
        let mass_power = p.pow(&p.eval(&mass, &values)?, &p.rational(&(-gamma_argument)));
        answer = p.add(
            &answer,
            &p.mul(
                &p.eval(&(coefficient * moment), &values)?,
                &p.mul(&gamma, &mass_power),
            ),
        );
    }
    let det_power = p.pow(&p.eval(&det, &values)?, &p.rational(&(-half)));
    answer = p.div(&p.mul(&answer, &det_power), &p.gamma_real(&p.real(power))?);
    Ok(if power % 2 == 0 {
        answer
    } else {
        p.neg(&answer)
    })
}

fn wick(
    loops: &[usize],
    links: &[(usize, usize)],
    external: &[Option<usize>],
    inverse: &[Vec<Atom>],
    gram: &[Vec<Atom>],
    dimension: &Atom,
) -> Result<Atom> {
    fn pairings(indices: &[usize]) -> Vec<Vec<(usize, usize)>> {
        if indices.is_empty() {
            return vec![vec![]];
        }
        let mut out = Vec::new();
        for j in 1..indices.len() {
            let rest = indices[1..]
                .iter()
                .copied()
                .filter(|&i| i != indices[j])
                .collect::<Vec<_>>();
            for mut pairs in pairings(&rest) {
                pairs.push((indices[0], indices[j]));
                out.push(pairs);
            }
        }
        out
    }
    let mut sum = Atom::new();
    let n = loops.len();
    for pairs in pairings(&(0..n).collect::<Vec<_>>()) {
        let mut components = linnet::union_find::UnionFind::new(
            external
                .iter()
                .map(|end| end.iter().copied().collect::<Vec<_>>())
                .collect(),
        );
        for &(i, j) in links.iter().chain(&pairs) {
            components.union(
                linnet::half_edge::involution::Hedge(i),
                linnet::half_edge::involution::Hedge(j),
                |mut left, mut right| {
                    left.append(&mut right);
                    left
                },
            );
        }
        let mut term = pairs.iter().fold(Atom::num(1), |a, &(i, j)| {
            -a * &inverse[loops[i]][loops[j]] / Atom::num(2)
        });
        for (_, ends) in components.iter_set_data() {
            match ends.as_slice() {
                [] => term *= dimension,
                [a, b] => term *= &gram[*a][*b],
                _ => return Err(Error::Numerical("invalid Wick contraction graph".into())),
            }
        }
        sum += term;
    }
    Ok(sum.together().cancel())
}
