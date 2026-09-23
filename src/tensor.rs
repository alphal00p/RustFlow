//! Exact isotropic tensor projection for vacuum boundary numerators.
use crate::algebra::inverse;
use crate::{Error, Result};
use std::collections::BTreeMap;
use symbolica::prelude::*;

type Pairing = Vec<(usize, usize)>;
type PairingOrbits = Vec<Vec<Pairing>>;
fn pairings(indices: &[usize]) -> Vec<Pairing> {
    if indices.is_empty() {
        return vec![Vec::new()];
    }
    let first = indices[0];
    let mut out = Vec::new();
    for j in 1..indices.len() {
        let remaining = indices[1..]
            .iter()
            .copied()
            .filter(|&i| i != indices[j])
            .collect::<Vec<_>>();
        for mut pairing in pairings(&remaining) {
            pairing.push((first, indices[j]));
            out.push(pairing);
        }
    }
    out
}
fn cycles(a: &Pairing, b: &Pairing, rank: usize) -> usize {
    let mut parent = (0..rank).collect::<Vec<_>>();
    fn root(parent: &[usize], mut i: usize) -> usize {
        while parent[i] != i {
            i = parent[i];
        }
        i
    }
    for &(i, j) in a.iter().chain(b) {
        let l = root(&parent, i);
        let r = root(&parent, j);
        parent[l] = r;
    }
    (0..rank).filter(|&i| root(&parent, i) == i).count()
}

fn vector_classes(gram: &[Vec<Atom>]) -> (Vec<usize>, Vec<usize>) {
    let mut representatives = Vec::<usize>::new();
    let classes = (0..gram.len())
        .map(|i| {
            if let Some(class) = representatives.iter().position(|&j| gram[i] == gram[j]) {
                class
            } else {
                representatives.push(i);
                representatives.len() - 1
            }
        })
        .collect();
    (classes, representatives)
}

// Wick's pairing sum, memoized by multiplicities of equal Gram rows. This
// avoids materializing (rank-1)!! labelled pairings for a repeated hard vector.
fn pairing_sum(gram: &[Vec<Atom>]) -> Result<Atom> {
    let (classes, representatives) = vector_classes(gram);
    let mut counts = vec![0; representatives.len()];
    for class in classes {
        counts[class] += 1;
    }
    fn sum(
        gram: &[Vec<Atom>],
        representatives: &[usize],
        counts: &mut [usize],
        memo: &mut BTreeMap<Vec<usize>, Atom>,
    ) -> Result<Atom> {
        if let Some(value) = memo.get(counts) {
            return Ok(value.clone());
        }
        let Some(first) = counts.iter().position(|&n| n != 0) else {
            return Ok(Atom::num(1));
        };
        if memo.len() >= 100_000 {
            return Err(Error::Limit(
                "tensor pairing sum exceeds 100000 states".into(),
            ));
        }
        let key = counts.to_vec();
        counts[first] -= 1;
        let mut value = Atom::new();
        for second in 0..counts.len() {
            let choices = counts[second];
            if choices == 0 || gram[representatives[first]][representatives[second]].is_zero() {
                continue;
            }
            counts[second] -= 1;
            value += Atom::num(choices as i64)
                * &gram[representatives[first]][representatives[second]]
                * sum(gram, representatives, counts, memo)?;
            counts[second] += 1;
        }
        counts[first] += 1;
        memo.insert(key, value.clone());
        Ok(value)
    }
    sum(gram, &representatives, &mut counts, &mut BTreeMap::new())
}

/// Project a product of scalar products (hard_vector[i] . external_vector[i]).
/// `hard_gram` and `external_gram` give all pairwise scalar products of the
/// labelled vectors; repeated vectors may appear as repeated Gram rows.
/// The returned expression is valid inside a rotationally invariant vacuum
/// integral. Odd rank vanishes. Repeated single-vector tensors support rank 32;
/// other tensors support rank 14 with at most 128 pairing orbits.
#[derive(Clone, Debug)]
pub struct TensorProjector {
    dimension: Atom,
    inverses: BTreeMap<Vec<usize>, (PairingOrbits, Vec<Vec<Atom>>)>,
}
impl TensorProjector {
    pub fn new(dimension: Atom) -> Self {
        Self {
            dimension,
            inverses: BTreeMap::new(),
        }
    }
    pub fn project(
        &mut self,
        hard_gram: &[Vec<Atom>],
        external_gram: &[Vec<Atom>],
    ) -> Result<Atom> {
        let rank = hard_gram.len();
        if external_gram.len() != rank
            || hard_gram
                .iter()
                .chain(external_gram)
                .any(|r| r.len() != rank)
        {
            return Err(Error::InvalidInput("tensor Gram dimensions".into()));
        }
        if rank % 2 == 1 {
            return Ok(Atom::new());
        }
        if rank == 0 {
            return Ok(Atom::num(1));
        }
        if rank > 32 {
            return Err(Error::Limit(format!(
                "vacuum tensor rank {rank} exceeds 32"
            )));
        }
        if (0..rank).any(|i| {
            (0..i).any(|j| {
                hard_gram[i][j] != hard_gram[j][i] || external_gram[i][j] != external_gram[j][i]
            })
        }) {
            return Err(Error::InvalidInput(
                "tensor Gram matrices must be symmetric".into(),
            ));
        }
        let (classes, representatives) = vector_classes(hard_gram);
        let class_count = representatives.len();
        if class_count == 1 {
            let denominator = (0..rank / 2)
                .fold(Atom::num(1), |value, k| {
                    value * (&self.dimension + Atom::num(2 * k as i64))
                })
                .together()
                .cancel();
            if denominator.is_zero() {
                return Err(Error::Numerical(
                    "singular isotropic tensor dimension".into(),
                ));
            }
            return Ok((hard_gram[0][0].clone().pow((rank / 2) as i64)
                * pairing_sum(external_gram)?
                / denominator)
                .together()
                .cancel());
        }
        if rank > 14 {
            return Err(Error::Limit(format!(
                "tensor rank {rank} with {class_count} distinct hard vectors exceeds the rank-14 orbit budget"
            )));
        }
        if !self.inverses.contains_key(&classes) {
            let pairings = pairings(&(0..rank).collect::<Vec<_>>());
            // Permuting identical hard vectors leaves the contraction vector
            // invariant. Pairings are in the same orbit exactly when they have
            // the same numbers of edges between vector classes. The inverse
            // problem therefore closes on this much smaller invariant space.
            let mut groups = BTreeMap::<Vec<usize>, Vec<Pairing>>::new();
            for pairing in pairings {
                let mut signature = vec![0; class_count * class_count];
                for &(i, j) in &pairing {
                    let (a, b) = (classes[i].min(classes[j]), classes[i].max(classes[j]));
                    signature[a * class_count + b] += 1;
                }
                groups.entry(signature).or_default().push(pairing);
                if groups.len() > 128 {
                    return Err(Error::Limit(
                        "tensor invariant basis exceeds 128 pairing orbits".into(),
                    ));
                }
            }
            let orbits = groups.into_values().collect::<PairingOrbits>();
            let gram = orbits
                .iter()
                .map(|a| {
                    orbits
                        .iter()
                        .map(|b| {
                            b.iter().fold(Atom::new(), |sum, pairing| {
                                sum + self
                                    .dimension
                                    .clone()
                                    .pow(cycles(&a[0], pairing, rank) as i64)
                            })
                        })
                        .collect()
                })
                .collect::<Vec<Vec<_>>>();
            self.inverses
                .insert(classes.clone(), (orbits, inverse(&gram)?));
        }
        let (orbits, inverse) = &self.inverses[&classes];
        let contract = |gram: &[Vec<Atom>], pairing: &Pairing| {
            pairing
                .iter()
                .fold(Atom::num(1), |a, &(i, j)| a * &gram[i][j])
        };
        let hard = orbits
            .iter()
            .map(|orbit| contract(hard_gram, &orbit[0]))
            .collect::<Vec<_>>();
        let external = orbits
            .iter()
            .map(|orbit| {
                orbit.iter().fold(Atom::new(), |sum, pairing| {
                    sum + contract(external_gram, pairing)
                })
            })
            .collect::<Vec<_>>();
        let mut result = Atom::new();
        for (i, row) in inverse.iter().enumerate() {
            for (j, c) in row.iter().enumerate() {
                result += c * &hard[j] * &external[i];
            }
        }
        Ok(result.together().cancel())
    }
}
