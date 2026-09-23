//! Exact isotropic tensor projection for vacuum boundary numerators.
use crate::algebra::inverse;
use crate::{Error, Result};
use std::collections::BTreeMap;
use symbolica::prelude::*;

type Pairing = Vec<(usize, usize)>;
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

/// Project a product of scalar products (hard_vector[i] . external_vector[i]).
/// `hard_gram` and `external_gram` give all pairwise scalar products of the
/// labelled vectors; repeated vectors may appear as repeated Gram rows.
/// The returned expression is valid inside a rotationally invariant vacuum
/// integral. Odd rank vanishes. The practical exact-algebra rank cap is eight.
#[derive(Clone, Debug)]
pub struct TensorProjector {
    dimension: Atom,
    inverses: BTreeMap<usize, (Vec<Pairing>, Vec<Vec<Atom>>)>,
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
        if rank > 8 {
            return Err(Error::Limit("vacuum tensor rank exceeds eight".into()));
        }
        if !self.inverses.contains_key(&rank) {
            let pairings = pairings(&(0..rank).collect::<Vec<_>>());
            let gram = pairings
                .iter()
                .map(|a| {
                    pairings
                        .iter()
                        .map(|b| self.dimension.clone().pow(cycles(a, b, rank) as i64))
                        .collect()
                })
                .collect::<Vec<Vec<_>>>();
            self.inverses.insert(rank, (pairings, inverse(&gram)?));
        }
        let (pairings, inverse) = &self.inverses[&rank];
        let contract = |gram: &[Vec<Atom>], pairing: &Pairing| {
            pairing
                .iter()
                .fold(Atom::num(1), |a, &(i, j)| a * &gram[i][j])
        };
        let hard = pairings
            .iter()
            .map(|pairing| contract(hard_gram, pairing))
            .collect::<Vec<_>>();
        let external = pairings
            .iter()
            .map(|pairing| contract(external_gram, pairing))
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
