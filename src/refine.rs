//! Optional D-factorizing basis refinement, following AMFlow 2.0's row swaps.
use crate::algebra::{inverse, matmul};
use crate::family::scalar_symbols;
use crate::reduction::{LinearCombination, ReducedSystem};
use crate::{DifferentialSystem, Error, Integral, Result};
use std::collections::BTreeSet;
use symbolica::prelude::*;

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct RefinementReport {
    pub basis_changes: usize,
    pub factorized: bool,
}

fn mixed_denominator(a: &Atom, epsilon: Symbol) -> Result<bool> {
    let rational: RationalPolynomial<IntegerRing, u16> = a
        .try_to_rational_polynomial(&Q, &Z, None)
        .map_err(|e| Error::InvalidInput(e.to_string()))?;
    let denominator = rational.denominator.to_expression().factor();
    let factors = if let AtomView::Mul(m) = denominator.as_view() {
        m.iter().map(|v| v.to_owned()).collect::<Vec<_>>()
    } else {
        vec![denominator]
    };
    for factor in factors {
        let mut symbols = BTreeSet::new();
        scalar_symbols(factor.as_view(), &mut symbols)?;
        if symbols.contains(&Atom::var(epsilon)) && symbols.len() > 1 {
            return Ok(true);
        }
    }
    Ok(false)
}

fn vector(terms: &LinearCombination, basis: &[Integral]) -> Vec<Atom> {
    basis
        .iter()
        .map(|i| terms.get(i).cloned().unwrap_or_default())
        .collect()
}
fn combination(row: Vec<Atom>, basis: &[Integral]) -> LinearCombination {
    basis
        .iter()
        .cloned()
        .zip(row)
        .filter(|(_, c)| !c.is_zero())
        .collect()
}

pub fn refine_basis(
    system: &mut ReducedSystem,
    variable: Symbol,
    epsilon: Symbol,
    max_changes: usize,
) -> Result<RefinementReport> {
    let n = system.basis.len();
    if n == 0 {
        return Ok(RefinementReport {
            basis_changes: 0,
            factorized: true,
        });
    }
    let mut visited = BTreeSet::new();
    let mut changes = 0;
    loop {
        if !visited.insert(system.basis.clone()) {
            return Ok(RefinementReport {
                basis_changes: changes,
                factorized: false,
            });
        }
        let mut chosen = None;
        'candidate: for (target, terms) in &system.candidates {
            if system.basis.contains(target) {
                continue;
            }
            for (column, integral) in system.basis.iter().enumerate() {
                if let Some(coefficient) = terms.get(integral)
                    && mixed_denominator(coefficient, epsilon)?
                {
                    chosen = Some((target.clone(), column, vector(terms, &system.basis)));
                    break 'candidate;
                }
            }
        }
        let Some((target, column, row)) = chosen else {
            let mut factorized = true;
            for terms in system.candidates.values() {
                for coefficient in terms.values() {
                    factorized &= !mixed_denominator(coefficient, epsilon)?;
                }
            }
            return Ok(RefinementReport {
                basis_changes: changes,
                factorized,
            });
        };
        if changes >= max_changes {
            return Ok(RefinementReport {
                basis_changes: changes,
                factorized: false,
            });
        }
        let mut forward = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| Atom::num(i64::from(i == j)))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        forward[column] = row;
        system
            .nonzero_conditions
            .push(forward[column][column].clone());
        let backward = inverse(&forward)?;
        let de = DifferentialSystem {
            variable,
            matrix: system.matrix.clone(),
        }
        .change_basis(&backward)?;
        let old_basis = system.basis.clone();
        system.basis[column] = target;
        system.matrix = de.matrix;
        system
            .transformations
            .push(crate::reduction::BasisTransformation {
                previous: old_basis.clone(),
                current: system.basis.clone(),
                matrix: backward.clone(),
            });
        for terms in &mut system.targets {
            *terms = combination(
                matmul(&[vector(terms, &old_basis)], &backward).remove(0),
                &system.basis,
            );
        }
        for terms in system.candidates.values_mut() {
            *terms = combination(
                matmul(&[vector(terms, &old_basis)], &backward).remove(0),
                &system.basis,
            );
        }
        changes += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;
    #[test]
    fn factorizing_swap_updates_derivative_and_targets() {
        let old = Integral(vec![1]);
        let new = Integral(vec![2]);
        let mut system = ReducedSystem {
            transformations: vec![],
            basis: vec![old.clone()],
            matrix: vec![vec![Atom::new()]],
            targets: vec![BTreeMap::from([(old.clone(), parse!("1/(eps+x)"))])],
            nonzero_conditions: vec![],
            candidates: BTreeMap::from([(
                new.clone(),
                BTreeMap::from([(old, parse!("1/(eps+x)"))]),
            )]),
        };
        let report = refine_basis(&mut system, symbol!("x"), symbol!("eps"), 10).unwrap();
        assert!(report.factorized);
        assert_eq!(report.basis_changes, 1);
        assert_eq!(system.transformations.len(), 1);
        assert!(
            (&system.transformations[0].matrix[0][0] - parse!("eps+x"))
                .expand()
                .is_zero()
        );
        assert_eq!(system.basis, vec![new.clone()]);
        assert_eq!(system.targets[0][&new], Atom::num(1));
        assert!(
            (&system.matrix[0][0] + parse!("1/(eps+x)"))
                .together()
                .cancel()
                .is_zero()
        );
    }
}
