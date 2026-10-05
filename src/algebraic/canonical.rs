//! Canonical algebraic differential forms without dense physical assembly.
use super::{
    AlgebraicKinematicSystem, AlgebraicSystem, SquareRoot, decoded, normalized_polynomial,
    quotient, rational,
};
use crate::diffexp::EpsilonSystem;
use crate::family::{imaginary_parameter, substitute};
use crate::kinematics::KinematicPath;
use crate::{DifferentialSystem, Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use symbolica::prelude::*;

/// A canonical form `dY = epsilon sum M_a dlog(L_a) Y` in physical invariants.
///
/// The exact letters, constant matrices and named square-root generators remain
/// separate until a path is supplied. Original letter/root domains are retained
/// even when a coefficient, Jacobian or sum of letters later cancels. Formal
/// root norms require invertibility on every sheet, and can therefore exclude a
/// point that is regular on a chosen sheet. No sheet-specific domain is inferred.
#[derive(Clone, Debug)]
pub struct CanonicalAlgebraicSystem {
    data: Arc<CanonicalData>,
}

#[derive(Debug)]
struct CanonicalData {
    epsilon: Symbol,
    variables: Vec<Symbol>,
    letters: Vec<Atom>,
    matrices: Vec<Vec<Vec<Atom>>>,
    roots: Vec<SquareRoot>,
    nonzero_conditions: Vec<Atom>,
}

struct RestrictedPath {
    rules: BTreeMap<Atom, Atom>,
    roots: Vec<SquareRoot>,
    guards: Vec<Atom>,
    allowed: BTreeSet<Atom>,
}

impl CanonicalAlgebraicSystem {
    pub fn new(
        epsilon: Symbol,
        variables: &[Symbol],
        letters: &[Atom],
        matrices: &[Vec<Vec<Atom>>],
        roots: Vec<SquareRoot>,
    ) -> Result<Self> {
        if variables.is_empty()
            || variables.iter().copied().collect::<BTreeSet<_>>().len() != variables.len()
            || letters.is_empty()
            || letters.len() != matrices.len()
        {
            return Err(Error::InvalidInput(
                "canonical data needs distinct variables and equal nonempty letters/matrices"
                    .into(),
            ));
        }
        let n = matrices[0].len();
        // Reuse the existing dimension and root-registry validation, without
        // constructing any physical logarithmic derivative or matrix sum.
        AlgebraicKinematicSystem {
            epsilon,
            derivatives: variables
                .iter()
                .map(|&v| (v, vec![vec![Atom::new(); n]; n]))
                .collect(),
            roots: roots.clone(),
        }
        .validate()?;
        let physical = variables
            .iter()
            .copied()
            .chain([imaginary_parameter()])
            .map(Atom::var)
            .collect::<BTreeSet<_>>();
        let mut allowed = physical.clone();
        allowed.extend(roots.iter().map(|r| Atom::var(r.symbol)));
        let constants = BTreeSet::from([Atom::var(imaginary_parameter())]);
        let root_atoms = roots
            .iter()
            .map(|r| Atom::var(r.symbol))
            .collect::<Vec<_>>();
        let mut guards = Vec::new();
        let physical_symbols = variables
            .iter()
            .copied()
            .chain([imaginary_parameter()])
            .collect::<BTreeSet<_>>();
        let algebraic_symbols = physical_symbols
            .iter()
            .copied()
            .chain(roots.iter().map(|r| r.symbol))
            .collect::<BTreeSet<_>>();
        for root in &roots {
            for condition in crate::physical_conditions::canonical_conditions(
                std::slice::from_ref(&root.radicand),
                &physical_symbols,
            )? {
                add_condition_parts(&condition, &physical, &mut guards)?;
            }
        }
        for (letter, matrix) in letters.iter().zip(matrices) {
            if matrix.len() != n || matrix.iter().any(|row| row.len() != n) {
                return Err(Error::InvalidInput(
                    "canonical matrices must have identical square dimensions".into(),
                ));
            }
            for coefficient in matrix.iter().flatten() {
                rational(coefficient, &constants)?;
            }
            rational(letter, &allowed)?;
            // Preserve raw negative-power bases inside each letter before native
            // rational conversion can cancel a removable source denominator.
            let conditions = crate::physical_conditions::canonical_conditions(
                std::slice::from_ref(letter),
                &algebraic_symbols,
            )?;
            for part in conditions {
                let normalized = normalized_polynomial(&part, &root_atoms, &roots)?;
                if normalized.is_empty() {
                    return Err(Error::InvalidInput(
                        "canonical letter vanishes or has zero denominator under root relations"
                            .into(),
                    ));
                }
                let inverse = quotient::invert_denominator(
                    &quotient::recompose(&normalized, &roots),
                    &roots,
                    &allowed,
                )?;
                for coefficient in inverse.values() {
                    let fraction = rational(coefficient, &allowed)?;
                    add_condition_parts(
                        &fraction.denominator.to_expression(),
                        &physical,
                        &mut guards,
                    )?;
                }
            }
        }
        canonicalize_conditions(&mut guards);
        Ok(Self {
            data: Arc::new(CanonicalData {
                epsilon,
                variables: variables.to_vec(),
                letters: letters.to_vec(),
                matrices: matrices.to_vec(),
                roots,
                nonzero_conditions: guards,
            }),
        })
    }

    pub fn epsilon(&self) -> Symbol {
        self.data.epsilon
    }
    pub fn variables(&self) -> &[Symbol] {
        &self.data.variables
    }
    pub fn letters(&self) -> &[Atom] {
        &self.data.letters
    }
    pub fn constant_matrices(&self) -> &[Vec<Vec<Atom>>] {
        &self.data.matrices
    }
    pub fn roots(&self) -> &[SquareRoot] {
        &self.data.roots
    }
    pub fn nonzero_conditions(&self) -> &[Atom] {
        &self.data.nonzero_conditions
    }
    pub fn dimension(&self) -> usize {
        self.data.matrices[0].len()
    }

    /// Restrict each exact letter before assembling a dense derivative matrix.
    /// All original domains are restricted first, including stationary variables.
    pub fn pullback(&self, path: &KinematicPath, order: usize) -> Result<AlgebraicSystem> {
        let restricted = self.restrict_path(path)?;
        let allowed = self.original_symbols();
        let letters = self
            .letters()
            .iter()
            .map(|letter| {
                restrict_rational(letter, &allowed, &restricted.allowed, &restricted.rules)
            })
            .collect::<Result<Vec<_>>>()?;
        // Reuse the existing total logarithmic derivative and validation. It now
        // operates only in one variable; no six-variable dense GCD is needed.
        let mut univariate = AlgebraicKinematicSystem::canonical_dlog(
            self.epsilon(),
            &[path.parameter],
            &letters,
            self.constant_matrices(),
            restricted.roots.clone(),
        )?;
        let matrix = univariate
            .derivatives
            .remove(&path.parameter)
            .ok_or_else(|| Error::InvalidInput("canonical path derivative is missing".into()))?;
        let result = AlgebraicSystem {
            system: EpsilonSystem::from_differential_system(
                &DifferentialSystem {
                    variable: path.parameter,
                    matrix,
                },
                self.epsilon(),
                order,
            )?,
            roots: restricted.roots,
            nonzero_conditions: restricted.guards,
        };
        result.validate()?;
        Ok(result)
    }

    fn original_symbols(&self) -> BTreeSet<Atom> {
        self.variables()
            .iter()
            .copied()
            .chain([imaginary_parameter()])
            .chain(self.roots().iter().map(|r| r.symbol))
            .map(Atom::var)
            .collect()
    }

    fn restrict_path(&self, path: &KinematicPath) -> Result<RestrictedPath> {
        path.validate()?;
        let variables = self.variables().iter().copied().collect::<BTreeSet<_>>();
        if path.parameter == self.epsilon()
            || path.parameter == imaginary_parameter()
            || variables.contains(&path.parameter)
            || self.roots().iter().any(|r| r.symbol == path.parameter)
            || variables != path.coordinates.keys().copied().collect::<BTreeSet<_>>()
        {
            return Err(Error::InvalidInput(
                "canonical path needs exactly the physical variables and a distinct parameter"
                    .into(),
            ));
        }
        let path_symbols =
            BTreeSet::from([Atom::var(path.parameter), Atom::var(imaginary_parameter())]);
        let mut allowed = path_symbols.clone();
        allowed.extend(self.roots().iter().map(|r| Atom::var(r.symbol)));
        let original = self.original_symbols();
        let rules = path
            .coordinates
            .iter()
            .map(|(&v, a)| (Atom::var(v), a.clone()))
            .collect::<BTreeMap<_, _>>();
        let mut guards = Vec::new();
        for coordinate in path.coordinates.values() {
            rational(coordinate, &path_symbols)?;
            // A coordinate may be zero; only its original denominator bases
            // restrict the path domain, including removable coordinate poles.
            for condition in crate::physical_conditions::rational_denominator_conditions(
                std::slice::from_ref(coordinate),
                &BTreeSet::from([path.parameter, imaginary_parameter()]),
            )? {
                add_condition_parts(&condition, &path_symbols, &mut guards)?;
            }
        }
        for guard in self.nonzero_conditions() {
            let value = restrict_rational(guard, &original, &allowed, &rules)?;
            add_condition_parts(&value, &path_symbols, &mut guards)?;
        }
        let roots = self
            .roots()
            .iter()
            .map(|r| {
                let radicand = restrict_rational(&r.radicand, &original, &allowed, &rules)?;
                add_condition_parts(&radicand, &path_symbols, &mut guards)?;
                Ok(SquareRoot {
                    symbol: r.symbol,
                    radicand,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        canonicalize_conditions(&mut guards);
        Ok(RestrictedPath {
            rules,
            roots,
            guards,
            allowed,
        })
    }
}

fn restrict_rational(
    expression: &Atom,
    original: &BTreeSet<Atom>,
    restricted: &BTreeSet<Atom>,
    rules: &BTreeMap<Atom, Atom>,
) -> Result<Atom> {
    let fraction = rational(expression, original)?;
    let denominator = decoded(&substitute(&fraction.denominator.to_expression(), rules));
    rational(&denominator, restricted)?;
    if denominator.is_zero() {
        return Err(Error::InvalidInput(
            "path lies identically on a source denominator".into(),
        ));
    }
    let numerator = decoded(&substitute(&fraction.numerator.to_expression(), rules));
    rational(&numerator, restricted)?;
    Ok(decoded(&(numerator / denominator)))
}

/// Record numerator and denominator individually. Cancellation against another
/// source condition or another logarithmic letter must never erase a hole.
fn add_condition_parts(
    expression: &Atom,
    allowed: &BTreeSet<Atom>,
    conditions: &mut Vec<Atom>,
) -> Result<()> {
    let fraction = rational(expression, allowed)?;
    for part in [
        fraction.numerator.to_expression(),
        fraction.denominator.to_expression(),
    ] {
        let value = decoded(&part);
        if value.is_zero() {
            return Err(Error::InvalidInput(
                "canonical source domain is identically singular".into(),
            ));
        }
        if !matches!(value.as_view(), AtomView::Num(_)) {
            conditions.push(value);
        }
    }
    Ok(())
}
fn canonicalize_conditions(conditions: &mut Vec<Atom>) {
    conditions.sort();
    conditions.dedup();
}
