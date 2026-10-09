use crate::algebra::{inverse, rref};
use crate::coefficient::{exact_coefficient_list, powers};
use crate::family::{LinearCombination, substitute};
use crate::{Error, Integral, Result};
use std::collections::BTreeMap;
use symbolica::prelude::*;

/// Completed affine inverse-propagator coordinates at the physical endpoint.
/// Physical slots retain their order, including dependent and zero-power slots.
/// Extra coordinates are numerator slots, never graph edges or required cuts.
#[derive(Clone, Debug)]
pub struct InversePropagatorBasis {
    coordinates: Vec<Atom>,
    slots: Vec<Atom>,
    labels: Vec<Atom>,
    inverse: BTreeMap<Atom, Atom>,
    physical_slots: usize,
}

fn symbol_of(a: &Atom) -> Result<Symbol> {
    if let AtomView::Var(v) = a.as_view() {
        Ok(v.get_symbol())
    } else {
        Err(Error::InvalidInput(
            "scalar coordinate is not a symbol".into(),
        ))
    }
}

impl InversePropagatorBasis {
    /// Exact affine completion; no statistical rank test or physical reduction.
    pub fn new(coordinates: Vec<Atom>, physical: Vec<Atom>, labels: Vec<Atom>) -> Result<Self> {
        let n = coordinates.len();
        if n == 0
            || coordinates
                .iter()
                .collect::<std::collections::BTreeSet<_>>()
                .len()
                != n
        {
            return Err(Error::InvalidInput(
                "empty or repeated scalar coordinates".into(),
            ));
        }
        for coordinate in &coordinates {
            symbol_of(coordinate)?;
        }
        let zero = coordinates
            .iter()
            .cloned()
            .map(|x| (x, Atom::new()))
            .collect();
        let mut slots = physical.clone();
        let mut selected = Vec::new();
        let mut rows = Vec::new();
        let mut constants = Vec::new();
        for (index, slot) in physical.iter().chain(coordinates.iter()).enumerate() {
            let mut row = Vec::new();
            let constant = substitute(slot, &zero).together().cancel();
            let mut residual = slot - &constant;
            for coordinate in &coordinates {
                let coefficient = slot.derivative(symbol_of(coordinate)?).together().cancel();
                for other in &coordinates {
                    if !coefficient.derivative(symbol_of(other)?).is_zero() {
                        return Err(Error::InvalidInput("nonaffine inverse propagator".into()));
                    }
                }
                residual -= &coefficient * coordinate;
                row.push(coefficient);
            }
            if !residual.together().cancel().is_zero() {
                return Err(Error::InvalidInput("nonaffine inverse propagator".into()));
            }
            // Rational coefficients ensure no hidden kinematic rank conditions.
            if row.iter().any(|a| Rational::try_from(a.as_view()).is_err()) {
                return Err(Error::Unsupported(
                    "basis routing coefficients must be rational".into(),
                ));
            }
            let mut candidate = rows.clone();
            candidate.push(row.clone());
            if rref(candidate).1.len() > rows.len() {
                let slot_index = if index < physical.len() {
                    index
                } else {
                    slots.push(slot.clone());
                    slots.len() - 1
                };
                selected.push(slot_index);
                rows.push(row);
                constants.push(constant);
            }
        }
        if rows.len() != n || labels.len() < slots.len() {
            return Err(Error::InvalidInput(
                "insufficient inverse-propagator labels".into(),
            ));
        }
        let labels = labels[..slots.len()].to_vec();
        if labels
            .iter()
            .collect::<std::collections::BTreeSet<_>>()
            .len()
            != labels.len()
        {
            return Err(Error::InvalidInput(
                "repeated inverse-propagator labels".into(),
            ));
        }
        for label in &labels {
            let symbol = symbol_of(label)?;
            if coordinates.contains(label) || slots.iter().any(|a| !a.derivative(symbol).is_zero())
            {
                return Err(Error::InvalidInput(
                    "inverse-propagator label collides with input".into(),
                ));
            }
        }
        let matrix = inverse(&rows)?;
        let mut inverse = BTreeMap::new();
        for (coordinate, row) in coordinates.iter().zip(matrix) {
            let value = row
                .iter()
                .enumerate()
                .fold(Atom::new(), |sum, (j, c)| {
                    sum + c * (&labels[selected[j]] - &constants[j])
                })
                .expand();
            inverse.insert(coordinate.clone(), value);
        }
        Ok(Self {
            coordinates,
            slots,
            labels,
            inverse,
            physical_slots: physical.len(),
        })
    }

    pub fn coordinates(&self) -> &[Atom] {
        &self.coordinates
    }
    pub fn slots(&self) -> &[Atom] {
        &self.slots
    }
    pub fn labels(&self) -> &[Atom] {
        &self.labels
    }
    pub fn inverse_map(&self) -> &BTreeMap<Atom, Atom> {
        &self.inverse
    }
    pub fn physical_slots(&self) -> usize {
        self.physical_slots
    }

    /// Convert once, before the auxiliary deformation. These indices and
    /// coefficients are then held fixed along that deformation.
    pub fn convert(
        &self,
        numerator: &Atom,
        denominator_powers: &[i16],
    ) -> Result<LinearCombination> {
        if denominator_powers.len() != self.physical_slots {
            return Err(Error::InvalidInput(
                "physical denominator power count".into(),
            ));
        }
        if self
            .labels
            .iter()
            .any(|label| !numerator.derivative(symbol_of(label).unwrap()).is_zero())
        {
            return Err(Error::InvalidInput(
                "input numerator contains a reserved inverse-propagator label".into(),
            ));
        }
        // Native polynomial conversion permits opaque coefficient atoms. Prove
        // every coefficient independent of all requested scalar coordinates.
        self.polynomial(numerator, &self.coordinates)?;
        let converted = substitute(numerator, &self.inverse).expand();
        let monomials = self.polynomial(&converted, &self.labels)?;
        let mut result = LinearCombination::new();
        for (degrees, coefficient) in monomials {
            let indices = degrees
                .iter()
                .enumerate()
                .map(|(j, &degree)| {
                    denominator_powers
                        .get(j)
                        .copied()
                        .unwrap_or(0)
                        .checked_sub(degree)
                        .ok_or_else(|| {
                            Error::Limit("finite-density numerator index overflow".into())
                        })
                })
                .collect::<Result<Vec<_>>>()?;
            *result.entry(Integral(indices)).or_default() += coefficient;
        }
        result.retain(|_, c| !c.is_zero());
        Ok(result)
    }

    fn polynomial(&self, expression: &Atom, variables: &[Atom]) -> Result<Vec<(Vec<i16>, Atom)>> {
        exact_coefficient_list(expression, variables)?
            .into_iter()
            .map(|(m, c)| {
                let degrees = powers(&m, variables)?;
                if degrees.iter().any(|&n| n < 0)
                    || variables
                        .iter()
                        .any(|v| symbol_of(v).is_ok_and(|s| !c.derivative(s).is_zero()))
                {
                    return Err(Error::InvalidInput(
                        "numerator must be polynomial in scalar and medium coordinates".into(),
                    ));
                }
                Ok((degrees, c))
            })
            .collect()
    }

    /// Differentiate the complete represented integrand with respect to an
    /// independent physical squared mass. Coefficient derivatives are essential
    /// to keep the original momentum numerator fixed.
    pub fn mass_derivative(
        &self,
        target: &LinearCombination,
        mass: Symbol,
    ) -> Result<LinearCombination> {
        let mut result = LinearCombination::new();
        for (integral, coefficient) in target {
            if integral.0.len() != self.slots.len() {
                return Err(Error::InvalidInput("basis index count".into()));
            }
            *result.entry(integral.clone()).or_default() += coefficient.derivative(mass);
            for (j, slot) in self.slots.iter().enumerate() {
                let derivative = slot.derivative(mass).together().cancel();
                if derivative.is_zero() || integral.0[j] == 0 {
                    continue;
                }
                if self
                    .coordinates
                    .iter()
                    .any(|x| !derivative.derivative(symbol_of(x).unwrap()).is_zero())
                {
                    return Err(Error::Unsupported(
                        "mass-dependent quadratic routing".into(),
                    ));
                }
                let mut raised = integral.clone();
                raised.0[j] = raised.0[j]
                    .checked_add(1)
                    .ok_or_else(|| Error::Limit("raised physical index overflow".into()))?;
                *result.entry(raised).or_default() -=
                    coefficient * Atom::num(integral.0[j]) * derivative;
            }
        }
        for coefficient in result.values_mut() {
            *coefficient = coefficient.together().cancel();
        }
        result.retain(|_, c| !c.is_zero());
        Ok(result)
    }

    pub fn reconstruct(&self, target: &LinearCombination) -> Result<Atom> {
        let mut expression = Atom::new();
        for (integral, coefficient) in target {
            if integral.0.len() != self.slots.len() {
                return Err(Error::InvalidInput("basis index count".into()));
            }
            expression += integral
                .0
                .iter()
                .zip(&self.slots)
                .fold(coefficient.clone(), |a, (&n, d)| {
                    a * d.clone().pow(-i64::from(n))
                });
        }
        Ok(expression.together().cancel())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn physical_mass_derivative_keeps_original_medium_numerator_fixed() {
        let g = parse!("fd_basis_g");
        let u = parse!("fd_basis_u");
        let m = symbol!("fd_basis_mass");
        let basis = InversePropagatorBasis::new(
            vec![g.clone(), u.clone()],
            vec![&g + Atom::var(m)],
            vec![parse!("fd_basis_r0"), parse!("fd_basis_r1")],
        )
        .unwrap();
        let numerator = &g * &u + u.clone().pow(2);
        let target = basis.convert(&numerator, &[1]).unwrap();
        let derivative = basis.mass_derivative(&target, m).unwrap();
        let expected = -numerator / (&g + Atom::var(m)).pow(2);
        assert!(
            (basis.reconstruct(&derivative).unwrap() - expected)
                .together()
                .cancel()
                .is_zero()
        );
        assert!(basis.convert(&(Atom::num(1) / &u), &[1]).is_err());
        assert!(basis.convert(&parse!("sin(fd_basis_u)"), &[1]).is_err());
        assert!(basis.convert(&parse!("fd_basis_r0"), &[1]).is_err());
    }
    #[test]
    fn redundant_physical_slots_are_preserved() {
        let g = parse!("fd_redundant_g");
        let basis = InversePropagatorBasis::new(
            vec![g.clone()],
            vec![&g + Atom::num(1), &g + Atom::num(2)],
            vec![parse!("fd_redundant_r0"), parse!("fd_redundant_r1")],
        )
        .unwrap();
        let target = basis.convert(&g, &[1, 0]).unwrap();
        assert!(target.keys().all(|i| i.0.len() == 2 && i.0[1] == 0));
        assert!(
            (basis.reconstruct(&target).unwrap() - &g / (&g + Atom::num(1)))
                .together()
                .cancel()
                .is_zero()
        );
    }
}
