//! Exact continuation of a certified input into oriented occupied cut factors.
//! This is a measure construction, not a contour-admission or evaluation certificate.
use super::PreparedDensityInput;
use super::guarded::IndexRole;
use super::measure::WeightedMeasure;
use crate::algebra::{determinant, inverse};
use crate::family::{LinearCombination, substitute};
use crate::{Error, Integral, IntegralFamily, Propagator, Result};
use std::collections::BTreeMap;
use symbolica::prelude::*;

/// One occupied energy with its independently variable physical shell mass.
#[derive(Clone, Debug)]
pub struct OccupiedShell {
    pub loop_index: usize,
    pub physical_slot: usize,
    pub upper_slot: usize,
    pub lower_slot: usize,
    pub mass_squared: Atom,
    pub chemical_potential: Rational,
}

/// Native factors in Gram coordinates, followed by contractions with u (u²=1).
/// The first `shells.len()` loops are future-oriented compact momenta. Remaining
/// loops use the ordinary native +i0 contour. Physical masses remain independent
/// until `at_physical_masses` is called. Target coefficients include the complete
/// Euclidean quadratic/medium phase, but no angular or routing measure factor.
#[derive(Clone, Debug)]
pub struct OccupiedCutFamily {
    coordinates: Vec<Atom>,
    factors: Vec<Atom>,
    roles: Vec<IndexRole>,
    targets: Vec<LinearCombination>,
    shells: Vec<OccupiedShell>,
    physical_slots: usize,
    input_slots: usize,
    loops: usize,
    name: String,
    /// old globally reversed Minkowski loops = inverse_routing * new loops.
    inverse_routing: Vec<Vec<Atom>>,
    routing_determinant: Rational,
    mass_assignments: BTreeMap<Atom, Atom>,
}

impl PreparedDensityInput {
    /// Construct a nonempty, connected cut complement admitted by the input
    /// graph. Actual continuation domains must still be established before
    /// reducing or evaluating its formal distribution products.
    pub fn occupied_cut(&self, slots: &[usize], budget: usize) -> Result<OccupiedCutFamily> {
        if slots.is_empty() {
            return Err(Error::InvalidInput(
                "an occupied cut must contain at least one charged edge".into(),
            ));
        }
        self.continued_cut(slots, budget)
    }

    /// The empty set is the ordinary vacuum contribution, with the same exact
    /// Wick conversion of the original polynomial and no occupied factors.
    pub(crate) fn continued_cut(
        &self,
        slots: &[usize],
        budget: usize,
    ) -> Result<OccupiedCutFamily> {
        if slots.windows(2).any(|s| s[0] >= s[1]) {
            return Err(Error::InvalidInput(
                "cut slots must be sorted and distinct".into(),
            ));
        }
        let mut routing = if slots.is_empty() {
            (0..self.input().loops)
                .map(|i| {
                    (0..self.input().loops)
                        .map(|j| Atom::num(i32::from(i == j)))
                        .collect()
                })
                .collect()
        } else {
            let certificate = self
                .cut_decomposition(budget)?
                .into_iter()
                .find(|c| {
                    c["cut_slots"].as_array().is_some_and(|s| {
                        s.len() == slots.len()
                            && s.iter()
                                .zip(slots)
                                .all(|(a, b)| a.as_u64() == Some(*b as u64))
                    })
                })
                .ok_or_else(|| {
                    Error::InvalidInput(
                        "occupied cut is not a connected admitted graph complement".into(),
                    )
                })?;
            let parse = |s: &str| {
                Atom::parse(s, "rustflow_density", Default::default())
                    .map_err(|e| Error::InvalidInput(e.to_string()))
            };
            certificate["loop_routing"]
                .as_array()
                .unwrap()
                .iter()
                .map(|row| {
                    row.as_array()
                        .unwrap()
                        .iter()
                        .map(|a| parse(a.as_str().unwrap()))
                        .collect::<Result<Vec<_>>>()
                })
                .collect::<Result<Vec<_>>>()?
        };
        let parse = |s: &str| {
            Atom::parse(s, "rustflow_density", Default::default())
                .map_err(|e| Error::InvalidInput(e.to_string()))
        };
        let input = self.input();
        let l = input.loops;
        let scalar_count = l * (l + 1) / 2;
        let potentials = input.chemical_potentials.iter().map(|s| {
            Rational::try_from(parse(s)?.as_view()).map_err(|_| Error::Unsupported(
                "occupied orientation requires exact real chemical potentials at the evaluation point".into()))
        }).collect::<Result<Vec<_>>>()?;
        let mut magnitudes = Vec::new();
        for (i, &slot) in slots.iter().enumerate() {
            let mu = input.edges[slot]
                .charges
                .iter()
                .zip(&potentials)
                .fold(Rational::zero(), |s, (&q, m)| s + m * &Rational::from(q));
            if mu.is_zero() {
                return Err(Error::InvalidInput(
                    "a zero chemical potential has no oriented occupied cut".into(),
                ));
            }
            let sign = if mu < Rational::zero() { -1 } else { 1 };
            for c in &mut routing[i] {
                *c *= Atom::num(sign);
            }
            magnitudes.push(mu * Rational::from(sign));
        }
        let det = Rational::try_from(determinant(routing.clone()).as_view())
            .map_err(|_| Error::InvalidInput("nonrational occupied routing determinant".into()))?;
        let inverse_routing = inverse(&routing)?;
        let coordinates = (0..scalar_count + l)
            .map(|i| Atom::var(symbol!(format!("rustflow_occupied::q_{i}"))))
            .collect::<Vec<_>>();
        for c in &coordinates {
            let AtomView::Var(v) = c.as_view() else {
                unreachable!()
            };
            if self
                .physical_masses()
                .iter()
                .chain(self.targets().iter().flat_map(|t| t.values()))
                .any(|a| !a.derivative(v.get_symbol()).is_zero())
            {
                return Err(Error::InvalidInput(
                    "input uses a reserved occupied coordinate".into(),
                ));
            }
        }
        let dot = |i: usize, j: usize| {
            let mut sum = Atom::zero();
            let mut n = 0;
            for a in 0..l {
                for b in a..l {
                    let c = if a == b {
                        &inverse_routing[i][a] * &inverse_routing[j][b]
                    } else {
                        &inverse_routing[i][a] * &inverse_routing[j][b]
                            + &inverse_routing[i][b] * &inverse_routing[j][a]
                    };
                    sum += c * &coordinates[n];
                    n += 1;
                }
            }
            sum.expand()
        };
        let mut wick = BTreeMap::new();
        let mut n = 0;
        for i in 0..l {
            for j in i..l {
                wick.insert(self.basis().coordinates()[n].clone(), -dot(i, j));
                n += 1;
            }
        }
        for i in 0..l {
            let energy = (0..l).fold(Atom::zero(), |s, j| {
                s + &inverse_routing[i][j] * &coordinates[scalar_count + j]
            });
            wick.insert(
                self.basis().coordinates()[scalar_count + i].clone(),
                Atom::i() * energy,
            );
        }
        let input_slots = self.basis().slots().len();
        let physical_slots = input.edges.len();
        let mut scales = Vec::new();
        let mut factors = Vec::new();
        for (i, slot) in self.basis().slots().iter().enumerate() {
            let medium =
                i >= physical_slots && self.basis().coordinates()[scalar_count..].contains(slot);
            let scale = if medium { Atom::i() } else { Atom::num(-1) };
            factors.push(
                (substitute(slot, &wick) / &scale)
                    .expand()
                    .together()
                    .cancel(),
            );
            scales.push(scale);
        }
        let mut roles = (0..input_slots)
            .map(|i| {
                if slots.contains(&i) {
                    IndexRole::RequiredCut
                } else {
                    IndexRole::Ordinary
                }
            })
            .collect::<Vec<_>>();
        let mut shells = Vec::new();
        for (i, (&slot, mu)) in slots.iter().zip(magnitudes).enumerate() {
            let upper_slot = factors.len();
            factors.push(Atom::num(mu.clone()) - &coordinates[scalar_count + i]);
            roles.push(IndexRole::Occupation);
            let lower_slot = factors.len();
            factors.push(coordinates[scalar_count + i].clone());
            roles.push(IndexRole::Occupation);
            shells.push(OccupiedShell {
                loop_index: i,
                physical_slot: slot,
                upper_slot,
                lower_slot,
                mass_squared: Atom::var(self.independent_masses()[slot]),
                chemical_potential: mu,
            });
        }
        let targets = self
            .targets()
            .iter()
            .map(|terms| {
                let mut out = LinearCombination::new();
                for (integral, coefficient) in terms {
                    if slots.iter().any(|&i| integral.0[i] <= 0) {
                        continue;
                    }
                    let mut coefficient = coefficient.clone();
                    for (&power, scale) in integral.0.iter().zip(&scales) {
                        coefficient *= scale.pow(-i64::from(power));
                    }
                    let mut index = integral.0.clone();
                    index.resize(factors.len(), 0);
                    *out.entry(Integral(index)).or_default() += coefficient;
                }
                out.retain(|_, c| !c.is_zero());
                out
            })
            .collect();
        Ok(OccupiedCutFamily {
            coordinates,
            factors,
            roles,
            targets,
            shells,
            physical_slots,
            input_slots,
            loops: l,
            name: format!("{} occupied {slots:?}", input.name),
            inverse_routing,
            routing_determinant: det,
            mass_assignments: self
                .independent_masses()
                .iter()
                .zip(self.physical_masses())
                .map(|(&s, m)| (Atom::var(s), m.clone()))
                .collect(),
        })
    }
}

impl OccupiedCutFamily {
    pub fn coordinates(&self) -> &[Atom] {
        &self.coordinates
    }
    pub fn factors(&self) -> &[Atom] {
        &self.factors
    }
    pub fn roles(&self) -> &[IndexRole] {
        &self.roles
    }
    pub fn targets(&self) -> &[LinearCombination] {
        &self.targets
    }
    pub fn shells(&self) -> &[OccupiedShell] {
        &self.shells
    }
    pub fn loops(&self) -> usize {
        self.loops
    }
    pub fn physical_slots(&self) -> usize {
        self.physical_slots
    }
    pub fn input_slots(&self) -> usize {
        self.input_slots
    }
    pub fn inverse_routing(&self) -> &[Vec<Atom>] {
        &self.inverse_routing
    }
    pub fn routing_determinant(&self) -> &Rational {
        &self.routing_determinant
    }

    /// Apply equal-mass/physical assignments only after all independent physical
    /// mass derivatives of the original target and distributions have been taken.
    pub fn at_physical_masses(&self) -> Self {
        let mut out = self.clone();
        for factor in &mut out.factors {
            *factor = substitute(factor, &self.mass_assignments);
        }
        for target in &mut out.targets {
            for coefficient in target.values_mut() {
                *coefficient = substitute(coefficient, &self.mass_assignments)
                    .together()
                    .cancel();
            }
        }
        for shell in &mut out.shells {
            shell.mass_squared = substitute(&shell.mass_squared, &self.mass_assignments);
        }
        out
    }

    /// Native D−eta deformation; shell and occupation factors, and all numerator
    /// completion slots, are fixed. Admissibility of its contour is independent
    /// of this exact algebraic construction.
    pub fn deformed_measure<const N: usize>(
        &self,
        eta: Symbol,
        shifted: &[usize],
    ) -> Result<WeightedMeasure<N>> {
        if shifted.is_empty()
            || shifted
                .iter()
                .any(|&i| i >= self.physical_slots || self.roles[i] != IndexRole::Ordinary)
            || shifted
                .iter()
                .collect::<std::collections::BTreeSet<_>>()
                .len()
                != shifted.len()
        {
            return Err(Error::InvalidInput(
                "occupied flow must shift distinct uncut physical factors".into(),
            ));
        }
        if self
            .factors
            .iter()
            .chain(self.targets.iter().flat_map(|t| t.values()))
            .any(|a| !a.derivative(eta).is_zero())
            || self.coordinates.iter().any(|a| *a == Atom::var(eta))
        {
            return Err(Error::InvalidInput(
                "auxiliary mass collides with the occupied input".into(),
            ));
        }
        let mut factors: [Atom; N] = self
            .factors
            .clone()
            .try_into()
            .map_err(|_| Error::InvalidInput("weighted slot arity mismatch".into()))?;
        let roles: [IndexRole; N] = self
            .roles
            .clone()
            .try_into()
            .map_err(|_| Error::InvalidInput("weighted slot arity mismatch".into()))?;
        for &i in shifted {
            factors[i] -= Atom::var(eta);
        }
        WeightedMeasure::new(self.coordinates.clone(), factors, roles)
    }

    /// Geometry for branch enumeration/ordinary hard projection. Occupation
    /// factors are deliberately absent: they are distributions, not denominators.
    pub fn region_family(&self, epsilon: Symbol, dimension: i64) -> Result<IntegralFamily> {
        let zero = self
            .coordinates
            .iter()
            .map(|c| (c.clone(), Atom::zero()))
            .collect();
        let propagators = self.factors[..self.input_slots]
            .iter()
            .map(|a| {
                Ok(Propagator {
                    constant: substitute(a, &zero),
                    scalar_products: self
                        .coordinates
                        .iter()
                        .map(|c| {
                            let AtomView::Var(v) = c.as_view() else {
                                unreachable!()
                            };
                            a.derivative(v.get_symbol())
                        })
                        .collect(),
                })
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(IntegralFamily {
            name: self.name.clone(),
            loops: (0..self.loops).map(|i| format!("q{i}")).collect(),
            external: vec!["u".into()],
            external_gram: vec![vec![Atom::one()]],
            propagators,
            physical_propagators: self.physical_slots,
            epsilon,
            dimension,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_density::DensityInput;

    fn input() -> DensityInput {
        serde_json::from_str(include_str!(
            "../../examples/finite_density/massive_two_loop_sunset.json"
        ))
        .unwrap()
    }

    #[test]
    fn occupied_routing_phases_support_and_fixed_shell_deformation() {
        let prepared = input().prepare().unwrap();
        let cut = prepared
            .occupied_cut(&[0, 1], 16)
            .unwrap()
            .at_physical_masses();
        let q = cut.coordinates();
        assert_eq!(cut.routing_determinant(), &Rational::one());
        assert!(
            (&cut.factors()[0] - &q[0] + Atom::num((1, 4)))
                .expand()
                .is_zero()
        );
        assert!(
            (&cut.factors()[1] - &q[2] + Atom::num((1, 4)))
                .expand()
                .is_zero()
        );
        assert!(
            (&cut.factors()[2] - (&q[0] - Atom::num(2) * &q[1] + &q[2] - Atom::one()))
                .expand()
                .is_zero()
        );
        assert_eq!(cut.factors()[5], Atom::one() - &q[3]);
        assert_eq!(cut.factors()[6], q[3]);
        assert_eq!(
            cut.targets()[0],
            BTreeMap::from([(Integral(vec![1, 1, 1, 0, 0, 0, 0, 0, 0]), Atom::num(-1))])
        );
        let eta = symbol!("occupied_geometry_test::eta");
        let measure = cut.deformed_measure::<9>(eta, &[2]).unwrap();
        assert_eq!(measure.roles()[0], IndexRole::RequiredCut);
        assert_eq!(measure.roles()[5], IndexRole::Occupation);
        assert!(cut.deformed_measure::<9>(eta, &[0]).is_err());
        assert!(cut.deformed_measure::<9>(eta, &[3]).is_err());
        assert!(cut.deformed_measure::<7>(eta, &[2]).is_err());
        let family = cut
            .region_family(symbol!("occupied_geometry_test::eps"), 4)
            .unwrap();
        assert_eq!(family.propagators.len(), 5);
        assert_eq!(family.physical_propagators, 3);
        family.validate().unwrap();
        let one = prepared
            .occupied_cut(&[1], 16)
            .unwrap()
            .at_physical_masses();
        assert_eq!(one.shells()[0].physical_slot, 1);
        assert!(
            (&one.factors()[1] - &one.coordinates()[0] + Atom::num((1, 4)))
                .expand()
                .is_zero()
        );
        one.deformed_measure::<7>(eta, &[0, 2]).unwrap();
    }

    #[test]
    fn negative_potentials_reverse_future_coordinates_and_missing_cuts_vanish() {
        let mut original = input();
        original.chemical_potentials = vec!["-1".into()];
        original.targets = vec![super::super::DensityTarget {
            powers: vec![1, 0, 1],
            numerator: "u1".into(),
        }];
        let prepared = original.prepare().unwrap();
        let cut = prepared
            .occupied_cut(&[0, 1], 16)
            .unwrap()
            .at_physical_masses();
        assert!(cut.targets()[0].is_empty());
        assert_eq!(cut.shells()[0].chemical_potential, Rational::one());
        assert_eq!(cut.factors()[3], -&cut.coordinates()[3]);
        assert!(prepared.occupied_cut(&[2], 16).is_err());
        assert!(prepared.occupied_cut(&[1, 0], 16).is_err());
        let single = prepared.occupied_cut(&[0], 16).unwrap();
        assert_eq!(single.targets()[0].values().next().unwrap(), &Atom::i());
    }
}
