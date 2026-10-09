//! Shared construction of occupied physical sources and their fixed-shell flow.
use std::collections::BTreeMap;

use symbolica::prelude::*;

use super::geometry::OccupiedCutFamily;
use super::guarded::{GuardedContext, GuardedMeasureIdentity, IndexBounds, IndexDomain};
use super::reduction::{AuxiliaryConvention, FixedShellDeformation};
use crate::{Error, Result};

/// Exact algebraic preparation. The measure identity must describe the actual
/// contour/support prescription; construction alone does not certify that the
/// products of distributions or continuation across thresholds are admissible.
pub struct PreparedWeightedSources<const N: usize> {
    pub context: GuardedContext<N>,
    pub deformation: FixedShellDeformation<N>,
    pub targets: Vec<BTreeMap<[i64; N], Atom>>,
}

impl OccupiedCutFamily {
    /// Derive Lorentz IBPs and distribution multiplication sources, with every
    /// numerator-completion power restricted to nonpositive integers. Physical
    /// masses must be assigned before this method, or included as parameters.
    pub fn guarded_sources<const N: usize>(
        &self,
        epsilon: Symbol,
        dimension: i64,
        eta: Symbol,
        shifted: &[usize],
        domain_budget: usize,
        mut parameters: Vec<Symbol>,
        mut identity: GuardedMeasureIdentity,
    ) -> Result<PreparedWeightedSources<N>> {
        let measure = self.deformed_measure::<N>(eta, shifted)?;
        let indices =
            std::array::from_fn(|slot| symbol!(format!("rustflow_occupied_indices::a_{slot}")));
        let mut identities = measure.lorentz_ibps(
            self.loops(),
            &(Atom::num(dimension) - Atom::num(2) * Atom::var(epsilon)),
            &indices,
            domain_budget,
        )?;
        let roles = *measure.roles();
        let mut bounds = *IndexDomain::for_roles(&roles).bounds();
        for bound in &mut bounds[self.physical_slots()..self.input_slots()] {
            *bound = IndexBounds::new(None, Some(0))
                .map_err(|error| Error::InvalidInput(error.to_string()))?;
        }
        let polynomial_domain =
            IndexDomain::new(bounds).map_err(|error| Error::InvalidInput(error.to_string()))?;
        identities = identities
            .into_iter()
            .filter_map(|mut source| {
                source.domain = source.domain.intersection(&polynomial_domain)?;
                Some(source)
            })
            .collect();
        // On a real occupied shell E^2-|p|^2=m^2>0, the lower endpoint E=0
        // has empty support, including every derivative of either delta. This
        // certificate uses assigned positive rational masses; no massless or
        // symbolic-mass limiting value is inferred. It is separate from IBP
        // algebra because the complex algebraic intersection is not empty.
        let mut zero_domains = Vec::new();
        for shell in self.shells() {
            let Ok(mass_squared) = Rational::try_from(shell.mass_squared.as_view()) else {
                continue;
            };
            if mass_squared <= Rational::zero() {
                continue;
            }
            let mut support = bounds;
            support[shell.physical_slot] = IndexBounds::new(Some(1), None)
                .map_err(|error| Error::InvalidInput(error.to_string()))?;
            support[shell.lower_slot] = IndexBounds::new(Some(1), None)
                .map_err(|error| Error::InvalidInput(error.to_string()))?;
            zero_domains.push(
                IndexDomain::new(support)
                    .map_err(|error| Error::InvalidInput(error.to_string()))?,
            );
            identity.support.push_str(&format!(
                "; certified real empty support: cut slot {}>=1 and lower-energy slot {}>=1, physical mass squared {}>0",
                shell.physical_slot, shell.lower_slot, mass_squared
            ));
        }
        for parameter in [eta, epsilon] {
            if !parameters.contains(&parameter) {
                parameters.push(parameter);
            }
        }
        let context = GuardedContext::new(identity, roles, indices, parameters, identities)?
            .with_measure_zero_domains(zero_domains)?;
        let mut selected = [false; N];
        for &slot in shifted {
            selected[slot] = true;
        }
        let deformation =
            FixedShellDeformation::new(eta, AuxiliaryConvention::NativeMinusEta, roles, selected)?
                .with_admitted_domain(polynomial_domain)?;
        let targets = self
            .targets()
            .iter()
            .map(|target| {
                target
                    .iter()
                    .map(|(integral, coefficient)| {
                        let indices = integral
                            .0
                            .iter()
                            .map(|&index| i64::from(index))
                            .collect::<Vec<_>>()
                            .try_into()
                            .map_err(|_| Error::InvalidInput("weighted target arity".into()))?;
                        Ok((indices, coefficient.clone()))
                    })
                    .collect::<Result<_>>()
            })
            .collect::<Result<_>>()?;
        Ok(PreparedWeightedSources {
            context,
            deformation,
            targets,
        })
    }
}
