//! Shared construction of occupied physical sources and their fixed-shell flow.
use std::collections::BTreeMap;

use symbolica::prelude::*;

use super::geometry::OccupiedCutFamily;
use super::guarded::{GuardedContext, GuardedMeasureIdentity, IndexBounds, IndexDomain};
use super::reduction::{AuxiliaryConvention, FixedShellDeformation};
use crate::{Error, Result};

/// Equivalent exact source presentations for bounded native discovery.
/// Selection controls source rows and their integer-domain partition, never
/// integral ordering, physical zeros, or the elimination backend.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum WeightedSourcePolicy {
    #[default]
    LegacyLorentz,
    ActiveLorentz,
    LorentzThenTangents,
    TangentsThenLorentz,
    NormalsThenLorentz,
    LorentzThenNormals,
}

impl WeightedSourcePolicy {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::LegacyLorentz => "legacy-lorentz",
            Self::ActiveLorentz => "active-lorentz",
            Self::LorentzThenTangents => "lorentz-then-tangents",
            Self::TangentsThenLorentz => "tangents-then-lorentz",
            Self::NormalsThenLorentz => "normals-then-lorentz",
            Self::LorentzThenNormals => "lorentz-then-normals",
        }
    }
}

impl std::str::FromStr for WeightedSourcePolicy {
    type Err = Error;

    fn from_str(value: &str) -> Result<Self> {
        match value {
            "legacy-lorentz" => Ok(Self::LegacyLorentz),
            "active-lorentz" => Ok(Self::ActiveLorentz),
            "lorentz-then-tangents" => Ok(Self::LorentzThenTangents),
            "tangents-then-lorentz" => Ok(Self::TangentsThenLorentz),
            "normals-then-lorentz" => Ok(Self::NormalsThenLorentz),
            "lorentz-then-normals" => Ok(Self::LorentzThenNormals),
            _ => Err(Error::InvalidInput(format!(
                "unknown weighted source presentation {value}"
            ))),
        }
    }
}

/// Physical source admission and presentation. Positive completion powers are
/// opt-in and require an exact massive compact-energy certificate per slot.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct WeightedSourceOptions {
    pub policy: WeightedSourcePolicy,
    pub positive_compact_energy_powers: bool,
}

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
        parameters: Vec<Symbol>,
        identity: GuardedMeasureIdentity,
    ) -> Result<PreparedWeightedSources<N>> {
        self.guarded_sources_with_policy(
            epsilon,
            dimension,
            eta,
            shifted,
            domain_budget,
            parameters,
            identity,
            WeightedSourcePolicy::default(),
        )
    }

    /// Select an exact source presentation explicitly. The default factory
    /// preserves the established Lorentz presentation; tangent alternatives are
    /// available for bounded native coverage diagnostics and further closure.
    pub fn guarded_sources_with_policy<const N: usize>(
        &self,
        epsilon: Symbol,
        dimension: i64,
        eta: Symbol,
        shifted: &[usize],
        domain_budget: usize,
        parameters: Vec<Symbol>,
        identity: GuardedMeasureIdentity,
        policy: WeightedSourcePolicy,
    ) -> Result<PreparedWeightedSources<N>> {
        self.guarded_sources_with_options(
            epsilon,
            dimension,
            eta,
            shifted,
            domain_budget,
            parameters,
            identity,
            WeightedSourceOptions {
                policy,
                ..Default::default()
            },
        )
    }

    /// Permit inverse powers of a completion only when it is certified equal
    /// to a nonzero rational multiple of one future energy with m²>0. Such a
    /// factor is nonsingular on every shell/surface distribution support. All
    /// other completions retain the polynomial domain; no numerical boundary
    /// support is inferred merely from this source admission.
    pub fn guarded_sources_with_options<const N: usize>(
        &self,
        epsilon: Symbol,
        dimension: i64,
        eta: Symbol,
        shifted: &[usize],
        domain_budget: usize,
        mut parameters: Vec<Symbol>,
        mut identity: GuardedMeasureIdentity,
        options: WeightedSourceOptions,
    ) -> Result<PreparedWeightedSources<N>> {
        let policy = options.policy;
        let measure = self.deformed_measure::<N>(eta, shifted)?;
        let indices =
            std::array::from_fn(|slot| symbol!(format!("rustflow_occupied_indices::a_{slot}")));
        let symbolic_dimension = Atom::num(dimension) - Atom::num(2) * Atom::var(epsilon);
        let lorentz =
            || measure.lorentz_ibps(self.loops(), &symbolic_dimension, &indices, domain_budget);
        let tangent = || {
            measure.compact_tangent_ibps(
                self.loops(),
                self.shells().len(),
                &symbolic_dimension,
                &indices,
                domain_budget,
            )
        };
        let legacy = || {
            measure.lorentz_ibps_legacy(self.loops(), &symbolic_dimension, &indices, domain_budget)
        };
        let normal = || {
            if !options.positive_compact_energy_powers {
                return Err(Error::InvalidInput(
                    "shell-normal sources require certified positive compact energy powers".into(),
                ));
            }
            let completions = (self.physical_slots()..self.input_slots())
                .filter_map(|slot| self.compact_energy_completion(slot))
                .map(|certificate| {
                    (
                        certificate.loop_index,
                        certificate.slot,
                        certificate.coefficient,
                    )
                })
                .collect::<Vec<_>>();
            if completions.is_empty() {
                return Err(Error::Unsupported(
                    "shell-normal sources have no certified massive energy completion".into(),
                ));
            }
            measure.compact_normal_ibps(self.loops(), &completions, &indices, domain_budget)
        };
        let mut identities = match policy {
            WeightedSourcePolicy::LegacyLorentz => legacy()?,
            WeightedSourcePolicy::ActiveLorentz => lorentz()?,
            WeightedSourcePolicy::LorentzThenTangents => {
                let mut sources = lorentz()?;
                sources.extend(tangent()?);
                sources
            }
            WeightedSourcePolicy::TangentsThenLorentz => {
                let mut sources = tangent()?;
                sources.extend(lorentz()?);
                sources
            }
            WeightedSourcePolicy::NormalsThenLorentz => {
                let mut sources = normal()?;
                sources.extend(legacy()?);
                sources
            }
            WeightedSourcePolicy::LorentzThenNormals => {
                let mut sources = legacy()?;
                sources.extend(normal()?);
                sources
            }
        };
        identity
            .measure
            .push_str(&format!("; source presentation={}", policy.as_str()));
        let roles = *measure.roles();
        let mut bounds = *IndexDomain::for_roles(&roles).bounds();
        identity.measure.push_str(&format!(
            "; positive compact energy completion powers={}",
            options.positive_compact_energy_powers,
        ));
        for slot in self.physical_slots()..self.input_slots() {
            if options.positive_compact_energy_powers {
                if let Some(certificate) = self.compact_energy_completion(slot) {
                    identity.support.push_str(&format!(
                        "; certified nonsingular completion slot {}={}*E_{} on future required shell with mass squared {}>0; all integer completion powers admitted",
                        certificate.slot, certificate.coefficient, certificate.loop_index, certificate.mass_squared,
                    ));
                    continue;
                }
            }
            bounds[slot] = IndexBounds::new(None, Some(0))
                .map_err(|error| Error::InvalidInput(error.to_string()))?;
        }
        let admitted_domain =
            IndexDomain::new(bounds).map_err(|error| Error::InvalidInput(error.to_string()))?;
        identities = identities
            .into_iter()
            .filter_map(|mut source| {
                source.domain = source.domain.intersection(&admitted_domain)?;
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
                .with_admitted_domain(admitted_domain)?;
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
