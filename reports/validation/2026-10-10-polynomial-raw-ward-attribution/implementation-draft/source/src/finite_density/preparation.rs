//! Shared construction of occupied physical sources and their fixed-shell flow.
use std::collections::BTreeMap;

use symbolica::prelude::*;

use super::geometry::OccupiedCutFamily;
use super::guarded::{GuardedContext, GuardedMeasureIdentity, IndexBounds, IndexDomain};
use super::massless_endpoint::MasslessFlowEvidence;
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
    /// Attach measure-aware polynomial zeros only through a bound sealed
    /// massless endpoint/origin proof. The default preserves existing sources.
    #[serde(default)]
    pub free_virtual_zero_sectors: bool,
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
        parameters: Vec<Symbol>,
        identity: GuardedMeasureIdentity,
        options: WeightedSourceOptions,
    ) -> Result<PreparedWeightedSources<N>> {
        self.guarded_sources_with_origin(
            epsilon,
            dimension,
            eta,
            shifted,
            domain_budget,
            parameters,
            identity,
            options,
            None,
            false,
        )
    }

    /// Use a sealed finite-positive-eta origin certificate for the proved
    /// massless flow classes. The resulting lower-contact zeros are defined by
    /// joint high-dimensional continuation of the complete smooth kernel; they
    /// are not inferred from massive empty support or a separately chosen PV.
    pub fn guarded_sources_with_massless_origin<const N: usize>(
        &self,
        epsilon: Symbol,
        dimension: i64,
        eta: Symbol,
        shifted: &[usize],
        domain_budget: usize,
        parameters: Vec<Symbol>,
        identity: GuardedMeasureIdentity,
        options: WeightedSourceOptions,
        origin: &MasslessFlowEvidence,
    ) -> Result<PreparedWeightedSources<N>> {
        self.guarded_sources_with_origin(
            epsilon,
            dimension,
            eta,
            shifted,
            domain_budget,
            parameters,
            identity,
            options,
            Some(origin),
            false,
        )
    }

    /// Explicit source-only presentation: prefix the polynomial raw Ward pair
    /// to the requested legacy/tangent policy. Requires a bound singleton-germ
    /// proof; other classes fail rather than silently applying another theorem.
    /// Existing factories, source options, zero domains and flow dispatch stay
    /// unchanged. This does not assert native closure or numerical accuracy.
    pub fn guarded_sources_with_massless_raw_ward<const N: usize>(
        &self,
        epsilon: Symbol,
        dimension: i64,
        eta: Symbol,
        shifted: &[usize],
        domain_budget: usize,
        parameters: Vec<Symbol>,
        identity: GuardedMeasureIdentity,
        options: WeightedSourceOptions,
        origin: &MasslessFlowEvidence,
    ) -> Result<PreparedWeightedSources<N>> {
        self.guarded_sources_with_origin(
            epsilon,
            dimension,
            eta,
            shifted,
            domain_budget,
            parameters,
            identity,
            options,
            Some(origin),
            true,
        )
    }

    fn guarded_sources_with_origin<const N: usize>(
        &self,
        epsilon: Symbol,
        dimension: i64,
        eta: Symbol,
        shifted: &[usize],
        domain_budget: usize,
        mut parameters: Vec<Symbol>,
        mut identity: GuardedMeasureIdentity,
        options: WeightedSourceOptions,
        origin: Option<&MasslessFlowEvidence>,
        raw_singleton_ward: bool,
    ) -> Result<PreparedWeightedSources<N>> {
        if options.free_virtual_zero_sectors && origin.is_none() {
            return Err(Error::InvalidInput(
                "free virtual zero sectors require bound sealed massless evidence".into(),
            ));
        }
        let massless_origin_loops = if let Some(origin) = origin {
            if origin.source_options() != options {
                return Err(Error::InvalidInput(
                    "massless source origin certificate uses different source options".into(),
                ));
            }
            if options.positive_compact_energy_powers {
                return Err(Error::Unsupported(
                    "massless source origin continuation requires polynomial completions".into(),
                ));
            }
            origin.validate_family(self, shifted)?;
            let loops = origin.certified_origin_loops(self, shifted)?;
            if loops.is_empty() {
                return Err(Error::InvalidInput(
                    "massless source origin certificate covers no occupied loop".into(),
                ));
            }
            identity.measure.push_str(&format!(
                "; massless source origin evidence={}",
                origin.source_identity()
            ));
            identity
                .branch
                .push_str(&format!("; {}", origin.origin_identity()));
            loops
        } else {
            Vec::new()
        };
        let free_virtual_zeros = if options.free_virtual_zero_sectors {
            origin
                .unwrap()
                .free_virtual_zero_supports(self, shifted, domain_budget)?
        } else {
            Vec::new()
        };
        let policy = options.policy;
        let measure = self.deformed_measure::<N>(eta, shifted)?;
        let physical_arity = measure.physical_arity();
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
        // An additive, explicitly requested presentation. The typed constructor
        // rejects all non-singleton classes and binds the exact deformed factors.
        // Defaults and existing source options keep their established semantics.
        let ward = if raw_singleton_ward {
            let evidence = origin
                .ok_or_else(|| {
                    Error::InvalidInput(
                        "raw singleton Ward sources require sealed massless origin evidence".into(),
                    )
                })?
                .singleton_raw_ward_evidence::<N>(self, shifted, eta)?;
            identity
                .measure
                .push_str(&format!("; raw Ward evidence={}", evidence.identity()));
            measure.raw_singleton_ward_sources(
                &evidence,
                &symbolic_dimension,
                &indices,
                domain_budget,
            )?
        } else {
            Vec::new()
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
        if !ward.is_empty() {
            // Raw (not upper-recentered) rows are the validated presentation.
            // These rows add no inverse-energy domain or support-zero premise.
            let mut prefixed = ward;
            prefixed.extend(identities);
            identities = prefixed;
        }
        identity
            .measure
            .push_str(&format!("; source presentation={}", policy.as_str()));
        let roles = *measure.roles();
        let mut bounds = *IndexDomain::for_roles(&roles).bounds();
        // Backend storage coordinates are absent from the physical measure.
        // Freeze them in source, support-zero, and deformation domains alike.
        bounds[physical_arity..].fill(IndexBounds::fixed(0));
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
        // This is a different physical statement from massive empty support.
        // At fixed eta>0 the sealed flow proof supplies a smooth complete
        // virtual/transfer kernel at each origin. With polynomial insertions,
        // every finite C_n H_lower,l label has vanishing required zero-energy
        // jets on a sufficiently high-Re(D) open domain. Joint meromorphic
        // continuation then gives zero pointwise in the integer-index family;
        // no one finite dimension is asserted to cover all unbounded labels.
        for loop_index in massless_origin_loops {
            let shell = self
                .shells()
                .iter()
                .find(|shell| shell.loop_index == loop_index)
                .ok_or_else(|| {
                    Error::InvalidInput(
                        "massless origin certificate refers to an absent compact loop".into(),
                    )
                })?;
            if !shell.mass_squared.is_zero() || shell.chemical_potential <= Rational::zero() {
                return Err(Error::InvalidInput("massless origin source zero needs mass zero and a positive separated upper endpoint".into()));
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
                "; certified joint dimensional origin zero jets: occupied loop {loop_index}, cut slot {}>=1 and lower-energy slot {}>=1; physical mass zero; fixed positive auxiliary mass; polynomial completions; upper endpoint {}>0",
                shell.physical_slot, shell.lower_slot, shell.chemical_potential,
            ));
        }
        // Only unrestricted virtual directions occur in these certificates.
        // Intersect the actual admitted box, preserving every completion,
        // occupation and frozen storage coordinate before exposing a zero.
        for proof in free_virtual_zeros {
            let mut support = bounds;
            for shell in self.shells() {
                support[shell.physical_slot] = IndexBounds::new(Some(1), None)
                    .map_err(|e| Error::InvalidInput(e.to_string()))?;
            }
            for &slot in proof.forced_nonpositive_slots() {
                if slot >= self.physical_slots()
                    || roles[slot] != super::guarded::IndexRole::Ordinary
                {
                    return Err(Error::InvalidInput(
                        "free-virtual zero certificate includes a compact or absent slot".into(),
                    ));
                }
                support[slot] = IndexBounds::new(None, Some(0))
                    .map_err(|e| Error::InvalidInput(e.to_string()))?;
            }
            let proposed =
                IndexDomain::new(support).map_err(|e| Error::InvalidInput(e.to_string()))?;
            if let Some(domain) = proposed.intersection(&admitted_domain) {
                zero_domains.push(domain);
            }
            identity.support.push_str(&format!(
                "; certified free virtual zero support: {}",
                proof.identity()
            ));
        }
        for parameter in [eta, epsilon] {
            if !parameters.contains(&parameter) {
                parameters.push(parameter);
            }
        }
        let context = GuardedContext::new_with_physical_arity(
            identity,
            roles,
            indices,
            parameters,
            identities,
            physical_arity,
        )?
        .with_measure_zero_domains(zero_domains)?;
        let mut selected = [false; N];
        for &slot in shifted {
            selected[slot] = true;
        }
        let deformation =
            FixedShellDeformation::new(eta, AuxiliaryConvention::NativeMinusEta, roles, selected)?
                .with_physical_arity(physical_arity)?
                .with_admitted_domain(admitted_domain)?;
        let targets = self
            .targets()
            .iter()
            .map(|target| {
                target
                    .iter()
                    .map(|(integral, coefficient)| {
                        if integral.0.len() != physical_arity {
                            return Err(Error::InvalidInput(
                                "weighted target physical arity".into(),
                            ));
                        }
                        // Target conversion has already been performed in the
                        // unpadded physical inverse-propagator basis.
                        let mut indices = [0_i64; N];
                        for (index, &physical) in indices.iter_mut().zip(&integral.0) {
                            *index = i64::from(physical);
                        }
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
