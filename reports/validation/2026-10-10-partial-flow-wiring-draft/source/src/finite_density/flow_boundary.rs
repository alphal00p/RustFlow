//! Automatic boundary matching for occupied, fixed-shell auxiliary mass flows.
//!
//! Every uncut physical quadratic is deformed. Compact momenta remain soft.
//! Each deformed quadratic denominator is purely hard or constant after
//! expansion; numerator-completion powers are nonpositive by default. An
//! explicit opt-in also admits inverse powers of individually certified
//! positive compact energies, which remain soft. Tensor projection therefore
//! leaves soft moments integrated against fixed shell and endpoint distributions. Any
//! unrestricted virtual soft loop is a scaleless polynomial integral.
//!
//! The hard vacuum factors have one large scale and use ordinary recursive
//! AMF. At generic symbolic epsilon their eta dependence is a pure power, so
//! the region coefficients have no explicit logarithms. Expanding the retained
//! epsilon-dependent exponents later can generate logarithms. Partial shifts
//! and retained nonpolynomial soft factors need separate weighted recursion.
use super::boundary::{
    IntegratedOccupiedBoundary, OccupiedBoundaryDistribution, OccupiedBoundaryLimits,
};
use super::compact::CompactShell;
use super::geometry::OccupiedCutFamily;
use crate::boundary::LogRegionBoundary;
use crate::frobenius::FrobeniusBasis;
use crate::{
    ComplexFloat as C, Error, FlowOptions, Integral, Precision, Progress, Result, RunContext,
    RustRedBackend, integrand, regions,
};
use ahash::HashMap;
use symbolica::prelude::*;

/// Region routing and the actual expansion depth used for each master.
#[derive(Clone, Debug)]
pub struct OccupiedRegionProvenance {
    pub hard: Vec<bool>,
    pub transformation: Vec<Vec<Atom>>,
    pub jacobian_determinant: Atom,
    /// Symbolic z=1/eta exponents, in the supplied master order.
    pub exponents: Vec<Atom>,
    /// None means the rank planner did not require this region coefficient.
    pub half_orders: Vec<Option<usize>>,
}

#[derive(Clone, Debug, Default)]
pub struct OccupiedBoundaryProvenance {
    pub regions: Vec<OccupiedRegionProvenance>,
    pub integrated_coefficients: usize,
    pub integrated_products: usize,
    pub scaleless_noncompact_products: usize,
    pub vanishing_required_cut_coefficients: usize,
    pub massless_origin: Option<&'static str>,
    pub partial_origin: Option<String>,
    pub certified_virtual_soft_products: usize,
    pub virtual_soft_certificates: Vec<serde_json::Value>,
}

#[derive(Clone, Debug)]
pub struct OccupiedBoundaryConstants {
    pub constants: Vec<C>,
    pub provenance: OccupiedBoundaryProvenance,
}

/// Boundary owner for a closed weighted master basis in native normalization.
/// Target Wick phases, the input routing Jacobian and the whole-amplitude
/// Euclidean normalization belong to the caller, and are not applied here.
pub struct OccupiedFlowBoundary<'a> {
    integrated: IntegratedOccupiedBoundary<'a>,
    options: &'a FlowOptions,
    context: &'a RunContext,
    epsilon_symbol: Symbol,
    max_half_order: usize,
    limits: OccupiedBoundaryLimits,
    positive_energy_powers: bool,
    massless_evidence: Option<super::massless_endpoint::MasslessFlowEvidence>,
    partial_continuation: Option<super::partial_origin::continuation::PartialContinuation>,
}

impl<'a> OccupiedFlowBoundary<'a> {
    pub fn new(
        backend: &'a RustRedBackend,
        options: &'a FlowOptions,
        context: &'a RunContext,
        epsilon_symbol: Symbol,
        max_half_order: usize,
        limits: OccupiedBoundaryLimits,
    ) -> Result<Self> {
        if max_half_order > 100 {
            return Err(Error::Limit(
                "occupied boundary half-order exceeds the integrand limit of 100".into(),
            ));
        }
        Ok(Self {
            integrated: IntegratedOccupiedBoundary::new(backend, options, context, limits)?,
            options,
            context,
            epsilon_symbol,
            max_half_order,
            limits,
            positive_energy_powers: false,
            massless_evidence: None,
            partial_continuation: None,
        })
    }

    /// Opt into positive completion indices only when geometry proves that
    /// the factor is a nonzero rational multiple of one future compact energy
    /// with assigned positive shell mass. The default remains polynomial.
    pub fn with_positive_compact_energy_powers(mut self, enabled: bool) -> Self {
        self.positive_energy_powers = enabled;
        self.integrated = self.integrated.with_positive_compact_energy_powers(enabled);
        self
    }

    /// Retain the sealed proof. `constants` rechecks its exact family and
    /// deformation and audits every supplied master before evaluating a seed.
    pub fn with_massless_evidence(
        mut self,
        evidence: &super::massless_endpoint::MasslessFlowEvidence,
    ) -> Result<Self> {
        if self.positive_energy_powers || self.partial_continuation.is_some() {
            return Err(Error::Unsupported(
                "massless boundary evidence excludes inverse compact energies".into(),
            ));
        }
        self.integrated = self.integrated.with_massless_flow_origin(evidence);
        self.massless_evidence = Some(evidence.clone());
        Ok(self)
    }

    pub(crate) fn with_partial_continuation(
        mut self,
        evidence: &super::partial_origin::continuation::PartialContinuation,
    ) -> Result<Self> {
        if self.positive_energy_powers || self.massless_evidence.is_some() {
            return Err(Error::Unsupported(
                "partial boundary requires its distinct polynomial-completion proof".into(),
            ));
        }
        self.integrated = self.integrated.with_partial_flow_origin(evidence);
        self.partial_continuation = Some(evidence.clone());
        Ok(self)
    }

    /// Derive the needed region coefficients and match the supplied complete
    /// Frobenius basis. No numerical specialization of eta exponents is used.
    ///
    /// `basis` contains all physical, numerator-completion and occupation
    /// indices. `shifted` may contain either the input slots or all weighted
    /// slots; any appended occupation shifts must be false. All uncut physical
    /// slots must be shifted, and every cut, occupation and completion slot
    /// must be fixed. Shell masses must have exact physical assignments. Positive
    /// completion indices require the explicit positive-energy opt-in and a
    /// matching exact geometry certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn constants(
        &self,
        family: &OccupiedCutFamily,
        basis: &[Integral],
        shifted: &[bool],
        epsilon: &Rational,
        p: Precision,
        solutions: &FrobeniusBasis,
    ) -> Result<OccupiedBoundaryConstants> {
        self.context.cancellation.check()?;
        let slots = family.factors().len();
        let input_slots = family.input_slots();
        if basis.is_empty()
            || basis.len() != solutions.columns.len()
            || p.bits != solutions.precision.bits
            || basis.iter().any(|i| i.0.len() != slots)
            || ![input_slots, slots].contains(&shifted.len())
        {
            return Err(Error::InvalidInput(
                "occupied boundary master or deformation dimensions".into(),
            ));
        }
        let mut required_cuts = vec![false; input_slots];
        let mut compact = vec![false; family.loops()];
        for shell in family.shells() {
            required_cuts[shell.physical_slot] = true;
            compact[shell.loop_index] = true;
        }
        if basis.iter().any(|master| {
            (family.physical_slots()..input_slots).any(|slot| {
                master.0[slot] > 0
                    && !(self.positive_energy_powers
                        && family.compact_energy_completion(slot).is_some())
            }) || family
                .shells()
                .iter()
                .any(|shell| master.0[shell.upper_slot] < 0 || master.0[shell.lower_slot] < 0)
        }) {
            return Err(Error::InvalidInput(
                "occupied boundary requires polynomial or explicitly certified positive-energy completion indices and nonnegative occupation indices".into(),
            ));
        }
        for (slot, &deformed) in shifted.iter().enumerate() {
            let expected = slot < family.physical_slots() && !required_cuts[slot];
            if deformed && !expected || self.partial_continuation.is_none() && deformed != expected
            {
                return Err(Error::Unsupported(
                    "occupied boundaries currently require every uncut physical quadratic, and only those quadratics, to be shifted".into(),
                ));
            }
        }
        if !shifted.iter().any(|&deformed| deformed) {
            return Err(Error::Unsupported(
                "an occupied family without uncut physical quadratics is a terminal moment, not an auxiliary mass flow".into(),
            ));
        }
        if let Some(evidence) = &self.massless_evidence {
            if self.positive_energy_powers {
                return Err(Error::Unsupported(
                    "massless boundary evidence excludes inverse compact energies".into(),
                ));
            }
            let selected = shifted
                .iter()
                .enumerate()
                .filter_map(|(slot, &selected)| selected.then_some(slot))
                .collect::<Vec<_>>();
            evidence.validate_family(family, &selected)?;
            evidence.validate_labels(family, basis)?;
        }
        if let Some(evidence) = &self.partial_continuation {
            if self.positive_energy_powers || self.massless_evidence.is_some() {
                return Err(Error::Unsupported(
                    "conflicting partial boundary authority".into(),
                ));
            }
            let selected = shifted
                .iter()
                .enumerate()
                .filter_map(|(slot, &on)| on.then_some(slot))
                .collect::<Vec<_>>();
            evidence.validate_family_binding(
                family,
                &selected,
                evidence.source_options(),
                self.context,
            )?;
            evidence.audit_labels(family, basis, self.context)?;
            crate::physical_conditions::validate_conditions_at(
                evidence.raw_nonzero_conditions(),
                self.epsilon_symbol,
                &std::collections::BTreeMap::from([(
                    self.epsilon_symbol,
                    Atom::num(epsilon.clone()),
                )]),
            )?;
        }
        let ordinary = family.region_family(self.epsilon_symbol, self.options.dimension)?;
        let shifted = &shifted[..input_slots];
        // Ordinary region expansion keeps its numerator-only completion
        // contract. Certified inverse energies are fixed compact soft factors:
        // remove them before expansion, and restore their exact factors below.
        let fixed_energies = basis
            .iter()
            .map(|master| {
                (family.physical_slots()..input_slots)
                    .filter(|&slot| master.0[slot] > 0)
                    .map(|slot| {
                        Ok((
                            family.compact_energy_completion(slot).ok_or_else(|| {
                                Error::InvalidInput("missing compact-energy certificate".into())
                            })?,
                            master.0[slot],
                        ))
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let targets = basis
            .iter()
            .map(|master| {
                Integral(
                    master.0[..input_slots]
                        .iter()
                        .enumerate()
                        .map(|(slot, &index)| {
                            if required_cuts[slot] || (slot >= family.physical_slots() && index > 0)
                            {
                                0
                            } else {
                                index
                            }
                        })
                        .collect(),
                )
            })
            .collect::<Vec<_>>();
        let distributions = basis
            .iter()
            .map(|master| {
                family.shells().iter().map(|shell| {
                    let mass_squared = Rational::try_from(shell.mass_squared.as_view()).map_err(|_| {
                        Error::Unsupported("occupied boundary shell masses need exact physical assignments".into())
                    })?;
                    Ok(OccupiedBoundaryDistribution {
                        source_loop_index: shell.loop_index,
                        shell: CompactShell { mass_squared, chemical_potential: shell.chemical_potential.clone() },
                        cut_index: master.0[shell.physical_slot],
                        upper_index: master.0[shell.upper_slot],
                        lower_index: master.0[shell.lower_slot],
                    })
                }).collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let regions = crate::cut_regions::enumerate_compact(&ordinary, &compact, self.context)?;
        // Occupied loops stay soft: every ordinary hard child has fewer loops
        // than this weighted parent. The existing ordinary owner then performs
        // its own massive-line recursion with the tadpole-only terminal policy.
        let transformed = regions
            .iter()
            .map(|region| ordinary.transform_loops(&region.transformation))
            .collect::<Result<Vec<_>>>()?;
        let energy_offset = family.loops() * (family.loops() + 1) / 2;
        for (region, transformed) in regions.iter().zip(&transformed) {
            for (certificate, _) in fixed_energies.iter().flatten() {
                let factor = &transformed.propagators[certificate.slot];
                if region.hard[certificate.loop_index]
                    || !factor.constant.is_zero()
                    || factor.scalar_products.iter().enumerate().any(|(i, c)| {
                        let expected = if i == energy_offset + certificate.loop_index {
                            Atom::num(certificate.coefficient.clone())
                        } else {
                            Atom::zero()
                        };
                        !(c - expected).together().cancel().is_zero()
                    })
                {
                    return Err(Error::Unsupported(
                        "certified compact energy does not remain fixed in the region".into(),
                    ));
                }
            }
        }
        // The checked fixed factors have eta degree zero, so stripping them
        // changes neither the symbolic exponent nor the required series depth.
        let powers = targets
            .iter()
            .map(|target| {
                regions
                    .iter()
                    .map(|region| {
                        self.context.cancellation.check()?;
                        Ok(
                            -regions::expand_region(&ordinary, target, shifted, region, 0)?
                                .eta_power,
                        )
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let orders = solutions.region_orders(&powers, self.max_half_order)?;
        let coefficient_count = orders
            .iter()
            .flatten()
            .flatten()
            .map(|order| order + 1)
            .sum::<usize>();
        if coefficient_count > self.limits.max_terms {
            return Err(Error::Limit(
                "occupied boundary coefficient budget exhausted".into(),
            ));
        }
        self.context.emit(Progress::BoundaryPlan {
            basis_size: basis.len(),
            region_series: orders
                .iter()
                .flatten()
                .filter(|order| order.is_some())
                .count(),
            coefficients: coefficient_count,
            max_half_order: orders
                .iter()
                .flatten()
                .flatten()
                .copied()
                .max()
                .unwrap_or(0),
        })?;
        let mut provenance = OccupiedBoundaryProvenance {
            massless_origin: self.massless_evidence.as_ref().map(|e| e.origin_identity()),
            partial_origin: self
                .partial_continuation
                .as_ref()
                .map(|e| e.origin_identity().to_owned()),
            regions: regions
                .iter()
                .enumerate()
                .map(|(r, region)| OccupiedRegionProvenance {
                    hard: region.hard.clone(),
                    transformation: region.transformation.clone(),
                    jacobian_determinant: region.jacobian_determinant.clone(),
                    exponents: powers.iter().map(|row| row[r].clone()).collect(),
                    half_orders: orders.iter().map(|row| row[r]).collect(),
                })
                .collect(),
            ..Default::default()
        };
        let parameters =
            HashMap::from_iter([(Atom::var(self.epsilon_symbol), p.rational(epsilon))]);
        let dimension = Rational::from(self.options.dimension) - epsilon * &Rational::from(2);
        let mut data = vec![Vec::new(); basis.len()];
        for (component, target) in targets.iter().enumerate() {
            for (r, region) in regions.iter().enumerate() {
                self.context.cancellation.check()?;
                let order = orders[component][r];
                for parity in 0..2 {
                    data[component].push(LogRegionBoundary {
                        exponent: (&powers[component][r] + Atom::num((parity as i64, 2)))
                            .together()
                            .cancel(),
                        coefficients: vec![
                            vec![None];
                            order.map_or(0, |order| (order + 2 - parity) / 2)
                        ],
                        max_log_power: Some(0),
                    });
                }
                let Some(order) = order else {
                    continue;
                };
                let expansion = regions::expand_region(&ordinary, target, shifted, region, order)?;
                let fixed_weight = fixed_energies[component].iter().fold(
                    Atom::one(),
                    |weight, (certificate, index)| {
                        weight
                            * (Atom::num(certificate.coefficient.clone())
                                * &expansion.coordinates[energy_offset + certificate.loop_index])
                                .pow(-i64::from(*index))
                    },
                );
                let determinant = p.eval(&expansion.jacobian_determinant, &parameters)?;
                let jacobian = p.pow(
                    &C::new(p.norm(&determinant), p.real(0)),
                    &p.rational(&dimension),
                );
                for (index, expression) in expansion.coefficients.iter().enumerate() {
                    self.context.cancellation.check()?;
                    let raw_expression = expression * &fixed_weight;
                    // Capture domains using the full parent coordinate system,
                    // before tensor projection or whole-expression cancellation.
                    let inherited = if self.partial_continuation.is_some() {
                        Some(super::boundary::partial_region_conditions(
                            &raw_expression,
                            &expansion.coordinates,
                            &compact,
                            transformed[r].external.len(),
                            self.limits,
                            self.context,
                        )?)
                    } else {
                        None
                    };
                    let projected = integrand::projected_factor_region(
                        &raw_expression.together().cancel(),
                        &expansion.coordinates,
                        &transformed[r],
                        &region.hard,
                        self.limits.max_terms,
                    )?;
                    let report = self
                        .integrated
                        .evaluate_projected_with_inherited_conditions(
                            &projected,
                            &distributions[component],
                            epsilon,
                            &parameters,
                            p,
                            inherited.as_deref(),
                        )?;
                    provenance.integrated_coefficients += 1;
                    provenance.integrated_products += report.integrated_products;
                    provenance.scaleless_noncompact_products +=
                        report.scaleless_noncompact_products;
                    provenance.vanishing_required_cut_coefficients +=
                        usize::from(report.vanishing_required_cut);
                    provenance.certified_virtual_soft_products +=
                        report.virtual_soft_certificates.len();
                    provenance
                        .virtual_soft_certificates
                        .extend(report.virtual_soft_certificates);
                    data[component][2 * r + index % 2].coefficients[index / 2][0] =
                        Some(&jacobian * &report.value);
                }
            }
        }
        Ok(OccupiedBoundaryConstants {
            constants: solutions.match_log_regions(&data)?,
            provenance,
        })
    }
}
