//! Direct occupied polynomial moments when no virtual pole survives.
//!
//! This terminal uses the same independent shell and endpoint distributions as
//! occupied AMF boundaries. No auxiliary connection or asymptotic matching is
//! introduced for a compact polynomial integral.
use super::PreparedDensityInput;
use super::boundary::{
    CompactOriginPrescription, IntegratedOccupiedBoundary, OccupiedBoundaryDistribution,
    OccupiedBoundaryLimits,
};
use super::compact::CompactShell;
use super::flow::OccupiedFlowEvaluation;
use super::flow_boundary::OccupiedBoundaryProvenance;
use super::geometry::OccupiedCutFamily;
use super::guarded::IndexRole;
use super::normalization::native_measure_to_euclidean;
use crate::{Error, FlowOptions, Precision, Prescription, Result, RunContext, RustRedBackend};
use ahash::HashMap;
use std::collections::{BTreeMap, BTreeSet};
use symbolica::prelude::*;

pub struct PreparedOccupiedTerminal {
    family: OccupiedCutFamily,
    dimension: i64,
    epsilon: Symbol,
    nonzero_conditions: Vec<Atom>,
    origin_prescription: CompactOriginPrescription,
}

impl PreparedOccupiedTerminal {
    /// Admit precisely the surviving cut terms whose noncut factors are
    /// polynomials. Other occupied sectors must use the weighted flow owner.
    pub fn prepare(
        input: &PreparedDensityInput,
        cuts: &[usize],
        options: &FlowOptions,
    ) -> Result<Option<Self>> {
        options.validate()?;
        if options.prescription != Prescription::PlusI0 {
            return Err(Error::Unsupported(
                "occupied polynomial terminals require the common +i0 convention".into(),
            ));
        }
        let family = input.occupied_cut(cuts, 65536)?.at_physical_masses();
        if family
            .targets()
            .iter()
            .flat_map(|target| target.keys())
            .any(|target| {
                (0..family.input_slots())
                    .any(|slot| family.roles()[slot] == IndexRole::Ordinary && target.0[slot] > 0)
            })
        {
            return Ok(None);
        }
        let mut origin_prescription = CompactOriginPrescription::Existing;
        for shell in family.shells() {
            let mass = Rational::try_from(shell.mass_squared.as_view()).map_err(|_| {
                Error::Unsupported(
                    "occupied terminal requires assigned rational shell masses".into(),
                )
            })?;
            if mass < 0 {
                return Err(Error::Unsupported(
                    "negative occupied shell mass squared".into(),
                ));
            }
            if mass.is_zero() {
                if shell.chemical_potential <= 0 {
                    return Err(Error::Unsupported("massless polynomial terminal requires mu>0; coincident zero-energy endpoints are not admitted".into()));
                }
                origin_prescription = CompactOriginPrescription::JointDimensionalMasslessTerminal;
            }
        }
        // Keep denominator domains even when the cut removes a target term or
        // a physical assignment cancels a coefficient. Include the original
        // polynomial spelling before its inverse-propagator conversion.
        let mut expressions = input
            .targets()
            .iter()
            .flat_map(|t| t.values().cloned())
            .collect::<Vec<_>>();
        for target in &input.input().targets {
            expressions.push(
                Atom::parse(&target.numerator, "rustflow_density", Default::default())
                    .map_err(|error| Error::InvalidInput(error.to_string()))?,
            );
        }
        let variables = expressions
            .iter()
            .flat_map(|a| a.get_all_symbols(true))
            .collect::<BTreeSet<_>>();
        let assignments = input
            .independent_masses()
            .iter()
            .zip(input.physical_masses())
            .map(|(&mass, value)| (Atom::var(mass), value.clone()))
            .collect::<BTreeMap<_, _>>();
        let nonzero_conditions =
            crate::physical_conditions::rational_denominator_conditions(&expressions, &variables)?
                .iter()
                .map(|condition| crate::family::substitute(condition, &assignments))
                .collect();
        Ok(Some(Self {
            family,
            dimension: options.dimension,
            epsilon: symbol!("rustflow_occupied::epsilon"),
            nonzero_conditions,
            origin_prescription,
        }))
    }

    pub fn family(&self) -> &OccupiedCutFamily {
        &self.family
    }

    /// Measure-origin identity retained in seed caches and evaluation reports.
    pub fn origin_prescription(&self) -> &'static str {
        self.origin_prescription.identity()
    }

    /// Euclidean contribution with target Wick coefficients and the routing
    /// determinant applied once. Threshold values and singular raised limits
    /// are decided by the compact distribution owner at the supplied dimension.
    pub fn evaluate(
        &self,
        epsilon: &Rational,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<OccupiedFlowEvaluation> {
        options.validate()?;
        context.cancellation.check()?;
        if options.dimension != self.dimension
            || epsilon.is_zero()
            || options.prescription != Prescription::PlusI0
        {
            return Err(Error::InvalidInput(
                "occupied terminal needs nonzero epsilon, the prepared dimension and common +i0 convention".into(),
            ));
        }
        crate::physical_conditions::validate_conditions_at(
            &self.nonzero_conditions,
            self.epsilon,
            &BTreeMap::from([(self.epsilon, Atom::num(epsilon.clone()))]),
        )?;
        let p = Precision::decimal(
            options
                .digits
                .checked_add(options.guard_digits)
                .ok_or_else(|| Error::Limit("occupied terminal precision overflow".into()))?,
        )?;
        let parameters = HashMap::from_iter([(Atom::var(self.epsilon), p.rational(epsilon))]);
        let backend = RustRedBackend {
            bubble_subloops: false,
            ..Default::default()
        };
        let limits = OccupiedBoundaryLimits::default();
        let boundary = IntegratedOccupiedBoundary::new(&backend, options, context, limits)?
            .with_terminal_origin(self.origin_prescription);
        let ordinary = self.family.region_family(self.epsilon, self.dimension)?;
        let hard = vec![false; self.family.loops()];
        let mut provenance = OccupiedBoundaryProvenance::default();
        let mut values = Vec::with_capacity(self.family.targets().len());
        for target in self.family.targets() {
            let mut value = p.zero();
            for (integral, coefficient) in target {
                context.cancellation.check()?;
                let polynomial = (0..self.family.input_slots())
                    .filter(|&slot| self.family.roles()[slot] == IndexRole::Ordinary)
                    .fold(Atom::one(), |product, slot| {
                        product * self.family.factors()[slot].pow(-i64::from(integral.0[slot]))
                    });
                let projected = crate::integrand::projected_factor_region(
                    &polynomial,
                    self.family.coordinates(),
                    &ordinary,
                    &hard,
                    limits.max_terms,
                )?;
                let distributions = self.family.shells().iter().map(|shell| {
                    Ok(OccupiedBoundaryDistribution {
                        source_loop_index: shell.loop_index,
                        shell: CompactShell {
                            mass_squared: Rational::try_from(shell.mass_squared.as_view())
                                .map_err(|_| Error::Unsupported("occupied terminal needs assigned rational shell masses".into()))?,
                            chemical_potential: shell.chemical_potential.clone(),
                        },
                        cut_index: integral.0[shell.physical_slot],
                        upper_index: integral.0[shell.upper_slot],
                        lower_index: integral.0[shell.lower_slot],
                    })
                }).collect::<Result<Vec<_>>>()?;
                let report = boundary.evaluate_projected_with_provenance(
                    &projected,
                    &distributions,
                    epsilon,
                    &parameters,
                    p,
                )?;
                provenance.integrated_coefficients += 1;
                provenance.integrated_products += report.integrated_products;
                provenance.scaleless_noncompact_products += report.scaleless_noncompact_products;
                provenance.vanishing_required_cut_coefficients +=
                    usize::from(report.vanishing_required_cut);
                value = p.add(
                    &value,
                    &p.mul(&p.eval(coefficient, &parameters)?, &report.value),
                );
            }
            values.push(value);
        }
        let dimension = Rational::from(self.dimension) - epsilon * &Rational::from(2);
        let determinant = self.family.routing_determinant();
        let absolute = if determinant < &Rational::zero() {
            -determinant.clone()
        } else {
            determinant.clone()
        };
        let measure = native_measure_to_euclidean(
            self.family.loops(),
            self.family.shells().len(),
            &Atom::num(dimension.clone()),
        )? * Atom::num(absolute).pow(Atom::num(-dimension));
        let factor = p.eval(&measure, &HashMap::default())?;
        Ok(OccupiedFlowEvaluation {
            values: values.iter().map(|value| p.mul(value, &factor)).collect(),
            boundary: provenance,
            nonzero_conditions: self.nonzero_conditions.clone(),
            contour_admission: format!(
                "independent occupied shells; polynomial remaining factors have no virtual poles; unrestricted polynomial loop integrals are scaleless in dimensional regularization; compact origin prescription={}; massive threshold products retain their thermal limits; massless zero-origin products are jointly continued from sufficiently large ReD before specializing dimension; no flowing massless endpoint admission",
                self.origin_prescription.identity()
            ),
            shifted_slots: vec![],
            basis_size: 0,
            construction: "compact_polynomial_moments",
            source_options: None,
        })
    }
}
