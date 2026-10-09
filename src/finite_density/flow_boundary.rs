//! Automatic boundary matching for occupied, fixed-shell auxiliary mass flows.
//!
//! Every uncut physical quadratic is deformed. Compact momenta remain soft,
//! and the admitted regions have at most one hard virtual loop. Their exact
//! integrand expansion is log-free at symbolic epsilon: the remaining hard
//! integrals have a single large scale, while the compact factors are
//! polynomials integrated against fixed shell and endpoint distributions.
//! More general soft denominators and multiloop hard regions need additional
//! weighted recursion and are rejected explicitly.
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
        })
    }

    /// Derive the needed region coefficients and match the supplied complete
    /// Frobenius basis. No numerical specialization of eta exponents is used.
    ///
    /// `basis` contains all physical, numerator-completion and occupation
    /// indices. `shifted` may contain either the input slots or all weighted
    /// slots; any appended occupation shifts must be false. All uncut physical
    /// slots must be shifted, and every cut, occupation and completion slot
    /// must be fixed. Shell masses must have exact physical assignments.
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
            master.0[family.physical_slots()..input_slots]
                .iter()
                .any(|&index| index > 0)
                || family
                    .shells()
                    .iter()
                    .any(|shell| master.0[shell.upper_slot] < 0 || master.0[shell.lower_slot] < 0)
        }) {
            return Err(Error::InvalidInput(
                "occupied boundary requires nonpositive numerator-completion indices and nonnegative occupation indices".into(),
            ));
        }
        for (slot, &deformed) in shifted.iter().enumerate() {
            let expected = slot < family.physical_slots() && !required_cuts[slot];
            if deformed != expected {
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
        let ordinary = family.region_family(self.epsilon_symbol, self.options.dimension)?;
        let shifted = &shifted[..input_slots];
        let targets = basis
            .iter()
            .map(|master| {
                Integral(
                    master.0[..input_slots]
                        .iter()
                        .enumerate()
                        .map(|(slot, &index)| if required_cuts[slot] { 0 } else { index })
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
        if regions
            .iter()
            .any(|region| region.hard.iter().filter(|&&hard| hard).count() > 1)
        {
            return Err(Error::Unsupported(
                "occupied flow boundaries currently admit at most one hard virtual loop per region"
                    .into(),
            ));
        }
        let transformed = regions
            .iter()
            .map(|region| ordinary.transform_loops(&region.transformation))
            .collect::<Result<Vec<_>>>()?;
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
                let determinant = p.eval(&expansion.jacobian_determinant, &parameters)?;
                let jacobian = p.pow(
                    &C::new(p.norm(&determinant), p.real(0)),
                    &p.rational(&dimension),
                );
                for (index, expression) in expansion.coefficients.iter().enumerate() {
                    self.context.cancellation.check()?;
                    let projected = integrand::projected_factor_region(
                        expression,
                        &expansion.coordinates,
                        &transformed[r],
                        &region.hard,
                        self.limits.max_terms,
                    )?;
                    let report = self.integrated.evaluate_projected_with_provenance(
                        &projected,
                        &distributions[component],
                        epsilon,
                        &parameters,
                        p,
                    )?;
                    provenance.integrated_coefficients += 1;
                    provenance.integrated_products += report.integrated_products;
                    provenance.scaleless_noncompact_products +=
                        report.scaleless_noncompact_products;
                    provenance.vanishing_required_cut_coefficients +=
                        usize::from(report.vanishing_required_cut);
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
