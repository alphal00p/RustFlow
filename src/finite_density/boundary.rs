//! Integrated coefficients of fixed-shell large-mass regions.
//!
//! Hard factors use ordinary recursive AMF with Gaussian/tadpole seeds. Soft
//! factors must be polynomials in the independent compact momenta, optionally
//! with inverse powers of individual positive-mass compact energies. HEPKit
//! projects their spatial angular dependence before the
//! one-shell distribution owner integrates radial and energy monomials.
//!
//! The returned normalization is native throughout: one virtual loop uses
//! `d^D k/(i*pi^(D/2))`, and one occupied loop uses
//! `pi^(-D/2) d^D q C_n(q^2-m^2) H_upper(mu-q0) H_lower(q0)`.
//! Input Wick phases and the whole-amplitude Euclidean measure are applied by
//! the caller. This owner neither moves a shell mass with eta nor discards a
//! retained virtual pole in a soft factor.
use super::compact::CompactShell;
use crate::coefficient::{exact_coefficient_list, powers};
use crate::family::substitute;
use crate::integrand::{ProjectedRegion, RegionFactorSpace};
use crate::recursive::{RecursiveBoundary, RecursiveTerminalPolicy};
use crate::{ComplexFloat as C, Error, FlowOptions, Precision, Result, RunContext, RustRedBackend};
use ahash::HashMap;
use std::collections::BTreeMap;
use std::sync::Mutex;
use symbolica::prelude::*;

/// This additional origin prescription is available only to the polynomial
/// terminal owner. It is not a continuation certificate for a flowing graph.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CompactOriginPrescription {
    Existing,
    JointDimensionalMasslessTerminal,
}

impl CompactOriginPrescription {
    pub(crate) fn identity(self) -> &'static str {
        match self {
            Self::Existing => "existing-compact-support-v1",
            Self::JointDimensionalMasslessTerminal => {
                "joint-high-D-massless-polynomial-origin-v1;mu>0;lower-contact-zero-jets;no-virtual-poles"
            }
        }
    }
}

/// The distribution carried by one compact loop in the region's supplied
/// routing. Indices refer to C_n and H_n, not ordinary denominator powers.
#[derive(Clone, Debug)]
pub struct OccupiedBoundaryDistribution {
    pub source_loop_index: usize,
    pub shell: CompactShell,
    pub cut_index: i16,
    pub upper_index: i16,
    pub lower_index: i16,
}

#[derive(Clone, Copy, Debug)]
pub struct OccupiedBoundaryLimits {
    pub max_terms: usize,
    pub max_compact_series_terms: usize,
    pub max_distribution_order: usize,
}

impl Default for OccupiedBoundaryLimits {
    fn default() -> Self {
        Self {
            max_terms: 10_000,
            max_compact_series_terms: 10_000,
            max_distribution_order: 32,
        }
    }
}

pub struct IntegratedOccupiedBoundary<'a> {
    ordinary: RecursiveBoundary<'a>,
    context: &'a RunContext,
    limits: OccupiedBoundaryLimits,
    moments: Mutex<BTreeMap<String, C>>,
    positive_energy_powers: bool,
    origin_prescription: CompactOriginPrescription,
}

/// Provenance of a single integrated region coefficient. The scaleless count
/// refers only to unrestricted, noncompact polynomial loop integrals.
#[derive(Clone, Debug)]
pub struct IntegratedOccupiedCoefficient {
    pub value: C,
    pub integrated_products: usize,
    pub scaleless_noncompact_products: usize,
    pub vanishing_required_cut: bool,
}

impl<'a> IntegratedOccupiedBoundary<'a> {
    pub fn new(
        backend: &'a RustRedBackend,
        options: &'a FlowOptions,
        context: &'a RunContext,
        limits: OccupiedBoundaryLimits,
    ) -> Result<Self> {
        if backend.bubble_subloops {
            return Err(Error::InvalidInput(
                "occupied boundaries require bubble_subloops=false".into(),
            ));
        }
        // Compact beta sums reserve eight working digits for accumulated
        // arithmetic. Their truncation target must still cover the requested
        // digits when this owner is called by the common flow driver.
        if options.guard_digits < 8 {
            return Err(Error::InvalidInput(
                "occupied boundaries require at least eight guard digits for compact moment summation"
                    .into(),
            ));
        }
        if limits.max_terms == 0
            || limits.max_compact_series_terms == 0
            || limits.max_distribution_order == 0
        {
            return Err(Error::InvalidInput(
                "zero occupied boundary work limit".into(),
            ));
        }
        Ok(Self {
            ordinary: RecursiveBoundary::new(backend, options, context)
                .with_terminal_policy(RecursiveTerminalPolicy::TadpolesOnly),
            context,
            limits,
            moments: Mutex::new(BTreeMap::new()),
            positive_energy_powers: false,
            origin_prescription: CompactOriginPrescription::Existing,
        })
    }

    /// Opt into Laurent monomials in individual future compact energies.
    /// Each negative power still requires the assigned positive-mass proof
    /// E>=m>0. Other soft denominators keep their existing rejection.
    pub fn with_positive_compact_energy_powers(mut self, enabled: bool) -> Self {
        self.positive_energy_powers = enabled;
        self
    }

    pub(crate) fn with_terminal_origin(mut self, prescription: CompactOriginPrescription) -> Self {
        self.origin_prescription = prescription;
        self
    }

    /// Integrate a scalar region coefficient after hard tensor projection.
    /// The soft space contains only compact loops; if it has one external
    /// vector, that vector is explicitly the normalized medium `u.u=1`.
    /// An additional noncompact soft loop is a scaleless zero only if the
    /// coefficient is polynomial in every coordinate involving it and its
    /// measure is unrestricted. Other surviving soft denominators need
    /// recursive weighted reduction and are rejected here.
    ///
    /// This integrates one coefficient only. The region's symbolic eta
    /// exponent and logarithmic order remain with the caller and are matched
    /// using `FrobeniusBasis::match_log_regions`, without epsilon specialization.
    pub fn evaluate_projected(
        &self,
        projected: &ProjectedRegion,
        distributions: &[OccupiedBoundaryDistribution],
        epsilon: &Rational,
        parameters: &HashMap<Atom, C>,
        p: Precision,
    ) -> Result<C> {
        Ok(self
            .evaluate_projected_with_provenance(projected, distributions, epsilon, parameters, p)?
            .value)
    }

    pub fn evaluate_projected_with_provenance(
        &self,
        projected: &ProjectedRegion,
        distributions: &[OccupiedBoundaryDistribution],
        epsilon: &Rational,
        parameters: &HashMap<Atom, C>,
        p: Precision,
    ) -> Result<IntegratedOccupiedCoefficient> {
        self.context.cancellation.check()?;
        if projected.terms.len() > self.limits.max_terms {
            return Err(Error::Limit(
                "occupied boundary product budget exhausted".into(),
            ));
        }
        if self.origin_prescription == CompactOriginPrescription::JointDimensionalMasslessTerminal
            && !projected.hard.template.loops.is_empty()
        {
            return Err(Error::Unsupported(
                "joint massless origin prescription admits polynomial terminals only, not virtual hard factors".into(),
            ));
        }
        let mut by_loop = BTreeMap::new();
        for distribution in distributions {
            if by_loop
                .insert(distribution.source_loop_index, distribution)
                .is_some()
            {
                return Err(Error::InvalidInput(
                    "repeated occupied boundary loop".into(),
                ));
            }
            validate_distribution(distribution, self.limits)?;
            if self.origin_prescription
                == CompactOriginPrescription::JointDimensionalMasslessTerminal
                && distribution.shell.mass_squared.is_zero()
                && distribution.shell.chemical_potential <= 0
            {
                return Err(Error::Unsupported(
                    "joint massless terminal origin requires positive chemical potential; coincident endpoints at mu=0 are not admitted".into(),
                ));
            }
        }
        if projected
            .hard
            .source_loop_indices
            .iter()
            .any(|i| by_loop.contains_key(i))
        {
            return Err(Error::InvalidInput(
                "occupied momentum cannot enter a hard region".into(),
            ));
        }
        let soft_distributions = projected
            .soft
            .source_loop_indices
            .iter()
            .filter_map(|i| by_loop.get(i).copied())
            .collect::<Vec<_>>();
        if soft_distributions.len() != by_loop.len() {
            return Err(Error::InvalidInput(
                "occupied loop missing from projected region".into(),
            ));
        }
        if distributions.iter().any(|d| d.cut_index <= 0) {
            return Ok(IntegratedOccupiedCoefficient {
                value: p.zero(),
                integrated_products: 0,
                scaleless_noncompact_products: 0,
                vanishing_required_cut: true,
            });
        }
        let noncompact_coordinates = noncompact_coordinates(&projected.soft, &by_loop)?;
        let dimension =
            Rational::from(projected.soft.template.dimension) - epsilon * &Rational::from(2);
        if projected.hard.template.dimension != projected.soft.template.dimension
            || projected.hard.template.epsilon != projected.soft.template.epsilon
        {
            return Err(Error::InvalidInput(
                "region sides have different dimensions".into(),
            ));
        }
        let mut parameters = parameters.clone();
        parameters.insert(
            Atom::var(projected.soft.template.epsilon),
            p.rational(epsilon),
        );
        let mut result = p.zero();
        let mut integrated_products = 0;
        let mut scaleless_noncompact_products = 0;
        for term in &projected.terms {
            self.context.cancellation.check()?;
            if !noncompact_coordinates.is_empty() {
                polynomial_terms(&term.soft, &noncompact_coordinates, self.limits.max_terms)
                    .map_err(|error| match error {
                        Error::Unsupported(_) => Error::Unsupported(
                            "retained nonpolynomial factor in a noncompact soft loop requires recursive weighted reduction".into(),
                        ),
                        other => other,
                    })?;
                // No distribution depends on this loop, and every occurrence
                // is polynomial. Its unrestricted DR integral is scaleless.
                // Occupied polynomial loops are never subject to this rule.
                scaleless_noncompact_products += 1;
                continue;
            }
            integrated_products += 1;
            let soft = self.integrate_polynomial(
                &term.soft,
                &projected.soft,
                &soft_distributions,
                &dimension,
                &parameters,
                p,
            )?;
            if soft == p.zero() {
                continue;
            }
            let hard = if projected.hard.template.loops.is_empty() {
                p.eval(&term.hard, &parameters)?
            } else {
                let terms = crate::integrand::to_integrals(
                    &term.hard,
                    &projected.hard.coordinates,
                    &projected.hard.template,
                    self.limits.max_terms,
                )?;
                let mut value = p.zero();
                for integral in terms {
                    self.context.cancellation.check()?;
                    let factor =
                        self.ordinary
                            .evaluate(&integral.family, &integral.integral, epsilon, p)?;
                    value = p.add(
                        &value,
                        &p.mul(&p.eval(&integral.coefficient, &parameters)?, &factor),
                    );
                }
                value
            };
            result = p.add(&result, &p.mul(&hard, &soft));
        }
        if !p.finite(&result) {
            return Err(Error::Numerical(
                "nonfinite integrated occupied boundary".into(),
            ));
        }
        Ok(IntegratedOccupiedCoefficient {
            value: result,
            integrated_products,
            scaleless_noncompact_products,
            vanishing_required_cut: false,
        })
    }

    fn integrate_polynomial(
        &self,
        expression: &Atom,
        space: &RegionFactorSpace,
        distributions: &[&OccupiedBoundaryDistribution],
        dimension: &Rational,
        parameters: &HashMap<Atom, C>,
        p: Precision,
    ) -> Result<C> {
        let loops = space.template.loops.len();
        if loops != distributions.len()
            || space.template.external.len() > 1
            || !space.template.external.is_empty()
                && space.template.external_gram != vec![vec![Atom::num(1)]]
        {
            return Err(Error::Unsupported(
                "occupied polynomial boundary needs compact loops and only a normalized medium vector".into(),
            ));
        }
        let mut pairs = Vec::new();
        for i in 0..loops {
            for j in i..loops {
                pairs.push((i, j));
            }
        }
        if space.coordinates.len() != pairs.len() + loops * space.template.external.len() {
            return Err(Error::InvalidInput(
                "occupied scalar coordinate count".into(),
            ));
        }
        // Only an individual occupied energy has the E>=m>0 certificate.
        // Negative Gram powers, inverse energy sums and hidden denominators
        // remain unsupported. Validate before introducing internal coordinates.
        let mut inverse_allowed = vec![false; space.coordinates.len()];
        if !space.template.external.is_empty() {
            for (i, distribution) in distributions.iter().enumerate() {
                inverse_allowed[pairs.len() + i] =
                    self.positive_energy_powers && distribution.shell.mass_squared > 0;
            }
        }
        laurent_terms(
            expression,
            &space.coordinates,
            &inverse_allowed,
            self.limits.max_terms,
        )?;
        if loops == 0 {
            return p.eval(expression, parameters);
        }
        let sources = std::iter::once(expression)
            .chain(space.coordinates.iter())
            .collect::<Vec<_>>();
        let energies = (0..loops)
            .map(|i| fresh_symbol(&format!("rustflow_occupied_boundary::energy_{i}"), &sources))
            .collect::<Result<Vec<_>>>()?;
        let spatial = pairs
            .iter()
            .map(|(i, j)| {
                fresh_symbol(
                    &format!("rustflow_occupied_boundary::spatial_{i}_{j}"),
                    &sources,
                )
            })
            .collect::<Result<Vec<_>>>()?;
        let mut replacements = BTreeMap::new();
        for (index, &(i, j)) in pairs.iter().enumerate() {
            replacements.insert(
                space.coordinates[index].clone(),
                &energies[i] * &energies[j] - &spatial[index],
            );
        }
        if !space.template.external.is_empty() {
            for (i, energy) in energies.iter().enumerate() {
                replacements.insert(space.coordinates[pairs.len() + i].clone(), energy.clone());
            }
        }
        let mut polynomial = substitute(expression, &replacements).expand();
        let d = dimension - &Rational::from(1);
        let mut projector = crate::tensor::TensorProjector::new(
            Atom::num(space.template.dimension - 1)
                - Atom::num(2) * Atom::var(space.template.epsilon),
        );
        let scalar = |i: usize, j: usize| {
            &spatial[pairs
                .iter()
                .position(|&pair| pair == (i.min(j), i.max(j)))
                .unwrap()]
        };
        // Each step removes every angular contraction involving one compact
        // vector. The count of vectors still carrying mixed contractions
        // strictly decreases; no recursive call can recreate an earlier one.
        for selected in 0..loops {
            self.context.cancellation.check()?;
            let mut next = Atom::new();
            for (degrees, coefficient) in
                polynomial_terms(&polynomial, &spatial, self.limits.max_terms)?
            {
                let mut partners = Vec::new();
                let mut unmixed = coefficient;
                for (index, (&power, &(a, b))) in degrees.iter().zip(&pairs).enumerate() {
                    if a == selected && b > selected {
                        if partners.len() + power as usize > 32 {
                            return Err(Error::Limit("compact tensor rank exceeds 32".into()));
                        }
                        partners.extend(std::iter::repeat_n(b, power as usize));
                    } else {
                        unmixed *= spatial[index].clone().pow(i64::from(power));
                    }
                }
                let hard_gram =
                    vec![vec![scalar(selected, selected).clone(); partners.len()]; partners.len()];
                let soft_gram = partners
                    .iter()
                    .map(|&i| partners.iter().map(|&j| scalar(i, j).clone()).collect())
                    .collect::<Vec<_>>();
                next += unmixed * projector.project(&hard_gram, &soft_gram)?;
            }
            polynomial = next.expand();
        }
        let variables = energies.iter().chain(&spatial).cloned().collect::<Vec<_>>();
        let inverse_allowed = distributions
            .iter()
            .map(|d| self.positive_energy_powers && d.shell.mass_squared > 0)
            .chain(std::iter::repeat_n(false, spatial.len()))
            .collect::<Vec<_>>();
        let mut result = p.zero();
        let digits = ((u64::from(p.bits.saturating_sub(16)) * 1000) / 3322) as u32;
        let digits = digits.saturating_sub(8).max(1);
        for (degrees, coefficient) in laurent_terms(
            &polynomial,
            &variables,
            &inverse_allowed,
            self.limits.max_terms,
        )? {
            let mut value = p.eval(&coefficient, parameters)?;
            for (i, distribution) in distributions.iter().enumerate() {
                let radial_index = pairs.iter().position(|&pair| pair == (i, i)).unwrap();
                let moment = self.distribution_moment(
                    distribution,
                    &d,
                    degrees[i],
                    degrees[loops + radial_index] as u16,
                    digits,
                    p,
                )?;
                value = p.mul(&value, &moment);
            }
            if pairs
                .iter()
                .enumerate()
                .any(|(j, &(a, b))| a != b && degrees[loops + j] != 0)
            {
                return Err(Error::Numerical(
                    "unprojected compact angular contraction".into(),
                ));
            }
            result = p.add(&result, &value);
        }
        Ok(result)
    }

    fn distribution_moment(
        &self,
        distribution: &OccupiedBoundaryDistribution,
        dimension: &Rational,
        energy_power: i16,
        radial_power: u16,
        digits: u32,
        p: Precision,
    ) -> Result<C> {
        // Include support, distribution orientation, normalization, dimension,
        // requested seed accuracy and arithmetic precision. Equal-mass shells
        // can share values, but distinct physical mass derivatives cannot.
        let key = format!(
            "raw-native-Cn-Hupper-Hlower-v2:{}:{}:{}:{}:{}:{}:{energy_power}:{radial_power}:{digits}:{}:{}",
            distribution.shell.mass_squared,
            distribution.shell.chemical_potential,
            distribution.cut_index,
            distribution.upper_index,
            distribution.lower_index,
            dimension,
            p.bits,
            self.origin_prescription.identity(),
        );
        if let Some(value) = self
            .moments
            .lock()
            .map_err(|_| Error::Numerical("occupied boundary moment memo poisoned".into()))?
            .get(&key)
            .cloned()
        {
            return Ok(value);
        }
        let value = distribution_moment(
            distribution,
            dimension,
            energy_power,
            radial_power,
            digits,
            self.limits,
            self.origin_prescription,
            p,
        )?;
        if !p.finite(&value) {
            return Err(Error::Numerical(
                "nonfinite occupied distribution moment".into(),
            ));
        }
        let mut memo = self
            .moments
            .lock()
            .map_err(|_| Error::Numerical("occupied boundary moment memo poisoned".into()))?;
        if memo.len() < self.limits.max_terms {
            memo.insert(key, value.clone());
        }
        Ok(value)
    }
}

fn noncompact_coordinates(
    space: &RegionFactorSpace,
    distributions: &BTreeMap<usize, &OccupiedBoundaryDistribution>,
) -> Result<Vec<Atom>> {
    let loops = space.template.loops.len();
    let external = space.template.external.len();
    if space.source_loop_indices.len() != loops
        || space.coordinates.len() != loops * (loops + 1) / 2 + loops * external
    {
        return Err(Error::InvalidInput(
            "projected soft coordinate dimensions".into(),
        ));
    }
    let free = space
        .source_loop_indices
        .iter()
        .map(|i| !distributions.contains_key(i))
        .collect::<Vec<_>>();
    let mut coordinates = Vec::new();
    let mut index = 0;
    for i in 0..loops {
        for j in i..loops {
            if free[i] || free[j] {
                coordinates.push(space.coordinates[index].clone());
            }
            index += 1;
        }
    }
    for is_free in free {
        for _ in 0..external {
            if is_free {
                coordinates.push(space.coordinates[index].clone());
            }
            index += 1;
        }
    }
    Ok(coordinates)
}

fn validate_distribution(
    d: &OccupiedBoundaryDistribution,
    limits: OccupiedBoundaryLimits,
) -> Result<()> {
    if d.upper_index < 0
        || d.lower_index < 0
        || d.shell.mass_squared < 0
        || d.shell.chemical_potential < 0
    {
        return Err(Error::InvalidInput(
            "invalid occupied boundary distribution or support".into(),
        ));
    }
    if usize::try_from(d.cut_index.max(0)).unwrap()
        + d.upper_index as usize
        + d.lower_index as usize
        > limits.max_distribution_order
    {
        return Err(Error::Limit(
            "occupied boundary distribution order limit".into(),
        ));
    }
    Ok(())
}

fn polynomial_terms(
    expression: &Atom,
    variables: &[Atom],
    budget: usize,
) -> Result<Vec<(Vec<i16>, Atom)>> {
    laurent_terms(expression, variables, &vec![false; variables.len()], budget)
}

fn laurent_terms(
    expression: &Atom,
    variables: &[Atom],
    inverse_allowed: &[bool],
    budget: usize,
) -> Result<Vec<(Vec<i16>, Atom)>> {
    if variables.len() != inverse_allowed.len() {
        return Err(Error::InvalidInput(
            "compact Laurent coordinate dimensions".into(),
        ));
    }
    let terms = exact_coefficient_list(expression, variables)?;
    if terms.len() > budget {
        return Err(Error::Limit(
            "occupied polynomial term budget exhausted".into(),
        ));
    }
    terms
        .into_iter()
        .map(|(monomial, coefficient)| {
            let degrees = powers(&monomial, variables)?;
            if degrees
                .iter()
                .zip(inverse_allowed)
                .any(|(&degree, &allowed)| degree < 0 && !allowed)
                || variables.iter().any(|variable| {
                    let AtomView::Var(v) = variable.as_view() else {
                        return true;
                    };
                    !coefficient.derivative(v.get_symbol()).is_zero()
                })
            {
                return Err(Error::Unsupported(
                    "nonpolynomial occupied soft factor requires recursive weighted reduction"
                        .into(),
                ));
            }
            Ok((degrees, coefficient))
        })
        .collect()
}

fn fresh_symbol(stem: &str, sources: &[&Atom]) -> Result<Atom> {
    for suffix in 0..10_000 {
        let label = Atom::var(symbol!(format!("{stem}_{suffix}")));
        if !sources
            .iter()
            .any(|source| source.contains(label.as_view()))
        {
            return Ok(label);
        }
    }
    Err(Error::Limit(
        "occupied boundary symbol allocation exhausted".into(),
    ))
}

fn distribution_moment(
    distribution: &OccupiedBoundaryDistribution,
    spatial_dimension: &Rational,
    energy_power: i16,
    radial_power: u16,
    digits: u32,
    limits: OccupiedBoundaryLimits,
    origin_prescription: CompactOriginPrescription,
    p: Precision,
) -> Result<C> {
    let shell = &distribution.shell;
    let cut = distribution.cut_index;
    if cut <= 0 {
        return Ok(p.zero());
    }
    let continued_massless = shell.mass_squared.is_zero()
        && origin_prescription == CompactOriginPrescription::JointDimensionalMasslessTerminal;
    if continued_massless
        && (shell.chemical_potential <= 0
            || spatial_dimension <= &Rational::zero()
            || energy_power < 0)
    {
        return Err(Error::Unsupported(
            "joint massless polynomial terminal requires mu>0, d>0 and nonnegative original energy powers".into(),
        ));
    }
    let gap = &shell.chemical_potential * &shell.chemical_potential - &shell.mass_squared;
    if gap < 0 {
        return Ok(p.zero());
    }
    if distribution.lower_index > 0 {
        if continued_massless {
            // For every fixed finite set of polynomial/distribution indices,
            // sufficiently large ReD makes all required energy jets at zero
            // vanish. Continue that jointly defined product meromorphically.
            // This is not a real empty-support assertion at a massless cone.
            return Ok(p.zero());
        }
        if shell.mass_squared > 0 || distribution.upper_index > 0 && shell.chemical_potential > 0 {
            return Ok(p.zero());
        }
        return Err(Error::Unsupported(
            "massless lower-energy endpoint distribution requires dimensional continuation".into(),
        ));
    }
    let value = if distribution.upper_index == 0 {
        // CompactShell returns the occupied correction: its n=1 seed is -W.
        // C_n=1/(n-1)! (d/dm^2)^(n-1) C_1, hence raw C_n=(-1)^n I_n^occ.
        let value = if continued_massless {
            shell.dimensionally_continued_massless_raised_moment(
                spatial_dimension,
                cut as u16,
                i32::from(energy_power),
                i32::from(radial_power),
                p,
            )?
        } else {
            shell.raised_laurent_moment(
                spatial_dimension,
                cut as u16,
                i32::from(energy_power),
                radial_power,
                digits,
                limits.max_compact_series_terms,
                p,
            )?
        };
        if cut % 2 == 0 { value } else { p.neg(&value) }
    } else {
        let a = spatial_dimension.clone() / Rational::from(2) + Rational::from(radial_power);
        let radial_order = i64::from(cut - 1);
        let upper_order = i64::from(distribution.upper_index - 1);
        if gap.is_zero() {
            let degree = &a - &Rational::from(cut) - Rational::from(upper_order);
            if shell.mass_squared > 0 && degree > 0 {
                return Ok(p.zero());
            }
            return Err(Error::Unsupported(
                "coincident shell and Fermi-surface threshold requires its thermal distributional limit".into(),
            ));
        }
        let e_symbol = symbol!("rustflow_occupied_boundary::surface_energy");
        let e = Atom::var(e_symbol);
        let mut coefficient = Rational::from(1);
        for j in 0..radial_order {
            coefficient *= -(a.clone() - Rational::from(j + 1)) / Rational::from(j + 1);
        }
        for j in 0..upper_order {
            coefficient /= Rational::from(-(j + 1));
        }
        let mut expression = Atom::num(coefficient)
            * e.clone().pow(i64::from(energy_power))
            * (e.clone().pow(2) - Atom::num(shell.mass_squared.clone()))
                .pow(Atom::num(&a - &Rational::from(cut)));
        for _ in 0..upper_order {
            expression = expression.derivative(e_symbol);
        }
        let pi = C::new(p.real(1).pi(), p.real(0));
        let angular_half = p.div(
            &p.pow(
                &pi,
                &p.rational(&(spatial_dimension.clone() / Rational::from(2))),
            ),
            &p.mul(
                &p.pow(&p.scale(&pi, 2, 1), &p.rational(spatial_dimension)),
                &p.gamma_real(
                    &p.rational(&(spatial_dimension.clone() / Rational::from(2)))
                        .re,
                )?,
            ),
        );
        p.mul(
            &angular_half,
            &p.eval(
                &expression,
                &HashMap::from_iter([(e, p.rational(&shell.chemical_potential))]),
            )?,
        )
    };
    let pi = C::new(p.real(1).pi(), p.real(0));
    let dimension = spatial_dimension + &Rational::from(1);
    let normalization = p.div(
        &p.pow(&p.scale(&pi, 2, 1), &p.rational(spatial_dimension)),
        &p.pow(&pi, &p.rational(&(dimension / Rational::from(2)))),
    );
    Ok(p.mul(&normalization, &value))
}

#[cfg(test)]
mod massless_terminal_origin_tests {
    use super::*;

    #[test]
    fn joint_massless_terminal_origin_retains_upper_contacts_and_zero_lower_jets() {
        let mut distribution = OccupiedBoundaryDistribution {
            source_loop_index: 0,
            shell: CompactShell {
                mass_squared: Rational::zero(),
                chemical_potential: Rational::one(),
            },
            cut_index: 2,
            upper_index: 0,
            lower_index: 1,
        };
        let p = Precision::decimal(60).unwrap();
        let limits = OccupiedBoundaryLimits::default();
        let value = |distribution: &OccupiedBoundaryDistribution, prescription| {
            distribution_moment(
                distribution,
                &Rational::from(3),
                0,
                0,
                45,
                limits,
                prescription,
                p,
            )
        };
        assert!(matches!(
            value(&distribution, CompactOriginPrescription::Existing),
            Err(Error::Unsupported(_))
        ));
        assert_eq!(
            value(
                &distribution,
                CompactOriginPrescription::JointDimensionalMasslessTerminal
            )
            .unwrap(),
            p.zero()
        );
        distribution.lower_index = 0;
        distribution.upper_index = 1;
        // Radially integrating raw C2 at d=3 gives -Ad/(4E), hence the
        // upper contact at E=mu=1 is -1/pi in native occupied normalization.
        let expected = p.eval(&parse!("-1/pi"), &HashMap::default()).unwrap();
        assert!(
            p.close(
                &value(
                    &distribution,
                    CompactOriginPrescription::JointDimensionalMasslessTerminal
                )
                .unwrap(),
                &expected,
                45
            )
        );
        distribution.lower_index = 1;
        assert_eq!(
            value(
                &distribution,
                CompactOriginPrescription::JointDimensionalMasslessTerminal
            )
            .unwrap(),
            p.zero()
        );
        distribution.shell.chemical_potential = Rational::zero();
        assert!(matches!(
            value(
                &distribution,
                CompactOriginPrescription::JointDimensionalMasslessTerminal
            ),
            Err(Error::Unsupported(_))
        ));
    }
}
