//! Native HEPKit graphs, kinematics and numerator algebra for RustFlow.
//!
//! DOT interpretation and momentum routing belong to HEPKit/Linnet. Spenso and
//! Idenso contract the numerator before HEPKit rewrites it in a complete
//! denominator basis. This module only translates that basis into RustRed's
//! scalar-product ordering and records exact linear combinations of integrals.

use crate::cuts::{
    CutDefinition, CutFamily, CutLine, CutMetadata, LoopPrescription, MomentumRouting,
};
use crate::reduction::LinearCombination;
use crate::{Error, Integral, IntegralFamily, KinematicPoint, Propagator, Result, RunContext};
use feynkit_graph::{DiagramEndpoint, EdgeId, FeynmanDiagram, IntegralFamily as NativeFamily};
use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use idenso::tensor::{AlgebraContraction, AlgebraSettings, SymbolicTensor};
use spenso::{network::parsing::AtomStructureExt, structure::TensorStructure};
use std::{collections::BTreeMap, sync::Arc};
use symbolica::domains::float::Complex;
use symbolica::prelude::*;

fn input_error(error: impl std::fmt::Display) -> Error {
    Error::InvalidInput(format!("HEPKit: {error}"))
}

/// A native physical graph and the integral powers absent from DOT metadata.
/// Native model floating-point default values are never used as exact inputs.
pub struct GraphIntegral {
    diagram: Arc<FeynmanDiagram>,
    kinematics: Kinematics,
    edges: Vec<EdgeId>,
    powers: Vec<i16>,
    algebra: AlgebraSettings,
}

impl GraphIntegral {
    pub fn new(diagram: Arc<FeynmanDiagram>, kinematics: &Kinematics) -> Result<Self> {
        diagram.validate().map_err(input_error)?;
        if !matches!(
            kinematics.dimension().to_symbolic().as_view(),
            AtomView::Var(_)
        ) {
            return Err(Error::InvalidInput("HEPKit input needs a symbolic tensor dimension, such as D; it is replaced by D0-2*epsilon after contraction".into()));
        }
        let edges = diagram
            .edges()
            .filter(|(_, ends, edge)| {
                ends.source.is_some()
                    && ends.target.is_some()
                    && edge.external.is_none()
                    && !edge.is_dummy
            })
            .map(|(id, _, _)| id)
            .collect::<Vec<_>>();
        let powers = vec![1; edges.len()];
        Ok(Self {
            diagram,
            kinematics: kinematics.clone(),
            edges,
            powers,
            algebra: AlgebraSettings {
                contract: AlgebraContraction::Minimal,
                ..AlgebraSettings::hep()
            },
        })
    }

    pub fn from_dot(model: Arc<Model>, dot: &str, kinematics: &Kinematics) -> Result<Self> {
        Self::new(
            Arc::new(FeynmanDiagram::from_dot(model, dot).map_err(input_error)?),
            kinematics,
        )
    }

    pub fn diagram(&self) -> &Arc<FeynmanDiagram> {
        &self.diagram
    }
    pub fn kinematics(&self) -> &Kinematics {
        &self.kinematics
    }
    pub fn propagator_edges(&self) -> &[EdgeId] {
        &self.edges
    }
    pub fn powers(&self) -> &[i16] {
        &self.powers
    }

    pub fn with_powers(mut self, powers: &BTreeMap<EdgeId, i16>) -> Result<Self> {
        for (edge, power) in powers {
            let slot = self.edges.iter().position(|id| id == edge).ok_or_else(|| {
                Error::InvalidInput(format!("edge {} is not an internal propagator", edge.0))
            })?;
            self.powers[slot] = *power;
        }
        Ok(self)
    }

    pub fn with_algebra_settings(mut self, settings: AlgebraSettings) -> Self {
        self.algebra = settings;
        self
    }

    fn admit_denominators(&self, point: &KinematicPoint, keep_cut_measure: bool) -> Result<()> {
        if !keep_cut_measure && !self.diagram.cuts().is_empty() {
            return Err(Error::Unsupported("native graph cuts require the cut-integral workflow; ordinary graph conversion does not remove cut metadata".into()));
        }
        let model = self.diagram.model();
        let momentum = function!(symbol!("UFO::P"), function!(symbol!("UFO::idx"), 1, 1));
        for (edge_id, _, edge) in self
            .diagram
            .edges()
            .filter(|(id, _, _)| self.edges.contains(id))
        {
            let particle = model.particle_by_id(edge.particle).map_err(input_error)?;
            let width = model.particle_width(edge.particle).map_err(input_error)?;
            let width_symbol = Atom::var(symbol!(&format!("UFO::{}", width.name)));
            let width_value = point
                .0
                .get(&width_symbol)
                .cloned()
                .or_else(|| width.expression.as_ref().map(|a| point.apply(a)));
            if width.name != "ZERO" && !width_value.is_some_and(|a| a.together().cancel().is_zero())
            {
                return Err(Error::Unsupported(format!(
                    "edge {} has a width not established to be exactly zero; the native quadratic builder does not include widths",
                    edge_id.0
                )));
            }
            if let Some((id, _)) = model.particle_propagator(particle).map_err(input_error)? {
                let propagator = model.propagator_by_id(id).map_err(input_error)?;
                let quadratic = crate::family::substitute(
                    &(momentum.clone().pow(2) - particle.symbolic_mass(model).pow(2)),
                    &BTreeMap::from([(Atom::var(symbol!("UFO::ZERO")), Atom::new())]),
                );
                if !(&propagator.denominator - quadratic)
                    .together()
                    .cancel()
                    .is_zero()
                {
                    return Err(Error::Unsupported(format!(
                        "edge {} has a custom model denominator not represented by the native quadratic builder",
                        edge_id.0
                    )));
                }
            }
        }
        Ok(())
    }

    /// Produce one exact weighted group per native partial-fraction family.
    /// Pass these groups to `solve_integral_combinations`, which multiplies the
    /// weights at nonzero epsilon before any Laurent fitting or truncation.
    pub fn integral_groups(
        &self,
        point: &KinematicPoint,
        epsilon: Symbol,
        dimension: i64,
        max_partial_fraction_states: usize,
        context: &RunContext,
    ) -> Result<Vec<(IntegralFamily, LinearCombination)>> {
        self.integral_groups_impl(
            point,
            epsilon,
            dimension,
            max_partial_fraction_states,
            context,
            false,
        )
    }

    /// Convert one selected native physical cut without changing denominator
    /// identities. Positive-energy momenta point from the native left amplitude
    /// to the right amplitude. The same native model and routed tensor numerator
    /// are retained; no model is reloaded or reinterpreted.
    ///
    /// All original denominator slots are kept, including inactive uncut slots.
    /// Dependent denominators currently require a cut-aware partial-fraction
    /// transformation and are rejected instead of losing their measure metadata.
    /// Loop prescriptions are explicit, in the native loop-basis order. Exactly
    /// one cut is selected by index; other cuts in the native inventory are not
    /// combined with it. `CutFamily::new` checks every oriented routing squared
    /// against the retained normalized denominator before returning the result.
    #[allow(clippy::too_many_arguments)] // Native cut selection and measure accompany ordinary conversion inputs.
    pub fn cut_integral_group(
        &self,
        cut_index: usize,
        point: &KinematicPoint,
        epsilon: Symbol,
        dimension: i64,
        loop_prescriptions: Vec<LoopPrescription>,
        context: &RunContext,
    ) -> Result<(CutFamily, LinearCombination)> {
        let cuts = self.cut_metadata(cut_index, loop_prescriptions)?;
        let mut groups = self.integral_groups_impl(point, epsilon, dimension, 0, context, true)?;
        if groups.len() != 1 {
            return Err(Error::Unsupported(
                "native cut decomposition must preserve one complete denominator family".into(),
            ));
        }
        let (family, terms) = groups.remove(0);
        let family = CutFamily::new(family, cuts)?;
        let mut retained = LinearCombination::new();
        for (integral, coefficient) in terms {
            if !family.is_cut_zero(&integral)? {
                retained.insert(integral, coefficient);
            }
        }
        Ok((family, retained))
    }

    fn cut_metadata(
        &self,
        cut_index: usize,
        loop_prescriptions: Vec<LoopPrescription>,
    ) -> Result<CutMetadata> {
        let cut =
            self.diagram.cuts().get(cut_index).ok_or_else(|| {
                Error::InvalidInput("native physical cut index out of range".into())
            })?;
        let basis = self.diagram.loop_momentum_basis();
        let mut lines = Vec::with_capacity(cut.cut.len());
        for half in &cut.cut {
            let propagator = self
                .edges
                .iter()
                .position(|edge| *edge == half.edge)
                .ok_or_else(|| {
                    Error::Unsupported(
                        "selected native cut is not an integrated physical propagator".into(),
                    )
                })?;
            let signature = basis.edge_signatures.get(&half.edge).ok_or_else(|| {
                Error::InvalidInput("native cut edge has no momentum signature".into())
            })?;
            let (loops, external) = signature.integer_coefficients();
            if loops.len() != basis.loop_edges.len() || external.len() != basis.external_edges.len()
            {
                return Err(Error::InvalidInput(
                    "native cut momentum signature dimensions".into(),
                ));
            }
            let sign = match half.endpoint {
                DiagramEndpoint::Source => 1,
                DiagramEndpoint::Target => -1,
            };
            let loops = loops
                .into_iter()
                .map(|n| Rational::from(sign * n as i64))
                .collect();
            let mut retained_external = Vec::new();
            for (edge, coefficient) in basis.external_edges.iter().zip(external) {
                if basis.dependent_externals.contains(edge) {
                    if coefficient != 0 {
                        return Err(Error::Unsupported(
                            "native cut routing retains a dependent external coordinate".into(),
                        ));
                    }
                } else {
                    retained_external.push(Rational::from(sign * coefficient as i64));
                }
            }
            lines.push(CutLine {
                propagator,
                definition: CutDefinition::PositiveEnergy {
                    momentum: MomentumRouting {
                        loops,
                        external: retained_external,
                    },
                },
            });
        }
        CutMetadata::new(lines, loop_prescriptions)
    }

    /// Prepare a cut numerator through native partial fractions. Original
    /// physical slots retain their EdgeId mapping across every independent child.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_cut_combination<'a>(
        &self,
        cut_index: usize,
        point: &KinematicPoint,
        epsilon: Symbol,
        channel: &crate::cuts::FutureTimelikeChannel,
        loop_prescriptions: Vec<LoopPrescription>,
        backend: &'a dyn crate::ReductionBackend,
        options: &crate::FlowOptions,
        max_partial_fraction_states: usize,
        context: &RunContext,
    ) -> Result<crate::PreparedCutCombination<'a>> {
        let cuts = self.cut_metadata(cut_index, loop_prescriptions)?;
        let (native, numerator) = self.native_integrand(point, epsilon, true, context)?;
        crate::PreparedCutCombination::from_hepkit(
            &native,
            &self.powers,
            &numerator,
            cuts,
            channel,
            epsilon,
            backend,
            options,
            max_partial_fraction_states,
            context,
        )
    }

    /// Prepare the selected cut and its full native numerator as one projection.
    /// Exact point substitutions are applied once by the native graph converter.
    /// The shared cut dispatcher then retains terminal normalization and IBP
    /// conditions and multiplies numerator weights before epsilon fitting.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_cut_projection<'a>(
        &self,
        cut_index: usize,
        point: &KinematicPoint,
        epsilon: Symbol,
        channel: &crate::cuts::FutureTimelikeChannel,
        loop_prescriptions: Vec<LoopPrescription>,
        backend: &'a dyn crate::ReductionBackend,
        options: &crate::FlowOptions,
        context: &RunContext,
    ) -> Result<crate::PreparedCutProjections<'a>> {
        let (family, weights) = self.cut_integral_group(
            cut_index,
            point,
            epsilon,
            options.dimension,
            loop_prescriptions,
            context,
        )?;
        crate::PreparedCutProjections::new(
            &family,
            channel,
            &[weights],
            &KinematicPoint::default(),
            backend,
            options,
            context,
        )
    }

    #[allow(clippy::too_many_arguments)] // Shared ordinary/cut conversion avoids a second tensor or routing implementation.
    fn integral_groups_impl(
        &self,
        point: &KinematicPoint,
        epsilon: Symbol,
        dimension: i64,
        max_partial_fraction_states: usize,
        context: &RunContext,
        keep_cut_measure: bool,
    ) -> Result<Vec<(IntegralFamily, LinearCombination)>> {
        let (native, numerator) =
            self.native_integrand(point, epsilon, keep_cut_measure, context)?;
        // Preserve the existing ordinary and singular-family conversion
        // semantics. The grouped cut path records raw domains separately.
        let numerator = numerator.together().cancel();
        let powers = self
            .powers
            .iter()
            .map(|&n| i32::from(n))
            .collect::<Vec<_>>();
        let parts = if keep_cut_measure {
            if !native.is_independent() {
                return Err(Error::Unsupported("dependent native cut denominators require a cut-aware partial-fraction transformation".into()));
            }
            vec![(Atom::one(), powers)]
        } else {
            native
                .partial_fraction(&powers, max_partial_fraction_states)
                .map_err(|error| match error {
                    feynkit_graph::IntegralFamilyError::PartialFractionLimit(_)
                    | feynkit_graph::IntegralFamilyError::PowerOverflow => {
                        Error::Limit(error.to_string())
                    }
                    _ => input_error(error),
                })?
        };
        let mut groups = Vec::new();
        for (index, (factor, powers)) in parts.into_iter().enumerate() {
            context.cancellation.check()?;
            let mut numerator = &numerator * factor;
            for (denominator, &power) in native.denominators().iter().zip(&powers) {
                if !keep_cut_measure && power < 0 {
                    numerator *= denominator.clone().pow(-i64::from(power));
                }
            }
            let physical_powers = powers
                .iter()
                .filter(|&&n| keep_cut_measure || n > 0)
                .copied()
                .collect::<Vec<_>>();
            let family = if keep_cut_measure {
                native.complete(&[]).map_err(input_error)?
            } else {
                native
                    .sector(&powers)
                    .map_err(input_error)?
                    .complete(&[])
                    .map_err(input_error)?
            };
            let (mut converted, terms) =
                native_projection(&family, &numerator, &physical_powers, epsilon, dimension)?;
            converted.name = format!("{}::part{index}", self.diagram.name());
            groups.push((converted, terms));
        }
        Ok(groups)
    }
    fn native_integrand(
        &self,
        point: &KinematicPoint,
        epsilon: Symbol,
        keep_cut_measure: bool,
        context: &RunContext,
    ) -> Result<(NativeFamily, Atom)> {
        context.cancellation.check()?;
        for (key, value) in &point.0 {
            if !matches!(key.as_view(), AtomView::Var(_)) {
                return Err(Error::InvalidInput("HEPKit parameter substitutions require scalar symbol keys; momentum assumptions belong in native Kinematics".into()));
            }
            crate::family::scalar_symbols(value.as_view(), &mut Default::default())?;
        }
        self.admit_denominators(point, keep_cut_measure)?;
        if point
            .0
            .get(&Atom::var(symbol!("UFO::ZERO")))
            .is_some_and(|a| !a.together().cancel().is_zero())
        {
            return Err(Error::InvalidInput(
                "the native model ZERO parameter cannot be overridden by a nonzero value".into(),
            ));
        }
        let tensor_dimension = self.kinematics.dimension().to_symbolic();
        if tensor_dimension == Atom::var(epsilon) || point.0.contains_key(&tensor_dimension) {
            return Err(Error::InvalidInput("the native tensor dimension must differ from epsilon and must not be substituted before contraction".into()));
        }
        // First ask the native graph to construct its momenta and propagators
        // in a neutral chart. Then specialize scalar kinematics exactly. The
        // formal-i bridge reuses the existing rational-field adapter so native
        // affine algebra can handle exact complex points as well as real ones.
        let neutral = Kinematics::in_dimension(&tensor_dimension).map_err(input_error)?;
        let raw = self
            .diagram
            .propagator_family(&neutral)
            .map_err(input_error)?;
        if raw.denominators().len() != self.powers.len() {
            return Err(Error::InvalidInput(
                "native edge/denominator ordering mismatch".into(),
            ));
        }
        for product in raw.scalar_products() {
            if self.kinematics.apply(product) != *product {
                return Err(Error::InvalidInput("native kinematic assumptions must not constrain integrated loop scalar products".into()));
            }
        }
        let bind = |a: &Atom| {
            crate::family::encode_complex(&point.apply(&self.kinematics.apply(a)))
                .together()
                .cancel()
        };
        let mut kinematics = neutral
            .with_momenta(
                raw.loop_momenta()
                    .iter()
                    .chain(raw.external_momenta())
                    .cloned(),
            )
            .map_err(input_error)?;
        for (i, left) in raw.external_momenta().iter().enumerate() {
            for right in &raw.external_momenta()[i..] {
                let value = bind(
                    &self
                        .kinematics
                        .scalar_product(left, right)
                        .map_err(input_error)?,
                );
                kinematics = kinematics
                    .with_scalar_product(left, right, value)
                    .map_err(input_error)?;
            }
        }
        let native = NativeFamily::new(
            raw.loop_momenta().to_vec(),
            raw.external_momenta().to_vec(),
            raw.denominators().iter().map(&bind).collect(),
            &kinematics,
        )
        .map_err(input_error)?;
        let projected = self.diagram.numerator()
            * self.diagram.projector()
            * self.diagram.numerator_prefactor()
            * self.diagram.overall_factor();
        // Dimension labels are input data. Native parsing can already turn a
        // literal four-dimensional metric trace into the scalar 4, so rewriting
        // only the surviving slots would give inconsistent dimensional algebra.
        if projected.with_lorentz_dimension(tensor_dimension.as_view()) != projected {
            return Err(Error::Unsupported(
                "numerator Lorentz indices must use the symbolic Kinematics dimension; fixed-dimensional tensor slots cannot be promoted after native parsing".into(),
            ));
        }
        let routed = self
            .diagram
            .loop_momentum_basis()
            .route_expression(&projected);
        let tensor = SymbolicTensor::infer(routed)
            .and_then(|t| t.simplify_algebra(&self.algebra))
            .and_then(|t| t.to_dots())
            .map_err(input_error)?;
        if !tensor.structure().canonical().is_scalar() {
            return Err(Error::Unsupported(
                "graph numerator retains free tensor indices; supply a scalar projector".into(),
            ));
        }
        // Keep represented numerator denominator domains until the cut adapter
        // has recorded them. Propagator specialization above remains canonical.
        let numerator = crate::family::encode_complex(
            &point.apply(&self.kinematics.apply(tensor.expression())),
        );
        Ok((native, numerator))
    }
}

fn native_scalar_products(family: &NativeFamily) -> Result<Vec<Atom>> {
    let mut products = Vec::new();
    for (i, left) in family.loop_momenta().iter().enumerate() {
        for right in &family.loop_momenta()[i..] {
            products.push(
                family
                    .kinematics()
                    .scalar_product(left, right)
                    .map_err(input_error)?,
            );
        }
    }
    for left in family.loop_momenta() {
        for right in family.external_momenta() {
            products.push(
                family
                    .kinematics()
                    .scalar_product(left, right)
                    .map_err(input_error)?,
            );
        }
    }
    Ok(products)
}

fn native_projection(
    family: &NativeFamily,
    numerator: &Atom,
    physical_powers: &[i32],
    epsilon: Symbol,
    dimension: i64,
) -> Result<(IntegralFamily, LinearCombination)> {
    let tensor_dimension = family.kinematics().dimension().to_symbolic();
    let converted = IntegralFamily::from_hepkit(family, physical_powers.len(), epsilon, dimension)?;
    let labels = (0..family.denominators().len())
        .map(|i| symbol!("symbolica_amflow::native_denominator").call(i))
        .collect::<Vec<_>>();
    let rewritten = family
        .rewrite_numerator(numerator, &labels)
        .map_err(input_error)?;
    let mut terms = LinearCombination::new();
    for (monomial, coefficient) in crate::coefficient::exact_coefficient_list(&rewritten, &labels)?
    {
        let numerator_powers = crate::integrand::powers(&monomial, &labels)?;
        if numerator_powers.iter().any(|&n| n < 0) {
            return Err(Error::Unsupported(
                "graph numerator must be polynomial in loop scalar products".into(),
            ));
        }
        let mut integral = Integral(vec![0; labels.len()]);
        for (slot, &power) in physical_powers.iter().enumerate() {
            integral.0[slot] = i16::try_from(power)
                .map_err(|_| Error::Limit("native propagator power exceeds i16".into()))?;
        }
        for (power, numerator_power) in integral.0.iter_mut().zip(numerator_powers) {
            *power = power
                .checked_sub(numerator_power)
                .ok_or_else(|| Error::Limit("native numerator power exceeds i16".into()))?;
        }
        let coefficient = decode(&coefficient, &tensor_dimension, epsilon, dimension);
        *terms.entry(integral).or_default() += coefficient;
    }
    for coefficient in terms.values_mut() {
        *coefficient = coefficient.together().cancel();
    }
    terms.retain(|_, coefficient| !coefficient.is_zero());
    Ok((converted, terms))
}

pub(crate) struct NativeCutGroup {
    pub family: CutFamily,
    pub weights: LinearCombination,
    pub original_slots: Vec<usize>,
}

pub(crate) struct NativeCutDecomposition {
    pub groups: Vec<NativeCutGroup>,
    pub coefficient_conditions: Vec<Atom>,
    pub coefficient_valuation: i32,
}

/// Preserve the represented coefficient domain before rewriting or dropping a
/// cut. HEPKit scalar products are native calls; temporary scalar names only
/// adapt them to the existing structural rational-domain validator.
fn cut_numerator_domain(
    native: &NativeFamily,
    numerator: &Atom,
    epsilon: Symbol,
    dimension: i64,
) -> Result<(Vec<Atom>, i32)> {
    use std::collections::BTreeSet;
    let coordinates = native
        .scalar_products()
        .iter()
        .enumerate()
        .map(|(i, product)| {
            (
                product.clone(),
                Atom::var(symbol!(format!(
                    "symbolica_amflow::cut_inventory_coordinate_{i}"
                ))),
            )
        })
        .collect::<BTreeMap<_, _>>();
    // Scalar-product identity includes its native tensor dimension. Recognize
    // those calls before replacing D, then decode only their scalar coefficients.
    let translated = decode_unsimplified(
        &crate::family::substitute(numerator, &coordinates),
        &native.kinematics().dimension().to_symbolic(),
        epsilon,
        dimension,
    );
    let mut allowed = coordinates.values().cloned().collect::<BTreeSet<_>>();
    allowed.extend([
        Atom::var(epsilon),
        Atom::var(crate::family::imaginary_parameter()),
    ]);
    let mut present = BTreeSet::new();
    crate::family::scalar_symbols(translated.as_view(), &mut present)?;
    if !present.is_subset(&allowed) {
        return Err(Error::InvalidInput(
            "cut numerator contains unsubstituted scalar coefficients".into(),
        ));
    }
    let variables = allowed
        .iter()
        .map(|atom| match atom.as_view() {
            AtomView::Var(v) => v.get_symbol(),
            _ => unreachable!(),
        })
        .collect();
    let guards =
        crate::physical_conditions::rational_denominator_conditions(&[translated], &variables)?;
    let allowed = BTreeSet::from([
        Atom::var(epsilon),
        Atom::var(crate::family::imaginary_parameter()),
    ]);
    for guard in &guards {
        let mut symbols = BTreeSet::new();
        crate::family::scalar_symbols(guard.as_view(), &mut symbols)?;
        if !symbols.is_subset(&allowed) {
            return Err(Error::Unsupported(
                "cut numerator must be polynomial in integrated scalar products".into(),
            ));
        }
    }
    let mut valuation = 0;
    for (monomial, coefficient) in
        crate::coefficient::exact_coefficient_list(numerator, native.scalar_products())?
    {
        if crate::integrand::powers(&monomial, native.scalar_products())?
            .iter()
            .any(|&n| n < 0)
        {
            return Err(Error::Unsupported(
                "cut numerator must be polynomial in integrated scalar products".into(),
            ));
        }
        valuation = valuation.min(crate::engine::projection_weight_valuation(
            &decode(
                &coefficient,
                &native.kinematics().dimension().to_symbolic(),
                epsilon,
                dimension,
            ),
            epsilon,
            &allowed,
        )?);
    }
    Ok((guards, valuation))
}

/// Exact decomposition only: HEPKit owns affine dependence, power transfer,
/// completion and numerator rewriting. The original measure is admitted first.
#[allow(clippy::too_many_arguments)]
pub(crate) fn decompose_cut_inventory(
    native: &NativeFamily,
    powers: &[i16],
    numerator: &Atom,
    cuts: &CutMetadata,
    channel: &crate::cuts::FutureTimelikeChannel,
    epsilon: Symbol,
    dimension: i64,
    max_states: usize,
    context: &RunContext,
) -> Result<NativeCutDecomposition> {
    context.cancellation.check()?;
    if powers.len() != native.denominators().len() {
        return Err(Error::InvalidInput(
            "native cut inventory power dimensions".into(),
        ));
    }
    let original = IntegralFamily::from_hepkit_inventory(native, powers.len(), epsilon, dimension)?;
    let inventory = crate::cuts::CutInventory::new(&original, cuts)?;
    crate::cut_flow::admit_inventory(inventory, channel, context)?;
    // Admission uses exact canonical coefficient values. Reconstruct the native
    // owner from those same values rather than trusting an earlier generic rank
    // or allowing uncancelled symbolic Add coefficients into its affine parser.
    let products = native_scalar_products(native)?;
    let canonical_denominators = original
        .propagators
        .iter()
        .map(|d| {
            products
                .iter()
                .zip(&d.scalar_products)
                .fold(d.constant.clone(), |sum, (x, coefficient)| {
                    sum + coefficient * x
                })
        })
        .collect();
    let mut kinematics = native.kinematics().clone();
    for (i, left) in native.external_momenta().iter().enumerate() {
        for (j, right) in native.external_momenta().iter().enumerate().take(i + 1) {
            kinematics = kinematics
                .with_scalar_product(left, right, original.external_gram[i][j].clone())
                .map_err(input_error)?;
        }
    }
    let canonical = NativeFamily::new(
        native.loop_momenta().to_vec(),
        native.external_momenta().to_vec(),
        canonical_denominators,
        &kinematics,
    )
    .map_err(input_error)?;
    let native = &canonical;
    let (coefficient_conditions, coefficient_valuation) =
        cut_numerator_domain(native, numerator, epsilon, dimension)?;
    let powers = powers.iter().map(|&p| i32::from(p)).collect::<Vec<_>>();
    let parts = native
        .partial_fraction(&powers, max_states)
        .map_err(|error| match error {
            feynkit_graph::IntegralFamilyError::PartialFractionLimit(_)
            | feynkit_graph::IntegralFamilyError::PowerOverflow => Error::Limit(error.to_string()),
            _ => input_error(error),
        })?;
    context.cancellation.check()?;
    let mut groups = BTreeMap::<Vec<usize>, NativeCutGroup>::new();
    for (coefficient, powers) in parts {
        context.cancellation.check()?;
        // For the retained normalized cut distribution D*C_n=C_(n-1), C_0=0.
        // The original inventory was checked even when this term now vanishes.
        if cuts.lines().keys().any(|&slot| powers[slot] <= 0) {
            continue;
        }
        let original_slots = powers
            .iter()
            .enumerate()
            .filter_map(|(slot, &power)| (power > 0).then_some(slot))
            .collect::<Vec<_>>();
        let physical_powers = original_slots
            .iter()
            .map(|&slot| powers[slot])
            .collect::<Vec<_>>();
        let mut numerator = numerator * coefficient;
        for (denominator, &power) in native.denominators().iter().zip(&powers) {
            if power < 0 {
                numerator *= denominator.clone().pow(-i64::from(power));
            }
        }
        let independent = native
            .sector(&powers)
            .and_then(|f| f.complete(&[]))
            .map_err(input_error)?;
        let (ordinary, weights) = native_projection(
            &independent,
            &numerator,
            &physical_powers,
            epsilon,
            dimension,
        )?;
        let metadata = CutMetadata::new(
            cuts.lines().iter().map(|(&slot, definition)| CutLine {
                propagator: original_slots.iter().position(|&old| old == slot).unwrap(),
                definition: definition.clone(),
            }),
            cuts.loop_prescriptions().to_vec(),
        )?;
        let family = CutFamily::new(ordinary, metadata)?;
        let weights = weights
            .into_iter()
            .filter_map(|(i, c)| match family.is_cut_zero(&i) {
                Ok(true) => None,
                Ok(false) => Some(Ok((i, c))),
                Err(e) => Some(Err(e)),
            })
            .collect::<Result<LinearCombination>>()?;
        if let Some(group) = groups.get_mut(&original_slots) {
            for (integral, coefficient) in weights {
                *group.weights.entry(integral).or_default() += coefficient;
            }
        } else {
            groups.insert(
                original_slots.clone(),
                NativeCutGroup {
                    family,
                    weights,
                    original_slots,
                },
            );
        }
    }
    for group in groups.values_mut() {
        context.cancellation.check()?;
        for coefficient in group.weights.values_mut() {
            *coefficient = coefficient.together().cancel();
        }
        group
            .weights
            .retain(|_, coefficient| !coefficient.is_zero());
    }
    context.cancellation.check()?;
    Ok(NativeCutDecomposition {
        groups: groups.into_values().collect(),
        coefficient_conditions,
        coefficient_valuation,
    })
}

fn decode(a: &Atom, native_dimension: &Atom, epsilon: Symbol, dimension: i64) -> Atom {
    decode_unsimplified(a, native_dimension, epsilon, dimension)
        .together()
        .cancel()
}

fn decode_unsimplified(a: &Atom, native_dimension: &Atom, epsilon: Symbol, dimension: i64) -> Atom {
    crate::family::substitute(
        a,
        &BTreeMap::from([
            (
                native_dimension.clone(),
                Atom::num(dimension) - Atom::num(2) * Atom::var(epsilon),
            ),
            (
                Atom::var(crate::family::imaginary_parameter()),
                Atom::num(Complex::new(Rational::from(0), Rational::from(1))),
            ),
        ]),
    )
}

impl IntegralFamily {
    /// Convert a complete native HEPKit family without changing its ordered
    /// denominator basis. Scalar-product coefficients are matched by Atom
    /// identity because native HEPKit and RustRed use different slot orderings.
    pub fn from_hepkit(
        family: &NativeFamily,
        physical_propagators: usize,
        epsilon: Symbol,
        dimension: i64,
    ) -> Result<Self> {
        if !family.is_complete() || !family.is_independent() {
            return Err(Error::InvalidInput(
                "native family must be independent and complete".into(),
            ));
        }
        let converted =
            Self::from_hepkit_inventory(family, physical_propagators, epsilon, dimension)?;
        converted.validate()?;
        Ok(converted)
    }

    /// Translate an authenticated native affine inventory without claiming that
    /// its denominator rows already form an independent reduction basis.
    pub(crate) fn from_hepkit_inventory(
        family: &NativeFamily,
        physical_propagators: usize,
        epsilon: Symbol,
        dimension: i64,
    ) -> Result<Self> {
        let native_dimension = family.kinematics().dimension().to_symbolic();
        if !matches!(native_dimension.as_view(), AtomView::Var(_))
            || native_dimension == Atom::var(epsilon)
        {
            return Err(Error::InvalidInput(
                "native tensor dimension must be a symbol distinct from epsilon".into(),
            ));
        }
        let bind = |a: &Atom| decode(a, &native_dimension, epsilon, dimension);
        let products = native_scalar_products(family)?;
        if products.len() != family.scalar_products().len()
            || products
                .iter()
                .any(|a| !family.scalar_products().contains(a))
        {
            return Err(Error::InvalidInput(
                "native scalar-product basis does not match RustRed coordinates".into(),
            ));
        }
        let mut propagators = Vec::new();
        for denominator in family.denominators() {
            let mut propagator = Propagator {
                constant: Atom::new(),
                scalar_products: vec![Atom::new(); products.len()],
            };
            for (monomial, coefficient) in
                crate::coefficient::exact_coefficient_list(denominator, &products)?
            {
                if monomial.is_one() {
                    propagator.constant += bind(&coefficient);
                } else if let Some(slot) = products.iter().position(|a| a == &monomial) {
                    propagator.scalar_products[slot] += bind(&coefficient);
                } else {
                    return Err(Error::Unsupported(
                        "native denominator is not affine in loop scalar products".into(),
                    ));
                }
            }
            propagators.push(propagator);
        }
        let external_gram = family
            .external_momenta()
            .iter()
            .map(|left| {
                family
                    .external_momenta()
                    .iter()
                    .map(|right| {
                        Ok(bind(
                            &family
                                .kinematics()
                                .scalar_product(left, right)
                                .map_err(input_error)?,
                        ))
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let converted = Self {
            name: "hepkit_integral_family".into(),
            loops: (0..family.loop_momenta().len())
                .map(|i| format!("k{i}"))
                .collect(),
            external: (0..family.external_momenta().len())
                .map(|i| format!("p{i}"))
                .collect(),
            external_gram,
            propagators,
            physical_propagators,
            epsilon,
            dimension,
        };
        Ok(converted)
    }
}

#[cfg(test)]
mod cut_inventory_tests {
    use super::*;

    #[test]
    fn cancelled_native_coefficients_are_canonical_before_partial_fractions() -> Result<()> {
        let d = Atom::var(symbol!("canonical_cut_inventory_test::D"));
        let r = Atom::var(symbol!("canonical_cut_inventory_test::r"));
        let p = Atom::var(symbol!("canonical_cut_inventory_test::p"));
        let a = Atom::var(symbol!("canonical_cut_inventory_test::a"));
        let b = Atom::var(symbol!("canonical_cut_inventory_test::b"));
        let tiny = Atom::num(10).pow(-1000);
        let zero: Atom = (&a + &b).pow(2) - a.clone().pow(2) - Atom::num(2) * &a * &b - b.pow(2);
        let coefficient = &tiny + zero;
        assert!(!matches!(coefficient.as_view(), AtomView::Num(_)));
        assert!((&coefficient - &tiny).together().cancel().is_zero());
        let kin = Kinematics::in_dimension(&d)
            .map_err(input_error)?
            .with_momenta([r.clone(), p.clone()])
            .map_err(input_error)?
            .with_scalar_product(&p, &p, Atom::one())
            .map_err(input_error)?;
        let r2 = kin.scalar_product(&r, &r).map_err(input_error)?;
        let pr = &p - &r;
        let native = NativeFamily::new(
            vec![r],
            vec![p],
            vec![
                r2.clone(),
                kin.scalar_product(&pr, &pr).map_err(input_error)?,
                coefficient * r2 - 2,
            ],
            &kin,
        )
        .map_err(input_error)?;
        let cuts = CutMetadata::new(
            [(0, 1, 0), (1, -1, 1)].map(|(propagator, q, external)| CutLine {
                propagator,
                definition: CutDefinition::PositiveEnergy {
                    momentum: MomentumRouting {
                        loops: vec![q.into()],
                        external: vec![external.into()],
                    },
                },
            }),
            vec![LoopPrescription::Insensitive],
        )?;
        let result = decompose_cut_inventory(
            &native,
            &[2, 1, 1],
            &Atom::one(),
            &cuts,
            &crate::cuts::FutureTimelikeChannel {
                external: vec![1.into()],
            },
            symbol!("canonical_cut_inventory_test::eps"),
            4,
            100,
            &RunContext::default(),
        )?;
        assert_eq!(result.groups.len(), 1);
        let group = &result.groups[0];
        assert_eq!(group.original_slots, vec![0, 1]);
        assert_eq!(group.weights.len(), 2);
        assert_eq!(group.weights[&Integral(vec![2, 1])], Atom::num((-1, 2)));
        assert!(
            (&group.weights[&Integral(vec![1, 1])] + tiny / 4)
                .together()
                .cancel()
                .is_zero()
        );
        Ok(())
    }
}
