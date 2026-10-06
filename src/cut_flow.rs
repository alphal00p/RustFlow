//! Auxiliary-mass flow with an explicit positive-energy cut measure.
//!
//! Compact final-state momenta stay soft. Partial mass placements retain soft
//! virtual families, which recursively deform all remaining uncut poles. Their
//! strictly smaller pole inventory bounds cut recursion; hard factors reuse
//! ordinary recursive vacuum boundaries.
//! This adapter reuses the ordinary flow, exact region algebra, native tensor
//! projection and cut-aware IBPs; it never takes the imaginary part of a mixed
//! graph to infer its cut.

use crate::boundary::{BoundaryProvider, LeadingBoundary, RegionBoundary};
use crate::cuts::{
    CutBackendView, CutDefinition, CutFamily, CutLine, CutMetadata, FutureTimelikeChannel,
    LoopPrescription, MomentumRouting,
};
use crate::frobenius::FrobeniusBasis;
use crate::integrand::IntegralTerm;
use crate::phase_space::{PreparedMasslessPhaseSpace, PreparedTwoBodyPhaseSpace, real_rational};
use crate::*;
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Mutex;
use symbolica::prelude::*;

/// A prepared positive-energy cut system with automatic virtual vacuum
/// boundaries. Current leaves are massive two-body and massless N-body phase
/// space. Virtual loops use +i0. Explicit partial placements require independent
/// surviving soft poles; dependent decompositions and algebraic-only cuts remain
/// explicitly unsupported.
pub struct PreparedCutFlow<'a> {
    family: CutFamily,
    channel: FutureTimelikeChannel,
    flow: Option<PreparedFlow>,
    backend: &'a dyn ReductionBackend,
    hard: Vec<bool>,
    shifted: Vec<bool>,
    targets: usize,
}

impl<'a> PreparedCutFlow<'a> {
    pub fn new(
        family: &CutFamily,
        channel: &FutureTimelikeChannel,
        targets: &[Integral],
        point: &KinematicPoint,
        backend: &'a dyn ReductionBackend,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Self> {
        options.validate()?;
        context.cancellation.check()?;
        validate_options(options)?;
        if targets.is_empty() {
            return Err(Error::InvalidInput("empty cut target list".into()));
        }
        let (mut ordinary, cuts) = family.at(point)?.into_parts();
        ordinary.dimension = options.dimension;
        let family = CutFamily::new(ordinary, cuts)?;
        for entry in family.family().external_gram.iter().flatten().chain(
            family
                .family()
                .propagators
                .iter()
                .flat_map(|d| std::iter::once(&d.constant).chain(&d.scalar_products)),
        ) {
            real_rational(entry, "cut-flow kinematics and denominator coefficient")?;
        }
        for target in targets {
            family.family().validate_integral(target)?;
        }
        let hard = family
            .cuts()
            .loop_prescriptions()
            .iter()
            .map(|prescription| match prescription {
                LoopPrescription::Insensitive => Ok(false),
                LoopPrescription::PlusI0 => Ok(true),
                LoopPrescription::MinusI0 => Err(Error::Unsupported(
                    "cut flow currently requires uniform +i0 virtual directions".into(),
                )),
            })
            .collect::<Result<Vec<_>>>()?;
        let boundary = CutProvenance::new(&family, &hard)?;
        let leaf = boundary.unit_leaf()?;
        let supported = crate::phase_space::supported_leaf(&leaf, channel)?;
        let shifted = placement(&family, options)?;
        family.validate_auxiliary_mask(&shifted)?;
        if supported {
            validate_uncut_domain(family.inventory(), &boundary, channel, true)?;
        }
        let flow = if supported && shifted.iter().any(|v| *v) {
            let backend = CutBackendView {
                backend,
                cuts: family.cuts(),
            };
            Some(PreparedFlow::new(
                family.family(),
                targets,
                &KinematicPoint::default(),
                &backend,
                &flow_options(options, &shifted),
                context,
            )?)
        } else {
            None
        };
        // An all-cut input should use its terminal API directly. Avoid an
        // implicit loss of target/ISP mapping when no deformation is needed.
        if supported && flow.is_none() {
            return Err(Error::Unsupported(
                "all-cut families use the existing phase-space terminal APIs".into(),
            ));
        }
        Ok(Self {
            family,
            channel: channel.clone(),
            flow,
            backend,
            hard,
            shifted,
            targets: targets.len(),
        })
    }

    pub fn family(&self) -> &CutFamily {
        &self.family
    }

    pub fn differential_system(&self) -> Option<&DifferentialSystem> {
        self.flow.as_ref().map(|flow| &flow.system)
    }

    /// Conditions of the prepared cut-aware reduction, checked at every sample.
    pub fn nonzero_conditions(&self) -> &[Atom] {
        self.flow
            .as_ref()
            .map_or(&[], |flow| flow.reduced.nonzero_conditions.as_slice())
    }

    pub fn evaluate(
        &self,
        epsilon: &Rational,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<ComplexFloat>> {
        options.validate()?;
        validate_options(options)?;
        validate_placement(options, &self.shifted)?;
        context.cancellation.check()?;
        if epsilon.is_zero() || options.dimension != self.family.family().dimension {
            return Err(Error::InvalidInput(
                "cut evaluation requires nonzero epsilon and the prepared dimension".into(),
            ));
        }
        let options = flow_options(options, &self.shifted);
        if self.flow.is_none() {
            return Ok(vec![
                Precision::decimal(
                    options.digits + options.guard_digits
                )?
                .zero();
                self.targets
            ]);
        }
        let ordinary = recursive::RecursiveBoundary::new(self.backend, &options, context);
        self.evaluate_recursive(
            epsilon,
            &options,
            context,
            &ordinary,
            &BTreeSet::new(),
            &Mutex::default(),
        )
    }

    fn evaluate_recursive(
        &self,
        epsilon: &Rational,
        options: &FlowOptions,
        context: &RunContext,
        ordinary: &recursive::RecursiveBoundary<'_>,
        ancestry: &BTreeSet<String>,
        memo: &Mutex<BTreeMap<String, ComplexFloat>>,
    ) -> Result<Vec<ComplexFloat>> {
        context.cancellation.check()?;
        let Some(flow) = &self.flow else {
            return Ok(vec![
                Precision::decimal(
                    options.digits + options.guard_digits
                )?
                .zero();
                self.targets
            ]);
        };
        let key = self.family.fingerprint()?;
        if ancestry.contains(&key) {
            return Err(Error::Unsupported("recursive cut boundary cycle".into()));
        }
        if ancestry.len() >= 32 {
            return Err(Error::Limit(
                "cut boundary recursion exceeds 32 levels".into(),
            ));
        }
        let mut ancestry = ancestry.clone();
        ancestry.insert(key);
        let boundary = CutBoundary {
            hard: &self.hard,
            channel: &self.channel,
            backend: self.backend,
            ordinary,
            options,
            context,
            cuts: &self.family,
            ancestry: &ancestry,
            memo,
        };
        flow.evaluate(epsilon, options, &boundary, context)
    }

    pub fn evaluate_samples(
        &self,
        samples: &[Rational],
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<Vec<ComplexFloat>>> {
        options.validate()?;
        let evaluate = |(index, epsilon): (usize, &Rational)| {
            context.emit(Progress::Sample {
                index,
                total: samples.len(),
            })?;
            self.evaluate(epsilon, options, context)
        };
        if options.workers == 1 {
            samples.iter().enumerate().map(evaluate).collect()
        } else {
            rayon::ThreadPoolBuilder::new()
                .num_threads(options.workers)
                .build()
                .map_err(|error| Error::InvalidInput(error.to_string()))?
                .install(|| samples.par_iter().enumerate().map(evaluate).collect())
        }
    }

    pub fn solve(
        &self,
        last: i32,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<LaurentExpansion>> {
        let loops = i32::try_from(self.family.family().loops.len())
            .map_err(|_| Error::Limit("cut Laurent loop count overflow".into()))?;
        let leading = loops
            .checked_mul(-2)
            .ok_or_else(|| Error::Limit("cut Laurent pole bound overflow".into()))?;
        engine::fit_samples_refined_leading(
            self.targets,
            leading,
            last,
            options,
            |samples, refined| self.evaluate_samples(samples, refined, context),
        )
    }
}

fn validate_options(options: &FlowOptions) -> Result<()> {
    if options.recursion != RecursionMode::Amf || options.prescription != Prescription::PlusI0 {
        return Err(Error::Unsupported(
            "cut flow currently uses AMF and +i0 virtual boundaries".into(),
        ));
    }
    Ok(())
}

pub(crate) fn placement(family: &CutFamily, options: &FlowOptions) -> Result<Vec<bool>> {
    let ordinary = family.family();
    let mut shifted = vec![false; ordinary.propagators.len()];
    match &options.mass_mode {
        MassMode::Auto | MassMode::All => {
            for (slot, value) in shifted
                .iter_mut()
                .enumerate()
                .take(ordinary.physical_propagators)
            {
                *value = !family.cuts().is_cut(slot);
            }
        }
        MassMode::Propagators(slots) => {
            if slots.is_empty()
                || slots
                    .iter()
                    .any(|&i| i >= ordinary.physical_propagators || family.cuts().is_cut(i))
            {
                return Err(Error::InvalidInput(
                    "cut deformation needs a nonempty subset of uncut physical denominators".into(),
                ));
            }
            for &slot in slots {
                shifted[slot] = true;
            }
        }
        _ => {
            return Err(Error::Unsupported(
                "cut flow supports All, Auto, or an explicit uncut propagator subset".into(),
            ));
        }
    }
    Ok(shifted)
}

fn validate_placement(options: &FlowOptions, shifted: &[bool]) -> Result<()> {
    match &options.mass_mode {
        MassMode::Auto | MassMode::All => Ok(()),
        MassMode::Propagators(slots)
            if (0..shifted.len()).all(|i| shifted[i] == slots.contains(&i))
                && slots.iter().all(|&i| i < shifted.len()) =>
        {
            Ok(())
        }
        _ => Err(Error::InvalidInput(
            "cut evaluation changed the prepared mass placement".into(),
        )),
    }
}

fn flow_options(options: &FlowOptions, shifted: &[bool]) -> FlowOptions {
    let mut options = options.clone();
    options.mass_mode = MassMode::Propagators(
        shifted
            .iter()
            .enumerate()
            .filter_map(|(slot, &shifted)| shifted.then_some(slot))
            .collect(),
    );
    options
}

/// Unscaled soft shells and their oriented physical momenta. Scalar-product
/// normalization in the ordinary converter must never redefine this measure.
struct CutProvenance {
    template: IntegralFamily,
    denominators: Vec<Propagator>,
    momenta: Vec<MomentumRouting>,
    coordinates: Vec<Atom>,
    real: Vec<usize>,
    uncut: Vec<Propagator>,
    prescriptions: Vec<LoopPrescription>,
}

impl CutProvenance {
    fn new(family: &CutFamily, hard: &[bool]) -> Result<Self> {
        let shifted = (0..family.family().propagators.len())
            .map(|i| !family.cuts().is_cut(i))
            .collect::<Vec<_>>();
        Self::project(family, hard, &shifted)
    }

    fn project(family: &CutFamily, hard: &[bool], shifted: &[bool]) -> Result<Self> {
        Self::project_inventory(family.inventory(), hard, shifted)
    }

    fn project_inventory(
        family: crate::cuts::CutInventory<'_>,
        hard: &[bool],
        shifted: &[bool],
    ) -> Result<Self> {
        let ordinary = family.family();
        let real = family
            .cuts()
            .loop_prescriptions()
            .iter()
            .enumerate()
            .filter_map(|(i, p)| (*p == LoopPrescription::Insensitive).then_some(i))
            .collect::<Vec<_>>();
        let kept = (0..hard.len()).filter(|&i| !hard[i]).collect::<Vec<_>>();
        if real.is_empty()
            || family.cuts().lines().len() != real.len() + 1
            || real.iter().any(|&i| hard[i])
        {
            return Err(Error::Unsupported(
                "cuts must form one complete soft final state in the insensitive loop subspace"
                    .into(),
            ));
        }
        let mut labels = Vec::new();
        for i in 0..hard.len() {
            for j in i..hard.len() {
                labels.push((i, j));
            }
        }
        for i in 0..hard.len() {
            for j in 0..ordinary.external.len() {
                labels.push((i, hard.len() + j));
            }
        }
        let selected = labels
            .iter()
            .enumerate()
            .filter_map(|(i, &(a, b))| {
                (kept.contains(&a) && (b >= hard.len() || kept.contains(&b))).then_some(i)
            })
            .collect::<Vec<_>>();
        let project = |d: &Propagator| Propagator {
            constant: d.constant.clone(),
            scalar_products: selected
                .iter()
                .map(|&i| d.scalar_products[i].clone())
                .collect(),
        };
        let mut denominators = Vec::new();
        let mut momenta = Vec::new();
        for (&slot, definition) in family.cuts().lines() {
            let CutDefinition::PositiveEnergy { momentum } = definition else {
                return Err(Error::Unsupported(
                    "distribution cuts do not specify physical positive-energy boundaries".into(),
                ));
            };
            if momentum
                .loops
                .iter()
                .zip(hard)
                .any(|(a, h)| *h && !a.is_zero())
            {
                return Err(Error::Unsupported(
                    "cut momenta must have no hard-direction component".into(),
                ));
            }
            momenta.push(MomentumRouting {
                loops: kept.iter().map(|&i| momentum.loops[i].clone()).collect(),
                external: momentum.external.clone(),
            });
            denominators.push(project(&ordinary.propagators[slot]));
        }
        let mut uncut = Vec::new();
        for (slot, d) in ordinary
            .propagators
            .iter()
            .enumerate()
            .take(ordinary.physical_propagators)
        {
            if family.cuts().is_cut(slot) || shifted[slot] {
                continue;
            }
            let branch = regions::branch(d, hard.len())?;
            if branch.iter().zip(hard).all(|(a, h)| !*h || a.is_zero()) {
                uncut.push(project(d));
            }
        }
        independent_soft_poles(denominators.iter().chain(&uncut))?;
        let template = IntegralFamily {
            loops: kept.iter().map(|&i| ordinary.loops[i].clone()).collect(),
            propagators: vec![],
            physical_propagators: 0,
            ..ordinary.clone()
        };
        let coordinates = (0..selected.len())
            .map(|i| Atom::var(symbol!(format!("symbolica_amflow::cut_soft_{i}"))))
            .collect();
        let prescriptions = kept
            .iter()
            .map(|&i| family.cuts().loop_prescriptions()[i])
            .collect();
        Ok(Self {
            template,
            denominators,
            momenta,
            coordinates,
            real,
            uncut,
            prescriptions,
        })
    }

    fn unit_leaf(&self) -> Result<CutFamily> {
        let expression = self.denominators.iter().fold(Atom::one(), |value, d| {
            value
                / self
                    .coordinates
                    .iter()
                    .zip(&d.scalar_products)
                    .fold(d.constant.clone(), |sum, (x, c)| sum + x * c)
        });
        let mut terms =
            integrand::to_integrals(&expression, &self.coordinates, &self.template, 10000)?;
        if terms.len() != 1 {
            return Err(Error::Unsupported(
                "unit cut shells do not define an independent leaf".into(),
            ));
        }
        self.restore(&mut terms.remove(0))?
            .ok_or_else(|| Error::InvalidInput("unit cut volume lost a required shell".into()))
    }

    fn restore(&self, term: &mut IntegralTerm) -> Result<Option<CutFamily>> {
        if term.family.loops != self.template.loops
            || term.family.external != self.template.external
        {
            return Err(Error::Unsupported(
                "boundary factor changed the retained cut coordinate space".into(),
            ));
        }
        let expected = self
            .denominators
            .iter()
            .chain(&self.uncut)
            .collect::<Vec<_>>();
        let mut matched = BTreeSet::new();
        let mut identities = Vec::new();
        for d in &expected {
            let pivot = d
                .scalar_products
                .iter()
                .position(|c| !c.is_zero())
                .ok_or_else(|| Error::Unsupported("constant soft pole".into()))?;
            let mut matches = Vec::new();
            for (slot, actual) in term.family.propagators.iter().enumerate() {
                let scale = (&actual.scalar_products[pivot] / &d.scalar_products[pivot])
                    .together()
                    .cancel();
                if !scale.is_zero()
                    && std::iter::once((&actual.constant, &d.constant))
                        .chain(actual.scalar_products.iter().zip(&d.scalar_products))
                        .all(|(a, b)| (a - &scale * b).together().cancel().is_zero())
                {
                    if real_rational(&scale, "soft pole normalization")? <= Rational::zero() {
                        return Err(Error::Unsupported(
                            "soft pole normalization must preserve the cut orientation and +i0"
                                .into(),
                        ));
                    }
                    matches.push((slot, scale));
                }
            }
            if matches.len() > 1 {
                return Err(Error::Unsupported(
                    "ambiguous soft pole identity in a boundary factor".into(),
                ));
            }
            let identity = matches.pop();
            if let Some((slot, _)) = &identity {
                matched.insert(*slot);
            }
            identities.push(identity);
        }
        if (0..term.family.physical_propagators).any(|i| !matched.contains(&i)) {
            return Err(Error::Unsupported(
                "an unauthenticated uncut pole survived in the soft boundary".into(),
            ));
        }
        // Validate every surviving pole before declaring a missing cut zero.
        if identities[..self.denominators.len()]
            .iter()
            .any(|identity| {
                identity
                    .as_ref()
                    .is_none_or(|(slot, _)| term.integral.0[*slot] <= 0)
            })
        {
            return Ok(None);
        }
        let mut lines = Vec::new();
        for (index, (expected, identity)) in expected.into_iter().zip(identities).enumerate() {
            let Some((slot, scale)) = identity else {
                continue;
            };
            term.coefficient *= scale.pow(-i64::from(term.integral.0[slot]));
            term.family.propagators[slot] = expected.clone();
            if index < self.denominators.len() {
                lines.push(CutLine {
                    propagator: slot,
                    definition: CutDefinition::PositiveEnergy {
                        momentum: self.momenta[index].clone(),
                    },
                });
            }
        }
        Ok(Some(CutFamily::new(
            term.family.clone(),
            CutMetadata::new(lines, self.prescriptions.clone())?,
        )?))
    }
}

/// Admit the complete original measure before partial fractions can remove a
/// pole or a cut. No dependent inventory is passed to the reduction backend.
pub(crate) fn admit_inventory(
    family: crate::cuts::CutInventory<'_>,
    channel: &FutureTimelikeChannel,
    context: &RunContext,
) -> Result<()> {
    context.cancellation.check()?;
    family.validate()?;
    let ordinary = family.family();
    for value in ordinary.external_gram.iter().flatten().chain(
        ordinary
            .propagators
            .iter()
            .flat_map(|d| std::iter::once(&d.constant).chain(&d.scalar_products)),
    ) {
        real_rational(value, "cut inventory kinematics and quadratic coefficient")?;
    }
    let hard = family
        .cuts()
        .loop_prescriptions()
        .iter()
        .map(|p| match p {
            LoopPrescription::Insensitive => Ok(false),
            LoopPrescription::PlusI0 => Ok(true),
            LoopPrescription::MinusI0 => Err(Error::Unsupported(
                "cut inventory requires uniform +i0 virtual directions".into(),
            )),
        })
        .collect::<Result<Vec<_>>>()?;
    let shifted = (0..ordinary.propagators.len())
        .map(|i| !family.cuts().is_cut(i))
        .collect::<Vec<_>>();
    let boundary = CutProvenance::project_inventory(family, &hard, &shifted)?;
    let leaf = boundary.unit_leaf()?;
    let supported = crate::phase_space::supported_leaf(&leaf, channel)?;
    validate_uncut_domain(family, &boundary, channel, supported)?;
    context.cancellation.check()
}

fn validate_uncut_domain(
    family: crate::cuts::CutInventory<'_>,
    boundary: &CutProvenance,
    channel: &FutureTimelikeChannel,
    supported: bool,
) -> Result<()> {
    let ordinary = family.family();
    for slot in 0..ordinary.physical_propagators {
        if family.cuts().is_cut(slot) {
            continue;
        }
        let propagator = &ordinary.propagators[slot];
        let branch = regions::branch(propagator, ordinary.loops.len())?;
        let pivot = branch.iter().position(|c| !c.is_zero()).unwrap();
        let offset = ordinary.loops.len() * (ordinary.loops.len() + 1) / 2;
        let diagonal = pivot * ordinary.loops.len() - pivot * pivot.saturating_sub(1) / 2;
        if real_rational(
            &propagator.scalar_products[diagonal],
            "uncut quadratic normalization",
        )? <= Rational::zero()
        {
            return Err(Error::Unsupported(
                "uncut quadratic normalization must preserve +i0".into(),
            ));
        }
        for (i, b) in branch.iter().enumerate() {
            for external in 0..ordinary.external.len() {
                if (&propagator.scalar_products[offset + i * ordinary.external.len() + external]
                    - b * &propagator.scalar_products
                        [offset + pivot * ordinary.external.len() + external])
                    .together()
                    .cancel()
                    .is_zero()
                {
                    continue;
                }
                return Err(Error::Unsupported(
                    "uncut denominator is not a physical quadratic momentum square".into(),
                ));
            }
        }
        if real_rational(&ordinary.mass_squared(slot)?, "uncut squared mass")? < Rational::zero() {
            return Err(Error::Unsupported(
                "cut flow currently requires nonnegative uncut squared masses".into(),
            ));
        }
        // A common imaginary displacement of virtual diagonal scalar products
        // preserves every affine relation. Its exact rate is read from the
        // original quadratic form, before any denominator normalization.
        let mut diagonal = 0;
        let mut rate = Rational::zero();
        for (loop_index, prescription) in family.cuts().loop_prescriptions().iter().enumerate() {
            if *prescription == LoopPrescription::PlusI0 {
                rate += real_rational(
                    &propagator.scalar_products[diagonal],
                    "virtual regulator rate",
                )?;
            }
            diagonal += ordinary.loops.len() - loop_index;
        }
        match family.propagator_loop_prescription(slot)? {
            LoopPrescription::PlusI0 if rate > Rational::zero() => {}
            LoopPrescription::Insensitive if rate.is_zero() => {
                if supported {
                    certify_soft_denominator(family, boundary, channel, slot)?;
                }
            }
            _ => {
                return Err(Error::Unsupported(
                    "no compatible positive virtual regulator for cut denominator inventory".into(),
                ));
            }
        }
    }
    Ok(())
}

/// Bound a real-only denominator using qi.qj >= 0 and sum_i qi=P. These exact
/// inequalities hold on the entire compact physical phase space. An inconclusive
/// bound is unsupported; endpoint checks cannot certify absence of an interior pole.
fn certify_soft_denominator(
    family: crate::cuts::CutInventory<'_>,
    boundary: &CutProvenance,
    channel: &FutureTimelikeChannel,
    slot: usize,
) -> Result<()> {
    let ordinary = family.family();
    if ordinary.external.len() != 1 || channel.external.len() != 1 {
        return Err(Error::Unsupported(
            "real-only cut propagator support bounds currently use one external channel".into(),
        ));
    }
    let n = boundary.real.len();
    let inverse = crate::algebra::inverse(
        &boundary.momenta[..n]
            .iter()
            .map(|q| q.loops.iter().cloned().map(Atom::num).collect())
            .collect::<Vec<_>>(),
    )?;
    let channel = &channel.external[0];
    let external = vec![Atom::num(channel.clone().inv()); n + 1];
    let vectors = inverse
        .iter()
        .map(|row| {
            (0..=n)
                .map(|particle| {
                    row.iter()
                        .enumerate()
                        .fold(Atom::zero(), |sum, (j, value)| {
                            sum + value
                                * (Atom::num(i64::from(j == particle))
                                    - Atom::num(boundary.momenta[j].external[0].clone())
                                        * &external[particle])
                        })
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut pairs = vec![vec![Atom::zero(); n + 1]; n + 1];
    let mut accumulate = |left: &[Atom], right: &[Atom], coefficient: &Atom| {
        for i in 0..=n {
            for j in i..=n {
                pairs[i][j] += coefficient
                    * if i == j {
                        &left[i] * &right[j]
                    } else {
                        &left[i] * &right[j] + &left[j] * &right[i]
                    };
            }
        }
    };
    let d = &ordinary.propagators[slot];
    let mut coordinate = 0;
    for i in 0..ordinary.loops.len() {
        for j in i..ordinary.loops.len() {
            if let (Some(a), Some(b)) = (
                boundary.real.iter().position(|&k| k == i),
                boundary.real.iter().position(|&k| k == j),
            ) {
                accumulate(&vectors[a], &vectors[b], &d.scalar_products[coordinate]);
            }
            coordinate += 1;
        }
    }
    for i in 0..ordinary.loops.len() {
        if let Some(a) = boundary.real.iter().position(|&k| k == i) {
            accumulate(&vectors[a], &external, &d.scalar_products[coordinate]);
        }
        coordinate += 1;
    }
    let masses = family
        .cuts()
        .lines()
        .keys()
        .map(|&cut| real_rational(&family.cut_mass_squared(cut)?, "cut squared mass"))
        .collect::<Result<Vec<_>>>()?;
    let mut constant = real_rational(&d.constant, "uncut constant")?;
    let mut coefficients = Vec::new();
    for i in 0..=n {
        constant += real_rational(&pairs[i][i], "cut support coefficient")? * &masses[i];
        for coefficient in &pairs[i][i + 1..] {
            coefficients.push(real_rational(coefficient, "cut support coefficient")?);
        }
    }
    let invariant =
        channel * channel * real_rational(&ordinary.external_gram[0][0], "channel invariant")?;
    let total = (invariant - masses.iter().sum::<Rational>()) / Rational::from(2);
    let lower = &constant + &(coefficients.iter().min().unwrap() * &total);
    let upper = &constant + &(coefficients.iter().max().unwrap() * &total);
    if lower <= Rational::zero() && upper >= Rational::zero() {
        return Err(Error::Unsupported(
            "cannot certify an uncut denominator away from zero on physical phase space".into(),
        ));
    }
    Ok(())
}

struct CutBoundary<'a, 'b> {
    hard: &'a [bool],
    channel: &'a FutureTimelikeChannel,
    backend: &'a dyn ReductionBackend,
    ordinary: &'a recursive::RecursiveBoundary<'b>,
    options: &'a FlowOptions,
    context: &'a RunContext,
    cuts: &'a CutFamily,
    ancestry: &'a BTreeSet<String>,
    memo: &'a Mutex<BTreeMap<String, ComplexFloat>>,
}

enum RequestFamily {
    Ordinary(IntegralFamily),
    Cut(CutFamily),
}
struct Request {
    family: RequestFamily,
    targets: BTreeSet<Integral>,
}
struct Product {
    coefficient: Atom,
    factors: Vec<(String, Integral)>,
}
struct Pending {
    component: usize,
    series: usize,
    order: usize,
    jacobian: ComplexFloat,
    products: Vec<Product>,
}

impl BoundaryProvider for CutBoundary<'_, '_> {
    fn leading(
        &self,
        _: &IntegralFamily,
        _: &[Integral],
        _: &Rational,
        _: Precision,
    ) -> Result<Vec<LeadingBoundary>> {
        Err(Error::Unsupported(
            "cut boundaries require the explicit deformation and region expansion".into(),
        ))
    }

    fn constants_with_deformation(
        &self,
        family: &IntegralFamily,
        basis: &[Integral],
        shifted: &[bool],
        epsilon: &Rational,
        p: Precision,
        solutions: &FrobeniusBasis,
    ) -> Result<Vec<ComplexFloat>> {
        let regions =
            if (0..family.physical_propagators).all(|i| self.cuts.cuts().is_cut(i) || shifted[i]) {
                vec![regions::LoopRegion {
                    transformation: (0..self.hard.len())
                        .map(|i| {
                            (0..self.hard.len())
                                .map(|j| Atom::num(i64::from(i == j)))
                                .collect()
                        })
                        .collect(),
                    hard: self.hard.to_vec(),
                    hard_branches: vec![],
                    jacobian_determinant: Atom::one(),
                }]
            } else {
                crate::cut_regions::enumerate(self.cuts, self.context)?
            };
        let mut provenance = Vec::new();
        for region in &regions {
            self.context.cancellation.check()?;
            let transformed = CutFamily::new(
                family.transform_loops(&region.transformation)?,
                self.cuts.cuts().clone(),
            )?;
            provenance.push(CutProvenance::project(&transformed, &region.hard, shifted)?);
        }
        let parent_poles = distinct_uncut_poles(self.cuts)?;
        let parameters =
            ahash::HashMap::from_iter([(Atom::var(family.epsilon), p.rational(epsilon))]);
        let powers = basis
            .iter()
            .map(|integral| {
                regions
                    .iter()
                    .map(|region| {
                        regions::expand_region(family, integral, shifted, region, 0)
                            .map(|expansion| -expansion.eta_power)
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let orders = solutions.region_orders(&powers, 32)?;
        self.context.emit(Progress::BoundaryPlan {
            basis_size: basis.len(),
            region_series: orders
                .iter()
                .flatten()
                .filter(|order| order.is_some())
                .count(),
            coefficients: orders
                .iter()
                .flatten()
                .flatten()
                .map(|order| order + 1)
                .sum(),
            max_half_order: orders
                .iter()
                .flatten()
                .flatten()
                .copied()
                .max()
                .unwrap_or(0),
        })?;
        let mut data = vec![Vec::new(); basis.len()];
        let mut requests = BTreeMap::<String, Request>::new();
        let mut pending = Vec::new();
        for (component, integral) in basis.iter().enumerate() {
            for (r, region) in regions.iter().enumerate() {
                self.context.cancellation.check()?;
                let order = orders[component][r];
                for parity in 0..2 {
                    data[component].push(RegionBoundary {
                        exponent: (&powers[component][r] + Atom::num((parity as i64, 2)))
                            .together()
                            .cancel(),
                        coefficients: vec![
                            p.zero();
                            order.map_or(0, |order| (order + 2 - parity) / 2)
                        ],
                    });
                }
                let Some(order) = order else {
                    continue;
                };
                let expansion = regions::expand_region(family, integral, shifted, region, order)?;
                let determinant = p.eval(&expansion.jacobian_determinant, &parameters)?;
                let jacobian = p.pow(
                    &ComplexFloat::new(p.norm(&determinant), p.real(0)),
                    &p.rational(&(Rational::from(family.dimension) - epsilon * &Rational::from(2))),
                );
                for (index, expression) in expansion.coefficients.iter().enumerate() {
                    self.context.cancellation.check()?;
                    let terms = integrand::factor_region(
                        expression,
                        &expansion.coordinates,
                        family,
                        &region.hard,
                        10000,
                    )?;
                    let mut products = Vec::new();
                    for term in terms {
                        self.context.cancellation.check()?;
                        let mut product = Product {
                            coefficient: term.coefficient,
                            factors: vec![],
                        };
                        let mut children = Vec::new();
                        let mut zero = false;
                        for mut factor in term.factors {
                            let child = if factor.family.external.is_empty() {
                                zero |= recursive::scaleless(&factor.family, &factor.integral)?;
                                RequestFamily::Ordinary(factor.family.clone())
                            } else {
                                let Some(cut) = provenance[r].restore(&mut factor)? else {
                                    zero = true;
                                    continue;
                                };
                                if distinct_uncut_poles(&cut)? >= parent_poles {
                                    return Err(Error::Unsupported("soft cut recursion did not strictly decrease distinct uncut poles".into()));
                                }
                                zero |= cut_scaleless(&cut, &factor.integral)?;
                                RequestFamily::Cut(cut)
                            };
                            product.coefficient *= &factor.coefficient;
                            let key = match &child {
                                RequestFamily::Ordinary(family) => {
                                    format!("ordinary:{}", family.convert()?.family.fingerprint())
                                }
                                RequestFamily::Cut(family) => {
                                    format!("cut:{}", family.fingerprint()?)
                                }
                            };
                            product.factors.push((key.clone(), factor.integral.clone()));
                            children.push((key, child, factor.integral));
                        }
                        if zero {
                            continue;
                        }
                        for (key, family, integral) in children {
                            requests
                                .entry(key)
                                .or_insert_with(|| Request {
                                    family,
                                    targets: BTreeSet::new(),
                                })
                                .targets
                                .insert(integral);
                            if requests.len() > 10000 {
                                return Err(Error::Limit(
                                    "cut boundary family budget exhausted".into(),
                                ));
                            }
                        }
                        products.push(product);
                    }
                    pending.push(Pending {
                        component,
                        series: 2 * r + index % 2,
                        order: index / 2,
                        jacobian: jacobian.clone(),
                        products,
                    });
                }
            }
        }
        let mut values = BTreeMap::new();
        for (key, request) in requests {
            self.context.cancellation.check()?;
            let targets = request.targets.into_iter().collect::<Vec<_>>();
            let evaluated = match request.family {
                RequestFamily::Ordinary(family) => {
                    self.ordinary.evaluate_many(&family, &targets, epsilon, p)?
                }
                RequestFamily::Cut(family) => self.evaluate_cut(&family, &targets, epsilon, p)?,
            };
            for (target, value) in targets.into_iter().zip(evaluated) {
                values.insert((key.clone(), target), value);
            }
        }
        for coefficient in pending {
            self.context.cancellation.check()?;
            let mut value = p.zero();
            for product in coefficient.products {
                let mut term = p.eval(&product.coefficient, &parameters)?;
                for key in product.factors {
                    let value = values.get(&key).ok_or_else(|| {
                        Error::Numerical("missing collected cut boundary factor".into())
                    })?;
                    term = p.mul(&term, value);
                }
                value = p.add(&value, &term);
            }
            data[coefficient.component][coefficient.series].coefficients[coefficient.order] =
                p.mul(&value, &coefficient.jacobian);
        }
        self.context.cancellation.check()?;
        solutions.match_regions(&data)
    }
}

/// Distinct affine pole shapes, including constants and excluding ISP slots.
fn distinct_uncut_poles(family: &CutFamily) -> Result<usize> {
    let mut shapes = BTreeSet::new();
    for (slot, d) in family
        .family()
        .propagators
        .iter()
        .enumerate()
        .take(family.family().physical_propagators)
    {
        if family.cuts().is_cut(slot) {
            continue;
        }
        let scale = d
            .scalar_products
            .iter()
            .find(|c| !c.is_zero())
            .ok_or_else(|| Error::Unsupported("constant uncut pole".into()))?;
        shapes.insert(
            std::iter::once(&d.constant)
                .chain(&d.scalar_products)
                .map(|a| (a / scale).together().cancel())
                .collect::<Vec<_>>(),
        );
    }
    Ok(shapes.len())
}

/// Reuse the native zero certificate with the required cuts when a partial
/// region leaves a polynomial soft virtual direction. Exclusion is not a proof.
fn cut_scaleless(family: &CutFamily, integral: &Integral) -> Result<bool> {
    if family
        .cuts()
        .loop_prescriptions()
        .iter()
        .all(|p| *p == LoopPrescription::Insensitive)
    {
        return Ok(false);
    }
    use rustred::sector::{
        Mask,
        zero::{Analyzer, Decision},
    };
    let converted = family.family().convert()?;
    let analyzer = Analyzer::try_new(&converted.family, family.native_restrictions()?)
        .map_err(|e| Error::Reduction(e.to_string()))?;
    let mask = Mask::try_new(integral.0.iter().map(|&n| n > 0))
        .map_err(|e| Error::InvalidInput(e.to_string()))?;
    Ok(matches!(
        analyzer
            .analyze(&mask)
            .map_err(|e| Error::Reduction(e.to_string()))?,
        Decision::ProvedZero(_)
    ))
}

impl CutBoundary<'_, '_> {
    fn evaluate_cut(
        &self,
        family: &CutFamily,
        targets: &[Integral],
        epsilon: &Rational,
        p: Precision,
    ) -> Result<Vec<ComplexFloat>> {
        self.context.cancellation.check()?;
        let prefix = format!("{}:{epsilon}:{}", family.fingerprint()?, p.bits);
        let keys = targets
            .iter()
            .map(|target| format!("{prefix}:{:?}", target.0))
            .collect::<Vec<_>>();
        let lock_error = |_| Error::Numerical("cut boundary memo poisoned".into());
        let missing = {
            let memo = self.memo.lock().map_err(lock_error)?;
            targets
                .iter()
                .zip(&keys)
                .filter(|(_, key)| !memo.contains_key(*key))
                .map(|(target, _)| target.clone())
                .collect::<Vec<_>>()
        };
        if !missing.is_empty() {
            let mut options = self.options.clone();
            options.mass_mode = MassMode::All;
            // Match the actual Frobenius working precision, including adaptive
            // retries. This remains working precision, not accuracy evidence.
            let working_digits =
                (u64::from(p.bits.saturating_sub(16)) * 1000).div_ceil(3322) as u32;
            let extra = working_digits.saturating_sub(options.digits + options.guard_digits);
            options.guard_digits = options
                .guard_digits
                .max(working_digits.saturating_sub(options.digits));
            options.series_order += extra as usize * 4 / 5;
            options.refine_rational_order(extra.div_ceil(20) as usize);
            let pure_cut =
                (0..family.family().physical_propagators).all(|i| family.cuts().is_cut(i));
            let evaluated = if pure_cut && family.family().loops.len() == 1 {
                PreparedTwoBodyPhaseSpace::new(
                    family,
                    self.channel,
                    &missing,
                    &KinematicPoint::default(),
                    self.backend,
                    &options,
                    self.context,
                )?
                .evaluate(epsilon, &options, self.context)?
            } else if pure_cut {
                PreparedMasslessPhaseSpace::new(
                    family,
                    self.channel,
                    &missing,
                    &KinematicPoint::default(),
                    self.backend,
                    &options,
                    self.context,
                )?
                .evaluate(epsilon, &options, self.context)?
            } else {
                let flow = PreparedCutFlow::new(
                    family,
                    self.channel,
                    &missing,
                    &KinematicPoint::default(),
                    self.backend,
                    &options,
                    self.context,
                )?;
                flow.evaluate_recursive(
                    epsilon,
                    &flow_options(&options, &flow.shifted),
                    self.context,
                    self.ordinary,
                    self.ancestry,
                    self.memo,
                )?
            };
            self.context.cancellation.check()?;
            let mut memo = self.memo.lock().map_err(lock_error)?;
            if memo
                .len()
                .checked_add(missing.len())
                .is_none_or(|size| size > 10000)
            {
                return Err(Error::Limit(
                    "cut boundary value memo budget exhausted".into(),
                ));
            }
            for (target, value) in missing.into_iter().zip(evaluated) {
                memo.insert(format!("{prefix}:{:?}", target.0), value);
            }
        }
        let memo = self.memo.lock().map_err(lock_error)?;
        keys.into_iter()
            .map(|key| {
                memo.get(&key)
                    .cloned()
                    .ok_or_else(|| Error::Numerical("missing evaluated cut boundary".into()))
            })
            .collect()
    }
}

fn independent_soft_poles<'a>(poles: impl Iterator<Item = &'a Propagator>) -> Result<()> {
    // Exclude the affine constant column: different masses on one momentum
    // are dependent scalar-product poles despite independent affine forms.
    let rows = poles.map(|d| d.scalar_products.clone()).collect::<Vec<_>>();
    let count = rows.len();
    if crate::algebra::rref(rows).1.len() != count {
        return Err(Error::Unsupported(
            "dependent soft poles need a cut-preserving decomposition before partial fractions"
                .into(),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distinct_masses_do_not_authorize_dependent_soft_partial_fractions() {
        let poles = [
            Propagator {
                constant: Atom::zero(),
                scalar_products: vec![Atom::one(), Atom::zero()],
            },
            Propagator {
                constant: Atom::num(-2),
                scalar_products: vec![Atom::zero(), Atom::one()],
            },
            Propagator {
                constant: Atom::num(-3),
                scalar_products: vec![Atom::zero(), Atom::one()],
            },
        ];
        assert!(
            matches!(independent_soft_poles(poles.iter()),Err(Error::Unsupported(message)) if message.contains("before partial fractions"))
        );
        assert!(independent_soft_poles(poles[..2].iter()).is_ok());
    }

    fn provenance() -> CutProvenance {
        let gram = vec![vec![Atom::num(25)]];
        let denominators = vec![
            Propagator::quadratic(&[2], &[0], Atom::one(), &gram).unwrap(),
            Propagator::quadratic(&[-2], &[1], Atom::num(4), &gram).unwrap(),
        ];
        CutProvenance {
            template: IntegralFamily {
                name: "cut_distribution_identity".into(),
                loops: vec!["r".into()],
                external: vec!["p".into()],
                external_gram: gram,
                propagators: vec![],
                physical_propagators: 0,
                epsilon: symbol!("cut_provenance_test::eps"),
                dimension: 4,
            },
            denominators,
            momenta: vec![
                MomentumRouting {
                    loops: vec![2.into()],
                    external: vec![0.into()],
                },
                MomentumRouting {
                    loops: vec![(-2).into()],
                    external: vec![1.into()],
                },
            ],
            coordinates: vec![
                Atom::var(symbol!("cut_provenance_test::rr")),
                Atom::var(symbol!("cut_provenance_test::rp")),
            ],
            real: vec![0],
            uncut: vec![],
            prescriptions: vec![LoopPrescription::Insensitive],
        }
    }

    #[test]
    fn boundary_cut_normalization_and_numerator_cancellation_are_distributional() -> Result<()> {
        let provenance = provenance();
        let d1 = Atom::num(4) * &provenance.coordinates[0] - Atom::one();
        let d2 = Atom::num(4) * &provenance.coordinates[0]
            - Atom::num(4) * &provenance.coordinates[1]
            + Atom::num(21);
        // The native converter normalizes both quadratic coefficients by four.
        // Restore a raised cut without changing either its shell or derivative.
        let mut terms = integrand::to_integrals(
            &(d1.clone().pow(-2) / &d2),
            &provenance.coordinates,
            &provenance.template,
            100,
        )?;
        assert_eq!(terms.len(), 1);
        let mut term = terms.remove(0);
        let cut = provenance.restore(&mut term)?.unwrap();
        assert!(term.coefficient.is_one());
        assert_eq!(term.integral.0.iter().sum::<i16>(), 3);
        assert_eq!(cut.cut_mass_squared(0)?, Atom::one());
        // D delta(D) = 0: cancellation cannot turn a cut into an ordinary line.
        let mut terms = integrand::to_integrals(
            &(&d1 / (&d1 * &d2)),
            &provenance.coordinates,
            &provenance.template,
            100,
        )?;
        assert_eq!(terms.len(), 1);
        assert!(provenance.restore(&mut terms.remove(0))?.is_none());
        // An unknown surviving pole is an error even if another cut is absent.
        let mut terms = integrand::to_integrals(
            &((&d1 + Atom::one()).pow(-1)),
            &provenance.coordinates,
            &provenance.template,
            100,
        )?;
        assert!(
            matches!(provenance.restore(&mut terms.remove(0)), Err(Error::Unsupported(message))
            if message.contains("unauthenticated uncut"))
        );
        Ok(())
    }
}
