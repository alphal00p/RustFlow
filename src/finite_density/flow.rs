//! Native occupied AMF: source discovery, integrated regions and shared transport.
//! The first numerical admission is a strictly massive, nonpinched sunset
//! domain. This is an explicit contour restriction, independent of input names.
use super::PreparedDensityInput;
use super::boundary::OccupiedBoundaryLimits;
use super::flow_boundary::OccupiedFlowBoundary;
use super::geometry::OccupiedCutFamily;
use super::guarded::GuardedMeasureIdentity;
use super::normalization::native_measure_to_euclidean;
use super::reduction::{
    WeightedClosureOptions, WeightedClosureOutcome, WeightedReducedSystem, prepare_weighted_system,
};
use crate::engine::connection::{ConnectionRequest, ConnectionTransport};
use crate::{
    ComplexFloat, DifferentialSystem, Error, FlowOptions, MassMode, Precision, Prescription,
    RecursionMode, Result, RunContext, RustRedBackend,
};
use symbolica::prelude::*;

/// A source-replayed closed occupied connection, with its own measure and
/// continuation evidence. Its matrix is never accepted as caller-supplied data.
pub struct PreparedOccupiedFlow<const N: usize> {
    family: OccupiedCutFamily,
    epsilon: Symbol,
    dimension: i64,
    shifted: Vec<bool>,
    closed: WeightedReducedSystem<N>,
    system: Option<DifferentialSystem>,
    transport: ConnectionTransport,
}

fn validate_options(options: &FlowOptions) -> Result<()> {
    options.validate()?;
    if options.prescription != Prescription::PlusI0 || options.recursion != RecursionMode::Amf {
        return Err(Error::Unsupported(
            "occupied flow requires AMF with the admitted common +i0 prescription".into(),
        ));
    }
    if !matches!(options.mass_mode, MassMode::All | MassMode::Auto) {
        return Err(Error::Unsupported(
            "integrated occupied boundaries currently require all uncut quadratic factors shifted"
                .into(),
        ));
    }
    Ok(())
}

/// Sufficient open physical-mass domain, verified from incidence and charges.
/// No integral values or topology names enter this admission. For its one-cut
/// bubble F=x mb²+(1−x)M²−x(1−x)ma²>0; for the two-cut term,
/// M²−(q1−q2)² >= M²−(m1−m2)²>0. Adding eta>=0 preserves both.
fn positive_sunset_domain(input: &PreparedDensityInput) -> Result<String> {
    let graph = input.input();
    if graph.loops != 2
        || graph.vertices != 2
        || graph.edges.len() != 3
        || graph.edges.iter().any(|e| e.vertices[0] == e.vertices[1])
    {
        return Err(Error::Unsupported("numerical occupied continuation is currently certified for a two-loop three-edge graph with positive masses".into()));
    }
    let charged = graph
        .edges
        .iter()
        .enumerate()
        .filter_map(|(i, e)| e.charges.iter().any(|&q| q != 0).then_some(i))
        .collect::<Vec<_>>();
    if charged.len() != 2 {
        return Err(Error::Unsupported(
            "massive sunset continuation requires exactly two charged edges".into(),
        ));
    }
    let neutral = (0..3).find(|i| !charged.contains(i)).unwrap();
    let masses = input
        .physical_masses()
        .iter()
        .map(|m| {
            Rational::try_from(m.as_view()).map_err(|_| {
                Error::Unsupported(
                    "continuation needs assigned positive rational mass squares".into(),
                )
            })
        })
        .collect::<Result<Vec<_>>>()?;
    if masses.iter().any(|m| m <= &Rational::zero()) {
        return Err(Error::Unsupported("the certified sunset domain requires strictly positive physical masses; a massless limiting contour is not inferred".into()));
    }
    let delta = &masses[charged[0]] + &masses[charged[1]] - &masses[neutral];
    if delta >= Rational::zero()
        && &delta * &delta >= Rational::from(4) * &masses[charged[0]] * &masses[charged[1]]
    {
        return Err(Error::Unsupported("sunset continuation requires |m1-m2|<M; a threshold or pinch needs a different contour admission".into()));
    }
    Ok(format!(
        "regulated thermal contour; positive-mass open sunset domain |m1-m2|<M; masses²={masses:?}; charged={charged:?}; fixed-shell eta>=0; meromorphic dimension before Laurent expansion"
    ))
}

impl<const N: usize> PreparedOccupiedFlow<N> {
    pub fn prepare(
        input: &PreparedDensityInput,
        cuts: &[usize],
        options: &FlowOptions,
        closure_options: WeightedClosureOptions,
        context: &RunContext,
    ) -> Result<Self> {
        validate_options(options)?;
        context.cancellation.check()?;
        let admission = positive_sunset_domain(input)?;
        let family = input.occupied_cut(cuts, 65536)?.at_physical_masses();
        let epsilon = symbol!("rustflow_occupied::epsilon");
        let eta = symbol!("rustflow_occupied::eta");
        let shifted_slots = (0..family.physical_slots())
            .filter(|&i| family.roles()[i] == super::guarded::IndexRole::Ordinary)
            .collect::<Vec<_>>();
        let shifted = (0..family.factors().len())
            .map(|i| shifted_slots.contains(&i))
            .collect::<Vec<_>>();
        let sources = family.guarded_sources::<N>(epsilon,options.dimension,eta,&shifted_slots,65536,vec![],GuardedMeasureIdentity {
            measure:format!("independent occupied Cn and Hn; input={}; factors={:?}",input.identity(),family.factors()),
            support:format!("real compact shell coordinates; {:?}",family.shells()),
            orientation:format!("future first {} loops; q=(-i PE0,-PEvec); inverse routing={:?}; determinant={}",cuts.len(),family.inverse_routing(),family.routing_determinant()),
            normalization:"virtual d^Dk/(i pi^(D/2)); occupied d^Dq/pi^(D/2); Euclidean target Wick phases already included".into(),
            branch:admission,
            deformation:format!("native Di-eta on {shifted_slots:?}; physical shells and upper/lower supports fixed; input indices and coefficients fixed"),
        })?;
        let closed = match prepare_weighted_system(
            &sources.context,
            &sources.targets,
            &sources.deformation,
            closure_options,
            context,
        )? {
            WeightedClosureOutcome::Closed(closed) => closed,
            WeightedClosureOutcome::Unresolved(failure) => {
                return Err(Error::IncompleteReduction(format!(
                    "occupied weighted connection unresolved: {}; rounds={}, provisional={}, unresolved={:?}",
                    failure.reason,
                    failure.diagnostics.rounds,
                    failure.provisional_frontier.len(),
                    failure.unresolved
                )));
            }
        };
        let system = closed.differential_system();
        Ok(Self {
            family,
            epsilon,
            dimension: options.dimension,
            shifted,
            closed,
            system,
            transport: ConnectionTransport::default(),
        })
    }

    pub fn differential_system(&self) -> Option<&DifferentialSystem> {
        self.system.as_ref()
    }
    pub fn reduced(&self) -> &crate::reduction::ReducedSystem {
        &self.closed.reduced
    }
    pub fn closure_diagnostics(&self) -> &super::reduction::WeightedClosureDiagnostics {
        &self.closed.diagnostics
    }

    /// Euclidean occupied contribution, including target phases and the exact
    /// loop-routing Jacobian. This is one cut contribution, not a full vacuum
    /// plus density assembly. The returned target weights are projected jointly
    /// at eta=0 by the same rational endpoint owner as ordinary AMF.
    pub fn evaluate(
        &self,
        epsilon: &Rational,
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<Vec<ComplexFloat>> {
        self.evaluate_with_start_scale(epsilon, options, context, 8)
    }

    /// Independently move the asymptotic starting point while retaining the
    /// prepared source program, boundary order and endpoint matching point.
    /// Accuracy is assessed by recomputation, not inferred from working digits.
    pub fn evaluate_with_start_scale(
        &self,
        epsilon: &Rational,
        options: &FlowOptions,
        context: &RunContext,
        start_scale: u32,
    ) -> Result<Vec<ComplexFloat>> {
        validate_options(options)?;
        if options.dimension != self.dimension || epsilon.is_zero() {
            return Err(Error::InvalidInput(
                "occupied sampling requires nonzero epsilon and the prepared dimension".into(),
            ));
        }
        let p = Precision::decimal(
            options
                .digits
                .checked_add(options.guard_digits)
                .ok_or_else(|| Error::Limit("occupied precision overflow".into()))?,
        )?;
        // Even an exact-zero target keeps the original reduction/target pole
        // conditions. The shared evaluator checks them before its empty-basis
        // branch, so do not return early here.
        let zero_system = DifferentialSystem {
            variable: self.closed.variable,
            matrix: vec![],
        };
        let system = self.system.as_ref().unwrap_or(&zero_system);
        let backend = RustRedBackend {
            bubble_subloops: false,
            ..Default::default()
        };
        let boundary = OccupiedFlowBoundary::new(
            &backend,
            options,
            context,
            self.epsilon,
            options.series_order.min(100),
            OccupiedBoundaryLimits::default(),
        )?;
        let values = self.transport.evaluate(
            ConnectionRequest {
                system,
                reduced: &self.closed.reduced,
                epsilon: self.epsilon,
                loops: self.family.loops(),
                positive_mass_contour: true,
                restrict_to_right_half_plane: true,
                start_scale,
            },
            epsilon,
            options,
            context,
            |_| Ok(()),
            |solutions, p| {
                Ok(boundary
                    .constants(
                        &self.family,
                        &self.closed.reduced.basis,
                        &self.shifted,
                        epsilon,
                        p,
                        solutions,
                    )?
                    .constants)
            },
        )?;
        let d = Rational::from(self.dimension) - epsilon * &Rational::from(2);
        let det = self.family.routing_determinant();
        let abs_det = if det < &Rational::zero() {
            -det.clone()
        } else {
            det.clone()
        };
        let measure = native_measure_to_euclidean(
            self.family.loops(),
            self.family.shells().len(),
            &Atom::num(d.clone()),
        )? * Atom::num(abs_det).pow(Atom::num(-d));
        let factor = p.eval(&measure, &ahash::HashMap::default())?;
        Ok(values.iter().map(|v| p.mul(v, &factor)).collect())
    }
}
