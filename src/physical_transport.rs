//! Cache-first transport between regular physical kinematic points.
//!
//! Every request selects a compatible boundary before forming its path. Accepted
//! intermediate physical points are independently checked and retained. The
//! ordinary interface supports regular straight paths. An explicitly bound
//! prescribed interface reuses the same engine for threshold detours. Explicit
//! singular endpoint requests retain regular matching anchors and separately
//! validated terminal evidence.
use crate::algebraic::{
    AlgebraicKinematicSystem, AlgebraicSystem, CanonicalAlgebraicSystem, CompiledAlgebraicSystem,
    PreparedAlgebraicSystem,
};
use crate::diffexp::{EpsilonSolution, EpsilonSystem, transport_epsilon};
use crate::kinematics::{KinematicPath, KinematicSystem};
use crate::transport_cache::{
    BoundaryAccuracy, BoundaryIdentity, BoundaryQuery, CachedBoundary, CachedPoint, EpsilonRange,
    PhysicalContinuation, PointKind, RootGerm, RootSheet, RustFlowCache, TransportCost,
};
use crate::{Error, FlowOptions, Precision, Prescription, Result, RunContext};
use std::{cell::RefCell, collections::BTreeMap, rc::Rc};
use symbolica::prelude::*;

/// A planner-generated affine physical route. Vertices are exact binary
/// rational images of the certified rounded contour, shared by all refinements.
#[derive(Clone, Debug)]
pub struct PhysicalRoute {
    pub path: KinematicPath,
    pub waypoints: Vec<Atom>,
    pub crossings: Vec<crate::contour::ContourCrossing>,
}

/// Geometry alone does not establish the integral's global monodromy. The
/// caller must admit each candidate route into the identity's declared domain.
pub trait HomotopyAdmission {
    fn admit(
        &self,
        source: &CachedBoundary,
        destination: &CachedPoint,
        route: &PhysicalRoute,
    ) -> Result<bool>;
}
impl<F> HomotopyAdmission for F
where
    F: Fn(&CachedBoundary, &CachedPoint, &PhysicalRoute) -> Result<bool>,
{
    fn admit(
        &self,
        source: &CachedBoundary,
        destination: &CachedPoint,
        route: &PhysicalRoute,
    ) -> Result<bool> {
        self(source, destination, route)
    }
}

#[derive(Clone, Copy)]
enum RouteMode<'a> {
    Straight,
    Prescribed(&'a dyn HomotopyAdmission),
}
impl RouteMode<'_> {
    fn is_prescribed(self) -> bool {
        matches!(self, Self::Prescribed(_))
    }
}

pub struct RustFlow<S = KinematicSystem> {
    system: S,
    identity: BoundaryIdentity,
}

/// Outcome of an actually attempted compatible cached source.
#[derive(Clone, Debug)]
pub enum BoundaryAttemptOutcome {
    Accepted,
    AccuracyRejected { message: String },
}

#[derive(Clone, Debug)]
pub struct BoundaryAttempt {
    pub starting_point: CachedPoint,
    /// The arbitrary-precision cost returned by the caller's policy.
    pub cost: Float,
    pub source_verified_digits: u32,
    pub outcome: BoundaryAttemptOutcome,
}

pub struct PhysicalResult {
    pub boundary: CachedBoundary,
    pub starting_point: CachedPoint,
    /// None for a compatible cache hit at precisely the requested coordinates.
    pub transport: Option<EpsilonSolution>,
    pub inserted_points: usize,
    /// Compatible sources tried in selection order, including an exact hit.
    pub boundary_attempts: Vec<BoundaryAttempt>,
}

macro_rules! endpoint_api {
    ($system:ty, $kind:ident) => {
        impl RustFlow<$system> {
            /// Evaluate a finite coefficientwise singular limit through a
            /// regular matching anchor. The caller explicitly admits the final
            /// chart approach; terminal limits never initialize regular transport.
            #[allow(clippy::too_many_arguments)]
            pub fn evaluate_endpoint(
                &self,
                cache: &mut RustFlowCache,
                request: &crate::singular_endpoint::EndpointRequest,
                options: &FlowOptions,
                context: &RunContext,
                policy: &dyn TransportCost,
                admission: &dyn crate::singular_endpoint::EndpointAdmission,
            ) -> Result<crate::singular_endpoint::EndpointResult> {
                evaluate_endpoint(
                    Connection::$kind(&self.system),
                    &self.identity,
                    cache,
                    request,
                    None,
                    options,
                    context,
                    policy,
                    RouteMode::Straight,
                    admission,
                )
            }

            /// Reach the regular matching anchor along an explicitly admitted
            /// prescribed route, then admit and evaluate its endpoint chart.
            #[allow(clippy::too_many_arguments)]
            pub fn evaluate_prescribed_endpoint(
                &self,
                cache: &mut RustFlowCache,
                request: &crate::singular_endpoint::EndpointRequest,
                options: &FlowOptions,
                context: &RunContext,
                policy: &dyn TransportCost,
                route_admission: &dyn HomotopyAdmission,
                endpoint_admission: &dyn crate::singular_endpoint::EndpointAdmission,
            ) -> Result<crate::singular_endpoint::EndpointResult> {
                evaluate_endpoint(
                    Connection::$kind(&self.system),
                    &self.identity,
                    cache,
                    request,
                    None,
                    options,
                    context,
                    policy,
                    RouteMode::Prescribed(route_admission),
                    endpoint_admission,
                )
            }
            /// Evaluate a finite coefficientwise singular limit through a
            /// regular matching anchor. The caller explicitly admits the final
            /// chart approach; terminal limits never initialize regular transport.
            #[allow(clippy::too_many_arguments)]
            pub fn evaluate_constrained_endpoint(
                &self,
                cache: &mut RustFlowCache,
                request: &crate::singular_endpoint::EndpointRequest,
                constraints: &crate::singular_endpoint::EndpointConstraints,
                options: &FlowOptions,
                context: &RunContext,
                policy: &dyn TransportCost,
                admission: &dyn crate::singular_endpoint::EndpointAdmission,
            ) -> Result<crate::singular_endpoint::EndpointResult> {
                evaluate_endpoint(
                    Connection::$kind(&self.system),
                    &self.identity,
                    cache,
                    request,
                    Some(constraints),
                    options,
                    context,
                    policy,
                    RouteMode::Straight,
                    admission,
                )
            }

            /// Reach the regular matching anchor along an explicitly admitted
            /// prescribed route, then admit and evaluate its endpoint chart.
            #[allow(clippy::too_many_arguments)]
            pub fn evaluate_prescribed_constrained_endpoint(
                &self,
                cache: &mut RustFlowCache,
                request: &crate::singular_endpoint::EndpointRequest,
                constraints: &crate::singular_endpoint::EndpointConstraints,
                options: &FlowOptions,
                context: &RunContext,
                policy: &dyn TransportCost,
                route_admission: &dyn HomotopyAdmission,
                endpoint_admission: &dyn crate::singular_endpoint::EndpointAdmission,
            ) -> Result<crate::singular_endpoint::EndpointResult> {
                evaluate_endpoint(
                    Connection::$kind(&self.system),
                    &self.identity,
                    cache,
                    request,
                    Some(constraints),
                    options,
                    context,
                    policy,
                    RouteMode::Prescribed(route_admission),
                    endpoint_admission,
                )
            }
        }
    };
}
endpoint_api!(KinematicSystem, Rational);
endpoint_api!(AlgebraicKinematicSystem, Algebraic);
endpoint_api!(CanonicalAlgebraicSystem, Canonical);

#[allow(clippy::too_many_arguments)]
fn evaluate_endpoint(
    system: Connection<'_>,
    identity: &BoundaryIdentity,
    cache: &mut RustFlowCache,
    request: &crate::singular_endpoint::EndpointRequest,
    constraints: Option<&crate::singular_endpoint::EndpointConstraints>,
    options: &FlowOptions,
    context: &RunContext,
    policy: &dyn TransportCost,
    mode: RouteMode<'_>,
    admission: &dyn crate::singular_endpoint::EndpointAdmission,
) -> Result<crate::singular_endpoint::EndpointResult> {
    use crate::singular_endpoint::{EndpointResult, PreparedEndpoint};
    options.validate()?;
    request.preflight(identity)?;
    if let Some(constraints) = constraints {
        constraints.preflight(identity, request.range, context)?;
    }
    request.chart.validate(identity)?;
    if mode.is_prescribed() != identity.physical_continuation().is_some() {
        return Err(Error::InvalidInput("endpoint matching route must use the identity's prescribed/ordinary transport interface".into()));
    }
    context.cancellation.check()?;
    if let Some(boundary) = cache.endpoint(
        identity,
        &request.chart,
        request.range,
        options.digits,
        constraints,
    )? {
        if !admission.admit(&boundary.matching_boundary, &request.chart)? {
            return Err(Error::InvalidInput(
                "endpoint approach was not admitted".into(),
            ));
        }
        context.cancellation.check()?;
        return Ok(EndpointResult {
            boundary,
            matching_transport: None,
            inserted_regular_points: 0,
            cache_hit: true,
            boundary_attempts: Vec::new(),
        });
    }
    let order = (i64::from(request.range.last) - i64::from(request.range.leading)) as usize;
    let source = match system.pullback(&request.chart.path, order)? {
        PulledConnection::Rational(system) => AlgebraicSystem {
            system,
            roots: Vec::new(),
            nonzero_conditions: Vec::new(),
        },
        PulledConnection::Algebraic(system) => system,
    };
    let prepared = PreparedEndpoint::new(source, identity, request, constraints, options, context)?;
    evaluate_physical_with(
        system,
        identity,
        cache,
        request.chart.matching_point()?,
        request.range,
        options,
        context,
        policy,
        mode,
        |regular| {
            if !admission.admit(&regular.boundary, &request.chart)? {
                return Err(Error::InvalidInput(
                    "endpoint approach was not admitted".into(),
                ));
            }
            let boundary = prepared.evaluate(regular.boundary, options, context)?;
            Ok((
                EndpointResult {
                    boundary: boundary.clone(),
                    matching_transport: regular.transport,
                    inserted_regular_points: regular.inserted_points,
                    cache_hit: false,
                    boundary_attempts: regular.boundary_attempts,
                },
                vec![boundary],
            ))
        },
    )
}

impl RustFlow {
    /// Regularize rational physical matrices with exact diagonal epsilon
    /// powers. The returned adapter preserves the original identity for
    /// checked boundary conversion and uses this same transport/cache engine.
    pub fn regularize_epsilon(&self, context: &RunContext) -> Result<crate::EpsilonShearedFlow> {
        let (system, shearing) = crate::EpsilonShearing::regularize(&self.system, context)?;
        let identity = self.identity.epsilon_sheared(&system, &shearing)?;
        Ok(crate::EpsilonShearedFlow {
            original_identity: self.identity.clone(),
            flow: Self { system, identity },
            shearing,
        })
    }

    /// Continue through planner-generated threshold detours. The identity must
    /// be bound to exact prescriptions and the caller must admit its homotopy.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_prescribed_to(
        &self,
        cache: &mut RustFlowCache,
        destination: &BTreeMap<Symbol, Atom>,
        range: EpsilonRange,
        options: &FlowOptions,
        context: &RunContext,
        policy: &dyn TransportCost,
        admission: &dyn HomotopyAdmission,
    ) -> Result<PhysicalResult> {
        let target = CachedPoint::Exact(destination.clone());
        evaluate_physical(
            Connection::Rational(&self.system),
            &self.identity,
            cache,
            target,
            range,
            options,
            context,
            policy,
            RouteMode::Prescribed(admission),
        )
    }

    pub fn new(
        system: KinematicSystem,
        basis: &[Atom],
        normalization: &Atom,
        prescription: Prescription,
        branch_domain: &str,
    ) -> Result<Self> {
        Self::with_conditions(
            system,
            basis,
            normalization,
            prescription,
            branch_domain,
            &[],
        )
    }

    /// Preserve all exact assumptions used in deriving the physical system.
    pub fn with_conditions(
        system: KinematicSystem,
        basis: &[Atom],
        normalization: &Atom,
        prescription: Prescription,
        branch_domain: &str,
        conditions: &[Atom],
    ) -> Result<Self> {
        let identity = BoundaryIdentity::with_conditions(
            &system,
            basis,
            normalization,
            prescription,
            branch_domain,
            conditions,
        )?;
        Ok(Self { system, identity })
    }

    /// Select a verified cached boundary and transport to this exact point.
    /// The cost policy must validate physical path/sheet compatibility. Nearby
    /// entries are starting values, never replacements for the requested point.
    /// The entire cache update is committed only after successful validation.
    pub fn evaluate_to(
        &self,
        cache: &mut RustFlowCache,
        destination: &BTreeMap<Symbol, Atom>,
        range: EpsilonRange,
        options: &FlowOptions,
        context: &RunContext,
        policy: &dyn TransportCost,
    ) -> Result<PhysicalResult> {
        evaluate_physical(
            Connection::Rational(&self.system),
            &self.identity,
            cache,
            CachedPoint::Exact(destination.clone()),
            range,
            options,
            context,
            policy,
            RouteMode::Straight,
        )
    }
}

impl<S> RustFlow<S> {
    /// Bind exact prescriptions into a distinct typed cache identity. A bound
    /// flow must use `evaluate_prescribed_to` with explicit homotopy admission.
    pub fn with_prescribed_continuation(
        mut self,
        continuation: PhysicalContinuation,
    ) -> Result<Self> {
        self.identity = self.identity.with_prescribed_continuation(continuation)?;
        Ok(self)
    }
    /// The exact common-basis physical connection.
    pub fn system(&self) -> &S {
        &self.system
    }
    pub fn identity(&self) -> &BoundaryIdentity {
        &self.identity
    }
}

impl RustFlow<AlgebraicKinematicSystem> {
    /// Apply the shared exact epsilon gauge after native root normalization.
    pub fn regularize_epsilon(
        &self,
        context: &RunContext,
    ) -> Result<crate::EpsilonShearedFlow<AlgebraicKinematicSystem>> {
        let (system, shearing) =
            crate::EpsilonShearing::regularize_algebraic(&self.system, context)?;
        let identity = self
            .identity
            .epsilon_sheared_algebraic(&system, &shearing)?;
        Ok(crate::EpsilonShearedFlow {
            original_identity: self.identity.clone(),
            flow: Self { system, identity },
            shearing,
        })
    }

    /// Continue through planner-generated threshold detours. The identity must
    /// be bound to exact prescriptions and the caller must admit its homotopy.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_prescribed_to(
        &self,
        cache: &mut RustFlowCache,
        destination: &BTreeMap<Symbol, Atom>,
        germ: &RootGerm,
        range: EpsilonRange,
        options: &FlowOptions,
        context: &RunContext,
        policy: &dyn TransportCost,
        admission: &dyn HomotopyAdmission,
    ) -> Result<PhysicalResult> {
        let target = CachedPoint::Exact(destination.clone()).with_root_germ(germ.clone())?;
        evaluate_physical(
            Connection::Algebraic(&self.system),
            &self.identity,
            cache,
            target,
            range,
            options,
            context,
            policy,
            RouteMode::Prescribed(admission),
        )
    }

    pub fn new_algebraic(
        system: AlgebraicKinematicSystem,
        basis: &[Atom],
        normalization: &Atom,
        prescription: Prescription,
        branch_domain: &str,
    ) -> Result<Self> {
        Self::with_algebraic_conditions(
            system,
            basis,
            normalization,
            prescription,
            branch_domain,
            &[],
        )
    }
    pub fn with_algebraic_conditions(
        system: AlgebraicKinematicSystem,
        basis: &[Atom],
        normalization: &Atom,
        prescription: Prescription,
        branch_domain: &str,
        conditions: &[Atom],
    ) -> Result<Self> {
        if system.roots.is_empty() {
            return Err(Error::InvalidInput(
                "use the rational RustFlow constructor for a system without registered roots"
                    .into(),
            ));
        }
        let identity = BoundaryIdentity::with_algebraic_conditions(
            &system,
            basis,
            normalization,
            prescription,
            branch_domain,
            conditions,
        )?;
        Ok(Self { system, identity })
    }
    /// Transport along a regular affine path with exact rational-complex coordinates.
    /// Every radicand must remain nonzero; the requested local root germ is explicit. The caller's
    /// admissibility policy additionally controls integral/logarithmic monodromy.
    #[allow(clippy::too_many_arguments)] // Ordinary physical inputs plus the required destination root germ.
    pub fn evaluate_to(
        &self,
        cache: &mut RustFlowCache,
        destination: &BTreeMap<Symbol, Atom>,
        germ: &RootGerm,
        range: EpsilonRange,
        options: &FlowOptions,
        context: &RunContext,
        policy: &dyn TransportCost,
    ) -> Result<PhysicalResult> {
        evaluate_physical(
            Connection::Algebraic(&self.system),
            &self.identity,
            cache,
            CachedPoint::Exact(destination.clone()).with_root_germ(germ.clone())?,
            range,
            options,
            context,
            policy,
            RouteMode::Straight,
        )
    }
}

impl RustFlow<CanonicalAlgebraicSystem> {
    /// Continue through planner-generated threshold detours. The identity must
    /// be bound to exact prescriptions and the caller must admit its homotopy.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_prescribed_to(
        &self,
        cache: &mut RustFlowCache,
        destination: &BTreeMap<Symbol, Atom>,
        germ: Option<&RootGerm>,
        range: EpsilonRange,
        options: &FlowOptions,
        context: &RunContext,
        policy: &dyn TransportCost,
        admission: &dyn HomotopyAdmission,
    ) -> Result<PhysicalResult> {
        let target = match (self.system.roots().is_empty(), germ) {
            (true, None) => CachedPoint::Exact(destination.clone()),
            (false, Some(germ)) => {
                CachedPoint::Exact(destination.clone()).with_root_germ(germ.clone())?
            }
            _ => {
                return Err(Error::InvalidInput(
                    "canonical root germ must be supplied exactly when roots are registered".into(),
                ));
            }
        };
        evaluate_physical(
            Connection::Canonical(&self.system),
            &self.identity,
            cache,
            target,
            range,
            options,
            context,
            policy,
            RouteMode::Prescribed(admission),
        )
    }

    pub fn new_canonical(
        system: CanonicalAlgebraicSystem,
        basis: &[Atom],
        normalization: &Atom,
        prescription: Prescription,
        branch_domain: &str,
    ) -> Result<Self> {
        Self::with_canonical_conditions(
            system,
            basis,
            normalization,
            prescription,
            branch_domain,
            &[],
        )
    }
    /// Keep the canonical form separate until a physical path is selected.
    pub fn with_canonical_conditions(
        system: CanonicalAlgebraicSystem,
        basis: &[Atom],
        normalization: &Atom,
        prescription: Prescription,
        branch_domain: &str,
        conditions: &[Atom],
    ) -> Result<Self> {
        let identity = BoundaryIdentity::with_canonical_conditions(
            &system,
            basis,
            normalization,
            prescription,
            branch_domain,
            conditions,
        )?;
        Ok(Self { system, identity })
    }
    /// Use an explicit germ when roots are registered, and None otherwise.
    /// The same regular affine path and caller monodromy contract as the dense
    /// algebraic interface applies; no dense physical matrix is assembled.
    #[allow(clippy::too_many_arguments)] // Physical transport inputs plus explicit optional root germ.
    pub fn evaluate_to(
        &self,
        cache: &mut RustFlowCache,
        destination: &BTreeMap<Symbol, Atom>,
        germ: Option<&RootGerm>,
        range: EpsilonRange,
        options: &FlowOptions,
        context: &RunContext,
        policy: &dyn TransportCost,
    ) -> Result<PhysicalResult> {
        let point = CachedPoint::Exact(destination.clone());
        let point = match (self.system.roots().is_empty(), germ) {
            (true, None) => point,
            (false, Some(germ)) => point.with_root_germ(germ.clone())?,
            _ => {
                return Err(Error::InvalidInput(
                    "canonical root germ must be supplied exactly when roots are registered".into(),
                ));
            }
        };
        evaluate_physical(
            Connection::Canonical(&self.system),
            &self.identity,
            cache,
            point,
            range,
            options,
            context,
            policy,
            RouteMode::Straight,
        )
    }
}

#[derive(Clone, Copy)]
enum Connection<'a> {
    Rational(&'a KinematicSystem),
    Algebraic(&'a AlgebraicKinematicSystem),
    Canonical(&'a CanonicalAlgebraicSystem),
}
impl Connection<'_> {
    fn pullback(&self, path: &KinematicPath, order: usize) -> Result<PulledConnection> {
        Ok(match self {
            Self::Rational(system) => {
                PulledConnection::Rational(EpsilonSystem::from_differential_system(
                    &system.pullback(path)?,
                    system.epsilon,
                    order,
                )?)
            }
            Self::Algebraic(system) => PulledConnection::Algebraic(system.pullback(path, order)?),
            Self::Canonical(system) => PulledConnection::Algebraic(system.pullback(path, order)?),
        })
    }
}
enum PulledConnection {
    Rational(EpsilonSystem),
    Algebraic(AlgebraicSystem),
}
impl PulledConnection {
    fn prepare(self, context: &RunContext) -> Result<PreparedConnection> {
        context.cancellation.check()?;
        Ok(match self {
            Self::Rational(system) => PreparedConnection::Rational(system),
            Self::Algebraic(system) => {
                PreparedConnection::Algebraic(system.prepare_with_context(context)?)
            }
        })
    }
}
enum PreparedConnection {
    Rational(EpsilonSystem),
    Algebraic(PreparedAlgebraicSystem),
}
impl PreparedConnection {
    fn compile(&self, p: Precision, context: &RunContext) -> Result<CompiledConnection> {
        context.cancellation.check()?;
        let compiled = match self {
            Self::Rational(system) => {
                CompiledConnection::Rational(system.compile(p, &Default::default())?)
            }
            Self::Algebraic(system) => {
                CompiledConnection::Algebraic(system.compile_with_context(p, context)?)
            }
        };
        context.cancellation.check()?;
        Ok(compiled)
    }
    fn transport(
        &self,
        source: &CachedBoundary,
        target: &CachedPoint,
        range: EpsilonRange,
        waypoints: &[Atom],
        options: &FlowOptions,
        context: &RunContext,
    ) -> Result<(EpsilonSolution, BTreeMap<usize, RootGerm>)> {
        match self {
            Self::Rational(system) => Ok((
                transport_epsilon(
                    system,
                    |p| source.as_epsilon_boundary(&Atom::new(), range, p),
                    waypoints,
                    options,
                    context,
                    true,
                )?,
                BTreeMap::new(),
            )),
            Self::Algebraic(system) => transport_algebraic_checked(
                system, source, target, range, waypoints, options, context,
            ),
        }
    }
}
enum CompiledConnection {
    Rational(crate::diffexp::CompiledEpsilonSystem),
    Algebraic(CompiledAlgebraicSystem),
}
impl CompiledConnection {
    fn poles(&self) -> &[crate::ComplexFloat] {
        match self {
            Self::Rational(c) => c.poles(),
            Self::Algebraic(c) => c.singularities(),
        }
    }
    fn singularity_polynomials(&self) -> &[Atom] {
        match self {
            Self::Rational(c) => c.singularity_polynomials(),
            Self::Algebraic(c) => c.singularity_polynomials(),
        }
    }
    fn epsilon_product_limit(&self) -> Option<usize> {
        match self {
            Self::Rational(c) => c.epsilon_product_limit(),
            Self::Algebraic(c) => c.epsilon_product_limit(),
        }
    }
    fn error_norm_integral_weighted(
        &self,
        start: &crate::ComplexFloat,
        end: &crate::ComplexFloat,
        weights: &[Float],
    ) -> Result<Float> {
        match self {
            Self::Rational(c) => c.error_norm_integral_weighted(start, end, weights),
            Self::Algebraic(c) => c.error_norm_integral_weighted(start, end, weights),
        }
    }
}

#[allow(clippy::too_many_arguments)] // Shared orchestration for both exact connection kinds.
fn evaluate_physical(
    system: Connection<'_>,
    identity: &BoundaryIdentity,
    cache: &mut RustFlowCache,
    target: CachedPoint,
    range: EpsilonRange,
    options: &FlowOptions,
    context: &RunContext,
    policy: &dyn TransportCost,
    mode: RouteMode<'_>,
) -> Result<PhysicalResult> {
    evaluate_physical_with(
        system,
        identity,
        cache,
        target,
        range,
        options,
        context,
        policy,
        mode,
        |result| Ok((result, Vec::new())),
    )
}

#[allow(clippy::too_many_arguments)]
fn evaluate_physical_with<T>(
    system: Connection<'_>,
    identity: &BoundaryIdentity,
    cache: &mut RustFlowCache,
    target: CachedPoint,
    range: EpsilonRange,
    options: &FlowOptions,
    context: &RunContext,
    policy: &dyn TransportCost,
    mode: RouteMode<'_>,
    mut complete: impl FnMut(
        PhysicalResult,
    ) -> Result<(T, Vec<crate::singular_endpoint::EndpointBoundary>)>,
) -> Result<T> {
    options.validate()?;
    if mode.is_prescribed() != identity.physical_continuation().is_some() {
        return Err(Error::InvalidInput("prescribed identities require evaluate_prescribed_to with explicit homotopy admission; ordinary identities require evaluate_to".into()));
    }
    context.cancellation.check()?;
    let p = Precision::decimal(options.digits + options.guard_digits)?;
    let query = BoundaryQuery::new(identity, &target, range, options.digits)?;
    // Reduction domains are enforced independently of the caller's sheet
    // policy. Reject unsafe sources before ranking, so another valid source
    // can be selected rather than failing after a nearest-source choice.
    let guarded_policy = GuardedCost {
        identity,
        policy,
        digits: options.digits,
        context,
        mode,
        system,
        range,
        planned: RefCell::new(BTreeMap::new()),
    };
    let mut excluded = std::collections::BTreeSet::new();
    let mut attempts = Vec::new();
    let mut last_accuracy = None;
    for _ in 0..options.max_boundary_attempts {
        context.cancellation.check()?;
        let Some((index, matched)) = cache.best_excluding_with_germ_policy(
            &query,
            &guarded_policy,
            p,
            &excluded,
            mode.is_prescribed() || !identity.roots().is_empty(),
        )?
        else {
            return Err(if let Some(message) = last_accuracy {
                Error::Accuracy(format!(
                    "none of {} compatible cached boundaries reached the target accuracy; last attempt: {message}",
                    attempts.len()
                ))
            } else {
                Error::IncompleteReduction("no compatible cached physical boundary with the requested epsilon range, verified accuracy and regular reduction domain".into())
            });
        };
        let source = matched.boundary.clone();
        let mut diagnostic = BoundaryAttempt {
            starting_point: source.point.clone(),
            cost: matched.cost,
            source_verified_digits: source.accuracy.verified_digits(),
            outcome: BoundaryAttemptOutcome::Accepted,
        };
        excluded.insert(index);
        context.cancellation.check()?;
        match evaluate_from_source(
            system,
            identity,
            target.clone(),
            source,
            range,
            options,
            context,
            guarded_policy
                .planned
                .borrow()
                .get(&diagnostic.starting_point.key()?)
                .cloned(),
        )
        .and_then(|(mut result, pending)| {
            result.boundary_attempts = attempts.clone();
            result.boundary_attempts.push(diagnostic.clone());
            let (output, endpoints) = complete(result)?;
            context.cancellation.check()?;
            cache.insert_batch(pending, endpoints)?;
            Ok(output)
        }) {
            Ok(result) => return Ok(result),
            Err(error @ (Error::Accuracy(_) | Error::InsufficientPrecision { .. })) => {
                // A different cached boundary can avoid the conditioning that
                // exhausted this route's bounded precision retries. Preserve
                // the attempt cap and never commit its failed trajectory.
                let message = match error {
                    Error::Accuracy(message) => message,
                    error => error.to_string(),
                };
                diagnostic.outcome = BoundaryAttemptOutcome::AccuracyRejected {
                    message: message.clone(),
                };
                attempts.push(diagnostic);
                last_accuracy = Some(message);
            }
            Err(error) => return Err(error),
        }
    }
    Err(Error::Limit(format!(
        "physical boundary attempt budget {} exhausted after accuracy failures; last attempt: {}",
        options.max_boundary_attempts,
        last_accuracy.unwrap_or_default()
    )))
}

#[allow(clippy::too_many_arguments)] // One already selected source and the shared physical transport inputs.
fn evaluate_from_source(
    system: Connection<'_>,
    identity: &BoundaryIdentity,
    target: CachedPoint,
    source: CachedBoundary,
    range: EpsilonRange,
    options: &FlowOptions,
    context: &RunContext,
    supplied_route: Option<Rc<PreparedRoute>>,
) -> Result<(PhysicalResult, Vec<CachedBoundary>)> {
    let p = Precision::decimal(options.digits + options.guard_digits)?;
    let destination = target.restart_coordinates()?;
    let coordinates = source.point.rounded_coordinates_as_exact()?;
    let count = (i64::from(range.last) - i64::from(range.leading) + 1) as usize;
    if source.point.root_germ() == target.root_germ() && source.point.same_coordinates(&target)? {
        let mut boundary = source.clone();
        boundary.range = range;
        boundary.coefficients.truncate(count);
        boundary.accuracy = source.accuracy.reindexed(
            source.accuracy.comparison_errors()[..count].to_vec(),
            "compatible exact-coordinate cache hit",
        )?;
        return Ok((
            PhysicalResult {
                boundary,
                starting_point: source.point,
                transport: None,
                inserted_points: 0,
                boundary_attempts: Vec::new(),
            },
            Vec::new(),
        ));
    }
    let route = if let Some(route) = supplied_route {
        route
    } else {
        let path = KinematicPath::straight_line(
            identity.path_parameter("physical_path_parameter")?,
            &coordinates,
            &destination,
        )?;
        if !identity.conditions_admit_straight_path_with_context(
            &source.point,
            &target,
            p,
            options.digits,
            context,
        )? {
            return Err(Error::Unsupported(
                "selected physical path leaves the reduction domain".into(),
            ));
        }
        let system = system.pullback(&path, count - 1)?.prepare(context)?;
        let prepared = system.compile(p, context)?;
        if prepared.poles().iter().any(|pole| {
            pole.re >= p.real(0)
                && pole.re <= p.real(1)
                && p.norm(&crate::ComplexFloat::new(p.real(0), pole.im.clone()))
                    <= p.tolerance(options.digits + 8)
        }) {
            return Err(Error::Unsupported("selected physical straight path meets a singularity; use explicit prescribed physical transport".into()));
        }
        Rc::new(PreparedRoute {
            physical: PhysicalRoute {
                path,
                waypoints: vec![Atom::one()],
                crossings: Vec::new(),
            },
            system,
        })
    };
    let path = &route.physical.path;
    let system = &route.system;
    let (mut solution, checkpoint_germs) = system.transport(
        &source,
        &target,
        range,
        &route.physical.waypoints,
        options,
        context,
    )?;
    let p = Precision {
        bits: solution.diagnostics.working_bits,
    };
    let prepared = system.compile(p, context)?;
    // A second transport from the same cached values checks integration
    // error only. Carry the supplied boundary uncertainty forward too.
    // Fix one scale per coefficient for the whole trajectory. A large
    // higher epsilon coefficient then contributes only through actual
    // couplings, rather than imposing its absolute error on every order.
    let weights = source.coefficients[..count]
        .iter()
        .flatten()
        .map(|value| {
            let norm = p.norm(value);
            if norm > p.real(1) { norm } else { p.real(1) }
        })
        .collect::<Vec<_>>();
    let mut input_relative_error = p.real(0);
    for (error, weight) in source.accuracy.comparison_errors()[..count]
        .iter()
        .flatten()
        .zip(&weights)
    {
        let relative = error.clone() / weight;
        if relative > input_relative_error {
            input_relative_error = relative;
        }
    }
    input_relative_error += p.tolerance(source.accuracy.verified_digits());
    let evidence_cap = source
        .accuracy
        .verified_digits()
        .min((options.digits + options.guard_digits).saturating_sub(10))
        .min(
            solution
                .diagnostics
                .conditioning_digits
                .unwrap_or(options.digits),
        );
    let mut integrated_norm = p.real(0);
    let mut amplification = p.real(1);
    let epsilon_product_limit = prepared.epsilon_product_limit();
    let mut accepted = Vec::new();
    for (index, segment) in solution.segments.iter().enumerate() {
        context.cancellation.check()?;
        integrated_norm +=
            prepared.error_norm_integral_weighted(&segment.center, &segment.end, &weights)?;
        amplification = crate::diffexp::amplification_from_integral(
            p,
            &integrated_norm,
            epsilon_product_limit,
        )?;
        if let Some(checkpoint) = solution.checkpoints.iter().find(|c| c.segment == index) {
            let mut checkpoint = checkpoint.clone();
            for (error, weight) in checkpoint
                .comparison_errors
                .iter_mut()
                .flatten()
                .zip(&weights)
            {
                *error += weight.clone() * &input_relative_error * &amplification;
            }
            if errors_meet(
                p,
                &checkpoint.coefficients,
                &checkpoint.comparison_errors,
                options.digits,
            ) {
                accepted.push(checkpoint);
            }
        }
    }
    for (error, weight) in solution
        .comparison_errors
        .iter_mut()
        .flatten()
        .zip(&weights)
    {
        *error += weight.clone() * &input_relative_error * &amplification;
    }
    if !errors_meet(
        p,
        &solution.coefficients,
        &solution.comparison_errors,
        options.digits,
    ) {
        return Err(Error::Accuracy("cached boundary uncertainty grows beyond the requested target accuracy; provide a more accurate or closer boundary".into()));
    }
    solution.verified_digits = Some(strongest_evidence(
        p,
        &solution.coefficients,
        &solution.comparison_errors,
        options.digits,
        evidence_cap,
    ));
    solution.checkpoints = accepted;
    let boundary = CachedBoundary {
        identity: identity.clone(),
        point: target,
        kind: PointKind::Physical,
        range,
        coefficients: solution.coefficients.clone(),
        accuracy: BoundaryAccuracy::from_solution(
            &solution,
            source.accuracy.verified_digits(),
            &format!(
                "independently refined physical transport including propagated cached input uncertainty; source: {}",
                source.accuracy.provenance(),
            ),
        )?,
    };
    boundary.validate()?;
    let mut pending = RustFlowCache::trajectory_boundaries_with_germs(
        identity,
        path,
        &solution,
        |index| Ok(checkpoint_germs.get(&index).cloned()),
        |index, segment| {
            if segment.end.im != p.real(0) {
                return Ok(None);
            }
            let Some(checkpoint) = solution.checkpoints.iter().find(|c| c.segment == index) else {
                return Ok(None);
            };
            let checked = strongest_evidence(
                p,
                &checkpoint.coefficients,
                &checkpoint.comparison_errors,
                options.digits,
                evidence_cap,
            );
            let accuracy = BoundaryAccuracy::supplied(
                checked,
                p.bits,
                checkpoint.comparison_errors.clone(),
                &format!(
                    "independent precision/order check at trajectory endpoint, with propagated cached input uncertainty; source: {}",
                    source.accuracy.provenance(),
                ),
            )?;
            Ok(Some((PointKind::Physical, accuracy)))
        },
    )?;
    let inserted = pending.len() + 1;
    pending.push(boundary.clone());
    Ok((
        PhysicalResult {
            boundary,
            starting_point: source.point,
            transport: Some(solution),
            inserted_points: inserted,
            boundary_attempts: Vec::new(),
        },
        pending,
    ))
}

pub(crate) fn errors_meet(
    p: Precision,
    values: &[Vec<crate::ComplexFloat>],
    errors: &[Vec<Float>],
    digits: u32,
) -> bool {
    values
        .iter()
        .flatten()
        .zip(errors.iter().flatten())
        .all(|(v, e)| {
            let norm = p.norm(v);
            let scale = if norm > p.real(1) { norm } else { p.real(1) };
            e.is_finite() && *e <= p.tolerance(digits) * scale
        })
}

// Stronger cache evidence requires explicitly passing tighter comparisons and
// propagated input errors. Arithmetic precision supplies only a conservative
// ceiling. Keep two decimal digits beyond any accuracy advertised to the bank.
pub(crate) fn strongest_evidence(
    p: Precision,
    values: &[Vec<crate::ComplexFloat>],
    errors: &[Vec<Float>],
    requested: u32,
    cap: u32,
) -> u32 {
    let mut checked = requested;
    for digits in requested.saturating_add(1)..=cap.saturating_sub(2) {
        if !errors_meet(p, values, errors, digits.saturating_add(2)) {
            break;
        }
        checked = digits;
    }
    checked
}

/// Compatibility name for the physical transport engine.
pub type PhysicalTransport = RustFlow;

struct GuardedCost<'a> {
    identity: &'a BoundaryIdentity,
    policy: &'a dyn TransportCost,
    digits: u32,
    context: &'a RunContext,
    mode: RouteMode<'a>,
    system: Connection<'a>,
    range: EpsilonRange,
    planned: RefCell<BTreeMap<String, Rc<PreparedRoute>>>,
}
impl TransportCost for GuardedCost<'_> {
    fn cost(
        &self,
        source: &CachedBoundary,
        target: &CachedPoint,
        p: Precision,
    ) -> Result<Option<Float>> {
        self.context.cancellation.check()?;
        match self.mode {
            RouteMode::Straight => {
                if !self.identity.conditions_admit_straight_path_with_context(
                    &source.point,
                    target,
                    p,
                    self.digits,
                    self.context,
                )? {
                    return Ok(None);
                }
            }
            RouteMode::Prescribed(admission) => {
                // A constant affine chart has no loop and cannot change sheet.
                // An exact compatible hit needs no homotopy admission, but the
                // caller's source policy still applies. Opposite sheets must
                // not let a zero distance outrank a valid source.
                if source.point.same_coordinates(target)? {
                    if source.point.root_germ() != target.root_germ() {
                        return Ok(None);
                    }
                    return self.policy.cost(source, target, p);
                }
                let route = prepare_prescribed_route(
                    self.system,
                    self.identity,
                    source,
                    target,
                    self.range,
                    p,
                    self.context,
                )?;
                if !admission.admit(source, target, &route.physical)? {
                    return Ok(None);
                }
                self.context.cancellation.check()?;
                self.planned
                    .borrow_mut()
                    .insert(source.point.key()?, Rc::new(route));
            }
        }
        self.policy.cost(source, target, p)
    }
    fn lower_bound(
        &self,
        source: &CachedBoundary,
        target: &CachedPoint,
        p: Precision,
    ) -> Result<Option<Float>> {
        self.context.cancellation.check()?;
        self.policy.lower_bound(source, target, p)
    }
    fn compare_tied_costs(
        &self,
        left: &CachedBoundary,
        right: &CachedBoundary,
        target: &CachedPoint,
        p: Precision,
    ) -> Result<std::cmp::Ordering> {
        self.context.cancellation.check()?;
        self.policy.compare_tied_costs(left, right, target, p)
    }
}

struct PreparedRoute {
    physical: PhysicalRoute,
    system: PreparedConnection,
}

#[allow(clippy::too_many_arguments)]
fn prepare_prescribed_route(
    system: Connection<'_>,
    identity: &BoundaryIdentity,
    source: &CachedBoundary,
    target: &CachedPoint,
    range: EpsilonRange,
    p: Precision,
    context: &RunContext,
) -> Result<PreparedRoute> {
    identity.validate_point(&source.point)?;
    identity.validate_point(target)?;
    let start = source.point.restart_coordinates()?;
    let finish = target.restart_coordinates()?;
    for value in start.values().chain(finish.values()) {
        crate::transport_cache::exact_real_rational(value)?;
    }
    let path = KinematicPath::straight_line(
        identity.path_parameter("prescribed_physical_path")?,
        &start,
        &finish,
    )?;
    let system = system
        .pullback(
            &path,
            (i64::from(range.last) - i64::from(range.leading)) as usize,
        )?
        .prepare(context)?;
    let compiled = system.compile(p, context)?;
    let mut polynomials = compiled.singularity_polynomials().to_vec();
    polynomials.extend(identity.prescribed_path_conditions(&path)?);
    let convention = identity
        .physical_continuation()
        .ok_or_else(|| Error::InvalidInput("missing prescribed identity".into()))?;
    let rules = path
        .coordinates
        .iter()
        .map(|(&s, a)| (Atom::var(s), a.clone()))
        .collect();
    let planner = crate::contour::PrescribedContour {
        variable: path.parameter,
        prescriptions: convention
            .prescriptions
            .iter()
            .map(|d| crate::contour::PolynomialPrescription {
                polynomial: crate::family::substitute(&d.polynomial, &rules)
                    .together()
                    .cancel(),
                prescription: d.prescription,
            })
            .collect(),
        unprescribed_side: convention.unprescribed_side,
    };
    let contour = planner.plan_with_context(p, &polynomials, &p.zero(), &p.i(1), context)?;
    let waypoints = contour
        .waypoints
        .iter()
        .map(|x| Atom::num(Complex::new(x.re.to_rational(), x.im.to_rational())))
        .collect();
    Ok(PreparedRoute {
        physical: PhysicalRoute {
            path,
            waypoints,
            crossings: contour.crossings,
        },
        system,
    })
}

struct BranchTrace {
    compiled: CompiledAlgebraicSystem,
    segments: Vec<crate::algebraic::BranchSegment>,
}
fn classify_germ(
    compiled: &CompiledAlgebraicSystem,
    state: &crate::algebraic::BranchState,
) -> Result<RootGerm> {
    let seeds = state
        .roots
        .iter()
        .map(|(s, _)| (*s, crate::algebraic::RootSeed::Principal))
        .collect();
    let principal = compiled.branch_state_at(&state.point, &seeds)?;
    let p = Precision {
        bits: state.point.re.as_raw().prec(),
    };
    if state.roots.len() != principal.roots.len() {
        return Err(Error::InvalidInput("root registry changed".into()));
    }
    let mut sheets = BTreeMap::new();
    for ((symbol, value), (expected_symbol, positive)) in state.roots.iter().zip(principal.roots) {
        if *symbol != expected_symbol || !p.finite(value) {
            return Err(Error::InvalidInput(
                "invalid transported root identity/value".into(),
            ));
        }
        let scale = p.norm(&positive);
        let plus = p.norm(&p.sub(value, &positive));
        let minus = p.norm(&p.add(value, &positive));
        let sheet = if plus * 4 < scale {
            RootSheet::Principal
        } else if minus * 4 < scale {
            RootSheet::Opposite
        } else {
            return Err(Error::Accuracy(
                "transported root sheet is ambiguous".into(),
            ));
        };
        sheets.insert(*symbol, sheet);
    }
    Ok(RootGerm { sheets })
}

#[allow(clippy::too_many_arguments)]
fn transport_algebraic_checked(
    system: &PreparedAlgebraicSystem,
    source: &CachedBoundary,
    target: &CachedPoint,
    range: EpsilonRange,
    waypoints: &[Atom],
    options: &FlowOptions,
    context: &RunContext,
) -> Result<(EpsilonSolution, BTreeMap<usize, RootGerm>)> {
    let seeds = source
        .point
        .root_germ()
        .map(RootGerm::seeds)
        .unwrap_or_default();
    let (solution, trace) = crate::diffexp::refine_epsilon_transport_with_metadata(
        options,
        context,
        true,
        |p, refined| {
            let compiled = system.compile_with_context(p, context)?;
            let points = waypoints
                .iter()
                .map(|a| p.eval(a, &Default::default()))
                .collect::<Result<Vec<_>>>()?;
            let boundary = source.as_epsilon_boundary(&Atom::new(), range, p)?;
            let coordinate_digits = crate::ode::conditioning::exact_waypoint_conditioning(
                p,
                &boundary.point,
                waypoints,
                &points,
                options.digits,
            )?;
            let mut result =
                compiled.transport(&boundary, &points, &seeds, refined, context, true)?;
            result.solution.diagnostics.conditioning_digits = Some(
                result
                    .solution
                    .diagnostics
                    .conditioning_digits
                    .map_or(coordinate_digits, |old| old.min(coordinate_digits)),
            );
            if !system.source().roots.is_empty()
                && Some(&classify_germ(&compiled, &result.branches)?) != target.root_germ()
            {
                return Err(Error::InvalidInput("planned continuation reaches a different root germ than requested; the route is not admitted to that destination sheet".into()));
            }
            Ok((
                result.solution,
                BranchTrace {
                    compiled,
                    segments: result.branch_segments,
                },
            ))
        },
        |old, old_trace, new, new_trace, index| {
            let state = &new_trace.segments[index].end;
            let germ = classify_germ(&new_trace.compiled, state)?;
            let point = &new.segments[index].end;
            let mut found = false;
            for (old_index, segment) in old_trace.segments.iter().enumerate() {
                if old.evaluate_segment(old_index, point).is_err() {
                    continue;
                }
                let branch = match old_trace.compiled.project_branch_segment(segment, point) {
                    Ok(branch) => branch,
                    Err(Error::Accuracy(_)) => return Ok(false),
                    Err(error) => return Err(error),
                };
                let old_germ = match classify_germ(&old_trace.compiled, &branch) {
                    Ok(germ) => germ,
                    Err(Error::Accuracy(_)) => return Ok(false),
                    Err(error) => return Err(error),
                };
                if old_germ != germ {
                    return Ok(false);
                }
                found = true;
            }
            Ok(found)
        },
    )?;
    let mut germs = BTreeMap::new();
    if !system.source().roots.is_empty() {
        for checkpoint in &solution.checkpoints {
            germs.insert(
                checkpoint.segment,
                classify_germ(&trace.compiled, &trace.segments[checkpoint.segment].end)?,
            );
        }
    }
    Ok((solution, germs))
}
