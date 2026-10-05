//! Physical boundary reuse for supplied differential systems.
//!
//! Exact coordinates and numerically chosen coordinates remain distinct. A
//! precision label is evidence supplied by the caller, never inferred from the
//! number of MPFR working bits. Path admissibility is required during selection:
//! distance alone cannot decide which continuation reaches the intended sheet.
use crate::algebraic::{AlgebraicKinematicSystem, CanonicalAlgebraicSystem, RootSeed, SquareRoot};
use crate::diffexp::{EpsilonBoundary, EpsilonSolution};
use crate::kinematics::{KinematicPath, KinematicSystem};
use crate::{ComplexFloat as C, Error, Precision, Prescription, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;
use std::sync::{
    Arc,
    atomic::{AtomicU64, Ordering},
};
use symbolica::coefficient::Coefficient;
use symbolica::prelude::*;

const VERSION: u32 = 5;
const MAGIC: &[u8] = b"AMFLOW-BOUNDARIES\0\x05";
const FILE: &str = "physical-boundaries.bin";
static NEXT_FILE: AtomicU64 = AtomicU64::new(0);

fn fingerprint<T: Serialize>(value: &T) -> Result<String> {
    Ok(
        blake3::hash(&serde_json::to_vec(value).map_err(|e| Error::Cache(e.to_string()))?)
            .to_hex()
            .to_string(),
    )
}

/// Inclusive epsilon powers. Reuse preserves the leading power because lower
/// coefficients can feed higher ones through the differential equations.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct EpsilonRange {
    pub leading: i32,
    pub last: i32,
}
impl EpsilonRange {
    pub fn new(leading: i32, last: i32) -> Result<Self> {
        let count = i64::from(last) - i64::from(leading) + 1;
        if !(1..=1024).contains(&count) {
            return Err(Error::InvalidInput(
                "boundary cache requires 1 through 1024 epsilon coefficients".into(),
            ));
        }
        Ok(Self { leading, last })
    }
    fn count(self) -> Result<usize> {
        Self::new(self.leading, self.last)?;
        Ok((i64::from(self.last) - i64::from(self.leading) + 1) as usize)
    }
    fn covers(self, requested: Self) -> bool {
        self.leading == requested.leading && self.last >= requested.last
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum PointKind {
    Physical,
    /// An explicitly retained numerical continuation detour; excluded by default.
    ContourDetour,
}

/// A local sign relative to the principal square root at the exact point.
/// This is discrete sheet information, not a precision claim about a root value.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum RootSheet {
    Principal,
    Opposite,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RootGerm {
    pub sheets: BTreeMap<Symbol, RootSheet>,
}
impl RootGerm {
    pub(crate) fn seeds(&self) -> BTreeMap<Symbol, RootSeed> {
        self.sheets
            .iter()
            .map(|(&s, sheet)| {
                (
                    s,
                    match sheet {
                        RootSheet::Principal => RootSeed::Principal,
                        RootSheet::Opposite => RootSeed::Opposite,
                    },
                )
            })
            .collect()
    }
    fn canonical(&self) -> Vec<(String, RootSheet)> {
        let mut out = self
            .sheets
            .iter()
            .map(|(&s, &sheet)| (Atom::var(s).to_canonical_string(), sheet))
            .collect::<Vec<_>>();
        out.sort_by(|a, b| a.0.cmp(&b.0));
        out
    }
}

/// Exact input, an exact path image of a numerical parameter, or independently
/// rounded numerical coordinates. The origin of a point remains explicit.
#[derive(Clone, Debug)]
pub enum CachedPoint {
    /// A point on the registered algebraic cover. The nested coordinate origin
    /// stays Exact, Derived, or Numerical; nested root wrappers are rejected.
    Algebraic {
        point: Box<CachedPoint>,
        germ: RootGerm,
    },
    Exact(BTreeMap<Symbol, Atom>),
    /// Exact coordinate image of an accepted, rounded parameter value. Keeping
    /// the image exact preserves multivariate path constraints on restart.
    Derived {
        coordinates: BTreeMap<Symbol, Atom>,
        working_bits: u32,
        provenance: String,
    },
    Numerical {
        coordinates: BTreeMap<Symbol, C>,
        working_bits: u32,
        provenance: String,
    },
}

fn exact_constant(a: AtomView<'_>) -> Result<()> {
    match a {
        AtomView::Num(n) if matches!(n.get_coeff_view().to_owned(), Coefficient::Complex(_)) => {
            Ok(())
        }
        AtomView::Add(v) => v.iter().try_for_each(exact_constant),
        AtomView::Mul(v) => v.iter().try_for_each(exact_constant),
        AtomView::Pow(v) => {
            let (base, exponent) = v.get_base_exp();
            exact_constant(base)?;
            if let AtomView::Num(n) = exponent
                && let Coefficient::Complex(c) = n.get_coeff_view().to_owned()
                && c.im.is_zero()
            {
                return Ok(());
            }
            Err(Error::InvalidInput(
                "exact point requires rational algebraic exponents".into(),
            ))
        }
        _ => Err(Error::InvalidInput(
            "exact point requires finite exact constants, not floats or free symbols".into(),
        )),
    }
}

pub(crate) fn exact_real_rational(a: &Atom) -> Result<()> {
    let reduced = a.together().cancel();
    if let AtomView::Num(n) = reduced.as_view()
        && let Coefficient::Complex(c) = n.get_coeff_view().to_owned()
        && c.im.is_zero()
    {
        return Ok(());
    }
    Err(Error::Unsupported(
        "prescribed real-contour planning requires exact real rational coordinates and radicands"
            .into(),
    ))
}

// The registered quotient-domain owner retains raw denominator occurrences.
fn algebraic_conditions(system: &AlgebraicKinematicSystem) -> Result<Vec<Atom>> {
    system.nonzero_conditions()
}

fn exact_metadata(a: AtomView<'_>) -> Result<()> {
    match a {
        AtomView::Num(n) if !matches!(n.get_coeff_view().to_owned(), Coefficient::Complex(_)) => {
            Err(Error::InvalidInput(
                "basis and normalization metadata must be exact".into(),
            ))
        }
        AtomView::Add(v) => v.iter().try_for_each(exact_metadata),
        AtomView::Mul(v) => v.iter().try_for_each(exact_metadata),
        AtomView::Pow(v) => v.iter().try_for_each(exact_metadata),
        AtomView::Fun(v) => v.iter().try_for_each(exact_metadata),
        _ => Ok(()),
    }
}

impl CachedPoint {
    pub fn with_root_germ(self, germ: RootGerm) -> Result<Self> {
        let result = Self::Algebraic {
            point: Box::new(self),
            germ,
        };
        result.validate()?;
        Ok(result)
    }
    pub fn root_germ(&self) -> Option<&RootGerm> {
        match self {
            Self::Algebraic { germ, .. } => Some(germ),
            _ => None,
        }
    }
    pub fn validate(&self) -> Result<()> {
        match self {
            Self::Algebraic { point, germ } => {
                if germ.sheets.is_empty() || point.root_germ().is_some() {
                    return Err(Error::InvalidInput(
                        "root germ must be nonempty and occur once per point".into(),
                    ));
                }
                point.validate()?;
            }

            Self::Exact(coordinates) | Self::Derived { coordinates, .. } => {
                if coordinates.is_empty() {
                    return Err(Error::InvalidInput("empty cached point".into()));
                }
                for value in coordinates.values() {
                    exact_constant(value.as_view())?;
                }
            }
            Self::Numerical {
                coordinates,
                working_bits,
                provenance,
            } => {
                if coordinates.is_empty() || *working_bits < 2 || provenance.trim().is_empty() {
                    return Err(Error::InvalidInput(
                        "numerical point requires coordinates, precision, and provenance".into(),
                    ));
                }
                for value in coordinates.values() {
                    if !value.re.is_finite()
                        || !value.im.is_finite()
                        || value.re.as_raw().prec() < *working_bits
                        || value.im.as_raw().prec() < *working_bits
                    {
                        return Err(Error::InvalidInput("numerical coordinate is nonfinite or lacks its declared working precision".into()));
                    }
                }
            }
        }
        if let Self::Derived {
            working_bits,
            provenance,
            ..
        } = self
            && (*working_bits < 2 || provenance.trim().is_empty())
        {
            return Err(Error::InvalidInput(
                "derived point requires parameter precision and provenance".into(),
            ));
        }
        Ok(())
    }

    /// Return retained exact coordinates for exact/derived points, or represent
    /// independently rounded machine coordinates as exact binary rationals. This
    /// does not undo rounding, improve accuracy, or change the point's origin.
    pub fn rounded_coordinates_as_exact(&self) -> Result<BTreeMap<Symbol, Atom>> {
        self.validate()?;
        match self {
            Self::Algebraic { point, .. } => point.rounded_coordinates_as_exact(),
            Self::Exact(coordinates) | Self::Derived { coordinates, .. } => Ok(coordinates.clone()),
            Self::Numerical { coordinates, .. } => Ok(coordinates
                .iter()
                .map(|(&symbol, value)| {
                    (
                        symbol,
                        Atom::num(symbolica::domains::float::Complex::new(
                            value.re.to_rational(),
                            value.im.to_rational(),
                        )),
                    )
                })
                .collect()),
        }
    }

    /// Coordinates for an exact restart chart, retaining a derived point's
    /// original path image rather than rounding its components independently.
    pub fn restart_coordinates(&self) -> Result<BTreeMap<Symbol, Atom>> {
        self.rounded_coordinates_as_exact()
    }

    pub fn evaluate(&self, p: Precision) -> Result<BTreeMap<Symbol, C>> {
        self.validate()?;
        match self {
            Self::Algebraic { point, .. } => point.evaluate(p),
            Self::Exact(coordinates) | Self::Derived { coordinates, .. } => coordinates
                .iter()
                .map(|(&v, a)| Ok((v, p.eval(a, &Default::default())?)))
                .collect(),
            Self::Numerical { coordinates, .. } => {
                Ok(coordinates.iter().map(|(&v, a)| (v, p.round(a))).collect())
            }
        }
    }
    fn coordinate_bits(&self) -> u32 {
        match self {
            Self::Algebraic { point, .. } => point.coordinate_bits(),
            Self::Exact(_) => u32::MAX,
            Self::Numerical { working_bits, .. } | Self::Derived { working_bits, .. } => {
                *working_bits
            }
        }
    }
    pub(crate) fn key(&self) -> Result<String> {
        match self {
            Self::Algebraic { point, germ } => {
                fingerprint(&("algebraic", point.key()?, germ.canonical()))
            }
            Self::Exact(coordinates) => {
                let mut pairs = coordinates
                    .iter()
                    .map(|(&v, a)| (Atom::var(v).to_canonical_string(), a.to_canonical_string()))
                    .collect::<Vec<_>>();
                pairs.sort();
                fingerprint(&("exact", pairs))
            }
            Self::Derived {
                coordinates,
                working_bits,
                provenance,
            } => {
                let mut pairs = coordinates
                    .iter()
                    .map(|(&v, a)| (Atom::var(v).to_canonical_string(), a.to_canonical_string()))
                    .collect::<Vec<_>>();
                pairs.sort();
                fingerprint(&("derived", working_bits, provenance, pairs))
            }
            Self::Numerical {
                coordinates,
                working_bits,
                provenance,
            } => {
                let mut pairs = coordinates
                    .iter()
                    .map(|(&v, a)| (Atom::var(v).to_canonical_string(), a))
                    .collect::<Vec<_>>();
                pairs.sort_by(|a, b| a.0.cmp(&b.0));
                fingerprint(&("numerical", working_bits, provenance, pairs))
            }
        }
    }
}

#[derive(Clone, Debug)]
enum IdentitySystem {
    Dense(KinematicSystem),
    Canonical(CanonicalAlgebraicSystem),
}
impl IdentitySystem {
    fn epsilon(&self) -> Symbol {
        match self {
            Self::Dense(system) => system.epsilon,
            Self::Canonical(system) => system.epsilon(),
        }
    }
    fn variables(&self) -> std::collections::BTreeSet<Symbol> {
        match self {
            Self::Dense(system) => system.derivatives.keys().copied().collect(),
            Self::Canonical(system) => system.variables().iter().copied().collect(),
        }
    }
    fn dimension(&self) -> usize {
        match self {
            Self::Dense(system) => system.derivatives.values().next().unwrap().len(),
            Self::Canonical(system) => system.dimension(),
        }
    }
    fn canonical_key(&self) -> Result<String> {
        let expression = AtomCore::to_canonical_string;
        Ok(match self {
            Self::Dense(system) => {
                let mut derivatives = system
                    .derivatives
                    .iter()
                    .map(|(&s, m)| {
                        (
                            Atom::var(s).to_canonical_string(),
                            m.iter()
                                .map(|row| row.iter().map(expression).collect::<Vec<_>>())
                                .collect::<Vec<_>>(),
                        )
                    })
                    .collect::<Vec<_>>();
                derivatives.sort_by(|a, b| a.0.cmp(&b.0));
                fingerprint(&("dense", derivatives))?
            }
            Self::Canonical(system) => fingerprint(&(
                "canonical-dlog",
                system.dimension(),
                system
                    .variables()
                    .iter()
                    .map(|&v| Atom::var(v).to_canonical_string())
                    .collect::<Vec<_>>(),
                system.letters().iter().map(expression).collect::<Vec<_>>(),
                system
                    .constant_matrices()
                    .iter()
                    .map(|matrix| {
                        matrix
                            .iter()
                            .map(|row| row.iter().map(expression).collect::<Vec<_>>())
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>(),
            ))?,
        })
    }
}

/// Exact convention for planner-generated affine physical continuations.
/// Geometry verifies the specified local sides; a caller must additionally
/// admit each route into this integral's global branch/homotopy domain.
#[derive(Clone, Debug)]
pub struct PhysicalContinuation {
    pub prescriptions: Vec<crate::contour::PolynomialPrescription>,
    pub unprescribed_side: Prescription,
    pub domain: String,
}
impl PhysicalContinuation {
    fn validate(&self, variables: &std::collections::BTreeSet<Symbol>) -> Result<()> {
        if self.domain.trim().is_empty() {
            return Err(Error::InvalidInput(
                "physical continuation needs an explicit homotopy domain".into(),
            ));
        }
        let allowed = variables.iter().copied().map(Atom::var).collect();
        for declaration in &self.prescriptions {
            let mut used = std::collections::BTreeSet::new();
            crate::family::scalar_symbols(declaration.polynomial.as_view(), &mut used)?;
            if !used.is_subset(&allowed) {
                return Err(Error::InvalidInput(
                    "physical prescription contains epsilon, roots, or undeclared variables".into(),
                ));
            }
            let fraction: RationalPolynomial<IntegerRing, u16> = declaration
                .polynomial
                .try_to_rational_polynomial(&Q, &Z, None)
                .map_err(|e| {
                    Error::Unsupported(format!(
                        "physical prescription must be an exact real polynomial: {e}"
                    ))
                })?;
            if fraction.numerator.is_zero()
                || !fraction.denominator.is_constant()
                || fraction
                    .numerator
                    .variables()
                    .iter()
                    .any(|v| !matches!(v, PolyVariable::Symbol(s) if variables.contains(s)))
            {
                return Err(Error::Unsupported(
                    "physical prescription must be a nonzero exact real polynomial".into(),
                ));
            }
            // Retain the raw expression's domain: a removable denominator is
            // not silently reinterpreted as a polynomial declaration.
            if !crate::physical_conditions::rational_denominator_conditions(
                std::slice::from_ref(&declaration.polynomial),
                variables,
            )?
            .is_empty()
            {
                return Err(Error::Unsupported(
                    "physical prescription contains an original nonconstant denominator".into(),
                ));
            }
        }
        Ok(())
    }
    fn canonical(&self) -> (String, bool, Vec<(String, bool)>) {
        (
            self.domain.clone(),
            matches!(self.unprescribed_side, Prescription::PlusI0),
            self.prescriptions
                .iter()
                .map(|p| {
                    (
                        p.polynomial.to_canonical_string(),
                        matches!(p.prescription, Prescription::PlusI0),
                    )
                })
                .collect(),
        )
    }
}

#[derive(Clone, Debug)]
struct IdentityData {
    key: String,
    system: IdentitySystem,
    variables: std::collections::BTreeSet<Symbol>,
    roots: Vec<SquareRoot>,
    basis: Vec<Atom>,
    normalization: Atom,
    prescription: Prescription,
    domain: String,
    conditions: Vec<Atom>,
    continuation: Option<PhysicalContinuation>,
}

/// Content identity for the ordered basis and its physical branch. `domain`
/// identifies a caller-defined homotopy/sheet domain, not just a kinematic region.
/// It must distinguish any additional convention that changes the values.
#[derive(Clone, Debug)]
pub struct BoundaryIdentity(Arc<IdentityData>);
impl BoundaryIdentity {
    pub fn new(
        system: &KinematicSystem,
        basis: &[Atom],
        normalization: &Atom,
        prescription: Prescription,
        domain: &str,
    ) -> Result<Self> {
        Self::with_conditions(system, basis, normalization, prescription, domain, &[])
    }

    /// Attach the exact nonzero assumptions used to derive the system. They
    /// remain part of cache identity even if matrix cancellation removes them.
    pub fn with_conditions(
        system: &KinematicSystem,
        basis: &[Atom],
        normalization: &Atom,
        prescription: Prescription,
        domain: &str,
        conditions: &[Atom],
    ) -> Result<Self> {
        Self::build(
            system,
            &[],
            basis,
            normalization,
            prescription,
            domain,
            conditions,
        )
    }

    pub fn with_algebraic_conditions(
        system: &AlgebraicKinematicSystem,
        basis: &[Atom],
        normalization: &Atom,
        prescription: Prescription,
        domain: &str,
        conditions: &[Atom],
    ) -> Result<Self> {
        system.validate()?;
        Self::build(
            &KinematicSystem {
                epsilon: system.epsilon,
                derivatives: system.derivatives.clone(),
            },
            &system.roots,
            basis,
            normalization,
            prescription,
            domain,
            conditions,
        )
    }

    /// Identity for an unassembled canonical connection. The ordered letters
    /// and constant matrices are distinct from a dense connection identity.
    /// Source guards were computed by the canonical constructor and are reused.
    pub fn with_canonical_conditions(
        system: &CanonicalAlgebraicSystem,
        basis: &[Atom],
        normalization: &Atom,
        prescription: Prescription,
        domain: &str,
        conditions: &[Atom],
    ) -> Result<Self> {
        let mut candidates = conditions.to_vec();
        candidates.extend_from_slice(system.nonzero_conditions());
        Self::from_parts(
            IdentitySystem::Canonical(system.clone()),
            system.roots(),
            basis,
            normalization,
            prescription,
            domain,
            &candidates,
        )
    }

    #[allow(clippy::too_many_arguments)] // Existing identity metadata plus its optional exact root registry.
    fn build(
        system: &KinematicSystem,
        roots: &[SquareRoot],
        basis: &[Atom],
        normalization: &Atom,
        prescription: Prescription,
        domain: &str,
        conditions: &[Atom],
    ) -> Result<Self> {
        if roots.is_empty() {
            system.validate()?;
        }
        let mut candidates = conditions.to_vec();
        if roots.is_empty() {
            candidates.extend(crate::physical_conditions::matrix_domain_conditions(
                system,
            )?);
        } else {
            candidates.extend(algebraic_conditions(&AlgebraicKinematicSystem {
                epsilon: system.epsilon,
                derivatives: system.derivatives.clone(),
                roots: roots.to_vec(),
            })?);
        }
        Self::from_parts(
            IdentitySystem::Dense(system.clone()),
            roots,
            basis,
            normalization,
            prescription,
            domain,
            &candidates,
        )
    }

    #[allow(clippy::too_many_arguments)] // Shared metadata for dense and canonical exact representations.
    fn from_parts(
        system: IdentitySystem,
        roots: &[SquareRoot],
        basis: &[Atom],
        normalization: &Atom,
        prescription: Prescription,
        domain: &str,
        candidates: &[Atom],
    ) -> Result<Self> {
        let physical_variables = system.variables();
        let mut variables = physical_variables.clone();
        variables.insert(system.epsilon());
        let mut conditions =
            crate::physical_conditions::canonical_conditions(candidates, &variables)?;
        if !roots.is_empty() {
            let leading = conditions
                .iter()
                .map(|condition| {
                    crate::physical_conditions::epsilon_leading_coefficient(
                        condition,
                        system.epsilon(),
                    )
                })
                .collect::<Result<Vec<_>>>()?;
            conditions.extend(leading);
            conditions = crate::physical_conditions::canonical_conditions(&conditions, &variables)?;
        }
        let n = system.dimension();
        if basis.len() != n || domain.trim().is_empty() || normalization.is_zero() {
            return Err(Error::InvalidInput("boundary identity needs an ordered basis, nonzero normalization, and homotopy domain".into()));
        }
        for a in basis.iter().chain(std::iter::once(normalization)) {
            exact_metadata(a.as_view())?;
        }
        let root_identity = roots
            .iter()
            .map(|r| {
                (
                    Atom::var(r.symbol).to_canonical_string(),
                    r.radicand.to_canonical_string(),
                )
            })
            .collect::<Vec<_>>();
        let key = fingerprint(&(
            VERSION,
            root_identity,
            env!("DEPENDENCY_SOURCE_DIGEST"),
            env!("PORT_SOURCE_DIGEST"),
            Atom::var(system.epsilon()).to_canonical_string(),
            system.canonical_key()?,
            basis
                .iter()
                .map(AtomCore::to_canonical_string)
                .collect::<Vec<_>>(),
            normalization.to_canonical_string(),
            matches!(prescription, Prescription::PlusI0),
            domain,
            conditions
                .iter()
                .map(AtomCore::to_canonical_string)
                .collect::<Vec<_>>(),
        ))?;
        Ok(Self(Arc::new(IdentityData {
            key,
            system,
            variables: physical_variables,
            roots: roots.to_vec(),
            basis: basis.to_vec(),
            normalization: normalization.clone(),
            prescription,
            domain: domain.into(),
            conditions,
            continuation: None,
        })))
    }
    /// Add an explicit routed representation; an ordinary identity cannot
    /// alias it even when the user-supplied branch-domain labels are equal.
    pub fn with_prescribed_continuation(&self, continuation: PhysicalContinuation) -> Result<Self> {
        if self.0.continuation.is_some() {
            return Err(Error::InvalidInput(
                "physical continuation is already bound to this identity".into(),
            ));
        }
        continuation.validate(&self.0.variables)?;
        let mut data = (*self.0).clone();
        data.key = fingerprint(&("prescribed-affine-v1", &data.key, continuation.canonical()))?;
        data.continuation = Some(continuation);
        Ok(Self(Arc::new(data)))
    }
    pub fn physical_continuation(&self) -> Option<&PhysicalContinuation> {
        self.0.continuation.as_ref()
    }
    pub(crate) fn prescribed_path_conditions(&self, path: &KinematicPath) -> Result<Vec<Atom>> {
        let rules = path
            .coordinates
            .iter()
            .map(|(&s, a)| (Atom::var(s), a.clone()))
            .collect();
        self.0
            .conditions
            .iter()
            .map(|condition| {
                let generic = crate::physical_conditions::epsilon_leading_coefficient(
                    condition,
                    self.0.system.epsilon(),
                )?;
                let restricted = crate::family::substitute(&generic, &rules)
                    .together()
                    .cancel();
                if restricted.is_zero() {
                    return Err(Error::Unsupported(
                        "physical path lies on an original nonzero condition".into(),
                    ));
                }
                Ok(restricted)
            })
            .collect()
    }
    pub fn key(&self) -> &str {
        &self.0.key
    }
    pub fn dimension(&self) -> usize {
        self.0.basis.len()
    }

    pub fn nonzero_conditions(&self) -> &[Atom] {
        &self.0.conditions
    }

    /// Check reduction assumptions and registered-root sheets along an exact affine path. The
    /// caller still owns its physical branch/homotopy admissibility policy.
    pub fn conditions_admit_straight_path(
        &self,
        source: &CachedPoint,
        target: &CachedPoint,
        p: Precision,
        digits: u32,
    ) -> Result<bool> {
        self.conditions_admit_straight_path_with_context(
            source,
            target,
            p,
            digits,
            &crate::RunContext::default(),
        )
    }

    pub(crate) fn conditions_admit_straight_path_with_context(
        &self,
        source: &CachedPoint,
        target: &CachedPoint,
        p: Precision,
        digits: u32,
        context: &crate::RunContext,
    ) -> Result<bool> {
        context.cancellation.check()?;
        self.validate_point(source)?;
        self.validate_point(target)?;
        if self.0.conditions.is_empty() && self.0.roots.is_empty() {
            return Ok(true);
        }
        let parameter = self.path_parameter("guard_path")?;
        let path = KinematicPath::straight_line(
            parameter,
            &source.rounded_coordinates_as_exact()?,
            &target.rounded_coordinates_as_exact()?,
        )?;
        if !self.0.roots.is_empty() {
            let rules = path
                .coordinates
                .iter()
                .map(|(&s, a)| (Atom::var(s), a.clone()))
                .collect();
            let source_germ = source.root_germ().unwrap();
            let target_germ = target.root_germ().unwrap();
            for root in &self.0.roots {
                let radicand = crate::family::substitute(&root.radicand, &rules);
                let Some(flip) =
                    crate::root_path::principal_flip(&radicand, path.parameter, p, context)?
                else {
                    return Ok(false);
                };
                if (source_germ.sheets[&root.symbol] != target_germ.sheets[&root.symbol]) != flip {
                    return Ok(false);
                }
            }
        }
        crate::physical_conditions::conditions_admit_path(
            &self.0.conditions,
            self.0.system.epsilon(),
            &path,
            p,
            digits,
        )
    }

    pub(crate) fn path_parameter(&self, stem: &str) -> Result<Symbol> {
        let mut index = 0usize;
        loop {
            let candidate = symbol!(&format!("symbolica_amflow::{stem}_{index}"));
            if candidate != self.0.system.epsilon()
                && !self.0.variables.contains(&candidate)
                && self.0.roots.iter().all(|r| r.symbol != candidate)
            {
                return Ok(candidate);
            }
            index = index
                .checked_add(1)
                .ok_or_else(|| Error::Limit("path symbol overflow".into()))?;
        }
    }

    pub fn roots(&self) -> &[SquareRoot] {
        &self.0.roots
    }

    fn validate_root_point(&self, point: &CachedPoint) -> Result<()> {
        if self.0.roots.is_empty() {
            if point.root_germ().is_some() {
                return Err(Error::InvalidInput(
                    "rational system cannot use a registered-root point".into(),
                ));
            }
            return Ok(());
        }
        let germ = point.root_germ().ok_or_else(|| {
            Error::InvalidInput("algebraic boundary requires an explicit root germ".into())
        })?;
        if germ.sheets.len() != self.0.roots.len()
            || self
                .0
                .roots
                .iter()
                .any(|r| !germ.sheets.contains_key(&r.symbol))
        {
            return Err(Error::InvalidInput(
                "cached root germ does not match the exact registry".into(),
            ));
        }
        let coordinates = point.restart_coordinates()?;
        for value in coordinates.values() {
            crate::root_path::exact_complex_rational(value)?;
        }
        let rules = coordinates
            .into_iter()
            .map(|(s, a)| (Atom::var(s), a))
            .collect();
        for root in &self.0.roots {
            let value = crate::family::substitute(&root.radicand, &rules)
                .together()
                .cancel();
            crate::root_path::exact_complex_rational(&value)?;
            if value.is_zero() {
                return Err(Error::InvalidInput(
                    "cached algebraic point is a root branch point".into(),
                ));
            }
        }
        Ok(())
    }

    pub(crate) fn validate_point(&self, point: &CachedPoint) -> Result<()> {
        point.validate()?;
        self.validate_root_point(point)?;
        let coordinates = point.rounded_coordinates_as_exact()?;
        if coordinates
            .keys()
            .copied()
            .collect::<std::collections::BTreeSet<_>>()
            != self.0.variables
        {
            return Err(Error::InvalidInput(
                "cached point coordinates do not match the system".into(),
            ));
        }
        crate::physical_conditions::validate_conditions_at(
            &self.0.conditions,
            self.0.system.epsilon(),
            &coordinates,
        )?;
        let rules = coordinates
            .iter()
            .map(|(&v, a)| (Atom::var(v), a.clone()))
            .collect();
        fn regular(a: AtomView<'_>, rules: &BTreeMap<Atom, Atom>) -> Result<()> {
            match a {
                AtomView::Add(v) => v.iter().try_for_each(|a| regular(a, rules)),
                AtomView::Mul(v) => v.iter().try_for_each(|a| regular(a, rules)),
                AtomView::Pow(v) => {
                    let (base, exponent) = v.get_base_exp();
                    if let AtomView::Num(n) = exponent
                        && let Coefficient::Complex(c) = n.get_coeff_view().to_owned()
                        && (c.re < 0 || !c.re.is_integer())
                        && crate::family::substitute(&base.to_owned(), rules)
                            .together()
                            .cancel()
                            .is_zero()
                    {
                        return Err(Error::InvalidInput(
                            "cached point is a matrix pole or algebraic branch point".into(),
                        ));
                    }
                    regular(base, rules)
                }
                _ => Ok(()),
            }
        }
        if let IdentitySystem::Dense(system) = &self.0.system {
            for a in system.derivatives.values().flatten().flatten() {
                regular(a.as_view(), &rules)?;
            }
        }
        Ok(())
    }
}

/// Caller-recorded accuracy evidence. Comparison errors are consistency estimates,
/// not rigorous error bounds. The provenance must explain how accuracy was obtained.
/// A verified digit count d uses the mixed coefficient scale
/// 10^(-d) * max(1, |c|), not d significant digits for arbitrarily small c.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BoundaryAccuracy {
    verified_digits: u32,
    working_bits: u32,
    input_verified_digits: u32,
    comparison_errors: Vec<Vec<Float>>,
    provenance: String,
}
impl BoundaryAccuracy {
    /// Record independently established input accuracy (e.g. an analytic value or
    /// upstream reference with precision metadata). This does not certify the claim.
    pub fn supplied(
        verified_digits: u32,
        working_bits: u32,
        comparison_errors: Vec<Vec<Float>>,
        provenance: &str,
    ) -> Result<Self> {
        let result = Self {
            verified_digits,
            working_bits,
            input_verified_digits: verified_digits,
            comparison_errors,
            provenance: provenance.into(),
        };
        result.validate()?;
        Ok(result)
    }
    /// Retain the weaker of the transport check and its initial boundary evidence.
    /// Repeating transport from the same low-accuracy boundary cannot improve this cap.
    pub fn from_solution(
        solution: &EpsilonSolution,
        input_verified_digits: u32,
        provenance: &str,
    ) -> Result<Self> {
        let checked = solution.verified_digits.ok_or_else(|| {
            Error::Accuracy("unverified transport is not reusable boundary evidence".into())
        })?;
        let result = Self {
            verified_digits: checked.min(input_verified_digits),
            working_bits: solution.diagnostics.working_bits,
            input_verified_digits,
            comparison_errors: solution.comparison_errors.clone(),
            provenance: provenance.into(),
        };
        result.validate()?;
        Ok(result)
    }
    pub fn verified_digits(&self) -> u32 {
        self.verified_digits
    }
    pub fn working_bits(&self) -> u32 {
        self.working_bits
    }
    pub fn comparison_errors(&self) -> &[Vec<Float>] {
        &self.comparison_errors
    }
    fn validate(&self) -> Result<()> {
        if self.verified_digits == 0
            || self.verified_digits > self.input_verified_digits
            || Precision::decimal(self.verified_digits)?.bits > self.working_bits
            || self.provenance.trim().is_empty()
            || self
                .comparison_errors
                .iter()
                .flatten()
                .any(|e| !e.is_finite() || *e < Float::with_val(32, 0))
        {
            return Err(Error::InvalidInput(
                "invalid boundary accuracy evidence".into(),
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub struct CachedBoundary {
    pub identity: BoundaryIdentity,
    pub point: CachedPoint,
    pub kind: PointKind,
    pub range: EpsilonRange,
    pub coefficients: Vec<Vec<C>>,
    pub accuracy: BoundaryAccuracy,
}
impl CachedBoundary {
    pub fn validate(&self) -> Result<()> {
        self.identity.validate_point(&self.point)?;
        self.accuracy.validate()?;
        let n = self.identity.dimension();
        let count = self.range.count()?;
        if self.coefficients.len() != count
            || self.coefficients.iter().any(|r| r.len() != n)
            || self.accuracy.comparison_errors.len() != count
            || self.accuracy.comparison_errors.iter().any(|r| r.len() != n)
        {
            return Err(Error::InvalidInput(
                "cached epsilon coefficients/evidence dimensions".into(),
            ));
        }
        let p = Precision {
            bits: self.accuracy.working_bits,
        };
        for (value, error) in self
            .coefficients
            .iter()
            .flatten()
            .zip(self.accuracy.comparison_errors.iter().flatten())
        {
            if !p.finite(value)
                || value.re.as_raw().prec() < p.bits
                || value.im.as_raw().prec() < p.bits
            {
                return Err(Error::InvalidInput(
                    "cached coefficient lacks declared working precision".into(),
                ));
            }
            let norm = p.norm(value);
            let scale = if norm > p.real(1) { norm } else { p.real(1) };
            if *error > p.tolerance(self.accuracy.verified_digits) * scale {
                return Err(Error::Accuracy(
                    "cached comparison error exceeds claimed accuracy".into(),
                ));
            }
        }
        Ok(())
    }

    /// Rebase the stored physical value to the new path's parameter coordinate.
    /// Coefficients are rounded, never upgraded in accuracy, at `precision`.
    pub fn as_epsilon_boundary(
        &self,
        parameter_point: &Atom,
        range: EpsilonRange,
        precision: Precision,
    ) -> Result<EpsilonBoundary> {
        self.validate()?;
        range.count()?;
        if precision.bits < 2 {
            return Err(Error::InvalidInput(
                "invalid boundary evaluation precision".into(),
            ));
        }
        if !self.range.covers(range) {
            return Err(Error::InvalidInput(
                "cached epsilon range does not cover requested leading/order".into(),
            ));
        }
        exact_constant(parameter_point.as_view())?;
        Ok(EpsilonBoundary {
            point: precision.eval(parameter_point, &Default::default())?,
            leading: range.leading,
            coefficients: self.coefficients[..range.count()?]
                .iter()
                .map(|r| r.iter().map(|v| precision.round(v)).collect())
                .collect(),
        })
    }
    fn key(&self) -> Result<String> {
        fingerprint(&(
            self.identity.key(),
            self.point.key()?,
            self.kind,
            self.range,
        ))
    }
}

/// Compatibility is checked before asking a policy to price a route.
pub trait TransportCost {
    /// `None` rejects an unsafe route, incompatible sheet, or other policy exclusion.
    fn cost(
        &self,
        source: &CachedBoundary,
        target: &CachedPoint,
        precision: Precision,
    ) -> Result<Option<Float>>;

    /// Optional cheap lower bound on the nonnegative, rounded value returned by
    /// `cost` at this precision. `None` means unknown, not an excluded route.
    /// Every returned bound must be finite, nonnegative, and no greater than
    /// `cost` when that route is accepted. A bound may also describe a route
    /// that `cost` rejects. Selection can skip `cost` when this bound exceeds
    /// an accepted cost; it therefore relies on this contract for skipped routes.
    ///
    /// Bounds must not depend on evaluation order or mutable policy state.
    /// Implementations that require cancellation checks while pricing a large
    /// bank should perform those checks here too. Errors propagate immediately.
    /// The default preserves exhaustive evaluation for existing policies.
    fn lower_bound(
        &self,
        _source: &CachedBoundary,
        _target: &CachedPoint,
        _precision: Precision,
    ) -> Result<Option<Float>> {
        Ok(None)
    }

    /// Resolve equal rounded costs, after both routes passed `cost`. The
    /// default prefers an exact coordinate hit; otherwise accuracy breaks the
    /// tie. Distance policies can distinguish costs beyond the initial working
    /// precision without imposing that ordering on unrelated cost policies.
    fn compare_tied_costs(
        &self,
        left: &CachedBoundary,
        right: &CachedBoundary,
        target: &CachedPoint,
        _precision: Precision,
    ) -> Result<std::cmp::Ordering> {
        compare_exact_hits(&left.point, &right.point, target)
    }
}

fn compare_exact_hits(
    left: &CachedPoint,
    right: &CachedPoint,
    target: &CachedPoint,
) -> Result<std::cmp::Ordering> {
    let target_germ = target.root_germ();
    let target = target.restart_coordinates()?;
    let is_hit = |point: &CachedPoint| -> Result<bool> {
        if point.root_germ() != target_germ {
            return Ok(false);
        }
        let point = point.restart_coordinates()?;
        Ok(point.keys().eq(target.keys())
            && point
                .iter()
                .all(|(s, value)| (value - &target[s]).together().cancel().is_zero()))
    };
    Ok(is_hit(right)?.cmp(&is_hit(left)?))
}

/// Squared Euclidean distance in scaled kinematic coordinates, computed in MPFR.
/// The mandatory callback checks path/sheet/pole admissibility before pricing.
pub struct ScaledDistance<F> {
    pub scales: BTreeMap<Symbol, Atom>,
    pub admissible: F,
}
impl<F> ScaledDistance<F> {
    fn differences(&self, source: &CachedPoint, target: &CachedPoint) -> Result<Vec<Atom>> {
        let source = source.restart_coordinates()?;
        let target = target.restart_coordinates()?;
        if !source.keys().eq(target.keys()) || self.scales.keys().any(|s| !source.contains_key(s)) {
            return Err(Error::InvalidInput(
                "distance coordinates/scales do not match".into(),
            ));
        }
        source
            .into_iter()
            .map(|(variable, value)| {
                let scale = self.scales.get(&variable).cloned().unwrap_or(Atom::num(1));
                exact_constant(scale.as_view())?;
                Ok(((value - &target[&variable]) / scale).together().cancel())
            })
            .collect()
    }

    fn checked_distance(
        &self,
        source: &CachedPoint,
        target: &CachedPoint,
        p: Precision,
    ) -> Result<Float> {
        for scale in self.scales.values() {
            exact_constant(scale.as_view())?;
            let scale = p.eval(scale, &Default::default())?;
            if !p.finite(&scale) || scale.im != p.real(0) || scale.re <= p.real(0) {
                return Err(Error::InvalidInput(
                    "coordinate scale must be finite and positive real".into(),
                ));
            }
        }
        // Cancel common offsets exactly before MPFR evaluation.
        self.distance(&self.differences(source, target)?, p)
    }

    fn distance(&self, differences: &[Atom], p: Precision) -> Result<Float> {
        let mut sum = p.real(0);
        for difference in differences {
            let difference = p.eval(difference, &Default::default())?;
            let term = p
                .mul(
                    &difference,
                    &C::new(difference.re.clone(), -difference.im.clone()),
                )
                .re;
            sum = Float::with_val(p.bits, sum.as_raw() + term.as_raw());
        }
        Ok(sum)
    }
}

fn rational_squared_distance(differences: &[Atom]) -> Option<Rational> {
    let mut sum = Rational::from(0);
    for difference in differences {
        let AtomView::Num(number) = difference.as_view() else {
            return None;
        };
        let Coefficient::Complex(value) = number.get_coeff_view().to_owned() else {
            return None;
        };
        sum += value.re.clone() * &value.re + value.im.clone() * &value.im;
    }
    Some(sum)
}

impl<F: Fn(&CachedBoundary, &CachedPoint) -> Result<bool>> TransportCost for ScaledDistance<F> {
    fn cost(
        &self,
        source: &CachedBoundary,
        target: &CachedPoint,
        p: Precision,
    ) -> Result<Option<Float>> {
        if !(self.admissible)(source, target)? {
            return Ok(None);
        }
        Ok(Some(self.checked_distance(&source.point, target, p)?))
    }

    fn lower_bound(
        &self,
        source: &CachedBoundary,
        target: &CachedPoint,
        p: Precision,
    ) -> Result<Option<Float>> {
        // Exactly the rounded cost that `cost` returns for an admissible route;
        // no binary64 screening or approximate nearest-neighbour cutoff.
        Ok(Some(self.checked_distance(&source.point, target, p)?))
    }

    fn compare_tied_costs(
        &self,
        left: &CachedBoundary,
        right: &CachedBoundary,
        target: &CachedPoint,
        p: Precision,
    ) -> Result<std::cmp::Ordering> {
        use std::cmp::Ordering;
        let hits = compare_exact_hits(&left.point, &right.point, target)?;
        if hits != Ordering::Equal {
            return Ok(hits);
        }
        let left = self.differences(&left.point, target)?;
        let right = self.differences(&right.point, target)?;
        if let (Some(left), Some(right)) = (
            rational_squared_distance(&left),
            rational_squared_distance(&right),
        ) {
            return Ok(left.cmp(&right));
        }
        // Algebraic distances need numerical comparison. Accept a refined
        // ordering only if two increased precisions agree; unresolved ties
        // retain the cache's normal accuracy preference.
        let mut previous = Ordering::Equal;
        for factor in [2, 4] {
            let refined = Precision {
                bits: p
                    .bits
                    .checked_mul(factor)
                    .ok_or_else(|| Error::Limit("distance comparison precision overflow".into()))?,
            };
            let left = self.distance(&left, refined)?;
            let right = self.distance(&right, refined)?;
            if !left.is_finite() || !right.is_finite() {
                return Err(Error::InvalidInput(
                    "nonfinite refined transport distance".into(),
                ));
            }
            let comparison = left.partial_cmp(&right).ok_or_else(|| {
                Error::InvalidInput("nonfinite refined transport distance".into())
            })?;
            if factor == 4 && comparison == previous {
                return Ok(comparison);
            }
            previous = comparison;
        }
        Ok(Ordering::Equal)
    }
}

pub struct BoundaryQuery<'a> {
    pub identity: &'a BoundaryIdentity,
    pub target: &'a CachedPoint,
    pub range: EpsilonRange,
    pub verified_digits: u32,
    pub minimum_coordinate_bits: u32,
    pub include_detours: bool,
}
impl<'a> BoundaryQuery<'a> {
    pub fn new(
        identity: &'a BoundaryIdentity,
        target: &'a CachedPoint,
        range: EpsilonRange,
        digits: u32,
    ) -> Result<Self> {
        range.count()?;
        if digits == 0 {
            return Err(Error::InvalidInput(
                "positive requested boundary accuracy required".into(),
            ));
        }
        let coordinate_digits = digits
            .checked_add(10)
            .ok_or_else(|| Error::Limit("coordinate precision overflow".into()))?;
        Ok(Self {
            identity,
            target,
            range,
            verified_digits: digits,
            minimum_coordinate_bits: Precision::decimal(coordinate_digits)?.bits,
            include_detours: false,
        })
    }
}

pub struct BoundaryMatch<'a> {
    pub boundary: &'a CachedBoundary,
    pub cost: Float,
}

#[derive(Clone, Debug, Default)]
pub struct RustFlowCache {
    entries: Vec<CachedBoundary>,
    index: BTreeMap<String, usize>,
}

/// Compatibility name for the native RustFlow physical boundary cache.
pub type BoundaryCache = RustFlowCache;

impl RustFlowCache {
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    pub fn entries(&self) -> &[CachedBoundary] {
        &self.entries
    }

    pub fn insert(&mut self, boundary: CachedBoundary) -> Result<()> {
        self.insert_many(vec![boundary])
    }

    /// Validate and key the complete batch before changing the bank. Existing
    /// entries are indexed once; growing the bank does not repeatedly serialize
    /// or hash every old point, and no full-bank clone is needed for rollback.
    pub fn insert_many(&mut self, boundaries: Vec<CachedBoundary>) -> Result<()> {
        let pending = boundaries
            .into_iter()
            .map(|boundary| {
                boundary.validate()?;
                Ok((boundary.key()?, boundary))
            })
            .collect::<Result<Vec<_>>>()?;
        for (key, boundary) in pending {
            if let Some(&index) = self.index.get(&key) {
                let old = &mut self.entries[index];
                if (
                    boundary.accuracy.verified_digits,
                    boundary.accuracy.working_bits,
                ) > (old.accuracy.verified_digits, old.accuracy.working_bits)
                {
                    *old = boundary;
                }
            } else {
                self.index.insert(key, self.entries.len());
                self.entries.push(boundary);
            }
        }
        Ok(())
    }

    pub fn best<'a>(
        &'a self,
        query: &BoundaryQuery<'_>,
        policy: &dyn TransportCost,
        p: Precision,
    ) -> Result<Option<BoundaryMatch<'a>>> {
        Ok(self
            .best_excluding(query, policy, p, &Default::default())?
            .map(|(_, matched)| matched))
    }

    // Failed physical transports are transactional, so indices remain stable
    // until a successful attempt inserts its validated result batch.
    pub(crate) fn best_excluding<'a>(
        &'a self,
        query: &BoundaryQuery<'_>,
        policy: &dyn TransportCost,
        p: Precision,
        excluded: &std::collections::BTreeSet<usize>,
    ) -> Result<Option<(usize, BoundaryMatch<'a>)>> {
        self.best_excluding_with_germ_policy(query, policy, p, excluded, false)
    }
    pub(crate) fn best_excluding_with_germ_policy<'a>(
        &'a self,
        query: &BoundaryQuery<'_>,
        policy: &dyn TransportCost,
        p: Precision,
        excluded: &std::collections::BTreeSet<usize>,
        allow_sheet_transition: bool,
    ) -> Result<Option<(usize, BoundaryMatch<'a>)>> {
        query.range.count()?;
        query.identity.validate_point(query.target)?;
        if query.verified_digits == 0
            || query.target.coordinate_bits() < query.minimum_coordinate_bits
        {
            return Err(Error::Accuracy(
                "target coordinates lack the required precision".into(),
            ));
        }
        if p.bits < 2 {
            return Err(Error::InvalidInput(
                "invalid transport cost precision".into(),
            ));
        }
        let validate_cost = |cost: &Float, name: &str| -> Result<()> {
            if !cost.is_finite() || cost < &p.real(0) {
                return Err(Error::InvalidInput(format!(
                    "transport {name} must be finite and nonnegative"
                )));
            }
            Ok(())
        };
        let mut candidates = Vec::new();
        for (index, boundary) in self.entries.iter().enumerate() {
            if excluded.contains(&index) {
                continue;
            }
            if boundary.identity.key() != query.identity.key()
                || (!allow_sheet_transition
                    && boundary.point.root_germ() != query.target.root_germ())
                || !boundary.range.covers(query.range)
                || boundary.accuracy.verified_digits < query.verified_digits
                || boundary.point.coordinate_bits() < query.minimum_coordinate_bits
                || (!query.include_detours && boundary.kind == PointKind::ContourDetour)
            {
                continue;
            }
            let lower_bound = policy.lower_bound(boundary, query.target, p)?;
            if let Some(bound) = &lower_bound {
                validate_cost(bound, "cost lower bound")?;
            }
            candidates.push((index, boundary, lower_bound));
        }
        // Unknown bounds sort first and are all evaluated. Known bounds sort
        // stably; original insertion order resolves otherwise identical ties.
        candidates.sort_by(|left, right| {
            match (&left.2, &right.2) {
                (None, None) => std::cmp::Ordering::Equal,
                (None, Some(_)) => std::cmp::Ordering::Less,
                (Some(_), None) => std::cmp::Ordering::Greater,
                (Some(left), Some(right)) => left.partial_cmp(right).unwrap(),
            }
            .then_with(|| left.0.cmp(&right.0))
        });
        let mut best: Option<(usize, BoundaryMatch<'a>)> = None;
        for (index, boundary, lower_bound) in candidates {
            if let (Some(bound), Some((_, old))) = (&lower_bound, &best)
                && bound > &old.cost
            {
                break;
            }
            if let Some(cost) = policy.cost(boundary, query.target, p)? {
                validate_cost(&cost, "cost")?;
                if lower_bound.as_ref().is_some_and(|bound| bound > &cost) {
                    return Err(Error::InvalidInput(
                        "transport cost is smaller than its advertised lower bound".into(),
                    ));
                }
                let replace = if let Some((old_index, old)) = &best {
                    if cost == old.cost {
                        let comparison =
                            policy.compare_tied_costs(boundary, old.boundary, query.target, p)?;
                        comparison.is_lt()
                            || (comparison.is_eq()
                                && (boundary.accuracy.verified_digits
                                    > old.boundary.accuracy.verified_digits
                                    || (boundary.accuracy.verified_digits
                                        == old.boundary.accuracy.verified_digits
                                        && index < *old_index)))
                    } else {
                        cost < old.cost
                    }
                } else {
                    true
                };
                if replace {
                    best = Some((index, BoundaryMatch { boundary, cost }));
                }
            }
        }
        Ok(best)
    }

    /// Retain accepted Taylor endpoints in physical coordinates. The callback
    /// must classify each endpoint and provide independently established evidence
    /// for that endpoint; `None` omits an unverified endpoint. Final-result accuracy
    /// must not be applied automatically to intermediate segments. Insertions are
    /// transactional if a later endpoint fails validation.
    pub fn insert_trajectory(
        &mut self,
        identity: &BoundaryIdentity,
        path: &KinematicPath,
        solution: &EpsilonSolution,
        evidence: impl FnMut(
            usize,
            &crate::ode::TaylorSegment,
        ) -> Result<Option<(PointKind, BoundaryAccuracy)>>,
    ) -> Result<usize> {
        let pending = Self::trajectory_boundaries(identity, path, solution, evidence)?;
        let count = pending.len();
        self.insert_many(pending)?;
        Ok(count)
    }

    /// Prepare validated checkpoints without mutating a bank. Append a final
    /// destination or other boundaries and commit them together with `insert_many`.
    pub fn trajectory_boundaries(
        identity: &BoundaryIdentity,
        path: &KinematicPath,
        solution: &EpsilonSolution,
        evidence: impl FnMut(
            usize,
            &crate::ode::TaylorSegment,
        ) -> Result<Option<(PointKind, BoundaryAccuracy)>>,
    ) -> Result<Vec<CachedBoundary>> {
        Self::trajectory_boundaries_with_germ(identity, path, solution, None, evidence)
    }

    pub(crate) fn trajectory_boundaries_with_germ(
        identity: &BoundaryIdentity,
        path: &KinematicPath,
        solution: &EpsilonSolution,
        germ: Option<&RootGerm>,
        evidence: impl FnMut(
            usize,
            &crate::ode::TaylorSegment,
        ) -> Result<Option<(PointKind, BoundaryAccuracy)>>,
    ) -> Result<Vec<CachedBoundary>> {
        Self::trajectory_boundaries_with_germs(
            identity,
            path,
            solution,
            |_| Ok(germ.cloned()),
            evidence,
        )
    }
    pub(crate) fn trajectory_boundaries_with_germs(
        identity: &BoundaryIdentity,
        path: &KinematicPath,
        solution: &EpsilonSolution,
        mut germ_at: impl FnMut(usize) -> Result<Option<RootGerm>>,
        mut evidence: impl FnMut(
            usize,
            &crate::ode::TaylorSegment,
        ) -> Result<Option<(PointKind, BoundaryAccuracy)>>,
    ) -> Result<Vec<CachedBoundary>> {
        path.validate()?;
        if path
            .coordinates
            .keys()
            .copied()
            .collect::<std::collections::BTreeSet<_>>()
            != identity.0.variables
            || path.parameter == identity.0.system.epsilon()
        {
            return Err(Error::InvalidInput(
                "trajectory path does not match boundary identity".into(),
            ));
        }
        if solution.coefficients.is_empty()
            || solution
                .coefficients
                .iter()
                .any(|r| r.len() != identity.dimension())
        {
            return Err(Error::InvalidInput("trajectory solution dimensions".into()));
        }
        let last = solution
            .leading
            .checked_add(
                i32::try_from(solution.coefficients.len())
                    .map_err(|_| Error::Limit("epsilon range overflow".into()))?
                    - 1,
            )
            .ok_or_else(|| Error::Limit("epsilon range overflow".into()))?;
        let range = EpsilonRange::new(solution.leading, last)?;
        let mut pending = Vec::new();
        for (index, segment) in solution.segments.iter().enumerate() {
            let Some((kind, accuracy)) = evidence(index, segment)? else {
                continue;
            };
            if !segment.end.re.is_finite() || !segment.end.im.is_finite() {
                return Err(Error::InvalidInput("nonfinite trajectory endpoint".into()));
            }
            let parameter = Atom::num(symbolica::domains::float::Complex::new(
                segment.end.re.to_rational(),
                segment.end.im.to_rational(),
            ));
            let parameters = BTreeMap::from([(Atom::var(path.parameter), parameter.clone())]);
            let coordinates = path
                .coordinates
                .iter()
                .map(|(&s, a)| {
                    (
                        s,
                        crate::family::substitute(a, &parameters)
                            .together()
                            .cancel(),
                    )
                })
                .collect();
            let point = CachedPoint::Derived {
                coordinates,
                working_bits: segment.working_bits,
                provenance: format!(
                    "exact image of accepted trajectory endpoint {index}; rounded parameter={}; path={}",
                    parameter.to_canonical_string(),
                    path.coordinates
                        .iter()
                        .map(|(&s, a)| format!(
                            "{}={}",
                            Atom::var(s).to_canonical_string(),
                            a.to_canonical_string()
                        ))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
            };
            let point = if let Some(germ) = germ_at(index)? {
                point.with_root_germ(germ)?
            } else {
                point
            };
            let boundary = CachedBoundary {
                identity: identity.clone(),
                point,
                kind,
                range,
                coefficients: solution.evaluate_segment(index, &segment.end)?,
                accuracy,
            };
            boundary.validate()?;
            pending.push(boundary);
        }
        Ok(pending)
    }
}

// Storage uses native Atom serialization and precision-bearing Symbolica Float
// serialization. Identity hashing instead uses canonical expressions, independent
// of process-local Symbol allocation order.
fn atom_bytes(a: &Atom) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    a.export(&mut out)?;
    Ok(out)
}
fn atom_read(bytes: &[u8]) -> Result<Atom> {
    let mut source = bytes;
    let result = Atom::import(&mut source, None).map_err(|e| Error::Cache(e.to_string()))?;
    if !source.is_empty() {
        return Err(Error::Cache("trailing cached Atom data".into()));
    }
    Ok(result)
}
fn symbol_bytes(s: Symbol) -> Result<Vec<u8>> {
    atom_bytes(&Atom::var(s))
}
fn symbol_read(bytes: &[u8]) -> Result<Symbol> {
    match atom_read(bytes)?.as_view() {
        AtomView::Var(v) => Ok(v.get_symbol()),
        _ => Err(Error::Cache(
            "cached coordinate name is not a symbol".into(),
        )),
    }
}

type StoredAtom = Vec<u8>;
type StoredMatrix = Vec<Vec<StoredAtom>>;

#[derive(Serialize, Deserialize)]
enum StoredSystem {
    Dense {
        derivatives: Vec<(StoredAtom, StoredMatrix)>,
    },
    Canonical {
        variables: Vec<StoredAtom>,
        letters: Vec<StoredAtom>,
        matrices: Vec<StoredMatrix>,
    },
}
#[derive(Serialize, Deserialize)]
enum StoredContinuation {
    PrescribedAffineV1 {
        domain: String,
        unprescribed_plus: bool,
        prescriptions: Vec<(StoredAtom, bool)>,
    },
}
impl StoredContinuation {
    fn encode(value: &PhysicalContinuation) -> Result<Self> {
        Ok(Self::PrescribedAffineV1 {
            domain: value.domain.clone(),
            unprescribed_plus: matches!(value.unprescribed_side, Prescription::PlusI0),
            prescriptions: value
                .prescriptions
                .iter()
                .map(|p| {
                    Ok((
                        atom_bytes(&p.polynomial)?,
                        matches!(p.prescription, Prescription::PlusI0),
                    ))
                })
                .collect::<Result<_>>()?,
        })
    }
    fn decode(self) -> Result<PhysicalContinuation> {
        let Self::PrescribedAffineV1 {
            domain,
            unprescribed_plus,
            prescriptions,
        } = self;
        let side = |plus| {
            if plus {
                Prescription::PlusI0
            } else {
                Prescription::MinusI0
            }
        };
        Ok(PhysicalContinuation {
            domain,
            unprescribed_side: side(unprescribed_plus),
            prescriptions: prescriptions
                .into_iter()
                .map(|(a, plus)| {
                    Ok(crate::contour::PolynomialPrescription {
                        polynomial: atom_read(&a)?,
                        prescription: side(plus),
                    })
                })
                .collect::<Result<_>>()?,
        })
    }
}
#[derive(Serialize, Deserialize)]
struct StoredIdentity {
    key: String,
    roots: Vec<(StoredAtom, StoredAtom)>,
    epsilon: StoredAtom,
    system: StoredSystem,
    basis: Vec<StoredAtom>,
    normalization: StoredAtom,
    plus_i0: bool,
    domain: String,
    conditions: Vec<StoredAtom>,
    continuation: Option<StoredContinuation>,
}
fn matrix_bytes(matrix: &[Vec<Atom>]) -> Result<StoredMatrix> {
    matrix
        .iter()
        .map(|r| r.iter().map(atom_bytes).collect())
        .collect()
}
fn matrix_read(matrix: StoredMatrix) -> Result<Vec<Vec<Atom>>> {
    matrix
        .into_iter()
        .map(|r| r.into_iter().map(|a| atom_read(&a)).collect())
        .collect()
}
impl StoredIdentity {
    fn encode(identity: &BoundaryIdentity) -> Result<Self> {
        let d = &identity.0;
        Ok(Self {
            key: d.key.clone(),
            roots: d
                .roots
                .iter()
                .map(|r| Ok((symbol_bytes(r.symbol)?, atom_bytes(&r.radicand)?)))
                .collect::<Result<_>>()?,
            epsilon: symbol_bytes(d.system.epsilon())?,
            system: match &d.system {
                IdentitySystem::Dense(system) => StoredSystem::Dense {
                    derivatives: system
                        .derivatives
                        .iter()
                        .map(|(&s, m)| Ok((symbol_bytes(s)?, matrix_bytes(m)?)))
                        .collect::<Result<_>>()?,
                },
                IdentitySystem::Canonical(system) => StoredSystem::Canonical {
                    variables: system
                        .variables()
                        .iter()
                        .map(|&s| symbol_bytes(s))
                        .collect::<Result<_>>()?,
                    letters: system
                        .letters()
                        .iter()
                        .map(atom_bytes)
                        .collect::<Result<_>>()?,
                    matrices: system
                        .constant_matrices()
                        .iter()
                        .map(|m| matrix_bytes(m))
                        .collect::<Result<_>>()?,
                },
            },
            basis: d.basis.iter().map(atom_bytes).collect::<Result<_>>()?,
            normalization: atom_bytes(&d.normalization)?,
            plus_i0: matches!(d.prescription, Prescription::PlusI0),
            domain: d.domain.clone(),
            conditions: d.conditions.iter().map(atom_bytes).collect::<Result<_>>()?,
            continuation: d
                .continuation
                .as_ref()
                .map(StoredContinuation::encode)
                .transpose()?,
        })
    }
    fn decode(self) -> Result<BoundaryIdentity> {
        let epsilon = symbol_read(&self.epsilon)?;
        let basis = self
            .basis
            .iter()
            .map(|a| atom_read(a))
            .collect::<Result<Vec<_>>>()?;
        let roots = self
            .roots
            .into_iter()
            .map(|(symbol, radicand)| {
                Ok(SquareRoot {
                    symbol: symbol_read(&symbol)?,
                    radicand: atom_read(&radicand)?,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let normalization = atom_read(&self.normalization)?;
        let prescription = if self.plus_i0 {
            Prescription::PlusI0
        } else {
            Prescription::MinusI0
        };
        let conditions = self
            .conditions
            .iter()
            .map(|a| atom_read(a))
            .collect::<Result<Vec<_>>>()?;
        let identity = match self.system {
            StoredSystem::Dense { derivatives } => {
                let system = KinematicSystem {
                    epsilon,
                    derivatives: derivatives
                        .into_iter()
                        .map(|(s, m)| Ok((symbol_read(&s)?, matrix_read(m)?)))
                        .collect::<Result<_>>()?,
                };
                BoundaryIdentity::build(
                    &system,
                    &roots,
                    &basis,
                    &normalization,
                    prescription,
                    &self.domain,
                    &conditions,
                )?
            }
            StoredSystem::Canonical {
                variables,
                letters,
                matrices,
            } => {
                let system = CanonicalAlgebraicSystem::new(
                    epsilon,
                    &variables
                        .iter()
                        .map(|s| symbol_read(s))
                        .collect::<Result<Vec<_>>>()?,
                    &letters
                        .iter()
                        .map(|a| atom_read(a))
                        .collect::<Result<Vec<_>>>()?,
                    &matrices
                        .into_iter()
                        .map(matrix_read)
                        .collect::<Result<Vec<_>>>()?,
                    roots,
                )?;
                BoundaryIdentity::with_canonical_conditions(
                    &system,
                    &basis,
                    &normalization,
                    prescription,
                    &self.domain,
                    &conditions,
                )?
            }
        };
        let identity = if let Some(continuation) = self.continuation {
            identity.with_prescribed_continuation(continuation.decode()?)?
        } else {
            identity
        };
        if identity.key() != self.key {
            return Err(Error::Cache("cached system/basis identity mismatch".into()));
        }
        Ok(identity)
    }
}
#[derive(Serialize, Deserialize)]
enum StoredPoint {
    Algebraic {
        point: Box<StoredPoint>,
        germ: Vec<(StoredAtom, RootSheet)>,
    },
    Exact(Vec<(Vec<u8>, Vec<u8>)>),
    Derived {
        coordinates: Vec<(Vec<u8>, Vec<u8>)>,
        working_bits: u32,
        provenance: String,
    },
    Numerical {
        coordinates: Vec<(Vec<u8>, C)>,
        working_bits: u32,
        provenance: String,
    },
}
impl StoredPoint {
    fn encode(point: &CachedPoint) -> Result<Self> {
        match point {
            CachedPoint::Algebraic { point, germ } => Ok(Self::Algebraic {
                point: Box::new(Self::encode(point)?),
                germ: germ
                    .sheets
                    .iter()
                    .map(|(&s, &sheet)| Ok((symbol_bytes(s)?, sheet)))
                    .collect::<Result<_>>()?,
            }),
            CachedPoint::Exact(c) => Ok(Self::Exact(
                c.iter()
                    .map(|(&s, a)| Ok((symbol_bytes(s)?, atom_bytes(a)?)))
                    .collect::<Result<_>>()?,
            )),
            CachedPoint::Derived {
                coordinates,
                working_bits,
                provenance,
            } => Ok(Self::Derived {
                coordinates: coordinates
                    .iter()
                    .map(|(&s, a)| Ok((symbol_bytes(s)?, atom_bytes(a)?)))
                    .collect::<Result<_>>()?,
                working_bits: *working_bits,
                provenance: provenance.clone(),
            }),
            CachedPoint::Numerical {
                coordinates,
                working_bits,
                provenance,
            } => Ok(Self::Numerical {
                coordinates: coordinates
                    .iter()
                    .map(|(&s, a)| Ok((symbol_bytes(s)?, a.clone())))
                    .collect::<Result<_>>()?,
                working_bits: *working_bits,
                provenance: provenance.clone(),
            }),
        }
    }
    fn decode(self) -> Result<CachedPoint> {
        match self {
            Self::Algebraic { point, germ } => point.decode()?.with_root_germ(RootGerm {
                sheets: germ
                    .into_iter()
                    .map(|(s, sheet)| Ok((symbol_read(&s)?, sheet)))
                    .collect::<Result<_>>()?,
            }),
            Self::Exact(c) => Ok(CachedPoint::Exact(
                c.into_iter()
                    .map(|(s, a)| Ok((symbol_read(&s)?, atom_read(&a)?)))
                    .collect::<Result<_>>()?,
            )),
            Self::Derived {
                coordinates,
                working_bits,
                provenance,
            } => Ok(CachedPoint::Derived {
                coordinates: coordinates
                    .into_iter()
                    .map(|(s, a)| Ok((symbol_read(&s)?, atom_read(&a)?)))
                    .collect::<Result<_>>()?,
                working_bits,
                provenance,
            }),
            Self::Numerical {
                coordinates,
                working_bits,
                provenance,
            } => Ok(CachedPoint::Numerical {
                coordinates: coordinates
                    .into_iter()
                    .map(|(s, a)| Ok((symbol_read(&s)?, a)))
                    .collect::<Result<_>>()?,
                working_bits,
                provenance,
            }),
        }
    }
}
#[derive(Serialize, Deserialize)]
struct StoredBoundary {
    identity: usize,
    point: StoredPoint,
    kind: PointKind,
    range: EpsilonRange,
    coefficients: Vec<Vec<C>>,
    accuracy: BoundaryAccuracy,
}
#[derive(Serialize, Deserialize)]
struct StoredCache {
    identities: Vec<StoredIdentity>,
    boundaries: Vec<StoredBoundary>,
}
#[derive(Serialize, Deserialize)]
struct Envelope {
    version: u32,
    implementation: String,
    dependencies: String,
    digest: String,
    payload: Vec<u8>,
}

impl RustFlowCache {
    /// Atomically save this snapshot. Concurrent snapshots may replace one another;
    /// callers sharing a directory should serialize writes or merge before saving.
    pub fn save(&self, directory: &Path) -> Result<()> {
        let mut identities = Vec::new();
        let mut keys = BTreeMap::new();
        let mut boundaries = Vec::new();
        for boundary in &self.entries {
            boundary.validate()?;
            let identity = if let Some(index) = keys.get(boundary.identity.key()) {
                *index
            } else {
                let index = identities.len();
                identities.push(StoredIdentity::encode(&boundary.identity)?);
                keys.insert(boundary.identity.key().to_owned(), index);
                index
            };
            boundaries.push(StoredBoundary {
                identity,
                point: StoredPoint::encode(&boundary.point)?,
                kind: boundary.kind,
                range: boundary.range,
                coefficients: boundary.coefficients.clone(),
                accuracy: boundary.accuracy.clone(),
            });
        }
        let payload = bincode::serde::encode_to_vec(
            &StoredCache {
                identities,
                boundaries,
            },
            bincode::config::standard(),
        )
        .map_err(|e| Error::Cache(e.to_string()))?;
        let envelope = Envelope {
            version: VERSION,
            implementation: env!("PORT_SOURCE_DIGEST").into(),
            dependencies: env!("DEPENDENCY_SOURCE_DIGEST").into(),
            digest: blake3::hash(&payload).to_hex().to_string(),
            payload,
        };
        std::fs::create_dir_all(directory)?;
        let temporary = directory.join(format!(
            ".physical-boundaries.{}.{}.tmp",
            std::process::id(),
            NEXT_FILE.fetch_add(1, Ordering::Relaxed)
        ));
        let write = (|| -> Result<()> {
            let mut file = std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&temporary)?;
            let mut writer = std::io::BufWriter::new(&mut file);
            writer.write_all(MAGIC)?;
            bincode::serde::encode_into_std_write(
                &envelope,
                &mut writer,
                bincode::config::standard(),
            )
            .map_err(|e| Error::Cache(e.to_string()))?;
            writer.flush()?;
            writer.get_ref().sync_all()?;
            std::fs::rename(&temporary, directory.join(FILE))?;
            Ok(())
        })();
        if write.is_err() {
            let _ = std::fs::remove_file(temporary);
        }
        write
    }

    /// Load and validate a snapshot; missing storage gives an empty cache. Corrupt
    /// or incompatible entries are errors, never silently treated as cache misses.
    pub fn load(directory: &Path) -> Result<Self> {
        let bytes = match std::fs::read(directory.join(FILE)) {
            Ok(bytes) => bytes,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(Self::default()),
            Err(e) => return Err(e.into()),
        };
        fn decode<T: serde::de::DeserializeOwned>(bytes: &[u8]) -> Result<T> {
            let (value, count) =
                bincode::serde::decode_from_slice(bytes, bincode::config::standard())
                    .map_err(|e| Error::Cache(e.to_string()))?;
            if count != bytes.len() {
                return Err(Error::Cache("trailing boundary-cache bytes".into()));
            }
            Ok(value)
        }
        let envelope: Envelope = decode(
            bytes
                .strip_prefix(MAGIC)
                .ok_or_else(|| Error::Cache("incompatible boundary cache codec".into()))?,
        )?;
        if envelope.version != VERSION
            || envelope.implementation != env!("PORT_SOURCE_DIGEST")
            || envelope.dependencies != env!("DEPENDENCY_SOURCE_DIGEST")
            || envelope.digest != blake3::hash(&envelope.payload).to_hex().as_str()
        {
            return Err(Error::Cache(
                "incompatible or corrupt boundary cache".into(),
            ));
        }
        let stored: StoredCache = decode(&envelope.payload)?;
        let identities = stored
            .identities
            .into_iter()
            .map(StoredIdentity::decode)
            .collect::<Result<Vec<_>>>()?;
        let mut cache = Self::default();
        for boundary in stored.boundaries {
            let identity = identities
                .get(boundary.identity)
                .ok_or_else(|| Error::Cache("unknown cached boundary identity".into()))?
                .clone();
            cache.insert(CachedBoundary {
                identity,
                point: boundary.point.decode()?,
                kind: boundary.kind,
                range: boundary.range,
                coefficients: boundary.coefficients,
                accuracy: boundary.accuracy,
            })?;
        }
        Ok(cache)
    }
}
