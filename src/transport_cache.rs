//! Physical boundary reuse for supplied differential systems.
//!
//! Exact coordinates and numerically chosen coordinates remain distinct. A
//! precision label is evidence supplied by the caller, never inferred from the
//! number of MPFR working bits. Path admissibility is required during selection:
//! distance alone cannot decide which continuation reaches the intended sheet.
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

const VERSION: u32 = 1;
const MAGIC: &[u8] = b"AMFLOW-BOUNDARIES\0\x01";
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

/// Exact input, an exact path image of a numerical parameter, or independently
/// rounded numerical coordinates. The origin of a point remains explicit.
#[derive(Clone, Debug)]
pub enum CachedPoint {
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
    pub fn validate(&self) -> Result<()> {
        match self {
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
            Self::Exact(_) => u32::MAX,
            Self::Numerical { working_bits, .. } | Self::Derived { working_bits, .. } => {
                *working_bits
            }
        }
    }
    fn key(&self) -> Result<String> {
        match self {
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
struct IdentityData {
    key: String,
    system: KinematicSystem,
    basis: Vec<Atom>,
    normalization: Atom,
    prescription: Prescription,
    domain: String,
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
        system.validate()?;
        let n = system.derivatives.values().next().unwrap().len();
        if basis.len() != n || domain.trim().is_empty() || normalization.is_zero() {
            return Err(Error::InvalidInput("boundary identity needs an ordered basis, nonzero normalization, and homotopy domain".into()));
        }
        for a in basis.iter().chain(std::iter::once(normalization)) {
            exact_metadata(a.as_view())?;
        }
        let mut derivatives = system
            .derivatives
            .iter()
            .map(|(&s, m)| {
                (
                    Atom::var(s).to_canonical_string(),
                    m.iter()
                        .map(|r| {
                            r.iter()
                                .map(AtomCore::to_canonical_string)
                                .collect::<Vec<_>>()
                        })
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<Vec<_>>();
        derivatives.sort_by(|a, b| a.0.cmp(&b.0));
        let key = fingerprint(&(
            VERSION,
            env!("DEPENDENCY_SOURCE_DIGEST"),
            env!("PORT_SOURCE_DIGEST"),
            Atom::var(system.epsilon).to_canonical_string(),
            derivatives,
            basis
                .iter()
                .map(AtomCore::to_canonical_string)
                .collect::<Vec<_>>(),
            normalization.to_canonical_string(),
            matches!(prescription, Prescription::PlusI0),
            domain,
        ))?;
        Ok(Self(Arc::new(IdentityData {
            key,
            system: system.clone(),
            basis: basis.to_vec(),
            normalization: normalization.clone(),
            prescription,
            domain: domain.into(),
        })))
    }
    pub fn key(&self) -> &str {
        &self.0.key
    }
    pub fn dimension(&self) -> usize {
        self.0.basis.len()
    }

    fn validate_point(&self, point: &CachedPoint) -> Result<()> {
        let coordinates = point.rounded_coordinates_as_exact()?;
        if !coordinates.keys().eq(self.0.system.derivatives.keys()) {
            return Err(Error::InvalidInput(
                "cached point coordinates do not match the system".into(),
            ));
        }
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
        for a in self.0.system.derivatives.values().flatten().flatten() {
            regular(a.as_view(), &rules)?;
        }
        Ok(())
    }
}

/// Caller-recorded accuracy evidence. Comparison errors are consistency estimates,
/// not rigorous error bounds. The provenance must explain how accuracy was obtained.
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
}

/// Squared Euclidean distance in scaled kinematic coordinates, computed in MPFR.
/// The mandatory callback checks path/sheet/pole admissibility before pricing.
pub struct ScaledDistance<F> {
    pub scales: BTreeMap<Symbol, Atom>,
    pub admissible: F,
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
        let source = source.point.evaluate(p)?;
        let target = target.evaluate(p)?;
        if !source.keys().eq(target.keys()) || self.scales.keys().any(|s| !source.contains_key(s)) {
            return Err(Error::InvalidInput(
                "distance coordinates/scales do not match".into(),
            ));
        }
        let mut sum = p.real(0);
        for (variable, value) in source {
            let scale = if let Some(a) = self.scales.get(&variable) {
                exact_constant(a.as_view())?;
                let scale = p.eval(a, &Default::default())?;
                if !p.finite(&scale) || scale.im != p.real(0) || scale.re <= p.real(0) {
                    return Err(Error::InvalidInput(
                        "coordinate scale must be finite and positive real".into(),
                    ));
                }
                scale
            } else {
                p.i(1)
            };
            let difference = p.div(&p.sub(&value, &target[&variable]), &scale);
            let term = p
                .mul(
                    &difference,
                    &C::new(difference.re.clone(), -difference.im.clone()),
                )
                .re;
            sum = Float::with_val(p.bits, sum.as_raw() + term.as_raw());
        }
        Ok(Some(sum))
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
        query.range.count()?;
        query.identity.validate_point(query.target)?;
        if query.verified_digits == 0
            || query.target.coordinate_bits() < query.minimum_coordinate_bits
        {
            return Err(Error::Accuracy(
                "target coordinates lack the required precision".into(),
            ));
        }
        let mut best: Option<BoundaryMatch<'a>> = None;
        for boundary in &self.entries {
            if boundary.identity.key() != query.identity.key()
                || !boundary.range.covers(query.range)
                || boundary.accuracy.verified_digits < query.verified_digits
                || boundary.point.coordinate_bits() < query.minimum_coordinate_bits
                || (!query.include_detours && boundary.kind == PointKind::ContourDetour)
            {
                continue;
            }
            if let Some(cost) = policy.cost(boundary, query.target, p)? {
                if !cost.is_finite() || cost < p.real(0) {
                    return Err(Error::InvalidInput(
                        "transport cost must be finite and nonnegative".into(),
                    ));
                }
                if best.as_ref().is_none_or(|old| {
                    cost < old.cost
                        || (cost == old.cost
                            && boundary.accuracy.verified_digits
                                > old.boundary.accuracy.verified_digits)
                }) {
                    best = Some(BoundaryMatch { boundary, cost });
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
        mut evidence: impl FnMut(
            usize,
            &crate::ode::TaylorSegment,
        ) -> Result<Option<(PointKind, BoundaryAccuracy)>>,
    ) -> Result<Vec<CachedBoundary>> {
        path.validate()?;
        if !path
            .coordinates
            .keys()
            .eq(identity.0.system.derivatives.keys())
            || path.parameter == identity.0.system.epsilon
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
            let boundary = CachedBoundary {
                identity: identity.clone(),
                point: CachedPoint::Derived {
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
                },
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
struct StoredIdentity {
    key: String,
    epsilon: Vec<u8>,
    derivatives: Vec<(StoredAtom, StoredMatrix)>,
    basis: Vec<Vec<u8>>,
    normalization: Vec<u8>,
    plus_i0: bool,
    domain: String,
}
impl StoredIdentity {
    fn encode(identity: &BoundaryIdentity) -> Result<Self> {
        let d = &identity.0;
        Ok(Self {
            key: d.key.clone(),
            epsilon: symbol_bytes(d.system.epsilon)?,
            derivatives: d
                .system
                .derivatives
                .iter()
                .map(|(&s, m)| {
                    Ok((
                        symbol_bytes(s)?,
                        m.iter()
                            .map(|r| r.iter().map(atom_bytes).collect())
                            .collect::<Result<_>>()?,
                    ))
                })
                .collect::<Result<_>>()?,
            basis: d.basis.iter().map(atom_bytes).collect::<Result<_>>()?,
            normalization: atom_bytes(&d.normalization)?,
            plus_i0: matches!(d.prescription, Prescription::PlusI0),
            domain: d.domain.clone(),
        })
    }
    fn decode(self) -> Result<BoundaryIdentity> {
        let system = KinematicSystem {
            epsilon: symbol_read(&self.epsilon)?,
            derivatives: self
                .derivatives
                .into_iter()
                .map(|(s, m)| {
                    Ok((
                        symbol_read(&s)?,
                        m.into_iter()
                            .map(|r| r.into_iter().map(|a| atom_read(&a)).collect())
                            .collect::<Result<_>>()?,
                    ))
                })
                .collect::<Result<_>>()?,
        };
        let basis = self
            .basis
            .iter()
            .map(|a| atom_read(a))
            .collect::<Result<Vec<_>>>()?;
        let identity = BoundaryIdentity::new(
            &system,
            &basis,
            &atom_read(&self.normalization)?,
            if self.plus_i0 {
                Prescription::PlusI0
            } else {
                Prescription::MinusI0
            },
            &self.domain,
        )?;
        if identity.key() != self.key {
            return Err(Error::Cache("cached system/basis identity mismatch".into()));
        }
        Ok(identity)
    }
}
#[derive(Serialize, Deserialize)]
enum StoredPoint {
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
