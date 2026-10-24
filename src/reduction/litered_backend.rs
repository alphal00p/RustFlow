//! RustRed LiteRed-lane reduction backend (RR-8 interim backend, UNCERTIFIED).
//!
//! Per exact symbolic-epsilon family this backend keeps one RustRed
//! [`DynamicLiteRedSession`]: candidate formulas are generated in process by
//! RustRed for the whole physical root downset, the RustRed sector-symmetry
//! census compresses the candidate terminals to a fixed basis, and every
//! request is reduced onto that basis through RustRed's modular application
//! DAG. Coefficients are exact: rational functions of `eta` at a sampled
//! epsilon (multi-prime Thiele/CRT/MQRR reconstruction, verified at a fresh
//! prime), or rationals when no free parameter remains.
//!
//! Contract with the differential closure: the basis is fixed per family
//! before the first request (it is the compressed family basis), so a
//! residual never later acquires a rule and earlier rules never change.
//! Requests that are basis integrals are returned as residuals.
//!
//! Limits: symbolic epsilon (`reduce`) needs bivariate reconstruction and is
//! not supported yet; cut families are not supported. The candidate formulas
//! carry no source-replay provenance and the basis is not proved minimal.
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

use rustred::family::IntegralKey;
use rustred::solver::litered::{
    ColumnOrder, DynamicLiteRedSession, LiteRedOptions, RationalFunction, SymmetryRows,
};
use symbolica::prelude::*;

use super::{LinearCombination, Reduction, ReductionBackend};
use crate::family::ConvertedFamily;
use crate::{Error, Integral, IntegralFamily, Progress, Result, RunContext};

/// Options of [`LiteRedRustRedBackend`].
#[derive(Clone, Copy, Debug)]
pub struct LiteRedBackendOptions {
    /// Threads for candidate generation and for concurrent modular samples.
    pub workers: usize,
    /// Census rows over the full representative groups (default) or only
    /// their generators.
    pub group_orbit_rows: bool,
    /// Basis order: `(lines, dots, rank, lex)` (default) or RustRed's order.
    pub lines_dots_rank: bool,
    pub seed: u64,
    pub max_primes: usize,
    /// Generate candidates with the coordinate priority "first ISP first,
    /// remaining ISPs last" (`rustred::solver::litered::isps_at_ends_priority`).
    /// On P14 spoke the default ISPs-last order stalls in the sector solver's
    /// exceptional geometry (sector 217); this priority generates in ~1 min.
    pub isps_at_ends: bool,
}

impl Default for LiteRedBackendOptions {
    fn default() -> Self {
        Self {
            workers: 8,
            group_orbit_rows: true,
            lines_dots_rank: true,
            seed: 0x5eed_1173,
            max_primes: 256,
            isps_at_ends: true,
        }
    }
}

struct FamilyLane {
    converted: ConvertedFamily,
    lane: DynamicLiteRedSession,
    /// Position of the epsilon parameter in the RustRed parameter order.
    epsilon: Option<usize>,
    /// RustFlow atom of every RustRed parameter, in parameter order.
    parameters: Vec<Atom>,
    basis: Vec<Integral>,
    nodes: BTreeMap<Integral, u32>,
}

/// UNCERTIFIED RustRed LiteRed-lane backend.
pub struct LiteRedRustRedBackend {
    pub options: LiteRedBackendOptions,
    families: Mutex<BTreeMap<String, Arc<Mutex<FamilyLane>>>>,
}

impl LiteRedRustRedBackend {
    pub fn new(options: LiteRedBackendOptions) -> Self {
        Self {
            options,
            families: Mutex::new(BTreeMap::new()),
        }
    }

    fn lane_options(&self) -> LiteRedOptions {
        LiteRedOptions {
            order: if self.options.lines_dots_rank {
                ColumnOrder::LinesDotsRank
            } else {
                ColumnOrder::Program
            },
            rows: if self.options.group_orbit_rows {
                SymmetryRows::GroupOrbit
            } else {
                SymmetryRows::Generators
            },
            seed: self.options.seed,
            max_primes: self.options.max_primes,
            workers: self.options.workers,
            ..Default::default()
        }
    }

    fn family_lane(
        &self,
        family: &IntegralFamily,
        context: &RunContext,
    ) -> Result<Arc<Mutex<FamilyLane>>> {
        let converted = family.convert()?;
        let key = converted.family.fingerprint().to_owned();
        if let Some(lane) = self.families.lock().expect("family map").get(&key) {
            return Ok(Arc::clone(lane));
        }
        context.cancellation.check()?;
        let root: Vec<bool> = (0..family.propagators.len())
            .map(|slot| slot < family.physical_propagators)
            .collect();
        // RustRed families are not Clone: convert once more for the lane.
        let rustred_family = Arc::new(family.convert()?.family);
        let priority = self
            .options
            .isps_at_ends
            .then(|| rustred::solver::litered::isps_at_ends_priority(&root))
            .flatten();
        let (mut lane, _report) = DynamicLiteRedSession::generate(
            rustred_family,
            &root,
            self.options.workers,
            true,
            self.lane_options(),
            priority.as_deref(),
        )
        .map_err(|e| Error::Reduction(format!("RustRed LiteRed lane: {e}")))?;
        let compression = lane
            .compress(&[])
            .map_err(|e| Error::Reduction(format!("RustRed LiteRed compression: {e}")))?;
        let context_parameters = converted.family.coefficient_context();
        let mut parameters = Vec::new();
        let mut epsilon = None;
        for (position, name) in lane.parameter_names().iter().enumerate() {
            let rustred_atom = context_parameters
                .parameter(name)
                .ok_or_else(|| Error::Reduction(format!("unknown RustRed parameter {name}")))?
                .to_expression();
            let original = converted
                .reverse
                .get(&rustred_atom)
                .cloned()
                .ok_or_else(|| Error::Reduction(format!("unmapped RustRed parameter {name}")))?;
            if original == Atom::var(family.epsilon) {
                epsilon = Some(position);
            }
            parameters.push(original);
        }
        let basis = compression
            .masters
            .iter()
            .map(|key| integral(key))
            .collect::<Result<Vec<_>>>()?;
        let state = Arc::new(Mutex::new(FamilyLane {
            converted,
            lane,
            epsilon,
            parameters,
            basis,
            nodes: BTreeMap::new(),
        }));
        let mut families = self.families.lock().expect("family map");
        Ok(Arc::clone(families.entry(key).or_insert(state)))
    }

    fn reduce_with(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        epsilon: Option<&Rational>,
        context: &RunContext,
    ) -> Result<Reduction> {
        context.emit(Progress::Reduction {
            integrals: targets.len(),
        })?;
        for target in targets {
            family.validate_integral(target)?;
        }
        let lane = self.family_lane(family, context)?;
        let mut lane = lane.lock().expect("family lane");
        let FamilyLane {
            converted: _,
            lane: session,
            epsilon: epsilon_position,
            parameters,
            basis,
            nodes,
        } = &mut *lane;
        let new: Vec<Integral> = targets
            .iter()
            .filter(|t| !nodes.contains_key(*t))
            .cloned()
            .collect();
        if !new.is_empty() {
            let keys = new.iter().map(key).collect::<Result<Vec<_>>>()?;
            let ids = session
                .add_targets(&keys)
                .map_err(|e| Error::Reduction(format!("RustRed LiteRed DAG: {e}")))?;
            nodes.extend(new.into_iter().zip(ids));
        }
        context.cancellation.check()?;
        let ids: Vec<u32> = targets.iter().map(|t| nodes[t]).collect();
        let mut fixed: Vec<Option<Rational>> = vec![None; parameters.len()];
        if let (Some(position), Some(value)) = (*epsilon_position, epsilon) {
            fixed[position] = Some(value.clone());
        }
        let free: Vec<usize> = (0..fixed.len()).filter(|&i| fixed[i].is_none()).collect();
        let coefficients: Vec<Vec<(Atom, Option<Atom>)>> = match free.as_slice() {
            [] => {
                let point: Vec<Rational> = fixed.into_iter().map(|v| v.unwrap()).collect();
                session
                    .reduce_exact_at_point(&point, &ids)
                    .map_err(|e| Error::Reduction(format!("RustRed LiteRed reduction: {e}")))?
                    .coefficients
                    .into_iter()
                    .map(|row| row.into_iter().map(|v| (Atom::num(v), None)).collect())
                    .collect()
            }
            [variable] => {
                let x = &parameters[*variable];
                session
                    .reduce_univariate(&fixed, &ids)
                    .map_err(|e| Error::Reduction(format!("RustRed LiteRed reduction: {e}")))?
                    .coefficients
                    .iter()
                    .map(|row| row.iter().map(|f| rational_function(f, x)).collect())
                    .collect()
            }
            _ => {
                return Err(Error::Unsupported(
                    "the RustRed LiteRed backend reconstructs at most one free parameter; \
                     use sampled epsilon (reduce_at_epsilon)"
                        .into(),
                ));
            }
        };
        let mut out = Reduction {
            residuals: basis.clone(),
            ..Default::default()
        };
        let mut conditions = BTreeMap::new();
        for (target, row) in targets.iter().zip(coefficients) {
            if basis.contains(target) {
                continue;
            }
            let mut combination = LinearCombination::new();
            for (master, (coefficient, denominator)) in basis.iter().zip(row) {
                if coefficient.is_zero() {
                    continue;
                }
                if let Some(denominator) = denominator {
                    conditions.insert(denominator.to_canonical_string(), denominator);
                }
                combination.insert(master.clone(), coefficient);
            }
            out.rules.insert(target.clone(), combination);
        }
        out.nonzero_conditions = conditions.into_values().collect();
        Ok(out)
    }
}

fn key(integral: &Integral) -> Result<IntegralKey> {
    IntegralKey::try_new(integral.0.iter().map(|&n| i64::from(n)))
        .map_err(|e| Error::InvalidInput(e.to_string()))
}

fn integral(key: &IntegralKey) -> Result<Integral> {
    key.powers()
        .iter()
        .map(|&n| i16::try_from(n).map_err(|_| Error::Limit("index exceeds i16".into())))
        .collect::<Result<Vec<_>>>()
        .map(Integral)
}

fn polynomial(coefficients: &[Rational], x: &Atom) -> Atom {
    let mut value = Atom::new();
    for (degree, coefficient) in coefficients.iter().enumerate() {
        if coefficient.is_zero() {
            continue;
        }
        let term = Atom::num(coefficient.clone()) * x.pow(Atom::num(degree as i64));
        value = value + term;
    }
    value
}

/// Exact coefficient and, when nonconstant, its denominator polynomial.
fn rational_function(function: &RationalFunction, x: &Atom) -> (Atom, Option<Atom>) {
    if function.is_zero() {
        return (Atom::new(), None);
    }
    let denominator = polynomial(&function.denominator, x);
    let coefficient = (polynomial(&function.numerator, x) / &denominator)
        .together()
        .cancel();
    let constant = function.denominator.iter().skip(1).all(Rational::is_zero);
    (coefficient, (!constant).then_some(denominator))
}

impl ReductionBackend for LiteRedRustRedBackend {
    fn identity(&self) -> String {
        format!(
            "rustred-litered:{}:{}:rows={}:order={}:seed={}:isps-at-ends={}:arities={:?}",
            env!("RUSTRED_SOURCE_DIGEST"),
            env!("DEPENDENCY_SOURCE_DIGEST"),
            if self.options.group_orbit_rows {
                "orbit"
            } else {
                "generators"
            },
            if self.options.lines_dots_rank {
                "lines-dots-rank"
            } else {
                "program"
            },
            self.options.seed,
            self.options.isps_at_ends,
            rustred::compiled_runtime_arities(),
        )
    }

    fn reduce_at_epsilon(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        epsilon: &Rational,
        context: &RunContext,
    ) -> Result<Reduction> {
        self.reduce_with(family, targets, Some(epsilon), context)
    }

    fn reduce(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        context: &RunContext,
    ) -> Result<Reduction> {
        self.reduce_with(family, targets, None, context)
    }
}
