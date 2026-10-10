//! Private class-preservation checker; physical identity validity belongs to the
//! typed source/origin factory. No per-application trace or uniform-D claim.
use super::guarded::{GuardedContext, IndexBounds, IndexDomain, IndexRole};
use super::preparation::WeightedSourceOptions;
use crate::{Error, Result, RunContext};
use rustred::algebra::CoefficientPolynomial;
use rustred::persistence::BinaryIoLimits;
use rustred::solver::Integral as NativeIntegral;
use rustred::solver::guarded::{GuardedProgram, GuardedSourceSystem};
use serde::Serialize;
use serde_json::{Value, json};
use std::collections::BTreeMap;
use std::sync::Arc;
use symbolica::prelude::*;
const VERSION: &str =
    "partial-source-image-class-v1;integer-boundary-zero;universal-finite-label;no-applied-trace";
fn invalid(s: impl Into<String>) -> Error {
    Error::InvalidInput(format!("source image class: {}", s.into()))
}
fn limit(s: &str) -> Error {
    Error::Limit(format!("source image class: {s}"))
}
fn unsupported(s: &str) -> Error {
    Error::Unsupported(format!("source image class: {s}"))
}

/// Constructed only by the bound partial bundle, after its checked physical
/// family/placement/options/origin validation. This payload grants no theorem.
#[derive(Clone, Debug)]
pub(crate) struct SourceImageBinding {
    pub(crate) input_identity: String,
    pub(crate) family_signature: Vec<String>,
    pub(crate) origin_identity: String,
    pub(crate) shifted_slots: Vec<usize>,
    pub(crate) source_options: WeightedSourceOptions,
    pub(crate) expected_roles: Vec<IndexRole>,
    pub(crate) physical_arity: usize,
    pub(crate) physical_slots: usize,
    pub(crate) input_slots: usize,
}
impl SourceImageBinding {
    fn report(&self) -> Value {
        json!({
            "input_identity":self.input_identity,"family_signature":self.family_signature,
            "origin_identity":self.origin_identity,"shifted_slots":self.shifted_slots,
            "source_options":self.source_options,
            "roles":self.expected_roles.iter().map(|r|format!("{r:?}")).collect::<Vec<_>>(),
            "physical_arity":self.physical_arity,"physical_slots":self.physical_slots,"input_slots":self.input_slots,
        })
    }
    fn digest(&self, budget: SourceImageBudget) -> Result<String> {
        let mut bytes = self
            .input_identity
            .len()
            .checked_add(self.origin_identity.len())
            .ok_or_else(|| limit("binding size overflow"))?;
        if self.family_signature.len() > budget.max_binding_fields
            || self.expected_roles.len() > budget.max_arity
            || self.shifted_slots.len() > budget.max_arity
        {
            return Err(limit("binding field count"));
        }
        for s in &self.family_signature {
            bytes = bytes
                .checked_add(s.len())
                .ok_or_else(|| limit("binding size overflow"))?;
        }
        if bytes > budget.max_binding_bytes {
            return Err(limit("binding text bytes"));
        }
        if self.input_identity.is_empty()
            || self.origin_identity.is_empty()
            || self.family_signature.is_empty()
        {
            return Err(invalid("missing typed physical binding"));
        }
        let encoded =
            serde_json::to_vec(&self.report()).map_err(|e| Error::Cache(e.to_string()))?;
        if encoded.len() > budget.max_binding_bytes {
            return Err(limit("encoded binding bytes"));
        }
        Ok(blake3::hash(&encoded).to_string())
    }
}
#[derive(Clone, Copy, Debug, Serialize)]
pub(crate) struct SourceImageBudget {
    pub(crate) max_arity: usize,
    pub(crate) max_sources: usize,
    pub(crate) max_terms: usize,
    pub(crate) max_conditions: usize,
    pub(crate) max_polynomial_terms: usize,
    pub(crate) max_polynomial_degree: u32,
    pub(crate) max_integer_bits: u32,
    pub(crate) max_face_specializations: usize,
    pub(crate) max_face_abs: u64,
    pub(crate) max_binding_fields: usize,
    pub(crate) max_binding_bytes: usize,
    pub(crate) max_source_bytes: usize,
}
impl Default for SourceImageBudget {
    fn default() -> Self {
        Self {
            max_arity: 32,
            max_sources: 16384,
            max_terms: 1_000_000,
            max_conditions: 1_000_000,
            max_polynomial_terms: 2_000_000,
            max_polynomial_degree: 128,
            max_integer_bits: 4096,
            max_face_specializations: 1_000_000,
            max_face_abs: 65536,
            max_binding_fields: 16384,
            max_binding_bytes: 4 * 1024 * 1024,
            max_source_bytes: 64 * 1024 * 1024,
        }
    }
}
#[derive(Clone, Debug, Serialize)]
struct Counters {
    sources: usize,
    raw_terms: usize,
    merged_terms: usize,
    source_conditions: usize,
    polynomial_terms: usize,
    checked_faces: usize,
    source_bytes: usize,
}
impl Counters {
    fn new() -> Self {
        Self {
            sources: 0,
            raw_terms: 0,
            merged_terms: 0,
            source_conditions: 0,
            polynomial_terms: 0,
            checked_faces: 0,
            source_bytes: 0,
        }
    }
}
/// Immutable evidence about the exact source corpus, not a new source theorem.
#[derive(Clone, Debug)]
pub(crate) struct SourceImageClassCertificate {
    binding_digest: String,
    source_digest: String,
    admitted: Vec<(Option<i64>, Option<i64>)>,
    physical_arity: usize,
    storage_capacity: usize,
    budget: SourceImageBudget,
    counts: Counters,
}
fn bounds<const N: usize>(a: &IndexDomain<N>) -> Vec<(Option<i64>, Option<i64>)> {
    a.bounds().iter().map(|b| (b.lower(), b.upper())).collect()
}
impl SourceImageClassCertificate {
    pub(crate) fn report(&self) -> Value {
        json!({"version":VERSION,"binding_blake3":self.binding_digest,"source_only_program_blake3":self.source_digest,"admitted_bounds":self.admitted,"physical_arity":self.physical_arity,"storage_capacity":self.storage_capacity,"budget":self.budget,"counts":self.counts,"coverage":"all nonzero source images pointwise on admitted integer guards; native seed/recenter/replay preserves this class","individual_application_trace":false,"one_dimension_for_unbounded_domain":false,"physical_identity_validity_certified_by_this_checker":false,"internal_cas_work_or_heap_bound":false})
    }
    pub(crate) fn validate_context<const N: usize>(
        &self,
        context: &GuardedContext<N>,
        admitted: &IndexDomain<N>,
        binding: &SourceImageBinding,
        run: &RunContext,
    ) -> Result<()> {
        run.cancellation.check()?;
        validate_layout(
            context.sources(),
            context.physical_arity(),
            admitted,
            binding,
            self.budget,
        )?;
        if self.physical_arity != context.physical_arity()
            || self.storage_capacity != N
            || self.admitted != bounds(admitted)
            || self.binding_digest != binding.digest(self.budget)?
        {
            return Err(invalid("certificate physical/domain binding mismatch"));
        }
        let (digest, _) = source_digest(context.sources(), self.budget, run)?;
        if digest != self.source_digest {
            return Err(invalid("certificate exact source context mismatch"));
        }
        run.cancellation.check()
    }
}
fn validate_layout<const N: usize>(
    sources: &Arc<GuardedSourceSystem<N>>,
    physical_arity: usize,
    admitted: &IndexDomain<N>,
    binding: &SourceImageBinding,
    budget: SourceImageBudget,
) -> Result<()> {
    if N == 0
        || N > budget.max_arity
        || physical_arity == 0
        || physical_arity > N
        || physical_arity != binding.physical_arity
        || binding.physical_slots > binding.input_slots
        || binding.input_slots > physical_arity
        || binding.expected_roles.len() != physical_arity
        || sources.roles()[..physical_arity] != binding.expected_roles
    {
        return Err(invalid("physical arity/role/layout mismatch"));
    }
    if binding.source_options.positive_compact_energy_powers {
        return Err(unsupported(
            "partial source class excludes inverse completions",
        ));
    }
    let mut expected = [IndexBounds::unbounded(); N];
    for i in 0..N {
        let role = sources.roles()[i];
        expected[i] = if i >= physical_arity {
            if role != IndexRole::Ordinary {
                return Err(invalid("nonordinary tail role"));
            }
            IndexBounds::fixed(0)
        } else if i >= binding.input_slots {
            if role != IndexRole::Occupation {
                return Err(invalid("nonoccupation endpoint role"));
            }
            IndexBounds::new(Some(0), None).map_err(|e| invalid(e.to_string()))?
        } else if i >= binding.physical_slots {
            if role != IndexRole::Ordinary {
                return Err(invalid("nonordinary completion role"));
            }
            IndexBounds::new(None, Some(0)).map_err(|e| invalid(e.to_string()))?
        } else {
            if role == IndexRole::Occupation {
                return Err(invalid("occupation in physical factor range"));
            }
            IndexBounds::unbounded()
        };
    }
    let a = IndexDomain::new(expected).map_err(|e| invalid(e.to_string()))?;
    if admitted != &a {
        return Err(invalid(
            "admitted box differs from polynomial occupied class",
        ));
    }
    let mut shifts = binding.shifted_slots.clone();
    shifts.sort_unstable();
    shifts.dedup();
    if shifts.is_empty()
        || shifts.len() != binding.shifted_slots.len()
        || shifts
            .iter()
            .any(|&s| s >= binding.physical_slots || sources.roles()[s] != IndexRole::Ordinary)
    {
        return Err(invalid("invalid fixed-shell placement"));
    }
    if sources
        .zero_domains()
        .iter()
        .any(|d| !d.is_subset_of(admitted))
    {
        return Err(invalid("zero box escapes admitted tuple class"));
    }
    Ok(())
}
fn inspect_poly(p: &CoefficientPolynomial, b: SourceImageBudget, count: &mut usize) -> Result<()> {
    *count = count
        .checked_add(p.nterms())
        .ok_or_else(|| limit("polynomial term counter overflow"))?;
    if *count > b.max_polynomial_terms {
        return Err(limit("polynomial term budget"));
    }
    for exponents in p.exponents_iter() {
        let degree = exponents
            .iter()
            .try_fold(0u32, |sum, &n| sum.checked_add(u32::from(n)))
            .ok_or_else(|| limit("polynomial degree overflow"))?;
        if degree > b.max_polynomial_degree {
            return Err(limit("polynomial degree budget"));
        }
    }
    if p.coefficients
        .iter()
        .any(|c| c.significant_bits() > u64::from(b.max_integer_bits))
    {
        return Err(limit("integer coefficient size"));
    }
    Ok(())
}
fn source_digest<const N: usize>(
    sources: &Arc<GuardedSourceSystem<N>>,
    budget: SourceImageBudget,
    run: &RunContext,
) -> Result<(String, usize)> {
    run.cancellation.check()?;
    let program = GuardedProgram::new(sources.clone(), vec![], [])
        .map_err(|e| Error::Reduction(e.to_string()))?;
    let bytes = program
        .encode_native(BinaryIoLimits {
            max_program_bytes: budget.max_source_bytes,
            max_state_bytes: budget.max_source_bytes,
            max_total_atom_bytes: budget.max_source_bytes,
            max_atom_bytes: budget.max_source_bytes,
            ..Default::default()
        })
        .map_err(|e| Error::Cache(e.to_string()))?;
    run.cancellation.check()?;
    Ok((blake3::hash(&bytes).to_string(), bytes.len()))
}
/// Bundle caller must validate its sealed physical/origin binding before this
/// function; this function certifies only the exact emitted index class.
pub(crate) fn certify_source_image_class<const N: usize>(
    context: &GuardedContext<N>,
    admitted: &IndexDomain<N>,
    binding: &SourceImageBinding,
    budget: SourceImageBudget,
    run: &RunContext,
) -> Result<SourceImageClassCertificate> {
    certify_native(
        context.sources(),
        context.physical_arity(),
        admitted,
        binding,
        budget,
        run,
    )
}
fn certify_native<const N: usize>(
    sources: &Arc<GuardedSourceSystem<N>>,
    physical_arity: usize,
    admitted: &IndexDomain<N>,
    binding: &SourceImageBinding,
    budget: SourceImageBudget,
    run: &RunContext,
) -> Result<SourceImageClassCertificate> {
    run.cancellation.check()?;
    validate_layout(sources, physical_arity, admitted, binding, budget)?;
    let binding_digest = binding.digest(budget)?;
    let system = sources.native_sources();
    if system.fixed().iter().any(Option::is_some) {
        return Err(unsupported(
            "prepared fixed raw-source coordinates are outside factory pattern",
        ));
    }
    if sources.sources().len() > budget.max_sources {
        return Err(limit("source row budget"));
    }
    if sources.sources().len() != system.rows().len() {
        return Err(invalid("source metadata/row mismatch"));
    }
    let mut counts = Counters::new();
    // Preflight complete corpus before any coefficient combination or encoding.
    for (info, row) in sources.sources().iter().zip(system.rows()) {
        run.cancellation.check()?;
        for d in info.domain.bounds() {
            if let (Some(lo), Some(hi)) = (d.lower(), d.upper()) {
                if lo == hi && lo.unsigned_abs() > budget.max_face_abs {
                    return Err(limit("fixed guard coordinate magnitude"));
                }
            }
        }
        if !info.domain.is_subset_of(admitted) {
            return Err(invalid(format!(
                "source {} base escapes admitted class",
                info.id
            )));
        }
        counts.raw_terms = counts
            .raw_terms
            .checked_add(row.len())
            .ok_or_else(|| limit("source term counter overflow"))?;
        counts.source_conditions = counts
            .source_conditions
            .checked_add(info.nonzero_conditions.len())
            .ok_or_else(|| limit("condition counter overflow"))?;
        if counts.raw_terms > budget.max_terms || counts.source_conditions > budget.max_conditions {
            return Err(limit("source term/condition budget"));
        }
        for term in row {
            if term.integral.powers().iter().any(|p| !p.is_symbolic()) {
                return Err(unsupported("nonsymbolic raw source pattern"));
            }
            inspect_poly(&term.coefficient, budget, &mut counts.polynomial_terms)?;
        }
        for p in &info.nonzero_conditions {
            inspect_poly(p, budget, &mut counts.polynomial_terms)?;
        }
    }
    for (info, row) in sources.sources().iter().zip(system.rows()) {
        run.cancellation.check()?;
        counts.sources += 1;
        let mut merged = BTreeMap::<NativeIntegral<N>, CoefficientPolynomial>::new();
        for t in row {
            run.cancellation.check()?;
            match merged.entry(t.integral) {
                std::collections::btree_map::Entry::Vacant(e) => {
                    e.insert(t.coefficient.clone());
                }
                std::collections::btree_map::Entry::Occupied(mut e) => {
                    let sum = e.get() + &t.coefficient;
                    *e.get_mut() = sum;
                }
            }
        }
        for (image, mut coefficient) in merged {
            run.cancellation.check()?;
            // Numeric guard faces substitute all fixed coordinates, even when
            // a different axis is the potentially invalid image coordinate.
            for (axis, domain) in info.domain.bounds().iter().enumerate() {
                if let (Some(lo), Some(hi)) = (domain.lower(), domain.upper()) {
                    if lo == hi {
                        coefficient =
                            coefficient.replace(system.index_variables()[axis], &Integer::from(lo));
                        run.cancellation.check()?;
                    }
                }
            }
            let mut retained_terms = 0;
            inspect_poly(&coefficient, budget, &mut retained_terms)?;
            if coefficient.is_zero() {
                continue;
            }
            counts.merged_terms += 1;
            for axis in 0..N {
                let shift = i64::from(image[axis].value());
                let g = info.domain.bounds()[axis];
                let a = admitted.bounds()[axis];
                let mut violations = Vec::new();
                if let Some(lo) = a.lower() {
                    let end = lo
                        .checked_sub(shift)
                        .and_then(|x| x.checked_sub(1))
                        .ok_or_else(|| limit("lower image boundary overflow"))?;
                    if g.lower().is_none_or(|x| x <= end) {
                        let hi = Some(g.upper().map_or(end, |x| x.min(end)));
                        violations.push((g.lower(), hi));
                    }
                }
                if let Some(hi) = a.upper() {
                    let start = hi
                        .checked_sub(shift)
                        .and_then(|x| x.checked_add(1))
                        .ok_or_else(|| limit("upper image boundary overflow"))?;
                    if g.upper().is_none_or(|x| x >= start) {
                        let lo = Some(g.lower().map_or(start, |x| x.max(start)));
                        violations.push((lo, g.upper()));
                    }
                }
                for (lo, hi) in violations {
                    let (Some(lo), Some(hi)) = (lo, hi) else {
                        return Err(unsupported("unbounded violating image faces"));
                    };
                    if lo > hi {
                        continue;
                    }
                    if lo.unsigned_abs() > budget.max_face_abs
                        || hi.unsigned_abs() > budget.max_face_abs
                    {
                        return Err(limit("face index magnitude"));
                    }
                    let count = usize::try_from(i128::from(hi) - i128::from(lo) + 1)
                        .map_err(|_| limit("face count overflow"))?;
                    let total = counts
                        .checked_faces
                        .checked_add(count)
                        .ok_or_else(|| limit("face counter overflow"))?;
                    if total > budget.max_face_specializations {
                        return Err(limit("face specialization budget"));
                    }
                    for value in lo..=hi {
                        run.cancellation.check()?;
                        counts.checked_faces += 1;
                        let specialized = coefficient
                            .replace(system.index_variables()[axis], &Integer::from(value));
                        run.cancellation.check()?;
                        let mut retained_terms = 0;
                        inspect_poly(&specialized, budget, &mut retained_terms)?;
                        if !specialized.is_zero() {
                            return Err(invalid(format!(
                                "source {} has nonzero inadmissible image at axis {axis}, base {value}, shift {shift}",
                                info.id
                            )));
                        }
                    }
                }
            }
        }
    }
    let (source_digest, bytes) = source_digest(sources, budget, run)?;
    counts.source_bytes = bytes;
    run.cancellation.check()?;
    Ok(SourceImageClassCertificate {
        binding_digest,
        source_digest,
        admitted: bounds(admitted),
        physical_arity,
        storage_capacity: N,
        budget,
        counts,
    })
}
#[cfg(test)]
#[path = "source_class/tests.rs"]
mod tests;
