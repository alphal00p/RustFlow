//! Complete unit-call memoization. This never substitutes a descendant tail.
use super::*;
use crate::algebra::ExactAlgebraLimits;
use crate::persistence::{
    BinaryIoError, BinaryIoLimits, CoefficientId, CoefficientTableBuilder, DecodedCoefficientTable,
};
use std::mem::size_of;
use std::sync::Mutex;

/// Limits on aggregate retained memo data, separate from native reduction work.
/// `max_owned_bytes` counts the optional memo/lock field, fixed boxed entry slots, boxed label/ID
/// records and boxed native coefficient/state bytes. Allocator overhead, the
/// immutable rule graph, serialization/decoding scratch and returned reductions
/// are excluded. No decoded coefficient objects are retained in the memo.
#[derive(Clone, Copy, Debug)]
pub struct GuardedUnitMemoLimits {
    pub max_entries: usize,
    pub max_output_terms: usize,
    pub max_conditions: usize,
    pub max_polynomial_terms: usize,
    pub max_owned_bytes: usize,
}
impl Default for GuardedUnitMemoLimits {
    fn default() -> Self {
        Self {
            max_entries: 64,
            max_output_terms: 4096,
            max_conditions: 32768,
            max_polynomial_terms: 262144,
            max_owned_bytes: 64 * 1024 * 1024,
        }
    }
}
#[derive(Debug)]
pub struct GuardedMemoizedReduction<const N: usize> {
    pub reduction: GuardedReduction<N>,
    /// Logical application counts remain in `reduction`; a hit performs none.
    pub cache_hit: bool,
}
#[derive(Clone, Copy, Debug, Default)]
pub struct GuardedUnitMemoUsage {
    pub entries: usize,
    pub output_terms: usize,
    pub conditions: usize,
    pub polynomial_terms: usize,
    pub owned_bytes: usize,
}
#[derive(Debug)]
struct Record<const N: usize> {
    integral: [i64; N],
    coefficient: CoefficientId,
    terminal: bool,
}
#[derive(Debug)]
struct Entry<const N: usize> {
    target: [i64; N],
    limits: GuardedReductionLimits,
    records: Box<[Record<N>]>,
    conditions: Box<[CoefficientId]>,
    state: Box<[u8]>,
    atoms: Box<[u8]>,
    applications: usize,
    polynomial_terms: usize,
}
impl<const N: usize> Entry<N> {
    fn heap_bytes(&self) -> Option<usize> {
        self.records
            .len()
            .checked_mul(size_of::<Record<N>>())?
            .checked_add(
                self.conditions
                    .len()
                    .checked_mul(size_of::<CoefficientId>())?,
            )?
            .checked_add(self.state.len())?
            .checked_add(self.atoms.len())
    }
    fn decode(&self) -> Result<GuardedReduction<N>, SolverError> {
        // The only producer is this immutable program's successful native result.
        // Preserve its normalized coefficient representation; do not redo GCDs.
        let table =
            DecodedCoefficientTable::import_generated(&self.state, &self.atoms, transport_limits())
                .map_err(io_error)?;
        let mut result = GuardedReduction {
            terms: BTreeMap::new(),
            unresolved: Vec::new(),
            nonzero_conditions: Vec::new(),
            rule_applications: self.applications,
        };
        for record in &self.records {
            let coefficient = table
                .coefficient(record.coefficient)
                .map_err(io_error)?
                .clone();
            if record.terminal {
                if result.terms.insert(record.integral, coefficient).is_some() {
                    return Err(invalid("duplicate terminal in owned unit memo"));
                }
            } else {
                result.unresolved.push(GuardedUnresolvedTerm {
                    integral: record.integral,
                    coefficient,
                    reason: GuardedApplicationFailure::NoApplicableRule,
                });
            }
        }
        for &id in &self.conditions {
            let coefficient = table.coefficient(id).map_err(io_error)?;
            if !coefficient.denominator.is_one() {
                return Err(invalid("nonpolynomial condition in owned unit memo"));
            }
            result
                .nonzero_conditions
                .push(coefficient.numerator.clone());
        }
        Ok(result)
    }
}
fn invalid(message: &str) -> SolverError {
    SolverError::InvalidInput(message.into())
}
fn io_error(error: BinaryIoError) -> SolverError {
    invalid(&format!("owned native unit memo: {error}"))
}
// These are trusted, ephemeral, internally generated records. Quotas below bound
// retained bytes. They deliberately do not pretend to bound encoder/CAS scratch.
fn transport_limits() -> BinaryIoLimits {
    BinaryIoLimits {
        max_program_bytes: usize::MAX,
        max_sections: usize::MAX,
        max_collection_entries: usize::MAX,
        max_state_bytes: usize::MAX,
        max_atom_bytes: usize::MAX,
        max_total_atom_bytes: usize::MAX,
        exact_algebra: ExactAlgebraLimits {
            max_exponent: u16::MAX,
            max_polynomial_terms: usize::MAX,
            max_term_operations: usize::MAX,
        },
    }
}
fn same_limits(a: GuardedReductionLimits, b: GuardedReductionLimits) -> bool {
    a.max_rule_applications == b.max_rule_applications
        && a.max_pending_integrals == b.max_pending_integrals
}
#[derive(Debug)]
pub(super) struct UnitMemo<const N: usize> {
    limits: GuardedUnitMemoLimits,
    slots: Box<[Option<Entry<N>>]>,
}
impl<const N: usize> UnitMemo<N> {
    fn new(limits: GuardedUnitMemoLimits) -> Result<Self, SolverError> {
        let base = size_of::<Option<Mutex<Self>>>()
            .checked_add(
                limits
                    .max_entries
                    .checked_mul(size_of::<Option<Entry<N>>>())
                    .ok_or_else(|| invalid("unit memo metadata size overflow"))?,
            )
            .ok_or_else(|| invalid("unit memo metadata size overflow"))?;
        if base > limits.max_owned_bytes {
            return Err(invalid("unit memo metadata exceeds owned-byte limit"));
        }
        let mut slots = Vec::new();
        slots
            .try_reserve_exact(limits.max_entries)
            .map_err(|_| invalid("cannot allocate bounded unit memo slots"))?;
        slots.resize_with(limits.max_entries, || None);
        Ok(Self {
            limits,
            slots: slots.into_boxed_slice(),
        })
    }
    fn usage(&self) -> GuardedUnitMemoUsage {
        let mut u = GuardedUnitMemoUsage {
            owned_bytes: size_of::<Option<Mutex<Self>>>()
                + self.slots.len() * size_of::<Option<Entry<N>>>(),
            ..Default::default()
        };
        for e in self.slots.iter().flatten() {
            u.entries += 1;
            u.output_terms += e.records.len();
            u.conditions += e.conditions.len();
            u.polynomial_terms += e.polynomial_terms;
            u.owned_bytes += e.heap_bytes().expect("admitted owned entry size");
        }
        u
    }
    fn hit(
        &mut self,
        target: [i64; N],
        limits: GuardedReductionLimits,
    ) -> Result<Option<GuardedReduction<N>>, SolverError> {
        let Some(i) = self.slots.iter().position(|e| {
            e.as_ref()
                .is_some_and(|e| e.target == target && same_limits(e.limits, limits))
        }) else {
            return Ok(None);
        };
        let result = self.slots[i].as_ref().expect("found entry").decode()?;
        self.slots[i..].rotate_left(1);
        Ok(Some(result))
    }
    fn insert(
        &mut self,
        target: [i64; N],
        limits: GuardedReductionLimits,
        result: &GuardedReduction<N>,
    ) -> Result<(), SolverError> {
        if self.slots.is_empty()
            || result
                .unresolved
                .iter()
                .any(|r| r.reason != GuardedApplicationFailure::NoApplicableRule)
        {
            return Ok(());
        }
        let Some(outputs) = result.terms.len().checked_add(result.unresolved.len()) else {
            return Ok(());
        };
        if outputs > self.limits.max_output_terms
            || result.nonzero_conditions.len() > self.limits.max_conditions
        {
            return Ok(());
        }
        let mut polynomial_terms = 0usize;
        let mut minimum_integer_bytes = 0usize;
        let mut admit_poly = |p: &CoefficientPolynomial| -> bool {
            let Some(t) = polynomial_terms.checked_add(p.nterms()) else {
                return false;
            };
            polynomial_terms = t;
            for c in &p.coefficients {
                let Some(b) = c
                    .significant_bits()
                    .checked_add(7)
                    .and_then(|v| usize::try_from(v / 8).ok())
                    .and_then(|v| minimum_integer_bytes.checked_add(v))
                else {
                    return false;
                };
                minimum_integer_bytes = b;
            }
            polynomial_terms <= self.limits.max_polynomial_terms
                && minimum_integer_bytes <= self.limits.max_owned_bytes
        };
        for c in result
            .terms
            .values()
            .chain(result.unresolved.iter().map(|r| &r.coefficient))
        {
            if !admit_poly(&c.numerator) || !admit_poly(&c.denominator) {
                return Ok(());
            }
        }
        for p in &result.nonzero_conditions {
            if !admit_poly(p) {
                return Ok(());
            }
        }
        let mut table = CoefficientTableBuilder::new(transport_limits());
        let mut records = Vec::new();
        for (integral, coefficient) in &result.terms {
            records.push(Record {
                integral: *integral,
                coefficient: table.intern(coefficient).map_err(io_error)?,
                terminal: true,
            });
        }
        for r in &result.unresolved {
            records.push(Record {
                integral: r.integral,
                coefficient: table.intern(&r.coefficient).map_err(io_error)?,
                terminal: false,
            });
        }
        let conditions = result
            .nonzero_conditions
            .iter()
            .map(|p| table.intern(&p.clone().into()).map_err(io_error))
            .collect::<Result<Vec<_>, _>>()?;
        let encoded = table.finish().map_err(io_error)?;
        let entry = Entry {
            target,
            limits,
            records: records.into_boxed_slice(),
            conditions: conditions.into_boxed_slice(),
            state: encoded.state.into_boxed_slice(),
            atoms: encoded.atoms.into_boxed_slice(),
            applications: result.rule_applications,
            polynomial_terms,
        };
        let Some(bytes) = entry.heap_bytes() else {
            return Ok(());
        };
        let base =
            size_of::<Option<Mutex<Self>>>() + self.slots.len() * size_of::<Option<Entry<N>>>();
        if bytes > self.limits.max_owned_bytes - base {
            return Ok(());
        }
        loop {
            let u = self.usage();
            if u.entries < self.slots.len()
                && u.output_terms <= self.limits.max_output_terms - outputs
                && u.conditions <= self.limits.max_conditions - entry.conditions.len()
                && u.polynomial_terms <= self.limits.max_polynomial_terms - polynomial_terms
                && u.owned_bytes <= self.limits.max_owned_bytes - bytes
            {
                break;
            }
            let i = self
                .slots
                .iter()
                .position(Option::is_some)
                .expect("one entry can fit after eviction");
            self.slots[i] = None;
        }
        let i = self
            .slots
            .iter()
            .position(Option::is_none)
            .expect("admitted slot");
        self.slots[i..].rotate_left(1);
        *self.slots.last_mut().expect("nonempty slots") = Some(entry);
        Ok(())
    }
}
impl<const N: usize> GuardedProgram<N> {
    /// Start a new serial whole-call memo scope. No entry survives this reset.
    /// Rule/source/order data remain immutable; plain `reduce` never uses it.
    pub fn with_unit_reduction_memo(
        mut self,
        limits: Option<GuardedUnitMemoLimits>,
    ) -> Result<Self, SolverError> {
        self.unit_memo = limits.map(UnitMemo::new).transpose()?.map(Mutex::new);
        Ok(self)
    }
    pub fn unit_reduction_memo_usage(&self) -> Result<Option<GuardedUnitMemoUsage>, SolverError> {
        self.unit_memo
            .as_ref()
            .map(|m| {
                m.lock()
                    .map(|m| m.usage())
                    .map_err(|_| invalid("native unit memo lock poisoned"))
            })
            .transpose()
    }
    pub fn reduce_memoized(
        &self,
        target: [i64; N],
        limits: GuardedReductionLimits,
    ) -> Result<GuardedMemoizedReduction<N>, SolverError> {
        let Some(memo) = &self.unit_memo else {
            return Ok(GuardedMemoizedReduction {
                reduction: self.reduce(target, limits)?,
                cache_hit: false,
            });
        };
        // Serialize cache calls; reductions never recurse through this method.
        let mut memo = memo
            .lock()
            .map_err(|_| invalid("native unit memo lock poisoned"))?;
        if let Some(reduction) = memo.hit(target, limits)? {
            return Ok(GuardedMemoizedReduction {
                reduction,
                cache_hit: true,
            });
        }
        let reduction = self.reduce(target, limits)?;
        memo.insert(target, limits, &reduction)?;
        Ok(GuardedMemoizedReduction {
            reduction,
            cache_hit: false,
        })
    }
    /// Add only old NoRule stopping points, replay every source rule, and retain
    /// complete unit-call maps by the exact NoRule-to-terminal equivalence.
    /// Ordinary unions/rebindings reset memo state instead.
    pub fn with_promoted_terminals_replayed(
        mut self,
        current_terminals: impl IntoIterator<Item = [i64; N]>,
        max_rules: usize,
    ) -> Result<Self, SolverError> {
        let terminals = Self::validate_terminals(&self.sources, current_terminals)?;
        if !self.terminals.is_subset(&terminals) {
            return Err(invalid(
                "unit memo terminal promotion cannot remove old terminals",
            ));
        }
        let added = terminals
            .difference(&self.terminals)
            .copied()
            .collect::<BTreeSet<_>>();
        for label in &added {
            let r = self.apply(label)?;
            if !matches!(
                r.status,
                GuardedApplicationStatus::Unresolved(GuardedApplicationFailure::NoApplicableRule)
            ) || !r.terms.is_empty()
                || !r.nonzero_conditions.is_empty()
            {
                return Err(invalid(
                    "unit memo promotion requires exact old NoApplicableRule labels",
                ));
            }
        }
        let memo = self.unit_memo.take();
        let mut result = self.with_terminals_replayed(terminals, max_rules)?;
        if let Some(memo) = memo {
            let mut memo = memo
                .into_inner()
                .map_err(|_| invalid("native unit memo lock poisoned"))?;
            for entry in memo.slots.iter_mut().flatten() {
                for record in &mut entry.records {
                    if !record.terminal && added.contains(&record.integral) {
                        record.terminal = true;
                    }
                }
            }
            result.unit_memo = Some(Mutex::new(memo));
        }
        Ok(result)
    }
}
#[cfg(test)]
mod tests;
