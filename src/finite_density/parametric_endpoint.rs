//! Query-only exact positive-parametric endpoint evidence.
//!
//! HEPKit owns the Symanzik construction. This owner resolves its positive
//! polynomials with complete coordinate blowups, retaining every chart and
//! original parameter map. It does not admit a thermal contour, evaluate a
//! period, or extend `MasslessFlowEvidence`.
use super::{PreparedDensityInput, geometry::OccupiedCutFamily};
use crate::algebra::{matmul, rref};
use crate::coefficient::{exact_coefficient_list, powers};
use crate::{Error, Result};
use feynkit_graph::IntegralFamily as NativeFamily;
use feynkit_kinematics::Kinematics;
use serde::Serialize;
use std::collections::BTreeSet;
use symbolica::prelude::*;

const VERSION: &str = "hepkit-positive-joint-eta-sectors-v2;full-off-null-identity;query-only;UV-meromorphic-before-high-D";

#[derive(Clone, Copy, Debug)]
pub struct ParametricEndpointBudget {
    pub parameters: usize,
    pub primary_charts: usize,
    pub sectors: usize,
    pub tree_nodes: usize,
    pub depth: usize,
    pub polynomial_terms: usize,
    pub exponent: u32,
    pub supports: usize,
    /// Bounds conservative integer-map and center-search work, not wall time.
    pub operations: usize,
}
impl Default for ParametricEndpointBudget {
    fn default() -> Self {
        Self {
            parameters: 8,
            primary_charts: 40320,
            sectors: 50000,
            tree_nodes: 200000,
            depth: 32,
            polynomial_terms: 20000,
            exponent: 1_000_000,
            supports: 256,
            operations: 200_000_000,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct PositiveMonomial {
    pub exponents: Vec<u32>,
    pub coefficient: String,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct EtaSectorLeaf {
    pub u: Vec<u32>,
    pub f: Vec<u32>,
    /// Pullback of d alpha_1 ... d alpha_(n-1) d eta/eta is z^(j-1) dz.
    pub jacobian_plus_one: Vec<u32>,
    pub eta: Vec<u32>,
    pub original_parameters: Vec<Vec<u32>>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub enum EtaSectorNode {
    Leaf(EtaSectorLeaf),
    /// The children contain exactly one chart per pivot in `axes`.
    Blowup {
        axes: Vec<usize>,
        children: Vec<(usize, Box<EtaSectorNode>)>,
    },
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct EtaPrimaryChart {
    /// alpha[order[0]]=1; alpha[order[k]]=z_0 ... z_(k-1).
    pub order: Vec<usize>,
    pub tree: EtaSectorNode,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct ParametricEtaBound {
    pub chart: usize,
    pub leaf: usize,
    pub axis: usize,
    pub eta_multiplicity: u32,
    /// Twice the D coefficient of the Mellin numerator, before dividing by m.
    pub twice_dimension_coefficient: i64,
    pub constant: i64,
    pub regular: bool,
}

/// Sealed construction, serializable for audit but deliberately not loadable as
/// an admission token. `replay` checks every retained chart against the native
/// polynomials. Input and exact family signatures participate in `identity`.
#[derive(Clone, Debug, Serialize)]
pub struct ParametricSupportEvidence {
    version: &'static str,
    input_identity: String,
    family_signature: Vec<String>,
    cuts: Vec<usize>,
    active_slots: Vec<usize>,
    routed_rows: Vec<Vec<String>>,
    virtual_loops: usize,
    active_virtual_rank: usize,
    pure_transfer_slots: Vec<(usize, String)>,
    u: Vec<PositiveMonomial>,
    v: Vec<PositiveMonomial>,
    charts: Vec<EtaPrimaryChart>,
    resolution_operations: usize,
    full_kinematics: Option<FullKinematicSymanzik>,
}
impl ParametricSupportEvidence {
    pub fn active_slots(&self) -> &[usize] {
        &self.active_slots
    }
    pub fn active_virtual_rank(&self) -> usize {
        self.active_virtual_rank
    }
    pub fn virtual_loops(&self) -> usize {
        self.virtual_loops
    }
    pub fn is_unrestricted_virtual_polynomial(&self) -> bool {
        self.active_virtual_rank < self.virtual_loops
    }
    pub fn u(&self) -> &[PositiveMonomial] {
        &self.u
    }
    pub fn v(&self) -> &[PositiveMonomial] {
        &self.v
    }
    pub fn charts(&self) -> &[EtaPrimaryChart] {
        &self.charts
    }
    pub fn pure_transfer_slots(&self) -> &[(usize, String)] {
        &self.pure_transfer_slots
    }
    pub fn full_kinematics(&self) -> Option<&FullKinematicSymanzik> {
        self.full_kinematics.as_ref()
    }
    pub fn label_envelope(
        &self,
        family: &OccupiedCutFamily,
        label: &crate::Integral,
    ) -> Result<ParametricLabelEnvelope> {
        if super::massless_endpoint::bound_family_signature(family) != self.family_signature {
            return Err(invalid("label envelope family binding mismatch"));
        }
        parametric_label_envelope(
            family,
            label,
            &self
                .pure_transfer_slots
                .iter()
                .map(|(s, _)| *s)
                .collect::<Vec<_>>(),
        )
    }
    /// Full finite-label query. Required-cut and lower-contact zeros still
    /// belong to their existing distribution owner and are not inferred here.
    pub fn bounds_for_label(
        &self,
        family: &OccupiedCutFamily,
        label: &crate::Integral,
    ) -> Result<ParametricLabelBoundQuery> {
        let envelope = self.label_envelope(family, label)?;
        validate_envelope_support(&envelope, &self.active_slots)?;
        let eta_bounds = self.bounds_for_polynomial_rank(
            &envelope.positive_virtual_powers,
            envelope.polynomial_momentum_rank,
            sum_u32(&envelope.shell_and_upper_jets)?,
        )?;
        let p = sum_i64(&envelope.positive_virtual_powers)?;
        let p0 = envelope
            .positive_transfer_powers
            .iter()
            .try_fold(0_i64, |s, (_, n)| checked_add(s, i64::from(*n)))?;
        let j = sum_i64(&envelope.shell_and_upper_jets)?;
        let mut compact_bounds = vec![CompactEndpointDegree {
            stratum: "two-vector angular measure".into(),
            twice_dimension_coefficient: 1,
            constant: -1,
        }];
        if !self.is_unrestricted_virtual_polynomial() {
            compact_bounds.push(CompactEndpointDegree {
                stratum: "combined virtual and deformed transfer endpoint jets".into(),
                twice_dimension_coefficient: self.virtual_loops as i64,
                constant: -checked_add(checked_add(p, p0)?, j)?,
            });
        }
        for (i, &ji) in envelope.shell_and_upper_jets.iter().enumerate() {
            compact_bounds.push(CompactEndpointDegree {
                stratum: format!("compact radial origin {i}"),
                twice_dimension_coefficient: 2,
                constant: -checked_add(2, checked_mul(2, i64::from(ji))?)?,
            });
        }
        Ok(ParametricLabelBoundQuery {
            envelope,
            eta_bounds,
            compact_bounds,
        })
    }
    pub fn identity(&self) -> Result<String> {
        let bytes = serde_json::to_vec(self).map_err(|e| Error::Cache(e.to_string()))?;
        Ok(format!("{VERSION}:{}", blake3::hash(&bytes).to_hex()))
    }
    /// Reconstruct the native HEPKit polynomials from this exact input/family
    /// and compare the whole construction, including the complete tree.
    pub fn replay_for(
        &self,
        input: &PreparedDensityInput,
        family: &OccupiedCutFamily,
        budget: ParametricEndpointBudget,
    ) -> Result<()> {
        let fresh = two_occupied_support_evidence(input, family, &self.active_slots, budget)?;
        if fresh.identity()? != self.identity()? {
            return Err(invalid("native family-bound evidence does not replay"));
        }
        self.replay(budget)
    }
    /// Replay all polynomial maps and every pivot branch. The native U/V
    /// coefficients are also rechecked as positive exact rationals.
    pub fn replay(&self, budget: ParametricEndpointBudget) -> Result<()> {
        validate_budget(budget)?;
        if self.is_unrestricted_virtual_polynomial() {
            if !self.u.is_empty()
                || !self.v.is_empty()
                || !self.charts.is_empty()
                || self.full_kinematics.is_some()
            {
                return Err(invalid("rank-deficient evidence carries parametric charts"));
            }
            return Ok(());
        }
        let full = self
            .full_kinematics
            .as_ref()
            .ok_or_else(|| invalid("missing full off-null Symanzik evidence"))?;
        if full.compact_loops != 2
            || full.virtual_loops != self.virtual_loops
            || full.parameter_slots != self.active_slots
            || full.u != self.u
            || full.pairs.len() != 1
            || full.pairs[0].pair != [0, 1]
            || full.pairs[0].terms != self.v
            || full.diagonal_masses.len() != 2
            || full.virtual_masses.len() != self.active_slots.len()
        {
            return Err(invalid(
                "off-null evidence differs from the retained null sector proof",
            ));
        }
        let n = self.active_slots.len();
        let orders = permutations(n, budget)?;
        if self.charts.len() != orders.len() {
            return Err(invalid("missing primary charts"));
        }
        let (u, f) = supports(&self.u, &self.v, n, budget)?;
        let mut counters = Counters::default();
        for (chart, order) in self.charts.iter().zip(orders) {
            if chart.order != order {
                return Err(invalid("primary order changed"));
            }
            let state = State::primary(&u, &f, &order, budget)?;
            replay_node(
                &chart.tree,
                &state,
                self.virtual_loops,
                0,
                budget,
                &mut counters,
            )?;
        }
        Ok(())
    }
    /// Sufficient Mellin inequalities for finite Gaussian insertions. The
    /// supplied rank and jets are explicit query inputs, not inferred physical
    /// admission. P includes the positive virtual powers only. A J-fold mass
    /// jet is charged once; an already differentiated representation uses J=0.
    pub fn bounds_for_polynomial_rank(
        &self,
        positive_powers: &[u32],
        polynomial_rank: u32,
        mass_jets: u32,
    ) -> Result<Vec<ParametricEtaBound>> {
        if positive_powers.len() != self.active_slots.len() || positive_powers.contains(&0) {
            return Err(invalid(
                "one strictly positive index is required per active slot",
            ));
        }
        if self.is_unrestricted_virtual_polynomial() {
            return Ok(Vec::new());
        }
        let p = positive_powers
            .iter()
            .try_fold(0_i64, |s, &n| s.checked_add(i64::from(n)))
            .ok_or_else(|| limit("positive index sum overflow"))?;
        let mut output = Vec::new();
        for (chart_index, chart) in self.charts.iter().enumerate() {
            let mut leaves = Vec::new();
            collect_leaves(&chart.tree, &mut leaves);
            for (leaf_index, leaf) in leaves.into_iter().enumerate() {
                for (axis, &m) in leaf.eta.iter().enumerate().filter(|(_, m)| **m > 0) {
                    let u = i64::from(leaf.u[axis]);
                    let f = i64::from(leaf.f[axis]);
                    let h = i64::try_from(self.virtual_loops)
                        .map_err(|_| limit("loop rank overflow"))?;
                    let d = checked_sub(checked_mul(h, f)?, checked_mul(h + 1, u)?)?;
                    let mut c = i64::from(leaf.jacobian_plus_one[axis]);
                    c = checked_add(c, checked_mul(p - i64::from(polynomial_rank), u)?)?;
                    c = checked_sub(c, checked_mul(checked_add(p, i64::from(mass_jets))?, f)?)?;
                    for (&power, map) in positive_powers.iter().zip(&leaf.original_parameters) {
                        c = checked_add(
                            c,
                            checked_mul(i64::from(power - 1), i64::from(map[axis]))?,
                        )?;
                    }
                    output.push(ParametricEtaBound {
                        chart: chart_index,
                        leaf: leaf_index,
                        axis,
                        eta_multiplicity: m,
                        twice_dimension_coefficient: d,
                        constant: c,
                        regular: d == 0,
                    });
                }
            }
        }
        Ok(output)
    }
}

fn invalid(message: &str) -> Error {
    Error::InvalidInput(format!("parametric endpoint: {message}"))
}
fn limit(message: &str) -> Error {
    Error::Limit(format!("parametric endpoint: {message}"))
}
fn unsupported(message: &str) -> Error {
    Error::Unsupported(format!("parametric endpoint: {message}"))
}
fn checked_add(a: i64, b: i64) -> Result<i64> {
    a.checked_add(b)
        .ok_or_else(|| limit("integer bound overflow"))
}
fn checked_sub(a: i64, b: i64) -> Result<i64> {
    a.checked_sub(b)
        .ok_or_else(|| limit("integer bound overflow"))
}
fn checked_mul(a: i64, b: i64) -> Result<i64> {
    a.checked_mul(b)
        .ok_or_else(|| limit("integer bound overflow"))
}
fn symbol(name: &str) -> Atom {
    Atom::var(symbol!(format!("rustflow_parametric_endpoint::{name}")))
}
fn native_error(e: impl std::fmt::Display) -> Error {
    unsupported(&format!("HEPKit Symanzik: {e}"))
}

/// Query an actual occupied family at one active positive virtual support.
/// Pure external factors must be nonzero multiples of the spacelike transfer;
/// they are retained separately and are not integrated as virtual parameters.
/// No source or physical flow admission is changed by this function.
pub fn two_occupied_support_evidence(
    input: &PreparedDensityInput,
    family: &OccupiedCutFamily,
    active_slots: &[usize],
    budget: ParametricEndpointBudget,
) -> Result<ParametricSupportEvidence> {
    validate_budget(budget)?;
    if family.shells().len() != 2 || family.loops() <= 2 || family.loops() > 16 {
        return Err(unsupported(
            "requires two occupied and one or more virtual loops (total <=16)",
        ));
    }
    if input.physical_masses().iter().any(|m| !m.is_zero())
        || family
            .shells()
            .iter()
            .any(|s| !s.mass_squared.is_zero() || s.chemical_potential <= 0)
    {
        return Err(unsupported(
            "requires all zero masses and strictly positive occupied endpoints",
        ));
    }
    let cuts = family
        .shells()
        .iter()
        .map(|s| s.physical_slot)
        .collect::<Vec<_>>();
    let expected = input.occupied_cut(&cuts, 1024)?.at_physical_masses();
    let signature = super::massless_endpoint::bound_family_signature(family);
    if signature != super::massless_endpoint::bound_family_signature(&expected) {
        return Err(invalid("family is not the exact assigned input family"));
    }
    if active_slots.windows(2).any(|w| w[0] >= w[1]) {
        return Err(invalid("active slots must be sorted and distinct"));
    }
    let routing = input
        .input()
        .edges
        .iter()
        .map(|e| {
            e.routing
                .iter()
                .map(|s| {
                    Atom::parse(s, "rustflow_density", Default::default())
                        .map_err(|e| invalid(&e.to_string()))
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    let routed = matmul(&routing, family.inverse_routing());
    let mut rational_rows = Vec::new();
    for row in &routed {
        rational_rows.push(
            row.iter()
                .map(|a| {
                    Rational::try_from(a.as_view())
                        .map_err(|_| unsupported("routing must be exact rational"))
                })
                .collect::<Result<Vec<_>>>()?,
        );
    }
    let mut virtual_slots = Vec::new();
    let mut pure = Vec::new();
    for (slot, row) in rational_rows
        .iter()
        .enumerate()
        .filter(|(s, _)| !cuts.contains(s))
    {
        if row[2..].iter().all(|r| r.is_zero()) {
            if row[0].is_zero() || !(&row[0] + &row[1]).is_zero() {
                return Err(unsupported(
                    "external-only factor is not a nonzero spacelike transfer",
                ));
            }
            pure.push((slot, (&row[0] * &row[0]).to_string()));
        } else {
            virtual_slots.push(slot);
        }
    }
    if active_slots.iter().any(|s| !virtual_slots.contains(s)) {
        return Err(invalid(
            "active support contains a cut, completion, or pure external slot",
        ));
    }
    if active_slots.len() > budget.parameters {
        return Err(limit("parameter count exceeds budget"));
    }
    let rank = rref(
        active_slots
            .iter()
            .map(|&s| routed[s][2..].to_vec())
            .collect(),
    )
    .1
    .len();
    let h = family.loops() - 2;
    let mut evidence = ParametricSupportEvidence {
        version: VERSION,
        input_identity: input.identity().to_owned(),
        family_signature: signature,
        cuts,
        active_slots: active_slots.to_vec(),
        routed_rows: routed
            .iter()
            .map(|r| r.iter().map(Atom::to_canonical_string).collect())
            .collect(),
        virtual_loops: h,
        active_virtual_rank: rank,
        pure_transfer_slots: pure,
        u: Vec::new(),
        v: Vec::new(),
        charts: Vec::new(),
        resolution_operations: 0,
        full_kinematics: None,
    };
    if rank < h {
        return Ok(evidence);
    }
    let orders = permutations(active_slots.len(), budget)?;
    let full = native_full_kinematics(&routed, active_slots, 2, h, budget)?;
    evidence.u = full.u.clone();
    evidence.v = full.pairs[0].terms.clone();
    evidence.full_kinematics = Some(full);
    if evidence.u.is_empty() {
        return Err(unsupported("full-rank support has zero native U"));
    }
    let (us, fs) = supports(&evidence.u, &evidence.v, active_slots.len(), budget)?;
    let mut counters = Counters::default();
    for order in orders {
        let state = State::primary(&us, &fs, &order, budget)?;
        let tree = resolve(state, h, 0, budget, &mut counters)?;
        evidence.charts.push(EtaPrimaryChart { order, tree });
    }
    evidence.resolution_operations = counters.operations;
    Ok(evidence)
}

/// Enumerate every active positive virtual support, including the empty and
/// rank-deficient supports. Exhausted budgets return an error, never a partial
/// coverage certificate. This remains a query-only operation.
pub fn two_occupied_all_supports(
    input: &PreparedDensityInput,
    family: &OccupiedCutFamily,
    budget: ParametricEndpointBudget,
) -> Result<Vec<ParametricSupportEvidence>> {
    // Perform binding, shape and physical-domain checks before any matrix
    // multiplication or row slicing, including for the empty active support.
    let empty = two_occupied_support_evidence(input, family, &[], budget)?;
    let rows = empty
        .routed_rows
        .iter()
        .map(|row| {
            row.iter()
                .map(|s| {
                    Atom::parse(s, "rustflow_density", Default::default())
                        .map_err(|e| invalid(&e.to_string()))
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    let slots = rows
        .iter()
        .enumerate()
        .filter(|(s, r)| {
            !family.shells().iter().any(|c| c.physical_slot == *s)
                && r[2..].iter().any(|a| !a.is_zero())
        })
        .map(|(s, _)| s)
        .collect::<Vec<_>>();
    let count = 1usize
        .checked_shl(u32::try_from(slots.len()).map_err(|_| limit("support exponent overflow"))?)
        .ok_or_else(|| limit("support count overflow"))?;
    if count > budget.supports {
        return Err(limit("active support count exceeds budget"));
    }
    let mut remaining = budget;
    let mut output = Vec::with_capacity(count);
    for mask in 0..count {
        let active = slots
            .iter()
            .enumerate()
            .filter(|(i, _)| mask & (1usize << i) != 0)
            .map(|(_, s)| *s)
            .collect::<Vec<_>>();
        let evidence = two_occupied_support_evidence(input, family, &active, remaining)?;
        let mut nodes = 0usize;
        let mut leaves = 0usize;
        for chart in &evidence.charts {
            tree_size(&chart.tree, &mut nodes, &mut leaves)?;
        }
        remaining.primary_charts = remaining
            .primary_charts
            .checked_sub(evidence.charts.len())
            .ok_or_else(|| limit("aggregate chart budget exhausted"))?;
        remaining.tree_nodes = remaining
            .tree_nodes
            .checked_sub(nodes)
            .ok_or_else(|| limit("aggregate tree budget exhausted"))?;
        remaining.sectors = remaining
            .sectors
            .checked_sub(leaves)
            .ok_or_else(|| limit("aggregate sector budget exhausted"))?;
        remaining.operations = remaining
            .operations
            .checked_sub(evidence.resolution_operations)
            .ok_or_else(|| limit("aggregate work budget exhausted"))?;
        output.push(evidence);
    }
    Ok(output)
}

fn validate_budget(budget: ParametricEndpointBudget) -> Result<()> {
    if budget.parameters > 16 || budget.depth > 128 {
        return Err(limit("hard coordinate/depth limit exceeded"));
    }
    Ok(())
}
fn tree_size(node: &EtaSectorNode, nodes: &mut usize, leaves: &mut usize) -> Result<()> {
    *nodes = nodes
        .checked_add(1)
        .ok_or_else(|| limit("tree size overflow"))?;
    match node {
        EtaSectorNode::Leaf(_) => {
            *leaves = leaves
                .checked_add(1)
                .ok_or_else(|| limit("leaf size overflow"))?
        }
        EtaSectorNode::Blowup { children, .. } => {
            for (_, child) in children {
                tree_size(child, nodes, leaves)?;
            }
        }
    }
    Ok(())
}

fn positive_polynomial(
    atom: &Atom,
    variables: &[Atom],
    budget: ParametricEndpointBudget,
) -> Result<Vec<PositiveMonomial>> {
    let terms = exact_coefficient_list(atom, variables)?;
    if terms.len() > budget.polynomial_terms {
        return Err(limit("polynomial term count exceeds budget"));
    }
    terms
        .into_iter()
        .map(|(monomial, coefficient)| {
            let coefficient = Rational::try_from(coefficient.as_view())
                .map_err(|_| unsupported("polynomial coefficient is not rational"))?;
            if coefficient <= 0 {
                return Err(unsupported(
                    "Symanzik polynomial has a nonpositive coefficient",
                ));
            }
            let exponents = powers(&monomial, variables)?
                .into_iter()
                .map(|p| u32::try_from(p).map_err(|_| unsupported("negative parameter power")))
                .collect::<Result<Vec<_>>>()?;
            if exponents.iter().any(|&p| p > budget.exponent) {
                return Err(limit("polynomial exponent exceeds budget"));
            }
            Ok(PositiveMonomial {
                exponents,
                coefficient: coefficient.to_string(),
            })
        })
        .collect()
}
fn supports(
    u: &[PositiveMonomial],
    v: &[PositiveMonomial],
    n: usize,
    budget: ParametricEndpointBudget,
) -> Result<(Vec<Vec<u32>>, Vec<Vec<u32>>)> {
    if u.len()
        .checked_mul(n)
        .and_then(|s| s.checked_add(v.len()))
        .is_none_or(|s| s > budget.polynomial_terms)
    {
        return Err(limit("joint polynomial construction exceeds term budget"));
    }
    for p in u.iter().chain(v) {
        let a = Atom::parse(&p.coefficient, "rustflow_density", Default::default())
            .map_err(|e| invalid(&e.to_string()))?;
        if p.exponents.len() != n
            || Rational::try_from(a.as_view())
                .map_err(|_| invalid("nonrational retained coefficient"))?
                <= 0
        {
            return Err(invalid("invalid retained positive polynomial"));
        }
    }
    let us = u
        .iter()
        .map(|m| {
            let mut p = m.exponents.clone();
            p.push(0);
            p
        })
        .collect::<Vec<_>>();
    let mut fs = v
        .iter()
        .map(|m| {
            let mut p = m.exponents.clone();
            p.push(0);
            p
        })
        .collect::<BTreeSet<_>>();
    for u in &us {
        for i in 0..n {
            let mut p = u.clone();
            p[i] = p[i]
                .checked_add(1)
                .ok_or_else(|| limit("exponent overflow"))?;
            p[n] = 1;
            fs.insert(p);
        }
    }
    if us.is_empty() || fs.is_empty() {
        return Err(invalid("empty full-rank polynomial support"));
    }
    if fs.len() > budget.polynomial_terms {
        return Err(limit("joint polynomial support exceeds budget"));
    }
    Ok((us, fs.into_iter().collect()))
}

#[derive(Clone)]
struct State {
    u: Vec<Vec<u32>>,
    f: Vec<Vec<u32>>,
    jac: Vec<u32>,
    eta: Vec<u32>,
    alpha: Vec<Vec<u32>>,
}
impl State {
    fn primary(
        u: &[Vec<u32>],
        f: &[Vec<u32>],
        order: &[usize],
        budget: ParametricEndpointBudget,
    ) -> Result<Self> {
        let n = order.len();
        if n == 0 {
            return Err(invalid("empty full-rank primary chart"));
        }
        let convert = |term: &Vec<u32>| -> Result<Vec<u32>> {
            let mut out = Vec::new();
            for j in 1..n {
                out.push(
                    order[j..]
                        .iter()
                        .try_fold(0u32, |s, &i| s.checked_add(term[i]))
                        .ok_or_else(|| limit("Hepp exponent overflow"))?,
                );
            }
            out.push(term[n]);
            if out.iter().any(|&exponent| exponent > budget.exponent) {
                return Err(limit("primary Hepp exponent exceeds budget"));
            }
            Ok(out)
        };
        let alpha = (0..n)
            .map(|i| {
                let mut unit = vec![0; n + 1];
                unit[i] = 1;
                convert(&unit)
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            u: u.iter().map(convert).collect::<Result<_>>()?,
            f: f.iter().map(convert).collect::<Result<_>>()?,
            jac: (0..n).map(|i| u32::try_from(n - 1 - i).unwrap()).collect(),
            eta: (0..n).map(|i| u32::from(i == n - 1)).collect(),
            alpha,
        })
    }
    fn change(
        &self,
        axes: &[usize],
        pivot: usize,
        budget: ParametricEndpointBudget,
    ) -> Result<Self> {
        let change = |v: &Vec<u32>| -> Result<Vec<u32>> {
            let mut r = v.clone();
            r[pivot] = axes
                .iter()
                .try_fold(0u32, |s, &i| s.checked_add(v[i]))
                .ok_or_else(|| limit("sector exponent overflow"))?;
            if r[pivot] > budget.exponent {
                return Err(limit("sector exponent exceeds budget"));
            }
            Ok(r)
        };
        let poly = |p: &Vec<Vec<u32>>| -> Result<Vec<Vec<u32>>> {
            Ok(p.iter()
                .map(change)
                .collect::<Result<BTreeSet<_>>>()?
                .into_iter()
                .collect())
        };
        Ok(Self {
            u: poly(&self.u)?,
            f: poly(&self.f)?,
            jac: change(&self.jac)?,
            eta: change(&self.eta)?,
            alpha: self.alpha.iter().map(change).collect::<Result<_>>()?,
        })
    }
}
fn valuation(poly: &[Vec<u32>]) -> (Vec<u32>, bool) {
    let n = poly[0].len();
    let v = (0..n)
        .map(|i| poly.iter().map(|p| p[i]).min().unwrap())
        .collect::<Vec<_>>();
    let unit = poly.contains(&v);
    (v, unit)
}
fn splitting_axes(
    poly: &[Vec<u32>],
    v: &[u32],
    budget: ParametricEndpointBudget,
    c: &mut Counters,
) -> Result<Vec<usize>> {
    let n = v.len();
    if n >= usize::BITS as usize {
        return Err(limit("blowup axis count overflow"));
    }
    charge(
        c,
        (1usize << n)
            .checked_mul(poly.len())
            .and_then(|s| s.checked_mul(n))
            .ok_or_else(|| limit("center-search cost overflow"))?,
        budget,
    )?;
    let mut best: Option<(usize, u64, Vec<usize>)> = None;
    for mask in 1usize..(1usize << n) {
        if mask.count_ones() < 2 {
            continue;
        }
        let axes = (0..n).filter(|&i| mask & (1 << i) != 0).collect::<Vec<_>>();
        if !poly.iter().all(|p| axes.iter().any(|&i| p[i] > v[i])) {
            continue;
        }
        let cost = poly
            .iter()
            .flat_map(|p| axes.iter().map(move |&i| u64::from(p[i] - v[i])))
            .try_fold(0_u64, |sum, x| sum.checked_add(x))
            .ok_or_else(|| limit("center-search score overflow"))?;
        let candidate = (axes.len(), cost, axes);
        if best.as_ref().is_none_or(|b| candidate < *b) {
            best = Some(candidate);
        }
    }
    best.map(|b| b.2)
        .ok_or_else(|| unsupported("positive-unit resolution found no complete coordinate center"))
}
#[derive(Default)]
struct Counters {
    nodes: usize,
    leaves: usize,
    operations: usize,
}
fn charge(c: &mut Counters, work: usize, budget: ParametricEndpointBudget) -> Result<()> {
    c.operations = c
        .operations
        .checked_add(work)
        .ok_or_else(|| limit("work counter overflow"))?;
    if c.operations > budget.operations {
        return Err(limit("sector arithmetic work budget exhausted"));
    }
    Ok(())
}
fn visit(c: &mut Counters, depth: usize, budget: ParametricEndpointBudget) -> Result<()> {
    c.nodes = c
        .nodes
        .checked_add(1)
        .ok_or_else(|| limit("tree node overflow"))?;
    if c.nodes > budget.tree_nodes || depth > budget.depth {
        return Err(limit("sector tree node/depth budget exhausted"));
    }
    Ok(())
}
fn leaf(
    state: &State,
    h: usize,
    c: &mut Counters,
    budget: ParametricEndpointBudget,
) -> Result<EtaSectorLeaf> {
    c.leaves = c
        .leaves
        .checked_add(1)
        .ok_or_else(|| limit("leaf count overflow"))?;
    if c.leaves > budget.sectors {
        return Err(limit("sector leaf budget exhausted"));
    }
    let (u, uu) = valuation(&state.u);
    let (f, ff) = valuation(&state.f);
    if !uu || !ff {
        return Err(invalid("leaf is not a pair of positive units"));
    }
    let mut regular = 0;
    for i in 0..state.eta.len() {
        let m = state.eta[i];
        if m == 0 {
            continue;
        }
        let d = checked_sub(
            checked_mul(h as i64, i64::from(f[i]))?,
            checked_mul(h as i64 + 1, i64::from(u[i]))?,
        )?;
        if d < 0 {
            return Err(unsupported("eta Mellin axis has a negative D slope"));
        }
        if d == 0 {
            if u[i] != 0
                || f[i] != 0
                || state.jac[i] != 0
                || m != 1
                || state.alpha.iter().any(|a| a[i] != 0)
            {
                return Err(unsupported(
                    "zero-D eta axis is not a regular parameter-independent axis",
                ));
            }
            regular += 1;
        }
        if state.eta[i].checked_add(u[i]).is_none_or(|x| x < f[i]) {
            return Err(invalid("eta Euler insertion has a negative valuation"));
        }
    }
    if regular > 1 {
        return Err(unsupported(
            "multiple regular eta axes could produce endpoint logarithms",
        ));
    }
    Ok(EtaSectorLeaf {
        u,
        f,
        jacobian_plus_one: state.jac.clone(),
        eta: state.eta.clone(),
        original_parameters: state.alpha.clone(),
    })
}
fn resolve(
    state: State,
    h: usize,
    depth: usize,
    budget: ParametricEndpointBudget,
    c: &mut Counters,
) -> Result<EtaSectorNode> {
    visit(c, depth, budget)?;
    charge(
        c,
        state
            .u
            .len()
            .checked_add(state.f.len())
            .and_then(|s| s.checked_add(state.alpha.len() + 2))
            .and_then(|s| s.checked_mul(state.eta.len()))
            .ok_or_else(|| limit("state work overflow"))?,
        budget,
    )?;
    let (f, ff) = valuation(&state.f);
    let (u, uu) = valuation(&state.u);
    if ff && uu {
        return Ok(EtaSectorNode::Leaf(leaf(&state, h, c, budget)?));
    }
    let axes = if !ff {
        splitting_axes(&state.f, &f, budget, c)?
    } else {
        splitting_axes(&state.u, &u, budget, c)?
    };
    let children = axes
        .iter()
        .map(|&p| {
            Ok((
                p,
                Box::new(resolve(
                    state.change(&axes, p, budget)?,
                    h,
                    depth + 1,
                    budget,
                    c,
                )?),
            ))
        })
        .collect::<Result<_>>()?;
    Ok(EtaSectorNode::Blowup { axes, children })
}
fn replay_node(
    node: &EtaSectorNode,
    state: &State,
    h: usize,
    depth: usize,
    budget: ParametricEndpointBudget,
    c: &mut Counters,
) -> Result<()> {
    visit(c, depth, budget)?;
    charge(
        c,
        state
            .u
            .len()
            .checked_add(state.f.len())
            .and_then(|s| s.checked_add(state.alpha.len() + 2))
            .and_then(|s| s.checked_mul(state.eta.len()))
            .ok_or_else(|| limit("replay work overflow"))?,
        budget,
    )?;
    match node {
        EtaSectorNode::Leaf(saved) => {
            if &leaf(state, h, c, budget)? != saved {
                return Err(invalid("retained leaf does not replay"));
            }
        }
        EtaSectorNode::Blowup { axes, children } => {
            if axes.len() < 2
                || axes.windows(2).any(|w| w[0] >= w[1])
                || axes.iter().any(|&i| i >= state.eta.len())
                || children.len() != axes.len()
            {
                return Err(invalid("incomplete or invalid blowup partition"));
            }
            for (&axis, (pivot, child)) in axes.iter().zip(children) {
                if axis != *pivot {
                    return Err(invalid("missing or repeated blowup pivot"));
                }
                replay_node(
                    child,
                    &state.change(axes, *pivot, budget)?,
                    h,
                    depth + 1,
                    budget,
                    c,
                )?;
            }
        }
    }
    Ok(())
}
fn collect_leaves<'a>(node: &'a EtaSectorNode, out: &mut Vec<&'a EtaSectorLeaf>) {
    match node {
        EtaSectorNode::Leaf(l) => out.push(l),
        EtaSectorNode::Blowup { children, .. } => {
            for (_, child) in children {
                collect_leaves(child, out);
            }
        }
    }
}
fn permutations(n: usize, budget: ParametricEndpointBudget) -> Result<Vec<Vec<usize>>> {
    if n == 0 || n > budget.parameters {
        return Err(limit("primary parameter count exceeds budget"));
    }
    let count = (1..=n)
        .try_fold(1usize, |s, i| s.checked_mul(i))
        .ok_or_else(|| limit("primary chart count overflow"))?;
    if count > budget.primary_charts {
        return Err(limit("primary chart count exceeds budget"));
    }
    fn next(current: &mut Vec<usize>, used: &mut [bool], out: &mut Vec<Vec<usize>>) {
        if current.len() == used.len() {
            out.push(current.clone());
            return;
        }
        for i in 0..used.len() {
            if !used[i] {
                used[i] = true;
                current.push(i);
                next(current, used, out);
                current.pop();
                used[i] = false;
            }
        }
    }
    let mut out = Vec::with_capacity(count);
    next(&mut Vec::new(), &mut vec![false; n], &mut out);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_density::DensityInput;

    fn prism() -> PreparedDensityInput {
        serde_json::from_str::<DensityInput>(include_str!(
            "../../examples/finite_density/triangular_prism.json"
        ))
        .unwrap()
        .prepare()
        .unwrap()
    }
    fn full_support(family: &OccupiedCutFamily, input: &PreparedDensityInput) -> Vec<usize> {
        let empty = two_occupied_support_evidence(input, family, &[], Default::default()).unwrap();
        (0..family.physical_slots())
            .filter(|s| {
                !empty.cuts.contains(s) && !empty.pure_transfer_slots.iter().any(|(p, _)| p == s)
            })
            .collect()
    }
    fn volume(node: &EtaSectorNode) -> Rational {
        match node {
            EtaSectorNode::Leaf(l) => l
                .jacobian_plus_one
                .iter()
                .zip(&l.eta)
                .fold(Rational::one(), |v, (&j, &m)| v / Rational::from(j + m)),
            EtaSectorNode::Blowup { children, .. } => children
                .iter()
                .fold(Rational::zero(), |s, (_, c)| s + volume(c)),
        }
    }
    #[test]
    fn actual_prism_all_virtual_supports_have_replayable_positive_sector_evidence() {
        let input = prism();
        for cuts in [[1, 5], [1, 7], [5, 7]] {
            let family = input
                .occupied_cut(&cuts, 1024)
                .unwrap()
                .at_physical_masses();
            let evidence = two_occupied_all_supports(&input, &family, Default::default()).unwrap();
            assert_eq!(evidence.len(), 64);
            assert!(evidence[0].is_unrestricted_virtual_polynomial());
            let mut deficient = 0;
            let mut massive_germ = 0;
            let mut regular = 0;
            for support in &evidence {
                support.replay(Default::default()).unwrap();
                if support.is_unrestricted_virtual_polynomial() {
                    deficient += 1;
                    continue;
                }
                massive_germ += usize::from(support.v().is_empty());
                let bounds = support
                    .bounds_for_polynomial_rank(&vec![1; support.active_slots.len()], 4, 1)
                    .unwrap();
                assert!(
                    bounds
                        .iter()
                        .all(|b| b.twice_dimension_coefficient > 0 || b.regular && b.constant == 0)
                );
                regular += usize::from(bounds.iter().any(|b| b.regular));
                // Independent Jacobian check: the complete maximum-gauge
                // projective charts have one unit cube per choice of maximum.
                let sum = support
                    .charts
                    .iter()
                    .fold(Rational::zero(), |s, c| s + volume(&c.tree));
                assert_eq!(sum, Rational::from(support.active_slots.len() as i64));
            }
            assert!(deficient > 0 && massive_germ > 0 && regular > 0);
            if let Some(directory) = std::env::var_os("RUSTFLOW_PARAMETRIC_EVIDENCE_REPORT") {
                let directory = std::path::PathBuf::from(directory);
                std::fs::create_dir_all(&directory).unwrap();
                let path =
                    directory.join(format!("prism-cuts-{}-{}-supports.json", cuts[0], cuts[1]));
                let file = std::io::BufWriter::new(std::fs::File::create(&path).unwrap());
                serde_json::to_writer(file, &serde_json::json!({
                    "scope": "query-only HEPKit exact Symanzik and complete sector evidence; no flow admission or period value",
                    "source": "src/finite_density/parametric_endpoint.rs",
                    "source_blake3": blake3::hash(include_bytes!("parametric_endpoint.rs")).to_hex().to_string(),
                    "input_fixture": "examples/finite_density/triangular_prism.json",
                    "input_fixture_blake3": blake3::hash(include_bytes!("../../examples/finite_density/triangular_prism.json")).to_hex().to_string(),
                    "cuts": cuts,
                    "supports": evidence,
                })).unwrap();
                eprintln!("Saved exact support evidence: {}", path.display());
            }
            eprintln!(
                "prism cuts {cuts:?}: 64 supports, {deficient} rank-deficient, {massive_germ} full-rank V=0, {regular} with regular eta branch"
            );
            let full = evidence.last().unwrap();
            assert_eq!(full.active_slots.len(), 6);
            assert_eq!(full.charts.len(), 720);
            let mut powers = vec![1; 6];
            powers[0] = 2;
            let bounds = full.bounds_for_polynomial_rank(&powers, 6, 0).unwrap();
            assert!(bounds.iter().all(|b| b.regular && b.constant == 0
                || b.twice_dimension_coefficient * 19 + 2 * b.constant > 0));
            // Query evidence deliberately does not widen the sealed permit.
            let shifted = (0..family.physical_slots())
                .filter(|s| !cuts.contains(s))
                .collect::<Vec<_>>();
            assert!(
                super::super::massless_endpoint::MasslessFlowEvidence::new(
                    &input,
                    &family,
                    &shifted,
                    Default::default()
                )
                .is_err()
            );
        }
    }

    #[test]
    fn native_sector_query_rejects_incomplete_trees_and_wrong_binding() {
        let input = prism();
        let family = input
            .occupied_cut(&[1, 5], 1024)
            .unwrap()
            .at_physical_masses();
        let active = full_support(&family, &input);
        let evidence =
            two_occupied_support_evidence(&input, &family, &active, Default::default()).unwrap();
        evidence
            .replay_for(&input, &family, Default::default())
            .unwrap();
        let original = evidence.identity().unwrap();
        let mut missing = evidence.clone();
        missing.charts.pop();
        assert_ne!(original, missing.identity().unwrap());
        assert!(missing.replay(Default::default()).is_err());
        let mut broken = evidence.clone();
        let node = broken
            .charts
            .iter_mut()
            .find_map(|c| {
                if let EtaSectorNode::Blowup { children, .. } = &mut c.tree {
                    Some(children)
                } else {
                    None
                }
            })
            .unwrap();
        node.pop();
        assert!(broken.replay(Default::default()).is_err());
        let mut wrong = evidence.clone();
        wrong.u[0].coefficient = "-1".into();
        assert!(wrong.replay(Default::default()).is_err());
        // Use an actual input of a different loop rank to exercise validation
        // before matmul rather than relying only on matching family shapes.
        let mut other: DensityInput = serde_json::from_str(include_str!(
            "../../examples/finite_density/massive_two_loop_sunset.json"
        ))
        .unwrap();
        for edge in &mut other.edges {
            edge.mass_squared = "0".into();
        }
        let other = other.prepare().unwrap();
        assert!(two_occupied_all_supports(&other, &family, Default::default()).is_err());
        let x = symbol("negative_test");
        assert!(positive_polynomial(&(-&x), &[x.clone()], Default::default()).is_err());
        assert!(positive_polynomial(&x.clone().pow(-1), &[x], Default::default()).is_err());
    }

    #[test]
    fn sector_query_budgets_fail_without_partial_support_coverage() {
        // Every original exponent is admissible; only their primary-chart
        // sum exceeds this limit. Catch it before allocating sector trees.
        assert!(matches!(State::primary(
            &[vec![0, 2, 2, 0]], &[vec![0, 2, 2, 0]], &[0, 1, 2],
            ParametricEndpointBudget { exponent: 3, ..Default::default() }
        ), Err(Error::Limit(ref message)) if message.contains("Hepp exponent")));
        let input = prism();
        let family = input
            .occupied_cut(&[1, 5], 1024)
            .unwrap()
            .at_physical_masses();
        let active = full_support(&family, &input);
        let b = ParametricEndpointBudget::default();
        for budget in [
            ParametricEndpointBudget { supports: 63, ..b },
            ParametricEndpointBudget { operations: 0, ..b },
            ParametricEndpointBudget { tree_nodes: 0, ..b },
            ParametricEndpointBudget { sectors: 0, ..b },
            ParametricEndpointBudget { depth: 0, ..b },
            ParametricEndpointBudget {
                primary_charts: 719,
                ..b
            },
            ParametricEndpointBudget {
                polynomial_terms: 1,
                ..b
            },
            ParametricEndpointBudget { exponent: 0, ..b },
        ] {
            let result = if budget.supports < 64 {
                two_occupied_all_supports(&input, &family, budget).map(|_| ())
            } else {
                two_occupied_support_evidence(&input, &family, &active, budget).map(|_| ())
            };
            assert!(
                matches!(result, Err(Error::Limit(_))),
                "unexpected budget result: {result:?}"
            );
        }
        assert!(two_occupied_support_evidence(&input, &family, &[1], b).is_err());
        assert!(two_occupied_support_evidence(&input, &family, &[0, 0], b).is_err());
    }
}
/// Signed real parameter coefficients are retained for the diagonal external
/// masses. Their sign is immaterial for same-sign pure-imaginary Feynman masses.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct RealMonomial {
    pub exponents: Vec<u32>,
    pub coefficient: String,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct PairSymanzikPolynomial {
    /// This multiplies 2 q_i.q_j, before null specialization.
    pub pair: [usize; 2],
    pub terms: Vec<PositiveMonomial>,
}
/// Exact native identity before any occupied mass is set to zero:
/// F=eta U sum(alpha)+sum_ij(2qi.qj)Vij+sum_i qi² Wi
///   +sum_j virtual_mass_j² U alpha_j.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct FullKinematicSymanzik {
    pub compact_loops: usize,
    pub virtual_loops: usize,
    pub parameter_slots: Vec<usize>,
    pub u: Vec<PositiveMonomial>,
    pub pairs: Vec<PairSymanzikPolynomial>,
    pub diagonal_masses: Vec<Vec<RealMonomial>>,
    pub virtual_masses: Vec<Vec<RealMonomial>>,
    pub regulator_prescription: &'static str,
}

fn real_polynomial(
    atom: &Atom,
    variables: &[Atom],
    budget: ParametricEndpointBudget,
) -> Result<Vec<RealMonomial>> {
    let terms = exact_coefficient_list(atom, variables)?;
    if terms.len() > budget.polynomial_terms {
        return Err(limit("real polynomial term count exceeds budget"));
    }
    terms
        .into_iter()
        .map(|(monomial, coefficient)| {
            let coefficient = Rational::try_from(coefficient.as_view())
                .map_err(|_| unsupported("off-null coefficient is not real rational"))?;
            let exponents = powers(&monomial, variables)?
                .into_iter()
                .map(|p| {
                    u32::try_from(p).map_err(|_| unsupported("negative off-null parameter power"))
                })
                .collect::<Result<Vec<_>>>()?;
            if exponents.iter().any(|&p| p > budget.exponent) {
                return Err(limit("off-null polynomial exponent exceeds budget"));
            }
            Ok(RealMonomial {
                exponents,
                coefficient: coefficient.to_string(),
            })
        })
        .collect()
}
fn derivative_in(atom: &Atom, variable: &Atom) -> Atom {
    let AtomView::Var(v) = variable.as_view() else {
        unreachable!()
    };
    atom.derivative(v.get_symbol()).expand()
}

fn native_full_kinematics(
    routed: &[Vec<Atom>],
    active_slots: &[usize],
    compact_loops: usize,
    virtual_loops: usize,
    budget: ParametricEndpointBudget,
) -> Result<FullKinematicSymanzik> {
    if compact_loops < 2
        || virtual_loops == 0
        || compact_loops + virtual_loops > 16
        || active_slots.is_empty()
        || active_slots.len() > budget.parameters
    {
        return Err(unsupported(
            "full kinematic query requires bounded nonempty virtual support",
        ));
    }
    let loops = (0..virtual_loops)
        .map(|i| symbol(&format!("k{i}")))
        .collect::<Vec<_>>();
    let external = (0..compact_loops)
        .map(|i| symbol(&format!("q{i}")))
        .collect::<Vec<_>>();
    let diagonals = (0..compact_loops)
        .map(|i| symbol(&format!("shell_mass_squared_{i}")))
        .collect::<Vec<_>>();
    let virtual_masses = (0..active_slots.len())
        .map(|i| symbol(&format!("virtual_mass_squared_{i}")))
        .collect::<Vec<_>>();
    let mut kin = Kinematics::new()
        .with_momenta(loops.iter().chain(&external).cloned())
        .map_err(native_error)?;
    for (q, m) in external.iter().zip(&diagonals) {
        kin = kin.with_mass_squared(q, m.clone()).map_err(native_error)?;
    }
    let mut pairs = Vec::new();
    for i in 0..compact_loops {
        for j in i + 1..compact_loops {
            let t = symbol(&format!("twice_dot_{i}_{j}"));
            kin = kin
                .with_scalar_product(&external[i], &external[j], &t / 2)
                .map_err(native_error)?;
            pairs.push(([i, j], t));
        }
    }
    let eta = symbol("eta");
    let momenta = external.iter().chain(&loops).cloned().collect::<Vec<_>>();
    let denominators = active_slots
        .iter()
        .zip(&virtual_masses)
        .map(|(&s, m)| {
            let row = routed
                .get(s)
                .ok_or_else(|| invalid("missing native routed row"))?;
            if row.len() != momenta.len() {
                return Err(invalid("native routed row dimension mismatch"));
            }
            let q = row.iter().zip(&momenta).map(|(a, k)| a * k).sum::<Atom>();
            Ok(kin.scalar_product(&q, &q).map_err(native_error)? - &eta - m)
        })
        .collect::<Result<Vec<_>>>()?;
    let native = NativeFamily::new(loops, external, denominators, &kin).map_err(native_error)?;
    let parameters = (0..active_slots.len())
        .map(|i| symbol(&format!("alpha{i}")))
        .collect::<Vec<_>>();
    let (u, f) = native.symanzik(&parameters).map_err(native_error)?;
    let mut rebuilt = &eta * &u * parameters.iter().cloned().sum::<Atom>();
    let mut pair_polynomials = Vec::new();
    let mut diagonal_polynomials = Vec::new();
    let mut mass_polynomials = Vec::new();
    let positive_u = positive_polynomial(&u, &parameters, budget)?;
    let mut count = positive_u.len();
    for (pair, t) in pairs {
        let coefficient = derivative_in(&f, &t);
        rebuilt += &t * &coefficient;
        let terms = positive_polynomial(&coefficient, &parameters, budget)?;
        count = count
            .checked_add(terms.len())
            .ok_or_else(|| limit("full kinematic term count overflow"))?;
        pair_polynomials.push(PairSymanzikPolynomial { pair, terms });
    }
    for a in &diagonals {
        let coefficient = derivative_in(&f, a);
        rebuilt += a * &coefficient;
        let terms = real_polynomial(&coefficient, &parameters, budget)?;
        count = count
            .checked_add(terms.len())
            .ok_or_else(|| limit("full kinematic term count overflow"))?;
        diagonal_polynomials.push(terms);
    }
    for (a, parameter) in virtual_masses.iter().zip(&parameters) {
        let coefficient = derivative_in(&f, a);
        if !(&coefficient - &u * parameter).expand().is_zero() {
            return Err(unsupported(
                "native virtual mass coefficient is not U alpha",
            ));
        }
        rebuilt += a * &coefficient;
        let terms = real_polynomial(&coefficient, &parameters, budget)?;
        count = count
            .checked_add(terms.len())
            .ok_or_else(|| limit("full kinematic term count overflow"))?;
        mass_polynomials.push(terms);
    }
    if count > budget.polynomial_terms {
        return Err(limit("full kinematic polynomial total exceeds budget"));
    }
    if !(f - rebuilt).expand().is_zero() {
        return Err(unsupported(
            "native F has unretained nonlinear kinematic dependence",
        ));
    }
    Ok(FullKinematicSymanzik {
        compact_loops,
        virtual_loops,
        parameter_slots: active_slots.to_vec(),
        u: positive_u,
        pairs: pair_polynomials,
        diagonal_masses: diagonal_polynomials,
        virtual_masses: mass_polynomials,
        regulator_prescription: "independent-native-mass-squared-minus-i-epsilon;real-spatial-future-sqrt-poles;epsilon-removed-at-fixed-positive-T-and-eta",
    })
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct PurePairTransfer {
    pub slot: usize,
    pub pair: [usize; 2],
    /// The actual deformed factor is eta + coefficient * h_ij.
    pub coefficient: String,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct CompactEndpointDegree {
    pub stratum: String,
    pub twice_dimension_coefficient: i64,
    pub constant: i64,
}
/// Query-only bounded one-Gaussian-block and multipair compact evidence.
/// This is not a `MasslessFlowEvidence` token and cannot enable source zeros.
#[derive(Clone, Debug, Serialize)]
pub struct OneVirtualCompactSupportEvidence {
    version: &'static str,
    input_identity: String,
    family_signature: Vec<String>,
    compact_loops: usize,
    cuts: Vec<usize>,
    active_slots: Vec<usize>,
    routed_rows: Vec<Vec<String>>,
    pure_transfers: Vec<PurePairTransfer>,
    full_kinematics: Option<FullKinematicSymanzik>,
    a_lower: Option<String>,
    a_upper: Option<String>,
    /// Bound for V_ij/U in the Gaussian completed-square convention, sum alpha=1.
    pair_q_upper: Vec<String>,
}
impl OneVirtualCompactSupportEvidence {
    pub fn label_envelope(
        &self,
        family: &OccupiedCutFamily,
        label: &crate::Integral,
    ) -> Result<ParametricLabelEnvelope> {
        if super::massless_endpoint::bound_family_signature(family) != self.family_signature {
            return Err(invalid("multipair label family binding mismatch"));
        }
        parametric_label_envelope(
            family,
            label,
            &self
                .pure_transfers
                .iter()
                .map(|t| t.slot)
                .collect::<Vec<_>>(),
        )
    }
    pub fn bounds_for_label(
        &self,
        family: &OccupiedCutFamily,
        label: &crate::Integral,
    ) -> Result<OneVirtualLabelBoundQuery> {
        let envelope = self.label_envelope(family, label)?;
        validate_envelope_support(&envelope, &self.active_slots)?;
        let powers = self
            .pure_transfers
            .iter()
            .map(|t| {
                envelope
                    .positive_transfer_powers
                    .iter()
                    .find(|(s, _)| *s == t.slot)
                    .map(|(_, n)| *n)
                    .unwrap_or(0)
            })
            .collect::<Vec<_>>();
        let compact_bounds = self.bounds_for_indices(
            &envelope.positive_virtual_powers,
            &powers,
            &envelope.shell_and_upper_jets,
        )?;
        Ok(OneVirtualLabelBoundQuery {
            envelope,
            compact_bounds,
        })
    }
    pub fn is_unrestricted_virtual_polynomial(&self) -> bool {
        self.active_slots.is_empty()
    }
    pub fn full_kinematics(&self) -> Option<&FullKinematicSymanzik> {
        self.full_kinematics.as_ref()
    }
    pub fn pure_transfers(&self) -> &[PurePairTransfer] {
        &self.pure_transfers
    }
    pub fn active_slots(&self) -> &[usize] {
        &self.active_slots
    }
    pub fn identity(&self) -> Result<String> {
        Ok(format!(
            "{}:{}",
            self.version,
            blake3::hash(&serde_json::to_vec(self).map_err(|e| Error::Cache(e.to_string()))?)
                .to_hex()
        ))
    }
    pub fn replay_for(
        &self,
        input: &PreparedDensityInput,
        family: &OccupiedCutFamily,
        budget: ParametricEndpointBudget,
    ) -> Result<()> {
        let fresh =
            one_virtual_compact_support_evidence(input, family, &self.active_slots, budget)?;
        if fresh.identity()? != self.identity()? {
            return Err(invalid("one-virtual native evidence does not replay"));
        }
        Ok(())
    }
    /// Conservative sufficient bounds for caller-supplied finite index data.
    /// No physical label or regulator admission is inferred by this query.
    pub fn bounds_for_indices(
        &self,
        positive_virtual_powers: &[u32],
        pure_transfer_powers: &[u32],
        jets: &[u32],
    ) -> Result<Vec<CompactEndpointDegree>> {
        if positive_virtual_powers.len() != self.active_slots.len()
            || positive_virtual_powers.contains(&0)
            || pure_transfer_powers.len() != self.pure_transfers.len()
            || jets.len() != self.compact_loops
        {
            return Err(invalid(
                "one-virtual bound data does not match its exact support",
            ));
        }
        let sum = |values: &[u32]| {
            values
                .iter()
                .try_fold(0_i64, |s, &n| checked_add(s, i64::from(n)))
        };
        let p = sum(positive_virtual_powers)?;
        let p0 = sum(pure_transfer_powers)?;
        let j = sum(jets)?;
        let transfer_and_jets = checked_add(p0, j)?;
        let mut out = vec![CompactEndpointDegree {
            stratum: "positive real-D compact Gram measure".into(),
            twice_dimension_coefficient: 2,
            constant: -(self.compact_loops as i64 + 1),
        }];
        if !self.is_unrestricted_virtual_polynomial() {
            out.push(CompactEndpointDegree {
                stratum: "bounded one-virtual Gaussian jets".into(),
                twice_dimension_coefficient: 1,
                constant: -checked_add(p, j)?,
            });
        }
        if transfer_and_jets > 0 {
            out.push(CompactEndpointDegree {
                stratum: "Holder bound for correlated compact pair angles".into(),
                twice_dimension_coefficient: 1,
                constant: -checked_add(1, transfer_and_jets)?,
            });
        }
        for (i, &ji) in jets.iter().enumerate() {
            let incident = self
                .pure_transfers
                .iter()
                .zip(pure_transfer_powers)
                .filter(|(t, _)| t.pair.contains(&i))
                .try_fold(0_i64, |s, (_, n)| checked_add(s, i64::from(*n)))?;
            out.push(CompactEndpointDegree {
                stratum: format!("compact radial origin {i}"),
                twice_dimension_coefficient: 2,
                constant: -checked_add(
                    2,
                    checked_add(checked_mul(2, i64::from(ji))?, checked_add(incident, j)?)?,
                )?,
            });
        }
        Ok(out)
    }
}

pub fn one_virtual_compact_support_evidence(
    input: &PreparedDensityInput,
    family: &OccupiedCutFamily,
    active_slots: &[usize],
    budget: ParametricEndpointBudget,
) -> Result<OneVirtualCompactSupportEvidence> {
    validate_budget(budget)?;
    let k = family.shells().len();
    if k < 2 || family.loops() != k + 1 || family.loops() > 16 {
        return Err(unsupported(
            "multipair Gaussian query requires exactly one virtual loop and at least two compacts",
        ));
    }
    if input.physical_masses().iter().any(|m| !m.is_zero())
        || family
            .shells()
            .iter()
            .any(|s| !s.mass_squared.is_zero() || s.chemical_potential <= 0)
    {
        return Err(unsupported(
            "multipair Gaussian query requires zero masses and positive chemical endpoints",
        ));
    }
    let cuts = family
        .shells()
        .iter()
        .map(|s| s.physical_slot)
        .collect::<Vec<_>>();
    let signature = super::massless_endpoint::bound_family_signature(family);
    if signature
        != super::massless_endpoint::bound_family_signature(
            &input.occupied_cut(&cuts, 1024)?.at_physical_masses(),
        )
    {
        return Err(invalid("multipair Gaussian family binding mismatch"));
    }
    if active_slots.windows(2).any(|w| w[0] >= w[1]) {
        return Err(invalid(
            "multipair active slots must be sorted and distinct",
        ));
    }
    if active_slots.len() > budget.parameters {
        return Err(limit("multipair parameter budget exhausted"));
    }
    let routing = input
        .input()
        .edges
        .iter()
        .map(|e| {
            e.routing
                .iter()
                .map(|s| {
                    Atom::parse(s, "rustflow_density", Default::default())
                        .map_err(|e| invalid(&e.to_string()))
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    let routed = matmul(&routing, family.inverse_routing());
    let rational_rows = routed
        .iter()
        .map(|row| {
            row.iter()
                .map(|a| {
                    Rational::try_from(a.as_view())
                        .map_err(|_| unsupported("multipair routing must be rational"))
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    let mut virtual_slots = Vec::new();
    let mut pure_transfers = Vec::new();
    for (s, row) in rational_rows
        .iter()
        .enumerate()
        .filter(|(s, _)| !cuts.contains(s))
    {
        if !row[k].is_zero() {
            virtual_slots.push(s);
            continue;
        }
        let nonzero = (0..k).filter(|&i| !row[i].is_zero()).collect::<Vec<_>>();
        if nonzero.len() != 2 || !(&row[nonzero[0]] + &row[nonzero[1]]).is_zero() {
            return Err(unsupported(
                "multipair external factor is not exactly one future difference",
            ));
        }
        pure_transfers.push(PurePairTransfer {
            slot: s,
            pair: [nonzero[0], nonzero[1]],
            coefficient: (&row[nonzero[0]] * &row[nonzero[0]]).to_string(),
        });
    }
    if active_slots.iter().any(|s| !virtual_slots.contains(s)) {
        return Err(invalid("multipair support contains a nonvirtual slot"));
    }
    let mut evidence = OneVirtualCompactSupportEvidence {
        version: "one-virtual-positive-multipair-Gaussian-Holder-v1;query-only;full-off-null-Feynman-mass-identity",
        input_identity: input.identity().into(),
        family_signature: signature,
        compact_loops: k,
        cuts,
        active_slots: active_slots.to_vec(),
        routed_rows: routed
            .iter()
            .map(|r| r.iter().map(Atom::to_canonical_string).collect())
            .collect(),
        pure_transfers,
        full_kinematics: None,
        a_lower: None,
        a_upper: None,
        pair_q_upper: Vec::new(),
    };
    if active_slots.is_empty() {
        return Ok(evidence);
    }
    let full = native_full_kinematics(&routed, active_slots, k, 1, budget)?;
    let scales = active_slots
        .iter()
        .map(|&s| &rational_rows[s][k] * &rational_rows[s][k])
        .collect::<Vec<_>>();
    let lower = scales.iter().min().unwrap().clone();
    let upper = scales.iter().max().unwrap().clone();
    let parameters = (0..active_slots.len())
        .map(|i| symbol(&format!("alpha{i}")))
        .collect::<Vec<_>>();
    let exact_a = scales
        .iter()
        .zip(&parameters)
        .map(|(a, x)| Atom::num(a.clone()) * x)
        .sum::<Atom>();
    if positive_polynomial(&exact_a, &parameters, budget)? != full.u {
        return Err(unsupported(
            "native one-loop U differs from exact sum alpha*a²",
        ));
    }
    for pair in &full.pairs {
        let sum = pair.terms.iter().try_fold(Rational::zero(), |s, t| {
            let a = Atom::parse(&t.coefficient, "rustflow_density", Default::default())
                .map_err(|e| invalid(&e.to_string()))?;
            Ok::<_, Error>(
                s + Rational::try_from(a.as_view())
                    .map_err(|_| invalid("nonrational retained pair coefficient"))?,
            )
        })?;
        evidence.pair_q_upper.push((sum / &lower).to_string());
    }
    evidence.a_lower = Some(lower.to_string());
    evidence.a_upper = Some(upper.to_string());
    evidence.full_kinematics = Some(full);
    Ok(evidence)
}

/// Complete active-support coverage for the one-virtual multipair query.
/// Empty support denotes an unrestricted virtual polynomial, not a massive
/// endpoint monomial. This operation never enables source-domain zeros.
/// The polynomial-term budget bounds the sum of all retained support polynomials.
pub fn one_virtual_compact_all_supports(
    input: &PreparedDensityInput,
    family: &OccupiedCutFamily,
    budget: ParametricEndpointBudget,
) -> Result<Vec<OneVirtualCompactSupportEvidence>> {
    let empty = one_virtual_compact_support_evidence(input, family, &[], budget)?;
    let slots = (0..empty.routed_rows.len())
        .filter(|slot| {
            !empty.cuts.contains(slot) && !empty.pure_transfers.iter().any(|t| t.slot == *slot)
        })
        .collect::<Vec<_>>();
    let count = 1usize
        .checked_shl(
            u32::try_from(slots.len()).map_err(|_| limit("multipair support count overflow"))?,
        )
        .ok_or_else(|| limit("multipair support count overflow"))?;
    if count > budget.supports {
        return Err(limit("multipair support budget exhausted"));
    }
    let mut out = Vec::with_capacity(count);
    let mut terms = 0_usize;
    for mask in 0..count {
        let active = slots
            .iter()
            .enumerate()
            .filter(|(i, _)| mask & (1usize << i) != 0)
            .map(|(_, s)| *s)
            .collect::<Vec<_>>();
        let evidence = one_virtual_compact_support_evidence(input, family, &active, budget)?;
        if let Some(full) = &evidence.full_kinematics {
            terms = terms
                .checked_add(full_kinematic_term_count(full)?)
                .ok_or_else(|| limit("multipair aggregate term count overflow"))?;
            if terms > budget.polynomial_terms {
                return Err(limit("multipair aggregate polynomial budget exhausted"));
            }
        }
        out.push(evidence);
    }
    Ok(out)
}

fn full_kinematic_term_count(full: &FullKinematicSymanzik) -> Result<usize> {
    full.pairs
        .iter()
        .map(|p| p.terms.len())
        .chain(full.diagonal_masses.iter().map(Vec::len))
        .chain(full.virtual_masses.iter().map(Vec::len))
        .try_fold(full.u.len(), |n, count| {
            n.checked_add(count)
                .ok_or_else(|| limit("full kinematic term count overflow"))
        })
}

#[cfg(test)]
mod full_kinematic_tests {
    use super::*;
    use crate::finite_density::DensityInput;
    #[test]
    fn full_native_symanzik_retains_off_null_and_independent_mass_terms() {
        let rows = vec![
            vec![Atom::one(), Atom::zero(), Atom::one()],
            vec![Atom::zero(), Atom::one(), Atom::one()],
        ];
        let full = native_full_kinematics(&rows, &[0, 1], 2, 1, Default::default()).unwrap();
        assert_eq!(full.pairs.len(), 1);
        assert_eq!(full.pairs[0].pair, [0, 1]);
        assert_eq!(
            full.pairs[0].terms,
            vec![PositiveMonomial {
                exponents: vec![1, 1],
                coefficient: "1".into()
            }]
        );
        let signed = vec![RealMonomial {
            exponents: vec![1, 1],
            coefficient: "-1".into(),
        }];
        assert_eq!(full.diagonal_masses, vec![signed.clone(), signed]);
        assert_eq!(full.virtual_masses.len(), 2);
        assert!(full.virtual_masses.iter().all(|p| p.len() == 2));
        assert!(
            full.regulator_prescription
                .contains("fixed-positive-T-and-eta")
        );
        let unsafe_sum = vec![
            vec![Atom::one(), Atom::one(), Atom::one()],
            vec![Atom::zero(), Atom::zero(), Atom::one()],
        ];
        assert!(
            matches!(native_full_kinematics(&unsafe_sum,&[0,1],2,1,Default::default()),Err(Error::Unsupported(ref m)) if m.contains("nonpositive coefficient"))
        );
        // Exactly 2 U +1 pair +2 diagonal +4 virtual-mass monomials.
        assert!(
            matches!(native_full_kinematics(&rows,&[0,1],2,1,ParametricEndpointBudget{polynomial_terms:8,..Default::default()}),Err(Error::Limit(ref m)) if m.contains("total"))
        );
    }
    #[test]
    fn generic_one_virtual_multipair_query_covers_every_prism_triple_cut_support() {
        let input = serde_json::from_str::<DensityInput>(include_str!(
            "../../examples/finite_density/triangular_prism.json"
        ))
        .unwrap()
        .prepare()
        .unwrap();
        let family = input
            .occupied_cut(&[1, 5, 7], 1024)
            .unwrap()
            .at_physical_masses();
        let all = one_virtual_compact_all_supports(&input, &family, Default::default()).unwrap();
        assert_eq!(all.len(), 8);
        assert!(all[0].is_unrestricted_virtual_polynomial());
        for evidence in &all {
            evidence
                .replay_for(&input, &family, Default::default())
                .unwrap();
            assert_eq!(evidence.pure_transfers.len(), 3);
            let bounds = evidence
                .bounds_for_indices(
                    &vec![1; evidence.active_slots.len()],
                    &[1, 1, 1],
                    &[1, 0, 0],
                )
                .unwrap();
            assert!(
                bounds
                    .iter()
                    .all(|b| b.twice_dimension_coefficient * 11 + 2 * b.constant > 0)
            );
            if let Some(full) = evidence.full_kinematics() {
                assert_eq!(full.pairs.len(), 3);
                assert_eq!(full.diagonal_masses.len(), 3);
                assert_eq!(full.virtual_masses.len(), evidence.active_slots.len());
                assert_eq!(evidence.pair_q_upper.len(), 3);
            }
        }
        if let Some(directory) = std::env::var_os("RUSTFLOW_PARAMETRIC_EVIDENCE_REPORT") {
            let directory = std::path::PathBuf::from(directory);
            std::fs::create_dir_all(&directory).unwrap();
            let path = directory.join("prism-cuts-1-5-7-one-virtual-supports.json");
            serde_json::to_writer(std::io::BufWriter::new(std::fs::File::create(&path).unwrap()),&serde_json::json!({
                "scope":"query-only full off-null one-virtual/multipair evidence; no physical admission or period value",
                "source_blake3":blake3::hash(include_bytes!("parametric_endpoint.rs")).to_hex().to_string(),
                "supports":all,
            })).unwrap();
            eprintln!("Saved multipair evidence: {}", path.display());
        }
        assert!(
            one_virtual_compact_all_supports(
                &input,
                &family,
                ParametricEndpointBudget {
                    supports: 7,
                    ..Default::default()
                }
            )
            .is_err()
        );
        let full = all.last().unwrap();
        let largest_support = all
            .iter()
            .filter_map(|e| e.full_kinematics())
            .map(|f| full_kinematic_term_count(f).unwrap())
            .max()
            .unwrap();
        assert!(matches!(
            one_virtual_compact_all_supports(
                &input,
                &family,
                ParametricEndpointBudget {
                    polynomial_terms: largest_support,
                    ..Default::default()
                }
            ),
            Err(Error::Limit(ref message)) if message.contains("aggregate polynomial")
        ));
        assert!(
            full.bounds_for_indices(&[1, 1, 1], &[1, 1], &[1, 0, 0])
                .is_err()
        );
        assert!(
            full.bounds_for_indices(&[1, 0, 1], &[1, 1, 1], &[1, 0, 0])
                .is_err()
        );
        let pair_family = input
            .occupied_cut(&[1, 5], 1024)
            .unwrap()
            .at_physical_masses();
        assert!(
            one_virtual_compact_support_evidence(&input, &pair_family, &[], Default::default())
                .is_err()
        );
        // The query still cannot act as a sealed physical flow permit.
        let shifted = (0..family.physical_slots())
            .filter(|s| ![1, 5, 7].contains(s))
            .collect::<Vec<_>>();
        assert!(
            super::super::massless_endpoint::MasslessFlowEvidence::new(
                &input,
                &family,
                &shifted,
                Default::default()
            )
            .is_err()
        );
    }
}
/// Exact finite inputs to a family-bound support query. This is not a sealed
/// physical flow permit and does not discard coefficient conditions.
#[derive(Clone, Debug, Serialize)]
pub struct ParametricLabelEnvelope {
    pub label: Vec<i16>,
    pub active_virtual_slots: Vec<usize>,
    pub positive_virtual_powers: Vec<u32>,
    pub positive_transfer_powers: Vec<(usize, u32)>,
    pub polynomial_momentum_rank: u32,
    pub shell_and_upper_jets: Vec<u32>,
    pub missing_required_cut: bool,
    pub has_lower_contact: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct ParametricLabelBoundQuery {
    pub envelope: ParametricLabelEnvelope,
    pub eta_bounds: Vec<ParametricEtaBound>,
    pub compact_bounds: Vec<CompactEndpointDegree>,
}
#[derive(Clone, Debug, Serialize)]
pub struct OneVirtualLabelBoundQuery {
    pub envelope: ParametricLabelEnvelope,
    pub compact_bounds: Vec<CompactEndpointDegree>,
}
fn sum_i64(values: &[u32]) -> Result<i64> {
    values
        .iter()
        .try_fold(0_i64, |s, &n| checked_add(s, i64::from(n)))
}
fn sum_u32(values: &[u32]) -> Result<u32> {
    values
        .iter()
        .try_fold(0_u32, |s, &n| s.checked_add(n))
        .ok_or_else(|| limit("finite jet count overflow"))
}
fn validate_envelope_support(envelope: &ParametricLabelEnvelope, active: &[usize]) -> Result<()> {
    if envelope.missing_required_cut || envelope.has_lower_contact {
        return Err(unsupported(
            "cut-vanishing and lower-contact labels need their separate distribution proof",
        ));
    }
    if envelope.active_virtual_slots != active {
        return Err(invalid(
            "finite label active support differs from the retained native proof",
        ));
    }
    Ok(())
}

fn exact_momentum_degree(factor: &Atom, family: &OccupiedCutFamily) -> Result<u32> {
    let coordinates = family.coordinates();
    let scalar_count = family
        .loops()
        .checked_mul(family.loops() + 1)
        .and_then(|x| x.checked_div(2))
        .ok_or_else(|| limit("Gram coordinate count overflow"))?;
    if coordinates.len() != scalar_count + family.loops() {
        return Err(invalid("occupied momentum-coordinate arity mismatch"));
    }
    let terms = exact_coefficient_list(factor, coordinates)?;
    let mut largest = 0_u32;
    for (monomial, coefficient) in terms {
        for coordinate in coordinates {
            let AtomView::Var(v) = coordinate.as_view() else {
                return Err(invalid("occupied coordinate is not a scalar symbol"));
            };
            if !coefficient.derivative(v.get_symbol()).is_zero() {
                return Err(unsupported(
                    "nonpolynomial dependence in an endpoint numerator factor",
                ));
            }
        }
        let exponents = powers(&monomial, coordinates)?;
        let mut degree = 0_u32;
        for (i, n) in exponents.into_iter().enumerate() {
            let n = u32::try_from(n).map_err(|_| {
                unsupported("inverse momentum monomial in an endpoint numerator factor")
            })?;
            degree = degree
                .checked_add(
                    n.checked_mul(if i < scalar_count { 2 } else { 1 })
                        .ok_or_else(|| limit("momentum degree overflow"))?,
                )
                .ok_or_else(|| limit("momentum degree overflow"))?;
        }
        largest = largest.max(degree);
    }
    Ok(largest)
}

/// `pure_transfer_slots` comes from a validated exact family-bound query, not
/// an arbitrary classification supplied by a flow caller. Before this envelope
/// can authorize anything, its constructor must remain private to that owner.
fn parametric_label_envelope(
    family: &OccupiedCutFamily,
    label: &crate::Integral,
    pure_transfer_slots: &[usize],
) -> Result<ParametricLabelEnvelope> {
    if label.0.len() != family.factors().len() {
        return Err(invalid("parametric endpoint label arity mismatch"));
    }
    // Validate every role first. An absent cut cannot hide an inverse
    // completion or negative occupation. Physical coefficient conditions
    // remain owned by the reduction/evaluation layer and are not erased here.
    if (family.physical_slots()..family.input_slots()).any(|s| label.0[s] > 0)
        || family
            .shells()
            .iter()
            .any(|s| label.0[s.upper_slot] < 0 || label.0[s.lower_slot] < 0)
    {
        return Err(unsupported(
            "endpoint labels require polynomial completions and nonnegative occupations",
        ));
    }
    let cuts = family
        .shells()
        .iter()
        .map(|s| s.physical_slot)
        .collect::<Vec<_>>();
    let mut active_virtual_slots = Vec::new();
    let mut positive_virtual_powers = Vec::new();
    let mut positive_transfer_powers = Vec::new();
    for slot in 0..family.physical_slots() {
        if cuts.contains(&slot) {
            continue;
        }
        let n = u32::from(label.0[slot].max(0) as u16);
        if pure_transfer_slots.contains(&slot) {
            positive_transfer_powers.push((slot, n));
        } else if n > 0 {
            active_virtual_slots.push(slot);
            positive_virtual_powers.push(n);
        }
    }
    let mut polynomial_momentum_rank = 0_u32;
    for (slot, &index) in label.0[..family.input_slots()].iter().enumerate() {
        if index >= 0 || cuts.contains(&slot) {
            continue;
        }
        let multiplicity = i32::from(index).unsigned_abs();
        let degree = exact_momentum_degree(&family.factors()[slot], family)?;
        polynomial_momentum_rank = polynomial_momentum_rank
            .checked_add(
                degree
                    .checked_mul(multiplicity)
                    .ok_or_else(|| limit("endpoint polynomial rank overflow"))?,
            )
            .ok_or_else(|| limit("endpoint polynomial rank overflow"))?;
    }
    let shell_and_upper_jets = family
        .shells()
        .iter()
        .map(|s| {
            let cut = u32::from(label.0[s.physical_slot].saturating_sub(1).max(0) as u16);
            let upper = u32::from(label.0[s.upper_slot].saturating_sub(1).max(0) as u16);
            cut.checked_add(upper)
                .ok_or_else(|| limit("shell/occupation jet count overflow"))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(ParametricLabelEnvelope {
        label: label.0.clone(),
        active_virtual_slots,
        positive_virtual_powers,
        positive_transfer_powers,
        polynomial_momentum_rank,
        shell_and_upper_jets,
        missing_required_cut: family
            .shells()
            .iter()
            .any(|s| label.0[s.physical_slot] <= 0),
        has_lower_contact: family.shells().iter().any(|s| label.0[s.lower_slot] > 0),
    })
}

#[cfg(test)]
mod label_envelope_tests {
    use super::*;
    use crate::finite_density::DensityInput;
    fn prism() -> PreparedDensityInput {
        serde_json::from_str::<DensityInput>(include_str!(
            "../../examples/finite_density/triangular_prism.json"
        ))
        .unwrap()
        .prepare()
        .unwrap()
    }
    fn bulk(family: &OccupiedCutFamily) -> crate::Integral {
        crate::Integral(
            (0..family.factors().len())
                .map(|i| i16::from(i < family.physical_slots()))
                .collect(),
        )
    }
    fn energy_slot(family: &OccupiedCutFamily) -> usize {
        (family.physical_slots()..family.input_slots())
            .find(|&s| exact_momentum_degree(&family.factors()[s], family).unwrap() == 1)
            .unwrap()
    }
    #[test]
    fn actual_label_degree_distinguishes_energy_quadratic_and_inactive_transfer_powers() {
        let input = prism();
        let family = input
            .occupied_cut(&[1, 5], 1024)
            .unwrap()
            .at_physical_masses();
        let empty =
            two_occupied_support_evidence(&input, &family, &[], Default::default()).unwrap();
        let active = (0..family.physical_slots())
            .filter(|s| {
                !empty.cuts.contains(s) && !empty.pure_transfer_slots.iter().any(|(p, _)| p == s)
            })
            .collect::<Vec<_>>();
        let proof =
            two_occupied_support_evidence(&input, &family, &active, Default::default()).unwrap();
        let mut label = bulk(&family);
        let transfer = proof.pure_transfer_slots[0].0;
        let energy = energy_slot(&family);
        label.0[transfer] = -2;
        label.0[energy] = -3;
        label.0[family.shells()[0].physical_slot] = 2;
        label.0[family.shells()[0].upper_slot] = 2;
        let query = proof.bounds_for_label(&family, &label).unwrap();
        assert_eq!(query.envelope.polynomial_momentum_rank, 7);
        assert_eq!(query.envelope.shell_and_upper_jets, vec![2, 0]);
        assert_eq!(query.envelope.positive_transfer_powers, vec![(transfer, 0)]);
        assert_eq!(query.envelope.positive_virtual_powers, vec![1; 6]);
        assert_eq!(
            query
                .compact_bounds
                .iter()
                .find(|b| b.stratum.starts_with("combined"))
                .unwrap()
                .constant,
            -8
        );
        assert!(
            query
                .eta_bounds
                .iter()
                .all(|b| b.regular && b.constant == 0 || b.twice_dimension_coefficient > 0)
        );
        let mut inactive = label.clone();
        inactive.0[active[0]] = 0;
        assert!(
            matches!(proof.bounds_for_label(&family,&inactive),Err(Error::InvalidInput(ref m))if m.contains("active support"))
        );
        let mut absent = label.clone();
        absent.0[family.shells()[0].physical_slot] = 0;
        assert!(
            proof
                .label_envelope(&family, &absent)
                .unwrap()
                .missing_required_cut
        );
        assert!(proof.bounds_for_label(&family, &absent).is_err());
        absent.0[energy] = 1;
        assert!(proof.label_envelope(&family, &absent).is_err());
        let mut negative = label.clone();
        negative.0[family.shells()[0].upper_slot] = -1;
        assert!(proof.label_envelope(&family, &negative).is_err());
        let mut lower = label;
        lower.0[family.shells()[0].lower_slot] = 1;
        assert!(
            proof
                .label_envelope(&family, &lower)
                .unwrap()
                .has_lower_contact
        );
        assert!(proof.bounds_for_label(&family, &lower).is_err());
        let other = input
            .occupied_cut(&[1, 7], 1024)
            .unwrap()
            .at_physical_masses();
        assert!(proof.label_envelope(&other, &bulk(&other)).is_err());
    }
    #[test]
    fn multipair_finite_label_query_retains_each_transfer_and_original_energy_degree() {
        let input = prism();
        let family = input
            .occupied_cut(&[1, 5, 7], 1024)
            .unwrap()
            .at_physical_masses();
        let proofs = one_virtual_compact_all_supports(&input, &family, Default::default()).unwrap();
        let proof = proofs.last().unwrap();
        let mut label = bulk(&family);
        label.0[proof.pure_transfers[0].slot] = 0;
        label.0[family.shells()[0].physical_slot] = 2;
        label.0[energy_slot(&family)] = -2;
        let query = proof.bounds_for_label(&family, &label).unwrap();
        assert_eq!(query.envelope.polynomial_momentum_rank, 2);
        assert_eq!(
            query
                .envelope
                .positive_transfer_powers
                .iter()
                .map(|(_, p)| p)
                .sum::<u32>(),
            2
        );
        assert_eq!(query.envelope.shell_and_upper_jets, vec![1, 0, 0]);
        assert!(
            query
                .compact_bounds
                .iter()
                .all(|b| b.twice_dimension_coefficient * 11 + 2 * b.constant > 0)
        );
        let angle = query
            .compact_bounds
            .iter()
            .find(|b| b.stratum.starts_with("Holder"))
            .unwrap();
        assert_eq!(angle.constant, -4);
    }
}
