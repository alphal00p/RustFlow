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

const VERSION: &str =
    "hepkit-positive-joint-eta-sectors-v1;query-only;UV-meromorphic-before-high-D";

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
            if !self.u.is_empty() || !self.v.is_empty() || !self.charts.is_empty() {
                return Err(invalid("rank-deficient evidence carries parametric charts"));
            }
            return Ok(());
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
    };
    if rank < h {
        return Ok(evidence);
    }
    let orders = permutations(active_slots.len(), budget)?;
    let loops = (0..h).map(|i| symbol(&format!("k{i}"))).collect::<Vec<_>>();
    let external = vec![symbol("q0"), symbol("q1")];
    let invariant = symbol("h");
    let eta = symbol("eta");
    let kin = Kinematics::new()
        .with_momenta(loops.iter().chain(&external).cloned())
        .map_err(native_error)?
        .with_mass_squared(&external[0], Atom::zero())
        .map_err(native_error)?
        .with_mass_squared(&external[1], Atom::zero())
        .map_err(native_error)?
        .with_scalar_product(&external[0], &external[1], &invariant / 2)
        .map_err(native_error)?;
    let momenta = external.iter().chain(&loops).cloned().collect::<Vec<_>>();
    let denoms = active_slots
        .iter()
        .map(|&s| {
            let q = routed[s]
                .iter()
                .zip(&momenta)
                .map(|(a, k)| a * k)
                .sum::<Atom>();
            Ok(kin.scalar_product(&q, &q).map_err(native_error)? - &eta)
        })
        .collect::<Result<Vec<_>>>()?;
    let native = NativeFamily::new(loops, external, denoms, &kin).map_err(native_error)?;
    let parameters = (0..active_slots.len())
        .map(|i| symbol(&format!("alpha{i}")))
        .collect::<Vec<_>>();
    let (u, f) = native.symanzik(&parameters).map_err(native_error)?;
    let AtomView::Var(iv) = invariant.as_view() else {
        unreachable!()
    };
    let v = f.derivative(iv.get_symbol()).expand();
    let sum = parameters.iter().cloned().sum::<Atom>();
    if !(f - &eta * &u * sum - &invariant * &v).expand().is_zero() {
        return Err(unsupported("native F is not eta U sum(alpha)+h V"));
    }
    evidence.u = positive_polynomial(&u, &parameters, budget)?;
    evidence.v = positive_polynomial(&v, &parameters, budget)?;
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
