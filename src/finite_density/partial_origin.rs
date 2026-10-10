//! Finite-positive-eta origin capability for virtual rank at most two.
//! No endpoint, infinity-boundary, integration or source-factory dispatch.
//! Limits bound adapter work and retained output, not internal CAS work.
//! Cooperative cancellation and conservative input envelopes precede native work.
//! Numerical admission belongs to a separately bound continuation consumer.
use super::guarded::{IndexBounds, IndexDomain, IndexRole};
use super::preparation::WeightedSourceOptions;
use super::{PreparedDensityInput, geometry::OccupiedCutFamily};
use crate::algebra::{matmul, rref};
use crate::coefficient::{exact_coefficient_list, powers};
use crate::{Error, Result, RunContext};
#[path = "partial_continuation.rs"]
pub(crate) mod continuation;
#[path = "partial_endpoint.rs"]
pub(crate) mod endpoint;
#[path = "partial_resources.rs"]
pub(crate) mod resources;
use feynkit_graph::IntegralFamily as NativeFamily;
use feynkit_kinematics::Kinematics;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use symbolica::prelude::*;

const VERSION: &str =
    "partial-origin-rank-two-v2;universal-finite-label-JN;no-endpoint-no-infinity";
const PRESCRIPTION: &str = "independent-native-mass-minus-i-epsilon;future-square-root-real-spatial;epsilon-removed-at-fixed-positive-eta-and-T;real-Fermi-T-limit;then-meromorphic-D;original-polynomial-fixed-in-mass-jets";
fn invalid(m: &str) -> Error {
    Error::InvalidInput(format!("partial origin: {m}"))
}
fn unsupported(m: &str) -> Error {
    Error::Unsupported(format!("partial origin: {m}"))
}
fn limit(m: &str) -> Error {
    Error::Limit(format!("partial origin: {m}"))
}
fn atom(name: &str) -> Atom {
    Atom::var(symbol!(format!("rustflow_partial_origin::{name}")))
}
fn add(a: u64, b: u64) -> Result<u64> {
    a.checked_add(b).ok_or_else(|| limit("bound overflow"))
}
fn mul(a: u64, b: u64) -> Result<u64> {
    a.checked_mul(b).ok_or_else(|| limit("bound overflow"))
}
fn bound(lo: Option<i64>, hi: Option<i64>) -> Result<IndexBounds> {
    IndexBounds::new(lo, hi).map_err(|e| invalid(&e.to_string()))
}
fn rat(a: &Atom) -> Result<Rational> {
    Rational::try_from(a.as_view()).map_err(|_| unsupported("nonrational row/coefficient"))
}
fn zero(a: Atom) -> bool {
    a.expand().together().cancel().is_zero()
}
fn derivative(a: &Atom, v: &Atom) -> Atom {
    let AtomView::Var(s) = v.as_view() else {
        unreachable!()
    };
    a.derivative(s.get_symbol()).expand()
}

#[derive(Clone, Copy, Debug, Serialize)]
pub(crate) struct PartialOriginBudget {
    pub physical_rows: usize,
    pub loops: usize,
    pub input_text_bytes: usize,
    pub geometry_terms: usize,
    pub supports: usize,
    pub charts: usize,
    /// Aggregate retained exact polynomial terms, including all support/chart replays.
    pub polynomial_terms: usize,
    /// Adapter events only, NOT an upper bound on internal HEPKit/Symbolica work.
    /// Structural preflight bounds rank/rows/input size; cancellation is cooperative.
    pub operations: usize,
    pub rational_digits: usize,
    pub degree: u32,
    pub expressions: resources::ExpressionLimits,
    pub aggregate_expression_bytes: usize,
    pub label_axes: usize,
    pub index_magnitude: u64,
}
impl Default for PartialOriginBudget {
    fn default() -> Self {
        Self {
            physical_rows: 8,
            loops: 8,
            input_text_bytes: 65536,
            geometry_terms: 1024,
            supports: 256,
            charts: 4096,
            polynomial_terms: 200000,
            operations: 2000000,
            rational_digits: 1024,
            degree: 1024,
            expressions: Default::default(),
            aggregate_expression_bytes: 16 * 1024 * 1024,
            label_axes: 128,
            index_magnitude: 65536,
        }
    }
}
struct Meter<'a> {
    run: &'a RunContext,
    budget: PartialOriginBudget,
    operations: usize,
    terms: usize,
    bytes: usize,
    charts: usize,
}
impl Meter<'_> {
    fn tick(&mut self, n: usize) -> Result<()> {
        self.run.cancellation.check()?;
        self.operations = self
            .operations
            .checked_add(n)
            .ok_or_else(|| limit("operation overflow"))?;
        if self.operations > self.budget.operations {
            return Err(limit("operation budget"));
        }
        Ok(())
    }
    fn preflight(&mut self, a: &Atom) -> Result<()> {
        self.bytes = self
            .bytes
            .checked_add(a.as_view().get_byte_size())
            .ok_or_else(|| limit("aggregate expression bytes overflow"))?;
        if self.bytes > self.budget.aggregate_expression_bytes {
            return Err(limit("aggregate expression bytes budget"));
        }
        let nodes = resources::expression(a, self.budget.expressions, self.run)?;
        self.tick(nodes)
    }
    fn terms(&mut self, n: usize) -> Result<()> {
        self.terms = self
            .terms
            .checked_add(n)
            .ok_or_else(|| limit("term overflow"))?;
        if self.terms > self.budget.polynomial_terms {
            return Err(limit("aggregate polynomial budget"));
        }
        self.tick(n)
    }
    fn rational(&mut self, r: &Rational) -> Result<()> {
        self.tick(1)?;
        if r.to_string().len() > self.budget.rational_digits {
            return Err(limit("rational size budget"));
        }
        Ok(())
    }
}
#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
struct Monomial {
    exponents: Vec<u32>,
    coefficient: String,
}
fn polynomial(a: &Atom, vars: &[Atom], positive: bool, m: &mut Meter) -> Result<Vec<Monomial>> {
    m.preflight(a)?;
    let list = exact_coefficient_list(a, vars)?;
    m.tick(1)?;
    m.terms(list.len())?;
    list.into_iter()
        .map(|(monomial, c)| {
            let c = rat(&c)?;
            m.rational(&c)?;
            if positive && c <= 0 {
                return Err(unsupported("nonpositive channel/unit coefficient"));
            }
            let exponents = powers(&monomial, vars)?
                .into_iter()
                .map(|e| u32::try_from(e).map_err(|_| unsupported("negative polynomial exponent")))
                .collect::<Result<Vec<_>>>()?;
            if exponents.iter().any(|e| *e > m.budget.degree) {
                return Err(limit("polynomial degree budget"));
            }
            Ok(Monomial {
                exponents,
                coefficient: c.to_string(),
            })
        })
        .collect()
}
fn strings(rows: &[Vec<Atom>]) -> Vec<Vec<String>> {
    rows.iter()
        .map(|r| r.iter().map(Atom::to_canonical_string).collect())
        .collect()
}
fn preflight_family(f: &OccupiedCutFamily, b: PartialOriginBudget, run: &RunContext) -> Result<()> {
    run.cancellation.check()?;
    if f.loops() > b.loops
        || f.physical_slots() > b.physical_rows
        || f.factors().len() > b.label_axes
    {
        return Err(limit("family preflight"));
    }
    let mut nodes = 0usize;
    let mut bytes = 0usize;
    for a in f
        .coordinates()
        .iter()
        .chain(f.factors())
        .chain(f.inverse_routing().iter().flatten())
        .chain(f.targets().iter().flat_map(|t| t.values()))
    {
        bytes = bytes
            .checked_add(a.as_view().get_byte_size())
            .ok_or_else(|| limit("family bytes overflow"))?;
        if bytes > b.aggregate_expression_bytes {
            return Err(limit("family bytes budget"));
        }
        nodes = nodes
            .checked_add(resources::expression(a, b.expressions, run)?)
            .ok_or_else(|| limit("family operation overflow"))?;
        if nodes > b.operations {
            return Err(limit("family operation preflight"));
        }
    }
    Ok(())
}
fn family_signature(f: &OccupiedCutFamily) -> Vec<String> {
    let mut s = vec![format!(
        "loops={};physical={};input={};det={}",
        f.loops(),
        f.physical_slots(),
        f.input_slots(),
        f.routing_determinant()
    )];
    s.extend(
        f.coordinates()
            .iter()
            .chain(f.factors())
            .map(Atom::to_canonical_string),
    );
    s.extend(f.roles().iter().map(|r| format!("{r:?}")));
    for row in f.inverse_routing() {
        s.extend(row.iter().map(Atom::to_canonical_string));
    }
    for x in f.shells() {
        s.push(format!(
            "shell={}:{}:{}:{}:{}:{}",
            x.loop_index,
            x.physical_slot,
            x.upper_slot,
            x.lower_slot,
            x.mass_squared,
            x.chemical_potential
        ));
    }
    for (t, terms) in f.targets().iter().enumerate() {
        for (i, c) in terms {
            s.push(format!("target={t}:{:?}:{}", i.0, c.to_canonical_string()));
        }
    }
    s
}
fn rank(rows: &[Vec<Atom>], h: usize) -> usize {
    if rows.is_empty() {
        0
    } else {
        rref(rows.to_vec()).1.into_iter().filter(|p| *p < h).count()
    }
}
fn null_vector(rows: &[Vec<Atom>], h: usize) -> Result<Vec<Atom>> {
    let (reduced, pivots) = rref(rows.to_vec());
    let free = (0..h)
        .find(|j| !pivots.contains(j))
        .ok_or_else(|| invalid("missing null direction"))?;
    let mut v = vec![Atom::zero(); h];
    v[free] = Atom::one();
    for (i, &p) in pivots.iter().enumerate() {
        if p < h {
            v[p] = -reduced[i][free].clone();
        }
    }
    if rows
        .iter()
        .any(|r| !zero(r.iter().zip(&v).map(|(a, b)| a * b).sum()))
    {
        return Err(invalid("null vector replay"));
    }
    Ok(v)
}
fn square(row: &[Atom], coordinates: &[Atom], loops: usize) -> Atom {
    let mut out = Atom::zero();
    let mut n = 0;
    for i in 0..loops {
        for j in i..loops {
            out += &row[i] * &row[j] * &coordinates[n] * if i == j { 1 } else { 2 };
            n += 1;
        }
    }
    out.expand()
}
fn affine(
    rows: &[Vec<Atom>],
    shifted: &BTreeSet<usize>,
    cuts: &BTreeSet<usize>,
    k: usize,
    h: usize,
    m: &mut Meter,
) -> Result<Vec<Vec<Atom>>> {
    let mut t = vec![vec![Atom::zero(); k]; h];
    for c in 0..k {
        let matrix = rows
            .iter()
            .enumerate()
            .filter(|(s, _)| !cuts.contains(s) && !shifted.contains(s))
            .map(|(_, r)| {
                let mut row = r[k..].to_vec();
                row.push(-r[c].clone());
                row
            })
            .collect::<Vec<_>>();
        m.tick(matrix.len() * h.max(1))?;
        let (reduced, pivots) = rref(matrix);
        m.tick(1)?;
        if pivots.contains(&h) {
            return Err(unsupported(
                "no common affine translation for unshifted rows",
            ));
        }
        for (i, &p) in pivots.iter().enumerate() {
            t[p][c] = reduced[i][h].clone();
        }
    }
    for (s, r) in rows
        .iter()
        .enumerate()
        .filter(|(s, _)| !cuts.contains(s) && !shifted.contains(s))
    {
        for c in 0..k {
            m.tick(h + 1)?;
            if !zero(&r[c] + r[k..].iter().zip(&t).map(|(a, b)| a * &b[c]).sum::<Atom>()) {
                return Err(invalid(&format!("affine replay slot {s}")));
            }
        }
    }
    Ok(t)
}

#[derive(Clone, Debug, Serialize)]
struct FullKinematics {
    u: Vec<Monomial>,
    pairs: Vec<([usize; 2], Vec<Monomial>)>,
    diagonal: Vec<Vec<Monomial>>,
    line_masses: Vec<Vec<Monomial>>,
    eta: Vec<Monomial>,
    native_u: String,
    native_f: String,
}
struct Gaussian {
    parameters: Vec<Atom>,
    u: Atom,
    means: Vec<Vec<Atom>>,
    v: Vec<Vec<Atom>>,
    native: FullKinematics,
}
/// Same native construction as the existing full-off-null query, with the actual partial eta map.
fn gaussian(
    rows: &[Vec<Atom>],
    active: &[usize],
    shifted: &BTreeSet<usize>,
    k: usize,
    h: usize,
    m: &mut Meter,
) -> Result<Gaussian> {
    m.tick(1)?;
    resources::dense_envelope(active.len(), 2 * (h + 1), m.budget.polynomial_terms)?;
    let loops = (0..h).map(|i| atom(&format!("k{i}"))).collect::<Vec<_>>();
    let external = (0..k).map(|i| atom(&format!("q{i}"))).collect::<Vec<_>>();
    let shell = (0..k)
        .map(|i| atom(&format!("shell{i}")))
        .collect::<Vec<_>>();
    let masses = (0..active.len())
        .map(|i| atom(&format!("mass{i}")))
        .collect::<Vec<_>>();
    let parameters = (0..active.len())
        .map(|i| atom(&format!("alpha{i}")))
        .collect::<Vec<_>>();
    let eta = atom("eta");
    let mut kin = Kinematics::new()
        .with_momenta(loops.iter().chain(&external).cloned())
        .map_err(|e| unsupported(&e.to_string()))?;
    for (q, a) in external.iter().zip(&shell) {
        kin = kin
            .with_mass_squared(q, a.clone())
            .map_err(|e| unsupported(&e.to_string()))?;
    }
    let mut pair_vars = Vec::new();
    for i in 0..k {
        for j in i + 1..k {
            let t = atom(&format!("dot_{i}_{j}"));
            kin = kin
                .with_scalar_product(&external[i], &external[j], &t / 2)
                .map_err(|e| unsupported(&e.to_string()))?;
            pair_vars.push(([i, j], t));
        }
    }
    let momenta = external.iter().chain(&loops).collect::<Vec<_>>();
    let factors = active
        .iter()
        .zip(&masses)
        .map(|(&slot, mass)| {
            let p = rows[slot]
                .iter()
                .zip(&momenta)
                .map(|(a, b)| a * (*b))
                .sum::<Atom>();
            Ok(kin
                .scalar_product(&p, &p)
                .map_err(|e| unsupported(&e.to_string()))?
                - mass
                - if shifted.contains(&slot) {
                    eta.clone()
                } else {
                    Atom::zero()
                })
        })
        .collect::<Result<Vec<_>>>()?;
    m.tick(active.len() * active.len() * h.max(1))?;
    for f in &factors {
        m.preflight(f)?;
    }
    let native = NativeFamily::new(loops, external, factors, &kin)
        .map_err(|e| unsupported(&e.to_string()))?;
    m.tick(1)?;
    let (u, f) = native
        .symanzik(&parameters)
        .map_err(|e| unsupported(&e.to_string()))?;
    m.preflight(&u)?;
    m.preflight(&f)?;
    let a = (0..h)
        .map(|i| {
            (0..h)
                .map(|j| {
                    active
                        .iter()
                        .zip(&parameters)
                        .map(|(&s, x)| x * &rows[s][k + i] * &rows[s][k + j])
                        .sum::<Atom>()
                        .expand()
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let b = (0..h)
        .map(|i| {
            (0..k)
                .map(|j| {
                    active
                        .iter()
                        .zip(&parameters)
                        .map(|(&s, x)| x * &rows[s][k + i] * &rows[s][j])
                        .sum::<Atom>()
                        .expand()
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let expected_u = if h == 1 {
        a[0][0].clone()
    } else {
        (&a[0][0] * &a[1][1] - &a[0][1] * &a[1][0]).expand()
    };
    if !zero(&u - &expected_u) {
        return Err(invalid(
            "HEPKit U differs from complete Gaussian determinant",
        ));
    }
    let adj = if h == 1 {
        vec![vec![Atom::one()]]
    } else {
        vec![
            vec![a[1][1].clone(), -a[0][1].clone()],
            vec![-a[1][0].clone(), a[0][0].clone()],
        ]
    };
    let means = matmul(&adj, &b);
    let mut v = vec![vec![Atom::zero(); k]; k];
    for i in 0..k {
        for j in 0..k {
            let c = active
                .iter()
                .zip(&parameters)
                .map(|(&s, x)| x * &rows[s][i] * &rows[s][j])
                .sum::<Atom>();
            v[i][j] = (&u * c - (0..h).map(|r| &b[r][i] * &means[r][j]).sum::<Atom>()).expand();
        }
    }
    let s = active
        .iter()
        .zip(&parameters)
        .filter(|(slot, _)| shifted.contains(slot))
        .map(|(_, x)| x.clone())
        .sum::<Atom>();
    let etac = (&u * &s).expand();
    if !zero(derivative(&f, &eta) - &etac) {
        return Err(invalid("HEPKit actual eta mass map mismatch"));
    }
    let mut rebuilt = &eta * &etac;
    let mut pairs = Vec::new();
    let mut diagonals = Vec::new();
    let mut line = Vec::new();
    for (pair, t) in pair_vars {
        let c = derivative(&f, &t);
        if !zero(&c + &v[pair[0]][pair[1]]) {
            return Err(invalid("off-null pair Gaussian replay"));
        }
        rebuilt += &t * &c;
        pairs.push((pair, polynomial(&c, &parameters, true, m)?));
    }
    for (i, t) in shell.iter().enumerate() {
        let c = derivative(&f, t);
        if !zero(&c + &v[i][i]) {
            return Err(invalid("off-null diagonal Gaussian replay"));
        }
        rebuilt += t * &c;
        diagonals.push(polynomial(&c, &parameters, false, m)?);
    }
    for (t, x) in masses.iter().zip(&parameters) {
        let c = derivative(&f, t);
        if !zero(&c - &u * x) {
            return Err(invalid("independent line regulator replay"));
        }
        rebuilt += t * &c;
        line.push(polynomial(&c, &parameters, true, m)?);
    }
    if !zero(&f - rebuilt) {
        return Err(invalid("HEPKit F has unretained nonlinear kinematics"));
    }
    Ok(Gaussian {
        parameters: parameters.clone(),
        u: u.clone(),
        means,
        v,
        native: FullKinematics {
            u: polynomial(&u, &parameters, true, m)?,
            pairs,
            diagonal: diagonals,
            line_masses: line,
            eta: polynomial(&etac, &parameters, true, m)?,
            native_u: u.to_canonical_string(),
            native_f: f.to_canonical_string(),
        },
    })
}

#[derive(Clone, Debug, Serialize)]
struct Chart {
    primary: usize,
    pivot: Option<usize>,
    parameter_map: Vec<String>,
    jacobian_x_power: usize,
    unit_lower_bound: String,
    u_unit: Vec<Monomial>,
    mean_units: Vec<Vec<Vec<Monomial>>>,
    v_units: Vec<Vec<Vec<Monomial>>>,
}
fn substitute(a: &Atom, map: &BTreeMap<Atom, Atom>) -> Atom {
    crate::family::substitute(a, map).expand()
}
fn charts(
    rows: &[Vec<Atom>],
    active: &[usize],
    k: usize,
    h: usize,
    g: &Gaussian,
    m: &mut Meter,
) -> Result<Vec<Chart>> {
    let mut out = Vec::new();
    for i in 0..active.len() {
        let pivots = if h == 1 {
            vec![None]
        } else {
            (0..active.len())
                .filter(|&j| {
                    !zero(
                        &rows[active[i]][k] * &rows[active[j]][k + 1]
                            - &rows[active[i]][k + 1] * &rows[active[j]][k],
                    )
                })
                .map(Some)
                .collect()
        };
        for pivot in pivots {
            m.tick(1)?;
            m.charts = m
                .charts
                .checked_add(1)
                .ok_or_else(|| limit("chart overflow"))?;
            if m.charts > m.budget.charts {
                return Err(limit("complete chart budget"));
            }
            let variables = (0..active.len() - 1)
                .map(|j| atom(&format!("z{j}")))
                .collect::<Vec<_>>();
            let mut map = BTreeMap::new();
            let mut position = if h == 1 { 0 } else { 1 };
            let mut transverse = 0usize;
            let mut parameter_map = Vec::new();
            for (e, &slot) in active.iter().enumerate() {
                let value = if e == i {
                    Atom::one()
                } else if Some(e) == pivot {
                    variables[0].clone()
                } else {
                    let value = variables[position].clone();
                    position += 1;
                    if h == 2
                        && !zero(
                            &rows[active[i]][k] * &rows[slot][k + 1]
                                - &rows[active[i]][k + 1] * &rows[slot][k],
                        )
                    {
                        &variables[0] * value
                    } else {
                        value
                    }
                };
                if h == 2
                    && !zero(
                        &rows[active[i]][k] * &rows[slot][k + 1]
                            - &rows[active[i]][k + 1] * &rows[slot][k],
                    )
                {
                    transverse += 1;
                }
                parameter_map.push(value.to_canonical_string());
                map.insert(g.parameters[e].clone(), value);
            }
            let divide = if h == 2 {
                variables[0].clone()
            } else {
                Atom::one()
            };
            let u = (substitute(&g.u, &map) / &divide)
                .together()
                .cancel()
                .expand();
            let units = polynomial(&u, &variables, true, m)?;
            let expected = if let Some(j) = pivot {
                (&rows[active[i]][k] * &rows[active[j]][k + 1]
                    - &rows[active[i]][k + 1] * &rows[active[j]][k])
                    .pow(2)
            } else {
                rows[active[i]][k].clone().pow(2)
            };
            let constant = units
                .iter()
                .find(|t| t.exponents.iter().all(|x| *x == 0))
                .ok_or_else(|| invalid("chart lacks unit constant"))?;
            if constant.coefficient != rat(&expected)?.to_string() {
                return Err(invalid("chart pivot constant mismatch"));
            }
            let mut means = Vec::new();
            for row in &g.means {
                let mut terms = Vec::new();
                for a in row {
                    terms.push(polynomial(
                        &(substitute(a, &map) / &divide).together().cancel().expand(),
                        &variables,
                        false,
                        m,
                    )?);
                }
                means.push(terms);
            }
            let mut v = Vec::new();
            for row in &g.v {
                let mut terms = Vec::new();
                for a in row {
                    terms.push(polynomial(
                        &(substitute(a, &map) / &divide).together().cancel().expand(),
                        &variables,
                        false,
                        m,
                    )?);
                }
                v.push(terms);
            }
            out.push(Chart {
                primary: active[i],
                pivot: pivot.map(|j| active[j]),
                parameter_map,
                jacobian_x_power: transverse.saturating_sub(1),
                unit_lower_bound: constant.coefficient.clone(),
                u_unit: units,
                mean_units: means,
                v_units: v,
            });
        }
    }
    Ok(out)
}

#[derive(Clone, Debug, Serialize)]
enum SupportClass {
    FreeVirtual {
        null_direction: Vec<String>,
    },
    RegulatedVacuumZero,
    FiniteEtaJets {
        shifted_positive: Vec<usize>,
        kinematics_vanish: bool,
    },
}
#[derive(Clone, Debug, Serialize)]
struct Support {
    active: Vec<usize>,
    rank: usize,
    class: SupportClass,
    full: Option<FullKinematics>,
    charts: Vec<Chart>,
}
#[derive(Clone, Debug, Serialize)]
pub(crate) struct PartialOriginCapability {
    version: &'static str,
    prescription: &'static str,
    input_identity: String,
    family_signature: Vec<String>,
    source_options: WeightedSourceOptions,
    cuts: Vec<usize>,
    shifted: Vec<usize>,
    virtual_slots: Vec<usize>,
    pure_compact_slots: Vec<usize>,
    routed_rows: Vec<Vec<String>>,
    translation: Vec<Vec<String>>,
    translated_rows: Vec<Vec<String>>,
    compact_loops: usize,
    virtual_loops: usize,
    schur_bound: String,
    auxiliary_coefficients: Vec<(usize, String)>,
    complete_supports: Vec<Support>,
    factor_degrees: Vec<u64>,
    operations: usize,
    polynomial_terms: usize,
}
#[derive(Clone, Debug, Serialize)]
pub(crate) struct FiniteLabelOriginAudit {
    pub label: Vec<i64>,
    pub active_virtual_slots: Vec<usize>,
    pub classification: String,
    pub positive_virtual_degree: u64,
    pub polynomial_degree: u64,
    pub jets_including_lower: u64,
    pub cut_vanishing: bool,
    pub dimension_witness: String,
    pub subtracted_terms: u64,
    pub remaining_mass_power: String,
}
impl PartialOriginCapability {
    pub(crate) fn new_with_context(
        input: &PreparedDensityInput,
        family: &OccupiedCutFamily,
        shifted: &[usize],
        options: WeightedSourceOptions,
        budget: PartialOriginBudget,
        run: &RunContext,
    ) -> Result<Self> {
        run.cancellation.check()?;
        if budget.physical_rows == 0
            || budget.loops == 0
            || budget.input_text_bytes == 0
            || budget.geometry_terms == 0
            || budget.supports == 0
            || budget.charts == 0
            || budget.polynomial_terms == 0
            || budget.operations == 0
            || budget.rational_digits == 0
            || budget.degree == 0
            || budget.label_axes == 0
            || budget.index_magnitude == 0
            || budget.aggregate_expression_bytes == 0
        {
            return Err(limit("zero proof budget"));
        }
        // Cheap structural/text preflight precedes geometry reconstruction and CAS expansion.
        let raw = input.input();
        if family.loops() != raw.loops
            || family.shells().is_empty()
            || family.physical_slots() != raw.edges.len()
            || raw.edges.iter().any(|e| e.routing.len() != raw.loops)
            || family.inverse_routing().len() != raw.loops
            || family
                .inverse_routing()
                .iter()
                .any(|r| r.len() != raw.loops)
        {
            return Err(invalid("input/family dimensions"));
        }
        if raw.loops > budget.loops || raw.edges.len() > budget.physical_rows {
            return Err(limit("input loop/physical row budget"));
        }
        let mut text_bytes = 0usize;
        for text in raw
            .edges
            .iter()
            .flat_map(|e| e.routing.iter().chain(std::iter::once(&e.mass_squared)))
            .chain(raw.chemical_potentials.iter())
            .chain(raw.targets.iter().map(|t| &t.numerator))
        {
            text_bytes = text_bytes
                .checked_add(text.len())
                .ok_or_else(|| limit("input text overflow"))?;
            if text_bytes > budget.input_text_bytes {
                return Err(limit("input text budget"));
            }
        }
        let mut m = Meter {
            run,
            budget,
            operations: 0,
            terms: 0,
            bytes: 0,
            charts: 0,
        };
        for a in family
            .coordinates()
            .iter()
            .chain(family.factors())
            .chain(family.targets().iter().flat_map(|t| t.values()))
        {
            m.preflight(a)?;
        }
        let routing = raw
            .edges
            .iter()
            .map(|e| {
                e.routing
                    .iter()
                    .map(|a| {
                        if a.len() > budget.rational_digits {
                            return Err(limit("routing input text size budget"));
                        }
                        let a = Atom::parse(a, "rustflow_density", Default::default())
                            .map_err(|e| invalid(&e.to_string()))?;
                        m.rational(&rat(&a)?)?;
                        Ok(a)
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        for a in family.inverse_routing().iter().flatten() {
            m.rational(&rat(a)?)?;
        }
        for s in family.shells() {
            m.rational(&s.chemical_potential)?;
        }

        let k = family.shells().len();
        let h = family
            .loops()
            .checked_sub(k)
            .ok_or_else(|| invalid("cut count"))?;
        if !(1..=2).contains(&h) {
            return Err(unsupported("requires one or two virtual coordinates"));
        }
        if options.positive_compact_energy_powers {
            return Err(unsupported("no inverse completion admission"));
        }
        if input.physical_masses().iter().any(|x| !x.is_zero())
            || family
                .shells()
                .iter()
                .any(|s| !s.mass_squared.is_zero() || s.chemical_potential <= 0)
        {
            return Err(unsupported(
                "requires zero physical masses and positive chemical magnitudes",
            ));
        }
        let cuts = family
            .shells()
            .iter()
            .map(|s| s.physical_slot)
            .collect::<Vec<_>>();
        m.tick(1)?;
        let expected = input
            .occupied_cut(&cuts, budget.geometry_terms)?
            .at_physical_masses();
        m.tick(1)?;
        let signature = family_signature(family);
        if signature != family_signature(&expected) {
            return Err(invalid("exact input/family binding mismatch"));
        }
        if shifted.is_empty()
            || shifted.windows(2).any(|p| p[0] >= p[1])
            || shifted
                .iter()
                .any(|s| *s >= family.physical_slots() || cuts.contains(s))
        {
            return Err(invalid("sorted distinct uncut physical placement required"));
        }
        if family.physical_slots() > budget.physical_rows {
            return Err(limit("physical row budget"));
        }
        m.tick(
            raw.edges
                .len()
                .checked_mul(raw.loops)
                .and_then(|x| x.checked_mul(raw.loops))
                .ok_or_else(|| limit("routing operation envelope"))?,
        )?;
        let rows = matmul(&routing, family.inverse_routing());
        m.tick(1)?;
        for (s, row) in rows.iter().enumerate() {
            m.tick(1)?;
            for a in row {
                m.rational(&rat(a)?)?;
            }
            if !zero(square(row, family.coordinates(), family.loops()) - &family.factors()[s]) {
                return Err(unsupported(
                    "physical factor is not its exact rational massless momentum square",
                ));
            }
        }
        let cutset = cuts.iter().copied().collect::<BTreeSet<_>>();
        let shiftset = shifted.iter().copied().collect::<BTreeSet<_>>();
        let t = affine(&rows, &shiftset, &cutset, k, h, &mut m)?;
        let mut translated = rows.clone();
        for r in &mut translated {
            for c in 0..k {
                r[c] =
                    (&r[c] + r[k..].iter().zip(&t).map(|(a, b)| a * &b[c]).sum::<Atom>()).expand();
            }
        }
        let mut virtual_slots = Vec::new();
        let mut pure = Vec::new();
        for (s, r) in translated
            .iter()
            .enumerate()
            .filter(|(s, _)| !cutset.contains(s))
        {
            if r[k..].iter().any(|a| !a.is_zero()) {
                virtual_slots.push(s);
            } else {
                if !shiftset.contains(&s) {
                    return Err(unsupported("unshifted pure compact physical factor"));
                }
                let nonzero = r[..k]
                    .iter()
                    .enumerate()
                    .filter(|(_, a)| !a.is_zero())
                    .collect::<Vec<_>>();
                if nonzero.len() != 2 || !zero(nonzero[0].1 + nonzero[1].1) {
                    return Err(unsupported(
                        "pure compact factor lacks spacelike difference gap",
                    ));
                }
                pure.push(s);
            }
        }
        let count = 1usize
            .checked_shl(
                u32::try_from(virtual_slots.len()).map_err(|_| limit("support count overflow"))?,
            )
            .ok_or_else(|| limit("support count overflow"))?;
        if count > budget.supports {
            return Err(limit("complete support budget"));
        }
        let mut degrees = Vec::new();
        let scalar = family.loops() * (family.loops() + 1) / 2;
        for (slot, factor) in family.factors().iter().enumerate() {
            if slot >= family.input_slots() {
                degrees.push(0);
                continue;
            }
            m.preflight(factor)?;
            let terms = exact_coefficient_list(factor, family.coordinates())?;
            m.tick(1)?;
            m.terms(terms.len())?;
            let mut degree = 0;
            for (mon, c) in terms {
                if family
                    .coordinates()
                    .iter()
                    .any(|v| !derivative(&c, v).is_zero())
                {
                    return Err(unsupported("nonpolynomial factor coefficient"));
                }
                let ps = powers(&mon, family.coordinates())?;
                let mut d = 0;
                for (i, p) in ps.into_iter().enumerate() {
                    if p < 0 {
                        return Err(unsupported("rational factor insertion"));
                    }
                    d = add(d, mul(p as u64, if i < scalar { 2 } else { 1 })?)?;
                }
                degree = degree.max(d);
            }
            degrees.push(degree);
        }
        let mut schur = Rational::zero();
        for slot in virtual_slots.iter().filter(|s| shiftset.contains(s)) {
            for c in &translated[*slot][..k] {
                let c = rat(c)?;
                schur += &c * &c;
                m.tick(2)?;
            }
        }
        m.rational(&schur)?;
        let mut supports = Vec::new();
        for mask in 0..count {
            m.tick(1)?;
            let active = virtual_slots
                .iter()
                .enumerate()
                .filter(|(i, _)| mask & (1usize << i) != 0)
                .map(|(_, s)| *s)
                .collect::<Vec<_>>();
            let a = active
                .iter()
                .map(|s| translated[*s][k..].to_vec())
                .collect::<Vec<_>>();
            m.tick(
                a.len()
                    .checked_mul(h)
                    .ok_or_else(|| limit("rank envelope"))?,
            )?;
            let r = rank(&a, h);
            m.tick(1)?;
            if r < h {
                supports.push(Support {
                    active,
                    rank: r,
                    class: SupportClass::FreeVirtual {
                        null_direction: null_vector(&a, h)?
                            .iter()
                            .map(Atom::to_canonical_string)
                            .collect(),
                    },
                    full: None,
                    charts: Vec::new(),
                });
                continue;
            }
            let g = gaussian(&translated, &active, &shiftset, k, h, &mut m)?;
            let chart = charts(&translated, &active, k, h, &g, &mut m)?;
            let shifted_positive = active
                .iter()
                .filter(|s| shiftset.contains(s))
                .copied()
                .collect::<Vec<_>>();
            let class = if shifted_positive.is_empty() {
                if g.v.iter().flatten().any(|a| !a.is_zero()) {
                    return Err(invalid("translated vacuum retains compact invariant"));
                }
                SupportClass::RegulatedVacuumZero
            } else {
                SupportClass::FiniteEtaJets {
                    shifted_positive,
                    kinematics_vanish: g.v.iter().flatten().all(Atom::is_zero),
                }
            };
            supports.push(Support {
                active,
                rank: r,
                class,
                full: Some(g.native),
                charts: chart,
            });
        }
        m.tick(1)?;
        Ok(Self {
            version: VERSION,
            prescription: PRESCRIPTION,
            input_identity: input.identity().into(),
            family_signature: signature,
            source_options: options,
            cuts,
            shifted: shifted.to_vec(),
            virtual_slots,
            pure_compact_slots: pure,
            routed_rows: strings(&rows),
            translation: strings(&t),
            translated_rows: strings(&translated),
            compact_loops: k,
            virtual_loops: h,
            schur_bound: schur.to_string(),
            auxiliary_coefficients: shifted.iter().map(|s| (*s, "1".into())).collect(),
            complete_supports: supports,
            factor_degrees: degrees,
            operations: m.operations,
            polynomial_terms: m.terms,
        })
    }
    pub(crate) fn identity(&self) -> Result<String> {
        let b = serde_json::to_vec(self).map_err(|e| Error::Cache(e.to_string()))?;
        Ok(format!("{VERSION}:{}", blake3::hash(&b)))
    }
    pub(crate) fn validate_binding(
        &self,
        input: &PreparedDensityInput,
        family: &OccupiedCutFamily,
        shifted: &[usize],
        options: WeightedSourceOptions,
    ) -> Result<()> {
        if self.input_identity != input.identity()
            || self.family_signature != family_signature(family)
            || self.shifted != shifted
            || self.source_options != options
        {
            return Err(invalid("capability binding mismatch"));
        }
        Ok(())
    }
    pub(crate) fn validate_binding_with_context(
        &self,
        input: &PreparedDensityInput,
        family: &OccupiedCutFamily,
        shifted: &[usize],
        options: WeightedSourceOptions,
        budget: PartialOriginBudget,
        run: &RunContext,
    ) -> Result<()> {
        preflight_family(family, budget, run)?;
        self.validate_binding(input, family, shifted, options)?;
        run.cancellation.check()?;
        Ok(())
    }
    pub(crate) fn replay_with_context(
        &self,
        input: &PreparedDensityInput,
        family: &OccupiedCutFamily,
        budget: PartialOriginBudget,
        run: &RunContext,
    ) -> Result<()> {
        run.cancellation.check()?;
        let other = Self::new_with_context(
            input,
            family,
            &self.shifted,
            self.source_options,
            budget,
            run,
        )?;
        if other.identity()? != self.identity()? {
            return Err(invalid("complete native/support/chart replay mismatch"));
        }
        Ok(())
    }
    pub(crate) fn audit_label_with_context(
        &self,
        family: &OccupiedCutFamily,
        label: &[i64],
        budget: PartialOriginBudget,
        run: &RunContext,
    ) -> Result<FiniteLabelOriginAudit> {
        run.cancellation.check()?;
        preflight_family(family, budget, run)?;
        if label.len() > budget.label_axes
            || label
                .iter()
                .any(|n| n.unsigned_abs() > budget.index_magnitude)
        {
            return Err(limit("finite label preflight"));
        }
        if self.family_signature != family_signature(family)
            || label.len() < family.factors().len()
            || label[family.factors().len()..].iter().any(|n| *n != 0)
        {
            return Err(invalid("finite label family/storage mismatch"));
        }
        self.audit_bound_label(family, label, budget, run)
    }
    fn audit_bound_label(
        &self,
        family: &OccupiedCutFamily,
        label: &[i64],
        budget: PartialOriginBudget,
        run: &RunContext,
    ) -> Result<FiniteLabelOriginAudit> {
        run.cancellation.check()?;
        if label.len() > budget.label_axes
            || label
                .iter()
                .any(|n| n.unsigned_abs() > budget.index_magnitude)
        {
            return Err(limit("finite label bound budget"));
        }
        if label.len() < family.factors().len()
            || label[family.factors().len()..].iter().any(|n| *n != 0)
        {
            return Err(invalid("finite label role/arity preflight"));
        }
        // All roles are validated before any algebraic zero classification.
        for (i, role) in family.roles().iter().enumerate() {
            if *role == IndexRole::Occupation && label[i] < 0 {
                return Err(invalid("negative occupation"));
            }
            if i >= family.physical_slots() && i < family.input_slots() && label[i] > 0 {
                return Err(unsupported("positive completion power"));
            }
        }
        let cut_vanishing = family.shells().iter().any(|s| label[s.physical_slot] <= 0);
        let active = self
            .virtual_slots
            .iter()
            .filter(|s| label[**s] > 0)
            .copied()
            .collect::<Vec<_>>();
        let support = self
            .complete_supports
            .iter()
            .find(|s| s.active == active)
            .ok_or_else(|| invalid("missing complete support"))?;
        let mut p = 0;
        for s in &active {
            p = add(p, label[*s] as u64)?;
        }
        let mut r = 0;
        for (i, &n) in label[..family.input_slots()].iter().enumerate() {
            if n < 0 {
                r = add(r, mul(n.unsigned_abs(), self.factor_degrees[i])?)?;
            }
        }
        let mut j = 0;
        let mut radial = 0;
        for s in family.shells() {
            run.cancellation.check()?;
            let n = label[s.physical_slot].max(1) as u64;
            let upper = label[s.upper_slot] as u64;
            let lower = label[s.lower_slot] as u64;
            j = add(
                j,
                add(
                    n - 1,
                    add(upper.saturating_sub(1), lower.saturating_sub(1))?,
                )?,
            )?;
            radial = radial.max(add(
                add(mul(2, n)?, lower)?,
                mul(2, upper.saturating_sub(1))?,
            )?);
        }
        let bound = mul(2, add(add(add(p, r)?, j)?, 2)?)?
            .max(radial)
            .max(self.compact_loops as u64 + 2);
        let b = i64::try_from(bound).map_err(|_| limit("dimension witness overflow"))?;
        let numerator = b
            .checked_mul(5)
            .and_then(|b| b.checked_add(1))
            .ok_or_else(|| limit("dimension witness overflow"))?;
        let d = Rational::from((numerator, 5));
        let n = bound / 2 + r + 1;
        let left = &d
            - &Rational::from(
                i64::try_from(add(add(p, j)?, n)?).map_err(|_| limit("jet witness overflow"))?,
            );
        // For h=1 there is no D-growing UV subtraction and lambda=D/2-P.
        let (n, left) = if self.virtual_loops == 1 {
            (
                0,
                &d / &Rational::from(2)
                    - Rational::from(i64::try_from(add(p, j)?).map_err(|_| limit("jet overflow"))?),
            )
        } else {
            (n, left)
        };
        if left <= 0 {
            return Err(invalid("finite-jet bound did not give strict margin"));
        }
        let classification = if cut_vanishing {
            "required-cut-zero".into()
        } else if family.shells().iter().any(|s| label[s.lower_slot] > 0) {
            "lower-origin-contact-zero".into()
        } else {
            match &support.class {
                SupportClass::FreeVirtual { .. } => "free-virtual-polynomial-zero".into(),
                SupportClass::RegulatedVacuumZero => "joint-regulator-vacuum-zero".into(),
                SupportClass::FiniteEtaJets { .. } => "finite-positive-eta-jets".into(),
            }
        };
        Ok(FiniteLabelOriginAudit {
            label: label.to_vec(),
            active_virtual_slots: active,
            classification,
            positive_virtual_degree: p,
            polynomial_degree: r,
            jets_including_lower: j,
            cut_vanishing,
            dimension_witness: d.to_string(),
            subtracted_terms: n,
            remaining_mass_power: left.to_string(),
        })
    }
    /// Every box is authorized by the universal finite-label theorem, not a finite target list.
    pub(crate) fn zero_domains_with_context<const N: usize>(
        &self,
        family: &OccupiedCutFamily,
        shifted: &[usize],
        options: WeightedSourceOptions,
        budget: PartialOriginBudget,
        run: &RunContext,
    ) -> Result<Vec<IndexDomain<N>>> {
        run.cancellation.check()?;
        preflight_family(family, budget, run)?;
        if N > budget.label_axes
            || self
                .complete_supports
                .len()
                .checked_add(family.shells().len())
                .is_none_or(|n| n > budget.supports)
        {
            return Err(limit("zero domain output budget"));
        }
        if self.family_signature != family_signature(family)
            || self.source_options != options
            || self.shifted != shifted
            || N < family.factors().len()
        {
            return Err(invalid("zero guard binding/arity"));
        }
        let mut base = [IndexBounds::unbounded(); N];
        for (i, b) in base.iter_mut().enumerate() {
            if i >= family.factors().len() {
                *b = IndexBounds::fixed(0);
            } else if i >= family.physical_slots() && i < family.input_slots() {
                *b = bound(None, Some(0))?;
            } else if family.roles()[i] == IndexRole::RequiredCut {
                *b = bound(Some(1), None)?;
            } else if family.roles()[i] == IndexRole::Occupation {
                *b = bound(Some(0), None)?;
            }
        }
        let mut out = Vec::new();
        for s in family.shells() {
            run.cancellation.check()?;
            let mut b = base;
            b[s.lower_slot] = bound(Some(1), None)?;
            out.push(IndexDomain::new(b).map_err(|e| invalid(&e.to_string()))?);
        }
        if options.free_virtual_zero_sectors {
            for s in &self.complete_supports {
                run.cancellation.check()?;
                if matches!(s.class, SupportClass::FiniteEtaJets { .. }) {
                    continue;
                }
                let mut b = base;
                for slot in &self.virtual_slots {
                    b[*slot] = if s.active.contains(slot) {
                        bound(Some(1), None)?
                    } else {
                        bound(None, Some(0))?
                    };
                }
                out.push(IndexDomain::new(b).map_err(|e| invalid(&e.to_string()))?);
            }
        }
        Ok(out)
    }
}

#[cfg(test)]
impl PartialOriginCapability {
    fn new(
        i: &PreparedDensityInput,
        f: &OccupiedCutFamily,
        s: &[usize],
        o: WeightedSourceOptions,
        b: PartialOriginBudget,
    ) -> Result<Self> {
        Self::new_with_context(i, f, s, o, b, &RunContext::default())
    }
    fn replay(
        &self,
        i: &PreparedDensityInput,
        f: &OccupiedCutFamily,
        b: PartialOriginBudget,
    ) -> Result<()> {
        self.replay_with_context(i, f, b, &RunContext::default())
    }
    fn audit_label(&self, f: &OccupiedCutFamily, l: &[i64]) -> Result<FiniteLabelOriginAudit> {
        self.audit_label_with_context(f, l, Default::default(), &RunContext::default())
    }
    fn zero_domains<const N: usize>(
        &self,
        f: &OccupiedCutFamily,
        s: &[usize],
        o: WeightedSourceOptions,
    ) -> Result<Vec<IndexDomain<N>>> {
        self.zero_domains_with_context(f, s, o, Default::default(), &RunContext::default())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn input() -> PreparedDensityInput {
        serde_json::from_str::<super::super::DensityInput>(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/examples/finite_density/massless_three_loop_chain.json"
        )))
        .unwrap()
        .prepare()
        .unwrap()
    }
    fn options() -> WeightedSourceOptions {
        WeightedSourceOptions {
            free_virtual_zero_sectors: true,
            ..Default::default()
        }
    }
    fn baseline(f: &OccupiedCutFamily) -> Vec<i64> {
        let mut a = vec![0; 16];
        for i in 0..f.physical_slots() {
            a[i] = 1;
        }
        a
    }
    fn proof(
        cuts: &[usize],
        shifted: &[usize],
    ) -> (
        PreparedDensityInput,
        OccupiedCutFamily,
        PartialOriginCapability,
    ) {
        let i = input();
        let f = i.occupied_cut(cuts, 1024).unwrap().at_physical_masses();
        let p =
            PartialOriginCapability::new(&i, &f, shifted, options(), Default::default()).unwrap();
        (i, f, p)
    }
    #[test]
    fn actual_singleton_and_both_double_masks_replay_all_native_supports() {
        for (cuts, shifted, expected) in [
            (vec![3], vec![0], 16),
            (vec![0, 3], vec![1, 2], 4),
            (vec![0, 3], vec![2, 4], 4),
        ] {
            let (i, f, p) = proof(&cuts, &shifted);
            assert_eq!(p.complete_supports.len(), expected);
            p.replay(&i, &f, Default::default()).unwrap();
            assert!(
                p.complete_supports
                    .iter()
                    .filter(|s| s.rank == p.virtual_loops)
                    .all(|s| !s.charts.is_empty() && s.full.is_some())
            );
            let a = baseline(&f);
            let audit = p.audit_label(&f, &a).unwrap();
            assert_eq!(audit.classification, "finite-positive-eta-jets");
            let zeros = p.zero_domains::<16>(&f, &shifted, options()).unwrap();
            assert!(
                !zeros
                    .iter()
                    .any(|d| d.contains(&a.clone().try_into().unwrap()))
            );
            for shell in f.shells() {
                let mut contact = a.clone();
                contact[shell.lower_slot] = 3;
                assert!(
                    zeros
                        .iter()
                        .any(|d| d.contains(&contact.clone().try_into().unwrap()))
                );
                assert_eq!(
                    p.audit_label(&f, &contact).unwrap().classification,
                    "lower-origin-contact-zero"
                );
            }
            if let Ok(dir) = std::env::var("PARTIAL_ORIGIN_REPORT") {
                std::fs::create_dir_all(&dir).unwrap();
                std::fs::write(
                    std::path::Path::new(&dir)
                        .join(format!("cuts-{:?}-shift-{:?}.json", cuts, shifted)),
                    serde_json::to_vec_pretty(&p).unwrap(),
                )
                .unwrap();
            }
        }
    }
    #[test]
    fn shifted_only_rank_counterexample_and_v_zero_are_not_finite_eta_zeros() {
        let (i, f, _) = proof(&[3], &[0]);
        let p =
            PartialOriginCapability::new(&i, &f, &[0, 1], options(), Default::default()).unwrap();
        let active = vec![1, 2, 4];
        let s = p
            .complete_supports
            .iter()
            .find(|s| s.active == active)
            .unwrap();
        assert_eq!(s.rank, 2);
        assert!(matches!(
            &s.class,
            SupportClass::FiniteEtaJets {
                kinematics_vanish: true,
                ..
            }
        ));
        let mut a = baseline(&f);
        a[0] = 0;
        assert_eq!(
            p.audit_label(&f, &a).unwrap().classification,
            "finite-positive-eta-jets"
        );
        assert!(
            !p.zero_domains::<16>(&f, &[0, 1], options())
                .unwrap()
                .iter()
                .any(|d| d.contains(&a.clone().try_into().unwrap()))
        );
    }
    #[test]
    fn all_finite_jets_and_negative_polynomial_powers_have_explicit_open_witness() {
        let (_, f, p) = proof(&[3], &[0]);
        for n in [1, 2, 9] {
            for upper in [0, 1, 7] {
                for lower in [0, 1, 11] {
                    let mut a = baseline(&f);
                    a[3] = n;
                    a[f.shells()[0].upper_slot] = upper;
                    a[f.shells()[0].lower_slot] = lower;
                    a[5] = -17;
                    let audit = p.audit_label(&f, &a).unwrap();
                    assert_eq!(
                        audit.jets_including_lower,
                        ((n - 1) + (upper - 1).max(0) + (lower - 1).max(0)) as u64
                    );
                    assert!(audit.subtracted_terms > 0);
                    assert!(audit.polynomial_degree >= 17);
                }
            }
        }
        let mut huge = baseline(&f);
        huge[0] = i64::MAX;
        assert!(matches!(p.audit_label(&f, &huge), Err(Error::Limit(_))));
    }
    #[test]
    fn binding_options_roles_and_storage_fail_before_zero() {
        let (i, f, p) = proof(&[3], &[0]);
        assert!(p.validate_binding(&i, &f, &[0, 1], options()).is_err());
        let mut opt = options();
        opt.free_virtual_zero_sectors = false;
        assert!(p.validate_binding(&i, &f, &[0], opt).is_err());
        assert!(p.zero_domains::<16>(&f, &[0, 1], options()).is_err());
        let other = i.occupied_cut(&[0], 1024).unwrap().at_physical_masses();
        assert!(p.audit_label(&other, &baseline(&other)).is_err());
        let mut a = baseline(&f);
        a[3] = 0;
        a[f.shells()[0].lower_slot] = -1;
        assert!(p.audit_label(&f, &a).is_err());
        a[f.shells()[0].lower_slot] = 0;
        a[5] = 1;
        assert!(p.audit_label(&f, &a).is_err());
        a[5] = 0;
        a[15] = 1;
        assert!(p.audit_label(&f, &a).is_err());
        let mut inv = options();
        inv.positive_compact_energy_powers = true;
        assert!(PartialOriginCapability::new(&i, &f, &[0], inv, Default::default()).is_err());
    }
    #[test]
    fn omitted_support_chart_and_native_off_null_data_cannot_replay() {
        let (i, f, p) = proof(&[3], &[0]);
        let mut q = p.clone();
        q.complete_supports.pop();
        assert!(q.replay(&i, &f, Default::default()).is_err());
        let mut q = p.clone();
        q.complete_supports
            .iter_mut()
            .find(|s| !s.charts.is_empty())
            .unwrap()
            .charts
            .pop();
        assert!(q.replay(&i, &f, Default::default()).is_err());
        let mut q = p.clone();
        q.complete_supports
            .iter_mut()
            .find_map(|s| s.full.as_mut())
            .unwrap()
            .native_f = "0".into();
        assert!(q.replay(&i, &f, Default::default()).is_err());
    }
    #[test]
    fn every_global_budget_rejects_without_partial_capability() {
        let i = input();
        let f = i.occupied_cut(&[3], 1024).unwrap().at_physical_masses();
        let b = PartialOriginBudget::default();
        for bad in [
            PartialOriginBudget { loops: 2, ..b },
            PartialOriginBudget {
                input_text_bytes: 1,
                ..b
            },
            PartialOriginBudget {
                geometry_terms: 0,
                ..b
            },
            PartialOriginBudget {
                physical_rows: 4,
                ..b
            },
            PartialOriginBudget { supports: 15, ..b },
            PartialOriginBudget { charts: 1, ..b },
            PartialOriginBudget {
                polynomial_terms: 1,
                ..b
            },
            PartialOriginBudget { operations: 1, ..b },
            PartialOriginBudget {
                rational_digits: 0,
                ..b
            },
            PartialOriginBudget { degree: 0, ..b },
        ] {
            assert!(matches!(
                PartialOriginCapability::new(&i, &f, &[0], options(), bad),
                Err(Error::Limit(_))
            ));
        }
    }
    #[test]
    fn invalid_partial_placement_and_physical_mass_reject() {
        let i = input();
        let f = i.occupied_cut(&[0, 3], 1024).unwrap().at_physical_masses();
        assert!(PartialOriginCapability::new(&i, &f, &[1], options(), Default::default()).is_err());
        assert!(PartialOriginCapability::new(&i, &f, &[2], options(), Default::default()).is_err());
        assert!(PartialOriginCapability::new(&i, &f, &[0], options(), Default::default()).is_err());
        let mut raw = i.input().clone();
        raw.edges[1].mass_squared = "1".into();
        let i = raw.prepare().unwrap();
        let f = i.occupied_cut(&[3], 1024).unwrap().at_physical_masses();
        assert!(PartialOriginCapability::new(&i, &f, &[0], options(), Default::default()).is_err());
    }
    #[test]
    fn actual_native_full_channel_rejects_timelike_sum() {
        let run = RunContext::default();
        let mut m = Meter {
            run: &run,
            bytes: 0,
            budget: Default::default(),
            operations: 0,
            terms: 0,
            charts: 0,
        };
        let rows = vec![
            vec![Atom::one(), Atom::zero(), Atom::one()],
            vec![Atom::zero(), Atom::num(-1), Atom::one()],
        ];
        assert!(
            matches!(gaussian(&rows,&[0,1],&BTreeSet::from([0,1]),2,1,&mut m),Err(Error::Unsupported(ref s)) if s.contains("nonpositive"))
        );
    }
    #[test]
    fn full_off_null_regulator_gap_and_same_phase_hold_on_closed_support_faces() {
        fn eval(p: &[Monomial], a: &[Rational]) -> Rational {
            p.iter().fold(Rational::zero(), |sum, t| {
                let mut c =
                    rat(
                        &Atom::parse(&t.coefficient, "partial_origin_test", Default::default())
                            .unwrap(),
                    )
                    .unwrap();
                for (x, e) in a.iter().zip(&t.exponents) {
                    for _ in 0..*e {
                        c = &c * x;
                    }
                }
                sum + c
            })
        }
        let energies = [(5i64, 3i64, 4i64), (13, 5, 12)];
        let temperature = 2i64;
        for (omega, gamma, r) in energies {
            assert_eq!(omega * omega - gamma * gamma, r * r);
            assert!(omega >= r);
            assert!(gamma < 3 * temperature);
        }
        let (wi, gi, ri) = energies[0];
        let (wj, gj, rj) = energies[1];
        let reverse = wi * wj - gi * gj;
        assert_eq!(
            reverse * reverse - ri * ri * rj * rj,
            (wi * gj - gi * wj).pow(2)
        );
        assert!(reverse >= ri * rj);
        let mut checked = 0;
        let mut zero_s_faces = 0;
        for (cuts, shifted) in [
            (vec![3], vec![0]),
            (vec![0, 3], vec![1, 2]),
            (vec![0, 3], vec![2, 4]),
        ] {
            let (_, _, p) = proof(&cuts, &shifted);
            for support in &p.complete_supports {
                let Some(full) = &support.full else {
                    continue;
                };
                for mask in 1usize..(1usize << support.active.len()) {
                    let alpha = (0..support.active.len())
                        .map(|i| Rational::from(if mask & (1 << i) != 0 { 1 } else { 0 }))
                        .collect::<Vec<_>>();
                    let u = eval(&full.u, &alpha);
                    if u.is_zero() {
                        continue;
                    }
                    let eta = Rational::from((1, 3));
                    let eta_coefficient = eval(&full.eta, &alpha);
                    let s = support
                        .active
                        .iter()
                        .zip(&alpha)
                        .filter(|(slot, _)| p.shifted.contains(slot))
                        .map(|(_, a)| a)
                        .sum::<Rational>();
                    assert_eq!(eta_coefficient, &u * &s);
                    for cosine in [-1, 0, 1] {
                        let mut real = &eta * &eta_coefficient;
                        let mut imaginary = Rational::zero();
                        for (pair, poly) in &full.pairs {
                            let c = eval(poly, &alpha);
                            assert!(c >= 0);
                            let (wi, gi, ri) = energies[pair[0]];
                            let (wj, gj, rj) = energies[pair[1]];
                            let real_dot = 2 * (wi * wj - gi * gj - ri * rj * cosine);
                            assert!(real_dot >= 2 * ri * rj * (1 - cosine));
                            real += &c * &Rational::from(real_dot);
                            imaginary += &c * &Rational::from(-2 * (wi * gj + wj * gi));
                        }
                        for (i, poly) in full.diagonal.iter().enumerate() {
                            let (w, g, _) = energies[i];
                            imaginary += eval(poly, &alpha) * Rational::from(-2 * w * g);
                        }
                        let mut expected_phase = Rational::zero();
                        for (j, poly) in full.line_masses.iter().enumerate() {
                            let eps = Rational::from(j as i64 + 1);
                            let c = eval(poly, &alpha);
                            assert_eq!(c, &u * &alpha[j]);
                            imaginary -= &c * &eps;
                            expected_phase -= &u * &alpha[j] * eps;
                        }
                        assert!(real >= &eta * &eta_coefficient);
                        if s.is_zero() {
                            zero_s_faces += 1;
                            assert!(full.pairs.iter().all(|(_, p)| eval(p, &alpha).is_zero()));
                            assert!(full.diagonal.iter().all(|p| eval(p, &alpha).is_zero()));
                            assert!(real.is_zero());
                            assert_eq!(imaginary, expected_phase);
                            assert!(imaginary < 0);
                        } else {
                            assert!(real > 0);
                        }
                        checked += 1;
                    }
                }
            }
        }
        assert!(checked > 100);
        assert!(zero_s_faces > 0);
        if let Ok(dir) = std::env::var("PARTIAL_ORIGIN_REPORT") {
            std::fs::write(std::path::Path::new(&dir).join("off-null-regulator-checks.json"),serde_json::to_vec_pretty(&serde_json::json!({"exact_face_angle_checks":checked,"S_zero_same_phase_checks":zero_s_faces,"pole_data_omega_gamma_radius":energies,"fermi_T":temperature,"strip_check":"gamma < 3*T < pi*T; same-sign line-mass regulators only","no_arbitrary_complex_null_scaling":true})).unwrap()).unwrap();
        }
    }
    #[test]
    fn fresh_actual_sources_bind_origin_zero_guards_and_reject_foreign_context() {
        use super::super::guarded::GuardedMeasureIdentity;
        use rustred::solver::guarded::{
            GuardedApplicationStatus, GuardedProgram, GuardedSource, GuardedSourceSystem,
        };
        use std::sync::Arc;
        let (i, f, p) = proof(&[3], &[0]);
        let epsilon = symbol!("partial_origin_source_epsilon");
        let eta = symbol!("partial_origin_source_eta");
        let id = GuardedMeasureIdentity {
            measure: format!("fresh exact family; origin={}", p.identity().unwrap()),
            support: "polynomial completions; massless future cuts; positive eta".into(),
            orientation: "existing occupied routing".into(),
            normalization: "existing C/H distribution normalization".into(),
            branch: PRESCRIPTION.into(),
            deformation: "only physical slot 0 minus eta; cuts fixed".into(),
        };
        // Existing formal factory supplies all original rows/guards. No endpoint token or old rule is imported.
        let formal = f
            .guarded_sources::<16>(epsilon, 4, eta, &[0], 4096, vec![], id)
            .unwrap();
        let original = formal.context.sources();
        assert!(!original.native_sources().rows().is_empty());
        assert!(original.zero_domains().is_empty());
        let fresh_rows = original
            .sources()
            .iter()
            .zip(original.native_sources().rows())
            .map(|(info, row)| {
                GuardedSource::new(info.id.clone(), row.clone(), info.domain.clone())
                    .with_nonzero_conditions(info.nonzero_conditions.clone())
            })
            .collect();
        p.validate_binding(&i, &f, &[0], options()).unwrap();
        let zero = p.zero_domains::<16>(&f, &[0], options()).unwrap();
        let sources = Arc::new(
            GuardedSourceSystem::new(
                format!(
                    "{};sealed-partial-origin={}",
                    original.measure_id(),
                    p.identity().unwrap()
                ),
                *original.roles(),
                *original.native_sources().index_variables(),
                fresh_rows,
            )
            .unwrap()
            .with_zero_domains(zero)
            .unwrap(),
        );
        assert_eq!(
            sources.native_sources().rows(),
            original.native_sources().rows()
        );
        assert_eq!(
            sources.native_sources().index_variables(),
            original.native_sources().index_variables()
        );
        for (a, b) in sources.sources().iter().zip(original.sources()) {
            assert_eq!(a.domain, b.domain);
            assert_eq!(a.nonzero_conditions, b.nonzero_conditions);
        }
        let program = GuardedProgram::new(sources.clone(), vec![], []).unwrap();
        let raw = program.encode_native(Default::default()).unwrap();
        let loaded =
            GuardedProgram::decode_generated(&raw, sources.clone(), Default::default()).unwrap();
        assert!(
            GuardedProgram::decode_generated(&raw, original.clone(), Default::default()).is_err()
        );
        let a: [i64; 16] = baseline(&f).try_into().unwrap();
        assert!(!matches!(
            loaded.apply(&a).unwrap().status,
            GuardedApplicationStatus::Zero
        ));
        let mut lower = a;
        lower[f.shells()[0].lower_slot] = 5;
        assert!(matches!(
            loaded.apply(&lower).unwrap().status,
            GuardedApplicationStatus::Zero
        ));
        let mut invalid = lower;
        invalid[5] = 1;
        assert!(!matches!(
            loaded.apply(&invalid).unwrap().status,
            GuardedApplicationStatus::Zero
        ));
        invalid = lower;
        invalid[15] = 1;
        assert!(!matches!(
            loaded.apply(&invalid).unwrap().status,
            GuardedApplicationStatus::Zero
        ));
        let mut vacuum = a;
        vacuum[0] = 0;
        assert!(matches!(
            loaded.apply(&vacuum).unwrap().status,
            GuardedApplicationStatus::Zero
        ));
        if let Ok(dir) = std::env::var("PARTIAL_ORIGIN_REPORT") {
            std::fs::write(std::path::Path::new(&dir).join("fresh-source-guard-binding.json"),serde_json::to_vec_pretty(&serde_json::json!({"fresh_rows":sources.sources().len(),"zero_domains":sources.zero_domains().len(),"measure_id":sources.measure_id(),"program_blake3":blake3::hash(&raw).to_string(),"program_bytes":raw.len(),"imported_rules":0,"all_original_rows_guards_conditions_unchanged":true})).unwrap()).unwrap();
        }
    }
}
