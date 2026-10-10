//! Exact endpoint evidence for two null occupied vectors and one spacelike
//! external momentum direction in every virtual denominator.
//!
//! The general multi-loop query supplies degree evidence only. A separate
//! sealed permit covers the proved singleton meromorphic germ and independent
//! rank-one virtual blocks, including the pure compact transfer, and binds the origin prescription to the physical family.
//! Native source closure, boundary integration, transport and endpoint
//! projection remain required. No virtual integral or period is evaluated here.
use super::{PreparedDensityInput, massless_contour::massless_channel_evidence};
use crate::algebra::{matmul, rref};
use crate::coefficient::{exact_coefficient_list, powers};
use crate::family::substitute;
use crate::{Error, Result};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use symbolica::prelude::*;

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct EndpointDegree {
    /// Twice the coefficient multiplying D.
    pub twice_dimension_coefficient: i32,
    pub constant: i32,
}
impl EndpointDegree {
    pub fn at(&self, dimension: &Rational) -> Rational {
        dimension * &Rational::from((self.twice_dimension_coefficient, 2))
            + &Rational::from(self.constant)
    }
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct EndpointInequality {
    pub stratum: String,
    /// The strict condition is degree(D)>0. Logarithms do not change a strict
    /// positive margin; equality is deliberately not certified.
    pub degree: EndpointDegree,
}

#[derive(Clone, Debug, Serialize)]
pub struct AffineSoftStratum {
    pub physical_slots: Vec<usize>,
    pub normal_rank: usize,
    pub positive_power_sum: i32,
    pub eta_degree: EndpointDegree,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq, PartialOrd, Ord)]
pub struct CompactMonomialEnvelope {
    /// E_1^r E_2^s h^(h_virtual*D/2 + h_power_offset).
    pub energy_powers: [i32; 2],
    pub h_power_offset: i32,
    /// H_0 is the upper occupation; H_s, s>0, is its normalized derivative.
    pub upper_indices: [u16; 2],
}

#[derive(Clone, Debug, Serialize)]
pub struct TargetEndpointEvidence {
    pub target_index: usize,
    pub original_cut_powers: [u16; 2],
    pub affine_soft_strata: Vec<AffineSoftStratum>,
    pub compact_envelope: Vec<CompactMonomialEnvelope>,
    pub inequalities: Vec<EndpointInequality>,
}
impl TargetEndpointEvidence {
    pub fn strict_margins(&self, dimension: &Rational) -> Vec<Rational> {
        self.inequalities
            .iter()
            .map(|x| x.degree.at(dimension))
            .collect()
    }
    pub fn inequalities_hold(&self, dimension: &Rational) -> bool {
        self.strict_margins(dimension)
            .iter()
            .all(|v| v > &Rational::zero())
    }
}

#[derive(Clone, Copy, Debug)]
pub struct EndpointEvidenceBudget {
    pub cut_sets: usize,
    pub vertex_partitions: usize,
    pub affine_subsets: usize,
    pub numerator_terms: usize,
    pub wick_branches: usize,
    pub mass_jet_states: usize,
}
impl Default for EndpointEvidenceBudget {
    fn default() -> Self {
        Self {
            cut_sets: 1024,
            vertex_partitions: 65536,
            affine_subsets: 65536,
            numerator_terms: 65536,
            wick_branches: 100000,
            mass_jet_states: 100000,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct SingleSpacelikeEndpointEvidence {
    pub input_identity: String,
    pub cut_slots: [usize; 2],
    pub shifted_physical_slots: Vec<usize>,
    pub virtual_loops: usize,
    pub virtual_uv_continuation_required: bool,
    /// k_old=k_new+shift*q_1, applied before Gaussian means k_new=M*p+ell.
    pub virtual_affine_shift: Vec<String>,
    /// Each uncut momentum is R*k_new+b*(q_1-q_2).
    pub denominator_rows: Vec<(usize, Vec<String>, String)>,
    pub targets: Vec<TargetEndpointEvidence>,
    pub prescription: &'static str,
    pub endpoint_admitted: bool,
    pub numerical_evaluation_admitted: bool,
}

fn parse(text: &str) -> Result<Atom> {
    Atom::parse(text, "rustflow_density", Default::default())
        .map_err(|e| Error::InvalidInput(e.to_string()))
}
fn scalar(name: &str) -> Atom {
    Atom::var(symbol!(format!("rustflow_massless_endpoint::{name}")))
}
fn rank(rows: Vec<Vec<Atom>>) -> usize {
    rref(rows).1.len()
}
fn rational(a: &Atom) -> Result<Rational> {
    Rational::try_from(a.as_view())
        .map_err(|_| Error::Unsupported("endpoint evidence requires exact rational routing".into()))
}
fn consume(left: &mut usize, context: &str) -> Result<()> {
    *left = left
        .checked_sub(1)
        .ok_or_else(|| Error::Limit(format!("massless endpoint {context} budget exhausted")))?;
    Ok(())
}

/// Solve A*x=b with exact native row reduction, free variables set to zero.
fn solve(a: &[Vec<Atom>], b: &[Atom], columns: usize) -> Result<Vec<Atom>> {
    let rows = a
        .iter()
        .zip(b)
        .map(|(r, b)| {
            let mut r = r.clone();
            r.push(b.clone());
            r
        })
        .collect::<Vec<_>>();
    let (rows, pivots) = rref(rows);
    if pivots.contains(&columns) {
        return Err(Error::Unsupported(
            "uncut denominator routing cannot be reduced to the single spacelike direction q1-q2"
                .into(),
        ));
    }
    let mut out = vec![Atom::zero(); columns];
    for (row, &pivot) in rows.iter().zip(&pivots) {
        out[pivot] = row[columns].clone();
    }
    Ok(out)
}

#[derive(Clone)]
struct Row {
    slot: usize,
    virtual_part: Vec<Atom>,
    external: Atom,
}

/// Find every nonempty closed consistent affine intersection at fixed p!=0.
/// A line belongs to the closure exactly when its augmented row belongs to the
/// span of the selected equations. Inconsistent intersections are absent.
fn affine_strata(rows: &[Row], powers: &[i16], budget: usize) -> Result<Vec<AffineSoftStratum>> {
    let virtual_rows = rows
        .iter()
        .filter(|r| r.virtual_part.iter().any(|x| !x.is_zero()))
        .collect::<Vec<_>>();
    let count = 1_usize
        .checked_shl(
            u32::try_from(virtual_rows.len())
                .map_err(|_| Error::Limit("affine subset exponent".into()))?,
        )
        .ok_or_else(|| Error::Limit("affine subset enumeration overflow".into()))?;
    if count - 1 > budget {
        return Err(Error::Limit(format!(
            "affine soft intersections need {} subsets; budget {budget}",
            count - 1
        )));
    }
    let mut seen = BTreeSet::new();
    let mut out = Vec::new();
    for mask in 1..count {
        let selected = virtual_rows
            .iter()
            .enumerate()
            .filter(|(j, _)| mask & (1_usize << j) != 0)
            .map(|(_, r)| *r)
            .collect::<Vec<_>>();
        let normal_rank = rank(selected.iter().map(|r| r.virtual_part.clone()).collect());
        let augmented = selected
            .iter()
            .map(|r| {
                let mut a = r.virtual_part.clone();
                a.push(r.external.clone());
                a
            })
            .collect::<Vec<_>>();
        if rank(augmented.clone()) != normal_rank {
            continue;
        }
        let closure = virtual_rows
            .iter()
            .filter(|r| {
                let mut a = augmented.clone();
                let mut row = r.virtual_part.clone();
                row.push(r.external.clone());
                a.push(row);
                rank(a) == normal_rank
            })
            .map(|r| r.slot)
            .collect::<Vec<_>>();
        if !seen.insert(closure.clone()) {
            continue;
        }
        let positive_power_sum = closure.iter().map(|&s| i32::from(powers[s])).sum();
        out.push(AffineSoftStratum {
            physical_slots: closure,
            normal_rank,
            positive_power_sum,
            eta_degree: EndpointDegree {
                twice_dimension_coefficient: i32::try_from(normal_rank)
                    .map_err(|_| Error::Limit("normal rank".into()))?,
                constant: -positive_power_sum,
            },
        });
    }
    Ok(out)
}

/// Formal centered Gaussian contraction. External vector labels 0,1,2 are
/// q1,q2,u. Others are centered virtual vectors ell_i. Covariance parameters
/// are deliberately independent; retaining impossible terms only weakens the
/// envelope and cannot certify a false improvement from a period cancellation.
fn wick(
    edges: &[(usize, usize)],
    dimension: &Atom,
    h: &Atom,
    external_gram: &[Vec<Atom>],
    remaining: &mut usize,
) -> Result<Atom> {
    consume(remaining, "Wick branch")?;
    let first = edges.iter().enumerate().find_map(|(e, &(a, b))| {
        if a >= 3 {
            Some((e, 0, a))
        } else if b >= 3 {
            Some((e, 1, b))
        } else {
            None
        }
    });
    let Some((edge, side, vector)) = first else {
        return Ok(edges
            .iter()
            .fold(Atom::num(1), |a, &(i, j)| a * &external_gram[i][j]));
    };
    let partner = if side == 0 {
        edges[edge].1
    } else {
        edges[edge].0
    };
    let mut total = Atom::zero();
    for (e, &(a, b)) in edges.iter().enumerate() {
        for (s, other) in [(0, a), (1, b)] {
            if other < 3 || (e == edge && s == side) {
                continue;
            }
            let lo = vector.min(other) - 3;
            let hi = vector.max(other) - 3;
            let covariance = scalar(&format!("cov_{lo}_{hi}")) * h;
            let mut next = edges
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != edge && *j != e)
                .map(|(_, pair)| *pair)
                .collect::<Vec<_>>();
            if e == edge {
                total +=
                    covariance * dimension * wick(&next, dimension, h, external_gram, remaining)?;
            } else {
                let other_partner = if s == 0 { b } else { a };
                next.push((partner, other_partner));
                total += covariance * wick(&next, dimension, h, external_gram, remaining)?;
            }
        }
    }
    Ok(total)
}

fn gaussian_envelope(
    input: &PreparedDensityInput,
    target_index: usize,
    inverse_routing: &[Vec<Atom>],
    shift: &[Atom],
    term_budget: usize,
    wick_budget: &mut usize,
) -> Result<Atom> {
    original_polynomial_rank(input, target_index)?;
    let loops = input.input().loops;
    let target = &input.input().targets[target_index];
    // Reconstruct the validated, convention-adjusted original numerator. The
    // independent masses introduced during basis conversion cancel exactly.
    let mut numerator = input.basis().reconstruct(&input.targets()[target_index])?;
    for (factor, &power) in input.basis().slots().iter().zip(&target.powers) {
        numerator *= factor.clone().pow(i64::from(power));
    }
    numerator = numerator.together().cancel().expand();
    for &mass in input.independent_masses() {
        if !numerator.derivative(mass).is_zero() {
            return Err(Error::InvalidInput(
                "reconstructed original numerator depends on a physical mass".into(),
            ));
        }
    }
    let vectors = loops + 1; // q1,q2,u, and loops-2 centered virtual vectors.
    let mut formal_gram = vec![vec![Atom::zero(); vectors]; vectors];
    let mut gram_variables = Vec::new();
    let mut gram_pairs = Vec::new();
    for i in 0..vectors {
        for j in i..vectors {
            let v = scalar(&format!("gram_{i}_{j}"));
            formal_gram[i][j] = v.clone();
            formal_gram[j][i] = v.clone();
            gram_variables.push(v);
            gram_pairs.push((i, j));
        }
    }
    let mut original_vectors = vec![vec![Atom::zero(); vectors]; loops];
    for i in 0..loops {
        original_vectors[i][0] = inverse_routing[i][0].clone();
        original_vectors[i][1] = inverse_routing[i][1].clone();
        for v in 0..loops - 2 {
            let c = &inverse_routing[i][v + 2];
            let mean = scalar(&format!("mean_{v}"));
            original_vectors[i][0] += c * (&shift[v] + &mean);
            original_vectors[i][1] -= c * mean;
            original_vectors[i][v + 3] = c.clone();
        }
    }
    let dot = |i: usize, j: usize| {
        let mut value = Atom::zero();
        for a in 0..vectors {
            for b in 0..vectors {
                value += &original_vectors[i][a] * &original_vectors[j][b] * &formal_gram[a][b];
            }
        }
        value
    };
    let mut replacements = BTreeMap::new();
    let mut index = 0;
    let coordinates = input.basis().coordinates();
    for i in 0..loops {
        for j in i..loops {
            replacements.insert(coordinates[index].clone(), -dot(i, j));
            index += 1;
        }
    }
    let imaginary = Atom::num(symbolica::domains::float::Complex::new(
        Rational::zero(),
        Rational::one(),
    ));
    for i in 0..loops {
        let energy = (0..vectors).fold(Atom::zero(), |s, j| {
            s + &original_vectors[i][j] * &formal_gram[j][2]
        });
        replacements.insert(coordinates[index + i].clone(), &imaginary * energy);
    }
    // Bound the actual distributive expansion before allocating it. The
    // conservative bound counts scalar mean coefficients too, so hidden large
    // sums in a routed vector cannot evade the Gram-monomial budget.
    let replacement_sizes = coordinates
        .iter()
        .map(|v| {
            let expanded = replacements[v].expand();
            match expanded.as_view() {
                AtomView::Add(a) => a.iter().count(),
                _ => 1,
            }
        })
        .collect::<Vec<_>>();
    let mut expansion_bound = 0_usize;
    for (monomial, _) in exact_coefficient_list(&numerator, coordinates)? {
        let exponents = powers(&monomial, coordinates)?;
        let mut term_size = 1_usize;
        for (&power, &size) in exponents.iter().zip(&replacement_sizes) {
            if power < 0 {
                return Err(Error::Unsupported(
                    "nonpolynomial original numerator".into(),
                ));
            }
            for _ in 0..power {
                term_size = term_size
                    .checked_mul(size)
                    .filter(|&v| v <= term_budget)
                    .ok_or_else(|| Error::Limit("Gaussian numerator expansion budget".into()))?;
            }
        }
        expansion_bound = expansion_bound
            .checked_add(term_size)
            .filter(|&v| v <= term_budget)
            .ok_or_else(|| Error::Limit("Gaussian numerator expansion budget".into()))?;
    }
    let numerator = substitute(&numerator, &replacements).expand();
    let terms = exact_coefficient_list(&numerator, &gram_variables)?;
    if terms.len() > term_budget {
        return Err(Error::Limit("Gaussian numerator term budget".into()));
    }
    let h = scalar("h");
    let a = scalar("mass_1");
    let b = scalar("mass_2");
    let e1 = scalar("energy_1");
    let e2 = scalar("energy_2");
    let d = scalar("dimension");
    let cross = (&a + &b + &h) / Atom::num(2);
    let external_gram = vec![
        vec![a, cross.clone(), e1.clone()],
        vec![cross, b, e2.clone()],
        vec![e1, e2, Atom::num(1)],
    ];
    let mut projected = Atom::zero();
    for (monomial, coefficient) in terms {
        let exponents = powers(&monomial, &gram_variables)?;
        let total_rank = exponents.iter().try_fold(0_i64, |s, &p| {
            s.checked_add(2 * i64::from(p))
                .ok_or_else(|| Error::Limit("Gaussian numerator rank overflow".into()))
        })?;
        if total_rank > 32 || exponents.iter().any(|&p| p < 0) {
            return Err(Error::Unsupported(
                "endpoint Gaussian numerator requires a polynomial of vector rank at most 32"
                    .into(),
            ));
        }
        let mut edges = Vec::new();
        for (&p, &edge) in exponents.iter().zip(&gram_pairs) {
            edges.extend(std::iter::repeat_n(edge, p as usize));
        }
        projected += coefficient * wick(&edges, &d, &h, &external_gram, wick_budget)?;
    }
    Ok(projected.expand())
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct JetState {
    energy: [i32; 2],
    mass: [i32; 2],
    h: i32,
    upper: [u16; 2],
}
fn add_state(
    states: &mut BTreeSet<JetState>,
    state: JetState,
    remaining: &mut usize,
) -> Result<()> {
    if states.insert(state) {
        consume(remaining, "mass-jet state")?;
    }
    Ok(())
}
fn mass_jet(
    states: BTreeSet<JetState>,
    loop_index: usize,
    virtual_loops: usize,
    remaining: &mut usize,
) -> Result<BTreeSet<JetState>> {
    let mut out = BTreeSet::new();
    for state in states {
        // Derivative at fixed spatial momentum of the original numerator,
        // shell factor 1/(2E), invariant h, and the moving upper occupation.
        // Factorials and nonzero signs cannot improve these degree bounds.
        if state.mass[loop_index] > 0 {
            let mut next = state.clone();
            next.mass[loop_index] -= 1;
            add_state(&mut out, next, remaining)?;
        }
        if state.energy[loop_index] != 1 {
            let mut next = state.clone();
            next.energy[loop_index] -= 2;
            add_state(&mut out, next, remaining)?;
        }
        // h has a positive D coefficient because there is a virtual loop, so
        // its exponent is not identically zero in the meromorphic family.
        if virtual_loops > 0 || state.h != 0 {
            let mut next = state.clone();
            next.h -= 1;
            add_state(&mut out, next, remaining)?;
            let mut next = state.clone();
            next.h -= 1;
            next.energy[loop_index] -= 1;
            next.energy[1 - loop_index] += 1;
            add_state(&mut out, next, remaining)?;
        }
        let mut next = state.clone();
        next.energy[loop_index] -= 1;
        next.upper[loop_index] = next.upper[loop_index]
            .checked_add(1)
            .ok_or_else(|| Error::Limit("upper occupation order overflow".into()))?;
        add_state(&mut out, next, remaining)?;
    }
    Ok(out)
}

fn compact_envelope(
    projected: &Atom,
    virtual_loops: usize,
    uncut_power_sum: i32,
    cut_powers: [u16; 2],
    remaining: &mut usize,
) -> Result<Vec<CompactMonomialEnvelope>> {
    let variables = [
        scalar("energy_1"),
        scalar("energy_2"),
        scalar("mass_1"),
        scalar("mass_2"),
        scalar("h"),
    ];
    let mut states = BTreeSet::new();
    for (monomial, coefficient) in exact_coefficient_list(projected, &variables)? {
        if coefficient.is_zero() {
            continue;
        }
        let p = powers(&monomial, &variables)?;
        if p.iter().any(|&p| p < 0) {
            return Err(Error::Unsupported("nonpolynomial Gaussian envelope".into()));
        }
        let p = p
            .into_iter()
            .map(|v| i32::try_from(v).map_err(|_| Error::Limit("compact degree overflow".into())))
            .collect::<Result<Vec<_>>>()?;
        add_state(
            &mut states,
            JetState {
                energy: [p[0], p[1]],
                mass: [p[2], p[3]],
                h: p[4] - uncut_power_sum,
                upper: [0, 0],
            },
            remaining,
        )?;
    }
    for (i, &n) in cut_powers.iter().enumerate() {
        for _ in 1..n {
            states = mass_jet(states, i, virtual_loops, remaining)?;
        }
    }
    Ok(states
        .into_iter()
        .filter(|s| s.mass == [0, 0])
        .map(|s| CompactMonomialEnvelope {
            energy_powers: s.energy,
            h_power_offset: s.h,
            upper_indices: s.upper,
        })
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect())
}

/// Generate source-bound evidence from graph/routing/powers/original numerators.
///
/// Scope: two independent future massless cuts, positive chemical magnitudes,
/// every uncut physical denominator present with positive integer power and
/// common positive eta, and an exact affine routing
/// to R*k+b*(q1-q2). All physical masses must vanish. The supported input
/// numerator is the already validated original polynomial, rank at most 32.
/// The resulting strict inequalities are sufficient *power-counting* conditions
/// for the specified virtual and compact strata; no admission flag is enabled.
pub fn single_spacelike_endpoint_evidence(
    input: &PreparedDensityInput,
    cuts: [usize; 2],
    shifted_physical_slots: &[usize],
    budget: EndpointEvidenceBudget,
) -> Result<SingleSpacelikeEndpointEvidence> {
    for target in &input.input().targets {
        if target.numerator.contains("rustflow_massless_endpoint::") {
            return Err(Error::InvalidInput(
                "numerator uses reserved endpoint evidence symbols".into(),
            ));
        }
    }
    let channel = massless_channel_evidence(
        input,
        &cuts,
        shifted_physical_slots,
        budget.cut_sets,
        budget.vertex_partitions,
    )?;
    if !channel.certified {
        return Err(Error::Unsupported(
            "common positive-eta channel evidence is incomplete".into(),
        ));
    }
    if input.physical_masses().iter().any(|m| !m.is_zero()) {
        return Err(Error::Unsupported(
            "single-spacelike endpoint evidence requires all physical masses zero".into(),
        ));
    }
    let family = input.occupied_cut(&cuts, budget.cut_sets)?;
    let loops = input.input().loops;
    let virtual_loops = loops.checked_sub(2).ok_or_else(|| {
        Error::Unsupported("two occupied directions require at least two loops".into())
    })?;
    if loops > 16 {
        return Err(Error::Limit(
            "single-spacelike endpoint loop rank exceeds 16".into(),
        ));
    }
    let input_rows = input
        .input()
        .edges
        .iter()
        .map(|edge| {
            edge.routing
                .iter()
                .map(|s| parse(s))
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    let routed = matmul(&input_rows, family.inverse_routing());
    for row in &routed {
        for entry in row {
            rational(entry)?;
        }
    }
    let uncut = shifted_physical_slots;
    let a = uncut
        .iter()
        .map(|&s| routed[s][2..].to_vec())
        .collect::<Vec<_>>();
    if rank(a.clone()) != virtual_loops {
        return Err(Error::Unsupported(
            "uncut virtual quadratic routing is rank deficient".into(),
        ));
    }
    let b = uncut
        .iter()
        .map(|&s| -(&routed[s][0] + &routed[s][1]))
        .collect::<Vec<_>>();
    let shift = solve(&a, &b, virtual_loops)?;
    let rows = uncut
        .iter()
        .map(|&slot| Row {
            slot,
            virtual_part: routed[slot][2..].to_vec(),
            external: -routed[slot][1].clone(),
        })
        .collect::<Vec<_>>();
    if rows
        .iter()
        .any(|r| r.external.is_zero() && r.virtual_part.iter().all(Atom::is_zero))
    {
        return Err(Error::Unsupported(
            "identically null uncut denominator".into(),
        ));
    }
    let mut targets = Vec::new();
    let mut wick_left = budget.wick_branches;
    let mut jet_left = budget.mass_jet_states;
    for (target_index, target) in input.input().targets.iter().enumerate() {
        if target.powers.iter().any(|&p| p <= 0) {
            return Err(Error::Unsupported(
                "endpoint evidence currently requires positive powers for every physical slot"
                    .into(),
            ));
        }
        let cut_powers = [
            u16::try_from(target.powers[cuts[0]]).unwrap(),
            u16::try_from(target.powers[cuts[1]]).unwrap(),
        ];
        if cut_powers.iter().any(|&p| p > 32) {
            return Err(Error::Limit("endpoint cut order exceeds 32".into()));
        }
        let strata = affine_strata(&rows, &target.powers, budget.affine_subsets)?;
        let numerator = gaussian_envelope(
            input,
            target_index,
            family.inverse_routing(),
            &shift,
            budget.numerator_terms,
            &mut wick_left,
        )?;
        let uncut_power_sum = uncut.iter().map(|&s| i32::from(target.powers[s])).sum();
        let compact = compact_envelope(
            &numerator,
            virtual_loops,
            uncut_power_sum,
            cut_powers,
            &mut jet_left,
        )?;
        let h = i32::try_from(virtual_loops).map_err(|_| Error::Limit("virtual rank".into()))?;
        let mut inequalities = strata
            .iter()
            .map(|s| EndpointInequality {
                stratum: format!("virtual affine soft intersection {:?}", s.physical_slots),
                degree: s.eta_degree.clone(),
            })
            .collect::<Vec<_>>();
        for (j, term) in compact.iter().enumerate() {
            for i in 0..2 {
                if term.upper_indices[i] == 0 {
                    inequalities.push(EndpointInequality {
                        stratum: format!("compact term {j}, radius {i} to zero"),
                        degree: EndpointDegree {
                            twice_dimension_coefficient: h + 2,
                            constant: -2 + term.energy_powers[i] + term.h_power_offset,
                        },
                    });
                }
            }
            inequalities.push(EndpointInequality {
                stratum: format!("compact term {j}, collinear angle"),
                degree: EndpointDegree {
                    twice_dimension_coefficient: h + 1,
                    constant: -1 + term.h_power_offset,
                },
            });
        }
        inequalities.push(EndpointInequality {
            stratum: "opposite angular endpoint".into(),
            degree: EndpointDegree {
                twice_dimension_coefficient: 1,
                constant: -1,
            },
        });
        targets.push(TargetEndpointEvidence {
            target_index,
            original_cut_powers: cut_powers,
            affine_soft_strata: strata,
            compact_envelope: compact,
            inequalities,
        });
    }
    Ok(SingleSpacelikeEndpointEvidence {
        input_identity: input.identity().to_owned(),
        cut_slots: cuts,
        shifted_physical_slots: uncut.to_vec(),
        virtual_loops,
        virtual_uv_continuation_required: virtual_loops > 0,
        virtual_affine_shift: shift.iter().map(Atom::to_canonical_string).collect(),
        denominator_rows: rows
            .iter()
            .map(|r| {
                (
                    r.slot,
                    r.virtual_part
                        .iter()
                        .map(Atom::to_canonical_string)
                        .collect(),
                    r.external.to_canonical_string(),
                )
            })
            .collect(),
        targets,
        prescription: "single-spacelike-degree-evidence-v1;original-mass-jets-before-zero;UV-meromorphic-forest-continuation-required;thermal-origin-and-uniform-remainder-proof-pending",
        endpoint_admitted: false,
        numerical_evaluation_admitted: false,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_density::DensityInput;
    fn e7() -> PreparedDensityInput {
        serde_json::from_str::<DensityInput>(include_str!(
            "../../examples/finite_density/chain_of_three_parallel_pairs.json"
        ))
        .unwrap()
        .prepare()
        .unwrap()
    }
    #[test]
    fn supplied_definitions_generate_strict_e7_degrees_without_period_values() {
        let input = e7();
        let e = single_spacelike_endpoint_evidence(
            &input,
            [0, 4],
            &[1, 2, 3, 5, 6],
            EndpointEvidenceBudget::default(),
        )
        .unwrap();
        assert_eq!(e.virtual_loops, 2);
        assert!(
            e.targets
                .iter()
                .all(|t| t.inequalities_hold(&Rational::from((15, 4))))
        );
        let raised = &e.targets[1];
        assert!(
            raised
                .compact_envelope
                .iter()
                .any(|m| m.energy_powers == [-2, 0] && m.h_power_offset == -3)
        );
        assert!(
            raised
                .compact_envelope
                .iter()
                .any(|m| m.energy_powers == [0, 0] && m.h_power_offset == -4)
        );
        assert!(
            raised
                .compact_envelope
                .iter()
                .any(|m| m.upper_indices == [1, 0])
        );
        assert!(raised.inequalities.iter().any(|x| x.degree
            == EndpointDegree {
                twice_dimension_coefficient: 4,
                constant: -7
            }));
        assert!(raised.inequalities.iter().any(|x| x.degree
            == EndpointDegree {
                twice_dimension_coefficient: 3,
                constant: -5
            }));
        assert!(!raised.inequalities_hold(&Rational::from((7, 2))));
        assert!(!e.endpoint_admitted && !e.numerical_evaluation_admitted);
    }
    #[test]
    fn pure_compact_transfer_needs_no_virtual_uv_continuation() {
        let mut definition = serde_json::from_str::<DensityInput>(include_str!(
            "../../examples/finite_density/massive_two_loop_sunset.json"
        ))
        .unwrap();
        for edge in &mut definition.edges {
            edge.mass_squared = "0".into();
        }
        let input = definition.prepare().unwrap();
        let e = single_spacelike_endpoint_evidence(
            &input,
            [0, 1],
            &[2],
            EndpointEvidenceBudget::default(),
        )
        .unwrap();
        assert_eq!(e.virtual_loops, 0);
        assert!(!e.virtual_uv_continuation_required);
        assert!(
            e.targets
                .iter()
                .all(|t| t.inequalities_hold(&Rational::from(7)))
        );
        assert!(e.targets.iter().all(|t| t.affine_soft_strata.is_empty()));
        assert!(!e.endpoint_admitted);
    }
    #[test]
    fn independent_compact_directions_are_rejected_by_exact_rank() {
        let input = serde_json::from_str::<DensityInput>(include_str!(
            "../../examples/finite_density/triangular_prism.json"
        ))
        .unwrap()
        .prepare()
        .unwrap();
        let error = single_spacelike_endpoint_evidence(
            &input,
            [1, 5],
            &[0, 2, 3, 4, 6, 7, 8],
            EndpointEvidenceBudget::default(),
        )
        .unwrap_err();
        assert!(
            matches!(error,Error::Unsupported(ref s) if s.contains("single spacelike direction"))
        );
    }
    #[test]
    fn orientation_and_affine_routing_do_not_supply_a_topology_shortcut() {
        let mut definition = e7().input().clone();
        // Reverse an occupied edge and a neutral edge together with incidence.
        for slot in [0, 1] {
            definition.edges[slot].vertices.swap(0, 1);
            for c in &mut definition.edges[slot].routing {
                *c = (-parse(c).unwrap()).to_string();
            }
            for q in &mut definition.edges[slot].charges {
                *q = -*q;
            }
        }
        // A genuine loop shear old P2=new P2+new P1 forces a nonzero
        // virtual affine shift; transform charge and the original polynomial.
        for edge in &mut definition.edges {
            edge.routing[0] =
                (&parse(&edge.routing[0]).unwrap() + &parse(&edge.routing[1]).unwrap()).to_string();
        }
        definition.loop_charges[1][0] = -1;
        let rules = BTreeMap::from([
            (parse("g1_2").unwrap(), parse("g1_2+g1_1").unwrap()),
            (parse("g2_4").unwrap(), parse("g2_4+g1_4").unwrap()),
        ]);
        for target in &mut definition.targets {
            target.numerator = substitute(&parse(&target.numerator).unwrap(), &rules).to_string();
        }
        definition.name = "unrelated_name".into();
        let input = definition.prepare().unwrap();
        let e = single_spacelike_endpoint_evidence(
            &input,
            [0, 4],
            &[1, 2, 3, 5, 6],
            EndpointEvidenceBudget::default(),
        )
        .unwrap();
        assert!(
            e.targets
                .iter()
                .all(|t| t.inequalities_hold(&Rational::from((15, 4))))
        );
        assert!(e.virtual_affine_shift.iter().any(|s| s != "0"));
        let mut budget = EndpointEvidenceBudget::default();
        budget.affine_subsets = 1;
        assert!(matches!(
            single_spacelike_endpoint_evidence(&input, [0, 4], &[1, 2, 3, 5, 6], budget),
            Err(Error::Limit(_))
        ));
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct OneLoopNullConeTarget {
    pub target_index: usize,
    pub uncut_positive_power_sum: i32,
    pub original_cut_mass_jet_order: u16,
    pub additional_energy_jet_order: u16,
    pub original_polynomial_vector_rank: u32,
    pub inequalities: Vec<EndpointInequality>,
}

#[derive(Clone, Debug, Serialize)]
pub struct OneLoopNullConeEvidence {
    pub input_identity: String,
    pub cut_slot: usize,
    pub shifted_physical_slots: Vec<usize>,
    /// Exact coefficients in k_e=a_e*k+b_e*q, with q future/null at the endpoint.
    pub denominator_rows: Vec<(usize, String, String)>,
    /// On the positive simplex A=sum a_e²*x_e >= this strictly positive bound.
    pub gaussian_a_lower_bound: String,
    /// 0<=Q=sum b_e²*x_e-(sum a_e*b_e*x_e)²/A <= this finite bound.
    pub external_q_upper_bound: String,
    pub chemical_magnitude: String,
    /// A nonempty proof domain witness, not a restriction on the eventual
    /// dimension reached by meromorphic continuation. Its fractional part also
    /// avoids every integer-shifted Gaussian Gamma pole in this finite family.
    pub high_dimension_witness: String,
    pub gaussian_gamma_poles_retained: bool,
    pub thermal_limit_order: &'static str,
    pub targets: Vec<OneLoopNullConeTarget>,
    pub prescription: &'static str,
    pub endpoint_admitted: bool,
    pub numerical_evaluation_admitted: bool,
}

fn original_polynomial_rank(input: &PreparedDensityInput, index: usize) -> Result<u32> {
    let target = &input.input().targets[index];
    let mut numerator = input.basis().reconstruct(&input.targets()[index])?;
    for (factor, &power) in input.basis().slots().iter().zip(&target.powers) {
        numerator *= factor.clone().pow(i64::from(power));
    }
    let numerator = numerator.together().cancel().expand();
    let scalar_count = input.input().loops * (input.input().loops + 1) / 2;
    let mut maximum = 0_u32;
    for (monomial, _) in exact_coefficient_list(&numerator, input.basis().coordinates())? {
        let p = powers(&monomial, input.basis().coordinates())?;
        if p.iter().any(|&n| n < 0) {
            return Err(Error::Unsupported(
                "null cone evidence needs the original polynomial numerator".into(),
            ));
        }
        let degree = p.iter().enumerate().try_fold(0_u32, |d, (i, &n)| {
            d.checked_add(u32::from(n as u16) * if i < scalar_count { 2 } else { 1 })
                .ok_or_else(|| Error::Limit("null cone numerator rank overflow".into()))
        })?;
        maximum = maximum.max(degree);
    }
    if maximum > 32 {
        return Err(Error::Limit(
            "null cone numerator vector rank exceeds 32".into(),
        ));
    }
    Ok(maximum)
}

/// Exact sufficient bounds for a single future null occupied momentum and one
/// virtual loop, with every uncut physical momentum depending on that loop.
/// Nonpositive ordinary indices are finite polynomial insertions. At least one
/// ordinary denominator must remain active; energy inverses are not admitted.
///
/// The polynomial numerator is retained, and no Gaussian period is evaluated.
/// The Gaussian parameter identity is first defined at UV-convergent complex D,
/// then continued meromorphically. Away from its Gamma poles, the finite tensor
/// terms have powers lambda=D/2-P+k, k>=0. The required J mass derivatives have
/// lambda-j, j<=J. Re D>2(P+J) makes all of these powers strictly positive.
/// Re D>2J+2 additionally controls the radial shell Jacobian derivatives.
/// These are target-level bounds; arbitrary generated source indices need their
/// own bounds and the shared massless origin prescription before flow admission.
pub fn one_loop_null_cone_evidence(
    input: &PreparedDensityInput,
    cut_slot: usize,
    shifted_physical_slots: &[usize],
    budget: EndpointEvidenceBudget,
) -> Result<OneLoopNullConeEvidence> {
    let channel = massless_channel_evidence(
        input,
        &[cut_slot],
        shifted_physical_slots,
        budget.cut_sets,
        budget.vertex_partitions,
    )?;
    if !channel.certified {
        return Err(Error::Unsupported(
            "null cone needs common positive-eta channel evidence".into(),
        ));
    }
    if input.input().loops != 2 || input.physical_masses().iter().any(|m| !m.is_zero()) {
        return Err(Error::Unsupported("one-loop null cone scope is exactly one occupied loop, one virtual loop and zero physical masses".into()));
    }
    let family = input.occupied_cut(&[cut_slot], budget.cut_sets)?;
    let routing = input
        .input()
        .edges
        .iter()
        .map(|edge| {
            edge.routing
                .iter()
                .map(|s| parse(s))
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    let rows = matmul(&routing, family.inverse_routing());
    let mut a_lower: Option<Rational> = None;
    let mut q_upper = Rational::zero();
    let mut denominator_rows = Vec::new();
    for &slot in shifted_physical_slots {
        let a = rational(&rows[slot][1])?;
        let b = rational(&rows[slot][0])?;
        if a.is_zero() {
            return Err(Error::Unsupported(
                "null cone bound excludes an uncut factor independent of its virtual loop".into(),
            ));
        }
        let a2 = &a * &a;
        let b2 = &b * &b;
        a_lower = Some(a_lower.map_or(a2.clone(), |old| old.min(a2)));
        q_upper = q_upper.max(b2);
        denominator_rows.push((slot, a.to_string(), b.to_string()));
    }
    let mut targets = Vec::new();
    for (target_index, target) in input.input().targets.iter().enumerate() {
        if target.powers[cut_slot] <= 0 {
            return Err(Error::Unsupported(
                "null cone target needs a positive occupied power".into(),
            ));
        }
        let j = u16::try_from(target.powers[cut_slot] - 1).unwrap();
        if j >= 32 {
            return Err(Error::Limit("null cone cut order exceeds 32".into()));
        }
        let p = shifted_physical_slots
            .iter()
            .map(|&s| i32::from(target.powers[s].max(0)))
            .sum::<i32>();
        if p == 0 {
            return Err(Error::Unsupported(
                "no active virtual denominator: use the polynomial terminal owner".into(),
            ));
        }
        let polynomial_insertions = shifted_physical_slots
            .iter()
            .map(|&s| 2 * u32::from(target.powers[s].min(0).unsigned_abs()))
            .sum::<u32>();
        let rank = original_polynomial_rank(input, target_index)?
            .checked_add(polynomial_insertions)
            .ok_or_else(|| Error::Limit("null cone polynomial insertion rank overflow".into()))?;
        if rank > 32 {
            return Err(Error::Limit(
                "null cone polynomial vector rank exceeds 32".into(),
            ));
        }
        targets.push(OneLoopNullConeTarget {
            target_index,
            uncut_positive_power_sum: p,
            original_cut_mass_jet_order: j,
            additional_energy_jet_order: 0,
            original_polynomial_vector_rank: rank,
            inequalities: vec![
                EndpointInequality {
                    stratum: "uniform Gaussian mass/thermal-cone jets".into(),
                    degree: EndpointDegree {
                        twice_dimension_coefficient: 1,
                        constant: -p - i32::from(j),
                    },
                },
                EndpointInequality {
                    stratum: "compact origin after shell Jacobian jets".into(),
                    degree: EndpointDegree {
                        twice_dimension_coefficient: 2,
                        constant: -2 - 2 * i32::from(j),
                    },
                },
            ],
        });
    }
    let largest = targets
        .iter()
        .map(|t| {
            t.uncut_positive_power_sum
                + i32::from(t.original_cut_mass_jet_order)
                + i32::from(t.additional_energy_jet_order)
        })
        .max()
        .unwrap();
    let witness = Rational::from((4 * largest + 1, 2));
    Ok(OneLoopNullConeEvidence {
        input_identity: input.identity().to_owned(),
        cut_slot,
        shifted_physical_slots: shifted_physical_slots.to_vec(),
        denominator_rows,
        gaussian_a_lower_bound: a_lower
            .ok_or_else(|| Error::Unsupported("null cone has no virtual factors".into()))?
            .to_string(),
        external_q_upper_bound: q_upper.to_string(),
        chemical_magnitude: family.shells()[0].chemical_potential.to_string(),
        high_dimension_witness: witness.to_string(),
        gaussian_gamma_poles_retained: true,
        thermal_limit_order: "remove-line-energy-regulator-at-fixed-positive-T-and-eta;then-real-Fermi-T0-and-eta0;then-meromorphic-D-continuation",
        targets,
        prescription: "one-loop-null-cone-v1;Gaussian-UV-meromorphic-continuation;common-eta-and-thermal-energy-cone;energy-regulator-removed-at-fixed-T;original-mass-jets-before-zero;source-origin-admission-separate",
        endpoint_admitted: false,
        numerical_evaluation_admitted: false,
    })
}

#[cfg(test)]
mod cone_tests {
    use super::*;
    use crate::finite_density::DensityInput;
    #[test]
    fn massless_small_graph_provides_a_uniform_null_cone() {
        let mut definition = serde_json::from_str::<DensityInput>(include_str!(
            "../../examples/finite_density/massive_two_loop_sunset.json"
        ))
        .unwrap();
        for edge in &mut definition.edges {
            edge.mass_squared = "0".into();
        }
        let input = definition.prepare().unwrap();
        for slot in [0, 1] {
            let shifted = (0..3).filter(|&i| i != slot).collect::<Vec<_>>();
            let e = one_loop_null_cone_evidence(
                &input,
                slot,
                &shifted,
                EndpointEvidenceBudget::default(),
            )
            .unwrap();
            assert_eq!(e.gaussian_a_lower_bound, "1");
            assert_eq!(e.external_q_upper_bound, "1");
            assert!(
                e.targets
                    .iter()
                    .flat_map(|t| &t.inequalities)
                    .all(|x| x.degree.at(&Rational::from((13, 2))) > Rational::zero())
            );
            assert!(!e.endpoint_admitted && !e.numerical_evaluation_admitted);
        }
    }
    #[test]
    fn nonpositive_ordinary_powers_only_enlarge_the_polynomial_insertions() {
        let mut definition = serde_json::from_str::<DensityInput>(include_str!(
            "../../examples/finite_density/massive_two_loop_sunset.json"
        ))
        .unwrap();
        for edge in &mut definition.edges {
            edge.mass_squared = "0".into();
        }
        definition.targets.truncate(1);
        definition.targets[0].powers = vec![2, -1, 1];
        let input = definition.prepare().unwrap();
        let e = one_loop_null_cone_evidence(&input, 0, &[1, 2], EndpointEvidenceBudget::default())
            .unwrap();
        assert_eq!(e.targets[0].uncut_positive_power_sum, 1);
        assert_eq!(e.targets[0].original_cut_mass_jet_order, 1);
        assert_eq!(e.targets[0].original_polynomial_vector_rank, 2);
        definition.targets[0].powers = vec![2, -1, 0];
        let input = definition.prepare().unwrap();
        assert!(
            matches!(one_loop_null_cone_evidence(&input,0,&[1,2],EndpointEvidenceBudget::default()),
            Err(Error::Unsupported(ref s)) if s.contains("no active virtual denominator"))
        );
    }
    #[test]
    fn extra_virtual_loops_cannot_inherit_the_one_loop_simplex_bound() {
        let input = serde_json::from_str::<DensityInput>(include_str!(
            "../../examples/finite_density/chain_of_three_parallel_pairs.json"
        ))
        .unwrap()
        .prepare()
        .unwrap();
        let error = one_loop_null_cone_evidence(
            &input,
            0,
            &[1, 2, 3, 4, 5, 6],
            EndpointEvidenceBudget::default(),
        )
        .unwrap_err();
        assert!(
            matches!(error,Error::Unsupported(ref s) if s.contains("exactly one occupied loop"))
        );
    }
}

/// One finite weighted label's endpoint/origin proof. The high-D witness is an
/// open-set witness before meromorphic continuation, not a sampling restriction.
#[derive(Clone, Debug, Serialize)]
pub struct MasslessLabelBound {
    pub indices: Vec<i16>,
    pub classification: &'static str,
    pub uncut_positive_power_sum: i32,
    pub active_virtual_rank: usize,
    pub virtual_block_positive_power_sums: Vec<i32>,
    pub pure_transfer_positive_power_sum: i32,
    pub mass_and_upper_jet_orders: Vec<i32>,
    pub inequalities: Vec<EndpointInequality>,
    pub high_dimension_witness: String,
}
#[derive(Clone, Debug, Serialize)]
pub struct MasslessLabelAudit {
    pub endpoint_structure: MasslessStructuralAudit,
    pub origin_identity: &'static str,
    pub labels: Vec<MasslessLabelBound>,
}

pub(crate) const MASSLESS_FLOW_ORIGIN_IDENTITY: &str = "joint-high-D-massless-flow-origin-v1;polynomial-completions;lower-contact-zero-jets;energy-regulator-removed-at-fixed-positive-T-and-eta;real-Fermi-T0;endpoint-then-meromorphic-D";

/// Exact routing evidence retained by the sealed endpoint permit. Bounds use
/// Gaussian parameters and no evaluated virtual period.
#[derive(Clone, Debug, Serialize)]
pub struct MasslessStructuralAudit {
    pub proof_version: &'static str,
    pub virtual_loops: usize,
    pub virtual_affine_shift: Vec<String>,
    /// Physical slot, native virtual row, external coefficient of q or q1-q2.
    pub denominator_rows: Vec<(usize, Vec<String>, String)>,
    /// None denotes a pure spacelike transfer denominator.
    pub block_assignments: Vec<(usize, Option<usize>)>,
    pub block_representatives: Vec<Vec<String>>,
    pub block_a_lower_bounds: Vec<String>,
    pub block_q_upper_bounds: Vec<String>,
}
const MASSLESS_ENDPOINT_PROOF_VERSION: &str = "positive-eta-gaussian-germ-and-rank-one-blocks-v1;fixed-T-regulator-order;UV-meromorphic-before-high-D-endpoint";

#[derive(Clone, Debug)]
struct CertifiedRow {
    slot: usize,
    virtual_part: Vec<Atom>,
    external: Atom,
}
#[derive(Clone, Debug)]
enum CertifiedEndpointClass {
    SingletonGerm {
        virtual_loops: usize,
        rows: Vec<CertifiedRow>,
        q_upper: Rational,
    },
    RankOneBlocks {
        virtual_loops: usize,
        rows: Vec<CertifiedRow>,
        affine_shift: Vec<Atom>,
        assignments: Vec<Option<usize>>,
        representatives: Vec<Vec<Atom>>,
        a_lower: Vec<Rational>,
        q_upper: Vec<Rational>,
    },
}
impl CertifiedEndpointClass {
    fn virtual_loops(&self) -> usize {
        match self {
            Self::SingletonGerm { virtual_loops, .. }
            | Self::RankOneBlocks { virtual_loops, .. } => *virtual_loops,
        }
    }
    fn rows(&self) -> &[CertifiedRow] {
        match self {
            Self::SingletonGerm { rows, .. } | Self::RankOneBlocks { rows, .. } => rows,
        }
    }
    fn structure(&self) -> MasslessStructuralAudit {
        let mut evidence = MasslessStructuralAudit {
            proof_version: MASSLESS_ENDPOINT_PROOF_VERSION,
            virtual_loops: self.virtual_loops(),
            virtual_affine_shift: Vec::new(),
            denominator_rows: self
                .rows()
                .iter()
                .map(|r| {
                    (
                        r.slot,
                        r.virtual_part
                            .iter()
                            .map(Atom::to_canonical_string)
                            .collect(),
                        r.external.to_canonical_string(),
                    )
                })
                .collect(),
            block_assignments: Vec::new(),
            block_representatives: Vec::new(),
            block_a_lower_bounds: Vec::new(),
            block_q_upper_bounds: Vec::new(),
        };
        match self {
            Self::SingletonGerm { q_upper, .. } => {
                evidence.block_q_upper_bounds.push(q_upper.to_string());
            }
            Self::RankOneBlocks {
                rows,
                affine_shift,
                assignments,
                representatives,
                a_lower,
                q_upper,
                ..
            } => {
                evidence.virtual_affine_shift =
                    affine_shift.iter().map(Atom::to_canonical_string).collect();
                evidence.block_assignments = rows
                    .iter()
                    .zip(assignments)
                    .map(|(r, b)| (r.slot, *b))
                    .collect();
                evidence.block_representatives = representatives
                    .iter()
                    .map(|row| row.iter().map(Atom::to_canonical_string).collect())
                    .collect();
                evidence.block_a_lower_bounds = a_lower.iter().map(ToString::to_string).collect();
                evidence.block_q_upper_bounds = q_upper.iter().map(ToString::to_string).collect();
            }
        }
        evidence
    }
}

/// Rational quadratic separation: after shifting k by an exact multiple of
/// q1, every uncut row depends on q1-q2 and at most one independent virtual
/// coordinate. The test also includes h=0 pure compact transfers.
fn rank_one_blocks(
    routed: &[Vec<Atom>],
    shifted: &[usize],
    h: usize,
) -> Result<CertifiedEndpointClass> {
    let a = shifted
        .iter()
        .map(|&s| routed[s][2..].to_vec())
        .collect::<Vec<_>>();
    if rank(a.clone()) != h {
        return Err(Error::Unsupported(
            "uncut rank-one block routing is rank deficient".into(),
        ));
    }
    let b = shifted
        .iter()
        .map(|&s| -(&routed[s][0] + &routed[s][1]))
        .collect::<Vec<_>>();
    let affine_shift = solve(&a, &b, h)?;
    let rows = shifted
        .iter()
        .map(|&slot| CertifiedRow {
            slot,
            virtual_part: routed[slot][2..].to_vec(),
            external: -routed[slot][1].clone(),
        })
        .collect::<Vec<_>>();
    let mut representatives: Vec<Vec<Atom>> = Vec::new();
    let mut assignments = Vec::new();
    let mut a_lower: Vec<Rational> = Vec::new();
    let mut q_upper: Vec<Rational> = Vec::new();
    for row in &rows {
        let Some(first) = row.virtual_part.iter().position(|x| !x.is_zero()) else {
            if row.external.is_zero() {
                return Err(Error::Unsupported(
                    "identically null uncut denominator".into(),
                ));
            }
            assignments.push(None);
            continue;
        };
        let scale = rational(&row.virtual_part[first])?;
        let normalized = row
            .virtual_part
            .iter()
            .map(|x| x / &row.virtual_part[first])
            .collect::<Vec<_>>();
        let a2 = &scale * &scale;
        let external = rational(&row.external)?;
        let b2 = &external * &external;
        let block = if let Some(i) = representatives.iter().position(|r| r == &normalized) {
            a_lower[i] = a_lower[i].clone().min(a2);
            q_upper[i] = q_upper[i].clone().max(b2);
            i
        } else {
            let i = representatives.len();
            representatives.push(normalized);
            a_lower.push(a2);
            q_upper.push(b2);
            i
        };
        assignments.push(Some(block));
    }
    if representatives.len() != h || rank(representatives.clone()) != h {
        return Err(Error::Unsupported(
            "virtual quadratics do not form independent rank-one blocks".into(),
        ));
    }
    Ok(CertifiedEndpointClass::RankOneBlocks {
        virtual_loops: h,
        rows,
        affine_shift,
        assignments,
        representatives,
        a_lower,
        q_upper,
    })
}

/// Measure-owner proof that an admitted positive ordinary support leaves an
/// unrestricted polynomial virtual direction. Fields are sealed; compact rows
/// never occur in the rank calculation.
#[derive(Clone, Debug, Serialize)]
pub struct FreeVirtualZeroSupport {
    proof_version: &'static str,
    virtual_rank: usize,
    allowed_positive_slots: Vec<usize>,
    forced_nonpositive_slots: Vec<usize>,
    /// Exact rational vector, normalized to first nonzero coefficient one.
    virtual_null_direction: Vec<String>,
}
impl FreeVirtualZeroSupport {
    pub fn forced_nonpositive_slots(&self) -> &[usize] {
        &self.forced_nonpositive_slots
    }
    pub fn identity(&self) -> String {
        format!(
            "{};virtual-rank={};allowed={:?};nonpositive={:?};virtual-null={:?}",
            self.proof_version,
            self.virtual_rank,
            self.allowed_positive_slots,
            self.forced_nonpositive_slots,
            self.virtual_null_direction
        )
    }
}
const FREE_VIRTUAL_ZERO_PROOF_VERSION: &str = "occupied-free-virtual-polynomial-zero-v1;fixed-eta;dimensional-virtual-measure;bound-origin;polynomial-completions";

fn checked_flat_subsets(n: usize, k: usize, budget: usize) -> Result<usize> {
    if k > n {
        return Ok(0);
    }
    let k = k.min(n - k);
    let mut count = 1usize;
    for i in 0..k {
        count = count
            .checked_mul(n - i)
            .ok_or_else(|| Error::Limit("free-virtual flat subset count overflow".into()))?
            / (i + 1);
        if count > budget {
            return Err(Error::Limit(format!(
                "free-virtual flat enumeration exceeds subset budget {budget}"
            )));
        }
    }
    if count > budget {
        return Err(Error::Limit(
            "free-virtual flat enumeration has zero budget".into(),
        ));
    }
    Ok(count)
}
fn free_virtual_flat(
    rows: &[CertifiedRow],
    h: usize,
    selected: Vec<Vec<Atom>>,
) -> Result<Option<FreeVirtualZeroSupport>> {
    let selected = if selected.is_empty() {
        vec![vec![Atom::zero(); h]]
    } else {
        selected
    };
    let (reduced, pivots) = rref(selected);
    if pivots.len() != h - 1 {
        return Ok(None);
    }
    let free = (0..h)
        .find(|j| !pivots.contains(j))
        .ok_or_else(|| Error::InvalidInput("free-virtual flat has no null direction".into()))?;
    let mut v = vec![Atom::zero(); h];
    v[free] = Atom::one();
    for (r, &pivot) in reduced.iter().zip(&pivots) {
        v[pivot] = -r[free].clone();
    }
    let normalization = v.iter().find(|x| !x.is_zero()).unwrap().clone();
    for x in &mut v {
        *x = (&*x / &normalization).together().cancel();
        rational(x)?;
    }
    let mut allowed = Vec::new();
    let mut forced = Vec::new();
    let mut flat_rows = Vec::new();
    for row in rows {
        let dot = row
            .virtual_part
            .iter()
            .zip(&v)
            .fold(Atom::zero(), |sum, (a, b)| sum + a * b)
            .together()
            .cancel();
        rational(&dot)?;
        if dot.is_zero() {
            allowed.push(row.slot);
            flat_rows.push(row.virtual_part.clone());
        } else {
            forced.push(row.slot);
        }
    }
    if forced.is_empty() || rank(flat_rows) != h - 1 {
        return Err(Error::InvalidInput(
            "free-virtual flat witness failed exact rank verification".into(),
        ));
    }
    allowed.sort_unstable();
    forced.sort_unstable();
    Ok(Some(FreeVirtualZeroSupport {
        proof_version: FREE_VIRTUAL_ZERO_PROOF_VERSION,
        virtual_rank: h,
        allowed_positive_slots: allowed,
        forced_nonpositive_slots: forced,
        virtual_null_direction: v.iter().map(Atom::to_canonical_string).collect(),
    }))
}

/// Sealed proof for singleton UV-meromorphic germs and independent rank-one
/// virtual blocks with complete high-D endpoint arguments. General single-spacelike *degree evidence* cannot construct this
/// permit. All fields are private, and every consumer rechecks family binding.
#[derive(Clone, Debug)]
pub struct MasslessFlowEvidence {
    input_identity: String,
    cuts: Vec<usize>,
    shifted: Vec<usize>,
    family_signature: Vec<String>,
    source_options: super::preparation::WeightedSourceOptions,
    class: CertifiedEndpointClass,
}

/// Proof data for one unraised massless occupied shell. Construction is only
/// through the bound SingletonGerm class, never a string classification.
/// This does not permit inverse energy powers. All fields remain private.
#[derive(Clone, Debug)]
pub(crate) struct SingletonRawWardEvidence {
    identity: String,
    coordinates: Vec<Atom>,
    factors: Vec<Atom>,
    roles: Vec<super::guarded::IndexRole>,
    cut: usize,
    upper: usize,
    lower: usize,
    physical_slots: usize,
    input_slots: usize,
    energy: Atom,
    energy_dependent_slots: Vec<usize>,
}

impl SingletonRawWardEvidence {
    pub(crate) fn identity(&self) -> &str {
        &self.identity
    }
    pub(crate) fn energy(&self) -> &Atom {
        &self.energy
    }
    pub(crate) fn cut(&self) -> usize {
        self.cut
    }
    pub(crate) fn upper(&self) -> usize {
        self.upper
    }
    pub(crate) fn lower(&self) -> usize {
        self.lower
    }
    pub(crate) fn physical_slots(&self) -> usize {
        self.physical_slots
    }
    pub(crate) fn input_slots(&self) -> usize {
        self.input_slots
    }
    pub(crate) fn energy_dependent_slots(&self) -> &[usize] {
        &self.energy_dependent_slots
    }
    pub(crate) fn validate_measure<const N: usize>(
        &self,
        coordinates: &[Atom],
        factors: &[Atom],
        roles: &[super::guarded::IndexRole; N],
    ) -> Result<()> {
        if self.factors.len() > N
            || coordinates != self.coordinates
            || factors != self.factors
            || roles[..self.roles.len()] != self.roles
            || roles[self.roles.len()..]
                .iter()
                .any(|r| *r != super::guarded::IndexRole::Ordinary)
        {
            return Err(Error::InvalidInput(
                "raw Ward certificate does not match the exact measure".into(),
            ));
        }
        Ok(())
    }
}

pub(crate) fn bound_family_signature(family: &super::geometry::OccupiedCutFamily) -> Vec<String> {
    let mut signature = vec![
        format!(
            "loops={};physical={};input={}",
            family.loops(),
            family.physical_slots(),
            family.input_slots()
        ),
        format!("det={}", family.routing_determinant()),
    ];
    signature.extend(family.coordinates().iter().map(Atom::to_canonical_string));
    signature.extend(family.factors().iter().map(Atom::to_canonical_string));
    signature.extend(family.roles().iter().map(|r| format!("role={r:?}")));
    for row in family.inverse_routing() {
        signature.extend(row.iter().map(Atom::to_canonical_string));
    }
    for s in family.shells() {
        signature.push(format!(
            "shell={}:{}:{}:{}:{}:{}",
            s.loop_index,
            s.physical_slot,
            s.upper_slot,
            s.lower_slot,
            s.mass_squared.to_canonical_string(),
            s.chemical_potential
        ));
    }
    for (target, terms) in family.targets().iter().enumerate() {
        for (label, coefficient) in terms {
            signature.push(format!(
                "target={target}:{:?}:{}",
                label.0,
                coefficient.to_canonical_string()
            ));
        }
    }
    signature
}

fn regular_dimension_witness(
    inequalities: &[EndpointInequality],
    virtual_loops: usize,
) -> Result<Rational> {
    let mut integer_bound = 0_i64;
    for inequality in inequalities {
        let a = i64::from(inequality.degree.twice_dimension_coefficient);
        if a <= 0 {
            return Err(Error::Unsupported(
                "endpoint proof has no verified high-D direction".into(),
            ));
        }
        let b = -2 * i64::from(inequality.degree.constant);
        integer_bound = integer_bound.max(if b > 0 { (b + a - 1) / a } else { 0 });
    }
    // kD/2 is nonintegral for every 1<=k<=h. Half-integer D would
    // land on the overall Gaussian poles for h=4,8,... .
    let denominator = 2 * i64::try_from(virtual_loops.max(1))
        .map_err(|_| Error::Limit("massless endpoint witness rank overflow".into()))?
        + 1;
    Ok(Rational::from((
        denominator * integer_bound + 1,
        denominator,
    )))
}

impl MasslessFlowEvidence {
    /// Bind the actual assigned family and all-uncut common-eta deformation.
    /// Numerator/source inverse energies are excluded even if an option asks
    /// for them; the high-D origin proof requires finite polynomial insertions.
    pub fn new(
        input: &PreparedDensityInput,
        family: &super::geometry::OccupiedCutFamily,
        shifted: &[usize],
        source_options: super::preparation::WeightedSourceOptions,
    ) -> Result<Self> {
        if input.input().loops < 2 || family.loops() != input.input().loops {
            return Err(Error::Unsupported(
                "massless flowing endpoint needs at least two loops and matching family rank"
                    .into(),
            ));
        }
        if family.loops() > 16 {
            return Err(Error::Limit(
                "massless endpoint loop rank exceeds 16".into(),
            ));
        }
        if source_options.positive_compact_energy_powers {
            return Err(Error::Unsupported(
                "massless flow origin excludes inverse compact-energy completions".into(),
            ));
        }
        if input.physical_masses().iter().any(|m| !m.is_zero()) {
            return Err(Error::Unsupported(
                "massless flow permit requires all physical masses zero".into(),
            ));
        }
        let cuts = family
            .shells()
            .iter()
            .map(|s| s.physical_slot)
            .collect::<Vec<_>>();
        let expected = input
            .occupied_cut(&cuts, EndpointEvidenceBudget::default().cut_sets)?
            .at_physical_masses();
        let signature = bound_family_signature(family);
        if signature != bound_family_signature(&expected) {
            return Err(Error::InvalidInput(
                "massless flow proof is not bound to the supplied physical family".into(),
            ));
        }
        let channels = massless_channel_evidence(input, &cuts, shifted, 1024, 65536)?;
        if !channels.certified {
            return Err(Error::Unsupported(
                "massless flow lacks common positive-eta channel evidence".into(),
            ));
        }
        if family
            .shells()
            .iter()
            .any(|s| !s.mass_squared.is_zero() || s.chemical_potential <= 0)
        {
            return Err(Error::Unsupported(
                "massless flow requires zero shell masses and positive chemical magnitudes".into(),
            ));
        }
        let routing = input
            .input()
            .edges
            .iter()
            .map(|edge| {
                edge.routing
                    .iter()
                    .map(|s| parse(s))
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let rows = matmul(&routing, family.inverse_routing());
        for row in &rows {
            for entry in row {
                rational(entry)?;
            }
        }
        let class = match cuts.len() {
            1 => {
                let h = family.loops() - 1;
                let virtual_rows = shifted
                    .iter()
                    .map(|&slot| rows[slot][1..].to_vec())
                    .collect::<Vec<_>>();
                if virtual_rows.iter().any(|r| r.iter().all(Atom::is_zero)) {
                    return Err(Error::Unsupported(
                        "singleton germ excludes external-only uncut physical factors".into(),
                    ));
                }
                if rank(virtual_rows.clone()) != h {
                    return Err(Error::Unsupported(
                        "singleton germ uncut virtual rows do not span the virtual space".into(),
                    ));
                }
                let mut q_upper = Rational::zero();
                for &slot in shifted {
                    let b = rational(&rows[slot][0])?;
                    q_upper = q_upper.max(&b * &b);
                }
                CertifiedEndpointClass::SingletonGerm {
                    virtual_loops: h,
                    rows: shifted
                        .iter()
                        .map(|&slot| CertifiedRow {
                            slot,
                            virtual_part: rows[slot][1..].to_vec(),
                            external: rows[slot][0].clone(),
                        })
                        .collect(),
                    q_upper,
                }
            }
            2 => rank_one_blocks(&rows, shifted, family.loops() - 2)?,
            _ => {
                return Err(Error::Unsupported(
                    "massless flow endpoint covers one or two occupied loops".into(),
                ));
            }
        };
        let evidence = Self {
            input_identity: input.identity().to_owned(),
            cuts,
            shifted: shifted.to_vec(),
            family_signature: signature,
            source_options,
            class,
        };
        let targets = family
            .targets()
            .iter()
            .flat_map(|t| t.keys().cloned())
            .collect::<Vec<_>>();
        evidence.validate_labels(family, &targets)?;
        Ok(evidence)
    }

    /// A narrowly guarded finite-eta source theorem, separate from the endpoint
    /// theorem. Its source bases have C1, lower H0, and zero indices for every
    /// medium-dependent factor. The row may contain polynomial energy images.
    pub(crate) fn singleton_raw_ward_evidence<const N: usize>(
        &self,
        family: &super::geometry::OccupiedCutFamily,
        shifted: &[usize],
        eta: Symbol,
    ) -> Result<SingletonRawWardEvidence> {
        use super::guarded::IndexRole;
        self.validate_family(family, shifted)?;
        if !matches!(&self.class, CertifiedEndpointClass::SingletonGerm { .. })
            || family.shells().len() != 1
            || self.source_options.positive_compact_energy_powers
        {
            return Err(Error::Unsupported(
                "raw Ward source needs the polynomial singleton-germ class".into(),
            ));
        }
        let shell = &family.shells()[0];
        if !shell.mass_squared.is_zero() || shell.chemical_potential <= Rational::zero() {
            return Err(Error::Unsupported(
                "raw Ward source needs mass zero and positive mu".into(),
            ));
        }
        if self.certified_origin_loops(family, shifted)? != vec![shell.loop_index] {
            return Err(Error::InvalidInput(
                "raw Ward source lacks its bound origin prescription".into(),
            ));
        }
        let loops = family.loops();
        let scalar_count = loops
            .checked_add(1)
            .and_then(|n| loops.checked_mul(n))
            .map(|n| n / 2)
            .ok_or_else(|| Error::Limit("raw Ward coordinate rank overflow".into()))?;
        let coordinate_count = scalar_count
            .checked_add(loops)
            .ok_or_else(|| Error::Limit("raw Ward coordinate count overflow".into()))?;
        if family.coordinates().len() != coordinate_count || shell.loop_index >= loops {
            return Err(Error::InvalidInput(
                "raw Ward coordinate layout mismatch".into(),
            ));
        }
        let energy = family.coordinates()[scalar_count + shell.loop_index].clone();
        let medium_symbols = family.coordinates()[scalar_count..]
            .iter()
            .map(|coordinate| {
                if let AtomView::Var(v) = coordinate.as_view() {
                    Ok(v.get_symbol())
                } else {
                    Err(Error::InvalidInput(
                        "raw Ward medium coordinate is not a symbol".into(),
                    ))
                }
            })
            .collect::<Result<Vec<_>>>()?;
        let diagonal_slot = (0..loops)
            .flat_map(|i| (i..loops).map(move |j| (i, j)))
            .position(|(i, j)| i == shell.loop_index && j == shell.loop_index)
            .ok_or_else(|| Error::InvalidInput("raw Ward shell Gram coordinate absent".into()))?;
        let measure = family.deformed_measure::<N>(eta, shifted)?;
        let factors = measure.factors();
        let roles = family.roles();
        if roles.len() != factors.len()
            || shell.physical_slot >= family.physical_slots()
            || family.physical_slots() > family.input_slots()
            || family.input_slots() > factors.len()
            || shell.upper_slot >= factors.len()
            || shell.lower_slot >= factors.len()
            || roles[shell.physical_slot] != IndexRole::RequiredCut
            || roles[shell.upper_slot] != IndexRole::Occupation
            || roles[shell.lower_slot] != IndexRole::Occupation
            || roles
                .iter()
                .filter(|r| **r == IndexRole::RequiredCut)
                .count()
                != 1
            || roles
                .iter()
                .filter(|r| **r == IndexRole::Occupation)
                .count()
                != 2
        {
            return Err(Error::InvalidInput(
                "raw Ward distribution role mismatch".into(),
            ));
        }
        let exact_zero = |a: Atom| a.expand().together().cancel().is_zero();
        if !exact_zero(&factors[shell.physical_slot] - &family.coordinates()[diagonal_slot])
            || !exact_zero(
                &factors[shell.upper_slot] - Atom::num(shell.chemical_potential.clone()) + &energy,
            )
            || !exact_zero(&factors[shell.lower_slot] - &energy)
        {
            return Err(Error::InvalidInput(
                "raw Ward shell or endpoint factor mismatch".into(),
            ));
        }
        let mut energy_dependent_slots = Vec::new();
        for (slot, factor) in factors.iter().enumerate() {
            let dependent = medium_symbols
                .iter()
                .any(|&s| !factor.derivative(s).together().cancel().is_zero());
            if slot < family.physical_slots() && dependent {
                return Err(Error::InvalidInput(
                    "raw Ward physical quadratic depends on the medium".into(),
                ));
            }
            if roles[slot] == IndexRole::Ordinary && dependent {
                energy_dependent_slots.push(slot);
            }
        }
        let identity = format!(
            "raw-polynomial-singleton-Ward-v1;C1;all-medium-dependent-base-indices-zero;lowerH0;upper0-or-positive;unrecentered;origin={};energy-dependent={:?};eta={}",
            self.source_identity(),
            energy_dependent_slots,
            Atom::var(eta).to_canonical_string(),
        );
        Ok(SingletonRawWardEvidence {
            identity,
            coordinates: family.coordinates().to_vec(),
            factors: factors.to_vec(),
            roles: roles.to_vec(),
            cut: shell.physical_slot,
            upper: shell.upper_slot,
            lower: shell.lower_slot,
            physical_slots: family.physical_slots(),
            input_slots: family.input_slots(),
            energy,
            energy_dependent_slots,
        })
    }

    /// Deliberately omit only the extra Ward theorem for other sealed classes.
    /// The family is still checked even when no Ward row is emitted.
    pub(crate) fn optional_singleton_raw_ward_evidence<const N: usize>(
        &self,
        family: &super::geometry::OccupiedCutFamily,
        shifted: &[usize],
        eta: Symbol,
    ) -> Result<Option<SingletonRawWardEvidence>> {
        self.validate_family(family, shifted)?;
        match &self.class {
            CertifiedEndpointClass::SingletonGerm { .. } => self
                .singleton_raw_ward_evidence::<N>(family, shifted, eta)
                .map(Some),
            CertifiedEndpointClass::RankOneBlocks { .. } => Ok(None),
        }
    }

    pub fn input_identity(&self) -> &str {
        &self.input_identity
    }
    pub fn source_options(&self) -> super::preparation::WeightedSourceOptions {
        self.source_options
    }

    pub fn cut_slots(&self) -> &[usize] {
        &self.cuts
    }
    pub fn shifted_slots(&self) -> &[usize] {
        &self.shifted
    }
    pub fn scope(&self) -> &'static str {
        match &self.class {
            CertifiedEndpointClass::SingletonGerm {
                virtual_loops: 1, ..
            } => "one-virtual-loop-null-cone",
            CertifiedEndpointClass::SingletonGerm { .. } => "singleton-UV-meromorphic-null-germ",
            CertifiedEndpointClass::RankOneBlocks {
                virtual_loops: 0, ..
            } => "two-compact-spacelike-transfer",
            CertifiedEndpointClass::RankOneBlocks { .. } => {
                "two-compact-independent-rank-one-blocks"
            }
        }
    }
    pub fn origin_identity(&self) -> &'static str {
        MASSLESS_FLOW_ORIGIN_IDENTITY
    }

    pub fn source_identity(&self) -> String {
        format!(
            "{};scope={};input={};cuts={:?};shifted={:?};options={:?};family={:?};proof={:?}",
            self.origin_identity(),
            self.scope(),
            self.input_identity,
            self.cuts,
            self.shifted,
            self.source_options,
            self.family_signature,
            self.class.structure()
        )
    }
    pub fn validate_family(
        &self,
        family: &super::geometry::OccupiedCutFamily,
        shifted: &[usize],
    ) -> Result<()> {
        if shifted != self.shifted || bound_family_signature(family) != self.family_signature {
            return Err(Error::InvalidInput(
                "massless flow evidence family or deformation mismatch".into(),
            ));
        }
        Ok(())
    }
    pub fn certified_origin_loops(
        &self,
        family: &super::geometry::OccupiedCutFamily,
        shifted: &[usize],
    ) -> Result<Vec<usize>> {
        self.validate_family(family, shifted)?;
        Ok(family.shells().iter().map(|s| s.loop_index).collect())
    }
    /// Enumerate complete maximal rank-deficient ordinary positive supports.
    /// This supplies proof data only; attaching its zero domains is explicitly
    /// opt-in at the guarded source factory. A budget failure supplies no zero.
    pub fn free_virtual_zero_supports(
        &self,
        family: &super::geometry::OccupiedCutFamily,
        shifted: &[usize],
        budget: usize,
    ) -> Result<Vec<FreeVirtualZeroSupport>> {
        self.validate_family(family, shifted)?;
        let h = self.class.virtual_loops();
        if h == 0 {
            return Ok(Vec::new());
        }
        let rows = self.class.rows();
        let mut flats = BTreeMap::new();
        if let CertifiedEndpointClass::RankOneBlocks {
            representatives, ..
        } = &self.class
        {
            if h > budget {
                return Err(Error::Limit(format!(
                    "free-virtual block flats need {h} domains; budget {budget}"
                )));
            }
            for omitted in 0..h {
                let selected = representatives
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i != omitted)
                    .map(|(_, r)| r.clone())
                    .collect();
                let flat = free_virtual_flat(rows, h, selected)?.ok_or_else(|| {
                    Error::InvalidInput("independent block hyperplane lost rank".into())
                })?;
                flats.insert(flat.allowed_positive_slots.clone(), flat);
            }
        } else {
            checked_flat_subsets(rows.len(), h - 1, budget)?;
            let k = h - 1;
            let mut selection = (0..k).collect::<Vec<_>>();
            loop {
                let selected = selection
                    .iter()
                    .map(|&i| rows[i].virtual_part.clone())
                    .collect();
                if let Some(flat) = free_virtual_flat(rows, h, selected)? {
                    flats.insert(flat.allowed_positive_slots.clone(), flat);
                    if flats.len() > budget {
                        return Err(Error::Limit(format!(
                            "free-virtual flats exceed domain budget {budget}"
                        )));
                    }
                }
                let Some(j) = (0..k).rev().find(|&j| selection[j] < rows.len() - k + j) else {
                    break;
                };
                selection[j] += 1;
                for l in j + 1..k {
                    selection[l] = selection[l - 1] + 1;
                }
            }
        }
        Ok(flats.into_values().collect())
    }

    /// Audit every basis, target and retained candidate label independently.
    /// A label's meromorphic proof domain can depend on its finite indices;
    /// taking the intersection for any finite requested set is harmless.
    pub fn validate_labels(
        &self,
        family: &super::geometry::OccupiedCutFamily,
        labels: &[crate::Integral],
    ) -> Result<MasslessLabelAudit> {
        self.validate_family(family, &self.shifted)?;
        let mut bounds = Vec::new();
        for label in labels {
            if label.0.len() != family.factors().len() {
                return Err(Error::InvalidInput(
                    "massless endpoint label arity mismatch".into(),
                ));
            }
            if (family.physical_slots()..family.input_slots()).any(|s| label.0[s] > 0)
                || family
                    .shells()
                    .iter()
                    .any(|s| label.0[s.upper_slot] < 0 || label.0[s.lower_slot] < 0)
            {
                return Err(Error::Unsupported("massless endpoint labels require polynomial completions and nonnegative occupation indices".into()));
            }
            let p = self.shifted.iter().try_fold(0_i32, |sum, &s| {
                sum.checked_add(i32::from(label.0[s].max(0)))
                    .ok_or_else(|| {
                        Error::Limit("massless endpoint positive power sum overflow".into())
                    })
            })?;
            let active_virtual_rank = rank(
                self.class
                    .rows()
                    .iter()
                    .filter(|r| label.0[r.slot] > 0)
                    .map(|r| r.virtual_part.clone())
                    .collect(),
            );
            let mut block_powers = vec![0_i32; self.class.virtual_loops()];
            let mut transfer_power = 0_i32;
            if let CertifiedEndpointClass::RankOneBlocks {
                rows, assignments, ..
            } = &self.class
            {
                for (row, block) in rows.iter().zip(assignments) {
                    let n = i32::from(label.0[row.slot].max(0));
                    let sum = match block {
                        Some(b) => &mut block_powers[*b],
                        None => &mut transfer_power,
                    };
                    *sum = sum.checked_add(n).ok_or_else(|| {
                        Error::Limit("massless block positive powers overflow".into())
                    })?;
                }
            }
            let jets = family
                .shells()
                .iter()
                .map(|s| {
                    i32::from(label.0[s.physical_slot].saturating_sub(1).max(0))
                        + i32::from(label.0[s.upper_slot].saturating_sub(1).max(0))
                })
                .collect::<Vec<_>>();
            let missing_cut = family
                .shells()
                .iter()
                .any(|s| label.0[s.physical_slot] <= 0);
            let lower = family.shells().iter().any(|s| label.0[s.lower_slot] > 0);
            let mut inequalities = Vec::new();
            let classification = if missing_cut {
                "vanishing-required-cut"
            } else if lower {
                for shell in family.shells().iter().filter(|s| label.0[s.lower_slot] > 0) {
                    inequalities.push(EndpointInequality {
                        stratum: format!(
                            "lower-origin zero jets at compact loop {}",
                            shell.loop_index
                        ),
                        degree: EndpointDegree {
                            twice_dimension_coefficient: 2,
                            constant: -2 * i32::from(label.0[shell.physical_slot])
                                - i32::from(label.0[shell.lower_slot])
                                - 2 * i32::from(label.0[shell.upper_slot].saturating_sub(1).max(0)),
                        },
                    });
                }
                "joint-dimensional-lower-contact-zero"
            } else {
                if active_virtual_rank < self.class.virtual_loops() {
                    "scaleless-unrestricted-virtual-polynomial"
                } else {
                    match &self.class {
                        CertifiedEndpointClass::SingletonGerm { virtual_loops, .. } => {
                            let j = jets[0];
                            inequalities.push(EndpointInequality {
                                stratum: "UV-meromorphic virtual null-germ endpoint jets".into(),
                                degree: EndpointDegree {
                                    twice_dimension_coefficient: i32::try_from(*virtual_loops)
                                        .unwrap(),
                                    constant: -p - j,
                                },
                            });
                            "uniform-null-cone-endpoint"
                        }
                        CertifiedEndpointClass::RankOneBlocks { virtual_loops, .. } => {
                            let total = jets.iter().sum::<i32>();
                            for (block, &power) in block_powers.iter().enumerate() {
                                inequalities.push(EndpointInequality {
                                    stratum: format!(
                                        "Gaussian virtual block {block} after finite jets"
                                    ),
                                    degree: EndpointDegree {
                                        twice_dimension_coefficient: 1,
                                        constant: -power - total,
                                    },
                                });
                            }
                            for (i, &j) in jets.iter().enumerate() {
                                inequalities.push(EndpointInequality {
                                    stratum: format!(
                                        "compact radius {i} after finite transfer jets"
                                    ),
                                    degree: EndpointDegree {
                                        twice_dimension_coefficient: 2,
                                        constant: -2 - 2 * j - transfer_power - total,
                                    },
                                });
                            }
                            inequalities.push(EndpointInequality {
                                stratum: "compact collinear transfer after finite jets".into(),
                                degree: EndpointDegree {
                                    twice_dimension_coefficient: 1,
                                    constant: -1 - transfer_power - total,
                                },
                            });
                            if *virtual_loops == 0 {
                                "pure-compact-transfer-endpoint"
                            } else {
                                "independent-rank-one-block-endpoint"
                            }
                        }
                    }
                }
            };
            // Even a factorized zero must have a concrete common high-D
            // domain for every remaining compact factor. This is particularly
            // relevant for P=0 virtual polynomials and a lower contact on only
            // one of two compact loops.
            if !missing_cut {
                for (i, &j) in jets.iter().enumerate() {
                    let degree = EndpointDegree {
                        twice_dimension_coefficient: 2,
                        constant: -2 - 2 * j,
                    };
                    if !inequalities.iter().any(|bound| bound.degree == degree) {
                        inequalities.push(EndpointInequality {
                            stratum: format!("common compact origin domain at loop {i}"),
                            degree,
                        });
                    }
                }
            }
            let witness = regular_dimension_witness(&inequalities, self.class.virtual_loops())?;
            bounds.push(MasslessLabelBound {
                indices: label.0.clone(),
                classification,
                uncut_positive_power_sum: p,
                active_virtual_rank,
                virtual_block_positive_power_sums: block_powers,
                pure_transfer_positive_power_sum: transfer_power,
                mass_and_upper_jet_orders: jets,
                inequalities,
                high_dimension_witness: witness.to_string(),
            });
        }
        Ok(MasslessLabelAudit {
            endpoint_structure: self.class.structure(),
            origin_identity: self.origin_identity(),
            labels: bounds,
        })
    }
}

#[cfg(test)]
mod raw_ward_draft_tests {
    use super::*;
    use crate::finite_density::{
        DensityInput,
        guarded::{GuardedMeasureIdentity, IndexRole},
        preparation::WeightedSourceOptions,
    };

    fn definition() -> DensityInput {
        let mut input: DensityInput = serde_json::from_str(include_str!(
            "../../examples/finite_density/massive_two_loop_sunset.json"
        ))
        .unwrap();
        for edge in &mut input.edges {
            edge.mass_squared = "0".into();
        }
        input
    }

    fn data(
        input: DensityInput,
        cuts: &[usize],
    ) -> (
        PreparedDensityInput,
        super::super::geometry::OccupiedCutFamily,
        Vec<usize>,
        MasslessFlowEvidence,
    ) {
        let input = input.prepare().unwrap();
        let family = input.occupied_cut(cuts, 16).unwrap().at_physical_masses();
        let shifted = (0..family.physical_slots())
            .filter(|s| family.roles()[*s] == IndexRole::Ordinary)
            .collect::<Vec<_>>();
        let proof =
            MasslessFlowEvidence::new(&input, &family, &shifted, WeightedSourceOptions::default())
                .unwrap();
        (input, family, shifted, proof)
    }

    fn identity() -> GuardedMeasureIdentity {
        GuardedMeasureIdentity {
            measure: "raw Ward draft: actual physical source factory".into(),
            support: "positive mu; joint origin proof supplied separately".into(),
            orientation: "future shell".into(),
            normalization: "source identity only".into(),
            branch: "fixed-T regulator order".into(),
            deformation: "all uncut D-eta".into(),
        }
    }

    #[test]
    fn raw_ward_maps_actual_energy_and_checks_every_base_guard() {
        let eta = symbol!("raw_ward_draft::eta");
        let d = Atom::var(symbol!("raw_ward_draft::D"));
        for mapping in 0..3 {
            let mut input = definition();
            if mapping == 1 {
                // Nonunit routing gives physical completion E_a/2. Charges
                // remain exactly conserved; no hand-edited proof object.
                input.edges[0].routing = vec!["2".into(), "0".into()];
                input.edges[1].routing = vec!["0".into(), "2".into()];
                input.edges[2].routing = vec!["2".into(), "-2".into()];
                input.edges[0].charges = vec![2];
                input.edges[1].charges = vec![2];
            } else if mapping == 2 {
                input.edges[0].routing = vec!["1".into(), "1".into()];
                input.edges[1].routing = vec!["0".into(), "1".into()];
                input.edges[2].routing = vec!["1".into(), "0".into()];
                input.loop_charges = vec![vec![0], vec![1]];
            }
            let (_, family, shifted, origin) = data(input, &[0]);
            let certificate = origin
                .singleton_raw_ward_evidence::<12>(&family, &shifted, eta)
                .unwrap();
            let measure = family.deformed_measure::<12>(eta, &shifted).unwrap();
            let indices =
                std::array::from_fn(|i| symbol!(format!("raw_ward_draft::a_{mapping}_{i}")));
            let rows = measure
                .raw_singleton_ward_sources(&certificate, &d, &indices, 2)
                .unwrap();
            assert_eq!(rows.len(), 2);
            assert!(rows.iter().all(|r| r.nonzero_conditions.is_empty()));
            for (branch, row) in rows.iter().enumerate() {
                let mut base = [0_i64; 12];
                base[certificate.cut] = 1;
                base[certificate.upper] = if branch == 0 { 0 } else { 3 };
                assert!(row.domain.contains(&base));
                for value in [0, 2, -1] {
                    let mut bad = base;
                    bad[certificate.cut] = value;
                    assert!(!row.domain.contains(&bad));
                }
                for &slot in &certificate.energy_dependent_slots {
                    for value in [-1, 1] {
                        let mut bad = base;
                        bad[slot] = value;
                        assert!(!row.domain.contains(&bad));
                    }
                }
                for slot in [certificate.lower, family.factors().len()] {
                    let mut bad = base;
                    bad[slot] = 1;
                    assert!(!row.domain.contains(&bad));
                }
                let mut bad = base;
                bad[certificate.upper] = -1;
                assert!(!row.domain.contains(&bad));
                for &slot in &shifted {
                    for n in [-3, 0, 4] {
                        let mut free = base;
                        free[slot] = n;
                        assert!(row.domain.contains(&free));
                    }
                }
                // Reconstruct the coefficient of H_(s+1) from actual factors.
                // This verifies cE and mixed-energy conversion without guessing slots.
                let mut reconstructed = Atom::zero();
                for term in &row.terms {
                    if term.shift[certificate.upper] == 0 {
                        assert!(term.shift.iter().all(|n| *n == 0));
                        assert!(
                            (&term.coefficient - &d + Atom::num(2))
                                .together()
                                .cancel()
                                .is_zero()
                        );
                        continue;
                    }
                    assert_eq!(term.shift[certificate.upper], 1);
                    let mut value = term.coefficient.clone();
                    for (slot, shift) in term.shift.iter().enumerate() {
                        if slot == certificate.upper {
                            continue;
                        }
                        assert!(*shift <= 0);
                        if *shift < 0 {
                            value *= family.factors()[slot].pow(Atom::num(-i32::from(*shift)));
                        }
                    }
                    reconstructed += value;
                }
                let expected = certificate.energy()
                    * if branch == 0 {
                        Atom::num(-1)
                    } else {
                        Atom::var(indices[certificate.upper])
                    };
                assert!(
                    (reconstructed - expected)
                        .expand()
                        .together()
                        .cancel()
                        .is_zero()
                );
            }
            assert!(
                measure
                    .raw_singleton_ward_sources(&certificate, &d, &indices, 1)
                    .is_err()
            );
            assert!(
                origin
                    .singleton_raw_ward_evidence::<12>(&family, &[], eta)
                    .is_err()
            );
            let changed = family
                .deformed_measure::<12>(symbol!("raw_ward_draft::different_eta"), &shifted)
                .unwrap();
            assert!(
                changed
                    .raw_singleton_ward_sources(&certificate, &d, &indices, 2)
                    .is_err()
            );
            let mut duplicate = indices;
            duplicate[1] = duplicate[0];
            assert!(
                measure
                    .raw_singleton_ward_sources(&certificate, &d, &duplicate, 2)
                    .is_err()
            );
        }
    }

    #[test]
    fn raw_ward_only_accepts_typed_singleton_and_preserves_factory_epoch() {
        let (input, family, shifted, origin) = data(definition(), &[0]);
        let epsilon = symbol!("raw_ward_draft::epsilon");
        let eta = symbol!("raw_ward_draft::eta_factory");
        let old_opts = WeightedSourceOptions::default();
        let opts = WeightedSourceOptions {
            policy: super::super::preparation::WeightedSourcePolicy::PolynomialClosure,
            ..old_opts
        };
        let ward_origin = MasslessFlowEvidence::new(&input, &family, &shifted, opts).unwrap();
        let ordinary = family
            .guarded_sources_with_massless_origin::<12>(
                epsilon,
                4,
                eta,
                &shifted,
                16,
                vec![],
                identity(),
                old_opts,
                &origin,
            )
            .unwrap();
        let ward = family
            .guarded_sources_with_massless_origin::<12>(
                epsilon,
                4,
                eta,
                &shifted,
                16,
                vec![],
                identity(),
                opts,
                &ward_origin,
            )
            .unwrap();
        assert_ne!(
            ordinary.context.sources().measure_id(),
            ward.context.sources().measure_id()
        );
        assert!(
            ward.context
                .sources()
                .measure_id()
                .contains("raw-polynomial-singleton-Ward-v1")
        );
        let mut label = [0; 12];
        label[family.shells()[0].physical_slot] = 1;
        assert!(!ward.context.sources().is_zero(&label));
        label[family.shells()[0].lower_slot] = 1;
        assert!(ward.context.sources().is_zero(&label));
        label[family.physical_slots()] = 1;
        assert!(!ward.context.sources().is_zero(&label));
        assert!(ward.deformation.derivative(label).is_err());
        let encoded = ordinary
            .context
            .discover(vec![], [], Default::default())
            .unwrap()
            .program
            .encode(Default::default())
            .unwrap();
        assert!(ward.context.decode(&encoded, Default::default()).is_err());
        let encoded = ward
            .context
            .discover(vec![], [], Default::default())
            .unwrap()
            .program
            .encode(Default::default())
            .unwrap();
        assert!(ward.context.decode(&encoded, Default::default()).is_ok());
        let different = WeightedSourceOptions {
            free_virtual_zero_sectors: true,
            ..opts
        };
        assert!(
            family
                .guarded_sources_with_massless_origin::<12>(
                    epsilon,
                    4,
                    eta,
                    &shifted,
                    16,
                    vec![],
                    identity(),
                    different,
                    &origin
                )
                .is_err()
        );
        let (_, double, shifted_double, proof_double) = data(definition(), &[0, 1]);
        assert!(
            proof_double
                .singleton_raw_ward_evidence::<12>(&double, &shifted_double, eta)
                .is_err()
        );
        assert!(
            origin
                .singleton_raw_ward_evidence::<12>(&double, &shifted_double, eta)
                .is_err()
        );
        let mut massive = definition();
        massive.edges[0].mass_squared = "1/4".into();
        let massive = massive.prepare().unwrap();
        let mass_family = massive.occupied_cut(&[0], 16).unwrap().at_physical_masses();
        assert!(MasslessFlowEvidence::new(&massive, &mass_family, &shifted, opts).is_err());
        let mut zero_mu = definition();
        zero_mu.chemical_potentials[0] = "0".into();
        let zero_mu = zero_mu.prepare().unwrap();
        if let Ok(zero_family) = zero_mu.occupied_cut(&[0], 16) {
            assert!(
                MasslessFlowEvidence::new(
                    &zero_mu,
                    &zero_family.at_physical_masses(),
                    &shifted,
                    opts
                )
                .is_err()
            );
        }
        assert!(
            MasslessFlowEvidence::new(
                &input,
                &family,
                &shifted,
                WeightedSourceOptions {
                    positive_compact_energy_powers: true,
                    ..opts
                }
            )
            .is_err()
        );
    }

    #[test]
    fn raw_ward_native_rule_replays_in_full_bound_source_context() {
        use crate::finite_density::guarded::{GuardedContext, IndexBounds, IndexDomain};
        let (_, family, shifted, origin) = data(definition(), &[0]);
        let eta = symbol!("raw_ward_draft::native_eta");
        let d = symbol!("raw_ward_draft::native_D");
        let measure = family.deformed_measure::<12>(eta, &shifted).unwrap();
        let proof = origin
            .singleton_raw_ward_evidence::<12>(&family, &shifted, eta)
            .unwrap();
        let indices = std::array::from_fn(|i| symbol!(format!("raw_ward_draft::native_a_{i}")));
        let rows = measure
            .raw_singleton_ward_sources(&proof, &Atom::var(d), &indices, 2)
            .unwrap();
        let mut description = identity();
        description.measure.push_str(proof.identity());
        let context = GuardedContext::new_with_physical_arity(
            description,
            *measure.roles(),
            indices,
            vec![d, eta],
            rows,
            measure.physical_arity(),
        )
        .unwrap();
        let energy_slot = family
            .factors()
            .iter()
            .position(|f| f == proof.energy())
            .unwrap();
        let mut target = [0_i64; 12];
        target[proof.cut] = 1;
        for slot in shifted {
            target[slot] = 1;
        }
        target[energy_slot] = -1;
        target[proof.upper] = 1;
        let mut bulk = target;
        bulk[energy_slot] = 0;
        bulk[proof.upper] = 0;
        let point = IndexDomain::new(target.map(IndexBounds::fixed)).unwrap();
        let program = context
            .discover(vec![point], [bulk], Default::default())
            .unwrap()
            .program;
        let restored = context
            .decode(
                &program.encode(Default::default()).unwrap(),
                Default::default(),
            )
            .unwrap();
        let reduction = restored.reduce(target, Default::default()).unwrap();
        assert!(reduction.unresolved.is_empty());
        assert_eq!(reduction.terms.len(), 1);
        assert!(
            (reduction.terms.get(&bulk).unwrap() - Atom::var(d) + Atom::num(2))
                .together()
                .cancel()
                .is_zero()
        );
        let mut raised = target;
        raised[proof.cut] = 2;
        let rejected = restored.reduce(raised, Default::default()).unwrap();
        assert!(!rejected.unresolved.is_empty());
        let mut extra_energy = target;
        extra_energy[energy_slot] = -2;
        assert!(
            !restored
                .reduce(extra_energy, Default::default())
                .unwrap()
                .unresolved
                .is_empty()
        );
    }

    #[test]
    fn raw_ward_upper_distribution_signs_match_direct_radial_moments() {
        // C1 shell measure after E integration is proportional to E^(D-3).
        // H_s=(-1)^(s-1) delta^(s-1)/(s-1)! for s>=1. This is an independent
        // exact evaluation of their action, not the source conversion algorithm.
        fn power(mu: &Rational, p: i64) -> Rational {
            if p < 0 {
                return Rational::one() / power(mu, -p);
            }
            (0..p).fold(Rational::one(), |x, _| x * mu)
        }
        fn moment(p: i64, s: usize, mu: &Rational) -> Rational {
            if s == 0 {
                return power(mu, p + 1) / Rational::from(p + 1);
            }
            let derivative = (0..s - 1).fold(Rational::one(), |x, k| {
                x * Rational::from(p - k as i64) / Rational::from((k + 1) as i64)
            });
            let sign = if s % 2 == 1 { 1 } else { -1 };
            Rational::from(sign) * derivative * power(mu, p - s as i64 + 1)
        }
        for d in [8_i64, 11, 14] {
            for mu in [Rational::from((2, 3)), Rational::from(2)] {
                assert_eq!(
                    Rational::from(d - 2) * moment(d - 3, 0, &mu) - moment(d - 2, 1, &mu),
                    Rational::zero()
                );
                for s in 1..=4 {
                    assert_eq!(
                        Rational::from(d - 2) * moment(d - 3, s, &mu)
                            + Rational::from(s as i64) * moment(d - 2, s + 1, &mu),
                        Rational::zero()
                    );
                    assert!(!moment(d - 2, s + 1, &mu).is_zero());
                }
            }
        }
    }
}

#[cfg(test)]
mod permit_tests {
    use super::*;
    use crate::finite_density::{DensityInput, preparation::WeightedSourceOptions};
    fn small() -> PreparedDensityInput {
        let mut definition = serde_json::from_str::<DensityInput>(include_str!(
            "../../examples/finite_density/massive_two_loop_sunset.json"
        ))
        .unwrap();
        for edge in &mut definition.edges {
            edge.mass_squared = "0".into();
        }
        definition.prepare().unwrap()
    }
    #[test]
    fn sealed_permit_audits_generated_zero_and_surface_labels() {
        let input = small();
        let family = input.occupied_cut(&[0], 16).unwrap().at_physical_masses();
        let proof =
            MasslessFlowEvidence::new(&input, &family, &[1, 2], WeightedSourceOptions::default())
                .unwrap();
        let mut label = crate::Integral(vec![0; family.factors().len()]);
        let shell = &family.shells()[0];
        label.0[shell.physical_slot] = 2;
        label.0[1] = 1;
        label.0[2] = 2;
        label.0[shell.upper_slot] = 3;
        let mut zero = label.clone();
        zero.0[shell.physical_slot] = 0;
        let mut polynomial = label.clone();
        polynomial.0[1] = -1;
        polynomial.0[2] = 0;
        let mut lower = label.clone();
        lower.0[shell.lower_slot] = 2;
        let audit = proof
            .validate_labels(&family, &[label.clone(), zero, polynomial, lower])
            .unwrap();
        assert_eq!(audit.labels[0].mass_and_upper_jet_orders, [3]);
        assert_eq!(audit.labels[1].classification, "vanishing-required-cut");
        assert_eq!(
            audit.labels[2].classification,
            "scaleless-unrestricted-virtual-polynomial"
        );
        assert_eq!(audit.labels[2].high_dimension_witness, "25/3");

        assert_eq!(
            audit.labels[3].classification,
            "joint-dimensional-lower-contact-zero"
        );
        let mut inverse = label;
        inverse.0[family.physical_slots()] = 1;
        assert!(proof.validate_labels(&family, &[inverse]).is_err());
        assert!(proof.validate_family(&family, &[1]).is_err());
        assert_eq!(proof.certified_origin_loops(&family, &[1, 2]).unwrap(), [0]);
    }
    #[test]
    fn e7_sealed_germ_and_block_permits_audit_generated_rank_loss() {
        let input = serde_json::from_str::<DensityInput>(include_str!(
            "../../examples/finite_density/chain_of_three_parallel_pairs.json"
        ))
        .unwrap()
        .prepare()
        .unwrap();
        for cuts in [vec![0], vec![4], vec![0, 4]] {
            let family = input.occupied_cut(&cuts, 16).unwrap().at_physical_masses();
            let shifted = (0..7).filter(|s| !cuts.contains(s)).collect::<Vec<_>>();
            let proof = MasslessFlowEvidence::new(
                &input,
                &family,
                &shifted,
                WeightedSourceOptions::default(),
            )
            .unwrap();
            let labels = family
                .targets()
                .iter()
                .flat_map(|t| t.keys().cloned())
                .collect::<Vec<_>>();
            let audit = proof.validate_labels(&family, &labels).unwrap();
            assert_eq!(audit.endpoint_structure.virtual_loops, 4 - cuts.len());
            for bound in &audit.labels {
                let d = rational(&parse(&bound.high_dimension_witness).unwrap()).unwrap();
                assert!(bound.inequalities.iter().all(|i| i.degree.at(&d) > 0));
                if bound.active_virtual_rank < 4 - cuts.len() {
                    assert_eq!(
                        bound.classification,
                        "scaleless-unrestricted-virtual-polynomial"
                    );
                }
            }
            let mut free = labels[0].clone();
            for &slot in &shifted {
                free.0[slot] = 0;
            }
            let zero = proof.validate_labels(&family, &[free.clone()]).unwrap();
            assert_eq!(
                zero.labels[0].classification,
                "scaleless-unrestricted-virtual-polynomial"
            );
            free.0[family.physical_slots()] = 1;
            assert!(proof.validate_labels(&family, &[free]).is_err());
            if cuts.len() == 2 {
                assert_eq!(proof.scope(), "two-compact-independent-rank-one-blocks");
                assert_eq!(audit.endpoint_structure.block_representatives.len(), 2);
                assert!(
                    audit
                        .endpoint_structure
                        .block_assignments
                        .iter()
                        .any(|(_, b)| b.is_none())
                );
            }
        }
    }
    #[test]
    fn rank_one_permit_rejects_genuinely_coupled_quadratics_and_prism() {
        let synthetic = [vec![1, -1, 1, 0], vec![1, -1, 0, 1], vec![1, -1, 1, 1]]
            .into_iter()
            .map(|r| r.into_iter().map(Atom::num).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        assert!(matches!(rank_one_blocks(&synthetic,&[0,1,2],2),
            Err(Error::Unsupported(ref s)) if s.contains("independent rank-one blocks")));
        let input = serde_json::from_str::<DensityInput>(include_str!(
            "../../examples/finite_density/triangular_prism.json"
        ))
        .unwrap()
        .prepare()
        .unwrap();
        let family = input
            .occupied_cut(&[1, 5], 32)
            .unwrap()
            .at_physical_masses();
        assert!(
            matches!(MasslessFlowEvidence::new(&input,&family,&[0,2,3,4,6,7,8],WeightedSourceOptions::default()),
            Err(Error::Unsupported(ref s)) if s.contains("single spacelike direction"))
        );
    }
    #[test]
    fn singleton_external_only_factors_remain_explicitly_unadmitted() {
        let mut definition = small().input().clone();
        definition.edges[1].routing = vec!["1".into(), "0".into()];
        definition.edges[2].vertices = [0, 0];
        definition.edges[2].routing = vec!["0".into(), "1".into()];
        definition.loop_charges[1][0] = 0;
        let input = definition.prepare().unwrap();
        let family = input.occupied_cut(&[0], 16).unwrap().at_physical_masses();
        assert!(
            matches!(MasslessFlowEvidence::new(&input,&family,&[1,2],WeightedSourceOptions::default()),
            Err(Error::Unsupported(ref s)) if s.contains("external-only"))
        );
    }
    #[test]
    fn sealed_germ_and_blocks_survive_routing_shear_and_orientation() {
        let mut definition = serde_json::from_str::<DensityInput>(include_str!(
            "../../examples/finite_density/chain_of_three_parallel_pairs.json"
        ))
        .unwrap();
        for slot in [0, 1] {
            definition.edges[slot].vertices.swap(0, 1);
            for c in &mut definition.edges[slot].routing {
                *c = (-parse(c).unwrap()).to_string();
            }
            for q in &mut definition.edges[slot].charges {
                *q = -*q;
            }
        }
        for edge in &mut definition.edges {
            edge.routing[0] =
                (&parse(&edge.routing[0]).unwrap() + &parse(&edge.routing[1]).unwrap()).to_string();
        }
        definition.loop_charges[1][0] = -1;
        let rules = BTreeMap::from([
            (parse("g1_2").unwrap(), parse("g1_2+g1_1").unwrap()),
            (parse("g2_4").unwrap(), parse("g2_4+g1_4").unwrap()),
        ]);
        for target in &mut definition.targets {
            target.numerator = substitute(&parse(&target.numerator).unwrap(), &rules).to_string();
        }
        definition.name = "renamed_and_sheared".into();
        let input = definition.prepare().unwrap();
        for cuts in [vec![0], vec![4], vec![0, 4]] {
            let family = input.occupied_cut(&cuts, 16).unwrap().at_physical_masses();
            let shifted = (0..7).filter(|s| !cuts.contains(s)).collect::<Vec<_>>();
            let proof = MasslessFlowEvidence::new(
                &input,
                &family,
                &shifted,
                WeightedSourceOptions::default(),
            )
            .unwrap();
            let mut label = family.targets()[0].keys().next().unwrap().clone();
            for row in proof.class.rows() {
                if !row.virtual_part.last().unwrap().is_zero() {
                    label.0[row.slot] = 0;
                }
            }
            let audit = proof.validate_labels(&family, &[label.clone()]).unwrap();
            assert!(audit.labels[0].active_virtual_rank < proof.class.virtual_loops());
            assert_eq!(
                audit.labels[0].classification,
                "scaleless-unrestricted-virtual-polynomial"
            );
            label.0[family.physical_slots()] = 1;
            assert!(proof.validate_labels(&family, &[label]).is_err());
            assert!(
                proof
                    .source_identity()
                    .contains(MASSLESS_ENDPOINT_PROOF_VERSION)
            );
            if cuts.len() == 2 {
                assert!(
                    proof
                        .class
                        .structure()
                        .virtual_affine_shift
                        .iter()
                        .any(|x| x != "0")
                );
            }
        }
    }
    #[test]
    fn generic_witness_avoids_all_integer_shifted_loop_poles() {
        let inequalities = [EndpointInequality {
            stratum: "test bound".into(),
            degree: EndpointDegree {
                twice_dimension_coefficient: 1,
                constant: -17,
            },
        }];
        for h in 1..=16 {
            let d = regular_dimension_witness(&inequalities, h).unwrap();
            assert!(inequalities[0].degree.at(&d) > 0);
            for k in 1..=h {
                let x = &d * &Rational::from((k as i64, 2));
                assert!(!x.is_integer());
            }
        }
    }
    #[test]
    fn free_virtual_flat_boxes_equal_all_rank_deficient_e7_supports() {
        let input = serde_json::from_str::<DensityInput>(include_str!(
            "../../examples/finite_density/chain_of_three_parallel_pairs.json"
        ))
        .unwrap()
        .prepare()
        .unwrap();
        for (cuts, count) in [(vec![0], 6), (vec![4], 6), (vec![0, 4], 2)] {
            let family = input.occupied_cut(&cuts, 16).unwrap().at_physical_masses();
            let shifted = (0..7).filter(|s| !cuts.contains(s)).collect::<Vec<_>>();
            let proof = MasslessFlowEvidence::new(
                &input,
                &family,
                &shifted,
                WeightedSourceOptions::default(),
            )
            .unwrap();
            let flats = proof
                .free_virtual_zero_supports(&family, &shifted, 16)
                .unwrap();
            assert_eq!(flats.len(), count);
            for mask in 0..1usize << shifted.len() {
                let active = shifted
                    .iter()
                    .enumerate()
                    .filter(|(j, _)| mask & (1 << j) != 0)
                    .map(|(_, s)| *s)
                    .collect::<BTreeSet<_>>();
                let r = rank(
                    proof
                        .class
                        .rows()
                        .iter()
                        .filter(|r| active.contains(&r.slot))
                        .map(|r| r.virtual_part.clone())
                        .collect(),
                );
                let covered = flats.iter().any(|f| {
                    f.forced_nonpositive_slots
                        .iter()
                        .all(|s| !active.contains(s))
                });
                assert_eq!(covered, r < proof.class.virtual_loops());
                for flat in &flats {
                    assert!(
                        flat.forced_nonpositive_slots
                            .iter()
                            .all(|s| !cuts.contains(s))
                    );
                    assert_eq!(
                        flat.virtual_null_direction.len(),
                        proof.class.virtual_loops()
                    );
                }
            }
            assert!(matches!(
                proof.free_virtual_zero_supports(&family, &shifted, 1),
                Err(Error::Limit(_))
            ));
            assert!(proof.free_virtual_zero_supports(&family, &[], 16).is_err());
        }
        assert!(matches!(
            checked_flat_subsets(usize::MAX, 8, usize::MAX),
            Err(Error::Limit(_))
        ));
    }
}
