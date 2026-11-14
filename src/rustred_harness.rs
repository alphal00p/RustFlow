//! Reader for RustRed harness solutions (`rustred-harness-solution.v1`).
//!
//! A harness solution is a reduction table at a FIXED rational dimension `d0`
//! of the eta-derivative targets of an eta-deformed family onto its masters:
//! `rules[k] = {target, terms: [{integral, coefficient}]}`, coefficients exact
//! rational functions of `eta` only. This module turns it into the input of
//! [`PreparedFlow::from_supplied_at_epsilon`]:
//!
//! * the family (from the harness family JSON: `denominators[i] = {coefficients
//!   over loop-loop coordinates, constant}`, Minkowski `q^2 - m^2`), with the
//!   deformed slot's `-eta` removed from its constant (RustFlow adds `D_i - eta`
//!   through the deformation mask; the reader checks that the harness constant
//!   is exactly `constant(0) - eta`, i.e. the same sign convention);
//! * the closed connection `d basis_i / d eta = sum_j A_ij basis_j`, built from
//!   `d/d eta I[a] = a_s I[a_s + 1]` (deformed slot `s`) and the rules;
//! * target maps from the rules (or identity rows when a target is a master).
//!
//! Masters with `a_s = 0` are eta-independent: their rows are zero, but they
//! appear as columns of the other rows and stay in the basis.
use crate::engine::{PreparedFlow, SuppliedAuxiliarySystem};
use crate::reduction::{LinearCombination, ReducedSystem};
use crate::*;
use serde::Deserialize;
use std::collections::BTreeMap;
use symbolica::prelude::*;

pub const SOLUTION_SCHEMA: &str = "rustred-harness-solution.v1";
/// Namespace in which harness coefficient strings are parsed (`eta` becomes `rustred_harness::eta`).
pub const NAMESPACE: &str = "rustred_harness";

pub fn eta_symbol() -> Symbol {
    symbol!("rustred_harness::eta")
}

#[derive(Clone, Debug, Deserialize)]
pub struct HarnessCoordinate {
    pub kind: String,
    pub left: usize,
    pub right: usize,
}

#[derive(Clone, Debug, Deserialize)]
pub struct HarnessDenominator {
    pub coefficients: Vec<String>,
    pub constant: String,
}

/// The harness family description (e.g. `classes/16359/family-eta-s0.json`).
#[derive(Clone, Debug, Deserialize)]
pub struct HarnessFamily {
    pub name: String,
    pub loop_momenta: Vec<String>,
    #[serde(default)]
    pub external_momenta: Vec<String>,
    pub coordinates: Vec<HarnessCoordinate>,
    pub denominators: Vec<HarnessDenominator>,
    pub parameters: Vec<String>,
    #[serde(default)]
    pub denominator_convention: Option<String>,
    #[serde(default)]
    pub sha256: Option<String>,
    /// Must all be zero: the reader has no power shifts.
    #[serde(default)]
    pub power_shifts: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Deserialize)]
pub struct HarnessTerm {
    pub integral: Vec<i16>,
    pub coefficient: String,
}

#[derive(Clone, Debug, Deserialize)]
pub struct HarnessRule {
    pub target: Vec<i16>,
    pub terms: Vec<HarnessTerm>,
}

#[derive(Clone, Debug, Deserialize)]
pub struct HarnessSolutionFamily {
    pub arity: usize,
    pub name: String,
    #[serde(default)]
    pub fingerprint: Option<String>,
}

#[derive(Clone, Debug, Deserialize)]
pub struct HarnessSolution {
    pub schema: String,
    pub family: HarnessSolutionFamily,
    pub parameters: Vec<String>,
    pub fixed: BTreeMap<String, String>,
    pub masters: Vec<Vec<i16>>,
    pub rules: Vec<HarnessRule>,
    #[serde(default)]
    pub provenance: serde_json::Value,
}

fn parse(text: &str) -> Result<Atom> {
    Atom::parse(text, NAMESPACE, Default::default())
        .map_err(|e| Error::InvalidInput(format!("harness expression {text:?}: {e}")))
}

fn rational_of(atom: &Atom) -> Result<Rational> {
    let atom = atom.clone().expand();
    if let AtomView::Num(n) = atom.as_view() {
        if let symbolica::coefficient::Coefficient::Complex(c) = n.get_coeff_view().to_owned() {
            if c.im.is_zero() {
                return Ok(c.re);
            }
        }
    }
    Err(Error::InvalidInput(format!(
        "not an exact rational number: {atom}"
    )))
}

impl HarnessSolution {
    pub fn from_slice(bytes: &[u8]) -> Result<Self> {
        let solution: Self = serde_json::from_slice(bytes)
            .map_err(|e| Error::InvalidInput(format!("harness solution JSON: {e}")))?;
        if solution.schema != SOLUTION_SCHEMA {
            return Err(Error::InvalidInput(format!(
                "unsupported harness schema {}",
                solution.schema
            )));
        }
        Ok(solution)
    }

    /// The fixed rational dimension `d0`.
    pub fn fixed_dimension(&self) -> Result<Rational> {
        let d = self
            .fixed
            .get("d")
            .ok_or_else(|| Error::InvalidInput("harness solution has no fixed d".into()))?;
        rational_of(&parse(d)?)
    }

    /// `epsilon = (dimension - d0) / 2` for an integer base dimension.
    pub fn epsilon_sample(&self, dimension: i64) -> Result<Rational> {
        let d0 = self.fixed_dimension()?;
        rational_of(&((Atom::num(dimension) - Atom::num(d0)) / Atom::num(2)))
    }
}

/// Undeformed RustFlow family plus deformation mask from a harness family.
/// `deformed` is the slot carrying `-eta` in its constant.
pub fn family_from_harness(
    harness: &HarnessFamily,
    epsilon: Symbol,
    dimension: i64,
) -> Result<(IntegralFamily, Vec<bool>)> {
    if !harness.external_momenta.is_empty() {
        return Err(Error::Unsupported(
            "harness reader: only vacuum (no external momenta) families".into(),
        ));
    }
    let loops = harness.loop_momenta.len();
    let mut expected = Vec::new();
    for i in 0..loops {
        for j in i..loops {
            expected.push((i, j));
        }
    }
    if harness.coordinates.len() != expected.len()
        || harness
            .coordinates
            .iter()
            .zip(&expected)
            .any(|(c, &(i, j))| c.kind != "loop_loop" || c.left != i || c.right != j)
    {
        return Err(Error::InvalidInput(
            "harness coordinates differ from RustFlow's loop-loop order (i <= j, row major)".into(),
        ));
    }
    if let Some(convention) = &harness.denominator_convention {
        if !convention.starts_with("Minkowski q^2-") {
            return Err(Error::Unsupported(format!(
                "harness denominator convention {convention:?} is not Minkowski q^2 - m^2"
            )));
        }
    }
    for (slot, shift) in harness.power_shifts.iter().enumerate() {
        if !rational_of(&parse(shift)?)?.is_zero() {
            return Err(Error::Unsupported(format!(
                "harness denominator {slot}: power shift {shift} (only unshifted families)"
            )));
        }
    }
    let eta = Atom::var(eta_symbol());
    let mut propagators = Vec::new();
    let mut mask = Vec::new();
    for (slot, denominator) in harness.denominators.iter().enumerate() {
        if denominator.coefficients.len() != expected.len() {
            return Err(Error::InvalidInput(format!(
                "harness denominator {slot}: coefficient count"
            )));
        }
        let scalar_products = denominator
            .coefficients
            .iter()
            .map(|c| parse(c))
            .collect::<Result<Vec<_>>>()?;
        let constant = parse(&denominator.constant)?.expand();
        let undeformed = constant.replace(eta.clone()).with(Atom::num(0)).expand();
        let slope = constant.derivative(eta_symbol()).expand();
        let shifted = if slope.is_zero() {
            false
        } else if (slope.clone() + Atom::num(1)).expand().is_zero() {
            true
        } else {
            return Err(Error::InvalidInput(format!(
                "harness denominator {slot}: constant {constant} is not c - eta (RustFlow deforms D_i - eta)"
            )));
        };
        rational_of(&undeformed)?;
        for value in &scalar_products {
            rational_of(value)?;
        }
        mask.push(shifted);
        propagators.push(Propagator {
            constant: undeformed,
            scalar_products,
        });
    }
    let family = IntegralFamily {
        name: harness.name.clone(),
        loops: harness.loop_momenta.clone(),
        external: vec![],
        external_gram: vec![],
        physical_propagators: propagators.len(),
        propagators,
        epsilon,
        dimension,
    };
    Ok((family, mask))
}

/// Exact connection and target maps of a harness solution.
pub fn reduced_system(
    solution: &HarnessSolution,
    mask: &[bool],
    targets: &[Integral],
) -> Result<ReducedSystem> {
    reduced_system_with(solution, mask, targets, false)
}

/// As [`reduced_system`]; with `homogeneous_eblock` the basis keeps only masters
/// with a positive index on a deformed slot and every eta-independent (e-free)
/// term is dropped. That is the homogeneous system satisfied by the
/// discontinuity across the eta cut (in-house dispersive route, prep5.py /
/// local_solutions.load_eblock); it is NOT the full AMF connection.
pub fn reduced_system_with(
    solution: &HarnessSolution,
    mask: &[bool],
    targets: &[Integral],
    homogeneous_eblock: bool,
) -> Result<ReducedSystem> {
    let deformed = |m: &[i16]| m.iter().zip(mask).any(|(&a, &q)| q && a > 0);
    if solution.family.arity != mask.len() {
        return Err(Error::InvalidInput(
            "harness arity differs from family".into(),
        ));
    }
    if solution.parameters.iter().any(|p| p != "d" && p != "eta") {
        return Err(Error::Unsupported(format!(
            "harness parameters {:?}: only d (fixed) and eta",
            solution.parameters
        )));
    }
    let basis = solution
        .masters
        .iter()
        .filter(|m| !homogeneous_eblock || deformed(m))
        .map(|m| Integral(m.clone()))
        .collect::<Vec<_>>();
    let all_masters = solution
        .masters
        .iter()
        .map(|m| Integral(m.clone()))
        .collect::<std::collections::BTreeSet<_>>();
    let index = basis
        .iter()
        .enumerate()
        .map(|(i, m)| (m.clone(), i))
        .collect::<BTreeMap<_, _>>();
    if index.len() != basis.len() {
        return Err(Error::InvalidInput("duplicate harness master".into()));
    }
    // Harness solutions may repeat a rule (seen for masters); repeated
    // targets must carry the same terms.
    let mut rules = BTreeMap::new();
    for rule in &solution.rules {
        if let Some(previous) = rules.insert(Integral(rule.target.clone()), &rule.terms) {
            if previous != &rule.terms {
                return Err(Error::InvalidInput(format!(
                    "harness has conflicting rules for {:?}",
                    rule.target
                )));
            }
        }
    }
    let map_of = |integral: &Integral| -> Result<LinearCombination> {
        if index.contains_key(integral) {
            return Ok(LinearCombination::from([(integral.clone(), Atom::one())]));
        }
        if homogeneous_eblock && all_masters.contains(integral) {
            return Ok(LinearCombination::new());
        }
        let terms = rules.get(integral).ok_or_else(|| {
            Error::IncompleteReduction(format!("harness has no rule for {:?}", integral.0))
        })?;
        let mut row = LinearCombination::new();
        for term in terms.iter() {
            let master = Integral(term.integral.clone());
            if homogeneous_eblock && !deformed(&term.integral) && all_masters.contains(&master) {
                continue;
            }
            if !index.contains_key(&master) {
                return Err(Error::IncompleteReduction(format!(
                    "harness rule term {:?} is not a master",
                    term.integral
                )));
            }
            let c = parse(&term.coefficient)?;
            if !c.is_zero() {
                let entry = row.entry(master).or_insert_with(Atom::new);
                *entry = (entry.clone() + c).together().cancel();
            }
        }
        row.retain(|_, c| !c.is_zero());
        Ok(row)
    };
    let n = basis.len();
    let mut matrix = vec![vec![Atom::new(); n]; n];
    for (i, master) in basis.iter().enumerate() {
        for (target, factor) in crate::family::eta_derivative(master, mask)? {
            for (column, coefficient) in map_of(&target)? {
                let j = index[&column];
                matrix[i][j] = (matrix[i][j].clone() + factor.clone() * coefficient)
                    .together()
                    .cancel();
            }
        }
    }
    let target_maps = targets.iter().map(map_of).collect::<Result<Vec<_>>>()?;
    let candidates = targets
        .iter()
        .cloned()
        .zip(target_maps.iter().cloned())
        .collect();
    Ok(ReducedSystem {
        basis,
        matrix,
        targets: target_maps,
        nonzero_conditions: vec![],
        candidates,
        transformations: vec![],
    })
}

/// Complete supplied declaration; `provenance_extra` should carry the file hash.
pub fn supplied_system(
    solution: &HarnessSolution,
    mask: &[bool],
    targets: &[Integral],
    provenance_extra: &str,
) -> Result<SuppliedAuxiliarySystem> {
    supplied_system_with(solution, mask, targets, provenance_extra, false)
}

pub fn supplied_system_with(
    solution: &HarnessSolution,
    mask: &[bool],
    targets: &[Integral],
    provenance_extra: &str,
    homogeneous_eblock: bool,
) -> Result<SuppliedAuxiliarySystem> {
    let reduced = reduced_system_with(solution, mask, targets, homogeneous_eblock)?;
    Ok(SuppliedAuxiliarySystem {
        variable: eta_symbol(),
        reduced,
        deformation_mask: mask.to_vec(),
        provenance: format!(
            "{SOLUTION_SCHEMA}{}; family {} ({}); fixed {:?}; {}; {}",
            if homogeneous_eblock {
                " (homogeneous e-block: e-free terms dropped)"
            } else {
                ""
            },
            solution.family.name,
            solution
                .family
                .fingerprint
                .as_deref()
                .map(|f| blake3::hash(f.as_bytes()).to_hex().to_string())
                .unwrap_or_default(),
            solution.fixed,
            solution.provenance,
            provenance_extra
        ),
    })
}

/// Load family + solution and prepare the flow at the solution's epsilon sample.
#[allow(clippy::too_many_arguments)]
pub fn prepare(
    harness_family: &HarnessFamily,
    solution: &HarnessSolution,
    targets: &[Integral],
    epsilon: Symbol,
    options: &FlowOptions,
    context: &RunContext,
    provenance_extra: &str,
    homogeneous_eblock: bool,
) -> Result<PreparedFlow> {
    let (family, mask) = family_from_harness(harness_family, epsilon, options.dimension)?;
    let sample = solution.epsilon_sample(options.dimension)?;
    let supplied = supplied_system_with(
        solution,
        &mask,
        targets,
        provenance_extra,
        homogeneous_eblock,
    )?;
    PreparedFlow::from_supplied_at_epsilon(
        &family,
        targets,
        &KinematicPoint::default(),
        supplied,
        options,
        context,
        &sample,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    // One-loop tadpole I[a] = int 1/(q^2 - 1 - eta)^a at fixed d0 = 4 - 2 eps:
    // I[2] = (d0/2 - 1)/(1 + eta) I[1]; with eps = 1/10, d0 = 19/5: I[2] = (9/10)/(1+eta) I[1].
    const FAMILY: &str = r#"{"name":"tad","loop_momenta":["k1"],"external_momenta":[],
        "coordinates":[{"kind":"loop_loop","left":0,"right":0}],
        "denominators":[{"coefficients":["1"],"constant":"-1-eta"}],
        "parameters":["d","eta"],"denominator_convention":"Minkowski q^2-mUV^2; stored unit mass squared"}"#;
    const SOLUTION: &str = r#"{"schema":"rustred-harness-solution.v1",
        "family":{"arity":1,"name":"tad","fingerprint":"x"},"parameters":["d","eta"],
        "fixed":{"d":"19/5"},"coefficient_field":{"kind":"Q"},"masters":[[1]],
        "rules":[{"target":[2],"terms":[{"integral":[1],"coefficient":"(9/10)/(1+eta)"}]},
                 {"target":[3],"terms":[{"integral":[1],"coefficient":"(9/10)*(-1/10)/(2*(1+eta)^2)"}]}],
        "stats":{},"provenance":{"lane":"test"}}"#;

    #[test]
    fn synthetic_tadpole_reads_and_prepares() -> Result<()> {
        let harness: HarnessFamily = serde_json::from_str(FAMILY).unwrap();
        let solution = HarnessSolution::from_slice(SOLUTION.as_bytes())?;
        assert_eq!(solution.epsilon_sample(4)?, Rational::from((1, 10)));
        let epsilon = symbol!("rustred_harness_test::eps");
        let (family, mask) = family_from_harness(&harness, epsilon, 4)?;
        assert_eq!(mask, vec![true]);
        assert_eq!(family.propagators[0].constant, Atom::num(-1));
        let targets = vec![Integral(vec![1]), Integral(vec![3])];
        let reduced = reduced_system(&solution, &mask, &targets)?;
        let expected = parse("(9/10)/(1+eta)")?;
        assert!(
            (reduced.matrix[0][0].clone() - expected)
                .together()
                .cancel()
                .is_zero()
        );
        assert_eq!(reduced.targets[0].len(), 1);
        let options = FlowOptions::default();
        let flow = prepare(
            &harness,
            &solution,
            &targets,
            epsilon,
            &options,
            &RunContext::default(),
            "test",
            false,
        )?;
        assert_eq!(flow.epsilon_sample(), Some(&Rational::from((1, 10))));
        assert_eq!(flow.deformation_mask(), &[true]);
        assert_eq!(flow.basis().len(), 1);
        Ok(())
    }

    #[test]
    fn homogeneous_eblock_drops_eta_free_columns() -> Result<()> {
        let solution = HarnessSolution::from_slice(
            SOLUTION
                .replace(r#""masters":[[1]]"#, r#""masters":[[1],[0]]"#)
                .replace(
                    r#""coefficient":"(9/10)/(1+eta)"}]"#,
                    r#""coefficient":"(9/10)/(1+eta)"},{"integral":[0],"coefficient":"1/eta"}]"#,
                )
                .as_bytes(),
        )?;
        let targets = vec![Integral(vec![1])];
        let full = reduced_system(&solution, &[true], &targets)?;
        assert_eq!(full.basis.len(), 2);
        assert!(!full.matrix[0][1].is_zero() && full.matrix[1][1].is_zero());
        let homogeneous = reduced_system_with(&solution, &[true], &targets, true)?;
        assert_eq!(homogeneous.basis, vec![Integral(vec![1])]);
        assert_eq!(homogeneous.matrix.len(), 1);
        Ok(())
    }

    #[test]
    fn power_shifts_and_conflicting_rules_are_refused() -> Result<()> {
        let epsilon = symbol!("rustred_harness_test::eps4");
        let shifted: HarnessFamily = serde_json::from_str(
            &FAMILY.replace(r#""parameters":"#, r#""power_shifts":["1"],"parameters":"#),
        )
        .unwrap();
        assert!(family_from_harness(&shifted, epsilon, 4).is_err());
        let unshifted: HarnessFamily = serde_json::from_str(
            &FAMILY.replace(r#""parameters":"#, r#""power_shifts":["0"],"parameters":"#),
        )
        .unwrap();
        assert!(family_from_harness(&unshifted, epsilon, 4).is_ok());
        let rule = r#"{"target":[2],"terms":[{"integral":[1],"coefficient":"(9/10)/(1+eta)"}]}"#;
        let repeated = SOLUTION.replace(r#""rules":["#, &format!(r#""rules":[{rule},"#));
        let targets = vec![Integral(vec![1])];
        reduced_system(
            &HarnessSolution::from_slice(repeated.as_bytes())?,
            &[true],
            &targets,
        )?;
        let conflicting = SOLUTION.replace(
            r#""rules":["#,
            &format!(r#""rules":[{},"#, rule.replace("(9/10)", "(8/10)")),
        );
        assert!(
            reduced_system(
                &HarnessSolution::from_slice(conflicting.as_bytes())?,
                &[true],
                &targets
            )
            .is_err()
        );
        Ok(())
    }

    #[test]
    fn wrong_deformation_sign_is_refused() {
        let harness: HarnessFamily =
            serde_json::from_str(&FAMILY.replace("-1-eta", "-1+eta")).unwrap();
        assert!(family_from_harness(&harness, symbol!("rustred_harness_test::eps2"), 4).is_err());
    }

    #[test]
    fn fixed_sample_rejects_epsilon_dependent_connection() -> Result<()> {
        let harness: HarnessFamily = serde_json::from_str(FAMILY).unwrap();
        let solution = HarnessSolution::from_slice(SOLUTION.as_bytes())?;
        let epsilon = symbol!("rustred_harness_test::eps3");
        let (family, mask) = family_from_harness(&harness, epsilon, 4)?;
        let targets = vec![Integral(vec![1])];
        let mut supplied = supplied_system(&solution, &mask, &targets, "test")?;
        supplied.reduced.matrix[0][0] =
            Atom::var(epsilon) / (Atom::num(1) + Atom::var(eta_symbol()));
        assert!(
            PreparedFlow::from_supplied_at_epsilon(
                &family,
                &targets,
                &KinematicPoint::default(),
                supplied,
                &FlowOptions::default(),
                &RunContext::default(),
                &Rational::from((1, 10)),
            )
            .is_err()
        );
        Ok(())
    }
}
