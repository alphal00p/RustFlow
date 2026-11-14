//! RustFlow adapter: the dispersive e-block from the RustRed harness reader (`rustred_harness::HarnessSolution`).
//! Coefficients in the RustRed "(num)/(den)" expanded form go through the fast exact parser; any other form
//! falls back to Symbolica (same namespace as the reader) and its exact rational-polynomial conversion.
use super::poly::QPoly;
use super::system::{parse_coefficient, EBlockSystem, RuleSet};
use crate::rustred_harness::{eta_symbol, HarnessSolution, NAMESPACE};
use std::collections::HashMap;
use symbolica::prelude::*;

fn poly_coefficients(p: &MultivariatePolynomial<IntegerRing, u16>) -> Result<QPoly, String> {
    let vars = p.variables();
    if vars.len() > 1 || (vars.len() == 1 && vars[0] != PolyVariable::Symbol(eta_symbol())) {
        return Err(format!("harness coefficient depends on {vars:?}, expected eta only"));
    }
    let mut out: Vec<rug::Rational> = Vec::new();
    for (c, ex) in p.coefficients.iter().zip(p.exponents_iter()) {
        let k = ex.first().copied().unwrap_or(0) as usize;
        if out.len() <= k {
            out.resize(k + 1, rug::Rational::new());
        }
        out[k] = c.to_string().parse::<rug::Rational>().map_err(|e| e.to_string())?;
    }
    if out.is_empty() {
        out.push(rug::Rational::new());
    }
    Ok(out)
}

/// exact (num, den) of a harness coefficient
pub fn parse_harness_coefficient(text: &str) -> Result<(QPoly, QPoly), String> {
    if let Ok(x) = parse_coefficient(text, "eta") {
        return Ok(x);
    }
    let atom = Atom::parse(text, NAMESPACE, Default::default()).map_err(|e| format!("harness expression {text:?}: {e}"))?;
    let rp: RationalPolynomial<IntegerRing, u16> = atom.try_to_rational_polynomial(&Q, &Z, None).map_err(|e| format!("{text:?}: {e}"))?;
    Ok((poly_coefficients(&rp.numerator)?, poly_coefficients(&rp.denominator)?))
}

/// the reader's solution as a rule set (fixed d0, masters, rules)
pub fn rule_set(sol: &HarnessSolution) -> crate::Result<RuleSet> {
    let d = sol.fixed.get("d").ok_or_else(|| crate::Error::InvalidInput("harness solution has no fixed d".into()))?;
    let d0 = d.trim().parse::<rug::Rational>().map_err(|e| crate::Error::InvalidInput(format!("fixed d {d:?}: {e}")))?;
    let masters = sol.masters.iter().map(|m| m.iter().map(|&x| x as i64).collect()).collect();
    let mut rules = HashMap::new();
    for r in &sol.rules {
        rules.insert(
            r.target.iter().map(|&x| x as i64).collect(),
            r.terms.iter().map(|t| (t.integral.iter().map(|&x| x as i64).collect(), t.coefficient.clone())).collect(),
        );
    }
    Ok(RuleSet { d0, masters, rules })
}

/// homogeneous e-block of a harness solution with deformed slot `e_slot` (same e-block rule as
/// `rustred_harness::reduced_system_with(.., homogeneous_eblock = true)`)
pub fn eblock(sol: &HarnessSolution, e_slot: usize) -> crate::Result<EBlockSystem> {
    let rs = rule_set(sol)?;
    EBlockSystem::from_rules_with(&rs, e_slot, &parse_harness_coefficient).map_err(crate::Error::InvalidInput)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rustred_harness::reduced_system_with;

    // two e-block masters [1, 0], [2, 0] (slot 0 deformed) and one e-free master [0, 1]
    const SOLUTION: &str = r#"{"schema":"rustred-harness-solution.v1",
        "family":{"arity":2,"name":"toy","fingerprint":"x"},"parameters":["d","eta"],
        "fixed":{"d":"19/5"},"coefficient_field":{"kind":"Q"},"masters":[[1,0],[2,0],[0,1]],
        "rules":[{"target":[3,0],"terms":[{"integral":[1,0],"coefficient":"(9/10)*(-1/10)/(2*(1+eta)^2)"},
                                          {"integral":[2,0],"coefficient":"(3-eta^2)/(1+eta)"},
                                          {"integral":[0,1],"coefficient":"(7)/(eta)"}]}],
        "stats":{},"provenance":{"lane":"test"}}"#;

    #[test]
    fn eblock_agrees_with_the_reader() {
        let sol = HarnessSolution::from_slice(SOLUTION.as_bytes()).unwrap();
        let sys = eblock(&sol, 0).unwrap();
        assert_eq!(sys.masters, vec![vec![1, 0], vec![2, 0]]);
        let reduced = reduced_system_with(&sol, &[true, false], &[], true).unwrap();
        assert_eq!(reduced.basis.len(), 2);
        // compare every entry at eta = 5/7 exactly
        let x = rug::Rational::from((5, 7));
        let ev = |p: &[rug::Rational]| p.iter().rev().fold(rug::Rational::new(), |acc, a| acc * &x + a);
        let eta = Atom::var(eta_symbol());
        for i in 0..2 {
            for j in 0..2 {
                let mine = sys.terms.iter().find(|t| t.i == i && t.j == j).map(|t| ev(&t.num.to_q()) / ev(&t.den.to_q())).unwrap_or_default();
                let theirs = reduced.matrix[i][j].replace(eta.clone()).with(Atom::num(Rational::from((5, 7)))).expand();
                let theirs = if theirs.is_zero() { rug::Rational::new() } else { theirs.to_string().parse::<rug::Rational>().unwrap() };
                assert_eq!(mine, theirs, "entry ({i}, {j})");
            }
        }
    }
}
