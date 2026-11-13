//! The homogeneous e-block of the eta-DE at a fixed rational d0: dD/dz = A(z) D, z = 1 + eta.
//!
//! Masters with a positive index on the deformed slot `e` form the e-block; d/deta I[a] = a_e I[a + 1_e]
//! is reduced with the RustRed rules and every term on an e-free master is dropped (the discontinuity
//! across the eta cut obeys the homogeneous system). Entries are combined EXACTLY per (i, j) as
//! n1/d1 + n2/d2 = (n1 d2 + n2 d1)/(d1 d2) (oracle: prep5.py `combine`, local_solutions.load_eblock;
//! same e-block rule as `rustred_harness::reduced_system_with(.., homogeneous_eblock = true)`).
use super::poly::{self, QPoly, ZPoly};
use rug::Rational;
use std::collections::{BTreeMap, HashMap};

/// one entry of A(z) (as a rational function of eta): A_ij = num(eta) / den(eta)
#[derive(Clone, Debug)]
pub struct Term {
    pub i: usize,
    pub j: usize,
    pub num: ZPoly,
    pub den: ZPoly,
}

#[derive(Clone, Debug)]
pub struct EBlockSystem {
    pub d0: Rational,
    pub e_slot: usize,
    pub masters: Vec<Vec<i64>>,
    pub terms: Vec<Term>,
}

/// minimal view of a RustRed harness solution (the RustFlow reader fills it, see harness.rs; tests build it by hand)
pub struct RuleSet {
    pub d0: Rational,
    pub masters: Vec<Vec<i64>>,
    /// target -> [(integral, coefficient text "(num)/(den)" in eta)]
    pub rules: HashMap<Vec<i64>, Vec<(Vec<i64>, String)>>,
}

impl RuleSet {
    /// parse a `rustred-harness-solution.v1` JSON directly (standalone tools; RustFlow uses the reader)
    pub fn from_json(v: &serde_json::Value) -> Result<RuleSet, String> {
        if v["schema"] != "rustred-harness-solution.v1" {
            return Err(format!("unsupported schema {}", v["schema"]));
        }
        let d0 = v["fixed"]["d"].as_str().ok_or("no fixed d")?.parse::<Rational>().map_err(|e| e.to_string())?;
        let vec_i = |x: &serde_json::Value| -> Vec<i64> { x.as_array().unwrap().iter().map(|a| a.as_i64().unwrap()).collect() };
        let masters = v["masters"].as_array().ok_or("no masters")?.iter().map(vec_i).collect();
        let mut rules = HashMap::new();
        for r in v["rules"].as_array().ok_or("no rules")? {
            let terms = r["terms"]
                .as_array()
                .unwrap()
                .iter()
                .map(|t| (vec_i(&t["integral"]), t["coefficient"].as_str().unwrap().to_string()))
                .collect();
            rules.insert(vec_i(&r["target"]), terms);
        }
        Ok(RuleSet { d0, masters, rules })
    }
}

// ---------------- coefficient parser: sums of c*eta^k, optionally "(P)/(Q)" ----------------
fn split_top_div(s: &str) -> (&str, &str) {
    let b = s.as_bytes();
    let mut depth = 0i32;
    for (i, &ch) in b.iter().enumerate() {
        match ch {
            b'(' => depth += 1,
            b')' => depth -= 1,
            b'/' if depth == 0 && i > 0 && b[i - 1] == b')' => return (&s[..i], &s[i + 1..]),
            _ => {}
        }
    }
    (s, "1")
}

fn strip_parens(mut s: &str) -> &str {
    loop {
        let t = s.trim();
        if t.starts_with('(') && t.ends_with(')') {
            // only if the outer parentheses match each other
            let b = t.as_bytes();
            let mut depth = 0;
            let mut outer = true;
            for (i, &ch) in b.iter().enumerate() {
                if ch == b'(' {
                    depth += 1;
                } else if ch == b')' {
                    depth -= 1;
                    if depth == 0 && i + 1 < b.len() {
                        outer = false;
                        break;
                    }
                }
            }
            if outer {
                s = &t[1..t.len() - 1];
                continue;
            }
        }
        return t;
    }
}

/// a sum of monomials [+-] [a[/b]] [*] [var[^k]] -> ascending coefficients
pub fn parse_sum(text: &str, var: &str) -> Result<QPoly, String> {
    let s: String = strip_parens(text).chars().filter(|c| !c.is_whitespace()).collect();
    let b = s.as_bytes();
    let mut out: BTreeMap<usize, Rational> = BTreeMap::new();
    let mut i = 0;
    if b.is_empty() {
        return Err(format!("empty polynomial in {text:?}"));
    }
    while i < b.len() {
        let mut sign = 1;
        while i < b.len() && (b[i] == b'+' || b[i] == b'-') {
            if b[i] == b'-' {
                sign = -sign;
            }
            i += 1;
        }
        // number
        let st = i;
        while i < b.len() && (b[i].is_ascii_digit() || b[i] == b'/') {
            i += 1;
        }
        let mut coef = if i > st {
            s[st..i].parse::<Rational>().map_err(|e| format!("coefficient {:?} in {text:?}: {e}", &s[st..i]))?
        } else {
            Rational::from(1)
        };
        if i < b.len() && b[i] == b'*' {
            i += 1;
        }
        let mut k = 0usize;
        if s[i..].starts_with(var) {
            i += var.len();
            k = 1;
            if i < b.len() && b[i] == b'^' {
                i += 1;
                let st2 = i;
                while i < b.len() && b[i].is_ascii_digit() {
                    i += 1;
                }
                k = s[st2..i].parse().map_err(|e| format!("exponent in {text:?}: {e}"))?;
            }
        } else if i == st {
            return Err(format!("cannot parse monomial at {:?} in {text:?}", &s[st..]));
        }
        if i < b.len() && b[i] == b'*' {
            // trailing "*eta" after "var^k" is not expected
            return Err(format!("unexpected '*' in {text:?}"));
        }
        if sign < 0 {
            coef = -coef;
        }
        *out.entry(k).or_insert_with(Rational::new) += coef;
    }
    let n = out.keys().max().copied().unwrap_or(0);
    let mut p = vec![Rational::new(); n + 1];
    for (k, c) in out {
        p[k] = c;
    }
    poly::trim(&mut p);
    Ok(p)
}

/// "(num)/(den)" -> (num, den) exact
pub fn parse_coefficient(text: &str, var: &str) -> Result<(QPoly, QPoly), String> {
    let (n, d) = split_top_div(text.trim());
    Ok((parse_sum(n, var)?, parse_sum(d, var)?))
}

impl EBlockSystem {
    pub fn from_rules(rs: &RuleSet, e_slot: usize) -> Result<EBlockSystem, String> {
        let all: HashMap<&Vec<i64>, usize> = rs.masters.iter().enumerate().map(|(i, m)| (m, i)).collect();
        let emast: Vec<Vec<i64>> = rs.masters.iter().filter(|m| m[e_slot] > 0).cloned().collect();
        let eidx: HashMap<&Vec<i64>, usize> = emast.iter().enumerate().map(|(i, m)| (m, i)).collect();
        let mut entries: BTreeMap<(usize, usize), Vec<(QPoly, QPoly)>> = BTreeMap::new();
        for m in &emast {
            let a = m[e_slot];
            let mut r = m.clone();
            r[e_slot] += 1;
            let im = eidx[m];
            if all.contains_key(&r) {
                let jr = *eidx.get(&r).ok_or("raised master is not in the e-block")?;
                entries.entry((im, jr)).or_default().push((vec![Rational::from(a)], vec![Rational::from(1)]));
                continue;
            }
            let terms = rs.rules.get(&r).ok_or_else(|| format!("missing derivative rule for {m:?}"))?;
            for (ti, c) in terms {
                if ti[e_slot] <= 0 {
                    continue; // e-free: no discontinuity
                }
                let j = *eidx.get(ti).ok_or_else(|| format!("rule term {ti:?} is not an e-block master"))?;
                let (n, d) = parse_coefficient(c, "eta")?;
                entries.entry((im, j)).or_default().push((poly::scale(&n, &Rational::from(a)), d));
            }
        }
        let mut terms = Vec::new();
        for ((i, j), lst) in entries {
            let mut n: QPoly = vec![Rational::new()];
            let mut d: QPoly = vec![Rational::from(1)];
            for (nn, dd) in lst {
                n = poly::add(&poly::mul(&n, &dd), &poly::mul(&nn, &d));
                d = poly::mul(&d, &dd);
            }
            if poly::is_zero(&n) {
                continue;
            }
            terms.push(Term { i, j, num: ZPoly::from_q(&n), den: ZPoly::from_q(&d) });
        }
        Ok(EBlockSystem { d0: rs.d0.clone(), e_slot, masters: emast, terms })
    }

    pub fn dim(&self) -> usize {
        self.masters.len()
    }

    pub fn epsilon(&self) -> Rational {
        (Rational::from(4) - &self.d0) / 2
    }

    /// distinct denominators (as primitive integer polynomials with a positive leading coefficient)
    pub fn distinct_denominators(&self) -> Vec<Vec<rug::Integer>> {
        let mut set = std::collections::BTreeSet::new();
        for t in &self.terms {
            set.insert(t.den.clone().prim_signed());
        }
        set.into_iter().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn parses_harness_coefficients() {
        let (n, d) = parse_coefficient("(-1/150000)/(1+eta)", "eta").unwrap();
        assert_eq!(n, vec![Rational::from((-1, 150000))]);
        assert_eq!(d, vec![Rational::from(1), Rational::from(1)]);
        let (n, d) = parse_coefficient("(2100232/165-17/3*eta+eta^2-eta^4)", "eta").unwrap();
        assert_eq!(d, vec![Rational::from(1)]);
        assert_eq!(n.len(), 5);
        assert_eq!(n[1], Rational::from((-17, 3)));
        assert_eq!(n[2], Rational::from(1));
        assert_eq!(n[4], Rational::from(-1));
        assert_eq!(parse_coefficient("1", "eta").unwrap().0, vec![Rational::from(1)]);
        assert_eq!(parse_coefficient("(eta-3/2*eta^2)/(2-eta)", "eta").unwrap().1, vec![Rational::from(2), Rational::from(-1)]);
    }
}
