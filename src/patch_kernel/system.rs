//! Exact rational-function systems for the patch kernel.
//!
//! `dy_i/dx = sum_j (N_ij(x) / D_ij(x)) y_j` with integer polynomials N, D.
//! Entries are kept in canonical (row, column) order and every distinct
//! polynomial is stored once, so the operation order of the kernel (and
//! therefore its rounding) depends only on the system, never on how it was
//! loaded.
use crate::{Error, Result};
use rug::Integer;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use symbolica::prelude::*;

pub const SYSTEM_SCHEMA: &str = "rustflow.patch-kernel.system.v1";

/// One nonzero matrix entry `scale * N / D` (indices into
/// [`RationalSystem::polys`]; N and D primitive, D with positive leading
/// coefficient, so proportional denominators are shared).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Entry {
    pub row: usize,
    pub col: usize,
    pub num: usize,
    pub den: usize,
    pub scale: rug::Rational,
}

#[derive(Clone, Debug)]
pub struct RationalSystem {
    pub dim: usize,
    /// Name of the variable `x` (information only).
    pub variable: String,
    /// Distinct integer polynomials, coefficients from degree 0 upwards,
    /// without trailing zeros (the zero polynomial is never stored).
    pub polys: Vec<Vec<Integer>>,
    /// Nonzero entries sorted by (row, col).
    pub entries: Vec<Entry>,
}

#[derive(Serialize, Deserialize)]
struct SystemFile {
    schema: String,
    dim: usize,
    variable: String,
    polys: Vec<Vec<String>>,
    entries: Vec<(usize, usize, usize, usize, String)>,
}

fn trim(p: &mut Vec<Integer>) {
    while p.len() > 1 && p.last().is_some_and(|c| *c == 0) {
        p.pop();
    }
}

fn parse_rational(text: &str) -> Result<rug::Rational> {
    text.trim()
        .parse::<rug::Rational>()
        .map_err(|e| Error::InvalidInput(format!("rational {text:?}: {e}")))
}

fn json_rational(v: &serde_json::Value) -> Result<rug::Rational> {
    match v {
        serde_json::Value::String(s) => parse_rational(s),
        serde_json::Value::Number(n) => n
            .as_i64()
            .map(rug::Rational::from)
            .ok_or_else(|| Error::InvalidInput(format!("non-integer number {n} in system"))),
        _ => Err(Error::InvalidInput(format!("bad rational {v}"))),
    }
}

/// Integer polynomials (N, D) with N/D = num/den for rational coefficient lists.
fn clear_denominators(
    num: &[rug::Rational],
    den: &[rug::Rational],
) -> (Vec<Integer>, Vec<Integer>) {
    let lcm = |p: &[rug::Rational]| p.iter().fold(Integer::from(1), |acc, c| acc.lcm(c.denom()));
    let (ln, ld) = (lcm(num), lcm(den));
    // N/D = (ln num / ln) / (ld den / ld) = (ld * (ln num)) / (ln * (ld den))
    let scale = |p: &[rug::Rational], l: &Integer, other: &Integer| {
        p.iter()
            .map(|c| (c.numer() * Integer::from(l / c.denom())) * other)
            .collect::<Vec<_>>()
    };
    (scale(num, &ln, &ld), scale(den, &ld, &ln))
}

impl RationalSystem {
    /// Build from raw entries `N/D`; normalises to `scale * N' / D'` with
    /// primitive N', D' (positive leading coefficient of D'), sorts and
    /// deduplicates polynomials.
    pub fn from_entries(
        dim: usize,
        variable: &str,
        raw: Vec<(usize, usize, Vec<Integer>, Vec<Integer>)>,
    ) -> Result<Self> {
        let mut raw = raw
            .into_iter()
            .filter_map(|(i, j, mut n, mut d)| {
                trim(&mut n);
                trim(&mut d);
                (!n.iter().all(|c| *c == 0)).then_some((i, j, n, d))
            })
            .collect::<Vec<_>>();
        raw.sort_by_key(|e| (e.0, e.1));
        let mut index: BTreeMap<Vec<Integer>, usize> = BTreeMap::new();
        let mut polys: Vec<Vec<Integer>> = Vec::new();
        let mut entries = Vec::with_capacity(raw.len());
        let mut last = None;
        for (i, j, mut n, mut d) in raw {
            if i >= dim || j >= dim {
                return Err(Error::InvalidInput(format!(
                    "entry ({i}, {j}) outside dimension {dim}"
                )));
            }
            if last == Some((i, j)) {
                return Err(Error::InvalidInput(format!("duplicate entry ({i}, {j})")));
            }
            last = Some((i, j));
            if d.iter().all(|c| *c == 0) {
                return Err(Error::InvalidInput(format!(
                    "zero denominator at ({i}, {j})"
                )));
            }
            let content = |p: &[Integer]| p.iter().fold(Integer::new(), |g, c| g.gcd(c));
            let mut cn = content(&n);
            let mut cd = content(&d);
            if d.last().is_some_and(|c| *c < 0) {
                cd = -cd;
            }
            if n.last().is_some_and(|c| *c < 0) {
                cn = -cn;
            }
            for c in n.iter_mut() {
                c.div_exact_mut(&cn);
            }
            for c in d.iter_mut() {
                c.div_exact_mut(&cd);
            }
            let scale = rug::Rational::from((cn, cd));
            let mut id = |p: Vec<Integer>| -> usize {
                if let Some(&k) = index.get(&p) {
                    return k;
                }
                polys.push(p.clone());
                index.insert(p, polys.len() - 1);
                polys.len() - 1
            };
            let num = id(n);
            let den = id(d);
            entries.push(Entry {
                row: i,
                col: j,
                num,
                den,
                scale,
            });
        }
        Ok(Self {
            dim,
            variable: variable.to_string(),
            polys,
            entries,
        })
    }

    /// From a square matrix of Symbolica rational functions in `variable`
    /// (for example the RustRed harness e-block:
    /// `rustred_harness::reduced_system_with(.., true).matrix`).
    pub fn from_matrix_atoms(matrix: &[Vec<Atom>], variable: Symbol) -> Result<Self> {
        let dim = matrix.len();
        if matrix.iter().any(|r| r.len() != dim) {
            return Err(Error::InvalidInput("square matrix required".into()));
        }
        let mut raw = Vec::new();
        for (i, row) in matrix.iter().enumerate() {
            for (j, a) in row.iter().enumerate() {
                if a.is_zero() {
                    continue;
                }
                let r: symbolica::domains::rational_polynomial::RationalPolynomial<
                    symbolica::domains::integer::IntegerRing,
                    u16,
                > = a
                    .try_to_rational_polynomial(&Q, &Z, None)
                    .map_err(|e| Error::InvalidInput(format!("entry ({i}, {j}): {e}")))?;
                raw.push((
                    i,
                    j,
                    univariate(&r.numerator, variable)?,
                    univariate(&r.denominator, variable)?,
                ));
            }
        }
        Self::from_entries(dim, &variable.to_string(), raw)
    }

    /// From the campaign's in-house dispersive input (`terms: [[i, j, num, den], ..]`,
    /// rational coefficient strings from degree 0, polynomials in eta).
    pub fn from_terms_json(dim: usize, terms: &serde_json::Value) -> Result<Self> {
        let list = terms
            .as_array()
            .ok_or_else(|| Error::InvalidInput("terms must be an array".into()))?;
        let mut raw = Vec::with_capacity(list.len());
        for t in list {
            let idx = |k: usize| {
                t[k].as_u64()
                    .map(|v| v as usize)
                    .ok_or_else(|| Error::InvalidInput(format!("bad term index in {t}")))
            };
            let poly = |k: usize| -> Result<Vec<rug::Rational>> {
                t[k].as_array()
                    .ok_or_else(|| Error::InvalidInput("term polynomial must be an array".into()))?
                    .iter()
                    .map(json_rational)
                    .collect()
            };
            let (n, d) = clear_denominators(&poly(2)?, &poly(3)?);
            raw.push((idx(0)?, idx(1)?, n, d));
        }
        Self::from_entries(dim, "eta", raw)
    }

    pub fn to_json(&self) -> serde_json::Value {
        serde_json::to_value(SystemFile {
            schema: SYSTEM_SCHEMA.into(),
            dim: self.dim,
            variable: self.variable.clone(),
            polys: self
                .polys
                .iter()
                .map(|p| p.iter().map(|c| c.to_string()).collect())
                .collect(),
            entries: self
                .entries
                .iter()
                .map(|e| (e.row, e.col, e.num, e.den, e.scale.to_string()))
                .collect(),
        })
        .expect("system serialises")
    }

    pub fn from_json(value: serde_json::Value) -> Result<Self> {
        let f: SystemFile = serde_json::from_value(value)
            .map_err(|e| Error::InvalidInput(format!("system file: {e}")))?;
        if f.schema != SYSTEM_SCHEMA {
            return Err(Error::InvalidInput(format!("system schema {}", f.schema)));
        }
        let polys = f
            .polys
            .iter()
            .map(|p| {
                p.iter()
                    .map(|c| {
                        c.parse::<Integer>()
                            .map_err(|e| Error::InvalidInput(e.to_string()))
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?;
        let raw = f
            .entries
            .iter()
            .map(|e| {
                let get = |k: usize| {
                    polys
                        .get(k)
                        .cloned()
                        .ok_or_else(|| Error::InvalidInput("poly index".into()))
                };
                let scale = parse_rational(&e.4)?;
                // N/D = scale * N' / D'
                let n = get(e.2)?.into_iter().map(|c| c * scale.numer()).collect();
                let d = get(e.3)?.into_iter().map(|c| c * scale.denom()).collect();
                Ok((e.0, e.1, n, d))
            })
            .collect::<Result<Vec<_>>>()?;
        let s = Self::from_entries(f.dim, &f.variable, raw)?;
        Ok(s)
    }

    /// blake3 of the canonical content (dimension, entries with their
    /// polynomials). Independent of the variable name and of poly numbering.
    pub fn fingerprint(&self) -> String {
        let mut h = blake3::Hasher::new();
        h.update(SYSTEM_SCHEMA.as_bytes());
        h.update(&(self.dim as u64).to_le_bytes());
        for e in &self.entries {
            h.update(&(e.row as u64).to_le_bytes());
            h.update(&(e.col as u64).to_le_bytes());
            h.update(e.scale.to_string().as_bytes());
            h.update(b";");
            for p in [&self.polys[e.num], &self.polys[e.den]] {
                h.update(&(p.len() as u64).to_le_bytes());
                for c in p {
                    h.update(c.to_string().as_bytes());
                    h.update(b",");
                }
            }
        }
        h.finalize().to_hex().to_string()
    }

    /// Multiply-adds per Taylor order of the recurrence (cost model).
    pub fn products_per_order(&self) -> usize {
        let mut groups = std::collections::BTreeSet::new();
        let mut cost = 0;
        for e in &self.entries {
            cost += self.polys[e.num].len();
            if groups.insert((e.row, e.den)) {
                cost += self.polys[e.den].len() - 1;
            }
        }
        cost
    }

    /// Exact comparison with another system on the same basis order up to a
    /// permutation `perm[i]` = index in `other` of this system's component i:
    /// returns the entries (i, j) where N1 D2 != N2 D1 or presence differs.
    pub fn differences(&self, other: &Self, perm: &[usize]) -> Vec<(usize, usize)> {
        let map = |s: &Self| -> BTreeMap<(usize, usize), (usize, usize, rug::Rational)> {
            s.entries
                .iter()
                .map(|e| ((e.row, e.col), (e.num, e.den, e.scale.clone())))
                .collect()
        };
        let (a, b) = (map(self), map(other));
        let mut out = Vec::new();
        let mut seen = std::collections::BTreeSet::new();
        for (&(i, j), (n1, d1, s1)) in &a {
            let key = (perm[i], perm[j]);
            seen.insert(key);
            match b.get(&key) {
                Some((n2, d2, s2)) => {
                    // s1 N1 D2 == s2 N2 D1
                    let f1 = Integer::from(s1.numer() * s2.denom());
                    let f2 = Integer::from(s2.numer() * s1.denom());
                    let scale = |p: Vec<Integer>, f: &Integer| {
                        p.into_iter().map(|c| c * f).collect::<Vec<_>>()
                    };
                    let lhs = scale(poly_mul(&self.polys[*n1], &other.polys[*d2]), &f1);
                    let rhs = scale(poly_mul(&other.polys[*n2], &self.polys[*d1]), &f2);
                    if lhs != rhs {
                        out.push((i, j));
                    }
                }
                None => out.push((i, j)),
            }
        }
        let inverse = {
            let mut inv = vec![usize::MAX; perm.len()];
            for (i, &p) in perm.iter().enumerate() {
                inv[p] = i;
            }
            inv
        };
        for &(i, j) in b.keys() {
            if !seen.contains(&(i, j)) {
                out.push((inverse[i], inverse[j]));
            }
        }
        out
    }
}

impl RationalSystem {
    /// From RustFlow's own differential system (exact rational rows).
    pub fn from_differential_system(system: &crate::DifferentialSystem) -> Result<Self> {
        Self::from_matrix_atoms(&system.matrix, system.variable)
    }

    /// Finite singular points of the system (roots of all denominators) in
    /// the system variable, as f64 pairs: every distinct denominator is
    /// reduced by the irreducible factors found so far (exact division) and
    /// only the remaining cofactor is factored (Symbolica).
    pub fn singular_points(&self) -> Result<Vec<(f64, f64)>> {
        let x = symbol!("rustflow_patch_kernel::x");
        let mut factors: Vec<Vec<Integer>> = Vec::new();
        let mut dens: Vec<usize> = self.entries.iter().map(|e| e.den).collect();
        dens.sort_unstable();
        dens.dedup();
        dens.sort_by_key(|&d| self.polys[d].len());
        for d in dens {
            let mut rest = self.polys[d].clone();
            for f in &factors {
                while rest.len() > f.len() - 1 && rest.len() > 1 {
                    match div_exact(&rest, f) {
                        Some(q) => rest = q,
                        None => break,
                    }
                }
            }
            if rest.len() <= 1 {
                continue;
            }
            let text = rest
                .iter()
                .enumerate()
                .map(|(k, c)| format!("({c})*rustflow_patch_kernel::x^{k}"))
                .collect::<Vec<_>>()
                .join("+");
            let atom = Atom::parse(&text, "rustflow_patch_kernel", Default::default())
                .map_err(|e| Error::InvalidInput(format!("denominator polynomial: {e}")))?;
            let r: symbolica::domains::rational_polynomial::RationalPolynomial<
                symbolica::domains::integer::IntegerRing,
                u16,
            > = atom
                .try_to_rational_polynomial(&Q, &Z, None)
                .map_err(|e| Error::InvalidInput(e.to_string()))?;
            for (f, _) in r.numerator.factor() {
                let mut c = univariate(&f, x)?;
                if c.len() <= 1 {
                    continue;
                }
                if c.last().is_some_and(|v| *v < 0) {
                    for v in c.iter_mut() {
                        *v = -std::mem::take(v);
                    }
                }
                if !factors.contains(&c) {
                    factors.push(c);
                }
            }
        }
        let p = crate::Precision::decimal(40)?;
        let mut out = Vec::new();
        for f in &factors {
            let coefficients: Vec<crate::ComplexFloat> = f
                .iter()
                .map(|c| {
                    crate::ComplexFloat::new(
                        symbolica::domains::float::Float::from(rug::Float::with_val(p.bits, c)),
                        symbolica::domains::float::Float::from(rug::Float::new(p.bits)),
                    )
                })
                .collect();
            for root in crate::ode::polynomial_roots(p, &coefficients, x)? {
                out.push((root.re.to_f64(), root.im.to_f64()));
            }
        }
        out.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
        Ok(out)
    }
}

/// Exact quotient `d / f` of integer polynomials if `f` divides `d` over Z.
fn div_exact(d: &[Integer], f: &[Integer]) -> Option<Vec<Integer>> {
    let (nd, nf) = (d.len() - 1, f.len() - 1);
    if nd < nf {
        return None;
    }
    let lc = &f[nf];
    let mut r = d.to_vec();
    let mut q = vec![Integer::new(); nd - nf + 1];
    for i in (0..=nd - nf).rev() {
        let c = &r[i + nf];
        if !c.is_divisible(lc) {
            return None;
        }
        let qi = Integer::from(c / lc);
        for (j, fj) in f.iter().enumerate() {
            r[i + j] -= Integer::from(&qi * fj);
        }
        q[i] = qi;
    }
    r.iter().all(|c| *c == 0).then_some(q)
}

fn poly_mul(a: &[Integer], b: &[Integer]) -> Vec<Integer> {
    let mut out = vec![Integer::new(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] += Integer::from(x * y);
        }
    }
    trim(&mut out);
    out
}

/// Dense integer coefficients of a polynomial in (at most) `variable`.
fn univariate(
    p: &symbolica::poly::polynomial::MultivariatePolynomial<
        symbolica::domains::integer::IntegerRing,
        u16,
    >,
    variable: Symbol,
) -> Result<Vec<Integer>> {
    let vars = p.get_vars_ref();
    let nv = vars.len();
    let mut position = None;
    for (k, v) in vars.iter().enumerate() {
        match v {
            symbolica::poly::PolyVariable::Symbol(s) if *s == variable => position = Some(k),
            other => {
                // Another variable may appear in the list only with zero exponents.
                let used = (0..p.nterms()).any(|t| p.exponents[t * nv + k] != 0);
                if used {
                    return Err(Error::InvalidInput(format!(
                        "coefficient depends on {other:?}, not only on {variable}"
                    )));
                }
            }
        }
    }
    let mut out: Vec<Integer> = vec![Integer::new()];
    for t in 0..p.nterms() {
        let deg = position
            .map(|k| p.exponents[t * nv + k] as usize)
            .unwrap_or(0);
        if out.len() <= deg {
            out.resize(deg + 1, Integer::new());
        }
        out[deg] += p.coefficients[t]
            .to_string()
            .parse::<Integer>()
            .map_err(|e| Error::InvalidInput(e.to_string()))?;
    }
    trim(&mut out);
    Ok(out)
}
