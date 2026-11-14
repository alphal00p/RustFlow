//! Soft boundary rows at z = 0 (oracle: prep5.py 'boundary values', soft.py / soft_gen.py 'solve').
//!
//! Soft region of the eta line (momentum p -> 0) of the e-block master j = (e^a, rest R):
//!   c_{j, n' - a} = (-1)^(sum j) e^(L gamma eps) Gamma(n' + D/2) Gamma(a - n' - D/2) / (Gamma(D/2) Gamma(a)) V_{R_j, n'}
//! (campaign normalisation e^(L gamma eps), Minkowski storage; D = d0), where V_{R,n'} (n' = 0, 1) are 4L all-massive
//! vacuum integrals (Euclidean, no e^(gamma eps)). Two sources for V:
//!   * READ: runs/boundary/V-a<att>-<idx>.json ('soft_V_euclid_nogamma', written by the in-house soft_gen.py solve);
//!   * NATIVE (`solve_v`): soft-inputs/plan.json (exact soft / endpoint polynomials over the 4L SUB family) + the RustRed
//!     reduction of those integrals at d0 (soft/solution-<d0>.json) + the 4L endpoint values of the RustFlow V2 sample
//!     (values/d5l/v9a/pr/samples/attempt-<a>/sample-<idx>.json): engine masters by exact elimination of the endpoint
//!     system, consistency of the redundant equations = gate, then the soft combinations.
use super::system::{parse_coefficient, EBlockSystem};
use rug::{Float, Rational};
use std::collections::{BTreeMap, HashMap};

#[derive(Clone, Debug)]
pub struct BoundaryRow {
    pub comp: usize,
    pub n: i64,
    pub value: Float,
}

/// V values keyed by "[m]|n'" (decimal strings) -> boundary rows c_{j, n' - a}
pub fn rows_from_v(sys: &EBlockSystem, v: &BTreeMap<String, String>, loops: u32, prec: u32) -> Result<Vec<BoundaryRow>, String> {
    let idx: HashMap<&Vec<i64>, usize> = sys.masters.iter().enumerate().map(|(i, m)| (m, i)).collect();
    let d = Float::with_val(prec, &sys.d0);
    let eps = Float::with_val(prec, (Float::with_val(prec, 4) - &d) / 2);
    let egam = Float::with_val(prec, Float::with_val(prec, rug::float::Constant::Euler) * &eps * loops).exp();
    let dh = Float::with_val(prec, &d / 2);
    let g_dh = dh.clone().gamma();
    let mut out = Vec::new();
    let mut unmatched = Vec::new();
    for (key, val) in v {
        let (mv, np) = key.split_once('|').ok_or("bad V key")?;
        let m: Vec<i64> = serde_json::from_str(mv).map_err(|e| e.to_string())?;
        let np: i64 = np.parse().map_err(|e: std::num::ParseIntError| e.to_string())?;
        let Some(&j) = idx.get(&m) else {
            unmatched.push(key.clone());
            continue;
        };
        let a = m[sys.e_slot];
        let g1 = Float::with_val(prec, &dh + np).gamma();
        let g2 = Float::with_val(prec, Float::with_val(prec, a - np) - &dh).gamma();
        let ga = Float::with_val(prec, a).gamma();
        let k = Float::with_val(prec, &g1 * &g2) / Float::with_val(prec, &g_dh * &ga);
        let vv = Float::parse(val).map(|x| Float::with_val(prec, x)).map_err(|e| e.to_string())?;
        let mut c = Float::with_val(prec, &egam * &k) * vv;
        if m.iter().sum::<i64>().rem_euclid(2) == 1 {
            c = -c;
        }
        out.push(BoundaryRow { comp: j, n: np - a, value: c });
    }
    if !unmatched.is_empty() {
        // REVIEW-B MAJOR 4: a V key that is not an e-block master means a wrong e slot or a V file of another system
        return Err(format!("{} V keys are not e-block masters of this system (wrong e slot or V file?), e.g. {}", unmatched.len(), unmatched[0]));
    }
    Ok(out)
}

pub fn read_v_file(path: &str) -> Result<(BTreeMap<String, String>, serde_json::Value), String> {
    let v: serde_json::Value = serde_json::from_slice(&std::fs::read(path).map_err(|e| format!("{path}: {e}"))?).map_err(|e| e.to_string())?;
    let mut out = BTreeMap::new();
    for (k, x) in v["soft_V_euclid_nogamma"].as_object().ok_or("no soft_V_euclid_nogamma")? {
        out.insert(k.clone(), x.as_str().ok_or("V value not a string")?.to_string());
    }
    let mut meta = v.clone();
    meta.as_object_mut().unwrap().remove("soft_V_euclid_nogamma");
    Ok((out, meta))
}

// ---------------- native V (port of soft_gen.py solve) ----------------
fn parse_plain_coef(text: &str) -> Result<Rational, String> {
    let (n, d) = parse_coefficient(text, "eta")?;
    if n.len() > 1 || d.len() > 1 {
        return Err(format!("eta-dependent coefficient in a plain reduction: {text}"));
    }
    Ok(Rational::from(&n[0] / &d[0]))
}

type Lin = BTreeMap<usize, Rational>;

/// V_{R,n'} from the 4L inputs. plan = soft-inputs/plan.json, sol = RustRed reduction of the plan's integrals at d0,
/// endpoint = [(family prefix + json target, target, value)] of the 4L sample, in the sample's order. Returns (V map, diagnostics).
pub fn solve_v(plan: &serde_json::Value, sol: &serde_json::Value, endpoint: &[(String, Vec<i64>, Float)], e_slot: usize, prec: u32) -> Result<(BTreeMap<String, String>, serde_json::Value), String> {
    let d0 = sol["fixed"]["d"].as_str().ok_or("no fixed d")?.parse::<Rational>().map_err(|e| e.to_string())?;
    let vec_i = |x: &serde_json::Value| -> Vec<i64> { x.as_array().unwrap().iter().map(|a| a.as_i64().unwrap()).collect() };
    let masters: Vec<Vec<i64>> = sol["masters"].as_array().ok_or("no masters")?.iter().map(vec_i).collect();
    let midx: HashMap<Vec<i64>, usize> = masters.iter().enumerate().map(|(i, m)| (m.clone(), i)).collect();
    let mut rules: HashMap<Vec<i64>, Vec<(Vec<i64>, Rational)>> = HashMap::new();
    for r in sol["rules"].as_array().ok_or("no rules")? {
        let terms = r["terms"].as_array().unwrap().iter().map(|t| Ok((vec_i(&t["integral"]), parse_plain_coef(t["coefficient"].as_str().unwrap())?))).collect::<Result<Vec<_>, String>>()?;
        rules.insert(vec_i(&r["target"]), terms);
    }
    let embed = |nu: &[i64]| -> Result<Vec<i64>, String> {
        let mut x = nu.to_vec();
        if x[e_slot] != 0 {
            return Err("plan integral uses the eta slot".into());
        }
        x[e_slot] = 1;
        Ok(x)
    };
    // E_sub(nu) of the embedded integral: sign (-1)^(sum mu - sum nu) per term (Euclidean <-> Minkowski storage)
    let red = |x: &[i64]| -> Result<Lin, String> {
        if let Some(&i) = midx.get(x) {
            return Ok(BTreeMap::from([(i, Rational::from(1))]));
        }
        let mut out = Lin::new();
        for (mu, c) in rules.get(x).ok_or_else(|| format!("no reduction rule for {x:?}"))? {
            let sgn = if (mu.iter().sum::<i64>() - x.iter().sum::<i64>()).rem_euclid(2) == 1 { -1 } else { 1 };
            *out.entry(midx[mu]).or_default() += Rational::from(c * sgn);
        }
        Ok(out)
    };
    let lin = |poly: &serde_json::Value| -> Result<Lin, String> {
        let mut out = Lin::new();
        for term in poly.as_array().ok_or("bad poly")? {
            let nu = vec_i(&term[0]);
            let dp = term[1].as_i64().unwrap();
            let mut c = term[2].as_str().unwrap().parse::<Rational>().map_err(|e| e.to_string())?;
            for _ in 0..dp {
                c /= &d0;
            }
            for (k, v) in red(&embed(&nu)?)? {
                *out.entry(k).or_default() += Rational::from(&c * &v);
            }
        }
        out.retain(|_, v| *v != 0);
        Ok(out)
    };
    // endpoint equations: rows (exact) = values (Float); the oracle orders them as the sample's endpoint list
    let mut piv: Vec<(usize, Lin, Float)> = Vec::new();
    let reduce_row = |piv: &[(usize, Lin, Float)], r: &Lin, val: &Float| -> (Lin, Float) {
        let mut r = r.clone();
        let mut val = val.clone();
        for (c, br, bv) in piv {
            if let Some(f) = r.get(c).cloned() {
                if f != 0 {
                    for (k, x) in br {
                        *r.entry(*k).or_default() -= Rational::from(&f * x);
                    }
                    val -= Float::with_val(prec, &f) * bv;
                    r.retain(|_, x| *x != 0);
                }
            }
        }
        (r, val)
    };
    let mut cons: Vec<(Float, Float)> = Vec::new();
    let mut scale = Float::new(prec);
    let eps_poly = &plan["endpoint_polys"];
    for (key, t, raw) in endpoint {
        let poly = eps_poly.get(key).ok_or_else(|| format!("endpoint target {key} not in plan"))?;
        let row = lin(poly)?;
        let mut v = raw.clone();
        if t.iter().sum::<i64>().rem_euclid(2) == 1 {
            v = -v;
        }
        if Float::with_val(prec, v.clone().abs()) > scale {
            scale = Float::with_val(prec, v.clone().abs());
        }
        let (rr, vv) = reduce_row(&piv, &row, &v);
        if rr.is_empty() {
            cons.push((vv, v));
            continue;
        }
        let c = *rr.keys().next().unwrap();
        let f = rr[&c].clone();
        let br: Lin = rr.iter().map(|(k, x)| (*k, Rational::from(x / &f))).collect();
        let bv = Float::with_val(prec, &vv / Float::with_val(prec, &f));
        let mut newp = Vec::with_capacity(piv.len() + 1);
        for (c2, br2, bv2) in piv.drain(..) {
            if let Some(g) = br2.get(&c).cloned() {
                let mut br2 = br2;
                for (k, x) in &br {
                    *br2.entry(*k).or_default() -= Rational::from(&g * x);
                }
                br2.retain(|_, x| *x != 0);
                let bv2 = bv2 - Float::with_val(prec, &g) * &bv;
                newp.push((c2, br2, bv2));
            } else {
                newp.push((c2, br2, bv2));
            }
        }
        newp.push((c, br, bv));
        piv = newp;
    }
    let floor = Float::with_val(prec, &scale * Float::with_val(prec, Float::i_exp(1, -200)));
    let mut worst = -1e9f64;
    for (vv, v) in &cons {
        let den = if Float::with_val(prec, v.clone().abs()) > floor { Float::with_val(prec, v.clone().abs()) } else { floor.clone() };
        let rel = Float::with_val(prec, vv.clone().abs() / den);
        if !rel.is_zero() {
            worst = worst.max(rel.log10().to_f64());
        }
    }
    let digits = (prec as f64 * 0.30103) as usize;
    let mut out = BTreeMap::new();
    let mut missing = Vec::new();
    for (key, poly) in plan["soft_polys"].as_object().ok_or("no soft_polys")? {
        let (rr, vv) = reduce_row(&piv, &lin(poly)?, &Float::new(prec));
        if !rr.is_empty() {
            missing.push(key.clone());
            continue;
        }
        out.insert(key.clone(), Float::with_val(prec, -vv).to_string_radix(10, Some(digits)));
    }
    let diag = serde_json::json!({"d0": d0.to_string(), "n_endpoint_eqs": endpoint.len(), "n_engine_masters": masters.len(),
        "rank_endpoint_system": piv.len(), "n_consistency_checks": cons.len(), "endpoint_consistency_worst_rel_log10": worst,
        "soft_missing": missing});
    Ok((out, diag))
}
