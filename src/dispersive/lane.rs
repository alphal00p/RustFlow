//! Lower-loop soft boundaries for a recursive use of the dispersive mode (REVIEW-B BLOCKER 1): the 4L prism lane's soft boundary
//! from the RustFlow P14 values only.
//!  * `values3l`: port of values/d5l/v9a/tools/extract3l.py: the soft system of a finished W4 class (3039) has rank 14 of 14 on the 40
//!    P14 endpoints; its engine masters T x M5, T x B4, T^2 x S3 (identified by the matroid components of their line momenta) give the
//!    3L masters M5, B4, S3 (T = Gamma(1 - D/2), the 1L tadpole E-value at unit mass).
//!  * `lane_v`: port of values/d5l/v9a/tools/soft4.py solve: lane soft V = sum c / D^dp x red(embedded integral) over the lane's soft
//!    reduction, whose soft masters are known products of T and the 3L masters (config soft_masters: {type, T}).
use super::boundary;
use super::softplan;
use super::system::parse_coefficient;
use rug::{Float, Rational};
use serde_json::{Value, json};
use std::collections::{BTreeMap, HashMap};

fn vi(x: &Value) -> Vec<i64> {
    x.as_array()
        .unwrap()
        .iter()
        .map(|a| a.as_i64().unwrap())
        .collect()
}

fn t_e(d0: &Rational, prec: u32) -> Float {
    Float::with_val(prec, 1 - Float::with_val(prec, d0) / 2u32).gamma()
}

/// 3L masters {M5, B4, S3} from the W4-class soft system and the P14 endpoints (only the non-'PR' rows are used)
pub fn values3l(
    fam5: &Value,
    plan: &Value,
    sol: &Value,
    sample: &Value,
    e_slot: usize,
    prec: u32,
) -> Result<(BTreeMap<String, Float>, Value), String> {
    let d0: Rational = sol["fixed"]["d"]
        .as_str()
        .ok_or("no fixed d")?
        .parse()
        .map_err(|e: rug::rational::ParseRationalError| e.to_string())?;
    let masters: Vec<Vec<i64>> = sol["masters"]
        .as_array()
        .ok_or("no masters")?
        .iter()
        .map(vi)
        .collect();
    // unit soft rows: V(m) = the E-value of engine master m (embed sets the e slot to 1: give nu with the e slot at 0)
    let mut unit = serde_json::Map::new();
    for m in &masters {
        if m[e_slot] != 1 {
            continue;
        }
        let mut nu = m.clone();
        nu[e_slot] = 0;
        unit.insert(serde_json::to_string(m).unwrap(), json!([[nu, 0, "1"]]));
    }
    let plan_u = json!({"endpoint_polys": plan["endpoint_polys"], "soft_polys": unit});
    let mut ep = Vec::new();
    for x in sample["endpoint"]
        .as_array()
        .ok_or("sample has no endpoints")?
    {
        if !x
            .get("family")
            .and_then(|f| f.as_str())
            .unwrap_or("")
            .is_empty()
        {
            continue; // P14 rows only
        }
        let t = vi(&x["target"]);
        let key = serde_json::to_string(&t).unwrap().replace(',', ", ");
        ep.push((
            key,
            t,
            Float::with_val(
                prec,
                Float::parse(x["value"]["re"].as_str().unwrap()).map_err(|e| e.to_string())?,
            ),
        ));
    }
    let (vals, diag) = boundary::solve_v(&plan_u, sol, &ep, e_slot, prec)?;
    let q = softplan::slot_momenta(fam5);
    let te = t_e(&d0, prec);
    let mut by_type: BTreeMap<String, Float> = BTreeMap::new();
    let mut found = serde_json::Map::new();
    for m in &masters {
        let Some(vs) = vals.get(&serde_json::to_string(m).unwrap()) else {
            continue;
        };
        if !m.iter().all(|&p| p == 0 || p == 1) {
            continue;
        }
        let lines: Vec<usize> = (0..m.len()).filter(|&s| m[s] > 0 && s != e_slot).collect();
        let vecs: Vec<Vec<Rational>> = lines
            .iter()
            .map(|&s| q[s].clone().ok_or("slot without momentum"))
            .collect::<Result<_, _>>()?;
        let comps = softplan::components(&vecs);
        let v = Float::with_val(prec, Float::parse(vs).map_err(|e| e.to_string())?);
        let (name, tpow) = match comps.as_slice() {
            [(1, 1), (5, 3)] => ("M5", 1),
            [(1, 1), (4, 3)] => ("B4", 1),
            [(1, 1), (1, 1), (3, 2)] => ("S3", 2),
            _ => continue,
        };
        let mut x = v;
        for _ in 0..tpow {
            x /= &te;
        }
        found.insert(name.into(), json!(m));
        by_type.insert(name.into(), x);
    }
    for t in ["M5", "B4", "S3"] {
        if !by_type.contains_key(t) {
            return Err(format!(
                "3L master {t} not found among the soft engine masters"
            ));
        }
    }
    Ok((
        by_type,
        json!({"endpoint_system": diag, "p14_rows": ep.len(), "masters_used": found}),
    ))
}

/// the lane's soft V from its soft plan, its soft reduction and the known 3L masters
pub fn lane_v(
    lane_plan: &Value,
    lane_config: &Value,
    sol: &Value,
    by_type: &BTreeMap<String, Float>,
    prec: u32,
) -> Result<BTreeMap<String, String>, String> {
    let e_slot = lane_config["e"].as_u64().ok_or("lane config has no e")? as usize;
    let d0: Rational = sol["fixed"]["d"]
        .as_str()
        .ok_or("no fixed d")?
        .parse()
        .map_err(|e: rug::rational::ParseRationalError| e.to_string())?;
    let te = t_e(&d0, prec);
    let masters: Vec<Vec<i64>> = sol["masters"]
        .as_array()
        .ok_or("no masters")?
        .iter()
        .map(vi)
        .collect();
    let mut mval: HashMap<Vec<i64>, Float> = HashMap::new();
    for m in &masters {
        let key = serde_json::to_string(m).unwrap().replace(',', ", ");
        let k = lane_config["soft_masters"]
            .get(&key)
            .ok_or_else(|| format!("lane soft master {key} has no known value"))?;
        let mut v = match k["type"].as_str() {
            Some("1") => Float::with_val(prec, 1),
            Some(t) => by_type
                .get(t)
                .ok_or_else(|| format!("no value of type {t}"))?
                .clone(),
            None => return Err("soft master without type".into()),
        };
        for _ in 0..k.get("T").and_then(|x| x.as_u64()).unwrap_or(0) {
            v *= &te;
        }
        mval.insert(m.clone(), v);
    }
    let mut rules: HashMap<Vec<i64>, Vec<(Vec<i64>, Rational)>> = HashMap::new();
    for r in sol["rules"].as_array().ok_or("no rules")? {
        let mut terms = Vec::new();
        for t in r["terms"].as_array().unwrap() {
            let (n, d) = parse_coefficient(t["coefficient"].as_str().unwrap(), "eta")?;
            if n.len() > 1 || d.len() > 1 {
                return Err("eta-dependent coefficient in a plain reduction".into());
            }
            terms.push((vi(&t["integral"]), Rational::from(&n[0] / &d[0])));
        }
        rules.insert(vi(&r["target"]), terms);
    }
    let red = |x: &[i64]| -> Result<Float, String> {
        if let Some(v) = mval.get(x) {
            return Ok(v.clone());
        }
        let mut tot = Float::new(prec);
        for (mu, c) in rules
            .get(x)
            .ok_or_else(|| format!("no reduction rule for {x:?}"))?
        {
            let sgn = if (mu.iter().sum::<i64>() - x.iter().sum::<i64>()).rem_euclid(2) == 1 {
                -1
            } else {
                1
            };
            tot += Float::with_val(prec, &mval[mu] * Float::with_val(prec, c)) * sgn;
        }
        Ok(tot)
    };
    let digits = (prec as f64 * super::DIGITS_PER_BIT) as usize;
    let mut out = BTreeMap::new();
    for (key, poly) in lane_plan["soft_polys"].as_object().ok_or("no soft_polys")? {
        let mut val = Float::new(prec);
        for term in poly.as_array().unwrap() {
            let mut nu = vi(&term[0]);
            let dp = term[1].as_i64().unwrap();
            let mut c: Rational = term[2]
                .as_str()
                .unwrap()
                .parse()
                .map_err(|e: rug::rational::ParseRationalError| e.to_string())?;
            for _ in 0..dp {
                c /= &d0;
            }
            if nu[e_slot] != 0 {
                return Err("soft integral uses the eta slot".into());
            }
            nu[e_slot] = 1;
            val += Float::with_val(prec, c) * red(&nu)?;
        }
        out.insert(key.clone(), val.to_string_radix(10, Some(digits)));
    }
    Ok(out)
}
