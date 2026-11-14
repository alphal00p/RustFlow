//! Recursive soft boundaries (REVIEW-B BLOCKER 1, W3c): the soft boundary of an L-loop LANE family (its eta line soft) built only
//! from endpoint values that RustFlow computed, plus exact RustRed reductions and exact Gamma functions.
//!
//! Chain (inputs named by a chain spec, see `lane_boundary`):
//!  1. LOWER system: the native soft plan (softplan::build) of a finished class whose soft graph is an 8-line 4L graph maps the P14
//!     endpoint targets of the 4L sample onto that graph's lines; with the class's RustRed soft reduction the P14 endpoint values fix
//!     the E-values of its soft engine masters (boundary::solve_v, exact elimination, consistency gate on the redundant rows).
//!  2. FACTORISATION: an engine master whose powers are all 0/1 is a unit-mass vacuum integral; the line momenta (eta line excluded)
//!     split into connected components of their linear matroid. A one-line component is a 1L tadpole T = Gamma(1 - D/2) (E-value,
//!     unit mass); every other component is a lower-loop vacuum master, identified by the ISOMORPHISM CLASS of its matroid
//!     (`canonical_matroid`): with every power 1 and unit masses the Feynman-parametric integrand depends only on U = sum over bases
//!     of prod x (unimodular routings) and on sum x, so isomorphic matroids give the same integral. A master with exactly one
//!     non-tadpole component gives that component's value = E-value / T^k.
//!  3. LANE: the soft plan of the lane family (softplan::build_soft_only), the lane's RustRed soft reduction, and its soft masters
//!     factorised the same way (T^k x prod of known components): every soft master must be fully determined, and a declared
//!     config "soft_masters" classification (in-house) is only CROSS-CHECKED, never used.
//! No number in this chain comes from a dispersive run: the inputs are the P14 endpoint values (RustFlow eta-flows + spectral gluing),
//! RustRed reductions, Gamma functions.
use super::boundary;
use super::softplan;
use super::system::parse_coefficient;
use rug::{Float, Rational};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet, HashMap};

type Vq = Vec<Rational>;

fn vi(x: &Value) -> Vec<i64> {
    x.as_array().unwrap().iter().map(|a| a.as_i64().unwrap()).collect()
}

fn rd(p: &str) -> Result<Value, String> {
    serde_json::from_slice(&std::fs::read(p).map_err(|e| format!("{p}: {e}"))?).map_err(|e| format!("{p}: {e}"))
}

/// exact row echelon: (pivot column, row) list of a set of vectors
fn echelon(vecs: &[&Vq]) -> Vec<(usize, Vq)> {
    let mut ech: Vec<(usize, Vq)> = Vec::new();
    for v in vecs {
        let mut v = (*v).clone();
        for (c, r) in &ech {
            if v[*c] != 0 {
                let f = Rational::from(&v[*c] / &r[*c]);
                for (x, y) in v.iter_mut().zip(r) {
                    *x -= Rational::from(&f * y);
                }
            }
        }
        if let Some(c) = v.iter().position(|x| *x != 0) {
            ech.push((c, v));
        }
    }
    ech
}

fn rank(vecs: &[&Vq]) -> usize {
    echelon(vecs).len()
}

/// connected components of the linear matroid of `vecs` (index groups, each sorted, groups sorted by first index):
/// union of the fundamental circuits of a greedy basis
pub fn component_groups(vecs: &[Vq]) -> Result<Vec<Vec<usize>>, String> {
    let n = vecs.len();
    if vecs.iter().any(|v| v.iter().all(|x| *x == 0)) {
        return Err("a line with zero momentum".into());
    }
    let mut basis: Vec<usize> = Vec::new();
    for i in 0..n {
        let mut s: Vec<&Vq> = basis.iter().map(|&b| &vecs[b]).collect();
        s.push(&vecs[i]);
        if rank(&s) == s.len() {
            basis.push(i);
        }
    }
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        let mut x = x;
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    for i in 0..n {
        if basis.contains(&i) {
            continue;
        }
        // basis element b is in the fundamental circuit of i iff removing b from the basis drops i out of the span
        for &b in &basis {
            let s: Vec<&Vq> = basis.iter().filter(|&&x| x != b).map(|&x| &vecs[x]).chain([&vecs[i]]).collect();
            if rank(&s) == basis.len() {
                let (x, y) = (find(&mut parent, i), find(&mut parent, b));
                parent[x] = y;
            }
        }
    }
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push(i);
    }
    let mut out: Vec<Vec<usize>> = groups.into_values().collect();
    out.sort();
    Ok(out)
}

/// canonical form of a (connected) linear matroid: "n<size>r<rank>:" + the lexicographically smallest sorted list of basis
/// bitmasks over all relabellings. Requires a unimodular-like set: every nonzero r x r minor (in coordinates of the span) has the
/// same absolute value, so that U = c x sum_bases prod x with one constant.
pub fn canonical_matroid(vecs: &[Vq]) -> Result<String, String> {
    let n = vecs.len();
    if n > 8 {
        return Err(format!("component with {n} lines: canonical form not implemented above 8"));
    }
    let all: Vec<&Vq> = vecs.iter().collect();
    let ech = echelon(&all);
    let r = ech.len();
    let cols: Vec<usize> = ech.iter().map(|(c, _)| *c).collect();
    // minors on the pivot coordinates (the projection onto them is an isomorphism on the span)
    let minor = |sub: &[usize]| -> Rational {
        let mut m: Vec<Vec<Rational>> = sub.iter().map(|&i| cols.iter().map(|&c| vecs[i][c].clone()).collect()).collect();
        let k = m.len();
        let mut d = Rational::from(1);
        for i in 0..k {
            let Some(p) = (i..k).find(|&rr| m[rr][i] != 0) else { return Rational::new() };
            if p != i {
                m.swap(i, p);
                d = -d;
            }
            let pv = m[i][i].clone();
            d *= &pv;
            for rr in (i + 1)..k {
                if m[rr][i] != 0 {
                    let f = Rational::from(&m[rr][i] / &pv);
                    let ri = m[i].clone();
                    for (x, y) in m[rr].iter_mut().zip(&ri) {
                        *x -= Rational::from(&f * y);
                    }
                }
            }
        }
        d
    };
    let mut bases: Vec<u32> = Vec::new();
    let mut absdet: Option<Rational> = None;
    for mask in 0u32..(1 << n) {
        if mask.count_ones() as usize != r {
            continue;
        }
        let sub: Vec<usize> = (0..n).filter(|i| mask >> i & 1 == 1).collect();
        let d = minor(&sub).abs();
        if d != 0 {
            match &absdet {
                None => absdet = Some(d),
                Some(a) if *a == d => {}
                Some(a) => return Err(format!("minors {a} and {d} differ: not a unimodular line set")),
            }
            bases.push(mask);
        }
    }
    // smallest relabelling (all permutations, n <= 8)
    let mut perm: Vec<usize> = (0..n).collect();
    let mut best: Option<Vec<u32>> = None;
    loop {
        let mut img: Vec<u32> = bases.iter().map(|&b| (0..n).filter(|i| b >> i & 1 == 1).map(|i| 1u32 << perm[i]).sum()).collect();
        img.sort();
        if best.as_ref().is_none_or(|b| img < *b) {
            best = Some(img);
        }
        // next permutation
        let Some(i) = (0..n.saturating_sub(1)).rev().find(|&i| perm[i] < perm[i + 1]) else { break };
        let j = (i + 1..n).rev().find(|&j| perm[j] > perm[i]).unwrap();
        perm.swap(i, j);
        perm[i + 1..].reverse();
    }
    Ok(format!("n{n}r{r}:{:?}", best.unwrap_or_default()))
}

/// label of a canonical form for reports: the usual names of the 1-3 loop unit-mass vacuum masters by (lines, loops)
pub fn label(canon: &str) -> String {
    let nr = canon.split(':').next().unwrap_or("");
    match nr {
        "n1r1" => "T".into(),
        "n3r2" => "S3".into(),
        "n4r3" => "B4".into(),
        "n5r3" => "M5".into(),
        _ => nr.into(),
    }
}

/// factorisation of a 0/1-power master (eta line excluded): (power of T, canonical forms of the other components)
pub fn factorise(m: &[i64], e_slot: usize, mom: &[Option<Vq>], allowed: &BTreeSet<usize>) -> Result<(u32, Vec<String>), String> {
    if m.iter().any(|&p| p != 0 && p != 1) {
        return Err(format!("{m:?}: dots or numerators"));
    }
    let lines: Vec<usize> = (0..m.len()).filter(|&s| m[s] == 1 && s != e_slot).collect();
    if let Some(s) = lines.iter().find(|s| !allowed.contains(s)) {
        return Err(format!("{m:?}: line {s} is not a slot of the soft family"));
    }
    let vecs: Vec<Vq> = lines.iter().map(|&s| mom[s].clone().ok_or_else(|| format!("slot {s} has no momentum"))).collect::<Result<_, _>>()?;
    let mut t = 0;
    let mut parts = Vec::new();
    for g in component_groups(&vecs)? {
        let sub: Vec<Vq> = g.iter().map(|&i| vecs[i].clone()).collect();
        if g.len() == 1 {
            t += 1;
        } else {
            parts.push(canonical_matroid(&sub)?);
        }
    }
    parts.sort();
    Ok((t, parts))
}

fn t_e(d0: &Rational, prec: u32) -> Float {
    Float::with_val(prec, 1 - Float::with_val(prec, d0) / 2u32).gamma()
}

fn digits(a: &Float, b: &Float) -> f64 {
    let p = a.prec().max(b.prec());
    if b.is_zero() {
        return if a.is_zero() { 99.0 } else { -Float::with_val(p, a.abs_ref()).log10().to_f64() };
    }
    let d = Float::with_val(p, Float::with_val(p, a - b) / b).abs();
    if d.is_zero() { 99.0 } else { -d.log10().to_f64() }
}

/// lower-loop vacuum masters (by canonical matroid) from the lower system's soft plan + RustRed soft reduction + the P14 endpoint rows
/// (rows whose 'family' is empty) of a 4L sample
pub fn lower_values(fam: &Value, plan: &Value, sol: &Value, sample: &Value, e_slot: usize, prec: u32) -> Result<(BTreeMap<String, Float>, Value), String> {
    let d0: Rational = sol["fixed"]["d"].as_str().ok_or("no fixed d")?.parse().map_err(|e: rug::rational::ParseRationalError| e.to_string())?;
    let masters: Vec<Vec<i64>> = sol["masters"].as_array().ok_or("no masters")?.iter().map(vi).collect();
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
    for x in sample["endpoint"].as_array().ok_or("sample has no endpoints")? {
        if x.get("family").and_then(|f| f.as_str()).unwrap_or("") != "" {
            continue;
        }
        let t = vi(&x["target"]);
        let key = serde_json::to_string(&t).unwrap().replace(',', ", ");
        ep.push((key, t, Float::with_val(prec, Float::parse(x["value"]["re"].as_str().unwrap()).map_err(|e| e.to_string())?)));
    }
    let (vals, diag) = boundary::solve_v(&plan_u, sol, &ep, e_slot, prec)?;
    let mom = softplan::slot_momenta(fam);
    let allowed: BTreeSet<usize> = plan["derived"]["sub"].as_array().ok_or("lower plan has no derived.sub")?.iter().map(|x| x.as_u64().unwrap() as usize).collect();
    let te = t_e(&d0, prec);
    let mut out: BTreeMap<String, Float> = BTreeMap::new();
    let mut used = Vec::new();
    for m in &masters {
        let Some(vs) = vals.get(&serde_json::to_string(m).unwrap()) else { continue };
        let Ok((t, parts)) = factorise(m, e_slot, &mom, &allowed) else { continue };
        if parts.len() != 1 {
            continue;
        }
        let mut x = Float::with_val(prec, Float::parse(vs).map_err(|e| e.to_string())?);
        for _ in 0..t {
            x /= &te;
        }
        let c = parts[0].clone();
        let mut rec = json!({"master": m, "T_power": t, "component": c, "label": label(&c)});
        if let Some(prev) = out.get(&c) {
            rec["agrees_with_previous_digits"] = json!((digits(&x, prev) * 10.0).round() / 10.0);
        } else {
            out.insert(c, x);
        }
        used.push(rec);
    }
    Ok((out, json!({"endpoint_system": diag, "p14_rows": ep.len(), "masters_used": used})))
}

/// the lane's soft V (soft_V_euclid_nogamma keys) from its soft plan, its RustRed soft reduction and the known lower-loop masters
pub fn lane_v(lane_plan: &Value, lane_fam: &Value, lane_config: &Value, sol: &Value, lower: &BTreeMap<String, Float>, prec: u32) -> Result<(BTreeMap<String, String>, Value), String> {
    let e_slot = lane_config["e"].as_u64().ok_or("lane config has no e")? as usize;
    let d0: Rational = sol["fixed"]["d"].as_str().ok_or("no fixed d")?.parse().map_err(|e: rug::rational::ParseRationalError| e.to_string())?;
    let te = t_e(&d0, prec);
    let mom = softplan::slot_momenta(lane_fam);
    let allowed: BTreeSet<usize> = lane_plan["derived"]["sub"].as_array().ok_or("lane plan has no derived.sub")?.iter().map(|x| x.as_u64().unwrap() as usize).collect();
    let masters: Vec<Vec<i64>> = sol["masters"].as_array().ok_or("no masters")?.iter().map(vi).collect();
    let mut mval: HashMap<Vec<i64>, Float> = HashMap::new();
    let mut cls = Vec::new();
    for m in &masters {
        if m[e_slot] != 1 {
            return Err(format!("lane soft master {m:?} has e-slot power {} (expected the 1L tadpole, power 1)", m[e_slot]));
        }
        let (t, parts) = factorise(m, e_slot, &mom, &allowed).map_err(|e| format!("lane soft master: {e}"))?;
        let mut v = Float::with_val(prec, 1);
        for _ in 0..t {
            v *= &te;
        }
        for c in &parts {
            v *= lower.get(c).ok_or_else(|| format!("lane soft master {m:?}: component {} ({c}) is not determined by the lower system", label(c)))?;
        }
        // cross-check of an in-house declared classification (never used for the value)
        let key = serde_json::to_string(m).unwrap().replace(',', ", ");
        let declared = lane_config.get("soft_masters").and_then(|s| s.get(&key));
        let check = match declared {
            None => json!("none declared"),
            Some(k) => {
                let ty = k["type"].as_str().unwrap_or("?");
                let tp = k.get("T").and_then(|x| x.as_u64()).unwrap_or(0) as u32;
                let mine: Vec<String> = parts.iter().map(|c| label(c)).collect();
                let ok = tp == t && ((ty == "1" && mine.is_empty()) || mine == vec![ty.to_string()]);
                if !ok {
                    return Err(format!("lane soft master {key}: native factorisation T^{t} x {mine:?} disagrees with the config ({ty}, T^{tp})"));
                }
                json!("agrees")
            }
        };
        cls.push(json!({"master": m, "T_power": t, "components": parts.iter().map(|c| label(c)).collect::<Vec<_>>(), "config_check": check}));
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
        for (mu, c) in rules.get(x).ok_or_else(|| format!("no reduction rule for {x:?}"))? {
            let sgn = if (mu.iter().sum::<i64>() - x.iter().sum::<i64>()).rem_euclid(2) == 1 { -1 } else { 1 };
            tot += Float::with_val(prec, &mval[mu] * Float::with_val(prec, c)) * sgn;
        }
        Ok(tot)
    };
    let dg = (prec as f64 * 0.30103) as usize;
    let mut out = BTreeMap::new();
    for (key, poly) in lane_plan["soft_polys"].as_object().ok_or("no soft_polys")? {
        let mut val = Float::new(prec);
        for term in poly.as_array().unwrap() {
            let mut nu = vi(&term[0]);
            let dp = term[1].as_i64().unwrap();
            let mut c: Rational = term[2].as_str().unwrap().parse().map_err(|e: rug::rational::ParseRationalError| e.to_string())?;
            for _ in 0..dp {
                c /= &d0;
            }
            if nu[e_slot] != 0 {
                return Err("soft integral uses the eta slot".into());
            }
            nu[e_slot] = 1;
            val += Float::with_val(prec, c) * red(&nu)?;
        }
        out.insert(key.clone(), val.to_string_radix(10, Some(dg)));
    }
    let n = out.len();
    Ok((out, json!({"soft_masters": cls, "n_soft_polys": n})))
}

fn d0_of(sol: &Value) -> Result<Rational, String> {
    sol["fixed"]["d"].as_str().ok_or("no fixed d")?.parse().map_err(|e: rug::rational::ParseRationalError| e.to_string())
}

/// The lane's soft V from a CHAIN SPEC and a 4L sample (only its P14 rows are read):
/// {"lower": {"family": 5L family JSON, "class_config": config of the lower class (e, j, top_powers),
///            "p14_family": native P14 family JSON, "soft_solution": RustRed soft reduction of the lower class at d0},
///  "lane":  {"family": lane family JSON, "soft_solution": the lane's RustRed soft reduction at d0}}
/// Paths may contain "{d0}" (replaced by NUM_DEN of the lane's d0). Returns (V map, provenance).
pub fn lane_boundary(spec_path: &str, sample_path: &str, lane_cfg: &Value, lane_masters: &[Vec<i64>], d0: &Rational, prec: u32) -> Result<(BTreeMap<String, String>, Value), String> {
    let spec = rd(spec_path)?;
    let tag = format!("{}_{}", d0.numer(), d0.denom());
    let path = |sect: &str, k: &str| -> Result<String, String> {
        Ok(spec[sect][k].as_str().ok_or_else(|| format!("chain spec: {sect}.{k} missing"))?.replace("{d0}", &tag))
    };
    let sample = rd(sample_path)?;
    let eps: Rational = sample["epsilon"].as_str().ok_or("sample has no epsilon")?.parse().map_err(|e: rug::rational::ParseRationalError| e.to_string())?;
    if Rational::from(4 - Rational::from(2 * &eps)) != *d0 {
        return Err(format!("4L sample epsilon {eps} does not match d0 {d0}"));
    }
    let t0 = std::time::Instant::now();
    // 1. lower system
    let fam5 = rd(&path("lower", "family")?)?;
    let cfg_lo = rd(&path("lower", "class_config")?)?;
    let e_lo = cfg_lo["e"].as_u64().ok_or("lower class config has no e")? as usize;
    let p14t: Vec<Vec<i64>> = sample["endpoint"].as_array().ok_or("sample has no endpoints")?.iter()
        .filter(|x| x.get("family").and_then(|f| f.as_str()).unwrap_or("") == "").map(|x| vi(&x["target"])).collect();
    let sp = softplan::build(&fam5, &cfg_lo, &[], &rd(&path("lower", "p14_family")?)?, &p14t, None);
    let sol_lo = rd(&path("lower", "soft_solution")?)?;
    if d0_of(&sol_lo)? != *d0 {
        return Err("lower soft solution is at another d0".into());
    }
    let (lower, d_lo) = lower_values(&fam5, &sp.plan, &sol_lo, &sample, e_lo, prec)?;
    // 2. lane
    let lfam = rd(&path("lane", "family")?)?;
    let lplan = softplan::build_soft_only(&lfam, lane_cfg, lane_masters);
    let sol_l = rd(&path("lane", "soft_solution")?)?;
    if d0_of(&sol_l)? != *d0 {
        return Err("lane soft solution is at another d0".into());
    }
    let (v, d_l) = lane_v(&lplan, &lfam, lane_cfg, &sol_l, &lower, prec)?;
    let prov = json!({"V": "native lane chain (W3c lower.rs): P14 endpoint values -> lower-class soft system -> lower-loop masters by matroid class -> lane soft V",
        "chain_spec": spec_path, "spec": spec, "p14_sample": sample_path, "p14_rows": p14t.len(),
        "lower_values": lower.iter().map(|(k, x)| json!({"component": k, "label": label(k), "value": x.to_string_radix(10, Some(60))})).collect::<Vec<_>>(),
        "lower_diag": d_lo, "lane_diag": d_l, "seconds": t0.elapsed().as_secs_f64()});
    Ok((v, prov))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(x: &[i64]) -> Vq {
        x.iter().map(|&a| Rational::from(a)).collect()
    }

    #[test]
    fn matroid_classes_of_small_vacuum_graphs() {
        // sunset k1, k2, k1 + k2 (2 loops) and a relabelled / re-routed copy
        let s3 = canonical_matroid(&[v(&[1, 0]), v(&[0, 1]), v(&[1, 1])]).unwrap();
        let s3b = canonical_matroid(&[v(&[1, -1]), v(&[1, 0]), v(&[0, 1])]).unwrap();
        assert_eq!(s3, s3b);
        assert_eq!(label(&s3), "S3");
        // 3L banana and M5 = Pi(q)^2 / (q^2 + 1); the dotted banana (two lines in series, parallel momenta) has the same size and
        // rank as M5 but another matroid
        let b4 = canonical_matroid(&[v(&[1, 0, 0]), v(&[0, 1, 0]), v(&[0, 0, 1]), v(&[1, 1, 1])]).unwrap();
        assert_eq!(label(&b4), "B4");
        let m5 = canonical_matroid(&[v(&[1, 0, 0]), v(&[0, 1, 0]), v(&[1, -1, 0]), v(&[0, 0, 1]), v(&[1, 0, -1])]).unwrap();
        let m5b = canonical_matroid(&[v(&[0, 0, 1]), v(&[1, 0, -1]), v(&[1, 0, 0]), v(&[1, -1, 0]), v(&[0, 1, 0])]).unwrap();
        let dotted = canonical_matroid(&[v(&[1, 0, 0]), v(&[0, 1, 0]), v(&[0, 0, 1]), v(&[0, 0, 1]), v(&[1, 1, 1])]).unwrap();
        assert_eq!(m5, m5b);
        assert_ne!(m5, dotted);
        assert_eq!(label(&m5), "M5");
        // non-unimodular line set is refused
        assert!(canonical_matroid(&[v(&[1, 0]), v(&[0, 1]), v(&[1, 2])]).is_err());
    }

    #[test]
    fn factorisation_into_tadpoles_and_components() {
        // slots: 0 = eta line (excluded), 1..3 = sunset, 4, 5 = two independent tadpoles; loop momenta (a, b, c, d, e)
        let mom: Vec<Option<Vq>> = vec![
            Some(v(&[0, 0, 0, 0, 1])),
            Some(v(&[1, 0, 0, 0, 0])),
            Some(v(&[0, 1, 0, 0, 0])),
            Some(v(&[1, 1, 0, 0, 0])),
            Some(v(&[0, 0, 1, 0, 0])),
            Some(v(&[0, 0, 0, 1, 0])),
        ];
        let allowed: BTreeSet<usize> = (1..6).collect();
        let (t, parts) = factorise(&[1, 1, 1, 1, 1, 1], 0, &mom, &allowed).unwrap();
        assert_eq!(t, 2);
        assert_eq!(parts.iter().map(|c| label(c)).collect::<Vec<_>>(), vec!["S3"]);
        let (t, parts) = factorise(&[1, 1, 0, 0, 1, 1], 0, &mom, &allowed).unwrap();
        assert_eq!((t, parts.len()), (3, 0));
        assert!(factorise(&[1, 2, 1, 1, 0, 0], 0, &mom, &allowed).is_err());
        assert!(factorise(&[1, 1, 1, 1, 0, 0], 0, &mom, &BTreeSet::from([1, 2])).is_err());
    }
}
