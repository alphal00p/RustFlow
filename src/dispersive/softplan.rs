//! Native soft-boundary plan (port of values/d5l/v7a/tools/soft_gen.py `demand`): the exact soft-region expansion of every
//! e-block master into 4L integrals of the SUB family (eta line soft: p = q_e -> 0), and the exact images of the 4L sample's
//! endpoint targets (P14 lane, prism PR lane) in the same SUB family through unimodular loop-momentum maps.
//!
//! Euclidean, Q = q^2 + 1, measure d^Dk/pi^(D/2), no e^(gamma eps). Loop momentum k_j is eliminated (k_j = c_j (p - sum_{i != j} c_i k_i),
//! q_e = sum c_i k_i, |c_j| = 1); slot s: q_s = a_s(L) + b_s p; SUB = the 10 slots with b_s = 0 (complete 4L family in L). For a
//! master (e^a, rest R):
//!   V_{R,0} = int R(p = 0),
//!   V_{R,1} = int R(0) [ sum_i b_i^2 (-a_i/Q_i + 2 a_i (a_i+1) K_i^2/(D Q_i^2)) + sum_{i<j} 4 a_i a_j b_i b_j K_i.K_j/(D Q_i Q_j) ].
//! Polynomials map (index vector over 15 slots, power of 1/D) -> rational coefficient. Output: the plan JSON of soft_gen.py
//! (endpoint_polys, soft_polys) and the RustRed demand (integrals embedded with a q_e tadpole).
use rug::Rational;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};

pub type Poly = BTreeMap<(Vec<i64>, i64), Rational>;
type Lin = BTreeMap<Option<usize>, Rational>;
type Vq = Vec<Rational>;

fn q(x: i64) -> Rational {
    Rational::from(x)
}

fn rat(v: &Value) -> Rational {
    match v {
        Value::String(s) => s.parse().unwrap(),
        Value::Number(n) => Rational::from(n.as_i64().unwrap()),
        _ => panic!("not a rational {v}"),
    }
}

fn invert(m: &[Vec<Rational>]) -> Option<Vec<Vec<Rational>>> {
    let n = m.len();
    let mut a: Vec<Vec<Rational>> = m.iter().enumerate().map(|(i, r)| r.iter().cloned().chain((0..n).map(|j| q((i == j) as i64))).collect()).collect();
    for i in 0..n {
        let piv = (i..n).find(|&r| a[r][i] != 0)?;
        a.swap(i, piv);
        let p = a[i][i].clone();
        for x in a[i].iter_mut() {
            *x /= &p;
        }
        for r in 0..n {
            if r != i && a[r][i] != 0 {
                let f = a[r][i].clone();
                let ri = a[i].clone();
                for (x, y) in a[r].iter_mut().zip(&ri) {
                    *x -= Rational::from(&f * y);
                }
            }
        }
    }
    Some(a.into_iter().map(|r| r[n..].to_vec()).collect())
}

fn det(m: &[Vec<Rational>]) -> Rational {
    let mut m = m.to_vec();
    let n = m.len();
    let mut d = q(1);
    for i in 0..n {
        let Some(piv) = (i..n).find(|&r| m[r][i] != 0) else { return q(0) };
        if piv != i {
            m.swap(i, piv);
            d = -d;
        }
        d *= &m[i][i];
        for r in i + 1..n {
            let f = Rational::from(&m[r][i] / &m[i][i]);
            let mi = m[i].clone();
            for (x, y) in m[r].iter_mut().zip(&mi) {
                *x -= Rational::from(&f * y);
            }
        }
    }
    d
}

/// slot -> 5x5 symmetric G with q^2 = k^T G k
fn n_loops(fam: &Value) -> usize {
    fam["loop_momenta"].as_array().map(|v| v.len()).unwrap_or_else(|| {
        1 + fam["coordinates"].as_array().unwrap().iter().map(|c| c["right"].as_u64().unwrap() as usize).max().unwrap()
    })
}

fn quad_forms(fam: &Value) -> Vec<Vec<Vec<Rational>>> {
    let nl = n_loops(fam);
    let co: Vec<(usize, usize)> = fam["coordinates"].as_array().unwrap().iter().map(|c| (c["left"].as_u64().unwrap() as usize, c["right"].as_u64().unwrap() as usize)).collect();
    fam["denominators"]
        .as_array()
        .unwrap()
        .iter()
        .map(|d| {
            let mut g = vec![vec![q(0); nl]; nl];
            for ((a, b), x) in co.iter().zip(d["coefficients"].as_array().unwrap()) {
                let x = rat(x);
                if a == b {
                    g[*a][*a] += x;
                } else {
                    let h = x / 2;
                    g[*a][*b] += &h;
                    g[*b][*a] += h;
                }
            }
            g
        })
        .collect()
}

fn momentum_any(g: &[Vec<Rational>]) -> Option<Vq> {
    let nl = g.len();
    let lead = (0..nl).find(|&a| g[a][a] != 0)?;
    if g[lead][lead] != 1 {
        return None;
    }
    let v: Vq = (0..nl).map(|b| g[lead][b].clone()).collect();
    for a in 0..nl {
        for b in 0..nl {
            if Rational::from(&v[a] * &v[b]) != g[a][b] {
                return None;
            }
        }
    }
    Some(v)
}

/// rank-1 momenta of a 4-loop family given by scalar-product rows over PAIRS (P14 native / prism family); None for non-squares
fn momenta4(rows: &[(Vec<Rational>, bool)]) -> Vec<Option<Vq>> {
    let pairs: Vec<(usize, usize)> = (0..4).flat_map(|a| (a..4).map(move |b| (a, b))).collect();
    rows.iter()
        .map(|(c, square)| {
            if !square {
                return None;
            }
            let cm: BTreeMap<(usize, usize), Rational> = pairs.iter().cloned().zip(c.iter().cloned()).collect();
            let lead = (0..4).find(|&a| cm[&(a, a)] != 0).unwrap();
            let mut v = vec![q(0); 4];
            v[lead] = q(1);
            for b in 0..4 {
                if b != lead {
                    v[b] = Rational::from(&cm[&(lead.min(b), lead.max(b))] / 2);
                }
            }
            for a in 0..4 {
                for b in a..4 {
                    let want = if a == b { Rational::from(&v[a] * &v[a]) } else { Rational::from(&v[a] * &v[b]) * 2 };
                    assert!(want == cm[&(a, b)], "not a rank-1 propagator");
                }
            }
            Some(v)
        })
        .collect()
}

struct Cfg {
    nl: usize,
    ns: usize,
    e: usize,
    loops: Vec<usize>,
    sub: Vec<usize>,
    pdep: Vec<usize>,
    kvec: BTreeMap<usize, Vq>,
    bco: BTreeMap<usize, Rational>,
    pd: BTreeMap<usize, usize>,
    qm: Vec<Vq>,
    inv: BTreeMap<(usize, usize), BTreeMap<usize, Rational>>,
}

impl Cfg {
    fn new(g: &[Vec<Vec<Rational>>], e: usize, j: usize) -> Cfg {
        let qm: Vec<Vq> = g.iter().map(|x| momentum_any(x).expect("slot without a rank-1 momentum")).collect();
        let qe = qm[e].clone();
        let cj = qe[j].clone();
        assert!(cj == 1 || cj == -1);
        let (nl, ns) = (qe.len(), qm.len());
        let loops: Vec<usize> = (0..nl).filter(|&i| i != j).collect();
        let sub: Vec<usize> = (0..ns).filter(|&s| qm[s][j] == 0).collect();
        assert_eq!(sub.len(), (nl - 1) * nl / 2, "SUB is not a complete (L-1)-loop family");
        let (mut pdep, mut kvec, mut bco, mut pd) = (Vec::new(), BTreeMap::new(), BTreeMap::new(), BTreeMap::new());
        for s in 0..ns {
            if s == e || qm[s][j] == 0 {
                continue;
            }
            let gq = &qm[s];
            let b = Rational::from(&gq[j] * &cj);
            let mut a = vec![q(0); nl];
            for &i in &loops {
                a[i] = gq[i].clone() - Rational::from(&gq[j] * &cj) * &qe[i];
            }
            let neg: Vq = a.iter().map(|x| Rational::from(-x)).collect();
            if let Some(&t) = sub.iter().find(|&&t| qm[t] == a || qm[t] == neg) {
                pd.insert(s, t);
            }
            pdep.push(s);
            kvec.insert(s, a);
            bco.insert(s, b);
        }
        // scalar products of L as linear forms in (Q_s - 1), s in SUB
        let pairs: Vec<(usize, usize)> = loops.iter().enumerate().flat_map(|(i, &a)| loops[i..].iter().map(move |&b| (a, b))).collect();
        let m: Vec<Vec<Rational>> = sub.iter().map(|&s| pairs.iter().map(|&(a, b)| if a == b { g[s][a][b].clone() } else { Rational::from(&g[s][a][b] * 2) }).collect()).collect();
        let minv = invert(&m).expect("SUB family is not complete");
        let mut inv = BTreeMap::new();
        for (k, ab) in pairs.iter().enumerate() {
            let row: BTreeMap<usize, Rational> = sub.iter().enumerate().filter(|(si, _)| minv[k][*si] != 0).map(|(si, &s)| (s, minv[k][si].clone())).collect();
            inv.insert(*ab, row);
        }
        Cfg { nl, ns, e, loops, sub, pdep, kvec, bco, pd, qm, inv }
    }

    /// u.v for 5-vectors supported on the loops -> linear form in Q_slot (None = constant)
    fn sp(&self, u: &[Rational], v: &[Rational]) -> Lin {
        let mut out = Lin::new();
        for &a in &self.loops {
            for &b in &self.loops {
                let c = Rational::from(&u[a] * &v[b]);
                if c == 0 {
                    continue;
                }
                let ab = if self.loops.iter().position(|&x| x == a) <= self.loops.iter().position(|&x| x == b) { (a, b) } else { (b, a) };
                for (s, x) in &self.inv[&ab] {
                    let cx = Rational::from(&c * x);
                    *out.entry(Some(*s)).or_default() += &cx;
                    *out.entry(None).or_default() -= cx;
                }
            }
        }
        out.retain(|_, v| *v != 0);
        out
    }

    fn p_qpow(&self, poly: &Poly, s: usize, m: i64) -> Poly {
        if let Some(&t) = self.pd.get(&s) {
            return p_shift(poly, t, m);
        }
        assert!(m <= 0, "positive power on a p-dependent slot without a SUB image");
        let k = &self.kvec[&s];
        let mut lin = self.sp(k, k);
        *lin.entry(None).or_default() += 1;
        lin.retain(|_, v| *v != 0);
        let mut p = poly.clone();
        for _ in 0..(-m) {
            p = p_mul_lin(&p, &lin);
        }
        p
    }

    fn soft_terms(&self, master: &[i64], n: u32) -> Poly {
        let mut nu = vec![0i64; self.ns];
        for &s in &self.sub {
            nu[s] = master[s];
        }
        let base: Poly = BTreeMap::from([((nu, 0), q(1))]);
        let act: Vec<(usize, i64)> = self.pdep.iter().filter(|&&s| master[s] != 0).map(|&s| (s, master[s])).collect();
        if n == 0 {
            let mut b0 = base;
            for &(s, a) in &act {
                b0 = self.p_qpow(&b0, s, a);
            }
            return b0;
        }
        assert_eq!(n, 1);
        let mut tot = Poly::new();
        for &(s, a) in &act {
            let b2 = Rational::from(&self.bco[&s] * &self.bco[&s]);
            let mut rest = base.clone();
            for &(s2, a2) in &act {
                if s2 != s {
                    rest = self.p_qpow(&rest, s2, a2);
                }
            }
            tot = p_add(&tot, &self.p_qpow(&rest, s, a + 1), &(Rational::from(-a) * &b2));
            if a * (a + 1) != 0 {
                let k = &self.kvec[&s];
                let t = p_mul_lin(&self.p_qpow(&rest, s, a + 2), &self.sp(k, k));
                tot = p_add(&tot, &p_dpow(&t, 1), &(Rational::from(2 * a * (a + 1)) * &b2));
            }
        }
        for x in 0..act.len() {
            for y in x + 1..act.len() {
                let ((s1, a1), (s2, a2)) = (act[x], act[y]);
                let mut rest = base.clone();
                for &(s3, a3) in &act {
                    if s3 != s1 && s3 != s2 {
                        rest = self.p_qpow(&rest, s3, a3);
                    }
                }
                let t = self.p_qpow(&self.p_qpow(&rest, s1, a1 + 1), s2, a2 + 1);
                let t = p_mul_lin(&t, &self.sp(&self.kvec[&s1], &self.kvec[&s2]));
                let f = Rational::from(4 * a1 * a2) * &self.bco[&s1] * &self.bco[&s2];
                tot = p_add(&tot, &p_dpow(&t, 1), &f);
            }
        }
        tot
    }

    fn embed(&self, nu: &[i64]) -> Vec<i64> {
        let mut x = nu.to_vec();
        assert!(self.pdep.iter().chain([&self.e]).all(|&s| x[s] == 0));
        x[self.e] = 1;
        x
    }
}

fn p_mul_lin(poly: &Poly, lin: &Lin) -> Poly {
    let mut out = Poly::new();
    for ((nu, dp), c) in poly {
        for (s, x) in lin {
            let key = match s {
                None => (nu.clone(), *dp),
                Some(s) => {
                    let mut l = nu.clone();
                    l[*s] -= 1;
                    (l, *dp)
                }
            };
            *out.entry(key).or_default() += Rational::from(c * x);
        }
    }
    out.retain(|_, v| *v != 0);
    out
}

fn p_shift(poly: &Poly, slot: usize, k: i64) -> Poly {
    poly.iter()
        .map(|((nu, dp), c)| {
            let mut l = nu.clone();
            l[slot] += k;
            ((l, *dp), c.clone())
        })
        .collect()
}

fn p_add(a: &Poly, b: &Poly, f: &Rational) -> Poly {
    let mut out = a.clone();
    for (k, v) in b {
        *out.entry(k.clone()).or_default() += Rational::from(f * v);
    }
    out.retain(|_, v| *v != 0);
    out
}

fn p_dpow(poly: &Poly, k: i64) -> Poly {
    poly.iter().map(|((nu, dp), c)| ((nu.clone(), dp + k), c.clone())).collect()
}

/// A (4x4, |det| 1) with lane momenta K = A L, and the propagator permutation lane -> W4 line
fn find_map(nat_q: &[Option<Vq>], sub_q: &[Vq]) -> (Vec<Vec<Rational>>, BTreeMap<usize, usize>) {
    let idx: Vec<usize> = (0..nat_q.len()).filter(|&i| nat_q[i].is_some()).collect();
    let mut basis = None;
    'b: for a in 0..idx.len() {
        for b in a + 1..idx.len() {
            for c in b + 1..idx.len() {
                for d in c + 1..idx.len() {
                    let cand = [idx[a], idx[b], idx[c], idx[d]];
                    let m: Vec<Vec<Rational>> = cand.iter().map(|&i| nat_q[i].clone().unwrap()).collect();
                    if det(&m) != 0 {
                        basis = Some(cand);
                        break 'b;
                    }
                }
            }
        }
    }
    let basis = basis.expect("no momentum basis");
    let bt: Vec<Vec<Rational>> = (0..4).map(|r| basis.iter().map(|&i| nat_q[i].as_ref().unwrap()[r].clone()).collect()).collect();
    let binv = invert(&bt).unwrap();
    let n = sub_q.len();
    // permutations of 4 out of n in lexicographic order (as itertools.permutations), signs in product order (1, -1)
    let mut perms = Vec::new();
    for a in 0..n {
        for b in 0..n {
            for c in 0..n {
                for d in 0..n {
                    let v = [a, b, c, d];
                    if (0..4).all(|x| (x + 1..4).all(|y| v[x] != v[y])) {
                        perms.push(v);
                    }
                }
            }
        }
    }
    for img in perms {
        for sm in 0..16 {
            let signs: Vec<i64> = (0..4).map(|k| if (sm >> (3 - k)) & 1 == 1 { -1 } else { 1 }).collect();
            let s: Vec<Vec<Rational>> = (0..4).map(|r| (0..4).map(|k| Rational::from(&sub_q[img[k]][r] * signs[k])).collect()).collect();
            let at: Vec<Vec<Rational>> = (0..4).map(|r| (0..4).map(|c| (0..4).fold(q(0), |acc, k| acc + Rational::from(&s[r][k] * &binv[k][c]))).collect()).collect();
            let dd = det(&at);
            if dd != 1 && dd != -1 {
                continue;
            }
            let mut perm = BTreeMap::new();
            let mut used = BTreeSet::new();
            let mut ok = true;
            for &i in &idx {
                let nq = nat_q[i].as_ref().unwrap();
                let v: Vq = (0..4).map(|r| (0..4).fold(q(0), |acc, c| acc + Rational::from(&at[r][c] * &nq[c]))).collect();
                let neg: Vq = v.iter().map(|x| Rational::from(-x)).collect();
                match (0..n).find(|j| !used.contains(j) && (sub_q[*j] == v || sub_q[*j] == neg)) {
                    Some(h) => {
                        used.insert(h);
                        perm.insert(i, h);
                    }
                    None => {
                        ok = false;
                        break;
                    }
                }
            }
            if ok {
                let a: Vec<Vec<Rational>> = (0..4).map(|r| (0..4).map(|c| at[c][r].clone()).collect()).collect();
                return (a, perm);
            }
        }
    }
    panic!("no lane map");
}

pub struct SoftPlan {
    pub plan: Value,
    pub demand: Value,
}

/// inputs: 5L family JSON, class config, e-block masters, P14 native family + its endpoint targets, prism family + targets (optional)
#[allow(clippy::too_many_arguments)]
pub fn build(fam: &Value, config: &Value, emast: &[Vec<i64>], p14: &Value, p14_targets: &[Vec<i64>], prism: Option<(&Value, &[Vec<i64>])>) -> SoftPlan {
    let g = quad_forms(fam);
    let e = config["e"].as_u64().unwrap() as usize;
    let j = config["j"].as_u64().unwrap() as usize;
    let cfg = Cfg::new(&g, e, j);
    let top: Vec<i64> = config["top_powers"].as_array().unwrap().iter().map(|x| x.as_i64().unwrap()).collect();
    let mut w4: BTreeSet<usize> = cfg.sub.iter().filter(|&&s| top[s] > 0).cloned().collect();
    for &s in &cfg.pdep {
        if top[s] > 0 {
            w4.insert(*cfg.pd.get(&s).expect("p-dependent top line without SUB image"));
        }
    }
    if let Some(x) = config.get("w4_extra").and_then(|v| v.as_u64()) {
        w4.insert(x as usize);
    }
    let w4: Vec<usize> = w4.into_iter().collect();
    // a 9-line soft graph is the prism lane (soft_gen.py: config 'prism' defaults to the PR endpoint targets)
    assert!(w4.len() == 8 || w4.len() == 9, "soft graph is not an 8- or 9-line 4L graph: {w4:?}");
    let prism = if w4.len() == 9 { Some(prism.expect("9-line soft graph needs the prism lane inputs")) } else { None };
    let sub_q: Vec<Vq> = w4.iter().map(|&s| cfg.loops.iter().map(|&l| cfg.qm[s][l].clone()).collect()).collect();
    // P14 lane
    let p14_rows: Vec<(Vec<Rational>, bool)> = p14["denominators"].as_array().unwrap().iter().map(|d| (d["scalar_products"].as_array().unwrap().iter().map(rat).collect(), d["constant"] == "-1")).collect();
    let p14_pairs: Vec<Vec<Rational>> = p14_rows.iter().map(|r| r.0.clone()).collect();
    let nq = momenta4(&p14_rows);
    let (a14, perm14) = find_map(&nq, &sub_q);
    let lift = |a: &[Vec<Rational>], row: usize| -> Vq {
        let mut u = vec![q(0); cfg.nl];
        for l in 0..4 {
            u[cfg.loops[l]] = a[row][l].clone();
        }
        u
    };
    let pairs4: Vec<(usize, usize)> = (0..4).flat_map(|a| (a..4).map(move |b| (a, b))).collect();
    let mut ep: BTreeMap<String, Poly> = BTreeMap::new();
    for t in p14_targets {
        let mut nu = vec![0i64; cfg.ns];
        let mut isps = Vec::new();
        for (i, &ti) in t.iter().enumerate() {
            if ti == 0 {
                continue;
            }
            if let Some(&h) = perm14.get(&i) {
                nu[w4[h]] += ti;
            } else {
                assert!(ti < 0);
                isps.push((i, -ti));
            }
        }
        let mut poly: Poly = BTreeMap::from([((nu, 0), q(1))]);
        for (i, k) in isps {
            let nz: Vec<(usize, usize)> = pairs4.iter().zip(&p14_pairs[i]).filter(|(_, x)| **x != 0).map(|(p, _)| *p).collect();
            assert!(nz.len() == 1 && nz[0].0 != nz[0].1);
            let (a, b) = nz[0];
            let lin = cfg.sp(&lift(&a14, a), &lift(&a14, b));
            for _ in 0..k {
                poly = p_mul_lin(&poly, &lin);
            }
        }
        ep.insert(serde_json::to_string(t).unwrap().replace(',', ", "), poly);
    }
    if let Some((prf, prt)) = prism {
        let rows: Vec<(Vec<Rational>, bool)> = prf["denominators"].as_array().unwrap().iter().map(|d| (d["coefficients"].as_array().unwrap().iter().map(rat).collect(), true)).collect();
        let mut prq = momenta4(&rows);
        let q9 = prq[9].clone().unwrap();
        prq[9] = None;
        let (apr, permpr) = find_map(&prq, &sub_q);
        for t in prt {
            let mut nu = vec![0i64; cfg.ns];
            for i in 0..9 {
                if t[i] != 0 {
                    nu[w4[permpr[&i]]] += t[i];
                }
            }
            let mut poly: Poly = BTreeMap::from([((nu, 0), q(1))]);
            if t[9] != 0 {
                assert!(t[9] < 0);
                let mut kk = vec![q(0); cfg.nl];
                for l in 0..4 {
                    kk[cfg.loops[l]] = (0..4).fold(q(0), |acc, a| acc + Rational::from(&q9[a] * &apr[a][l]));
                }
                let mut lin = cfg.sp(&kk, &kk);
                *lin.entry(None).or_default() += 1;
                lin.retain(|_, v| *v != 0);
                for _ in 0..(-t[9]) {
                    poly = p_mul_lin(&poly, &lin);
                }
            }
            ep.insert(format!("PR{}", serde_json::to_string(t).unwrap().replace(',', ", ")), poly);
        }
    }
    let n1 = config.get("n1").cloned().unwrap_or(json!("a2"));
    let mut soft: BTreeMap<String, Poly> = BTreeMap::new();
    for m in emast {
        let key = serde_json::to_string(m).unwrap().replace(',', ", ");
        soft.insert(format!("{key}|0"), cfg.soft_terms(m, 0));
        let take1 = match &n1 {
            Value::String(s) if s == "all" => true,
            Value::String(s) if s == "a2" => m[e] >= 2,
            Value::Array(list) => list.iter().any(|x| x.as_array().map(|v| v.iter().map(|y| y.as_i64().unwrap()).collect::<Vec<_>>()) == Some(m.clone())),
            _ => false,
        };
        if take1 {
            soft.insert(format!("{key}|1"), cfg.soft_terms(m, 1));
        }
    }
    let mut ints = BTreeSet::new();
    for poly in ep.values().chain(soft.values()) {
        for (nu, _) in poly.keys() {
            ints.insert(cfg.embed(nu));
        }
    }
    let ser = |p: &Poly| -> Value { json!(p.iter().map(|((nu, dp), c)| json!([nu, dp, c.to_string()])).collect::<Vec<_>>()) };
    let plan = json!({"config": config, "derived": {"loops": cfg.loops, "sub": cfg.sub, "pdep": cfg.pdep, "w4lines": w4},
        "A_P14_from_L": a14.iter().map(|r| r.iter().map(|x| x.to_string()).collect::<Vec<_>>()).collect::<Vec<_>>(),
        "perm_P14_to_W4line": perm14.iter().map(|(k, v)| (k.to_string(), json!(v))).collect::<serde_json::Map<_, _>>(),
        "endpoint_polys": ep.iter().map(|(k, v)| (k.clone(), ser(v))).collect::<serde_json::Map<_, _>>(),
        "soft_polys": soft.iter().map(|(k, v)| (k.clone(), ser(v))).collect::<serde_json::Map<_, _>>(),
        "emast": emast, "n_integrals": ints.len()});
    let topmask: u64 = w4.iter().chain([&e]).map(|&s| 1u64 << s).sum();
    let demand = json!({"schema": "c2.t2-demand.v1", "job": format!("w3a-{}-soft", config["class"]), "class_rep": config["class"],
        "slot_top": topmask, "targets": ints.iter().collect::<Vec<_>>(), "expected": {}});
    SoftPlan { plan, demand }
}

/// soft polynomials only (soft4.py demand: the soft family's masters are known lower-loop values, no endpoint map)
pub fn build_soft_only(fam: &Value, config: &Value, emast: &[Vec<i64>]) -> Value {
    let g = quad_forms(fam);
    let e = config["e"].as_u64().unwrap() as usize;
    let j = config["j"].as_u64().unwrap() as usize;
    let cfg = Cfg::new(&g, e, j);
    let n1 = config.get("n1").cloned().unwrap_or(json!("a2"));
    let mut soft = serde_json::Map::new();
    let ser = |p: &Poly| -> Value { json!(p.iter().map(|((nu, dp), c)| json!([nu, dp, c.to_string()])).collect::<Vec<_>>()) };
    for m in emast {
        let key = serde_json::to_string(m).unwrap().replace(',', ", ");
        soft.insert(format!("{key}|0"), ser(&cfg.soft_terms(m, 0)));
        let take1 = match &n1 {
            Value::String(s) if s == "all" => true,
            Value::String(s) if s == "a2" => m[e] >= 2,
            _ => false,
        };
        if take1 {
            soft.insert(format!("{key}|1"), ser(&cfg.soft_terms(m, 1)));
        }
    }
    json!({"config": config, "derived": {"loops": cfg.loops, "sub": cfg.sub, "pdep": cfg.pdep}, "endpoint_polys": {}, "soft_polys": soft})
}

/// matroid components of a set of vectors (exact): list of (size, rank) sorted, via fundamental circuits of a greedy basis
pub fn components(vecs: &[Vq]) -> Vec<(usize, usize)> {
    let n = vecs.len();
    // greedy basis with exact elimination
    let mut basis: Vec<usize> = Vec::new();
    let mut ech: Vec<(usize, Vq)> = Vec::new(); // (pivot column, reduced row) of the basis vectors
    let reduce = |ech: &[(usize, Vq)], v: &Vq| -> Vq {
        let mut v = v.clone();
        for (c, r) in ech {
            if v[*c] != 0 {
                let f = Rational::from(&v[*c] / &r[*c]);
                for (x, y) in v.iter_mut().zip(r) {
                    *x -= Rational::from(&f * y);
                }
            }
        }
        v
    };
    for (i, v) in vecs.iter().enumerate() {
        let r = reduce(&ech, v);
        if let Some(c) = r.iter().position(|x| *x != 0) {
            ech.push((c, r));
            basis.push(i);
        }
    }
    let rank = basis.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut Vec<usize>, x: usize) -> usize {
        let mut x = x;
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    // fundamental circuit of a non-basis vector: the basis vectors with a nonzero coefficient in its expansion
    for i in 0..n {
        if basis.contains(&i) {
            continue;
        }
        // solve sum_k x_k b_k = v by exact Gauss on the rank x dim system
        let dimv = vecs[i].len();
        let mut a: Vec<Vec<Rational>> = (0..dimv).map(|r| basis.iter().map(|&b| vecs[b][r].clone()).chain([vecs[i][r].clone()]).collect()).collect();
        let mut row = 0;
        let mut pivc = Vec::new();
        for c in 0..rank {
            let Some(p) = (row..dimv).find(|&r| a[r][c] != 0) else { continue };
            a.swap(row, p);
            let pv = a[row][c].clone();
            for x in a[row].iter_mut() {
                *x /= &pv;
            }
            for r in 0..dimv {
                if r != row && a[r][c] != 0 {
                    let f = a[r][c].clone();
                    let rr = a[row].clone();
                    for (x, y) in a[r].iter_mut().zip(&rr) {
                        *x -= Rational::from(&f * y);
                    }
                }
            }
            pivc.push(c);
            row += 1;
        }
        for (r, &c) in pivc.iter().enumerate() {
            if a[r][rank] != 0 {
                let (x, y) = (find(&mut parent, i), find(&mut parent, basis[c]));
                parent[x] = y;
            }
        }
    }
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push(i);
    }
    let mut out: Vec<(usize, usize)> = groups
        .values()
        .map(|g| {
            let mut e2: Vec<(usize, Vq)> = Vec::new();
            for &i in g {
                let r = reduce(&e2, &vecs[i]);
                if let Some(c) = r.iter().position(|x| *x != 0) {
                    e2.push((c, r));
                }
            }
            (g.len(), e2.len())
        })
        .collect();
    out.sort();
    out
}

/// the rank-1 momentum (over the loop momenta) of every slot of a family
pub fn slot_momenta(fam: &Value) -> Vec<Option<Vq>> {
    quad_forms(fam).iter().map(|g| momentum_any(g)).collect()
}

/// exact comparison of two plans' polynomial maps: (number of keys compared, mismatching keys)
pub fn compare(a: &Value, b: &Value) -> (usize, Vec<String>) {
    let norm = |v: &Value| -> BTreeMap<(Vec<i64>, i64), Rational> {
        v.as_array().unwrap().iter().map(|t| ((t[0].as_array().unwrap().iter().map(|x| x.as_i64().unwrap()).collect(), t[1].as_i64().unwrap()), rat(&t[2]))).collect()
    };
    let mut n = 0;
    let mut bad = Vec::new();
    for sect in ["endpoint_polys", "soft_polys"] {
        let empty = serde_json::Map::new();
        let (ma, mb) = (a[sect].as_object().unwrap_or(&empty), b[sect].as_object().unwrap_or(&empty));
        for k in ma.keys().chain(mb.keys()).collect::<BTreeSet<_>>() {
            n += 1;
            match (ma.get(k), mb.get(k)) {
                (Some(x), Some(y)) if norm(x) == norm(y) => {}
                _ => bad.push(format!("{sect}:{k}")),
            }
        }
    }
    (n, bad)
}
