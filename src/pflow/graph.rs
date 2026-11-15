//! Graph-level region data of a one-mass massless propagator family (E-normalisation: Euclidean, measure d^dk / pi^(d/2) per
//! loop, p^2 = 1). Line indices are exact rationals (h = d0/2 enters through one-loop insertions); values are MPFR.
//!
//! Regions of J(z) = family integral with 1/(k_e^2 + z)^(a_e) on line e:
//!  * z -> 0, soft set S (contains e): hard graph H = G \ S (a p-integral), soft vacuum graph = S with H's vertices identified;
//!    S \ e is a massless two-point function C (k_e^2)^-alpha, so the soft part is C Vm(a_e, alpha) z^(h - a_e - alpha).
//!    Exponent class: L_s h (L_s = loops of the soft graph). Valid iff H connects the external vertices, H + (external edge) is
//!    bridgeless and the soft graph is bridgeless.
//!  * z -> infinity, hard subgraph gamma (contains e, connected, bridgeless): gamma \ e at zero external momenta gives
//!    C Vm(a_e, alpha), times the co-graph G / gamma (a p-integral). Class L_gamma h.
//! Only the LOWEST order of each class is a clean boundary row (higher orders mix in subleading terms of other regions);
//! `exact_series` gives all orders for the simple masters (pure power; one-loop massive bubble after massless reductions).
//! The massless two-point evaluator uses series / parallel (G functions) / p-carrying bridges and the IBP triangle rule
//! (integer indices); a graph it cannot reduce is looked up in the `Provided` table (values from a lower-loop run of the same
//! mode, keyed by `canonical`), otherwise reported (`None`), never guessed.
use rug::{Float, Rational};
use std::collections::HashMap;

/// massless two-point values supplied from outside (E-normalisation, P^2 = 1), keyed by `canonical`
#[derive(Debug, Default)]
pub struct Provided {
    pub values: HashMap<String, Float>,
    /// canonical keys of irreducible graphs that were not in the table (diagnostic)
    pub missing: std::sync::Mutex<Vec<String>>,
}

impl Clone for Provided {
    fn clone(&self) -> Provided {
        Provided { values: self.values.clone(), missing: std::sync::Mutex::new(self.missing.lock().unwrap().clone()) }
    }
}

impl Provided {
    /// JSON list of {"edges": [[u, v, a], ...], "ext": [s, t], "value": "..."} (a integer or "p/q")
    pub fn from_json(v: &serde_json::Value, prec: u32) -> Result<Provided, String> {
        let mut values = HashMap::new();
        for x in v.as_array().ok_or("provided: a list")? {
            let edges: Vec<(usize, usize, Rational)> = x["edges"].as_array().ok_or("provided: edges")?.iter().map(|e| {
                let a = match &e[2] { serde_json::Value::String(s) => s.parse::<Rational>().map_err(|e| e.to_string()), n => Ok(Rational::from(n.as_i64().unwrap_or(0))) };
                Ok((e[0].as_u64().ok_or("u")? as usize, e[1].as_u64().ok_or("v")? as usize, a?))
            }).collect::<Result<_, String>>()?;
            let ext = (x["ext"][0].as_u64().ok_or("ext")? as usize, x["ext"][1].as_u64().ok_or("ext")? as usize);
            let val = Float::with_val(prec, Float::parse(x["value"].as_str().ok_or("provided: value string")?).map_err(|e| e.to_string())?);
            values.insert(canonical(&edges, ext.0, ext.1), val);
        }
        Ok(Provided { values, missing: Default::default() })
    }
}

fn connected_plain(edges: &[(usize, usize, Rational)], a: usize, b: usize, skip: usize) -> bool {
    let mut seen = vec![a];
    let mut st = vec![a];
    while let Some(x) = st.pop() {
        for (i, e) in edges.iter().enumerate() {
            if i == skip {
                continue;
            }
            let y = if e.0 == x { e.1 } else if e.1 == x { e.0 } else { continue };
            if !seen.contains(&y) {
                seen.push(y);
                st.push(y);
            }
        }
    }
    seen.contains(&b)
}

/// canonical key of a massless two-point graph (same as s1/proto/prov.py): contract index-0 lines; repeat {contract a bridge
/// (it carries the full external momentum: factor 1 at P^2 = 1); merge two lines in series at an internal degree-2 vertex}; then
/// the lexicographic minimum over relabellings (externals -> 0, 1 in both orders, internal vertices -> 2..) of the sorted list
/// (min(u,v), max(u,v), index), printed as JSON
pub fn canonical(edges: &[(usize, usize, Rational)], s: usize, t: usize) -> String {
    let (mut s, mut t) = (s, t);
    let mut e: Vec<(usize, usize, Rational)> = edges.to_vec();
    let relabel = |e: &mut Vec<(usize, usize, Rational)>, from: usize, to: usize| {
        for x in e.iter_mut() {
            if x.0 == from { x.0 = to; }
            if x.1 == from { x.1 = to; }
        }
    };
    while let Some(i) = e.iter().position(|x| x.2 == 0) {
        let (u, v, _) = e.remove(i);
        if u != v {
            relabel(&mut e, v, u);
            if s == v { s = u; }
            if t == v { t = u; }
        }
    }
    loop {
        if let Some(i) = (0..e.len()).find(|&i| !connected_plain(&e, e[i].0, e[i].1, i)) {
            let (u, v, _) = e.remove(i);
            relabel(&mut e, v, u);
            if s == v { s = u; }
            if t == v { t = u; }
            continue;
        }
        let mut deg: Vec<(usize, usize)> = Vec::new();
        for x in &e {
            for y in [x.0, x.1] {
                match deg.iter_mut().find(|d| d.0 == y) { Some(d) => d.1 += 1, None => deg.push((y, 1)) }
            }
        }
        // python iterates a Counter in insertion order; the merge order does not change the result
        let Some(&(v, _)) = deg.iter().find(|d| d.1 == 2 && d.0 != s && d.0 != t && !e.iter().any(|x| x.0 == d.0 && x.1 == d.0)) else { break };
        let idx: Vec<usize> = (0..e.len()).filter(|&i| e[i].0 == v || e[i].1 == v).collect();
        let (e1, e2) = (e[idx[0]].clone(), e[idx[1]].clone());
        let u = if e1.1 == v { e1.0 } else { e1.1 };
        let w = if e2.1 == v { e2.0 } else { e2.1 };
        e.remove(idx[1]);
        e.remove(idx[0]);
        e.push((u, w, Rational::from(&e1.2 + &e2.2)));
    }
    let mut verts: Vec<usize> = e.iter().flat_map(|x| [x.0, x.1]).filter(|&x| x != s && x != t).collect();
    verts.sort();
    verts.dedup();
    let n = verts.len();
    let mut best: Option<Vec<(usize, usize, String)>> = None;
    let mut perm: Vec<usize> = (0..n).collect();
    let go = |perm: &[usize], best: &mut Option<Vec<(usize, usize, String)>>| {
        for (a, b) in [(s, t), (t, s)] {
            let map = |x: usize| if x == a { 0 } else if x == b { 1 } else { 2 + perm[verts.iter().position(|&y| y == x).unwrap()] };
            let mut key: Vec<(usize, usize, String)> = e.iter().map(|x| { let (p, q) = (map(x.0), map(x.1)); (p.min(q), p.max(q), x.2.to_string()) }).collect();
            key.sort();
            if best.as_ref().map_or(true, |b| key < *b) { *best = Some(key); }
        }
    };
    // Heap's algorithm over all relabellings of the internal vertices
    let mut c = vec![0usize; n];
    go(&perm, &mut best);
    let mut i = 0;
    while i < n {
        if c[i] < i {
            if i % 2 == 0 { perm.swap(0, i) } else { perm.swap(c[i], i) }
            go(&perm, &mut best);
            c[i] += 1;
            i = 0;
        } else {
            c[i] = 0;
            i += 1;
        }
    }
    let best = best.unwrap_or_default();
    format!("[{}]", best.iter().map(|(u, v, a)| format!("[{u}, {v}, \"{a}\"]")).collect::<Vec<_>>().join(", "))
}

pub type Edge = (usize, usize);

#[derive(Clone, Debug)]
pub struct Graph {
    pub edges: Vec<Edge>,
    pub ext: (usize, usize),
    pub loops: usize,
}

/// 1/Gamma(x) for an exact rational x (zero at non-positive integers)
pub fn rgamma(x: &Rational, prec: u32) -> Float {
    if x.denom() == &1 && *x <= 0 {
        return Float::new(prec);
    }
    let g = Float::with_val(prec, x).gamma();
    Float::with_val(prec, 1 / g)
}
pub fn gamma(x: &Rational, prec: u32) -> Option<Float> {
    if x.denom() == &1 && *x <= 0 {
        return None;
    }
    Some(Float::with_val(prec, x).gamma())
}
fn rat_i64(x: &Rational) -> Option<i64> {
    if x.denom() == &1 { x.numer().to_i64() } else { None }
}
fn q(n: i64) -> Rational {
    Rational::from(n)
}
/// massless one-loop bubble int 1/(k^2)^a ((k-p)^2)^b, p^2 = 1
pub fn g_bubble(a: &Rational, b: &Rational, h: &Rational, prec: u32) -> Option<Float> {
    let g1 = gamma(&(Rational::from(a + b) - h.clone()), prec)?;
    let g2 = gamma(&Rational::from(h - a), prec)?;
    let g3 = gamma(&Rational::from(h - b), prec)?;
    let mut v = Float::with_val(prec, &g1 * &g2);
    v *= &g3;
    v *= rgamma(a, prec);
    v *= rgamma(b, prec);
    v *= rgamma(&(Rational::from(h * 2u32) - a.clone() - b.clone()), prec);
    Some(v)
}
/// int (k^2)^-al (k^2 + 1)^-a  (coefficient of (m^2)^(h - a - al))
pub fn vm(a: &Rational, al: &Rational, h: &Rational, prec: u32) -> Option<Float> {
    let g1 = gamma(&Rational::from(h - al), prec)?;
    let g2 = gamma(&(Rational::from(a + al) - h.clone()), prec)?;
    let mut v = Float::with_val(prec, &g1 * &g2);
    v *= rgamma(h, prec);
    v *= rgamma(a, prec);
    Some(v)
}
pub fn tad(a: &Rational, h: &Rational, prec: u32) -> Option<Float> {
    let g = gamma(&(a.clone() - h.clone()), prec)?;
    Some(Float::with_val(prec, g * rgamma(a, prec)))
}

#[derive(Clone, Debug)]
struct E2 {
    u: usize,
    v: usize,
    a: Rational,
}

/// outcome of a massless evaluation
enum Ev {
    Zero,
    Val(Float, Rational),
}

fn connected(edges: &[E2], a: usize, b: usize, skip: Option<usize>) -> bool {
    if a == b {
        return true;
    }
    let mut seen = vec![a];
    let mut st = vec![a];
    while let Some(x) = st.pop() {
        for (i, e) in edges.iter().enumerate() {
            if Some(i) == skip {
                continue;
            }
            let y = if e.u == x {
                e.v
            } else if e.v == x {
                e.u
            } else {
                continue;
            };
            if !seen.contains(&y) {
                seen.push(y);
                st.push(y);
            }
        }
    }
    seen.contains(&b)
}

/// massless two-point function between s and t: value = C (P^2)^-alpha. None = not reducible by the implemented rules.
fn eval2pt(mut edges: Vec<E2>, mut s: usize, mut t: usize, h: &Rational, prec: u32, depth: usize, pv: &Provided) -> Option<Ev> {
    if depth > 400 {
        return None;
    }
    // index-0 lines are absent propagators: contract
    while let Some(i) = edges.iter().position(|e| e.a == 0) {
        let e = edges.remove(i);
        if e.u != e.v {
            for x in edges.iter_mut() {
                if x.u == e.v {
                    x.u = e.u;
                }
                if x.v == e.v {
                    x.v = e.u;
                }
            }
            if s == e.v {
                s = e.u;
            }
            if t == e.v {
                t = e.u;
            }
        }
    }
    if s == t {
        return Some(if edges.is_empty() { Ev::Val(Float::with_val(prec, 1), q(0)) } else { Ev::Zero });
    }
    if !pv.values.is_empty() && edges.iter().all(|e| e.a.denom() == &1) {
        let ce: Vec<(usize, usize, Rational)> = edges.iter().map(|e| (e.u, e.v, e.a.clone())).collect();
        if let Some(v) = pv.values.get(&canonical(&ce, s, t)) {
            // C (P^2)^-alpha with alpha = sum a - L h
            let mut vs: Vec<usize> = edges.iter().flat_map(|e| [e.u, e.v]).collect();
            vs.push(s);
            vs.push(t);
            vs.sort();
            vs.dedup();
            let l = edges.len() as i64 - vs.len() as i64 + 1;
            let mut alpha = q(0);
            for e in &edges {
                alpha += &e.a;
            }
            alpha -= Rational::from(h * l);
            return Some(Ev::Val(Float::with_val(prec, v), alpha));
        }
    }
    let mut c = Float::with_val(prec, 1);
    let mut alpha = q(0);
    loop {
        if edges.iter().any(|e| e.u == e.v) {
            return Some(Ev::Zero);
        }
        // parallel
        let mut did = false;
        'par: for i in 0..edges.len() {
            for j in 0..i {
                let (a, b) = (&edges[i], &edges[j]);
                if (a.u == b.u && a.v == b.v) || (a.u == b.v && a.v == b.u) {
                    let gv = g_bubble(&b.a, &a.a, h, prec)?;
                    c *= gv;
                    let na = Rational::from(&b.a + &a.a) - h.clone();
                    edges[j].a = na;
                    edges.remove(i);
                    did = true;
                    break 'par;
                }
            }
        }
        if did {
            continue;
        }
        // degrees
        let mut deg = std::collections::BTreeMap::new();
        for e in &edges {
            *deg.entry(e.u).or_insert(0usize) += 1;
            *deg.entry(e.v).or_insert(0usize) += 1;
        }
        if deg.iter().any(|(v, d)| *d == 1 && *v != s && *v != t) {
            return Some(Ev::Zero);
        }
        if let Some((&v, _)) = deg.iter().find(|(v, d)| **d == 2 && **v != s && **v != t) {
            let idx: Vec<usize> = (0..edges.len()).filter(|&i| edges[i].u == v || edges[i].v == v).collect();
            let (e1, e2) = (edges[idx[0]].clone(), edges[idx[1]].clone());
            let u = if e1.v == v { e1.u } else { e1.v };
            let w = if e2.v == v { e2.u } else { e2.v };
            edges.remove(idx[1]);
            edges.remove(idx[0]);
            edges.push(E2 { u, v: w, a: Rational::from(&e1.a + &e2.a) });
            continue;
        }
        // bridges
        let mut bridge = None;
        for i in 0..edges.len() {
            if !connected(&edges, edges[i].u, edges[i].v, Some(i)) {
                bridge = Some(i);
                break;
            }
        }
        if let Some(i) = bridge {
            let rest: Vec<E2> = edges.iter().enumerate().filter(|(k, _)| *k != i).map(|(_, e)| e.clone()).collect();
            if connected(&rest, s, t, None) {
                return Some(Ev::Zero); // a zero-momentum line
            }
            let e = edges.remove(i);
            alpha += &e.a;
            for x in edges.iter_mut() {
                if x.u == e.v {
                    x.u = e.u;
                }
                if x.v == e.v {
                    x.v = e.u;
                }
            }
            if s == e.v {
                s = e.u;
            }
            if t == e.v {
                t = e.u;
            }
            continue;
        }
        break;
    }
    if edges.is_empty() {
        return Some(if s == t { Ev::Val(c, alpha) } else { Ev::Zero });
    }
    if edges.len() == 1 && ((edges[0].u == s && edges[0].v == t) || (edges[0].u == t && edges[0].v == s)) {
        alpha += &edges[0].a;
        return Some(Ev::Val(c, alpha));
    }
    let tr = triangle(&edges, s, t, h, prec, depth, pv);
    if tr.is_none() && depth == 0 {
        let ce: Vec<(usize, usize, Rational)> = edges.iter().map(|e| (e.u, e.v, e.a.clone())).collect();
        let k = canonical(&ce, s, t);
        let mut m = pv.missing.lock().unwrap();
        if !m.contains(&k) {
            m.push(k);
        }
    }
    match tr? {
        Ev::Zero => Some(Ev::Zero),
        Ev::Val(c2, a2) => {
            c *= c2;
            alpha += a2;
            Some(Ev::Val(c, alpha))
        }
    }
}

fn is_pos_int(a: &Rational) -> bool {
    a.denom() == &1 && *a > 0
}

/// IBP triangle rule: (d - 2 a1 - a2 - a3) I = a2 [I(2+,1-) - I(2+,fY-)] + a3 [I(3+,1-) - I(3+,fX-)]
fn triangle(edges: &[E2], s: usize, t: usize, h: &Rational, prec: u32, depth: usize, pv: &Provided) -> Option<Ev> {
    let inc = |x: usize| -> Vec<usize> { (0..edges.len()).filter(|&i| edges[i].u == x || edges[i].v == x).collect() };
    for i1 in 0..edges.len() {
        if !is_pos_int(&edges[i1].a) {
            continue;
        }
        for (x, y) in [(edges[i1].u, edges[i1].v), (edges[i1].v, edges[i1].u)] {
            if x == s || x == t || y == s || y == t {
                continue;
            }
            for &i3 in &inc(x) {
                if i3 == i1 {
                    continue;
                }
                let w = if edges[i3].u == x { edges[i3].v } else { edges[i3].u };
                for &i2 in &inc(y) {
                    if i2 == i1 || i2 == i3 {
                        continue;
                    }
                    let w2 = if edges[i2].u == y { edges[i2].v } else { edges[i2].u };
                    if w2 != w {
                        continue;
                    }
                    let fx: Vec<usize> = inc(x).into_iter().filter(|&i| i != i1 && i != i3).collect();
                    let fy: Vec<usize> = inc(y).into_iter().filter(|&i| i != i1 && i != i2).collect();
                    if fx.len() != 1 || fy.len() != 1 || !is_pos_int(&edges[fx[0]].a) || !is_pos_int(&edges[fy[0]].a) {
                        continue;
                    }
                    let (a1, a2, a3) = (edges[i1].a.clone(), edges[i2].a.clone(), edges[i3].a.clone());
                    let den = Rational::from(h * 2u32) - a1.clone() * 2u32 - a2.clone() - a3.clone();
                    if den == 0 {
                        continue;
                    }
                    let modi = |ch: &[(usize, i64)]| -> Vec<E2> {
                        let mut n = edges.to_vec();
                        for &(i, d) in ch {
                            n[i].a += d;
                        }
                        n
                    };
                    let terms = [
                        (a2.clone(), modi(&[(i2, 1), (i1, -1)])),
                        (Rational::from(-&a2), modi(&[(i2, 1), (fy[0], -1)])),
                        (a3.clone(), modi(&[(i3, 1), (i1, -1)])),
                        (Rational::from(-&a3), modi(&[(i3, 1), (fx[0], -1)])),
                    ];
                    let mut tot = Float::new(prec);
                    let mut al = None;
                    for (cf, g) in terms {
                        match eval2pt(g, s, t, h, prec, depth + 1, pv)? {
                            Ev::Zero => {}
                            Ev::Val(cg, ag) => {
                                tot += Float::with_val(prec, &cg * &cf);
                                al = Some(ag);
                            }
                        }
                    }
                    return Some(match al {
                        None => Ev::Zero,
                        Some(a) => Ev::Val(Float::with_val(prec, &tot / &den), a),
                    });
                }
            }
        }
    }
    None
}

fn find(par: &mut Vec<usize>, x: usize) -> usize {
    let mut r = x;
    while par[r] != r {
        r = par[r];
    }
    let mut y = x;
    while par[y] != r {
        let n = par[y];
        par[y] = r;
        y = n;
    }
    r
}

/// contract index-0 lines; returns (edges, ext, indices, position of e) or Err for numerators
fn contract_zero(g: &Graph, a: &[i64], e: usize) -> Result<(Vec<Edge>, (usize, usize), Vec<i64>, usize), String> {
    let nl = g.edges.len();
    if a[..nl].iter().any(|&x| x < 0) || a[nl..].iter().any(|&x| x != 0) {
        return Err("numerators (negative indices or ISPs) are not handled by the graph rows".into());
    }
    let nv = g.edges.iter().map(|&(u, v)| u.max(v)).max().unwrap_or(0).max(g.ext.0).max(g.ext.1) + 1;
    let mut par: Vec<usize> = (0..nv).collect();
    for i in 0..nl {
        if a[i] == 0 {
            let (u, v) = (find(&mut par, g.edges[i].0), find(&mut par, g.edges[i].1));
            if u != v {
                par[u] = v;
            }
        }
    }
    let keep: Vec<usize> = (0..nl).filter(|&i| a[i] != 0).collect();
    let edges = keep.iter().map(|&i| (find(&mut par, g.edges[i].0), find(&mut par, g.edges[i].1))).collect();
    let ext = (find(&mut par, g.ext.0), find(&mut par, g.ext.1));
    let pos = keep.iter().position(|&i| i == e).ok_or("the massive line has index 0")?;
    Ok((edges, ext, keep.iter().map(|&i| a[i]).collect(), pos))
}

fn bridgeless(edges: &[Edge]) -> bool {
    let e2: Vec<E2> = edges.iter().map(|&(u, v)| E2 { u, v, a: q(1) }).collect();
    (0..edges.len()).all(|i| connected(&e2, edges[i].0, edges[i].1, Some(i)))
}

/// massless value of the family integral at indices a (no numerators), None if not reducible here
pub fn massless(g: &Graph, a: &[i64], h: &Rational, prec: u32, pv: &Provided) -> Option<Float> {
    let nl = g.edges.len();
    if a[..nl].iter().any(|&x| x < 0) || a[nl..].iter().any(|&x| x != 0) {
        return None;
    }
    let nv = g.edges.iter().map(|&(u, v)| u.max(v)).max().unwrap_or(0).max(g.ext.0).max(g.ext.1) + 1;
    let mut par: Vec<usize> = (0..nv).collect();
    for i in 0..nl {
        if a[i] == 0 {
            let (u, v) = (find(&mut par, g.edges[i].0), find(&mut par, g.edges[i].1));
            if u != v {
                par[u] = v;
            }
        }
    }
    let edges: Vec<E2> = (0..nl).filter(|&i| a[i] != 0).map(|i| E2 { u: find(&mut par, g.edges[i].0), v: find(&mut par, g.edges[i].1), a: q(a[i]) }).collect();
    let (s, t) = (find(&mut par, g.ext.0), find(&mut par, g.ext.1));
    if s == t {
        return Some(Float::new(prec));
    }
    let mut vs: Vec<usize> = edges.iter().flat_map(|e| [e.u, e.v]).collect();
    vs.push(s);
    vs.push(t);
    vs.sort();
    vs.dedup();
    if (edges.len() as i64) - (vs.len() as i64) + 1 < g.loops as i64 {
        return Some(Float::new(prec)); // an empty loop: scaleless
    }
    match eval2pt(edges, s, t, h, prec, 0, pv)? {
        Ev::Zero => Some(Float::new(prec)),
        Ev::Val(c, _) => Some(c),
    }
}

/// leading coefficient of every region class: side 0 -> {class k: [(n, value)]} with exponent z^(k h + n) (n = -sum_S a);
/// side infinity -> x = 1/z exponent x^(n - k h) (n = sum_gamma a). None if a needed evaluation is not reducible.
pub fn regions(g: &Graph, e: usize, a: &[i64], h: &Rational, prec: u32, at_zero: bool, pv: &Provided) -> Result<Option<Vec<(i64, i64, Float)>>, String> {
    let (edges, (ea, eb), aa, ee) = contract_zero(g, a, e)?;
    let nl = edges.len();
    let mut nv = edges.iter().map(|&(u, v)| u.max(v)).max().unwrap_or(0).max(ea).max(eb) + 1;
    nv = nv.max(1);
    let mut out: Vec<(i64, i64, Float)> = Vec::new();
    let others: Vec<usize> = (0..nl).filter(|&i| i != ee).collect();
    for mask in 0u64..(1u64 << others.len()) {
        let mut sset = vec![ee];
        for (k, &i) in others.iter().enumerate() {
            if mask >> k & 1 == 1 {
                sset.push(i);
            }
        }
        let inset = |i: usize| sset.contains(&i);
        let val;
        let class;
        let nexp;
        if at_zero {
            let hard: Vec<usize> = (0..nl).filter(|&i| !inset(i)).collect();
            let hedges: Vec<E2> = hard.iter().map(|&i| E2 { u: edges[i].0, v: edges[i].1, a: q(aa[i]) }).collect();
            if !connected(&hedges, ea, eb, None) {
                continue;
            }
            let mut he: Vec<Edge> = hard.iter().map(|&i| edges[i]).collect();
            he.push((ea, eb));
            if !bridgeless(&he) {
                continue;
            }
            // identify H's vertices (H connected: one representative)
            let mut par: Vec<usize> = (0..nv).collect();
            for &i in &hard {
                let (u, v) = (find(&mut par, edges[i].0), find(&mut par, edges[i].1));
                if u != v {
                    par[u] = v;
                }
            }
            let sedges: Vec<Edge> = sset.iter().map(|&i| (find(&mut par, edges[i].0), find(&mut par, edges[i].1))).collect();
            if !bridgeless(&sedges) {
                continue;
            }
            let mut sv: Vec<usize> = sedges.iter().flat_map(|&(u, v)| [u, v]).collect();
            sv.sort();
            sv.dedup();
            let ls = sedges.len() as i64 - sv.len() as i64 + 1;
            let ch = match eval2pt(hedges, ea, eb, h, prec, 0, pv) {
                None => return Ok(None),
                Some(Ev::Zero) => continue,
                Some(Ev::Val(c, _)) => c,
            };
            let (u, v) = sedges[0];
            let x: Vec<E2> = (1..sset.len()).map(|j| E2 { u: sedges[j].0, v: sedges[j].1, a: q(aa[sset[j]]) }).collect();
            let (cx, alx) = if u == v {
                if !x.is_empty() {
                    continue;
                }
                (Float::with_val(prec, 1), q(0))
            } else {
                match eval2pt(x, u, v, h, prec, 0, pv) {
                    None => return Ok(None),
                    Some(Ev::Zero) => continue,
                    Some(Ev::Val(c, a)) => (c, a),
                }
            };
            let Some(vmv) = vm(&q(aa[ee]), &alx, h, prec) else { continue };
            val = Float::with_val(prec, &ch * &cx) * vmv;
            class = ls;
            nexp = -sset.iter().map(|&i| aa[i]).sum::<i64>();
        } else {
            let gedges: Vec<Edge> = sset.iter().map(|&i| edges[i]).collect();
            let mut gv: Vec<usize> = gedges.iter().flat_map(|&(u, v)| [u, v]).collect();
            gv.sort();
            gv.dedup();
            let ge2: Vec<E2> = gedges.iter().map(|&(u, v)| E2 { u, v, a: q(1) }).collect();
            if !gv.iter().all(|&x| connected(&ge2, gv[0], x, None)) || !bridgeless(&gedges) {
                continue;
            }
            let lg = gedges.len() as i64 - gv.len() as i64 + 1;
            let mut par: Vec<usize> = (0..nv).collect();
            for &(u0, v0) in &gedges {
                let (u, v) = (find(&mut par, u0), find(&mut par, v0));
                if u != v {
                    par[u] = v;
                }
            }
            let co: Vec<E2> = (0..nl).filter(|&i| !inset(i)).map(|i| E2 { u: find(&mut par, edges[i].0), v: find(&mut par, edges[i].1), a: q(aa[i]) }).collect();
            let (u, v) = edges[ee];
            let x: Vec<E2> = sset[1..].iter().map(|&i| E2 { u: edges[i].0, v: edges[i].1, a: q(aa[i]) }).collect();
            let (cx, alx) = match eval2pt(x, u, v, h, prec, 0, pv) {
                None => return Ok(None),
                Some(Ev::Zero) => continue,
                Some(Ev::Val(c, a)) => (c, a),
            };
            let (fa, fb) = (find(&mut par, ea), find(&mut par, eb));
            let cc = if fa == fb {
                if !co.is_empty() {
                    continue;
                }
                Float::with_val(prec, 1)
            } else {
                match eval2pt(co, fa, fb, h, prec, 0, pv) {
                    None => return Ok(None),
                    Some(Ev::Zero) => continue,
                    Some(Ev::Val(c, _)) => c,
                }
            };
            let Some(vmv) = vm(&q(aa[ee]), &alx, h, prec) else { continue };
            val = Float::with_val(prec, &cx * &vmv) * cc;
            class = lg;
            nexp = sset.iter().map(|&i| aa[i]).sum::<i64>();
        }
        match out.iter().position(|(k, n, _)| *k == class && *n == nexp) {
            Some(i) => out[i].2 += &val,
            None => out.push((class, nexp, val)),
        }
    }
    Ok(Some(out))
}

fn poch(x: &Rational, n: usize) -> Rational {
    let mut r = q(1);
    for i in 0..n {
        r *= Rational::from(x + i as u32);
    }
    r
}

/// (kp, m) with expo = kp h - m, m integer
fn classify(expo: &Rational, h: &Rational) -> Option<(i64, i64)> {
    for kp in 0..12i64 {
        let m = Rational::from(h * kp) - expo.clone();
        if m.denom() == &1 {
            return Some((kp, m.numer().to_i64()?));
        }
    }
    None
}

/// simple masters: ('power', C, expo) or ('bubble', C, ae, beta, b) after massless series/parallel reductions; see graph.py
enum Simple {
    Power(Float, Rational),
    Bubble(Float, i64, Rational, Rational),
}

/// split at a cut vertex v (v != a, b; removing v disconnects a from b): returns (edges of the part with the massive line, edges of
/// the massless part, externals of the massive part). None if there is no such vertex, if a third (vacuum) part hangs at v, or if the
/// massive line is not in one of the two parts.
fn series_split(edges: &[(usize, usize, Rational, bool)], a: usize, b: usize) -> Option<(Vec<(usize, usize, Rational, bool)>, Vec<(usize, usize, Rational, bool)>, usize, usize)> {
    let mut verts: Vec<usize> = edges.iter().flat_map(|x| [x.0, x.1]).collect();
    verts.sort();
    verts.dedup();
    let reach = |from: usize, cut: usize| -> Vec<usize> {
        let mut seen = vec![from];
        let mut st = vec![from];
        while let Some(x) = st.pop() {
            for e in edges {
                let y = if e.0 == x { e.1 } else if e.1 == x { e.0 } else { continue };
                if y != cut && !seen.contains(&y) {
                    seen.push(y);
                    st.push(y);
                }
            }
        }
        seen
    };
    for &v in verts.iter().filter(|&&v| v != a && v != b) {
        let ra = reach(a, v);
        if ra.contains(&b) {
            continue;
        }
        let rb = reach(b, v);
        // every edge must belong to the a-side or the b-side (no vacuum blob hanging at v)
        let side = |x: &(usize, usize, Rational, bool)| -> Option<bool> {
            let w = if x.0 == v { x.1 } else { x.0 };
            if ra.contains(&w) { Some(true) } else if rb.contains(&w) { Some(false) } else { None }
        };
        let mut pa = Vec::new();
        let mut pb = Vec::new();
        for x in edges {
            if x.0 == v && x.1 == v {
                return None;
            }
            match side(x)? {
                true => pa.push(x.clone()),
                false => pb.push(x.clone()),
            }
        }
        let ea_side = pa.iter().any(|x| x.3);
        return if ea_side { Some((pa, pb, a, v)) } else { Some((pb, pa, v, b)) };
    }
    None
}

fn simple_form(g: &Graph, e: usize, a: &[i64], h: &Rational, prec: u32) -> Option<Simple> {
    let (edges0, (mut ea, mut eb), aa, ee) = contract_zero(g, a, e).ok()?;
    if ea == eb {
        return None;
    }
    // (u, v, index, massive)
    let mut edges: Vec<(usize, usize, Rational, bool)> = edges0.iter().enumerate().map(|(i, &(u, v))| (u, v, q(aa[i]), i == ee)).collect();
    let mut c = Float::with_val(prec, 1);
    loop {
        let em = edges.iter().find(|x| x.3).unwrap().clone();
        if em.0 == em.1 {
            let rest: Vec<E2> = edges.iter().filter(|x| !x.3).map(|x| E2 { u: x.0, v: x.1, a: x.2.clone() }).collect();
            return match eval2pt(rest, ea, eb, h, prec, 0, &Provided::default())? {
                Ev::Zero => Some(Simple::Power(Float::new(prec), q(0))),
                Ev::Val(cr, _) => {
                    let t = tad(&em.2, h, prec)?;
                    Some(Simple::Power(Float::with_val(prec, &c * &cr) * t, Rational::from(h - &em.2)))
                }
            };
        }
        let mut did = false;
        'p: for i in 0..edges.len() {
            if edges[i].3 {
                continue;
            }
            if edges[i].0 == edges[i].1 {
                return Some(Simple::Power(Float::new(prec), q(0)));
            }
            for j in 0..i {
                if edges[j].3 {
                    continue;
                }
                let same = (edges[i].0 == edges[j].0 && edges[i].1 == edges[j].1) || (edges[i].0 == edges[j].1 && edges[i].1 == edges[j].0);
                if same {
                    c *= g_bubble(&edges[j].2, &edges[i].2, h, prec)?;
                    edges[j].2 = Rational::from(&edges[j].2 + &edges[i].2) - h.clone();
                    edges.remove(i);
                    did = true;
                    break 'p;
                }
            }
        }
        if did {
            continue;
        }
        let mut deg = std::collections::BTreeMap::new();
        for x in &edges {
            *deg.entry(x.0).or_insert(0usize) += 1;
            *deg.entry(x.1).or_insert(0usize) += 1;
        }
        let mut merged = false;
        for (&v, &dv) in &deg {
            if dv == 2 && v != ea && v != eb {
                let idx: Vec<usize> = (0..edges.len()).filter(|&i| edges[i].0 == v || edges[i].1 == v).collect();
                if edges[idx[0]].3 || edges[idx[1]].3 {
                    continue;
                }
                let (x1, x2) = (edges[idx[0]].clone(), edges[idx[1]].clone());
                let u = if x1.1 == v { x1.0 } else { x1.1 };
                let w = if x2.1 == v { x2.0 } else { x2.1 };
                edges.remove(idx[1]);
                edges.remove(idx[0]);
                edges.push((u, w, Rational::from(&x1.2 + &x2.2), false));
                merged = true;
                break;
            }
        }
        if merged {
            continue;
        }
        // a cut vertex v splits the graph into two two-point parts in series (both carry p): the massless part is a factor
        // (its value at p^2 = 1), the part with the massive line stays; hanging vacuum blobs are scaleless
        if let Some((part_e, part_m, na, nb)) = series_split(&edges, ea, eb) {
            let rest: Vec<E2> = part_m.iter().map(|x| E2 { u: x.0, v: x.1, a: x.2.clone() }).collect();
            let (ma, mb) = if na == ea { (nb, eb) } else { (ea, na) };
            match eval2pt(rest, ma, mb, h, prec, 0, &Provided::default())? {
                Ev::Zero => return Some(Simple::Power(Float::new(prec), q(0))),
                Ev::Val(cr, _) => c *= cr,
            }
            edges = part_e;
            (ea, eb) = (na, nb);
            continue;
        }
        break;
    }
    let em = edges.iter().find(|x| x.3).unwrap().clone();
    let rest: Vec<_> = edges.iter().filter(|x| !x.3).cloned().collect();
    let ab = |x: &(usize, usize, Rational, bool)| (x.0 == ea && x.1 == eb) || (x.0 == eb && x.1 == ea);
    if ab(&em) && rest.len() == 1 && ab(&rest[0]) {
        return Some(Simple::Bubble(c, rat_i64(&em.2)?, q(0), rest[0].2.clone()));
    }
    if rest.len() == 2 {
        for k in 0..2 {
            let (s_, o) = (&rest[k], &rest[1 - k]);
            if !ab(o) {
                continue;
            }
            let ends_e = [em.0, em.1];
            let ends_s = [s_.0, s_.1];
            let mid: Vec<usize> = ends_e.iter().filter(|x| ends_s.contains(x)).cloned().collect();
            if mid.len() == 1 && mid[0] != ea && mid[0] != eb {
                let other_e = if em.0 == mid[0] { em.1 } else { em.0 };
                let other_s = if s_.0 == mid[0] { s_.1 } else { s_.0 };
                if (other_e == ea && other_s == eb) || (other_e == eb && other_s == ea) {
                    return Some(Simple::Bubble(c, rat_i64(&em.2)?, s_.2.clone(), o.2.clone()));
                }
            }
        }
    }
    None
}

/// all-order rows [(n, value)] of class k for simple masters (None = not simple); same n convention as `regions`
pub fn exact_series(g: &Graph, e: usize, a: &[i64], h: &Rational, at_zero: bool, k: i64, nmax: usize, prec: u32) -> Option<Vec<(i64, Float)>> {
    match simple_form(g, e, a, h, prec)? {
        Simple::Power(c, expo) => {
            if c.is_zero() {
                return Some(vec![]);
            }
            let (kp, m) = classify(&expo, h)?;
            if kp != k {
                return Some(vec![]);
            }
            Some(vec![(if at_zero { -m } else { m }, c)])
        }
        Simple::Bubble(c, ae, beta, b) => {
            let aeq = q(ae);
            let mut out = Vec::new();
            let ang = |n: usize| -> Rational {
                // (b)_n (b - h + 1)_n / (n! (h)_n)
                let num = poch(&b, n) * poch(&(Rational::from(&b - h) + 1u32), n);
                let den = poch(&q(1), n) * poch(h, n);
                num / den
            };
            if at_zero {
                let (kp, m) = classify(&(Rational::from(h - &aeq) - beta.clone()), h)?;
                if kp != k {
                    return Some(vec![]);
                }
                for n in 0..nmax {
                    let g1 = gamma(&(Rational::from(h + n as u32) - beta.clone()), prec)?;
                    let g2 = match gamma(&(Rational::from(&aeq + &beta) - h.clone() - n as u32), prec) {
                        Some(x) => x,
                        None => continue,
                    };
                    let mut v = Float::with_val(prec, &c * &ang(n));
                    v *= g1;
                    v *= g2;
                    v *= rgamma(h, prec);
                    v *= rgamma(&aeq, prec);
                    out.push((-m + n as i64, v));
                }
            } else {
                let (kp, m) = classify(&(Rational::from(h - &aeq) - beta.clone() - b.clone()), h)?;
                if kp != k {
                    return Some(vec![]);
                }
                for n in 0..nmax {
                    let Some(vmv) = vm(&aeq, &(Rational::from(&beta + &b) + n as u32), h, prec) else { continue };
                    out.push((m + n as i64, Float::with_val(prec, &c * &ang(n)) * vmv));
                }
            }
            Some(out)
        }
    }
}
