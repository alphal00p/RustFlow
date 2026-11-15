//! PFLOW tests: seconds-scale at one loop (hand IBP DE) and graph-level checks; two- and three-loop RustRed DEs as #[ignore].
use super::engine::{PConfig, run_b};
use super::graph::{self, Graph};
use crate::dispersive::system::{EBlockSystem, RuleSet, parse_coefficient};
use rug::{Float, Rational};

fn graph_of(name: &str) -> Graph {
    let p = format!("{}/fixtures/pflow/graph-{name}.json", env!("CARGO_MANIFEST_DIR"));
    let f: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap();
    let edges = f["lines"].as_array().unwrap().iter().map(|l| (l["endpoints"][0].as_u64().unwrap() as usize, l["endpoints"][1].as_u64().unwrap() as usize)).collect();
    Graph { edges, ext: (f["external"]["in"].as_u64().unwrap() as usize, f["external"]["out"].as_u64().unwrap() as usize), loops: f["loops"].as_u64().unwrap() as usize }
}

fn digits(a: &Float, b: &Float) -> f64 {
    let d = Float::with_val(a.prec(), a - b);
    if d.is_zero() { 999.0 } else { -(d.abs() / Float::with_val(a.prec(), b.abs_ref())).log10().to_f64() }
}

/// one-loop bubble, mass on line 0 (Minkowski RustRed conventions): masters I(1,1), I(1,0);
/// E(2,1) = [(1-h) T/z - (d-3) B]/(z+1) (E-normalisation) written as a harness rule in I = (-1)^(sum a) E
fn one_loop_system(eps: Rational) -> EBlockSystem {
    let d = Rational::from(4) - eps * 2u32;
    let h = Rational::from(&d / 2u32);
    let cb = Rational::from(&d - 3u32);
    let ct = Rational::from(1 - h.clone());
    let sol = serde_json::json!({"schema": "rustred-harness-solution.v1", "fixed": {"d": d.to_string()}, "masters": [[1,1],[1,0]],
        "rules": [{"target": [2,1], "terms": [{"integral": [1,1], "coefficient": format!("({cb})/(1+z)")}, {"integral": [1,0], "coefficient": format!("({ct})/(z+z^2)")}]},
                  {"target": [2,0], "terms": [{"integral": [1,0], "coefficient": format!("({})/(z)", Rational::from(-&ct))}]}]});
    let rs = RuleSet::from_json(&sol).unwrap();
    EBlockSystem::from_rules_with(&rs, 0, &|c: &str| parse_coefficient(c, "z")).unwrap()
}

#[test]
fn one_loop_bubble_route_b() {
    for eps in [Rational::from((1, 100)), Rational::from((-3, 10))] {
        let sys = one_loop_system(eps.clone());
        let g = graph_of("P1L");
        let mut cfg = PConfig::new(40, 50, 60, 1);
        cfg.moments = vec![1, 2];
        let r = run_b(&sys, &g, 0, &cfg).unwrap();
        let h = Rational::from(&sys.d0 / 2u32);
        let prec = r.moments[0].1[0].prec();
        // I(1,1) = E(1,1); moment m = 2 gives E(2,1) = M(-1) (a_e = 1), and I(2,1) = -E(2,1) is not a master: compare E
        let ex1 = graph::g_bubble(&Rational::from(1), &Rational::from(1), &h, prec).unwrap();
        let ex2 = graph::g_bubble(&Rational::from(2), &Rational::from(1), &h, prec).unwrap();
        let d1 = digits(&r.moments[0].1[0], &ex1);
        let d2 = digits(&r.moments[1].1[0], &ex2);
        assert!(d1 > 30.0 && d2 > 28.0, "eps {eps}: {d1} {d2} {}", r.diag);
        // the massive tadpole has no analytic part
        assert!(r.moments[0].1[1].clone().abs() < Float::with_val(prec, 1e-30));
    }
}

#[test]
fn triangle_rule_kite() {
    let g = graph_of("P2L");
    let h = Rational::from(2) - Rational::from((1, 100));
    let prec = 200;
    let k = graph::massless(&g, &[1, 1, 1, 1, 1], &h, prec).unwrap();
    // 2 G(1,1) [G(2,1) - G(2, 3-h)] / (d - 4)
    let o = Rational::from(1);
    let t = Rational::from(2);
    let g11 = graph::g_bubble(&o, &o, &h, prec).unwrap();
    let g21 = graph::g_bubble(&t, &o, &h, prec).unwrap();
    let g2 = graph::g_bubble(&t, &(Rational::from(3) - h.clone()), &h, prec).unwrap();
    let ex = Float::with_val(prec, &g21 - &g2) * g11 * 2u32 / Float::with_val(prec, Rational::from(&h * 2u32) - 4u32);
    assert!(digits(&k, &ex) > 55.0);
}

fn run_fixture(file: &str, fam: &str, w: u32, m: i64) -> (EBlockSystem, super::engine::PResult, Graph) {
    let p = format!("{}/fixtures/pflow/{file}", env!("CARGO_MANIFEST_DIR"));
    let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap();
    let rs = RuleSet::from_json(&v).unwrap();
    let sys = EBlockSystem::from_rules_with(&rs, 0, &|c: &str| parse_coefficient(c, "z")).unwrap();
    let g = graph_of(fam);
    let mut cfg = PConfig::new(w, w + 20, m, g.loops as u32);
    cfg.exact_rows = m as usize;
    let r = run_b(&sys, &g, 0, &cfg).unwrap();
    (sys, r, g)
}

fn check_exact(sys: &EBlockSystem, r: &super::engine::PResult, g: &Graph, min_digits: f64) -> usize {
    let h = Rational::from(&sys.d0 / 2u32);
    let mut n = 0;
    for (j, a) in sys.masters.iter().enumerate() {
        let prec = r.moments[0].1[j].prec();
        let Some(ex) = graph::massless(g, a, &h, prec) else { continue };
        let sgn = if a.iter().sum::<i64>() % 2 != 0 { -1 } else { 1 };
        let v = Float::with_val(prec, &r.moments[0].1[j] * sgn);
        if ex.is_zero() {
            assert!(v.clone().abs() < Float::with_val(prec, 1e-60), "{a:?}: {v} should vanish");
        } else {
            let d = digits(&v, &ex);
            assert!(d > min_digits, "{a:?}: {d} digits, diag {}", r.diag);
        }
        n += 1;
    }
    n
}

#[test]
#[ignore = "two-loop kite family at 80+ digits (about a minute in release)"]
fn kite_family_80_digits() {
    let (sys, r, g) = run_fixture("P2L_m0-solution-199_50.json", "P2L", 100, 170);
    assert!(check_exact(&sys, &r, &g, 80.0) >= 6);
}

#[test]
#[ignore = "three-loop NO family at 80+ digits (minutes in release)"]
fn no_family_80_digits() {
    let (sys, r, g) = run_fixture("P3L_NO_m0-solution-62811_15625.json", "P3L_NO", 100, 170);
    assert!(check_exact(&sys, &r, &g, 80.0) >= 14);
}
