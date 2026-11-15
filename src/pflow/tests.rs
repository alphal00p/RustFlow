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
    let k = graph::massless(&g, &[1, 1, 1, 1, 1], &h, prec, &graph::Provided::default()).unwrap();
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
        let Some(ex) = graph::massless(g, a, &h, prec, &graph::Provided::default()) else { continue };
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

/// bug 1 (source fix): the singular-point finder works in eta = z - 1 and deflates eta = -1, i.e. z = 0 ... the threshold z = -1
/// must survive (it was dropped when the denominators were handed over in z)
#[test]
fn singular_points_keep_the_threshold() {
    let sys = one_loop_system(Rational::from((1, 100)));
    let pts = super::engine::singular_points(&sys);
    assert!(pts.iter().any(|p| (p.0 + 1.0).abs() < 1e-12 && p.1.abs() < 1e-12), "z = -1 missing: {pts:?}");
    assert!(pts.iter().any(|p| p.0.abs() < 1e-12 && p.1.abs() < 1e-12), "z = 0 missing: {pts:?}");
}

/// bug 2 (source fix): rank-deficient least squares returns the minimum-pivot solution, its rank and the null directions, no NaN
#[test]
fn lstsq_rank_deficient_is_finite() {
    let prec = 200;
    let f = |x: f64| Float::with_val(prec, x);
    let a = vec![vec![f(1.0), f(1.0), f(0.0)], vec![f(2.0), f(2.0), f(0.0)], vec![f(0.0), f(0.0), f(3.0)]];
    let b = vec![f(1.0), f(2.0), f(6.0)];
    let (x, rank, null) = super::engine::lstsq_rr_null(&a, &b, prec);
    assert_eq!(rank, 2);
    assert_eq!(null.len(), 1);
    assert!(x.iter().all(|v| v.is_finite()), "{x:?}");
    for (row, bi) in a.iter().zip(&b) {
        let v = Float::with_val(prec, Float::dot(row.iter().zip(&x))) - bi;
        assert!(v.abs() < 1e-50, "residual");
    }
    // the null direction is (1, -1, 0) up to scale
    let n = &null[0];
    assert!(Float::with_val(prec, &n[0] + &n[1]).abs() < 1e-50 && n[2].clone().abs() < 1e-50 && n[0].clone().abs() > 1e-10);
}

/// bug 3 (source fix): split points follow spurious singular points at O(eps) from 0 and -1 and at O(1/eps^2)
#[test]
fn split_points_follow_spurious_points() {
    let cfg = PConfig::new(60, 80, 100, 4);
    let pts = [(0.0, 0.0), (-1.0, 0.0), (-4e-5, 0.0), (-1.0 - 2e-5, 0.0), (-5.1e10, 0.0), (1.0, 0.0)];
    let sp = super::engine::split_points(&cfg, &pts);
    assert!(sp.s0.to_f64() <= 0.4 * 4e-5 * (1.0 + 1e-6) && sp.s0.to_f64() > 0.0);
    assert!(sp.r.to_f64() <= 0.4 * 2e-5 * (1.0 + 1e-6) && sp.r.to_f64() > 0.0);
    assert!(sp.s_end.to_f64() >= 2.5 * 5.1e10 * (1.0 - 1e-6));
    // no spurious points: the configured splits stay
    let sp = super::engine::split_points(&cfg, &[(0.0, 0.0), (-1.0, 0.0)]);
    assert_eq!((sp.s0, sp.r, sp.s_end), (cfg.s0.clone(), cfg.r.clone(), cfg.s_end.clone()));
}

/// row precision budget: a spurious point at distance delta costs log10(radius/delta) digits per order
#[test]
fn order_budget_from_spurious_distance() {
    let (nmax, loss) = super::engine::order_budget(-2, 220, 60, 1.0, 2.2e-5);
    assert!((loss - 4.657).abs() < 1e-3);
    assert_eq!(nmax, -2 + 32);
    assert_eq!(super::engine::order_budget(0, 80, 60, 1.0, 1.0).0, i64::MAX);
}

/// canonical keys of the provided lower-loop values agree with s1/proto/prov.py (relabelling, index-0 lines, bridges, series)
#[test]
fn canonical_key_matches_python() {
    let r = |x: i64| Rational::from(x);
    let kite = vec![(1, 3, r(1)), (1, 4, r(1)), (3, 4, r(1)), (3, 2, r(1)), (4, 2, r(2))];
    let want = r#"[[0, 2, "1"], [0, 3, "1"], [1, 2, "1"], [1, 3, "2"], [2, 3, "1"]]"#;
    assert_eq!(graph::canonical(&kite, 1, 2), want);
    // relabelled, with an index-0 line and a bridge on the external vertex
    let other = vec![(7, 5, r(1)), (9, 7, r(2)), (5, 9, r(1)), (9, 8, r(1)), (5, 8, r(1)), (8, 6, r(0)), (6, 1, r(1))];
    assert_eq!(graph::canonical(&other, 1, 7), want);
    // a three-loop hard part of the four-loop NO1 family (prov-1990761.json, part 0)
    let e: Vec<(usize, usize, Rational)> = [[4, 7, 1], [5, 8, 1], [6, 8, 1], [5, 2, 1], [6, 7, 1], [8, 4, 1], [7, 3, 1], [3, 5, 1], [1, 4, 2], [3, 1, 1]]
        .iter().map(|x| (x[0] as usize, x[1] as usize, r(x[2]))).collect();
    assert_eq!(graph::canonical(&e, 1, 2), r#"[[0, 2, "1"], [0, 3, "1"], [1, 2, "1"], [1, 4, "2"], [2, 5, "1"], [3, 4, "1"], [3, 5, "2"], [4, 5, "1"]]"#);
}

/// one driver: the propagator class config gives the same run as the flags (two-loop kite family at fixed d0)
#[test]
fn propagator_class_config_equals_flags() {
    let dir = format!("{}/fixtures/pflow", env!("CARGO_MANIFEST_DIR"));
    let de = format!("{dir}/P2L_m0-solution-199_50.json");
    let cfgp = format!("{dir}/config-P2L.json");
    assert!(super::cli::is_propagator_config(&cfgp));
    let s = |x: &str| x.to_string();
    let a = super::cli::run(&[s("--de"), de.clone(), s("--class-config"), cfgp]).unwrap();
    let b = super::cli::run(&[s("--de"), de, s("--graph"), format!("{dir}/graph-P2L.json"), s("--e"), s("0"), s("--w"), s("40"), s("--wf"), s("60"), s("--m"), s("80"), s("--moments"), s("1,2")]).unwrap();
    assert_eq!(a["values"], b["values"]);
    assert!(!a["values"].as_array().unwrap().is_empty());
}

/// simple masters through a cut vertex: (massive bubble) x (massless two-point part) in series; the all-order rows exist and their
/// lowest rows equal the region evaluator's leading rows on both sides (four-loop NO1, sector {0,4,5,6,7,8,9})
#[test]
fn simple_form_through_a_cut_vertex() {
    let g = graph_of("P4L_NO1");
    let h = Rational::from((1990761, 1000000));
    let prec = 200;
    for a in [[1i64, 0, 0, 0, 1, 1, 1, 1, 2, 1, 0, 0, 0, 0], [1, 0, 0, 0, 1, 1, 1, 1, 1, 2, 0, 0, 0, 0]] {
        for at_zero in [true, false] {
            let reg = graph::regions(&g, 0, &a, &h, prec, at_zero, &graph::Provided::default()).unwrap().unwrap();
            let k = reg[0].0;
            let mine: Vec<_> = reg.iter().filter(|x| x.0 == k).collect();
            let nmin = mine.iter().map(|x| x.1).min().unwrap();
            let mut lead = Float::new(prec);
            for x in mine.iter().filter(|x| x.1 == nmin) {
                lead += &x.2;
            }
            let ex = graph::exact_series(&g, 0, &a, &h, at_zero, k, 12, prec).expect("simple through the cut vertex");
            let (n0, v0) = ex.iter().min_by_key(|x| x.0).unwrap();
            assert_eq!(*n0, nmin, "{a:?} at_zero {at_zero}");
            assert!(digits(v0, &lead) > 50.0, "{a:?} at_zero {at_zero}: {v0} vs {lead}");
        }
    }
}
