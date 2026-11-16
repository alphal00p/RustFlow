//! PFLOW tests: seconds-scale at one loop (hand IBP DE) and graph-level checks; two- and three-loop RustRed DEs as #[ignore].
use super::engine::{PConfig, run_b};
use crate::dispersive::run;
use super::graph::{self, Graph};
use crate::dispersive::system::{EBlockSystem, RuleSet, parse_coefficient};
use rug::ops::Pow;
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

/// every family with nonzero boundary rows is fitted with rank > 0, and no physical row is dropped by the budget (H0)
fn check_ranks(r: &super::engine::PResult) {
    for f in r.diag["families"].as_array().unwrap() {
        if f["nonzero_rows"].as_u64().unwrap() > 0 {
            assert!(f["rank"].as_u64().unwrap() > 0, "family fitted with rank 0: {f}");
            assert!(f["rows_dropped_budget"].as_u64().unwrap() < f["rows_offered"].as_u64().unwrap(), "every row dropped: {f}");
        }
    }
    assert!(r.accepted, "not accepted: {}", r.diag["acceptance"]["failures"]);
}

#[test]
#[ignore = "two-loop kite family at 80+ digits (about a minute in release)"]
fn kite_family_80_digits() {
    let (sys, r, g) = run_fixture("P2L_m0-solution-199_50.json", "P2L", 100, 170);
    check_ranks(&r);
    assert!(check_exact(&sys, &r, &g, 80.0) >= 6);
}

#[test]
#[ignore = "three-loop NO family at 80+ digits (minutes in release)"]
fn no_family_80_digits() {
    let (sys, r, g) = run_fixture("P3L_NO_m0-solution-62811_15625.json", "P3L_NO", 100, 170);
    check_ranks(&r);
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

/// canonical keys of the provided lower-loop values agree with the reference strings of the boundary request tools (relabelling, index-0 lines, bridges, series)
#[test]
fn canonical_key_matches_python() {
    let r = |x: i64| Rational::from(x);
    let kite = vec![(1, 3, r(1)), (1, 4, r(1)), (3, 4, r(1)), (3, 2, r(1)), (4, 2, r(2))];
    let want = r#"[[0, 2, "1"], [0, 3, "1"], [1, 2, "1"], [1, 3, "2"], [2, 3, "1"]]"#;
    assert_eq!(graph::canonical(&kite, 1, 2), want);
    // relabelled, with an index-0 line and a bridge on the external vertex
    let other = vec![(7, 5, r(1)), (9, 7, r(2)), (5, 9, r(1)), (9, 8, r(1)), (5, 8, r(1)), (8, 6, r(0)), (6, 1, r(1))];
    assert_eq!(graph::canonical(&other, 1, 7), want);
    // a three-loop hard part of the four-loop NO1 family (a hard part of the NO1 boundary)
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

/// per-cut DE blocks: the two-loop kite DE split into two blocks gives the same run; a block with a different rule is refused
#[test]
fn de_blocks_merge_and_conflicts() {
    let dir = format!("{}/fixtures/pflow", env!("CARGO_MANIFEST_DIR"));
    let de = format!("{dir}/P2L_m0-solution-199_50.json");
    let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&de).unwrap()).unwrap();
    let rules = v["rules"].as_array().unwrap().clone();
    let tmp = std::env::temp_dir().join(format!("pflow-blocks-{}", std::process::id()));
    std::fs::create_dir_all(&tmp).unwrap();
    let mut a = v.clone();
    let mut b = v.clone();
    a["rules"] = serde_json::Value::Array(rules.iter().step_by(2).cloned().collect());
    b["rules"] = serde_json::Value::Array(rules.iter().skip(1).step_by(2).cloned().collect());
    let (pa, pb) = (tmp.join("a.json"), tmp.join("b.json"));
    std::fs::write(&pa, a.to_string()).unwrap();
    std::fs::write(&pb, b.to_string()).unwrap();
    let s = |x: &str| x.to_string();
    let whole = super::cli::run(&[s("--de"), de.clone(), s("--graph"), format!("{dir}/graph-P2L.json"), s("--e"), s("0"), s("--w"), s("30"), s("--wf"), s("45"), s("--m"), s("60")]).unwrap();
    let split = super::cli::run(&[s("--de"), format!("{},{}", pa.display(), pb.display()), s("--graph"), format!("{dir}/graph-P2L.json"), s("--e"), s("0"), s("--w"), s("30"), s("--wf"), s("45"), s("--m"), s("60")]).unwrap();
    assert_eq!(whole["values"], split["values"]);
    // a conflicting rule for the same target
    let mut c = v.clone();
    let mut r0 = rules[0].clone();
    r0["terms"][0]["coefficient"] = serde_json::Value::String("(1)/(1+z)".into());
    c["rules"] = serde_json::Value::Array(vec![r0]);
    let pc = tmp.join("c.json");
    std::fs::write(&pc, c.to_string()).unwrap();
    assert!(matches!(super::cli::load_blocks(&format!("{},{}", de, pc.display())), Err(e) if e.contains("differs")));
    let _ = std::fs::remove_dir_all(&tmp);
}

/// checkpoints: a second run in the same work directory resumes every transport and reproduces the values exactly
#[test]
fn workdir_resume_reproduces() {
    let sys = one_loop_system(Rational::from((1, 100)));
    let g = graph_of("P1L");
    let tmp = std::env::temp_dir().join(format!("pflow-wd-{}", std::process::id()));
    let mut cfg = PConfig::new(30, 40, 50, 1);
    cfg.moments = vec![1, 2];
    cfg.workdir = Some(tmp.to_string_lossy().into_owned());
    std::fs::create_dir_all(&tmp).unwrap();
    let r1 = run_b(&sys, &g, 0, &cfg).unwrap();
    assert!(std::fs::read_dir(&tmp).unwrap().any(|f| f.unwrap().file_name().to_string_lossy().starts_with("tr-")));
    let r2 = run_b(&sys, &g, 0, &cfg).unwrap();
    for (a, b) in r1.moments.iter().zip(&r2.moments) {
        assert_eq!(a.1, b.1);
    }
    let _ = std::fs::remove_dir_all(&tmp);
}


/// the non-planar NO top master (not reducible by the graph evaluator) against its literature eps-expansion (G-scheme
/// E_G = E / (eps G(1,1))^3 = 20 z5 + (68 z3^2 + 50 z6 - 80 z5) eps + (80 z5 - 272 z3^2 - 200 z6 + 204 z3 z4 + 450 z7) eps^2 + O(eps^3))
/// at eps = 1/160000 (truncation ~ 1e-14 relative)
#[test]
#[ignore = "three-loop NO family at w60, eps = 1/160000 (about 2 minutes in release)"]
fn no_top_master_against_its_expansion() {
    let (sys, r, _g) = run_fixture("P3L_NO_m0-solution-319999_80000.json", "P3L_NO", 60, 110);
    check_ranks(&r);
    let prec = 300;
    let z = |k: u32| Float::with_val(prec, Float::with_val(prec, k).zeta());
    let eps = Float::with_val(prec, Rational::from((1, 160000)));
    let c0 = Float::with_val(prec, &z(5) * 20u32);
    let c1 = Float::with_val(prec, Float::with_val(prec, z(3).square() * 68u32) + Float::with_val(prec, &z(6) * 50u32) - Float::with_val(prec, &z(5) * 80u32));
    let c2 = Float::with_val(prec, &z(5) * 80u32) - Float::with_val(prec, z(3).square() * 272u32) - Float::with_val(prec, &z(6) * 200u32)
        + Float::with_val(prec, Float::with_val(prec, &z(3) * &z(4)) * 204u32) + Float::with_val(prec, &z(7) * 450u32);
    let ser = Float::with_val(prec, &c0 + Float::with_val(prec, &c1 * &eps)) + Float::with_val(prec, &c2 * Float::with_val(prec, eps.square_ref()));
    // eps G(1,1) = Gamma(1+eps) Gamma(1-eps)^2 / Gamma(2-2eps)
    let one = Float::with_val(prec, 1);
    let g11 = Float::with_val(prec, Float::with_val(prec, &one + &eps).gamma() * Float::with_val(prec, &one - &eps).gamma().square()) / Float::with_val(prec, Float::with_val(prec, 2u32) - Float::with_val(prec, &eps * 2u32)).gamma();
    let j = sys.masters.iter().position(|m| m == &vec![1, 1, 1, 1, 1, 1, 1, 1, 0]).expect("NO corner is a master");
    let e_val = Float::with_val(prec, &r.moments[0].1[j]); // sum of indices even: E = I
    let eg = Float::with_val(prec, &e_val / Float::with_val(prec, g11.pow(3u32)));
    let d = digits(&eg, &ser);
    assert!(d > 12.0, "NO corner E_G {eg} vs its expansion {ser}: {d} digits");
}

fn kite_cfg(w: u32, moments: Vec<i32>) -> (EBlockSystem, Graph, PConfig, std::collections::HashMap<Vec<i64>, Vec<(Vec<i64>, String)>>) {
    let p = format!("{}/fixtures/pflow/P2L_m0-solution-199_50.json", env!("CARGO_MANIFEST_DIR"));
    let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap();
    let rs = RuleSet::from_json(&v).unwrap();
    let sys = EBlockSystem::from_rules_with(&rs, 0, &|c: &str| parse_coefficient(c, "z")).unwrap();
    let mut cfg = PConfig::new(w, w + 20, (w as i64 * 17) / 10, 2);
    cfg.moments = moments;
    (sys, graph_of("P2L"), cfg, rs.rules)
}

/// m = 2 is (a_e)_1 E[a + 1_e] (checked with a_e = 2), m = 0 is E[a_e -> 0] for a_e = 1, m = 0 with a_e = 2 vanishes, and a
/// non-master target from the moments (extra_value), all against the graph evaluator
#[test]
fn moments_m2_m0_and_extra_value() {
    let (sys, g, cfg, rules) = kite_cfg(40, vec![1, 2, 0]);
    let r = run_b(&sys, &g, 0, &cfg).unwrap();
    assert!(r.accepted, "{}", r.diag["acceptance"]["failures"]);
    let h = Rational::from(&sys.d0 / 2u32);
    let prec = r.moments[0].1[0].prec();
    let none = graph::Provided::default();
    let ev = |a: &[i64]| graph::massless(&g, a, &h, prec, &none).unwrap();
    let e_of = |m: i32, a: &[i64]| -> Float {
        let k = r.moments.iter().position(|x| x.0 == m).unwrap();
        let j = sys.masters.iter().position(|x| x == a).unwrap();
        let sgn = if a.iter().sum::<i64>() % 2 != 0 { -1 } else { 1 };
        Float::with_val(prec, &r.moments[k].1[j] * sgn)
    };
    // m = 2, a_e = 2: value = (2)_1 / 1! E[(3,1,0,0,1)]
    let v = Float::with_val(prec, e_of(2, &[2, 1, 0, 0, 1]) / 2u32);
    assert!(digits(&v, &ev(&[3, 1, 0, 0, 1])) > 30.0);
    // m = 2, a_e = 1: E[(2,1,1,1,0)]
    assert!(digits(&e_of(2, &[1, 1, 1, 1, 0]), &ev(&[2, 1, 1, 1, 0])) > 30.0);
    // m = 0, a_e = 1: (-1)^0 (1)_0 / 0! E[a_e -> 0]
    assert!(digits(&e_of(0, &[1, 1, 1, 1, 0]), &ev(&[0, 1, 1, 1, 0])) > 30.0);
    // m = 0, a_e = 2: n = -1, no integer term: vanishes
    assert!(e_of(0, &[2, 1, 0, 0, 1]).abs() < Float::with_val(prec, 1e-30) * ev(&[2, 1, 0, 0, 1]).abs());
    // a non-master target through its rule and the moments
    let efree: std::collections::HashMap<Vec<i64>, Float> = sys_efree(&sys, &rules, &g, &h, prec);
    let (x, _) = super::engine::extra_value(&rules, &sys, &[2, 0, 1, 1, 1], &r.moments, &efree, "z", prec).unwrap();
    assert!(digits(&x, &ev(&[2, 0, 1, 1, 1])) > 30.0);
}

/// massless values of the z-free masters appearing in the rules (graph evaluator)
fn sys_efree(sys: &EBlockSystem, rules: &std::collections::HashMap<Vec<i64>, Vec<(Vec<i64>, String)>>, g: &Graph, h: &Rational, prec: u32) -> std::collections::HashMap<Vec<i64>, Float> {
    let mut out = std::collections::HashMap::new();
    for terms in rules.values() {
        for (t, _) in terms {
            if t[sys.e_slot] <= 0 {
                if let Some(v) = graph::massless(g, t, h, prec, &graph::Provided::default()) {
                    out.insert(t.clone(), v);
                }
            }
        }
    }
    out
}

/// the provided table: a correct entry for a graph the evaluator reduces is accepted (cross-check), a wrong one fails the
/// acceptance; parsing refuses a float index, a wrong d0 and conflicting duplicates
#[test]
fn provided_table_checks() {
    // three-loop NO fixture: the one-loop bubble between the external vertices occurs in its regions
    let p = format!("{}/fixtures/pflow/P3L_NO_m0-solution-62811_15625.json", env!("CARGO_MANIFEST_DIR"));
    let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(p).unwrap()).unwrap();
    let rs = RuleSet::from_json(&v).unwrap();
    let sys = EBlockSystem::from_rules_with(&rs, 0, &|c: &str| parse_coefficient(c, "z")).unwrap();
    let g = graph_of("P3L_NO");
    let mut cfg = PConfig::new(30, 50, 60, 3);
    let prec = run::bits(cfg.wf);
    let h = Rational::from(&sys.d0 / 2u32);
    let one = Rational::from(1);
    let bubble = graph::g_bubble(&one, &one, &h, prec).unwrap();
    let tab = |val: &str, idx: serde_json::Value, d0: String| serde_json::json!({"schema": "pflow-provided.v1", "normalisation": "E, P^2=1", "d0": d0, "digits": 60,
        "values": [{"edges": [[0, 1, idx.clone()], [0, 1, 1]], "ext": [0, 1], "value": val}]});
    let good = tab(&bubble.to_string_radix(10, Some(80)), serde_json::json!(1), sys.d0.to_string());
    cfg.provided = graph::Provided::from_json(&good, prec, &sys.d0).unwrap();
    let r = run_b(&sys, &g, 0, &cfg).unwrap();
    assert!(r.accepted, "{}", r.diag["acceptance"]["failures"]);
    let bad = tab("1.0", serde_json::json!(1), sys.d0.to_string());
    cfg.provided = graph::Provided::from_json(&bad, prec, &sys.d0).unwrap();
    let r = run_b(&sys, &g, 0, &cfg).unwrap();
    assert!(!r.accepted, "a wrong provided value must fail the acceptance");
    assert!(r.diag["acceptance"]["failures"].to_string().contains("provided values = evaluator"));
    assert!(graph::Provided::from_json(&tab("1.0", serde_json::json!(1.0), sys.d0.to_string()), prec, &sys.d0).is_err(), "float index");
    assert!(graph::Provided::from_json(&tab("1.0", serde_json::json!(1), "3".into()), prec, &sys.d0).is_err(), "wrong d0");
    let mut dup = good.clone();
    dup["values"] = serde_json::json!([good["values"][0], {"edges": [[0, 1, 1], [0, 1, 1]], "ext": [0, 1], "value": "2.0"}]);
    assert!(graph::Provided::from_json(&dup, prec, &sys.d0).is_err(), "conflicting duplicate");
}

/// n_lo follows the rows: a configured nlo above the lowest row is lowered, the result is unchanged
#[test]
fn nlo_follows_the_rows() {
    let (sys, g, cfg, _) = kite_cfg(30, vec![1]);
    let r1 = run_b(&sys, &g, 0, &cfg).unwrap();
    let mut c2 = cfg.clone();
    c2.nlo = 5;
    let r2 = run_b(&sys, &g, 0, &c2).unwrap();
    assert!(r1.accepted && r2.accepted);
    for f in r2.diag["families"].as_array().unwrap() {
        assert!(f["n_lo"].as_i64().unwrap() < 5, "{f}");
    }
    for (a, b) in r1.moments[0].1.iter().zip(&r2.moments[0].1) {
        if !a.is_zero() {
            assert!(digits(a, b) > 25.0);
        }
    }
}

/// the family checkpoint fingerprint follows the DE content
#[test]
fn system_fingerprint_follows_the_content() {
    let (sys, _, _, _) = kite_cfg(30, vec![1]);
    let mut s2 = sys.clone();
    s2.terms[0].num.content += Rational::from(1);
    assert_ne!(super::engine::system_fingerprint(&sys), super::engine::system_fingerprint(&s2));
}

/// the CLI on the repo's three-loop DE fixture (it contains RustRed's identity rule for a master), and --two-precision together
/// with --workdir (one checkpoint subdirectory per pass)
#[test]
fn cli_on_the_3l_fixture_and_two_precision_workdir() {
    let dir = format!("{}/fixtures/pflow", env!("CARGO_MANIFEST_DIR"));
    let s = |x: &str| x.to_string();
    let v = super::cli::run(&[s("--de"), format!("{dir}/P3L_NO_m0-solution-62811_15625.json"), s("--graph"), format!("{dir}/graph-P3L_NO.json"),
        s("--e"), s("0"), s("--w"), s("30"), s("--wf"), s("50"), s("--m"), s("60")]).unwrap();
    assert_eq!(v["accepted"], true, "{}", v["diag"]["acceptance"]["failures"]);
    let tmp = std::env::temp_dir().join(format!("pflow-tp-{}", std::process::id()));
    let v = super::cli::run(&[s("--de"), format!("{dir}/P2L_m0-solution-199_50.json"), s("--graph"), format!("{dir}/graph-P2L.json"),
        s("--e"), s("0"), s("--w"), s("30"), s("--wf"), s("45"), s("--m"), s("60"), s("--two-precision"), s("5"), s("--workdir"), tmp.to_string_lossy().into_owned()]).unwrap();
    assert_eq!(v["accepted"], true);
    assert!(v["values"][0]["digits_two_precision"].as_f64().unwrap() > 15.0);
    let _ = std::fs::remove_dir_all(&tmp);
}

/// four-loop regression (NO1 family, eps = 0.009239, boundary = p-integral's exact reductions with t105 from the 3L mode):
/// accepted, and all 73 masters with a closed form agree to >= 40 digits at w45
#[test]
#[ignore = "four-loop NO1 family (about 30-60 minutes at 2 threads, 2 GB)"]
fn four_loop_no1_regression() {
    let dir = format!("{}/fixtures/pflow", env!("CARGO_MANIFEST_DIR"));
    let v: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(format!("{dir}/P4L_NO1_m0-solution-1990761_500000.json")).unwrap()).unwrap();
    let rs = RuleSet::from_json(&v).unwrap();
    let sys = EBlockSystem::from_rules_with(&rs, 0, &|c: &str| parse_coefficient(c, "z")).unwrap();
    let g = graph_of("P4L_NO1");
    let mut cfg = PConfig::new(45, 200, 90, 4);
    cfg.threads = 2;
    let pv: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(format!("{dir}/provided-P4L_NO1-1990761_500000.json")).unwrap()).unwrap();
    cfg.provided = graph::Provided::from_json(&pv, run::bits(cfg.wf), &sys.d0).unwrap();
    let r = run_b(&sys, &g, 0, &cfg).unwrap();
    check_ranks(&r);
    let ic = &r.diag["independent_checks"];
    assert!(ic["closed_forms"].as_u64().unwrap() >= 73, "{ic}");
    assert!(ic["worst_closed_form_log10"].as_f64().unwrap() < -40.0, "{ic}");
}
