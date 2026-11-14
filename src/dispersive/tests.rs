//! Tests of the dispersive pipeline (seconds).
//!
//! (1) End to end on a 2-master system with a closed-form dispersion integral (detour, resume, thread-count identity, stale
//!     checkpoint refusal):
//!       D1(z) = z^lam (1 - z)^(-3),  d1(s) = s^lam (1 + s)^(-3);   D2 = (z + 1/2) D1,  d2(s) = (1/2 - s) d1(s)
//!       A11 = lam/z - 3/(z - 1), A22 = A11 + 1/(z + 1/2)  (axis singular point z = -1/2, non-branching for this solution)
//!       I1 = B(lam + 1, 3 - lam),  I2 = B(lam + 1, 3 - lam)/2 - B(lam + 2, 2 - lam).
//! (2) Pivoting on a NON-diagonal 3-master system with a near-resonant pair (REVIEW-B MAJOR 5): A11 as above, A21 = kappa/z,
//!     A22 = mu/z, mu = lam + 2 - 1/1000 (the lam + Z family of D2 has the small denominator n - 2 + 1/1000 at n = 2), A32 = 1/(z + 1/2),
//!     A33 = lam/z + 1/(1 - z). Production pivoting (complete, threshold) must solve the recursion to rounding and agree with each other
//!     and across two precisions.
//! (3) Soft boundary rows against the oracle (prep5 rows of the in-house W165 input at node a1 idx 022), J prefactor.
use super::boundary::{self, BoundaryRow};
use super::dispersion::{self, FamilyNum};
use super::field::{ModP, Scalar};
use super::frobenius::{self, Pivoting};
use super::laurent::{self, Point};
use super::poly::{self, QPoly, ZPoly};
use super::run::{self, Config};
use super::system::{EBlockSystem, Term};
use super::transport::TaylorRecurrence;
use rug::{Float, Rational};
use std::collections::BTreeMap;

fn beta(a: &Float, b: &Float, prec: u32) -> Float {
    let ga = a.clone().gamma();
    let gb = b.clone().gamma();
    let gab = Float::with_val(prec, a + b).gamma();
    Float::with_val(prec, &ga * &gb) / gab
}

fn term(i: usize, j: usize, n: &QPoly, d: &QPoly) -> Term {
    Term { i, j, num: ZPoly::from_q(n), den: ZPoly::from_q(d) }
}

/// eta-polynomials: A11 = lam/(1 + eta) - 3/eta = (-3 + (lam - 3) eta) / (eta + eta^2)
fn a11(lam: &Rational) -> (QPoly, QPoly) {
    (vec![Rational::from(-3), Rational::from(lam - Rational::from(3))], vec![Rational::new(), Rational::from(1), Rational::from(1)])
}

fn two_master(d0: &Rational) -> EBlockSystem {
    let lam = Rational::from(d0 / 2);
    let (n11, d11) = a11(&lam);
    let d3 = vec![Rational::from((3, 2)), Rational::from(1)];
    let n22 = poly::add(&poly::mul(&n11, &d3), &d11);
    let d22 = poly::mul(&d11, &d3);
    EBlockSystem { d0: d0.clone(), e_slot: 0, masters: vec![vec![1, 0], vec![1, 1]], terms: vec![term(0, 0, &n11, &d11), term(1, 1, &n22, &d22)] }
}

fn three_master(d0: &Rational) -> EBlockSystem {
    let lam = Rational::from(d0 / 2);
    let (n11, d11) = a11(&lam);
    let z = vec![Rational::from(1), Rational::from(1)]; // z = 1 + eta
    let kappa = Rational::from((7, 3));
    let mu = Rational::from(&lam + Rational::from(2)) - Rational::from((1, 1000));
    // A33 = lam/z + 1/(1 - z) = lam/(1 + eta) - 1/eta = (-1 + (lam - 1) eta) / (eta + eta^2)
    let n33 = vec![Rational::from(-1), Rational::from(&lam - Rational::from(1))];
    EBlockSystem {
        d0: d0.clone(),
        e_slot: 0,
        masters: vec![vec![1, 0], vec![1, 1], vec![1, 2]],
        terms: vec![
            term(0, 0, &n11, &d11),
            term(1, 0, &vec![kappa], &z),
            term(1, 1, &vec![mu], &z),
            term(2, 1, &vec![Rational::from(1)], &vec![Rational::from((3, 2)), Rational::from(1)]),
            term(2, 2, &n33, &d11),
        ],
    }
}

fn workdir(tag: &str) -> String {
    let d = std::env::temp_dir().join(format!("w3a-test-{tag}-{}", std::process::id()));
    std::fs::remove_dir_all(&d).ok();
    std::fs::create_dir_all(&d).unwrap();
    d.to_string_lossy().into()
}

#[test]
fn two_master_dispersion_integral() {
    let d0 = Rational::from((299999, 75000));
    let sys = two_master(&d0);
    let mut cfg = Config::new(30, 70, 24, 50, 1);
    cfg.nthreads = 2;
    cfg.workdir = workdir("a");
    let plan = run::plan(&sys, &cfg).unwrap();
    assert_eq!(plan.zero.dim, 2, "soft family dimension");
    assert_eq!(plan.inf[0].dim, 2, "infinity family dimension");
    assert!(plan.patches.iter().any(|p| !p.is_real()), "a detour around z = -1/2 is planned");
    let pf = run::bits(cfg.wf);
    let bnd = vec![BoundaryRow { comp: 0, n: 0, value: Float::with_val(pf, 1) }, BoundaryRow { comp: 1, n: 0, value: Float::with_val(pf, 0.5) }];
    let kernel = TaylorRecurrence::new(&sys, 2);
    let out = run::run_sample(&sys, &cfg, &plan, &bnd, &kernel, serde_json::json!({"boundary": "test"})).unwrap();
    let lam = Float::with_val(pf, &d0) / 2u32;
    let b1 = beta(&Float::with_val(pf, &lam + 1u32), &Float::with_val(pf, 3u32 - &lam), pf);
    let b2 = beta(&Float::with_val(pf, &lam + 2u32), &Float::with_val(pf, 2u32 - &lam), pf);
    let want = [b1.clone(), Float::with_val(pf, &b1 / 2u32) - b2];
    for (c, w) in want.iter().enumerate() {
        let got = Float::with_val(pf, Float::parse(out["I"][c]["I"].as_str().unwrap()).unwrap());
        let digits = -Float::with_val(pf, (got - w) / w).abs().log10().to_f64();
        assert!(digits > 25.0, "master {c}: {digits:.1} digits; out {}", out);
    }
    assert!(out["matching"]["worst_rel_log10"].as_f64().unwrap() < -25.0, "{}", out["matching"]);
    assert!(out["boundary_fit"]["worst_rel_log10"].as_f64().unwrap() < -50.0, "{}", out["boundary_fit"]);
    assert_eq!(out["transport"]["detour_check"], "pass");
    assert!(out["zero_diag"]["de_residual_log10"].as_f64().unwrap() < -50.0, "{}", out["zero_diag"]);
    // resume from the checkpoints: identical
    let again = run::run_sample(&sys, &cfg, &plan, &bnd, &kernel, serde_json::json!({"boundary": "test"})).unwrap();
    assert_eq!(again["I"], out["I"]);
    // another thread count, fresh workdir: identical output
    let mut cfg1 = cfg.clone();
    cfg1.nthreads = 1;
    cfg1.workdir = workdir("b");
    let k1 = TaylorRecurrence::new(&sys, 1);
    let one = run::run_sample(&sys, &cfg1, &plan, &bnd, &k1, serde_json::json!({"boundary": "test"})).unwrap();
    assert_eq!(one["I"], out["I"], "thread-count identity");
    // a changed boundary value in the SAME workdir must be refused (REVIEW-B BLOCKER 2), never silently resumed
    let mut bnd2 = bnd.clone();
    bnd2[1].value = Float::with_val(pf, 0.25);
    let err = run::run_sample(&sys, &cfg, &plan, &bnd2, &kernel, serde_json::json!({"boundary": "test"})).unwrap_err();
    assert!(err.contains("refusing to resume"), "{err}");
    // a changed pivoting mode too
    let mut cfg2 = cfg.clone();
    cfg2.pivoting = "complete".into();
    assert!(run::run_sample(&sys, &cfg2, &plan, &bnd, &kernel, serde_json::json!({})).is_err());
    for d in [&cfg.workdir, &cfg1.workdir] {
        std::fs::remove_dir_all(d).ok();
    }
}

#[test]
fn pivoting_on_a_coupled_near_resonant_system() {
    let d0 = Rational::from((299999, 75000));
    let sys = three_master(&d0);
    let e = sys.dim();
    let lam = Rational::from(&d0 / 2);
    let (n_lo, mm) = (-2i64, 16i64);
    let b = laurent::laurent::<ModP>(&sys, Point::Zero, mm + 8, 0, 1);
    let out = frobenius::frobenius::<ModP>(&b, e, &ModP::from_rational(&lam, 0), n_lo, mm, None, 0, 1);
    let spec = frobenius::family_spec(&out, e, n_lo, mm, &lam, &Rational::from((1, 4)), "soft").unwrap();
    assert_eq!(spec.dim, 2, "lam + Z family: the D1-driven chain and the homogeneous D3");
    let run_mode = |p: Pivoting, prec: u32| -> (FamilyNum, f64) {
        let bf = laurent::scaled(&laurent::laurent::<Float>(&sys, Point::Zero, mm + 8, prec, 1), &spec.rho, prec);
        let o = frobenius::frobenius_with::<Float>(&bf, e, &Float::with_val(prec, &lam), n_lo, mm, Some(&spec.ranks), Some(&spec.pivots), p, prec, 2);
        let fam = FamilyNum { spec: spec.clone(), basis: frobenius::stable_basis(&o, &spec, e, prec).unwrap() };
        let r = dispersion::family_residual(&bf, &fam, e, 8, prec);
        (fam, r)
    };
    let (fc, rc) = run_mode(Pivoting::Complete, 300);
    let (ft, rt) = run_mode(Pivoting::Threshold(30), 300);
    let (fh, _) = run_mode(Pivoting::Complete, 400);
    assert!(rc < -70.0 && rt < -70.0, "recursion residuals complete 1e{rc:.1} threshold 1e{rt:.1}");
    let agree = |a: &FamilyNum, b: &FamilyNum| -> f64 {
        let mut worst = -1e9f64;
        for (r1, r2) in a.basis.iter().zip(&b.basis) {
            for (x, y) in r1.iter().zip(r2) {
                if !x.is_zero() || !y.is_zero() {
                    let d = Float::with_val(400, Float::with_val(400, x) - y).abs();
                    let s = Float::with_val(400, x.clone().abs()).max(&Float::with_val(400, y.clone().abs()));
                    if !d.is_zero() {
                        worst = worst.max((d / s).log10().to_f64());
                    }
                }
            }
        }
        worst
    };
    assert!(agree(&fc, &ft) < -60.0, "complete vs threshold 1e{:.1}", agree(&fc, &ft));
    assert!(agree(&fc, &fh) < -60.0, "two precisions 1e{:.1}", agree(&fc, &fh));
}

#[test]
fn soft_rows_match_the_oracle() {
    // V (soft_V_euclid_nogamma of classes/16359/runs/boundary/V-a1-022.json, 90 digits) and the prep5 row c_{j, n'-a} of the in-house
    // W165 input $L/s/a1-p165r-022/inp.json (90 digits)
    let cases = [
        ("[1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 3, 1]|0", "0.0028967578278820641768824189925127861775046204516482970525577643962998058519443122100816", -1, "-434.5232592873785439743148962773236218552225249224974236726061715987876007413579595070805"),
        ("[3, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1]|1", "-0.006023254980581328016227594624163673733755135696591400022409429841653184564720840949014", -2, "-903.4991424219234734435814908911383873427396610370298371414463399669205660564697486737314"),
        ("[1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 2, 1, 0, 0]|0", "17578521268.791285028898803534141717119899678326069761047825945154173814569814634957033695", -1, "-2636836356027850.690967837004754549573941221583996742359396735871273472154012216948601424"),
    ];
    let prec = 400;
    let masters: Vec<Vec<i64>> = cases.iter().map(|c| serde_json::from_str(c.0.split('|').next().unwrap()).unwrap()).collect();
    let sys = EBlockSystem { d0: Rational::from((299999, 75000)), e_slot: 0, masters, terms: vec![] };
    let v: BTreeMap<String, String> = cases.iter().map(|c| (c.0.to_string(), c.1.to_string())).collect();
    let rows = boundary::rows_from_v(&sys, &v, 5, prec).unwrap();
    for (k, c) in cases.iter().enumerate() {
        let m: Vec<i64> = serde_json::from_str(c.0.split('|').next().unwrap()).unwrap();
        let r = rows.iter().find(|r| sys.masters[r.comp] == m).unwrap();
        assert_eq!(r.n, c.2, "case {k}: order");
        let want = Float::with_val(prec, Float::parse(c.3).unwrap());
        let d = -Float::with_val(prec, (r.value.clone() - &want) / &want).abs().log10().to_f64();
        assert!(d > 80.0, "case {k}: {d:.1} digits");
    }
    // an unknown key is an error, not a silent skip
    let mut bad = v.clone();
    bad.insert("[9, 9, 9, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1]|0".into(), "1".into());
    assert!(boundary::rows_from_v(&sys, &bad, 5, prec).is_err());
    // J = -(sin(pi d0/2)/pi) I: d0 = 3 -> 1/pi
    let p = dispersion::prefactor(&Rational::from(3), 200);
    let want = Float::with_val(200, 1) / Float::with_val(200, rug::float::Constant::Pi);
    assert!(Float::with_val(200, p - want).abs() < 1e-55);
}

#[test]
fn modp_and_mpfr_ranks_agree_on_a_small_system() {
    let d0 = Rational::from((299999, 75000));
    let sys = two_master(&d0);
    let lam = Rational::from(&d0 / 2);
    let b = laurent::laurent::<ModP>(&sys, Point::Zero, 20, 0, 1);
    let out = frobenius::frobenius::<ModP>(&b, 2, &ModP::from_rational(&lam, 0), -2, 12, None, 0, 1);
    let spec = frobenius::family_spec(&out, 2, -2, 12, &lam, &Rational::from((1, 4)), "soft").unwrap();
    assert_eq!((spec.dim, spec.first_order), (2, 0));
    let bf = laurent::laurent::<Float>(&sys, Point::Zero, 20, 200, 1);
    let bs = laurent::scaled(&bf, &spec.rho, 200);
    let of = frobenius::frobenius::<Float>(&bs, 2, &Float::with_val(200, &lam), -2, 12, Some(&spec.ranks), 200, 1);
    assert!(of.diag.min_rank_gap_log10 > 40.0, "rank gap {}", of.diag.min_rank_gap_log10);
    let basis = frobenius::stable_basis(&of, &spec, 2, 200).unwrap();
    // D1 = z^lam (1 - z)^(-3): c_1 / c_0 = 3 (scaled by rho = 1/4)
    let pv = spec.pivot_rows.iter().position(|&(c, n)| c == 0 && n == 0).unwrap();
    let r = Float::with_val(200, &basis[((1 - spec.n_lo) as usize) * 2][pv] / &basis[((0 - spec.n_lo) as usize) * 2][pv]);
    assert!(Float::with_val(200, r - Float::with_val(200, 3) / 4u32).abs() < 1e-50);
}

#[test]
fn pivot_modes_parse() {
    use frobenius::PivotMode;
    assert_eq!(PivotMode::parse(frobenius::DEFAULT_PIVOTING).unwrap(), PivotMode::Float(Pivoting::Threshold(30)));
    assert_eq!(PivotMode::parse("shadow").unwrap(), PivotMode::Shadow(Pivoting::Complete));
    assert_eq!(PivotMode::parse("shadow:threshold:20").unwrap(), PivotMode::Shadow(Pivoting::Threshold(20)));
    assert_eq!(PivotMode::parse("shadow:first").unwrap().label(), "shadow:first");
    assert!(PivotMode::parse("bogus").is_err());
}
