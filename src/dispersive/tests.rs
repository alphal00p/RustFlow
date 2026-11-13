//! End-to-end tests of the dispersive pipeline on small systems with a closed-form dispersion integral.
//!
//! Two masters, lam = d0/2 (near-integer, eps = 1/150000 as at node a1 idx 022):
//!   D1(z) = z^lam (1 - z)^(-3),        d1(s) = s^lam (1 + s)^(-3)
//!   D2(z) = (z + 1/2) D1(z),           d2(s) = (1/2 - s) d1(s)
//! as the diagonal system A11 = lam/z - 3/(z - 1), A22 = A11 + 1/(z + 1/2) (singular point z = -1/2 on the cut,
//! non-branching for this solution: a detour is planned and the dropped imaginary part must stay at rounding level).
//! Soft family at 0: x^lam, dim 2; family at infinity: x^(-lam) with x = 1/z, i.e. z^lam z^(-n) (L = 1), as the
//! eta^(k d0/2 - n) families of the real system.
//!   I1 = int_0^inf d1/(s + 1) = B(lam + 1, 3 - lam),  I2 = B(lam + 1, 3 - lam)/2 - B(lam + 2, 2 - lam).
use super::boundary::BoundaryRow;
use super::poly::{self, QPoly, ZPoly};
use super::run::{self, Config};
use super::system::{EBlockSystem, Term};
use super::transport::TaylorRecurrence;
use rug::{Float, Rational};

fn q(v: &[Rational]) -> QPoly {
    v.to_vec()
}

fn beta(a: &Float, b: &Float, prec: u32) -> Float {
    let ga = a.clone().gamma();
    let gb = b.clone().gamma();
    let gab = Float::with_val(prec, a + b).gamma();
    Float::with_val(prec, &ga * &gb) / gab
}

fn two_master(d0: &Rational) -> EBlockSystem {
    let lam = Rational::from(d0 / 2);
    // A11 = lam/(1 + eta) - 3/eta = (-3 + (lam - 3) eta) / (eta + eta^2)
    let n11 = q(&[Rational::from(-3), Rational::from(&lam - Rational::from(3))]);
    let d11 = q(&[Rational::new(), Rational::from(1), Rational::from(1)]);
    // A22 = A11 + 1/(eta + 3/2)
    let d3 = q(&[Rational::from((3, 2)), Rational::from(1)]);
    let n22 = poly::add(&poly::mul(&n11, &d3), &d11);
    let d22 = poly::mul(&d11, &d3);
    EBlockSystem {
        d0: d0.clone(),
        e_slot: 0,
        masters: vec![vec![1, 0], vec![1, 1]],
        terms: vec![
            Term { i: 0, j: 0, num: ZPoly::from_q(&n11), den: ZPoly::from_q(&d11) },
            Term { i: 1, j: 1, num: ZPoly::from_q(&n22), den: ZPoly::from_q(&d22) },
        ],
    }
}

#[test]
fn two_master_dispersion_integral() {
    let d0 = Rational::from((299999, 75000));
    let sys = two_master(&d0);
    let dir = std::env::temp_dir().join(format!("w3a-test-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let cfg = Config { w: 30, wf: 70, loops: 1, nthreads: 2, workdir: dir.to_string_lossy().into(), ..Default::default() };
    let plan = run::plan(&sys, &cfg).unwrap();
    assert_eq!(plan.zero.dim, 2, "soft family dimension");
    assert_eq!(plan.inf[0].dim, 2, "infinity family dimension");
    assert!(plan.patches.iter().any(|p| !p.is_real()), "a detour around z = -1/2 is planned");
    let pf = run::bits(cfg.wf);
    let bnd = vec![
        BoundaryRow { comp: 0, n: 0, value: Float::with_val(pf, 1) },
        BoundaryRow { comp: 1, n: 0, value: Float::with_val(pf, 0.5) },
    ];
    let kernel = TaylorRecurrence::new(&sys, 2);
    let out = run::run_sample(&sys, &cfg, &plan, &bnd, &kernel, serde_json::json!({})).unwrap();
    let lam = Float::with_val(pf, &d0) / 2u32;
    let b1 = beta(&Float::with_val(pf, &lam + 1u32), &Float::with_val(pf, 3u32 - &lam), pf);
    let b2 = beta(&Float::with_val(pf, &lam + 2u32), &Float::with_val(pf, 2u32 - &lam), pf);
    let want = [b1.clone(), Float::with_val(pf, &b1 / 2u32) - b2];
    for (c, w) in want.iter().enumerate() {
        let got = Float::with_val(pf, Float::parse(out["I"][c]["I"].as_str().unwrap()).unwrap());
        let rel = Float::with_val(pf, (got - w) / w).abs();
        let digits = -rel.log10().to_f64();
        assert!(digits > 25.0, "master {c}: {digits:.1} digits; out {}", out);
    }
    assert!(out["matching"]["worst_rel_log10"].as_f64().unwrap() < -25.0, "{}", out["matching"]);
    assert!(out["boundary_fit"]["worst_rel_log10"].as_f64().unwrap() < -50.0, "{}", out["boundary_fit"]);
    assert!(out["transport"]["max_log10_im_over_re_after_detour"].as_f64().unwrap() < -25.0, "{}", out["transport"]);
    // resume: the checkpoints of the zero and transport stages give the identical result
    let again = run::run_sample(&sys, &cfg, &plan, &bnd, &kernel, serde_json::json!({})).unwrap();
    assert_eq!(again["I"], out["I"]);
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn modp_and_mpfr_ranks_agree_on_a_small_system() {
    use super::field::{ModP, Scalar};
    use super::frobenius;
    use super::laurent::{self, Point};
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
    let row_c0 = ((0 - spec.n_lo) as usize) * 2;
    let row_c1 = ((1 - spec.n_lo) as usize) * 2;
    let pv = spec.pivot_rows.iter().position(|&(c, n)| c == 0 && n == 0).unwrap();
    let r = Float::with_val(200, &basis[row_c1][pv] / &basis[row_c0][pv]);
    let want = Float::with_val(200, 3) / 4u32;
    assert!(Float::with_val(200, r - want).abs() < 1e-50);
}
