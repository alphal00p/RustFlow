//! Unit tests of the patch kernel on small systems with known solutions.
use super::checkpoint::CheckpointPolicy;
use super::path::{PatchSpec, PlannedPath};
use super::*;
use crate::options::CancellationToken;
use rug::{Complex, Float, Integer, Rational};
use symbolica::prelude::*;

fn ints(v: &[i64]) -> Vec<Integer> {
    v.iter().map(|&x| Integer::from(x)).collect()
}

fn q(s: &str) -> Rational {
    s.parse().unwrap()
}

/// Straight path from `a` to `b` (complex rationals) in `n` equal patches.
fn line(a: (Rational, Rational), b: (Rational, Rational), n: u32, order: usize) -> Vec<PatchSpec> {
    let step = (
        Rational::from(&b.0 - &a.0) / n,
        Rational::from(&b.1 - &a.1) / n,
    );
    (0..n)
        .map(|k| PatchSpec {
            centre: (
                Rational::from(&a.0 + Rational::from(&step.0 * k)),
                Rational::from(&a.1 + Rational::from(&step.1 * k)),
            ),
            step: step.clone(),
            order,
        })
        .collect()
}

fn path(segments: Vec<Vec<PatchSpec>>) -> PlannedPath {
    PlannedPath { offset: Rational::new(), patches: segments.into_iter().flatten().collect() }
}

fn real(v: &[f64], bits: u32) -> Vec<Complex> {
    v.iter().map(|&x| Complex::with_val(bits, (x, 0))).collect()
}

fn run_full(system: &RationalSystem, p: &PlannedPath, config: TransportConfig, y0: &[Complex]) -> TransportState {
    let t = Transport::new(system, p, config).unwrap();
    let mut s = t.initial_state(y0).unwrap();
    let out = t.run(&mut s, usize::MAX, y0, &CancellationToken::default(), None, |_| {}).unwrap();
    assert_eq!(out, TransportOutcome::Completed);
    s
}

fn digits(a: &Complex, b: &Complex) -> f64 {
    accuracy::agreement_digits(a, b, 1000.0)
}

/// y' = 2x/(1 + x^2) y has the polynomial solution y = 1 + x^2.
fn quadratic() -> RationalSystem {
    RationalSystem::from_entries(1, "x", vec![(0, 0, ints(&[0, 2]), ints(&[1, 0, 1]))]).unwrap()
}

#[test]
fn polynomial_solution_real_and_complex_paths() {
    let sys = quadratic();
    let bits = bits_for_digits(40);
    // real path 0 -> 3 (radius limited by +-i)
    let p = path(vec![line((q("0"), q("0")), (q("3"), q("0")), 12, 60)]);
    let s = run_full(&sys, &p, TransportConfig::new(bits), &real(&[1.0], bits));
    assert!(digits(&s.y[0], &Complex::with_val(bits, (10, 0))) > 38.0);
    // complex path 0 -> 2 + i/2 -> 3
    let p = path(vec![
        line((q("0"), q("0")), (q("2"), q("1/2")), 10, 70),
        line((q("2"), q("1/2")), (q("3"), q("0")), 5, 70),
    ]);
    let s = run_full(&sys, &p, TransportConfig::new(bits), &real(&[1.0], bits));
    assert!(digits(&s.y[0], &Complex::with_val(bits, (10, 0))) > 38.0);
}

#[test]
fn rotation_system_gives_cos_and_sin() {
    // y1' = y2, y2' = -y1 (constant entries)
    let sys = RationalSystem::from_entries(
        2,
        "x",
        vec![(0, 1, ints(&[1]), ints(&[1])), (1, 0, ints(&[-1]), ints(&[1]))],
    )
    .unwrap();
    let bits = bits_for_digits(50);
    let p = path(vec![
        line((q("0"), q("0")), (q("1/2"), q("1/3")), 3, 60),
        line((q("1/2"), q("1/3")), (q("1"), q("0")), 3, 60),
    ]);
    let s = run_full(&sys, &p, TransportConfig::new(bits), &real(&[1.0, 0.0], bits));
    let one = Float::with_val(bits, 1);
    let c = Complex::with_val(bits, (one.clone().cos(), 0));
    let sn = Complex::with_val(bits, (-one.sin(), 0));
    assert!(digits(&s.y[0], &c) > 48.0, "{}", digits(&s.y[0], &c));
    assert!(digits(&s.y[1], &sn) > 48.0);
}

#[test]
fn branch_point_detour_gives_the_monodromy_phase() {
    // y' = a/(x - 1) y, y(0) = 1: y = (1 - x)^a; above x = 1, y(2) = exp(-i pi a)
    let a = q("1/3");
    let sys = RationalSystem::from_entries(
        1,
        "x",
        vec![(0, 0, vec![a.numer().clone()], vec![Integer::from(-a.denom()), a.denom().clone()])],
    )
    .unwrap();
    let bits = bits_for_digits(40);
    let p = path(vec![
        line((q("0"), q("0")), (q("1/2"), q("0")), 4, 80),
        line((q("1/2"), q("0")), (q("1/2"), q("1/2")), 4, 80),
        line((q("1/2"), q("1/2")), (q("3/2"), q("1/2")), 8, 80),
        line((q("3/2"), q("1/2")), (q("3/2"), q("0")), 4, 80),
        line((q("3/2"), q("0")), (q("2"), q("0")), 4, 80),
    ]);
    let s = run_full(&sys, &p, TransportConfig::new(bits), &real(&[1.0], bits));
    let pi = Float::with_val(bits, rug::float::Constant::Pi);
    let phase = Float::with_val(bits, &pi * Float::with_val(bits, &a));
    let expect = Complex::with_val(bits, (phase.clone().cos(), -phase.sin()));
    assert!(digits(&s.y[0], &expect) > 38.0, "{}", digits(&s.y[0], &expect));
}

#[test]
fn cauchy_weight_integral_matches_closed_form() {
    // y = 1 + x^2, weight 1/(a - x) with a = -1:
    // int_0^X (1 + x^2)/(-1 - x) dx = -(X^2/2 - X + 2 log(1 + X))
    let sys = quadratic();
    let bits = bits_for_digits(40);
    let p = path(vec![line((q("0"), q("0")), (q("2"), q("0")), 16, 70)]);
    let mut c = TransportConfig::new(bits);
    c.cauchy_pole = Some(q("-1"));
    let s = run_full(&sys, &p, c, &real(&[1.0], bits));
    let x = Float::with_val(bits, 2);
    let exact = -(Float::with_val(bits, &x * &x) / 2u32 - &x + Float::with_val(bits, 3).ln() * 2u32);
    let got = &s.integral.as_ref().unwrap()[0];
    assert!(digits(got, &Complex::with_val(bits, (exact, 0))) > 38.0);
}

/// A dense-ish 9 x 9 system with several shared denominators and numerators
/// of various degrees (singular points at -3, -5/2, 4 and 2 +- i).
fn mixed() -> RationalSystem {
    let dens = [ints(&[3, 1]), ints(&[5, 2]), ints(&[-4, 1]), ints(&[5, -4, 1]), ints(&[15, 11, 2])];
    let mut raw = Vec::new();
    for i in 0..9i64 {
        for j in 0..9i64 {
            if (i * 7 + j * 3) % 4 == 0 && i != j {
                continue;
            }
            let n = ints(&[(i - j) % 5 + 1, (i * j) % 3 - 1, if (i + j) % 3 == 0 { 2 } else { 0 }]);
            let d = dens[((i + 2 * j) % 5) as usize].clone();
            raw.push((i as usize, j as usize, n, d));
        }
    }
    RationalSystem::from_entries(9, "x", raw).unwrap()
}

fn mixed_path(order: usize) -> PlannedPath {
    path(vec![
        line((q("0"), q("0")), (q("1/2"), q("0")), 2, order),
        line((q("1/2"), q("0")), (q("1"), q("1/4")), 2, order),
        line((q("1"), q("1/4")), (q("3/2"), q("0")), 2, order),
        line((q("3/2"), q("0")), (q("2"), q("0")), 2, order),
    ])
}

fn ramp(n: usize, bits: u32) -> Vec<Complex> {
    (0..n).map(|j| Complex::with_val(bits, (Rational::from((j as i64 + 3, 7)), 0))).collect()
}

#[test]
fn output_is_identical_for_any_thread_count() {
    let sys = mixed();
    let bits = bits_for_digits(60);
    let p = mixed_path(90);
    let y0 = ramp(9, bits);
    let mut hashes = Vec::new();
    for threads in [1, 2, 4, 8] {
        let mut c = TransportConfig::new(bits);
        c.threads = threads;
        c.cauchy_pole = Some(q("-7"));
        let s = run_full(&sys, &p, c, &y0);
        hashes.push(s.identity_hash());
    }
    assert!(hashes.windows(2).all(|w| w[0] == w[1]), "{hashes:?}");
}

#[test]
fn kill_and_resume_is_identical_also_with_other_threads() {
    let sys = mixed();
    let bits = bits_for_digits(60);
    let p = mixed_path(90);
    let y0 = ramp(9, bits);
    let mut c = TransportConfig::new(bits);
    c.cauchy_pole = Some(q("-7"));
    c.threads = 3;
    let reference = run_full(&sys, &p, c.clone(), &y0);
    let dir = std::env::temp_dir().join(format!("rf-pk-{}-{}", std::process::id(), line!()));
    std::fs::create_dir_all(&dir).unwrap();
    let file = dir.join("state.json");
    let policy = CheckpointPolicy { file: file.clone(), every_patches: 2, every_seconds: 0.0 };
    // first run: cancelled after the third patch (cancel at a patch boundary)
    let t = Transport::new(&sys, &p, c.clone()).unwrap();
    let mut s = t.initial_state(&y0).unwrap();
    let token = CancellationToken::default();
    let tk = token.clone();
    let out = t
        .run(&mut s, usize::MAX, &y0, &token, Some(&policy), |st| {
            if st.next_patch == 3 {
                tk.cancel();
            }
        })
        .unwrap();
    assert_eq!(out, TransportOutcome::Cancelled);
    assert_eq!(s.next_patch, 3);
    // resume with other thread counts, cancelling again mid-way (inside a patch)
    for threads in [1, 8] {
        let mut c2 = c.clone();
        c2.threads = threads;
        let t2 = Transport::new(&sys, &p, c2).unwrap();
        let mut s2 = t2.resume(&file, &y0, true).unwrap();
        assert!(s2.next_patch >= 3);
        let token = CancellationToken::default();
        let out = t2.run(&mut s2, usize::MAX, &y0, &token, Some(&policy), |_| {}).unwrap();
        assert_eq!(out, TransportOutcome::Completed);
        assert_eq!(s2.identity_hash(), reference.identity_hash());
        // reset the checkpoint to patch 3 for the next thread count
        let mut s3 = t2.initial_state(&y0).unwrap();
        let tk = CancellationToken::default();
        let tk2 = tk.clone();
        t2.run(&mut s3, usize::MAX, &y0, &tk, Some(&policy), |st| {
            if st.next_patch == 3 {
                tk2.cancel();
            }
        })
        .unwrap();
    }
    // a cancellation from another thread at an arbitrary time also resumes exactly
    let t3 = Transport::new(&sys, &p, c.clone()).unwrap();
    let mut s4 = t3.initial_state(&y0).unwrap();
    let token = CancellationToken::default();
    let tk = token.clone();
    let killer = std::thread::spawn(move || {
        std::thread::sleep(std::time::Duration::from_millis(30));
        tk.cancel();
    });
    let _ = t3.run(&mut s4, usize::MAX, &y0, &token, Some(&policy), |_| {}).unwrap();
    killer.join().unwrap();
    let mut s5 = t3.resume(&file, &y0, true).unwrap();
    t3.run(&mut s5, usize::MAX, &y0, &CancellationToken::default(), Some(&policy), |_| {}).unwrap();
    assert_eq!(s5.identity_hash(), reference.identity_hash());
    // a different precision must refuse the checkpoint
    let mut c4 = c.clone();
    c4.bits += 1;
    let t4 = Transport::new(&sys, &p, c4).unwrap();
    assert!(t4.resume(&file, &y0, true).is_err());
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn two_precision_comparison_and_overlap_diagnostics() {
    let sys = mixed();
    let lo = bits_for_digits(30);
    let hi = bits_for_digits(60);
    // the low run uses the lower order (its own truncation)
    let a = run_full(&sys, &mixed_path(45), TransportConfig::new(lo), &ramp(9, lo));
    let b = run_full(&sys, &mixed_path(90), TransportConfig::new(hi), &ramp(9, hi));
    let r = DigitsReport::new(&a, 30.0, &b, 60.0);
    assert!(r.min_agreement > 20.0, "{:?}", r.agreement);
    // overlap consistency is bounded by the same accuracy
    assert!(r.min_overlap_digits > 20.0, "{:?}", r.overlap_digits);
    let r2 = DigitsReport::new(&b, 60.0, &b, 60.0);
    assert!(r2.min_agreement >= 60.0);
}

#[test]
fn tracked_mode_agrees_with_fixed_mode() {
    let sys = mixed();
    let bits = bits_for_digits(40);
    let p = mixed_path(70);
    let y0 = ramp(9, bits);
    let f = run_full(&sys, &p, TransportConfig::new(bits), &y0);
    let mut c = TransportConfig::new(bits);
    c.mode = ArithmeticMode::Tracked;
    let t = run_full(&sys, &p, c, &y0);
    for (a, b) in f.y.iter().zip(&t.y) {
        let d = digits(a, b);
        let claimed = accuracy::tracked_digits(b);
        assert!(d > 30.0, "{d}");
        assert!(claimed > 10.0 && claimed <= 60.0, "{claimed}");
    }
}

#[test]
fn system_from_atoms_equals_system_from_entries() {
    let x = symbol!("pk_test::pk_test_x");
    let e = |s: &str| Atom::parse(s, "pk_test", Default::default()).unwrap();
    let m = vec![
        vec![e("2*pk_test_x/(1+pk_test_x^2)"), Atom::new()],
        vec![e("(3/7)/(2*pk_test_x - 5)"), e("(pk_test_x^2-1/2)/((pk_test_x+3)*(2*pk_test_x - 5))")],
    ];
    let a = RationalSystem::from_matrix_atoms(&m, x).unwrap();
    let b = RationalSystem::from_entries(
        2,
        "x",
        vec![
            (0, 0, ints(&[0, 2]), ints(&[1, 0, 1])),
            (1, 0, ints(&[3]), ints(&[-35, 14])),
            (1, 1, ints(&[-1, 0, 2]), ints(&[-30, 2, 4])),
        ],
    )
    .unwrap();
    assert!(a.differences(&b, &[0, 1]).is_empty());
    assert_eq!(a.fingerprint(), RationalSystem::from_json(a.to_json()).unwrap().fingerprint());
    let wrong = RationalSystem::from_entries(2, "x", vec![(0, 0, ints(&[0, 2]), ints(&[1, 0, 1]))]).unwrap();
    assert_eq!(a.differences(&wrong, &[0, 1]).len(), 2);
}

#[test]
fn path_json_round_trip_and_contiguity() {
    let p = mixed_path(20);
    let back = PlannedPath::from_json(p.to_json(), &Rational::new()).unwrap();
    assert_eq!(back, p);
    assert_eq!(back.fingerprint(), p.fingerprint());
    let mut bad = p.clone();
    bad.patches.swap(1, 2);
    assert!(bad.check_contiguous().is_err());
}
