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
    PlannedPath {
        offset: Rational::new(),
        patches: segments.into_iter().flatten().collect(),
    }
}

fn real(v: &[f64], bits: u32) -> Vec<Complex> {
    v.iter().map(|&x| Complex::with_val(bits, (x, 0))).collect()
}

fn run_full(
    system: &RationalSystem,
    p: &PlannedPath,
    config: TransportConfig,
    y0: &[Complex],
) -> TransportState {
    let t = Transport::new(system, p, config).unwrap();
    let mut s = t.initial_state(y0).unwrap();
    let out = t
        .run(
            &mut s,
            usize::MAX,
            y0,
            &CancellationToken::default(),
            None,
            |_| {},
        )
        .unwrap();
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
    let s = run_full(
        &sys,
        &p,
        TransportConfig::generic(bits),
        &real(&[1.0], bits),
    );
    assert!(digits(&s.y[0], &Complex::with_val(bits, (10, 0))) > 38.0);
    // complex path 0 -> 2 + i/2 -> 3
    let p = path(vec![
        line((q("0"), q("0")), (q("2"), q("1/2")), 10, 70),
        line((q("2"), q("1/2")), (q("3"), q("0")), 5, 70),
    ]);
    let s = run_full(
        &sys,
        &p,
        TransportConfig::generic(bits),
        &real(&[1.0], bits),
    );
    assert!(digits(&s.y[0], &Complex::with_val(bits, (10, 0))) > 38.0);
}

#[test]
fn rotation_system_gives_cos_and_sin() {
    // y1' = y2, y2' = -y1 (constant entries)
    let sys = RationalSystem::from_entries(
        2,
        "x",
        vec![
            (0, 1, ints(&[1]), ints(&[1])),
            (1, 0, ints(&[-1]), ints(&[1])),
        ],
    )
    .unwrap();
    let bits = bits_for_digits(50);
    let p = path(vec![
        line((q("0"), q("0")), (q("1/2"), q("1/3")), 3, 60),
        line((q("1/2"), q("1/3")), (q("1"), q("0")), 3, 60),
    ]);
    let s = run_full(
        &sys,
        &p,
        TransportConfig::generic(bits),
        &real(&[1.0, 0.0], bits),
    );
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
        vec![(
            0,
            0,
            vec![a.numer().clone()],
            vec![Integer::from(-a.denom()), a.denom().clone()],
        )],
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
    let s = run_full(
        &sys,
        &p,
        TransportConfig::generic(bits),
        &real(&[1.0], bits),
    );
    let pi = Float::with_val(bits, rug::float::Constant::Pi);
    let phase = Float::with_val(bits, &pi * Float::with_val(bits, &a));
    let expect = Complex::with_val(bits, (phase.clone().cos(), -phase.sin()));
    assert!(
        digits(&s.y[0], &expect) > 38.0,
        "{}",
        digits(&s.y[0], &expect)
    );
}

#[test]
fn cauchy_weight_integral_matches_closed_form() {
    // y = 1 + x^2, weight 1/(a - x) with a = -1:
    // int_0^X (1 + x^2)/(-1 - x) dx = -(X^2/2 - X + 2 log(1 + X))
    let sys = quadratic();
    let bits = bits_for_digits(40);
    let p = path(vec![line((q("0"), q("0")), (q("2"), q("0")), 16, 70)]);
    let mut c = TransportConfig::generic(bits);
    c.cauchy_pole = Some(q("-1"));
    let s = run_full(&sys, &p, c, &real(&[1.0], bits));
    let x = Float::with_val(bits, 2);
    let exact =
        -(Float::with_val(bits, &x * &x) / 2u32 - &x + Float::with_val(bits, 3).ln() * 2u32);
    let got = &s.integral.as_ref().unwrap()[0];
    assert!(digits(got, &Complex::with_val(bits, (exact, 0))) > 38.0);
}

/// A dense-ish 9 x 9 system with several shared denominators and numerators
/// of various degrees (singular points at -3, -5/2, 4 and 2 +- i).
fn mixed() -> RationalSystem {
    let dens = [
        ints(&[3, 1]),
        ints(&[5, 2]),
        ints(&[-4, 1]),
        ints(&[5, -4, 1]),
        ints(&[15, 11, 2]),
    ];
    let mut raw = Vec::new();
    for i in 0..9i64 {
        for j in 0..9i64 {
            if (i * 7 + j * 3) % 4 == 0 && i != j {
                continue;
            }
            let n = ints(&[
                (i - j) % 5 + 1,
                (i * j) % 3 - 1,
                if (i + j) % 3 == 0 { 2 } else { 0 },
            ]);
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
    (0..n)
        .map(|j| Complex::with_val(bits, (Rational::from((j as i64 + 3, 7)), 0)))
        .collect()
}

#[test]
fn output_is_identical_for_any_thread_count() {
    let sys = mixed();
    let bits = bits_for_digits(60);
    let p = mixed_path(90);
    let y0 = ramp(9, bits);
    let mut hashes = Vec::new();
    for threads in [1, 2, 4, 8] {
        let mut c = TransportConfig::generic(bits);
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
    let mut c = TransportConfig::generic(bits);
    c.cauchy_pole = Some(q("-7"));
    c.threads = 3;
    let reference = run_full(&sys, &p, c.clone(), &y0);
    let dir = std::env::temp_dir().join(format!("rf-pk-{}-{}", std::process::id(), line!()));
    std::fs::create_dir_all(&dir).unwrap();
    let file = dir.join("state.json");
    let mut policy = CheckpointPolicy::new(file.clone());
    policy.every_patches = 2;
    // leg 1 (3 threads): cancelled at a patch boundary after patch 3
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
    // leg 2: resume with 1 and with 8 threads, each from the patch-3 checkpoint
    let saved = std::fs::read(&file).unwrap();
    for threads in [1, 8] {
        std::fs::write(&file, &saved).unwrap();
        let mut c2 = c.clone();
        c2.threads = threads;
        let t2 = Transport::new(&sys, &p, c2).unwrap();
        let mut s2 = t2.resume(&file, &y0, true).unwrap();
        assert_eq!(s2.next_patch, 3);
        let out = t2
            .run(
                &mut s2,
                usize::MAX,
                &y0,
                &CancellationToken::default(),
                Some(&policy),
                |_| {},
            )
            .unwrap();
        assert_eq!(out, TransportOutcome::Completed);
        assert_eq!(s2.identity_hash(), reference.identity_hash());
    }
    // deterministic cancellation INSIDE patch 5 (after order 16), then resume
    std::fs::remove_file(&file).unwrap();
    let mut c3 = c.clone();
    c3.debug_cancel_at = Some((5, 16));
    let t3 = Transport::new(&sys, &p, c3).unwrap();
    let mut s3 = t3.initial_state(&y0).unwrap();
    let out = t3
        .run(
            &mut s3,
            usize::MAX,
            &y0,
            &CancellationToken::default(),
            Some(&policy),
            |_| {},
        )
        .unwrap();
    assert_eq!(out, TransportOutcome::Cancelled);
    assert_eq!(s3.next_patch, 5);
    let mut c4 = c.clone();
    c4.threads = 2;
    let t4 = Transport::new(&sys, &p, c4).unwrap();
    let mut s4 = t4.resume(&file, &y0, true).unwrap();
    assert_eq!(s4.next_patch, 5);
    t4.run(
        &mut s4,
        usize::MAX,
        &y0,
        &CancellationToken::default(),
        Some(&policy),
        |_| {},
    )
    .unwrap();
    assert_eq!(s4.identity_hash(), reference.identity_hash());
    // a cancellation from another thread at an arbitrary time also resumes exactly
    let t5 = Transport::new(&sys, &p, c.clone()).unwrap();
    let mut s5 = t5.initial_state(&y0).unwrap();
    let token = CancellationToken::default();
    let tk = token.clone();
    let killer = std::thread::spawn(move || {
        std::thread::sleep(std::time::Duration::from_millis(30));
        tk.cancel();
    });
    let mut fresh = policy.clone();
    fresh.fresh = true;
    let _ = t5
        .run(&mut s5, usize::MAX, &y0, &token, Some(&fresh), |_| {})
        .unwrap();
    killer.join().unwrap();
    let mut s6 = t5.resume(&file, &y0, true).unwrap();
    t5.run(
        &mut s6,
        usize::MAX,
        &y0,
        &CancellationToken::default(),
        Some(&policy),
        |_| {},
    )
    .unwrap();
    assert_eq!(s6.identity_hash(), reference.identity_hash());
    // a different precision must refuse the checkpoint; a garbled one too
    let mut c6 = c.clone();
    c6.bits += 1;
    let t6 = Transport::new(&sys, &p, c6).unwrap();
    assert!(t6.resume(&file, &y0, true).is_err());
    std::fs::write(&file, b"{\"schema\": \"x\"").unwrap();
    assert!(t5.resume(&file, &y0, true).is_err());
    std::fs::remove_dir_all(&dir).unwrap();
}

/// Review A, B1: a fresh start never overwrites an existing checkpoint.
#[test]
fn fresh_start_does_not_overwrite_a_checkpoint() {
    let sys = mixed();
    let bits = bits_for_digits(30);
    let p = mixed_path(45);
    let y0 = ramp(9, bits);
    let dir = std::env::temp_dir().join(format!("rf-pk-{}-{}", std::process::id(), line!()));
    std::fs::create_dir_all(&dir).unwrap();
    let file = dir.join("c.json");
    let policy = CheckpointPolicy::new(file.clone());
    let t = Transport::new(&sys, &p, TransportConfig::generic(bits)).unwrap();
    let mut s = t.initial_state(&y0).unwrap();
    let out = t
        .run(
            &mut s,
            3,
            &y0,
            &CancellationToken::default(),
            Some(&policy),
            |_| {},
        )
        .unwrap();
    assert_eq!(out, TransportOutcome::Partial);
    let before = std::fs::read(&file).unwrap();
    // same run started again from patch 0: refused, file untouched
    let mut s0 = t.initial_state(&y0).unwrap();
    assert!(
        t.run(
            &mut s0,
            1,
            &y0,
            &CancellationToken::default(),
            Some(&policy),
            |_| {}
        )
        .is_err()
    );
    assert_eq!(std::fs::read(&file).unwrap(), before);
    // another run (other initial vector) on the same file: refused
    let other = ramp(9, bits).into_iter().rev().collect::<Vec<_>>();
    let mut so = t.initial_state(&other).unwrap();
    assert!(
        t.run(
            &mut so,
            1,
            &other,
            &CancellationToken::default(),
            Some(&policy),
            |_| {}
        )
        .is_err()
    );
    assert_eq!(std::fs::read(&file).unwrap(), before);
    // explicit fresh start is allowed
    let mut fresh = policy.clone();
    fresh.fresh = true;
    let mut s1 = t.initial_state(&y0).unwrap();
    t.run(
        &mut s1,
        1,
        &y0,
        &CancellationToken::default(),
        Some(&fresh),
        |_| {},
    )
    .unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
}
fn fingerprinted(
    sys: &RationalSystem,
    p: &PlannedPath,
    c: TransportConfig,
    y0: &[Complex],
) -> (TransportState, checkpoint::Fingerprint) {
    let t = Transport::new(sys, p, c).unwrap();
    let mut s = t.initial_state(y0).unwrap();
    t.run(
        &mut s,
        usize::MAX,
        y0,
        &CancellationToken::default(),
        None,
        |_| {},
    )
    .unwrap();
    let fp = t.fingerprint(y0);
    (s, fp)
}

#[test]
fn two_precision_comparison_and_overlap_diagnostics() {
    let sys = mixed();
    let lo = bits_for_digits(30);
    let hi = bits_for_digits(60);
    // the low run uses the lower order (its own truncation); same points
    let (a, fa) = fingerprinted(
        &sys,
        &mixed_path(45),
        TransportConfig::generic(lo),
        &ramp(9, lo),
    );
    let (b, fb) = fingerprinted(
        &sys,
        &mixed_path(90),
        TransportConfig::generic(hi),
        &ramp(9, hi),
    );
    let r = DigitsReport::new(&a, &fa, &b, &fb).unwrap();
    assert!(r.min_agreement > 20.0, "{:?}", r.agreement);
    // overlap consistency is bounded by the same accuracy
    assert!(r.min_overlap_digits > 20.0, "{:?}", r.overlap_digits);
    // review A, M3: not comparable -> refused
    assert!(DigitsReport::new(&b, &fb, &b, &fb).is_err(), "same bits");
    let (c, fc) = fingerprinted(
        &sys,
        &mixed_path(90),
        TransportConfig::generic(hi),
        &ramp(9, hi).into_iter().rev().collect::<Vec<_>>(),
    );
    assert!(
        DigitsReport::new(&a, &fa, &c, &fc).is_err(),
        "initial vector"
    );
    let mut partial = b.clone();
    partial.next_patch -= 1;
    assert!(
        DigitsReport::new(&a, &fa, &partial, &fb).is_err(),
        "incomplete"
    );
    let other_points = path(vec![line((q("0"), q("0")), (q("2"), q("0")), 8, 90)]);
    let (d, fd) = fingerprinted(
        &sys,
        &other_points,
        TransportConfig::generic(hi),
        &ramp(9, hi),
    );
    assert!(DigitsReport::new(&a, &fa, &d, &fd).is_err(), "path points");
}

/// Review A, M2: NaN is 0 digits; a pole on the integral weight fails fast.
#[test]
fn nan_counts_as_zero_digits_and_bad_weights_fail() {
    let bits = 200;
    let nan = Complex::with_val(bits, (Float::with_val(bits, rug::float::Special::Nan), 0));
    let one = Complex::with_val(bits, (1, 0));
    assert_eq!(accuracy::agreement_digits(&nan, &one, 70.0), 0.0);
    assert_eq!(accuracy::agreement_digits(&one, &one, 70.0), 70.0);
    let inf = Complex::with_val(
        bits,
        (Float::with_val(bits, rug::float::Special::Infinity), 0),
    );
    assert_eq!(accuracy::agreement_digits(&inf, &inf, 70.0), 0.0);
    // cauchy pole at a patch centre (review R2b): refused instead of an infinite integral
    let sys = quadratic();
    let p = path(vec![line((q("0"), q("0")), (q("2"), q("0")), 16, 40)]);
    let mut c = TransportConfig::generic(bits);
    c.cauchy_pole = Some(q("1"));
    let t = Transport::new(&sys, &p, c).unwrap();
    let mut s = t.initial_state(&[one.clone()]).unwrap();
    let before = s.identity_hash();
    assert!(
        t.run(
            &mut s,
            usize::MAX,
            &[one.clone()],
            &CancellationToken::default(),
            None,
            |_| {}
        )
        .is_err()
    );
    assert!(s.next_patch < 16);
    let _ = before;
}

/// Review A, M2(b): a patch whose disc contains a singular point diverges
/// and is refused (state unchanged).
#[test]
fn divergent_patch_fails_fast() {
    // y' = (1/3)/(x - 1) y: non-terminating series of radius 1 at x = 0
    let a = q("1/3");
    let sys = RationalSystem::from_entries(
        1,
        "x",
        vec![(
            0,
            0,
            vec![a.numer().clone()],
            vec![Integer::from(-a.denom()), a.denom().clone()],
        )],
    )
    .unwrap();
    let bits = bits_for_digits(40);
    let p = path(vec![line((q("0"), q("0")), (q("9/5"), q("0")), 1, 80)]);
    let t = Transport::new(&sys, &p, TransportConfig::generic(bits)).unwrap();
    let y0 = real(&[1.0], bits);
    let mut s = t.initial_state(&y0).unwrap();
    let h = s.identity_hash();
    let e = t.run(
        &mut s,
        usize::MAX,
        &y0,
        &CancellationToken::default(),
        None,
        |_| {},
    );
    assert!(e.is_err());
    assert_eq!(s.identity_hash(), h);
    assert_eq!(s.next_patch, 0);
}

/// Review A, M1 (R1): under real projection the overlap check is skipped at
/// the detour exit instead of raising a false alarm; elsewhere it holds.
#[test]
fn overlap_check_at_detour_exit_under_projection() {
    let a = q("1/3");
    let sys = RationalSystem::from_entries(
        1,
        "x",
        vec![(
            0,
            0,
            vec![a.numer().clone()],
            vec![Integer::from(-a.denom()), a.denom().clone()],
        )],
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
    let y0 = real(&[1.0], bits);
    for project in [false, true] {
        let mut c = TransportConfig::generic(bits);
        c.project_real_axis = project;
        let s = run_full(&sys, &p, c, &y0);
        let exit = &s.records[20];
        assert_eq!(exit.overlap_skipped, project);
        assert!(
            s.overlap_log10[0] < -35.0,
            "project {project}: {}",
            s.overlap_log10[0]
        );
        if project {
            // the real part of exp(-i pi/3) survives the projection
            let pi = Float::with_val(bits, rug::float::Constant::Pi);
            let re = Float::with_val(bits, &pi / 3u32).cos();
            assert!(digits(&s.y[0], &Complex::with_val(bits, (re, 0))) > 38.0);
            assert!(s.max_dropped_imag_log10() > -1.0);
        }
    }
}

/// Review A, m7: real arithmetic refuses complex inputs.
#[test]
fn real_scalars_refuse_complex_inputs() {
    let sys = quadratic();
    let plan = patch::KernelPlan::new(&sys);
    let bits = 100;
    let cancel = CancellationToken::default();
    let (z, h, hi) = (q("0"), q("1/8"), q("1/8"));
    let job = patch::PatchJob {
        system: &sys,
        plan: &plan,
        centre: (&z, &z),
        step: (&h, &hi),
        order: 20,
        bits,
        threads: 1,
        cauchy_u: None,
        probes: vec![],
        cancel: &cancel,
        pin_threads: false,
        debug_cancel_at_order: None,
    };
    assert!(patch::run_patch::<Float>(&job, &real(&[1.0], bits)).is_err());
    assert!(patch::run_patch::<Complex>(&job, &real(&[1.0], bits)).is_ok());
}

#[test]
fn singular_points_from_the_denominators() {
    let pts = mixed().singular_points().unwrap();
    let expect = [
        (-3.0, 0.0),
        (-2.5, 0.0),
        (2.0, -1.0),
        (2.0, 1.0),
        (4.0, 0.0),
    ];
    assert_eq!(pts.len(), expect.len(), "{pts:?}");
    for (a, b) in pts.iter().zip(expect) {
        assert!(
            (a.0 - b.0).abs() < 1e-12 && (a.1 - b.1).abs() < 1e-12,
            "{pts:?}"
        );
    }
}

#[test]
fn dispersive_configuration_is_explicit() {
    let c = TransportConfig::dispersive(300);
    assert!(c.project_real_axis);
    assert_eq!(c.cauchy_pole, Some(Rational::from(1)));
    let g = TransportConfig::generic(300);
    assert!(!g.project_real_axis && g.cauchy_pole.is_none());
}

#[test]
fn tracked_mode_agrees_with_fixed_mode() {
    let sys = mixed();
    let bits = bits_for_digits(40);
    let p = mixed_path(70);
    let y0 = ramp(9, bits);
    let f = run_full(&sys, &p, TransportConfig::generic(bits), &y0);
    let mut c = TransportConfig::generic(bits);
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
        vec![
            e("(3/7)/(2*pk_test_x - 5)"),
            e("(pk_test_x^2-1/2)/((pk_test_x+3)*(2*pk_test_x - 5))"),
        ],
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
    assert_eq!(
        a.fingerprint(),
        RationalSystem::from_json(a.to_json())
            .unwrap()
            .fingerprint()
    );
    let wrong = RationalSystem::from_entries(2, "x", vec![(0, 0, ints(&[0, 2]), ints(&[1, 0, 1]))])
        .unwrap();
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
