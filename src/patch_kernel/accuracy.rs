//! Accuracy assessment of a fixed-precision transport by recomputation.
//!
//! * Two working precisions: the same path points, each run with its own
//!   truncation orders (`PlannedPath::with_orders`), compared entry by entry.
//!   The agreement ESTIMATES the error of the lower-precision run: rounding
//!   and truncation differ between the runs, but anything common to both
//!   (system, initial vector, path points, the projection) is invisible to
//!   it. It is not a rigorous bound.
//! * Patch overlap: consecutive patch series agree at the midpoint of the
//!   earlier patch (local truncation + rounding, per component; skipped at
//!   detour exits under real projection).
//! * Series tail: relative size of the last Taylor terms at the step.
//! * Tracked mode only: the precision Symbolica's tracking float has left.
//!
//! NaN never counts as agreement: it is reported as 0 digits.
use super::checkpoint::Fingerprint;
use super::transport::{LOG10_NONE, TransportState};
use crate::{Error, Result};
use rug::{Complex, Float};
use serde::{Deserialize, Serialize};

fn digits_of(diff: Float, scale: Float, cap: f64) -> f64 {
    if diff.is_nan() || scale.is_nan() {
        return 0.0;
    }
    if diff.is_zero() {
        return cap;
    }
    if scale.is_zero() {
        return 0.0;
    }
    let d = -(diff / scale).log10().to_f64();
    if d.is_nan() { 0.0 } else { d.clamp(0.0, cap) }
}

fn prec_of(a: &Complex, b: &Complex) -> u32 {
    a.prec()
        .0
        .max(b.prec().0)
        .max(a.prec().1)
        .max(b.prec().1)
        .max(64)
}

/// Decimal digits of agreement of `a` and `b` relative to the larger modulus
/// (capped at `cap`; NaN or infinity -> 0).
pub fn agreement_digits(a: &Complex, b: &Complex, cap: f64) -> f64 {
    let d = Complex::with_val(prec_of(a, b), a - b);
    let an = Float::with_val(64, a.abs_ref());
    let bn = Float::with_val(64, b.abs_ref());
    if an.is_infinite() || bn.is_infinite() {
        return 0.0;
    }
    let scale = if an > bn { an } else { bn };
    digits_of(Float::with_val(64, d.abs_ref()), scale, cap)
}

/// Agreement of the real and of the imaginary parts separately, each
/// relative to its own size (a part that is zero in both gets `cap`).
pub fn part_digits(a: &Complex, b: &Complex, cap: f64) -> [f64; 2] {
    let p = prec_of(a, b);
    let one = |x: &Float, y: &Float| {
        if x.is_infinite() || y.is_infinite() {
            return 0.0;
        }
        let d = Float::with_val(p, x - y).abs();
        let s = Float::with_val(64, x.clone().abs()).max(&Float::with_val(64, y.clone().abs()));
        digits_of(Float::with_val(64, &d), s, cap)
    };
    [one(a.real(), b.real()), one(a.imag(), b.imag())]
}

/// Per-entry agreement of two vectors.
pub fn compare_digits(a: &[Complex], b: &[Complex], cap: f64) -> Vec<f64> {
    a.iter()
        .zip(b)
        .map(|(x, y)| agreement_digits(x, y, cap))
        .collect()
}

/// Digits Symbolica's tracking float claims for a value (its precision).
pub fn tracked_digits(z: &Complex) -> f64 {
    let (pr, pi) = z.prec();
    f64::from(pr.min(pi)) * std::f64::consts::LOG10_2
}

/// Working decimal digits of a run of `bits` (inverse of `bits_for_digits`).
pub fn working_digits(bits: u32) -> f64 {
    f64::from(bits) * std::f64::consts::LOG10_2 - 10.0
}

/// Per-entry digits report of a run against a higher-precision reference.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DigitsReport {
    /// Working digits of the run under assessment and of the reference.
    pub digits: f64,
    pub reference_digits: f64,
    /// Agreement with the reference run, per entry (final vector).
    pub agreement: Vec<f64>,
    /// Agreement of real and imaginary parts separately, per entry.
    pub agreement_parts: Vec<[f64; 2]>,
    /// Agreement of the accumulated integral (if any).
    pub integral_agreement: Option<Vec<f64>>,
    /// -max overlap discrepancy per entry (digits of local consistency).
    pub overlap_digits: Vec<f64>,
    /// -max relative series tail per entry.
    pub tail_digits: Vec<f64>,
    /// Tracked mode: the precision left in each entry (digits).
    pub tracked_digits: Option<Vec<f64>>,
    /// Largest relative imaginary part dropped by the real projection (log10).
    pub dropped_imag_log10: f64,
    /// Patches whose overlap check was skipped (detour exits, projection).
    pub overlap_skipped: usize,
    pub min_agreement: f64,
    pub median_agreement: f64,
    pub min_overlap_digits: f64,
    pub min_integral_agreement: Option<f64>,
}

fn neg_log(v: f64) -> f64 {
    if v.is_nan() {
        0.0
    } else if v <= LOG10_NONE + 1.0 {
        f64::INFINITY
    } else {
        -v
    }
}

/// Minimum where NaN counts as 0 digits.
fn min_of(v: &[f64]) -> f64 {
    v.iter()
        .map(|x| if x.is_nan() { 0.0 } else { *x })
        .fold(f64::INFINITY, f64::min)
}

impl DigitsReport {
    /// Assess `run` against `reference` (higher working precision). Refuses
    /// runs that are not comparable: different system, initial vector,
    /// Cauchy pole, projection or path points, an incomplete run, or equal
    /// working precision.
    pub fn new(
        run: &TransportState,
        run_fp: &Fingerprint,
        reference: &TransportState,
        reference_fp: &Fingerprint,
    ) -> Result<Self> {
        let mut why = Vec::new();
        if run_fp.system != reference_fp.system {
            why.push("system");
        }
        if run_fp.initial_64 != reference_fp.initial_64 {
            why.push("initial vector");
        }
        if run_fp.cauchy_pole != reference_fp.cauchy_pole {
            why.push("cauchy pole");
        }
        if run_fp.project_real_axis != reference_fp.project_real_axis {
            why.push("real projection");
        }
        if run_fp.path_points != reference_fp.path_points
            || run_fp.path_len != reference_fp.path_len
        {
            why.push("path points");
        }
        if run.next_patch != run_fp.path_len || reference.next_patch != reference_fp.path_len {
            why.push("incomplete run");
        }
        if run_fp.bits == reference_fp.bits {
            why.push("same working precision");
        }
        if run.y.len() != reference.y.len() {
            why.push("dimension");
        }
        if !why.is_empty() {
            return Err(Error::InvalidInput(format!(
                "runs are not comparable: {}",
                why.join(", ")
            )));
        }
        let digits = working_digits(run_fp.bits);
        let reference_digits = working_digits(reference_fp.bits);
        let cap = reference_digits + 10.0;
        let agreement = compare_digits(&run.y, &reference.y, cap);
        let agreement_parts = run
            .y
            .iter()
            .zip(&reference.y)
            .map(|(a, b)| part_digits(a, b, cap))
            .collect();
        let integral_agreement = match (&run.integral, &reference.integral) {
            (Some(a), Some(b)) => Some(compare_digits(a, b, cap)),
            _ => None,
        };
        let overlap_digits: Vec<f64> = run.overlap_log10.iter().map(|&v| neg_log(v)).collect();
        let tail_digits: Vec<f64> = run.tail_log10.iter().map(|&v| neg_log(v)).collect();
        let mut sorted: Vec<f64> = agreement
            .iter()
            .map(|x| if x.is_nan() { 0.0 } else { *x })
            .collect();
        sorted.sort_by(f64::total_cmp);
        let median = if sorted.is_empty() {
            0.0
        } else {
            sorted[sorted.len() / 2]
        };
        Ok(Self {
            digits,
            reference_digits,
            min_agreement: min_of(&agreement),
            median_agreement: median,
            min_overlap_digits: min_of(&overlap_digits),
            min_integral_agreement: integral_agreement.as_deref().map(min_of),
            tracked_digits: (run_fp.mode == "Tracked")
                .then(|| run.y.iter().map(tracked_digits).collect()),
            dropped_imag_log10: run.max_dropped_imag_log10(),
            overlap_skipped: run.records.iter().filter(|r| r.overlap_skipped).count(),
            agreement,
            agreement_parts,
            integral_agreement,
            overlap_digits,
            tail_digits,
        })
    }
}
