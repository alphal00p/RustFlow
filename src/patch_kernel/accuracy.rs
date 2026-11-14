//! Accuracy assessment of a fixed-precision transport by recomputation.
//!
//! * Two working precisions: the same path points, each run with its own
//!   truncation orders (`PlannedPath::with_orders`), compared entry by entry.
//!   Digits of agreement are a measured lower bound for the LOWER-precision
//!   run (rounding and truncation both differ between the runs).
//! * Patch overlap: consecutive patch series agree at the midpoint of the
//!   earlier patch (local truncation + rounding, per component).
//! * Series tail: relative size of the last Taylor term at the step.
//! * Tracked mode only: the precision Symbolica's tracking float has left.
use super::transport::{LOG10_NONE, TransportState};
use rug::{Complex, Float};
use serde::{Deserialize, Serialize};

/// Decimal digits of agreement of `a` and `b` (capped at `cap`).
pub fn agreement_digits(a: &Complex, b: &Complex, cap: f64) -> f64 {
    let p = a
        .prec()
        .0
        .max(b.prec().0)
        .max(a.prec().1)
        .max(b.prec().1)
        .max(64);
    let d = Complex::with_val(p, a - b);
    let dn = Float::with_val(64, d.abs_ref());
    let an = Float::with_val(64, a.abs_ref());
    let bn = Float::with_val(64, b.abs_ref());
    let scale = if an > bn { an } else { bn };
    if dn.is_zero() {
        return cap;
    }
    if scale.is_zero() {
        return 0.0;
    }
    (-(dn / scale).log10().to_f64()).clamp(0.0, cap)
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

/// Per-entry digits report.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DigitsReport {
    /// Working digits of the run under assessment and of the reference.
    pub digits: f64,
    pub reference_digits: f64,
    /// Agreement with the reference run, per entry (final vector).
    pub agreement: Vec<f64>,
    /// Agreement of the accumulated integral (if any).
    pub integral_agreement: Option<Vec<f64>>,
    /// -max overlap discrepancy per entry (digits of local consistency).
    pub overlap_digits: Vec<f64>,
    /// -max relative series tail per entry.
    pub tail_digits: Vec<f64>,
    pub min_agreement: f64,
    pub median_agreement: f64,
    pub min_overlap_digits: f64,
    pub min_integral_agreement: Option<f64>,
}

fn neg_log(v: f64) -> f64 {
    if v <= LOG10_NONE + 1.0 {
        f64::INFINITY
    } else {
        -v
    }
}

fn min_of(v: &[f64]) -> f64 {
    v.iter().copied().fold(f64::INFINITY, f64::min)
}

impl DigitsReport {
    /// Assess `run` (working digits `digits`) against `reference` (more digits).
    pub fn new(
        run: &TransportState,
        digits: f64,
        reference: &TransportState,
        reference_digits: f64,
    ) -> Self {
        let cap = reference_digits + 10.0;
        let agreement = compare_digits(&run.y, &reference.y, cap);
        let integral_agreement = match (&run.integral, &reference.integral) {
            (Some(a), Some(b)) => Some(compare_digits(a, b, cap)),
            _ => None,
        };
        let overlap_digits: Vec<f64> = run.overlap_log10.iter().map(|&v| neg_log(v)).collect();
        let tail_digits: Vec<f64> = run.tail_log10.iter().map(|&v| neg_log(v)).collect();
        let mut sorted = agreement.clone();
        sorted.sort_by(f64::total_cmp);
        let median = if sorted.is_empty() {
            f64::NAN
        } else {
            sorted[sorted.len() / 2]
        };
        Self {
            digits,
            reference_digits,
            min_agreement: min_of(&agreement),
            median_agreement: median,
            min_overlap_digits: min_of(&overlap_digits),
            min_integral_agreement: integral_agreement.as_deref().map(min_of),
            agreement,
            integral_agreement,
            overlap_digits,
            tail_digits,
        }
    }
}
