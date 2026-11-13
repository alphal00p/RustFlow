//! Patch-by-patch transport along a planned path, with exact state,
//! diagnostics and checkpoints.
use super::checkpoint::{
    self, CHECKPOINT_SCHEMA, CheckpointPolicy, ExactComplex, Fingerprint, vector_from_exact,
    vector_to_exact,
};
use super::patch::{KernelPlan, PatchJob, PatchOutput, log10_abs, run_patch};
use super::path::PlannedPath;
use super::scalar::{ArithmeticMode, TrackedComplex, TrackedReal};
use super::system::RationalSystem;
use crate::options::CancellationToken;
use crate::{Error, Result};
use rug::{Complex, Float, Rational};
use serde::{Deserialize, Serialize};
use std::time::Instant;

/// Sentinel for "no difference / no information" in log10 diagnostics.
pub const LOG10_NONE: f64 = -1e9;

#[derive(Clone, Debug)]
pub struct TransportConfig {
    /// Working precision in bits (see [`super::bits_for_digits`]).
    pub bits: u32,
    pub mode: ArithmeticMode,
    /// Worker threads inside one patch (results do not depend on it).
    pub threads: usize,
    /// Accumulate `int y(z) / (a - z) dz` along the path (pole `a` in the
    /// path variable `z`); the dispersion integrand of W3 is
    /// `d(s)/(s + 1) ds = -d(z)/(1 - z) dz` on `z = -s`, i.e. `a = 1` with the
    /// opposite sign.
    pub cauchy_pole: Option<Rational>,
    /// On real patches keep the state real: imaginary parts left by complex
    /// detours are dropped and their largest relative size is recorded.
    pub project_real_axis: bool,
    /// Compare the series of consecutive patches at the midpoint of the
    /// earlier one (patch-overlap consistency).
    pub overlap_check: bool,
}

impl TransportConfig {
    pub fn new(bits: u32) -> Self {
        Self {
            bits,
            mode: ArithmeticMode::Fixed,
            threads: 1,
            cauchy_pole: None,
            project_real_axis: false,
            overlap_check: true,
        }
    }
}

/// Deterministic per-patch record (seconds are informational only).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PatchRecord {
    pub index: usize,
    pub order: usize,
    pub complex: bool,
    /// max_j log10 |overlap discrepancy_j| relative to the values.
    pub overlap_log10: f64,
    /// max_j relative series tail log10 |y_{N-1} h^{N-1}| / |y|.
    pub tail_log10: f64,
    /// max_j log10 |Im y_j| / |Re y_j| dropped before this (real) patch.
    pub dropped_imag_log10: f64,
    pub seconds_setup: f64,
    pub seconds_recurrence: f64,
    pub seconds_finish: f64,
}

#[derive(Clone, Debug)]
pub struct Probe {
    /// Absolute point (path variable).
    pub point: (Rational, Rational),
    pub values: Vec<Complex>,
}

/// State after `next_patch` completed patches. Everything except the
/// `seconds_*` fields of the records is bit-reproducible.
#[derive(Clone, Debug)]
pub struct TransportState {
    pub next_patch: usize,
    pub y: Vec<Complex>,
    pub integral: Option<Vec<Complex>>,
    pub probe: Option<Probe>,
    pub records: Vec<PatchRecord>,
    /// Per component: max over patches of the overlap discrepancy (log10).
    pub overlap_log10: Vec<f64>,
    /// Per component: max over patches of the relative series tail (log10).
    pub tail_log10: Vec<f64>,
}

#[derive(Serialize, Deserialize)]
struct StateFile {
    schema: String,
    fingerprint: Fingerprint,
    next_patch: usize,
    y: Vec<ExactComplex>,
    integral: Option<Vec<ExactComplex>>,
    probe: Option<((String, String), Vec<ExactComplex>)>,
    records: Vec<PatchRecord>,
    overlap_log10: Vec<f64>,
    tail_log10: Vec<f64>,
}

impl TransportState {
    pub fn to_json(&self, fingerprint: &Fingerprint) -> serde_json::Value {
        serde_json::to_value(StateFile {
            schema: CHECKPOINT_SCHEMA.into(),
            fingerprint: fingerprint.clone(),
            next_patch: self.next_patch,
            y: vector_to_exact(&self.y),
            integral: self.integral.as_deref().map(vector_to_exact),
            probe: self.probe.as_ref().map(|p| {
                ((p.point.0.to_string(), p.point.1.to_string()), vector_to_exact(&p.values))
            }),
            records: self.records.clone(),
            overlap_log10: self.overlap_log10.clone(),
            tail_log10: self.tail_log10.clone(),
        })
        .expect("state serialises")
    }

    /// Parse a checkpoint; returns the state and its fingerprint.
    pub fn from_json(value: serde_json::Value) -> Result<(Self, Fingerprint)> {
        let f: StateFile = serde_json::from_value(value)
            .map_err(|e| Error::InvalidInput(format!("checkpoint: {e}")))?;
        if f.schema != CHECKPOINT_SCHEMA {
            return Err(Error::InvalidInput(format!("checkpoint schema {}", f.schema)));
        }
        let q = |s: &str| {
            s.parse::<Rational>().map_err(|e| Error::InvalidInput(format!("rational {s:?}: {e}")))
        };
        let probe = match &f.probe {
            Some(((a, b), v)) => Some(Probe { point: (q(a)?, q(b)?), values: vector_from_exact(v)? }),
            None => None,
        };
        Ok((
            Self {
                next_patch: f.next_patch,
                y: vector_from_exact(&f.y)?,
                integral: f.integral.as_deref().map(vector_from_exact).transpose()?,
                probe,
                records: f.records,
                overlap_log10: f.overlap_log10,
                tail_log10: f.tail_log10,
            },
            f.fingerprint,
        ))
    }

    /// Exact identity of the numerical state (vector, integral, probe and the
    /// deterministic diagnostics), for thread/resume identity checks.
    pub fn identity_hash(&self) -> String {
        let mut h = blake3::Hasher::new();
        h.update(checkpoint::vector_hash(&self.y).as_bytes());
        if let Some(i) = &self.integral {
            h.update(checkpoint::vector_hash(i).as_bytes());
        }
        if let Some(p) = &self.probe {
            h.update(format!("{}+{}i", p.point.0, p.point.1).as_bytes());
            h.update(checkpoint::vector_hash(&p.values).as_bytes());
        }
        for r in &self.records {
            h.update(
                format!(
                    "{}:{}:{}:{:e}:{:e}:{:e}|",
                    r.index, r.order, r.complex, r.overlap_log10, r.tail_log10, r.dropped_imag_log10
                )
                .as_bytes(),
            );
        }
        for v in self.overlap_log10.iter().chain(&self.tail_log10) {
            h.update(&v.to_le_bytes());
        }
        h.finalize().to_hex().to_string()
    }
}

/// How a run ended.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TransportOutcome {
    /// All requested patches done.
    Completed,
    /// Cancelled (signal); the state is the last completed patch and was
    /// flushed to the checkpoint file if one is configured.
    Cancelled,
}

pub struct Transport<'a> {
    pub system: &'a RationalSystem,
    pub path: &'a PlannedPath,
    pub config: TransportConfig,
    plan: KernelPlan,
    system_fingerprint: String,
    binary: String,
}

impl<'a> Transport<'a> {
    pub fn new(system: &'a RationalSystem, path: &'a PlannedPath, config: TransportConfig) -> Result<Self> {
        path.check_contiguous()?;
        if config.bits < 16 {
            return Err(Error::InvalidInput("working precision too small".into()));
        }
        Ok(Self {
            plan: KernelPlan::new(system),
            system_fingerprint: system.fingerprint(),
            binary: checkpoint::binary_hash(),
            system,
            path,
            config,
        })
    }

    /// Use a precomputed system fingerprint (saves hashing a large system).
    pub fn with_system_fingerprint(mut self, fingerprint: String) -> Self {
        self.system_fingerprint = fingerprint;
        self
    }

    pub fn fingerprint(&self, initial: &[Complex]) -> Fingerprint {
        Fingerprint {
            kernel_version: super::KERNEL_VERSION.into(),
            system: self.system_fingerprint.clone(),
            path: self.path.fingerprint(),
            bits: self.config.bits,
            mode: format!("{:?}", self.config.mode),
            cauchy_pole: self.config.cauchy_pole.as_ref().map(|a| a.to_string()),
            project_real_axis: self.config.project_real_axis,
            overlap_check: self.config.overlap_check,
            initial: checkpoint::vector_hash(initial),
            binary: self.binary.clone(),
        }
    }

    /// Initial state from exact values (rounded to the working precision).
    pub fn initial_state(&self, y0: &[Complex]) -> Result<TransportState> {
        if y0.len() != self.system.dim {
            return Err(Error::InvalidInput("initial vector dimension".into()));
        }
        let bits = self.config.bits;
        Ok(TransportState {
            next_patch: 0,
            y: y0.iter().map(|z| Complex::with_val(bits, z)).collect(),
            integral: self
                .config
                .cauchy_pole
                .as_ref()
                .map(|_| vec![Complex::new(bits); self.system.dim]),
            probe: None,
            records: Vec::new(),
            overlap_log10: vec![LOG10_NONE; self.system.dim],
            tail_log10: vec![LOG10_NONE; self.system.dim],
        })
    }

    /// Load a checkpoint and check it against this transport and `initial`.
    pub fn resume(&self, file: &std::path::Path, initial: &[Complex], check_binary: bool) -> Result<TransportState> {
        let bytes = std::fs::read(file)?;
        let value: serde_json::Value = serde_json::from_slice(&bytes)
            .map_err(|e| Error::InvalidInput(format!("checkpoint JSON: {e}")))?;
        let (state, fp) = TransportState::from_json(value)?;
        let mine = self.fingerprint(initial);
        let bad = mine.mismatches(&fp, check_binary);
        if !bad.is_empty() {
            return Err(Error::InvalidInput(format!("checkpoint fingerprint mismatch: {bad:?}")));
        }
        if state.next_patch > self.path.patches.len() || state.y.len() != self.system.dim {
            return Err(Error::InvalidInput("checkpoint does not fit the path/system".into()));
        }
        Ok(state)
    }

    pub fn write_checkpoint(&self, file: &std::path::Path, state: &TransportState, initial: &[Complex]) -> Result<()> {
        let json = state.to_json(&self.fingerprint(initial));
        checkpoint::write_atomic(file, serde_json::to_string(&json).expect("json").as_bytes())
    }

    /// Transport patches `state.next_patch .. until`. Checkpoints according to
    /// `policy` (and on cancellation); `on_patch` sees every completed patch.
    pub fn run(
        &self,
        state: &mut TransportState,
        until: usize,
        initial: &[Complex],
        cancel: &CancellationToken,
        policy: Option<&CheckpointPolicy>,
        mut on_patch: impl FnMut(&TransportState),
    ) -> Result<TransportOutcome> {
        let until = until.min(self.path.patches.len());
        let mut last_write = Instant::now();
        let mut since = 0usize;
        while state.next_patch < until {
            if cancel.is_cancelled() {
                break;
            }
            match self.step(state, cancel) {
                Ok(()) => {}
                Err(Error::Cancelled) => break,
                Err(e) => return Err(e),
            }
            since += 1;
            on_patch(state);
            if let Some(p) = policy {
                let due = (p.every_patches > 0 && since >= p.every_patches)
                    || (p.every_seconds > 0.0 && last_write.elapsed().as_secs_f64() >= p.every_seconds)
                    || state.next_patch == until;
                if due {
                    self.write_checkpoint(&p.file, state, initial)?;
                    last_write = Instant::now();
                    since = 0;
                }
            }
        }
        if state.next_patch < until {
            if let Some(p) = policy {
                self.write_checkpoint(&p.file, state, initial)?;
            }
            return Ok(TransportOutcome::Cancelled);
        }
        Ok(TransportOutcome::Completed)
    }

    /// One patch; on error (including cancellation) `state` is unchanged.
    pub fn step(&self, state: &mut TransportState, cancel: &CancellationToken) -> Result<()> {
        let k = state.next_patch;
        let spec = &self.path.patches[k];
        let bits = self.config.bits;
        let offset = &self.path.offset;
        let centre_x = (Rational::from(&spec.centre.0 - offset), spec.centre.1.clone());
        let has_imag = state.y.iter().any(|z| !z.imag().is_zero());
        let real_patch = spec.is_real();
        let complex = !real_patch || (has_imag && !self.config.project_real_axis);
        // state handed to the patch
        let mut dropped = LOG10_NONE;
        let y0: Vec<Complex> = if complex {
            state.y.clone()
        } else {
            state
                .y
                .iter()
                .map(|z| {
                    if !z.imag().is_zero() && !z.real().is_zero() {
                        let r = log10_abs(&Complex::with_val(64, (z.imag(), 0)))
                            - log10_abs(&Complex::with_val(64, (z.real(), 0)));
                        dropped = dropped.max(r);
                    }
                    let (pr, pi) = z.prec();
                    Complex::with_val((pr, pi), (z.real(), 0))
                })
                .collect()
        };
        let mut probes = Vec::new();
        let check = self.config.overlap_check && state.probe.is_some();
        if check {
            let p = state.probe.as_ref().unwrap();
            probes.push((
                Rational::from(&p.point.0 - &spec.centre.0),
                Rational::from(&p.point.1 - &spec.centre.1),
            ));
        }
        let half = (Rational::from(&spec.step.0 / 2u32), Rational::from(&spec.step.1 / 2u32));
        if self.config.overlap_check {
            probes.push(half.clone());
        }
        let job = PatchJob {
            system: self.system,
            plan: &self.plan,
            centre: (&centre_x.0, &centre_x.1),
            step: (&spec.step.0, &spec.step.1),
            order: spec.order,
            bits,
            threads: self.config.threads,
            cauchy_u: self.config.cauchy_pole.as_ref().map(|a| {
                (Rational::from(a - &spec.centre.0), Rational::from(-&spec.centre.1))
            }),
            probes,
            cancel,
        };
        let out: PatchOutput = match (self.config.mode, complex) {
            (ArithmeticMode::Fixed, false) => run_patch::<Float>(&job, &y0)?,
            (ArithmeticMode::Fixed, true) => run_patch::<Complex>(&job, &y0)?,
            (ArithmeticMode::Tracked, false) => run_patch::<TrackedReal>(&job, &y0)?,
            (ArithmeticMode::Tracked, true) => run_patch::<TrackedComplex>(&job, &y0)?,
        };
        // ---- commit ----
        let mut overlap_max = LOG10_NONE;
        let mut probe_vals = out.probes.into_iter();
        if check {
            let mine = probe_vals.next().expect("overlap probe");
            let prev = &state.probe.as_ref().unwrap().values;
            for (j, (a, b)) in mine.iter().zip(prev).enumerate() {
                let p = a.prec().0.max(b.prec().0).max(a.prec().1).max(b.prec().1);
                let d = Complex::with_val(p, a - b);
                let scale = log10_abs(a).max(log10_abs(b));
                let r = if d.real().is_zero() && d.imag().is_zero() || scale < -1e8 {
                    LOG10_NONE
                } else {
                    log10_abs(&d) - scale
                };
                overlap_max = overlap_max.max(r);
                state.overlap_log10[j] = state.overlap_log10[j].max(r);
            }
        }
        state.probe = if self.config.overlap_check {
            Some(Probe {
                point: (
                    Rational::from(&spec.centre.0 + &half.0),
                    Rational::from(&spec.centre.1 + &half.1),
                ),
                values: probe_vals.next().expect("midpoint probe"),
            })
        } else {
            None
        };
        let mut tail_max = LOG10_NONE;
        for (j, t) in out.tail_log10.iter().enumerate() {
            tail_max = tail_max.max(*t);
            state.tail_log10[j] = state.tail_log10[j].max(*t);
        }
        if let (Some(acc), Some(add)) = (state.integral.as_mut(), out.integral.as_ref()) {
            for (a, b) in acc.iter_mut().zip(add) {
                let p = a.prec();
                let sum = Complex::with_val(p, &*a + b);
                *a = sum;
            }
        }
        state.y = out.value;
        state.records.push(PatchRecord {
            index: k,
            order: spec.order,
            complex,
            overlap_log10: overlap_max,
            tail_log10: tail_max,
            dropped_imag_log10: dropped,
            seconds_setup: out.seconds_setup,
            seconds_recurrence: out.seconds_recurrence,
            seconds_finish: out.seconds_finish,
        });
        state.next_patch = k + 1;
        Ok(())
    }
}
