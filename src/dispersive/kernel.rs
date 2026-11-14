//! The dispersive mode's PRODUCTION transport: RustFlow's fixed-precision patch kernel ([`crate::patch_kernel`]).
//!
//! The e-block system ([`EBlockSystem`], exact Q[eta] numerators/denominators from the RustRed DE) enters the kernel as a
//! [`RationalSystem`] in `x = eta`; the planned patches ([`Patch`], path variable `z = 1 + eta`) as a [`PlannedPath`] with offset 1.
//! The kernel runs the dispersive configuration: real projection on real patches (the discontinuity is real; the dropped imaginary
//! part is recorded), Cauchy pole at `z = 1` (integrand `d(s)/(s + 1) ds` on `z = -s`), row-parallel recurrence bit-identical for
//! any thread count, exact checkpoints after every patch (resumed with any thread count).
//!
//! Sign: the kernel accumulates `+int y(z) / (1 - z) dz` per patch; the dispersive interface ([`TransportOutput::weighted_integral`])
//! is `-int_0^h y(z_c + t) / (1 - z_c - t) dt`, i.e. the opposite sign.
//!
//! The stand-in [`super::transport::TaylorRecurrence`] stays as the test reference (same operations, no checkpoints).
use super::path::Patch;
use super::system::EBlockSystem;
use super::transport::{PatchDiag, Transport, TransportOutput};
use crate::options::CancellationToken;
use crate::patch_kernel::checkpoint::CheckpointPolicy;
use crate::patch_kernel::{PatchSpec, PlannedPath, RationalSystem, TransportConfig, TransportOutcome};
use rug::{Complex, Integer, Rational};
use std::path::{Path, PathBuf};

/// How a checkpointed transport ended.
pub enum Checkpointed {
    Completed(TransportOutput),
    /// cancelled (signal): the kernel checkpoint holds the last completed patch
    Cancelled { next_patch: usize },
}

/// Control of a checkpointed transport stage.
pub struct StageControl {
    pub checkpoint: PathBuf,
    pub cancel: CancellationToken,
    /// write the checkpoint at least every this many seconds (and after the last patch / on cancel)
    pub every_seconds: f64,
    /// ... and after every this many patches (0 = by time only)
    pub every_patches: usize,
    /// stop (as if cancelled) once this patch index is reached: short-segment runs and kill tests
    pub stop_at_patch: Option<usize>,
}

/// e-block -> kernel system (exact; same canonical fingerprint as the harness-reader build of W2)
pub fn rational_system(sys: &EBlockSystem) -> Result<RationalSystem, String> {
    let ints = |p: &super::poly::ZPoly, mul: &Integer| -> Vec<Integer> { p.prim.iter().map(|c| Integer::from(c * mul)).collect() };
    let raw: Vec<(usize, usize, Vec<Integer>, Vec<Integer>)> = sys
        .terms
        .iter()
        .map(|t| {
            // (cn Pn) / (cd Pd) = (num(cn) den(cd) Pn) / (den(cn) num(cd) Pd)
            let a = Integer::from(t.num.content.numer() * t.den.content.denom());
            let b = Integer::from(t.num.content.denom() * t.den.content.numer());
            (t.i, t.j, ints(&t.num, &a), ints(&t.den, &b))
        })
        .collect();
    RationalSystem::from_entries(sys.dim(), "eta", raw).map_err(|e| e.to_string())
}

pub fn planned_path(patches: &[Patch]) -> PlannedPath {
    PlannedPath {
        offset: Rational::from(1),
        patches: patches.iter().map(|p| PatchSpec { centre: p.centre.clone(), step: p.step.clone(), order: p.order }).collect(),
    }
}

pub struct PatchKernel {
    pub system: RationalSystem,
    pub system_fingerprint: String,
    pub threads: usize,
    pub pin_threads: bool,
}

impl PatchKernel {
    pub fn new(sys: &EBlockSystem, threads: usize) -> Result<PatchKernel, String> {
        let system = rational_system(sys)?;
        let system_fingerprint = system.fingerprint();
        Ok(PatchKernel { system, system_fingerprint, threads: threads.max(1), pin_threads: true })
    }

    fn config(&self, prec: u32) -> TransportConfig {
        let mut c = TransportConfig::new(prec);
        // the dispersive configuration (REVIEW-A M5): real projection, Cauchy pole at z = 1
        c.cauchy_pole = Some(Rational::from(1));
        c.project_real_axis = true;
        c.threads = self.threads;
        c.pin_threads = self.pin_threads;
        c
    }

    fn output(&self, st: &crate::patch_kernel::TransportState) -> TransportOutput {
        let neg: Vec<Complex> = st.integral.as_ref().expect("dispersive transport has an integral").iter().map(|z| Complex::with_val(z.prec(), -z)).collect();
        let patches = st
            .records
            .iter()
            .map(|r| PatchDiag {
                seconds: r.seconds_setup + r.seconds_recurrence + r.seconds_finish,
                setup_seconds: r.seconds_setup,
                dropped_im_over_re_log10: r.dropped_imag_log10,
            })
            .collect();
        TransportOutput { end: st.y.clone(), weighted_integral: neg, patches }
    }

    /// kernel-level diagnostics of a finished transport (for out.json)
    pub fn summary(st: &crate::patch_kernel::TransportState) -> serde_json::Value {
        let worst = |it: &mut dyn Iterator<Item = f64>| it.fold(-1e9f64, |m, x| if x.is_nan() || m.is_nan() { f64::NAN } else { m.max(x) });
        let ov = worst(&mut st.records.iter().map(|r| r.overlap_log10));
        let tl = worst(&mut st.records.iter().map(|r| r.tail_log10));
        let dr = worst(&mut st.records.iter().map(|r| r.dropped_imag_log10));
        let setup: f64 = st.records.iter().map(|r| r.seconds_setup).sum();
        let rec: f64 = st.records.iter().map(|r| r.seconds_recurrence).sum();
        let fin: f64 = st.records.iter().map(|r| r.seconds_finish).sum();
        let nz = |x: f64| if x.is_finite() { serde_json::json!(x) } else { serde_json::json!(format!("{x}")) };
        serde_json::json!({"max_overlap_log10": nz(ov), "max_tail_log10": nz(tl), "max_dropped_imag_log10": nz(dr),
            "seconds_setup": setup, "seconds_recurrence": rec, "seconds_finish": fin, "identity_hash": st.identity_hash(),
            "patches_done": st.next_patch})
    }

    /// The full path with a checkpoint file: resumes the file when it exists (never starts over it: a file of another run is an
    /// error), writes it after every patch / `every_seconds`, flushes on cancel.
    pub fn run_checkpointed(&self, y0: &[Complex], patches: &[Patch], prec: u32, ctl: &StageControl) -> Result<(Checkpointed, serde_json::Value), String> {
        let path = planned_path(patches);
        let tr = crate::patch_kernel::Transport::new(&self.system, &path, self.config(prec))
            .map_err(|e| e.to_string())?
            .with_system_fingerprint(self.system_fingerprint.clone());
        let mut state = if ctl.checkpoint.exists() {
            let st = tr.resume(&ctl.checkpoint, y0, false).map_err(|e| format!("{}: {e} (refusing to start over it; use a fresh workdir)", ctl.checkpoint.display()))?;
            eprintln!("[dispersive] transport kernel: resumed {} at patch {} of {}", ctl.checkpoint.display(), st.next_patch, patches.len());
            st
        } else {
            tr.initial_state(y0).map_err(|e| e.to_string())?
        };
        let until = ctl.stop_at_patch.unwrap_or(patches.len()).min(patches.len());
        let policy = checkpoint_policy(&ctl.checkpoint, ctl.every_patches, ctl.every_seconds);
        let total = patches.len();
        let t0 = std::time::Instant::now();
        let outcome = tr
            .run(&mut state, until, y0, &ctl.cancel, Some(&policy), |s| {
                let r = s.records.last().unwrap();
                eprintln!(
                    "[dispersive] patch {}/{} order {} {} setup {:.1} s rec {:.1} s tail {:.1} overlap {:.1} ({:.0} s)",
                    r.index + 1,
                    total,
                    r.order,
                    if r.complex { "C" } else { "R" },
                    r.seconds_setup,
                    r.seconds_recurrence,
                    r.tail_log10,
                    r.overlap_log10,
                    t0.elapsed().as_secs_f64()
                );
            })
            .map_err(|e| e.to_string())?;
        let summary = Self::summary(&state);
        if state.next_patch == patches.len() {
            return Ok((Checkpointed::Completed(self.output(&state)), summary));
        }
        // stopped early (signal or stop_at_patch): make sure the checkpoint holds this state
        let _ = outcome == TransportOutcome::Completed;
        if !ctl.checkpoint.exists() {
            tr.write_checkpoint(&ctl.checkpoint, &state, y0).map_err(|e| e.to_string())?;
        }
        Ok((Checkpointed::Cancelled { next_patch: state.next_patch }, summary))
    }
}

/// the only place that builds the kernel's checkpoint policy (never `fresh`: a run resumes or refuses, REVIEW-A B1)
fn checkpoint_policy(file: &Path, every_patches: usize, every_seconds: f64) -> CheckpointPolicy {
    CheckpointPolicy { file: file.to_path_buf(), every_patches, every_seconds }
}

impl Transport for PatchKernel {
    fn id(&self) -> String {
        format!("{} (patch_kernel, dispersive config: real projection, Cauchy pole 1)", crate::patch_kernel::KERNEL_VERSION)
    }

    /// one-shot transport (no checkpoint): segments, tests
    fn transport(&self, _sys: &EBlockSystem, y0: &[Complex], patches: &[Patch], prec: u32) -> Result<TransportOutput, String> {
        let path = planned_path(patches);
        let tr = crate::patch_kernel::Transport::new(&self.system, &path, self.config(prec))
            .map_err(|e| e.to_string())?
            .with_system_fingerprint(self.system_fingerprint.clone());
        let mut st = tr.initial_state(y0).map_err(|e| e.to_string())?;
        let cancel = CancellationToken::default();
        tr.run(&mut st, patches.len(), y0, &cancel, None, |_| {}).map_err(|e| e.to_string())?;
        Ok(self.output(&st))
    }

    fn as_patch_kernel(&self) -> Option<&PatchKernel> {
        Some(self)
    }
}
