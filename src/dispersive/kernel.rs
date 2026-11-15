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
    /// Override of the kernel's series-tail bound (log10, per component relative tail); `None` = the kernel default
    /// -(working decimal digits)/2. Not part of any fingerprint: it only decides whether a patch is accepted.
    pub max_tail_log10: Option<f64>,
    /// Tail heal ([`TAIL_HEAL_POLICY`]): a patch that fails the tail bound is retried at higher Taylor orders (same step, same
    /// precision) instead of failing the run. Default on.
    pub tail_heal: bool,
}

/// Self-healing of the series-tail guard (W7-TAIL, 15 Nov 2026). In the last real stretch of the dispersive path (|z| >> the
/// singular points, steps |z|/3) the components of the discontinuity spread over ~110 orders of magnitude; the per-component
/// relative tail of the tiniest components is a TRUNCATION tail driven by the dominant ones (it follows the Taylor order, ~0.8
/// digits per term at node 022, not the precision). The planned order (W + 8 digits for the nearest singular point) then
/// leaves a relative tail above -(digits)/2 in components 1e-50..1e-80 of the vector at the lower precision of a two-precision
/// run. Remedy: retry that ONE patch with more Taylor terms (same centre, step and precision; the weight and probe series as in
/// the kernel) and accept the first retry whose largest per-component relative tail is within the same bound. Patches that pass
/// at the planned order are untouched (bit-identical); a healed patch is recorded in the kernel state with its order (record
/// order != planned order), reported in out.json (`tail_heals`) and replayed identically on resume. If no retry passes, the run
/// fails as before (the guard stays).
pub const TAIL_HEAL_POLICY: &str = "tail-heal-1: retry a failing patch at orders ceil(5N/4), ceil(3N/2), 2N (same step and precision); accept the first within the tail bound";

/// retry orders of the tail heal for a patch of planned order `n`
pub fn tail_heal_orders(n: usize) -> [usize; 3] {
    [(5 * n).div_ceil(4), (3 * n).div_ceil(2), 2 * n]
}

fn is_tail_failure(e: &crate::Error) -> bool {
    matches!(e, crate::Error::Numerical(m) if m.contains("series tail"))
}

impl PatchKernel {
    pub fn new(sys: &EBlockSystem, threads: usize) -> Result<PatchKernel, String> {
        let system = rational_system(sys)?;
        let system_fingerprint = system.fingerprint();
        Ok(PatchKernel { system, system_fingerprint, threads: threads.max(1), pin_threads: true, max_tail_log10: None, tail_heal: true })
    }

    fn config(&self, prec: u32) -> TransportConfig {
        // the kernel's named dispersive configuration (REVIEW-A M5): real projection, Cauchy pole at z = 1
        let mut c = TransportConfig::dispersive(prec);
        c.threads = self.threads;
        c.pin_threads = self.pin_threads;
        if let Some(b) = self.max_tail_log10 {
            c.max_tail_log10 = b;
        }
        c
    }

    /// Tail heal of patch `state.next_patch` (see [`TAIL_HEAL_POLICY`]): Ok(true) = `state` advanced by the healed patch,
    /// Ok(false) = no retry passed (`state` unchanged), Err(Cancelled) = cancelled (`state` unchanged).
    fn heal_patch(
        &self,
        path: &PlannedPath,
        state: &mut crate::patch_kernel::TransportState,
        prec: u32,
        cancel: &CancellationToken,
    ) -> Result<bool, crate::Error> {
        let k = state.next_patch;
        let planned = path.patches[k].order;
        for order in tail_heal_orders(planned) {
            let mut spec = path.patches[k].clone();
            spec.order = order;
            let one = PlannedPath { offset: path.offset.clone(), patches: vec![spec] };
            let tr = crate::patch_kernel::Transport::new(&self.system, &one, self.config(prec))?.with_system_fingerprint(self.system_fingerprint.clone());
            // the same kernel step on the same state, with the patch as patch 0 of a one-patch path
            let mut st = state.clone();
            st.next_patch = 0;
            let t = std::time::Instant::now();
            match tr.step(&mut st, cancel) {
                Ok(()) => {
                    st.records.last_mut().expect("record").index = k;
                    st.next_patch = k + 1;
                    let r = st.records.last().unwrap();
                    eprintln!(
                        "[dispersive] tail heal: patch {} order {planned} -> {order}: tail {:.1} overlap {:.1} within the bound {:.1} ({:.0} s): accepted",
                        k + 1,
                        r.tail_log10,
                        r.overlap_log10,
                        tr.config.max_tail_log10,
                        t.elapsed().as_secs_f64()
                    );
                    *state = st;
                    return Ok(true);
                }
                Err(e) if is_tail_failure(&e) => {
                    eprintln!("[dispersive] tail heal: patch {} order {planned} -> {order}: {e} ({:.0} s)", k + 1, t.elapsed().as_secs_f64());
                }
                Err(e) => return Err(e),
            }
        }
        Ok(false)
    }

    /// healed patches of a kernel state (record order != planned order), for out.json
    pub fn tail_heals(st: &crate::patch_kernel::TransportState, patches: &[Patch]) -> serde_json::Value {
        serde_json::Value::Array(
            st.records
                .iter()
                .filter(|r| patches.get(r.index).is_some_and(|p| p.order != r.order))
                .map(|r| serde_json::json!({"patch": r.index, "planned_order": patches[r.index].order, "order": r.order, "tail_log10": r.tail_log10, "overlap_log10": r.overlap_log10}))
                .collect(),
        )
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
        let show = |s: &crate::patch_kernel::TransportState| {
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
        };
        let outcome = loop {
            match tr.run(&mut state, until, y0, &ctl.cancel, Some(&policy), show) {
                Ok(o) => break o,
                // the failed patch left `state` at the last good patch; heal it, checkpoint it, continue with the planned path
                Err(e) if self.tail_heal && is_tail_failure(&e) && state.next_patch < until => {
                    let k = state.next_patch;
                    eprintln!("[dispersive] patch {}/{total}: {e}; tail heal: retrying at orders {:?}", k + 1, tail_heal_orders(patches[k].order));
                    match self.heal_patch(&path, &mut state, prec, &ctl.cancel) {
                        Ok(true) => {
                            let _lock = crate::patch_kernel::checkpoint::CheckpointLock::acquire(&ctl.checkpoint).map_err(|e| e.to_string())?;
                            tr.write_checkpoint(&ctl.checkpoint, &state, y0).map_err(|e| e.to_string())?;
                            show(&state);
                        }
                        Ok(false) => return Err(format!("{e}; tail heal failed at orders {:?} (TAIL_HEAL_POLICY {TAIL_HEAL_POLICY})", tail_heal_orders(patches[k].order))),
                        // cancelled inside the heal: `state` unchanged; the next run() sees the token, flushes and returns Cancelled
                        Err(crate::Error::Cancelled) => {}
                        Err(e2) => return Err(format!("{e}; tail heal: {e2}")),
                    }
                }
                Err(e) => return Err(e.to_string()),
            }
        };
        let mut summary = Self::summary(&state);
        summary["tail_bound_log10"] = serde_json::json!(tr.config.max_tail_log10);
        summary["tail_heal"] = serde_json::json!(if self.tail_heal { TAIL_HEAL_POLICY } else { "off" });
        summary["tail_heals"] = Self::tail_heals(&state, patches);
        if state.next_patch == patches.len() {
            return Ok((Checkpointed::Completed(self.output(&state)), summary));
        }
        // stopped early (signal or stop_at_patch): make sure the checkpoint holds this state
        let _ = outcome == TransportOutcome::Partial;
        if !ctl.checkpoint.exists() {
            tr.write_checkpoint(&ctl.checkpoint, &state, y0).map_err(|e| e.to_string())?;
        }
        Ok((Checkpointed::Cancelled { next_patch: state.next_patch }, summary))
    }
}

/// the only place that builds the kernel's checkpoint policy (never `fresh`: a run resumes or refuses, REVIEW-A B1)
fn checkpoint_policy(file: &Path, every_patches: usize, every_seconds: f64) -> CheckpointPolicy {
    CheckpointPolicy { every_patches, every_seconds, fresh: false, ..CheckpointPolicy::new(file) }
}

impl Transport for PatchKernel {
    fn id(&self) -> String {
        format!(
            "{} algorithm {} (patch_kernel, dispersive config: real projection, Cauchy pole 1)",
            crate::patch_kernel::KERNEL_VERSION,
            crate::patch_kernel::KERNEL_ALGORITHM
        )
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
