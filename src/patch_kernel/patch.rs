//! One Taylor patch: exact shift, row-parallel recurrence, finish.
//!
//! `y' = sum_t (N_t / D_t) y_j(t)`; per group g = (row i, denominator D),
//! `u_g = sum_{t in g} (N_t / D) y_j(t)` obeys `D u_g = sum_t N_t y_j(t)`, so
//! `u_g[n] = sum_t sum_k nv_t[k] y_j[n - k] - sum_{k >= 1} dv[k] u_g[n - k]`
//! with `nv = N / D(0)`, `dv = D / D(0)` (coefficients at the patch centre), and
//! `y_i[n + 1] = (sum_{g of row i} u_g[n]) / (n + 1)`.
//!
//! Determinism: per order, phase A computes every group value `u_g[n]`
//! (fixed operation sequence per group) into the group's own ring, phase B
//! gives row `i` of order `n + 1`: zero, one addition of `u_g[n]` per group in
//! canonical order, one division by `n + 1`. Work is taken dynamically in
//! fixed chunks (fast threads take more, so a contended CPU does not stall the
//! others), but which thread computes a value never changes its operations.
//! Phases are separated by barriers, so every read sees completed values.
use super::scalar::KernelScalar;
use super::shift::{Centre, shift_guarded};
use super::system::RationalSystem;
use crate::options::CancellationToken;
use crate::{Error, Result};
use rug::{Complex, Float, Rational};
use std::sync::Barrier;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::time::Instant;

/// Guard bits of the exact-shift conversion (each normalised coefficient is
/// the correctly rounded quotient of two values accurate to `bits + GUARD`).
pub const SHIFT_GUARD_BITS: u32 = 64;

/// Patch-independent structure of the recurrence for one system.
#[derive(Clone, Debug)]
pub struct KernelPlan {
    /// (row, den poly, entry indices) in canonical order.
    pub groups: Vec<(usize, usize, Vec<usize>)>,
    /// Distinct denominators in first-use order.
    pub dens: Vec<usize>,
    /// den poly -> index into `dens`.
    den_slot: std::collections::BTreeMap<usize, usize>,
    /// Cost of each row per order at full depth (for the worker assignment).
    pub row_cost: Vec<usize>,
    /// Groups of row i: `row_groups[i].0 .. row_groups[i].1` (groups are row-sorted).
    pub row_groups: Vec<(usize, usize)>,
    /// Work chunks of phase A (group ranges of similar cost) and phase B (row ranges).
    pub group_chunks: Vec<(usize, usize)>,
    pub row_chunks: Vec<(usize, usize)>,
}

impl KernelPlan {
    pub fn new(system: &RationalSystem) -> Self {
        let mut gidx = std::collections::BTreeMap::new();
        let mut groups: Vec<(usize, usize, Vec<usize>)> = Vec::new();
        let mut dens = Vec::new();
        let mut den_slot = std::collections::BTreeMap::new();
        for (t, e) in system.entries.iter().enumerate() {
            let g = *gidx.entry((e.row, e.den)).or_insert_with(|| {
                groups.push((e.row, e.den, Vec::new()));
                groups.len() - 1
            });
            groups[g].2.push(t);
            den_slot.entry(e.den).or_insert_with(|| {
                dens.push(e.den);
                dens.len() - 1
            });
        }
        // Entries are sorted by row, so groups are already in row order;
        // keep a stable sort for safety.
        groups.sort_by_key(|g| g.0);
        let mut row_cost = vec![0usize; system.dim];
        for (row, den, ts) in &groups {
            row_cost[*row] += ts
                .iter()
                .map(|&t| system.polys[system.entries[t].num].len())
                .sum::<usize>()
                + system.polys[*den].len();
        }
        let mut row_groups = vec![(0usize, 0usize); system.dim];
        let mut g = 0;
        for (i, rg) in row_groups.iter_mut().enumerate() {
            let start = g;
            while g < groups.len() && groups[g].0 == i {
                g += 1;
            }
            *rg = (start, g);
        }
        let gcost: Vec<usize> = groups
            .iter()
            .map(|(_, den, ts)| {
                ts.iter()
                    .map(|&t| system.polys[system.entries[t].num].len())
                    .sum::<usize>()
                    + system.polys[*den].len()
            })
            .collect();
        let total: usize = gcost.iter().sum();
        let target = (total / 2048).max(1);
        let mut group_chunks = Vec::new();
        let (mut start, mut acc) = (0usize, 0usize);
        for (k, c) in gcost.iter().enumerate() {
            acc += c;
            if acc >= target {
                group_chunks.push((start, k + 1));
                start = k + 1;
                acc = 0;
            }
        }
        if start < groups.len() {
            group_chunks.push((start, groups.len()));
        }
        let row_chunks = (0..system.dim)
            .step_by(8)
            .map(|i| (i, (i + 8).min(system.dim)))
            .collect();
        Self {
            groups,
            dens,
            den_slot,
            row_cost,
            row_groups,
            group_chunks,
            row_chunks,
        }
    }
}

/// CPUs this process may run on (Linux; empty elsewhere), increasing.
pub fn allowed_cpus() -> Vec<usize> {
    #[cfg(target_os = "linux")]
    {
        // SAFETY: plain libc calls on a zero-initialised cpu_set_t.
        unsafe {
            let mut set: libc::cpu_set_t = std::mem::zeroed();
            if libc::sched_getaffinity(0, std::mem::size_of::<libc::cpu_set_t>(), &mut set) != 0 {
                return Vec::new();
            }
            (0..libc::CPU_SETSIZE as usize)
                .filter(|&c| libc::CPU_ISSET(c, &set))
                .collect()
        }
    }
    #[cfg(not(target_os = "linux"))]
    {
        Vec::new()
    }
}

/// Pin the calling thread to `cpu` (best effort; results never depend on it).
/// The kernel's NUMA-imbalance allowance otherwise keeps a few threads of a
/// job whose CPUs span NUMA nodes on ONE CPU.
pub fn pin_current_thread(cpu: usize) {
    #[cfg(target_os = "linux")]
    // SAFETY: plain libc calls; 0 = the calling thread.
    unsafe {
        let mut set: libc::cpu_set_t = std::mem::zeroed();
        libc::CPU_SET(cpu, &mut set);
        let _ = libc::sched_setaffinity(0, std::mem::size_of::<libc::cpu_set_t>(), &set);
    }
    #[cfg(not(target_os = "linux"))]
    let _ = cpu;
}

/// `[f(0), .., f(n - 1)]` computed by up to `threads` scoped threads (thread k
/// pinned to `pins[k]` if given); the result is in index order and `f(i)`
/// itself does not depend on the threads.
pub fn par_map<R: Send, F: Fn(usize) -> R + Sync>(
    n: usize,
    threads: usize,
    pins: &[usize],
    f: F,
) -> Vec<R> {
    if threads <= 1 || n < 2 {
        return (0..n).map(&f).collect();
    }
    let blocks = (threads * 16).min(n);
    let chunk = n.div_ceil(blocks);
    let counter = AtomicUsize::new(0);
    let (fr, cr) = (&f, &counter);
    let mut parts: Vec<(usize, Vec<R>)> = std::thread::scope(|sc| {
        let handles: Vec<_> = (0..threads.min(blocks))
            .map(|k| {
                let pin = pins.get(k).copied();
                sc.spawn(move || {
                    if let Some(c) = pin {
                        pin_current_thread(c);
                    }
                    let mut out = Vec::new();
                    loop {
                        let lo = cr.fetch_add(1, Ordering::Relaxed) * chunk;
                        if lo >= n {
                            break;
                        }
                        let hi = (lo + chunk).min(n);
                        out.push((lo, (lo..hi).map(fr).collect::<Vec<R>>()));
                    }
                    out
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
            .collect()
    });
    parts.sort_by_key(|p| p.0);
    parts.into_iter().flat_map(|p| p.1).collect()
}

/// Inputs of one patch. Points are in the system variable `x`.
pub struct PatchJob<'a> {
    pub system: &'a RationalSystem,
    pub plan: &'a KernelPlan,
    pub centre: (&'a Rational, &'a Rational),
    pub step: (&'a Rational, &'a Rational),
    /// Number of Taylor coefficients (orders 0 .. order - 1).
    pub order: usize,
    pub bits: u32,
    pub threads: usize,
    /// Integral weight `1 / (u - t)` with `u = a - centre` (if any).
    pub cauchy_u: Option<(Rational, Rational)>,
    /// Extra evaluation offsets from the centre (overlap probes).
    pub probes: Vec<(Rational, Rational)>,
    pub cancel: &'a CancellationToken,
    /// Pin worker threads to distinct allowed CPUs.
    pub pin_threads: bool,
}

/// Result of one patch (state values exact as `rug::Complex`).
pub struct PatchOutput {
    pub value: Vec<Complex>,
    pub integral: Option<Vec<Complex>>,
    pub probes: Vec<Vec<Complex>>,
    /// Per component log10(|y_{N-1} h^{N-1}| / max(|y(0)|, |y(h)|)).
    pub tail_log10: Vec<f64>,
    pub seconds_setup: f64,
    pub seconds_recurrence: f64,
    pub seconds_finish: f64,
    /// Seconds each recurrence worker spent computing (outside barriers);
    /// informational, not part of any identity.
    pub worker_busy: Vec<f64>,
}

struct SPtr<T>(*mut T);
// SAFETY: see `recurrence`: disjoint writes, barrier-ordered reads.
unsafe impl<T: Send> Send for SPtr<T> {}
unsafe impl<T: Send + Sync> Sync for SPtr<T> {}

/// log10 |x| (finite, -1e9 for zero).
pub fn log10_abs(x: &Complex) -> f64 {
    let a = Float::with_val(64, x.abs_ref());
    if a.is_zero() {
        -1e9
    } else {
        a.log10().to_f64()
    }
}

/// Run one patch on scalar type `S` from the exact state `y0`.
pub fn run_patch<S: KernelScalar>(job: &PatchJob, y0: &[Complex]) -> Result<PatchOutput> {
    let system = job.system;
    let plan = job.plan;
    let bits = job.bits;
    let dim = system.dim;
    let nord = job.order.max(2);
    let threads = job.threads.max(1);
    if y0.len() != dim {
        return Err(Error::InvalidInput("state dimension".into()));
    }
    // ---- setup: exact shift of every polynomial to the centre ----
    let t_setup = Instant::now();
    let centre = Centre::new(job.centre.0, job.centre.1);
    if !S::COMPLEX && !centre.is_real() {
        return Err(Error::InvalidInput(
            "real arithmetic at a complex centre".into(),
        ));
    }
    let gbits = bits + SHIFT_GUARD_BITS;
    let cpus = if job.pin_threads && threads > 1 {
        allowed_cpus()
    } else {
        Vec::new()
    };
    let pins: &[usize] = if cpus.len() >= threads {
        &cpus[..threads]
    } else {
        &[]
    };
    // denominators: D(0) (guarded) and dv[k] = D[k + 1] / D[0]
    let den_data: Vec<(Complex, Vec<S>)> = par_map(plan.dens.len(), threads, pins, |q| {
        let s = shift_guarded(&system.polys[plan.dens[q]], &centre, gbits);
        let dv = s[1..].iter().map(|c| S::quotient(c, &s[0], bits)).collect();
        (s[0].clone(), dv)
    });
    for (d0, _) in &den_data {
        if d0.real().is_zero() && d0.imag().is_zero() {
            return Err(Error::Numerical(
                "patch centre is a pole of the system".into(),
            ));
        }
    }
    job.cancel.check()?;
    let nv: Vec<Vec<S>> = par_map(system.entries.len(), threads, pins, |t| {
        let e = &system.entries[t];
        let s = shift_guarded(&system.polys[e.num], &centre, gbits);
        let d0 = &den_data[plan.den_slot[&e.den]].0;
        if e.scale == 1 {
            s.iter().map(|c| S::quotient(c, d0, bits)).collect()
        } else {
            // scale * N[k] / D(0): the guarded scale enters the guarded numerator
            let f = rug::Float::with_val(gbits, &e.scale);
            s.iter()
                .map(|c| S::quotient(&Complex::with_val(gbits, c * &f), d0, bits))
                .collect()
        }
    });
    let seconds_setup = t_setup.elapsed().as_secs_f64();
    job.cancel.check()?;
    // ---- recurrence ----
    // Two phases per order, both dynamically scheduled over fixed chunks:
    // (A) every group g computes u_g[n] into its own ring slot; (B) every row i
    // sums u_g[n] over its groups in canonical order and divides by n + 1.
    // Each value receives the same operations in the same order whichever
    // thread computes it, so the result does not depend on the threads.
    let t_rec = Instant::now();
    let workers = threads.min(plan.group_chunks.len()).max(1);
    // Component-major storage y[j][m]: the inner loops walk the orders of one
    // component, whose values (and limbs) are allocated contiguously.
    let mut y: Vec<Vec<S>> = y0
        .iter()
        .map(|v| {
            let mut col = Vec::with_capacity(nord);
            col.push(S::from_state(v, bits));
            for _ in 1..nord {
                col.push(S::zero(bits));
            }
            col
        })
        .collect();
    // ring of the D-recurrence per group; after phase A of order n, slot
    // n % dd (slot 0 if dd = 0) holds u_g[n]
    let mut rings: Vec<Vec<S>> = plan
        .groups
        .iter()
        .map(|(_, den, _)| vec![S::zero(bits); den_data[plan.den_slot[den]].1.len().max(1)])
        .collect();
    let (cancelled, worker_busy) = {
        let ptrs: Vec<SPtr<S>> = y.iter_mut().map(|v| SPtr(v.as_mut_ptr())).collect();
        let rptrs: Vec<SPtr<S>> = rings.iter_mut().map(|v| SPtr(v.as_mut_ptr())).collect();
        let counters_a: Vec<AtomicUsize> = (0..nord).map(|_| AtomicUsize::new(0)).collect();
        let counters_b: Vec<AtomicUsize> = (0..nord).map(|_| AtomicUsize::new(0)).collect();
        let barrier = Barrier::new(workers);
        let stop = AtomicBool::new(false);
        let (ptrs, rptrs, ca, cb, barrier, stop, nv, den_data) = (
            &ptrs, &rptrs, &counters_a, &counters_b, &barrier, &stop, &nv, &den_data,
        );
        std::thread::scope(|sc| {
            let handles: Vec<_> = (0..workers)
                .map(|w| {
                    let pin = pins.get(w).copied();
                    sc.spawn(move || -> (bool, f64) {
                        if let Some(c) = pin {
                            pin_current_thread(c);
                        }
                        let mut busy = 0.0f64;
                        let mut t_busy = Instant::now();
                        let mut a = S::zero(bits);
                        let mut b = S::zero(bits);
                        let mut tmp = S::zero(bits);
                        for n in 0..nord - 1 {
                            // SAFETY (phase A): y orders m <= n are complete and not
                            // written in this phase; ring g is touched only by the
                            // thread that took the chunk containing g.
                            let yv = |m: usize, j: usize| -> &S { unsafe { &*ptrs[j].0.add(m) } };
                            loop {
                                let c = ca[n].fetch_add(1, Ordering::Relaxed);
                                let Some(&(g0, g1)) = plan.group_chunks.get(c) else {
                                    break;
                                };
                                for g in g0..g1 {
                                    let (_, den, ts) = &plan.groups[g];
                                    let dv = &den_data[plan.den_slot[den]].1[..];
                                    let dd = dv.len();
                                    let ring = unsafe {
                                        std::slice::from_raw_parts_mut(rptrs[g].0, dd.max(1))
                                    };
                                    a.set_zero();
                                    for &t in ts.iter() {
                                        let c = &nv[t];
                                        let j = system.entries[t].col;
                                        for k in 0..c.len().min(n + 1) {
                                            a.add_prod(&c[k], yv(n - k, j), &mut tmp);
                                        }
                                    }
                                    if dd > 0 && n > 0 {
                                        b.set_zero();
                                        for k in 1..=dd.min(n) {
                                            b.add_prod(&dv[k - 1], &ring[(n - k) % dd], &mut tmp);
                                        }
                                        a.sub_ref(&b);
                                    }
                                    let slot = if dd > 0 { n % dd } else { 0 };
                                    std::mem::swap(&mut ring[slot], &mut a);
                                }
                            }
                            busy += t_busy.elapsed().as_secs_f64();
                            barrier.wait();
                            t_busy = Instant::now();
                            // SAFETY (phase B): rings are only read; y[i][n + 1] is
                            // written only by the thread that took row i's chunk.
                            loop {
                                let c = cb[n].fetch_add(1, Ordering::Relaxed);
                                let Some(&(r0, r1)) = plan.row_chunks.get(c) else {
                                    break;
                                };
                                for i in r0..r1 {
                                    let target = unsafe { &mut *ptrs[i].0.add(n + 1) };
                                    let (g0, g1) = plan.row_groups[i];
                                    for g in g0..g1 {
                                        let dd = den_data[plan.den_slot[&plan.groups[g].1]].1.len();
                                        let slot = if dd > 0 { n % dd } else { 0 };
                                        target.add_ref(unsafe { &*rptrs[g].0.add(slot) });
                                    }
                                    target.div_u32((n + 1) as u32);
                                }
                            }
                            busy += t_busy.elapsed().as_secs_f64();
                            let leader = barrier.wait().is_leader();
                            if (n + 1) % 8 == 0 {
                                // One consistent stop decision: written by the leader
                                // between two barriers, read by all after the second.
                                if leader {
                                    stop.store(job.cancel.is_cancelled(), Ordering::Relaxed);
                                }
                                barrier.wait();
                                if stop.load(Ordering::Relaxed) {
                                    return (true, busy);
                                }
                            }
                            t_busy = Instant::now();
                        }
                        (false, busy)
                    })
                })
                .collect();
            let results: Vec<(bool, f64)> = handles
                .into_iter()
                .map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
                .collect();
            (
                results.iter().any(|r| r.0),
                results.iter().map(|r| r.1).collect::<Vec<f64>>(),
            )
        })
    };
    if cancelled {
        return Err(Error::Cancelled);
    }
    let seconds_recurrence = t_rec.elapsed().as_secs_f64();
    // ---- finish: per component, independent ----
    let t_fin = Instant::now();
    let h = S::from_rational(job.step.0, job.step.1, bits);
    let u = job
        .cauchy_u
        .as_ref()
        .map(|(r, i)| S::from_rational(r, i, bits));
    let probes: Vec<S> = job
        .probes
        .iter()
        .map(|(r, i)| S::from_rational(r, i, bits))
        .collect();
    let habs = Float::with_val(
        64,
        Complex::with_val(64, (job.step.0, job.step.1)).abs_ref(),
    );
    let horner = |j: usize, x: &S| -> S {
        let mut acc = S::zero(bits);
        for n in (0..nord).rev() {
            acc.mul_ref(x);
            acc.add_ref(&y[j][n]);
        }
        acc
    };
    let per: Vec<(Complex, Option<Complex>, Vec<Complex>, f64)> =
        par_map(dim, threads, pins, |j| {
            let val = horner(j, &h);
            let integral = u.as_ref().map(|u| {
                let mut p: Vec<S> = Vec::with_capacity(nord);
                let mut cur = S::zero(bits);
                for n in 0..nord {
                    cur.add_ref(&y[j][n]);
                    cur.div_ref(u);
                    p.push(cur.clone());
                }
                let mut acc = S::zero(bits);
                for n in (0..nord).rev() {
                    let mut q = p[n].clone();
                    q.div_u32((n + 1) as u32);
                    acc.add_ref(&q);
                    acc.mul_ref(&h);
                }
                acc.to_state()
            });
            let pv = probes.iter().map(|x| horner(j, x).to_state()).collect();
            let val_c = val.to_state();
            let scale = log10_abs(&val_c).max(log10_abs(&y0[j]));
            let last = log10_abs(&y[j][nord - 1].to_state());
            let tail = if last < -1e8 || scale < -1e8 {
                -1e9
            } else {
                last + (nord - 1) as f64 * habs.to_f64().log10() - scale
            };
            (val_c, integral, pv, tail)
        });
    let seconds_finish = t_fin.elapsed().as_secs_f64();
    let mut value = Vec::with_capacity(dim);
    let mut integral = job.cauchy_u.as_ref().map(|_| Vec::with_capacity(dim));
    let mut probe_vals: Vec<Vec<Complex>> = vec![Vec::with_capacity(dim); job.probes.len()];
    let mut tail_log10 = Vec::with_capacity(dim);
    for (v, i, pv, t) in per {
        value.push(v);
        if let (Some(acc), Some(i)) = (integral.as_mut(), i) {
            acc.push(i);
        }
        for (k, p) in pv.into_iter().enumerate() {
            probe_vals[k].push(p);
        }
        tail_log10.push(t);
    }
    Ok(PatchOutput {
        value,
        integral,
        probes: probe_vals,
        tail_log10,
        seconds_setup,
        seconds_recurrence,
        seconds_finish,
        worker_busy,
    })
}
