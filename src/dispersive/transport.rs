//! Transport of the discontinuity along the planned patches: the INTERFACE to the fixed-precision kernel (W2) and a
//! stand-in implementation (`TaylorRecurrence`).
//!
//! Semantics (oracle disp3y-pt rpatch): patch p has centre z_c, step h, Taylor order N. y(z_c + t) = sum_n y_n t^n solves
//! y' = A(z) y. Output per patch: y(z_c + h) and -int_0^h y(z_c + t) / (1 - z_c - t) dt (on z = -s this is the piece of
//! int d(s)/(s + 1) ds). Real patches (centre and step real) take the REAL part of the incoming vector (the
//! discontinuity is real; the dropped imaginary part is recorded as a diagnostic); complex patches run in Complex.
use super::path::Patch;
use super::poly::{self, ZPoly};
use super::system::EBlockSystem;
use rug::{Assign, Complex, Float, Rational};
use std::collections::BTreeMap;
use std::time::Instant;

#[derive(Clone, Debug)]
pub struct PatchDiag {
    pub seconds: f64,
    pub setup_seconds: f64,
    /// log10 max_j |Im y_j| / |Re y_j| of the incoming vector of a real patch (dropped imaginary part)
    pub dropped_im_over_re_log10: f64,
}

pub struct TransportOutput {
    pub end: Vec<Complex>,
    pub weighted_integral: Vec<Complex>,
    pub patches: Vec<PatchDiag>,
}

/// the kernel interface (W2 implements it with checkpoints and row parallelism)
pub trait Transport {
    fn transport(&self, sys: &EBlockSystem, y0: &[Complex], patches: &[Patch], prec: u32) -> Result<TransportOutput, String>;
    /// kernel identity (algorithm + version), part of the checkpoint fingerprint
    fn id(&self) -> String;
    /// the production kernel (checkpointed per patch); the stage driver then runs the whole path through it
    fn as_patch_kernel(&self) -> Option<&super::kernel::PatchKernel> {
        None
    }
}

// ---------------- scalar for the stand-in (Float on the real axis, Complex off it) ----------------
pub trait Ts: Clone + Send + Sync {
    fn zero(prec: u32) -> Self;
    fn from_q(re: &Rational, im: &Rational, prec: u32) -> Self;
    fn add_prod(&mut self, x: &Self, y: &Self, tmp: &mut Self);
    fn sub_prod(&mut self, x: &Self, y: &Self, tmp: &mut Self);
    fn addr(&mut self, o: &Self);
    fn subr(&mut self, o: &Self);
    fn mulr(&mut self, o: &Self);
    fn divr(&mut self, o: &Self);
    fn divu(&mut self, u: u32);
    fn set_zero(&mut self);
    fn to_c(&self) -> Complex;
}
impl Ts for Float {
    fn zero(prec: u32) -> Self {
        Float::new(prec)
    }
    fn from_q(re: &Rational, _im: &Rational, prec: u32) -> Self {
        Float::with_val(prec, re)
    }
    #[inline]
    fn add_prod(&mut self, x: &Self, y: &Self, tmp: &mut Self) {
        tmp.assign(x * y);
        *self += &*tmp;
    }
    #[inline]
    fn sub_prod(&mut self, x: &Self, y: &Self, tmp: &mut Self) {
        tmp.assign(x * y);
        *self -= &*tmp;
    }
    fn addr(&mut self, o: &Self) {
        *self += o;
    }
    fn subr(&mut self, o: &Self) {
        *self -= o;
    }
    fn mulr(&mut self, o: &Self) {
        *self *= o;
    }
    fn divr(&mut self, o: &Self) {
        *self /= o;
    }
    fn divu(&mut self, u: u32) {
        *self /= u;
    }
    fn set_zero(&mut self) {
        self.assign(rug::float::Special::Zero);
    }
    fn to_c(&self) -> Complex {
        Complex::with_val(self.prec(), (self, 0))
    }
}
impl Ts for Complex {
    fn zero(prec: u32) -> Self {
        Complex::new(prec)
    }
    fn from_q(re: &Rational, im: &Rational, prec: u32) -> Self {
        Complex::with_val(prec, (re, im))
    }
    #[inline]
    fn add_prod(&mut self, x: &Self, y: &Self, tmp: &mut Self) {
        tmp.assign(x * y);
        *self += &*tmp;
    }
    #[inline]
    fn sub_prod(&mut self, x: &Self, y: &Self, tmp: &mut Self) {
        tmp.assign(x * y);
        *self -= &*tmp;
    }
    fn addr(&mut self, o: &Self) {
        *self += o;
    }
    fn subr(&mut self, o: &Self) {
        *self -= o;
    }
    fn mulr(&mut self, o: &Self) {
        *self *= o;
    }
    fn divr(&mut self, o: &Self) {
        *self /= o;
    }
    fn divu(&mut self, u: u32) {
        *self /= u;
    }
    fn set_zero(&mut self) {
        self.mut_real().assign(rug::float::Special::Zero);
        self.mut_imag().assign(rug::float::Special::Zero);
    }
    fn to_c(&self) -> Complex {
        self.clone()
    }
}

/// raw row pointer shared by the patch workers (disjoint writes, barrier-ordered reads)
struct SPtr<T>(*mut T);
unsafe impl<T: Send> Send for SPtr<T> {}
unsafe impl<T: Send + Sync> Sync for SPtr<T> {}

/// [f(0), .., f(n - 1)] by up to nt scoped threads (result in index order)
pub fn par_map<R: Send, F: Fn(usize) -> R + Sync>(n: usize, nt: usize, f: F) -> Vec<R> {
    if nt <= 1 || n < 2 {
        return (0..n).map(&f).collect();
    }
    let blocks = (nt * 8).min(n);
    let chunk = n.div_ceil(blocks);
    let ctr = std::sync::atomic::AtomicUsize::new(0);
    let (fr, cr) = (&f, &ctr);
    let mut parts: Vec<(usize, Vec<R>)> = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..nt.min(blocks))
            .map(|_| {
                sc.spawn(move || {
                    let mut out = Vec::new();
                    loop {
                        let lo = cr.fetch_add(1, std::sync::atomic::Ordering::Relaxed) * chunk;
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
        hs.into_iter().flat_map(|h| h.join().unwrap()).collect()
    });
    parts.sort_by_key(|p| p.0);
    parts.into_iter().flat_map(|p| p.1).collect()
}

/// stand-in kernel: exact dyadic Taylor shifts of every distinct polynomial, then the row-grouped polynomial
/// recurrence of disp3y-pt (groups = (row, denominator): D u = sum_t N_t y_j).
pub struct TaylorRecurrence {
    pub nthreads: usize,
    polys: Vec<ZPoly>,
    tn: Vec<usize>,
    td: Vec<usize>,
    groups: Vec<(usize, usize, Vec<usize>)>, // (row, den poly id, terms)
}

impl TaylorRecurrence {
    pub fn new(sys: &EBlockSystem, nthreads: usize) -> TaylorRecurrence {
        let mut idx: BTreeMap<ZPoly, usize> = BTreeMap::new();
        let mut polys = Vec::new();
        let mut get = |p: &ZPoly| -> usize {
            if let Some(i) = idx.get(p) {
                return *i;
            }
            polys.push(p.clone());
            idx.insert(p.clone(), polys.len() - 1);
            polys.len() - 1
        };
        let tn: Vec<usize> = sys.terms.iter().map(|t| get(&t.num)).collect();
        let td: Vec<usize> = sys.terms.iter().map(|t| get(&t.den)).collect();
        let mut gidx: BTreeMap<(usize, usize), usize> = BTreeMap::new();
        let mut groups: Vec<(usize, usize, Vec<usize>)> = Vec::new();
        for (t, term) in sys.terms.iter().enumerate() {
            let g = *gidx.entry((term.i, td[t])).or_insert_with(|| {
                groups.push((term.i, td[t], Vec::new()));
                groups.len() - 1
            });
            groups[g].2.push(t);
        }
        groups.sort_by_key(|g| g.0); // stable: rows in order, groups of a row in first-appearance order
        TaylorRecurrence { nthreads, polys, tn, td, groups }
    }

    fn patch<T: Ts>(&self, sys: &EBlockSystem, y0: &[T], p: &Patch, prec: u32) -> Result<(Vec<T>, Vec<T>, f64), String> {
        let ts = Instant::now();
        let e = sys.dim();
        let nord = p.order;
        let one = Rational::from(1);
        let ecr = Rational::from(&p.centre.0 - &one);
        let ci = &p.centre.1;
        let is_den: Vec<bool> = {
            let mut v = vec![false; self.polys.len()];
            for &d in &self.td {
                v[d] = true;
            }
            v
        };
        let shifted: Vec<(bool, Vec<T>)> = par_map(self.polys.len(), self.nthreads, |q| {
            let (re, im) = poly::shift_complex(&self.polys[q], &ecr, ci);
            // oracle: assert(patch centre on a pole): an exactly vanishing denominator at the centre
            let pole = is_den[q] && re[0] == 0 && im[0] == 0;
            (pole, re.iter().zip(&im).map(|(a, b)| T::from_q(a, b, prec)).collect())
        });
        if shifted.iter().any(|x| x.0) {
            return Err(format!("patch centre ({}, {}) is a pole of the system", p.centre.0, p.centre.1));
        }
        let shifted: Vec<Vec<T>> = shifted.into_iter().map(|x| x.1).collect();
        // normalised numerators N/D[0] per term and denominators D[1..]/D[0] per group
        let nv: Vec<Vec<T>> = par_map(sys.terms.len(), self.nthreads, |t| {
            let d0 = &shifted[self.td[t]][0];
            shifted[self.tn[t]].iter().map(|x| {
                let mut v = x.clone();
                v.divr(d0);
                v
            }).collect()
        });
        let dn: Vec<Vec<T>> = self
            .groups
            .iter()
            .map(|g| {
                let d = &shifted[g.1];
                d[1..].iter().map(|x| {
                    let mut v = x.clone();
                    v.divr(&d[0]);
                    v
                }).collect()
            })
            .collect();
        let setup = ts.elapsed().as_secs_f64();
        // rows split into contiguous chunks of groups (groups sorted by row)
        let nt = self.nthreads.max(1);
        let mut bounds = vec![0usize];
        {
            let cost: Vec<usize> = self.groups.iter().zip(&dn).map(|(g, d)| g.2.iter().map(|&t| nv[t].len()).sum::<usize>() + d.len() + 1).collect();
            let total: usize = cost.iter().sum();
            let mut acc = 0;
            for (k, g) in self.groups.iter().enumerate() {
                acc += cost[k];
                let last_of_row = k + 1 == self.groups.len() || self.groups[k + 1].0 != g.0;
                if last_of_row && acc * nt >= total * bounds.len() && bounds.len() < nt {
                    bounds.push(k + 1);
                }
            }
            if *bounds.last().unwrap() != self.groups.len() {
                bounds.push(self.groups.len());
            }
        }
        // PERSISTENT workers for the whole patch (one barrier per order): scoped threads spawned per order are too short-lived for the
        // scheduler to spread them over far-apart allotted CPUs (measured: 4 workers packed on one core). Every row is owned by one
        // worker (rows never cross a bound) and receives its groups' contributions in group order -> identical for any thread count.
        let mut rings: Vec<Vec<T>> = self.groups.iter().zip(&dn).map(|(_, d)| vec![T::zero(prec); d.len().max(1)]).collect();
        let mut y: Vec<Vec<T>> = Vec::with_capacity(nord);
        y.push(y0.to_vec());
        for _ in 1..nord {
            y.push(vec![T::zero(prec); e]);
        }
        let rows_of: Vec<Vec<usize>> = bounds
            .windows(2)
            .map(|w| {
                let mut r: Vec<usize> = self.groups[w[0]..w[1]].iter().map(|g| g.0).collect();
                r.dedup();
                r
            })
            .collect();
        {
            let yp: Vec<SPtr<T>> = y.iter_mut().map(|v| SPtr(v.as_mut_ptr())).collect();
            let nchunks = bounds.len() - 1;
            let barrier = std::sync::Barrier::new(nchunks);
            let mut ring_chunks: Vec<&mut [Vec<T>]> = Vec::new();
            let mut rest: &mut [Vec<T>] = &mut rings;
            for w in bounds.windows(2) {
                let (a, b) = rest.split_at_mut(w[1] - w[0]);
                ring_chunks.push(a);
                rest = b;
            }
            let (ypr, bar, groups, nvr, dnr, terms) = (&yp, &barrier, &self.groups, &nv, &dn, &sys.terms);
            std::thread::scope(|sc| {
                for ((w, rg), rows) in bounds.windows(2).zip(ring_chunks).zip(&rows_of) {
                    let (lo, hi) = (w[0], w[1]);
                    sc.spawn(move || {
                        let mut a = T::zero(prec);
                        let mut b = T::zero(prec);
                        let mut tmp = T::zero(prec);
                        for n in 0..nord - 1 {
                            // SAFETY: during order n rows m <= n are complete (written before the previous barrier) and only read;
                            // entry (n + 1, i) is written only by the owner of row i; the barrier orders the orders
                            let yv = |m: usize, j: usize| -> &T { unsafe { &*ypr[m].0.add(j) } };
                            for g in lo..hi {
                                let (row, _, ts) = &groups[g];
                                a.set_zero();
                                for &t in ts {
                                    let j = terms[t].j;
                                    let nt_ = &nvr[t];
                                    for k in 0..nt_.len().min(n + 1) {
                                        a.add_prod(&nt_[k], yv(n - k, j), &mut tmp);
                                    }
                                }
                                let dt = &dnr[g];
                                let dd = dt.len();
                                let ring = &mut rg[g - lo];
                                if dd > 0 && n > 0 {
                                    b.set_zero();
                                    for k in 1..=dd.min(n) {
                                        b.add_prod(&dt[k - 1], &ring[(n - k) % dd], &mut tmp);
                                    }
                                    a.subr(&b);
                                }
                                unsafe { (*ypr[n + 1].0.add(*row)).addr(&a) };
                                if dd > 0 {
                                    ring[n % dd].clone_from(&a);
                                }
                            }
                            for &i in rows {
                                unsafe { (*ypr[n + 1].0.add(i)).divu((n + 1) as u32) };
                            }
                            bar.wait();
                        }
                    });
                }
            });
        }
        let hh = T::from_q(&p.step.0, &p.step.1, prec);
        let u = T::from_q(&Rational::from(&one - &p.centre.0), &Rational::from(-ci), prec);
        let vi: Vec<(T, T)> = par_map(e, self.nthreads, |j| {
            let mut acc = T::zero(prec);
            for n in (0..nord).rev() {
                acc.mulr(&hh);
                acc.addr(&y[n][j]);
            }
            let val = acc;
            let mut pv: Vec<T> = Vec::with_capacity(nord);
            let mut cur = T::zero(prec);
            for n in 0..nord {
                cur.addr(&y[n][j]);
                cur.divr(&u);
                pv.push(cur.clone());
            }
            let mut acc = T::zero(prec);
            for n in (0..nord).rev() {
                let mut q = pv[n].clone();
                q.divu((n + 1) as u32);
                acc.addr(&q);
                acc.mulr(&hh);
            }
            let mut neg = T::zero(prec);
            neg.subr(&acc);
            (val, neg)
        });
        let (val, integ): (Vec<T>, Vec<T>) = vi.into_iter().unzip();
        Ok((val, integ, setup))
    }
}

impl Transport for TaylorRecurrence {
    fn id(&self) -> String {
        "w3a-TaylorRecurrence-v2 (exact dyadic shifts, row-grouped recurrence, persistent workers)".into()
    }
    fn transport(&self, sys: &EBlockSystem, y0: &[Complex], patches: &[Patch], prec: u32) -> Result<TransportOutput, String> {
        let e = sys.dim();
        let mut ycur: Vec<Complex> = y0.to_vec();
        let mut mid = vec![Complex::new(prec); e];
        let mut diags = Vec::with_capacity(patches.len());
        for p in patches {
            let tp = Instant::now();
            let mut dropped = -1e9f64;
            let (setup, end, ig): (f64, Vec<Complex>, Vec<Complex>) = if p.is_real() {
                let y0r: Vec<Float> = ycur
                    .iter()
                    .map(|z| {
                        let re = Float::with_val(prec, z.real());
                        if !re.is_zero() && !z.imag().is_zero() {
                            let r = z.imag().clone().abs().log10().to_f64() - re.clone().abs().log10().to_f64();
                            dropped = dropped.max(r);
                        }
                        re
                    })
                    .collect();
                let (v, ig, s) = self.patch::<Float>(sys, &y0r, p, prec)?;
                (s, v.iter().map(|x| x.to_c()).collect(), ig.iter().map(|x| x.to_c()).collect())
            } else {
                let (v, ig, s) = self.patch::<Complex>(sys, &ycur, p, prec)?;
                (s, v, ig)
            };
            ycur = end;
            for (m, x) in mid.iter_mut().zip(&ig) {
                *m += x;
            }
            diags.push(PatchDiag { seconds: tp.elapsed().as_secs_f64(), setup_seconds: setup, dropped_im_over_re_log10: dropped });
        }
        Ok(TransportOutput { end: ycur, weighted_integral: mid, patches: diags })
    }
}
