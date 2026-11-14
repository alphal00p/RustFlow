//! Linear-algebra Frobenius method for formal local solutions x^lam sum_{n >= n_lo} c_n x^n of theta c = B(x) c
//! (non-Fuchsian Laporta basis; oracle: local_solutions.local_dim / prep5.modp_family + family_spec / disp frobenius
//! + stable_basis). The truncated solution space is an explicit parametrisation T (rows = coefficients
//! (n - n_lo) * E + component, columns = free parameters), refined order by order by incremental null spaces of
//!   G = (lam + m) T_m - sum_k B_k T_{m-k}.
//! The SAME code runs exactly mod p (ranks are measured; the family structure: rank plateau = dimension, k_stable,
//! pivot coordinates) and in MPFR (ranks imposed from the mod-p run, complete pivoting, cancellation noise snapped
//! to zero as the oracle). T is stored as sparse rows.
use super::field::{ModP, Scalar};
use super::laurent::Lmat;
use rug::{Float, Rational};
use std::time::Instant;

pub type SparseRow<S> = Vec<(u32, S)>;

#[derive(Clone, Debug, Default)]
pub struct FrobDiag {
    pub pole_rank: i64,
    pub steps: usize,
    pub q_final: usize,
    pub top: i64,
    /// min over steps of log10(last pivot) - log10(largest rejected entry): the rank gap (MPFR runs)
    pub min_rank_gap_log10: f64,
    /// smallest accepted pivot (log10, MPFR runs)
    pub min_pivot_log10: f64,
    pub t_nnz_max: usize,
    pub seconds: f64,
    /// structural pivoting: smallest pivot relative to the largest entry of its row, and the sum of these (log10)
    pub min_structural_pivot_rel_log10: f64,
    pub sum_structural_pivot_rel_log10: f64,
    /// structural pivoting fell back to complete pivoting from this step on (pivot unresolvable at the working precision)
    pub structural_fallback_step: Option<usize>,
}

pub struct FrobOut<S> {
    pub t: Vec<SparseRow<S>>,
    pub q: usize,
    pub ranks: Vec<usize>,
    /// pivot (row of G, column of T) sequence of every step
    pub pivots: Vec<Vec<(u32, u32)>>,
    pub diag: FrobDiag,
}

/// pivot strategy of the MPFR run
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pivoting {
    /// largest entry (oracle default)
    Complete,
    /// oracle DISP_MARK: sparsest among the entries within 2^-tau of the largest
    Threshold(i32),
    /// the exact run's pivot sequence: T keeps the exact sparsity; small (near-resonant) pivots cost rounding digits, paid by guard bits
    Structural,
}

impl Pivoting {
    pub fn from_env() -> Pivoting {
        match std::env::var("DISPERSIVE_PIVOTING").as_deref() {
            Ok("complete") => Pivoting::Complete,
            Ok("structural") | Err(_) => Pivoting::Structural,
            Ok(x) if x.starts_with("shadow") => Pivoting::Structural, // shadow modes: handled by run::family_num
            Ok(x) if x.starts_with("threshold:") => Pivoting::Threshold(x[10..].parse().expect("threshold:<bits>")),
            Ok(x) => panic!("DISPERSIVE_PIVOTING={x}: complete | structural | threshold:<bits>"),
        }
    }
}

/// run the Frobenius recursion; `ranks = None` measures the ranks (exact types only)
#[allow(clippy::too_many_arguments)]
pub fn frobenius<S: Scalar>(b: &Lmat<S>, e: usize, lam: &S, n_lo: i64, mm: i64, ranks: Option<&[usize]>, prec: u32, nthreads: usize) -> FrobOut<S> {
    let pv = match Pivoting::from_env() {
        Pivoting::Structural => Pivoting::Complete, // no pivot sequence given
        x => x,
    };
    frobenius_with(b, e, lam, n_lo, mm, ranks, None, pv, prec, nthreads)
}

/// as `frobenius`; `pivots` = the exact run's pivot sequence (needed for Pivoting::Structural)
#[allow(clippy::too_many_arguments)]
pub fn frobenius_with<S: Scalar>(b: &Lmat<S>, e: usize, lam: &S, n_lo: i64, mm: i64, ranks: Option<&[usize]>, pivots: Option<&[Vec<(u32, u32)>]>, pivoting: Pivoting, prec: u32, nthreads: usize) -> FrobOut<S> {
    assert!(ranks.is_some() || S::EXACT, "MPFR Frobenius needs the ranks of the exact run");
    // ModP: first nonzero. Shadowed: the requested strategy on the float magnitudes among EXACT nonzeros (Structural = first nonzero).
    let mut pivoting = if S::EXACT && !S::SNAP { Pivoting::Complete } else { pivoting };
    let first_nonzero = S::EXACT && (!S::SNAP || pivoting == Pivoting::Structural);
    assert!(pivoting != Pivoting::Structural || pivots.is_some() || S::EXACT, "structural pivoting needs the exact pivot sequence");
    let mut pivots_out: Vec<Vec<(u32, u32)>> = Vec::new();
    let clock = Instant::now();
    let r = -*b.keys().next().unwrap();
    let keep = prec as i32 - 40;
    let mut t: Vec<SparseRow<S>> = Vec::new();
    let mut q = 0usize;
    let mut top = n_lo - 1;
    let mut step = 0usize;
    let mut ranks_out = Vec::new();
    let mut diag = FrobDiag { pole_rank: r, min_rank_gap_log10: f64::INFINITY, min_pivot_log10: f64::INFINITY, min_structural_pivot_rel_log10: f64::INFINITY, ..Default::default() };
    let nt = nthreads.max(1);
    for m in (n_lo - r)..(n_lo + mm) {
        let nb = m + r;
        if nb >= n_lo {
            for c in 0..e {
                t.push(vec![((q + c) as u32, S::one(prec))]);
            }
            q += e;
            top = nb;
        }
        diag.t_nnz_max = diag.t_nnz_max.max(t.iter().map(|x| x.len()).sum());
        // ---- G (E x q, dense rows) ----
        let lamm = {
            let mut l = lam.clone();
            l.add_assign_s(&S::from_rational(&Rational::from(m), prec));
            l
        };
        let mut pairs: Vec<Vec<(&S, usize)>> = vec![Vec::new(); e];
        for (k, bk) in b {
            let jb = m - k;
            if jb < n_lo || jb > top {
                continue;
            }
            let base = ((jb - n_lo) as usize) * e;
            for (i, j, v) in bk {
                pairs[*i].push((v, base + *j));
            }
        }
        let mut g: Vec<Vec<S>> = vec![Vec::new(); e];
        {
            let chunk = e.div_ceil(nt);
            let (tr, pr, lr) = (&t, &pairs, &lamm);
            std::thread::scope(|sc| {
                for (ci, gch) in g.chunks_mut(chunk).enumerate() {
                    sc.spawn(move || {
                        let mut tmp = S::zero(prec);
                        let mut emax = vec![i32::MIN; q];
                        for (ii, grow) in gch.iter_mut().enumerate() {
                            let i = ci * chunk + ii;
                            *grow = vec![S::zero(prec); q];
                            if m >= n_lo {
                                for (c, v) in &tr[((m - n_lo) as usize) * e + i] {
                                    grow[*c as usize] = v.mul_s(lr, prec);
                                }
                            }
                            if pr[i].is_empty() {
                                continue;
                            }
                            if S::SNAP {
                                for c in 0..q {
                                    emax[c] = grow[c].exponent();
                                }
                            }
                            for (v, rw) in &pr[i] {
                                for (c, x) in &tr[*rw] {
                                    let c = *c as usize;
                                    let ex = grow[c].sub_mul(v, x, &mut tmp);
                                    if ex > emax[c] {
                                        emax[c] = ex;
                                    }
                                }
                            }
                            if S::SNAP {
                                for c in 0..q {
                                    grow[c].snap(emax[c], keep);
                                }
                            }
                        }
                    });
                }
            });
        }
        // ---- elimination: complete pivoting (MPFR) / first nonzero (mod p) ----
        let target = ranks.map(|rk| rk[step]);
        step += 1;
        let mut used_r = vec![false; e];
        let mut used_c = vec![false; q];
        let mut prow: Vec<usize> = Vec::new();
        let mut pcol: Vec<usize> = Vec::new();
        let mut tmp = S::zero(prec);
        let mut last = f64::INFINITY;
        let markowitz: Option<i32> = if let Pivoting::Threshold(t) = pivoting { Some(t) } else { None };
        let tcol: Vec<u64> = if markowitz.is_some() {
            let mut v = vec![0u64; q];
            for row in &t {
                for (c, _) in row {
                    v[*c as usize] += 1;
                }
            }
            v
        } else {
            Vec::new()
        };
        loop {
            if let Some(rk) = target {
                if prow.len() == rk {
                    break;
                }
            }
            let (mut bi, mut bc) = (usize::MAX, usize::MAX);
            if pivoting == Pivoting::Structural && !S::EXACT {
                let (r, c) = pivots.unwrap()[step - 1][prow.len()];
                (bi, bc) = (r as usize, c as usize);
                // pivot size relative to the largest entry of its G row: the rounding loss this pivot can cause
                let mut rowmax = S::zero(prec);
                for c in 0..q {
                    if !used_c[c] && g[bi][c].cmp_mag(&rowmax) == std::cmp::Ordering::Greater {
                        rowmax = g[bi][c].clone();
                    }
                }
                let rel = g[bi][bc].log10_abs() - rowmax.log10_abs();
                if g[bi][bc].is_zero() || rel < -(prec as f64 * 0.30103) + 20.0 {
                    // unresolvable at this precision: the remaining pivots of this family are chosen by complete pivoting (the
                    // structural sequence no longer applies once a pivot differs); recorded in the diagnostics
                    diag.structural_fallback_step = Some(step);
                    pivoting = Pivoting::Complete;
                    bi = usize::MAX;
                    bc = usize::MAX;
                } else {
                    diag.min_structural_pivot_rel_log10 = diag.min_structural_pivot_rel_log10.min(rel);
                    diag.sum_structural_pivot_rel_log10 += rel.min(0.0);
                }
            }
            if bi != usize::MAX {
            } else if first_nonzero {
                'outer: for i in 0..e {
                    if used_r[i] {
                        continue;
                    }
                    for c in 0..q {
                        if !used_c[c] && !g[i][c].is_zero() {
                            bi = i;
                            bc = c;
                            break 'outer;
                        }
                    }
                }
                if S::SNAP && bi != usize::MAX {
                    // shadowed run: the exact pivot's size relative to its row (rounding growth diagnostic)
                    let mut rowmax = S::zero(prec);
                    for c in 0..q {
                        if !used_c[c] && g[bi][c].cmp_mag(&rowmax) == std::cmp::Ordering::Greater {
                            rowmax = g[bi][c].clone();
                        }
                    }
                    let rel = g[bi][bc].log10_abs() - rowmax.log10_abs();
                    diag.min_structural_pivot_rel_log10 = diag.min_structural_pivot_rel_log10.min(rel);
                    diag.sum_structural_pivot_rel_log10 += rel.min(0.0);
                }
            } else {
                // largest |G| among unused rows/columns; ties -> first in row-major order (thread-count independent)
                let chunk = e.div_ceil(nt);
                let (gr, ur, uc) = (&g, &used_r, &used_c);
                let scan = move |lo: usize| -> Option<(usize, usize)> {
                    let mut best: Option<(usize, usize)> = None;
                    for i in lo..(lo + chunk).min(e) {
                        if ur[i] {
                            continue;
                        }
                        for c in 0..q {
                            if !uc[c] && !gr[i][c].is_zero() && best.is_none_or(|(x, y)| gr[i][c].cmp_mag(&gr[x][y]) == std::cmp::Ordering::Greater) {
                                best = Some((i, c));
                            }
                        }
                    }
                    best
                };
                let cands: Vec<Option<(usize, usize)>> = if nt > 1 {
                    std::thread::scope(|sc| {
                        let hs: Vec<_> = (0..e).step_by(chunk).map(|lo| sc.spawn(move || scan(lo))).collect();
                        hs.into_iter().map(|h| h.join().unwrap()).collect()
                    })
                } else {
                    vec![scan(0)]
                };
                for (x, y) in cands.into_iter().flatten() {
                    if bi == usize::MAX || g[x][y].cmp_mag(&g[bi][bc]) == std::cmp::Ordering::Greater {
                        bi = x;
                        bc = y;
                    }
                }
                if let (Some(tau), true) = (markowitz, bi != usize::MAX) {
                    // oracle DISP_MARK: threshold pivoting keeps T sparse (complete pivoting fills it); the basis is invariant
                    let lim = g[bi][bc].exponent().saturating_sub(tau);
                    let rc: Vec<u64> = (0..e).map(|i| if used_r[i] { 0 } else { (0..q).filter(|&c| !used_c[c] && !g[i][c].is_zero()).count() as u64 }).collect();
                    let mut cc = vec![0u64; q];
                    for i in 0..e {
                        if !used_r[i] {
                            for c in 0..q {
                                if !used_c[c] && !g[i][c].is_zero() {
                                    cc[c] += 1;
                                }
                            }
                        }
                    }
                    let mut bcost = u64::MAX;
                    for i in 0..e {
                        if used_r[i] {
                            continue;
                        }
                        for c in 0..q {
                            if !used_c[c] && !g[i][c].is_zero() && g[i][c].exponent() >= lim {
                                let cost = (rc[i] - 1) * (tcol[c] + cc[c] - 1);
                                if cost < bcost {
                                    bcost = cost;
                                    bi = i;
                                    bc = c;
                                }
                            }
                        }
                    }
                }
            }
            if bi == usize::MAX {
                assert!(target.is_none(), "rank deficit at step {step}: {} of {:?} pivots", prow.len(), target);
                break;
            }
            last = g[bi][bc].log10_abs();
            assert!(last > -1e8, "pivot ({bi}, {bc}) of step {step} is zero in floating point: raise the precision");
            used_r[bi] = true;
            used_c[bc] = true;
            let pv = g[bi][bc].clone();
            let mut prow_v: Vec<(usize, S)> = Vec::new();
            for c in 0..q {
                if !g[bi][c].is_zero() {
                    g[bi][c] = g[bi][c].div_s(&pv, prec);
                    prow_v.push((c, g[bi][c].clone()));
                }
            }
            let elim = |lo: usize, gch: &mut [Vec<S>], tmp: &mut S| {
                for (ii, grow) in gch.iter_mut().enumerate() {
                    if lo + ii == bi || grow[bc].is_zero() {
                        continue;
                    }
                    let f = grow[bc].clone();
                    for (c, pvv) in &prow_v {
                        let e0 = grow[*c].exponent();
                        let ex = grow[*c].sub_mul(pvv, &f, tmp).max(e0);
                        grow[*c].snap(ex, keep);
                    }
                    grow[bc].set_zero();
                }
            };
            if nt > 1 && S::SNAP {
                let chunk = e.div_ceil(nt);
                let el = &elim;
                std::thread::scope(|sc| {
                    for (ci, gch) in g.chunks_mut(chunk).enumerate() {
                        sc.spawn(move || {
                            let mut t2 = S::zero(prec);
                            el(ci * chunk, gch, &mut t2)
                        });
                    }
                });
            } else {
                elim(0, &mut g, &mut tmp);
            }
            prow.push(bi);
            pcol.push(bc);
        }
        ranks_out.push(prow.len());
        pivots_out.push(prow.iter().zip(&pcol).map(|(&r, &c)| (r as u32, c as u32)).collect());
        if !S::EXACT && !prow.is_empty() {
            let mut left = -1e9f64;
            for i in 0..e {
                if used_r[i] {
                    continue;
                }
                for c in 0..q {
                    if !used_c[c] {
                        left = left.max(g[i][c].log10_abs());
                    }
                }
            }
            diag.min_rank_gap_log10 = diag.min_rank_gap_log10.min(last - left);
            diag.min_pivot_log10 = diag.min_pivot_log10.min(last);
        }
        if prow.is_empty() {
            continue;
        }
        // ---- T <- T N: new column a = old free column f_a minus sum_p (old pivot column c_p) G[r_p][f_a] ----
        let mut newidx = vec![u32::MAX; q];
        let mut pidx = vec![usize::MAX; q];
        let mut nf = 0u32;
        for c in 0..q {
            if !used_c[c] {
                newidx[c] = nf;
                nf += 1;
            }
        }
        for (p, &c) in pcol.iter().enumerate() {
            pidx[c] = p;
        }
        let rcoef: Vec<Vec<(u32, S)>> = prow
            .iter()
            .map(|&ri| (0..q).filter(|&c| !used_c[c] && !g[ri][c].is_zero()).map(|c| (newidx[c], g[ri][c].clone())).collect())
            .collect();
        let nfu = nf as usize;
        {
            let chunk = t.len().div_ceil(nt).max(1);
            let (ni, pi, rc) = (&newidx, &pidx, &rcoef);
            std::thread::scope(|sc| {
                for tch in t.chunks_mut(chunk) {
                    sc.spawn(move || {
                        let mut acc: Vec<S> = vec![S::zero(prec); nfu];
                        let mut emax = vec![i32::MIN; nfu];
                        let mut touched: Vec<u32> = Vec::new();
                        let mut mark = vec![false; nfu];
                        let mut tmp = S::zero(prec);
                        for row in tch.iter_mut() {
                            for (c, x) in row.iter() {
                                let c = *c as usize;
                                if ni[c] != u32::MAX {
                                    let a = ni[c] as usize;
                                    acc[a].add_assign_s(x);
                                    emax[a] = emax[a].max(acc[a].exponent());
                                    if !mark[a] {
                                        mark[a] = true;
                                        touched.push(a as u32);
                                    }
                                }
                            }
                            for (c, x) in row.iter() {
                                let p = pi[*c as usize];
                                if p == usize::MAX {
                                    continue;
                                }
                                for (a, gv) in &rc[p] {
                                    let a = *a as usize;
                                    let e0 = acc[a].exponent();
                                    let ex = acc[a].sub_mul(x, gv, &mut tmp);
                                    emax[a] = emax[a].max(ex).max(e0);
                                    if !mark[a] {
                                        mark[a] = true;
                                        touched.push(a as u32);
                                    }
                                }
                            }
                            touched.sort_unstable();
                            let mut nr: SparseRow<S> = Vec::with_capacity(touched.len());
                            for &a in &touched {
                                let a = a as usize;
                                acc[a].snap(emax[a], keep);
                                if !acc[a].is_zero() {
                                    nr.push((a as u32, acc[a].clone()));
                                }
                                acc[a].set_zero();
                                emax[a] = i32::MIN;
                                mark[a] = false;
                            }
                            touched.clear();
                            *row = nr;
                        }
                    });
                }
            });
        }
        q = nfu;
        if S::SNAP && std::env::var_os("DISPERSIVE_VERBOSE").is_some() {
            eprintln!("  frobenius step {step} m {m}: rank {} q {q} T nnz {} ({:.1} s)", prow.len(), t.iter().map(|x| x.len()).sum::<usize>(), clock.elapsed().as_secs_f64());
        }
    }
    diag.steps = step;
    diag.q_final = q;
    diag.top = top;
    diag.seconds = clock.elapsed().as_secs_f64();
    FrobOut { t, q, ranks: ranks_out, pivots: pivots_out, diag }
}

/// the MPFR values of a shadowed run
pub fn shadowed_to_float(out: FrobOut<super::field::Shadowed>) -> FrobOut<Float> {
    FrobOut {
        t: out.t.into_iter().map(|row| row.into_iter().map(|(c, x)| (c, x.f)).collect()).collect(),
        q: out.q,
        ranks: out.ranks,
        pivots: out.pivots,
        diag: out.diag,
    }
}

/// structure of a family from the exact run (oracle prep5 family_spec + greedy)
#[derive(Clone, Debug)]
pub struct FamilySpec {
    pub name: String,
    pub lam: Rational,
    pub n_lo: i64,
    pub mm: i64,
    pub ranks: Vec<usize>,
    pub rho: Rational,
    pub dim: usize,
    pub k_stable: i64,
    pub first_order: i64,
    pub pivot_rows: Vec<(usize, i64)>,
    pub pole_rank: i64,
    /// pivot sequence of the exact run (structural pivoting of the MPFR run)
    pub pivots: Vec<Vec<(u32, u32)>>,
}

impl FamilySpec {
    pub fn to_json(&self) -> serde_json::Value {
        serde_json::json!({"name": self.name, "lam": self.lam.to_string(), "n_lo": self.n_lo, "M": self.mm, "ranks": self.ranks,
            "rho": self.rho.to_string(), "dim": self.dim, "k_stable": self.k_stable, "first_order": self.first_order,
            "pivot_rows": self.pivot_rows.iter().map(|(c, n)| serde_json::json!([c, n])).collect::<Vec<_>>(), "pole_rank": self.pole_rank})
    }
}

/// greedy pivot coordinates over the rows of the exact T (lowest orders first), rank by order, plateau
pub fn family_spec(out: &FrobOut<ModP>, e: usize, n_lo: i64, mm: i64, lam: &Rational, rho: &Rational, name: &str) -> Result<FamilySpec, String> {
    let q = out.q;
    let mut basis: Vec<(usize, Vec<u32>)> = Vec::new(); // (pivot column, fully reduced dense row)
    let mut picked: Vec<(usize, i64)> = Vec::new();
    let mut rbo: Vec<(i64, usize)> = Vec::new();
    let p = super::field::PRIME;
    let mut dim: Option<usize> = None;
    let nrows = out.t.len();
    for i in 0..nrows {
        let n = n_lo + (i / e) as i64;
        let mut v = vec![0u32; q];
        for (c, x) in &out.t[i] {
            v[*c as usize] = x.0;
        }
        for (c, b) in &basis {
            let f = v[*c] as u64;
            if f != 0 {
                for k in 0..q {
                    if b[k] != 0 {
                        v[k] = ((v[k] as u64 + p - f * b[k] as u64 % p) % p) as u32;
                    }
                }
            }
        }
        if let Some(c) = v.iter().position(|&x| x != 0) {
            let inv = ModP(v[c]).inv().0 as u64;
            for x in v.iter_mut() {
                *x = (*x as u64 * inv % p) as u32;
            }
            for (_, b) in basis.iter_mut() {
                let f = b[c] as u64;
                if f != 0 {
                    for k in 0..q {
                        if v[k] != 0 {
                            b[k] = ((b[k] as u64 + p - f * v[k] as u64 % p) % p) as u32;
                        }
                    }
                }
            }
            basis.push((c, v));
            picked.push((i % e, n));
        }
        if i % e == e - 1 {
            rbo.push((n, picked.len()));
            // plateau: 6 consecutive orders with the same positive rank
            if dim.is_none() {
                for w in rbo.windows(6) {
                    if w.iter().all(|x| x.1 == w[0].1) && w[0].1 > 0 {
                        dim = Some(w[0].1);
                        break;
                    }
                }
            }
            if let Some(d) = dim {
                if picked.len() > d {
                    break;
                }
            }
        }
    }
    let dim = dim.ok_or_else(|| format!("{name}: no rank plateau {rbo:?}"))?;
    let kst = rbo.iter().filter(|x| x.1 == dim).map(|x| x.0).max().unwrap();
    let first = rbo.iter().filter(|x| x.1 > 0).map(|x| x.0).min().unwrap();
    let piv: Vec<(usize, i64)> = picked.iter().filter(|x| x.1 <= kst).take(dim).cloned().collect();
    Ok(FamilySpec {
        name: name.to_string(),
        lam: lam.clone(),
        n_lo,
        mm,
        ranks: out.ranks.clone(),
        rho: rho.clone(),
        dim,
        k_stable: kst,
        first_order: first,
        pivot_rows: piv,
        pole_rank: out.diag.pole_rank,
        pivots: out.pivots.clone(),
    })
}

/// basis rows (n - n_lo) * E + comp, n <= k_stable: unit vectors on the pivot coordinates (oracle disp stable_basis)
pub fn stable_basis(out: &FrobOut<Float>, spec: &FamilySpec, e: usize, prec: u32) -> Result<Vec<Vec<Float>>, String> {
    let q = out.q;
    let dim = spec.pivot_rows.len();
    let mut u: Vec<Vec<Float>> = spec
        .pivot_rows
        .iter()
        .map(|&(c, n)| {
            let mut row = vec![Float::new(prec); q];
            for (k, x) in &out.t[((n - spec.n_lo) as usize) * e + c] {
                row[*k as usize] = x.clone();
            }
            row
        })
        .collect();
    let mut aug: Vec<Vec<Float>> = (0..dim).map(|i| (0..dim).map(|k| Float::with_val(prec, (i == k) as u32)).collect()).collect();
    let mut used_r = vec![false; dim];
    let mut used_c = vec![false; q];
    let mut sel_col = vec![0usize; dim];
    let mut tmp = Float::new(prec);
    for _ in 0..dim {
        let (mut bi, mut bc) = (usize::MAX, usize::MAX);
        let mut best = Float::new(prec);
        for i in 0..dim {
            if used_r[i] {
                continue;
            }
            for c in 0..q {
                if !used_c[c] && u[i][c].cmp_abs(&best) == Some(std::cmp::Ordering::Greater) {
                    best.clone_from(&u[i][c]);
                    bi = i;
                    bc = c;
                }
            }
        }
        if bi == usize::MAX || best.is_zero() {
            return Err(format!("{}: stable basis: singular pivot coordinates", spec.name));
        }
        used_r[bi] = true;
        used_c[bc] = true;
        sel_col[bi] = bc;
        let pv = u[bi][bc].clone();
        for c in 0..q {
            u[bi][c] /= &pv;
        }
        for c in 0..dim {
            aug[bi][c] /= &pv;
        }
        let (urow, arow) = (u[bi].clone(), aug[bi].clone());
        for i in 0..dim {
            if i == bi || u[i][bc].is_zero() {
                continue;
            }
            let f = u[i][bc].clone();
            for c in 0..q {
                if !urow[c].is_zero() {
                    use rug::Assign;
                    tmp.assign(&urow[c] * &f);
                    u[i][c] -= &tmp;
                }
            }
            for c in 0..dim {
                use rug::Assign;
                tmp.assign(&arow[c] * &f);
                aug[i][c] -= &tmp;
            }
        }
    }
    let mut colsel = vec![usize::MAX; q];
    for (i, &c) in sel_col.iter().enumerate() {
        colsel[c] = i;
    }
    let nrows = (((spec.k_stable - spec.n_lo + 1) as usize) * e).min(out.t.len());
    let rows: Vec<Vec<Float>> = out.t[..nrows]
        .iter()
        .map(|row| {
            let mut v = vec![Float::new(prec); dim];
            for (c, x) in row {
                let i = colsel[*c as usize];
                if i == usize::MAX {
                    continue;
                }
                for a in 0..dim {
                    v[a] += Float::with_val(prec, x * &aug[i][a]);
                }
            }
            v
        })
        .collect();
    Ok(rows)
}
