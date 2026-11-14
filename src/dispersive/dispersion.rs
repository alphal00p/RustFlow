//! Boundary fit at z = 0, head [0, s0], matching at S to the families at infinity, tail [S, inf), assembly of
//! I_j = int_0^inf d_j(s)/(s + 1) ds (oracle: disp3y-pt main()).
use super::boundary::BoundaryRow;
use super::frobenius::FamilySpec;
use rug::{Assign, Float, Rational, ops::Pow};

/// numerical family: basis rows (n - n_lo) * E + comp for n <= k_stable, columns = dim (unit on the pivot coordinates)
pub struct FamilyNum {
    pub spec: FamilySpec,
    pub basis: Vec<Vec<Float>>,
}

pub fn log10abs(x: &Float) -> f64 {
    if x.is_zero() {
        -1e9
    } else {
        x.clone().abs().log10().to_f64()
    }
}

/// least squares min ||A x - b|| by modified Gram-Schmidt QR with one reorthogonalisation; (x, residual norm)
pub fn lstsq(a: &[Vec<Float>], b: &[Float], prec: u32) -> (Vec<Float>, Float) {
    let m = a.len();
    let n = a[0].len();
    let mut q: Vec<Vec<Float>> = (0..n)
        .map(|j| (0..m).map(|i| a[i][j].clone()).collect())
        .collect();
    let mut rr = vec![vec![Float::new(prec); n]; n];
    let mut tmp = Float::new(prec);
    for j in 0..n {
        for _pass in 0..2 {
            for k in 0..j {
                tmp.assign(Float::dot(q[k].iter().zip(&q[j])));
                rr[k][j] += &tmp;
                let (qk, qj) = if k < j {
                    let (lo, hi) = q.split_at_mut(j);
                    (&lo[k], &mut hi[0])
                } else {
                    unreachable!()
                };
                for i in 0..m {
                    let t2 = Float::with_val(prec, &qk[i] * &tmp);
                    qj[i] -= &t2;
                }
            }
        }
        let nrm = Float::with_val(prec, Float::dot(q[j].iter().zip(&q[j]))).sqrt();
        rr[j][j] = nrm.clone();
        for i in 0..m {
            q[j][i] /= &nrm;
        }
    }
    let qtb: Vec<Float> = (0..n)
        .map(|j| Float::with_val(prec, Float::dot(q[j].iter().zip(b))))
        .collect();
    let mut x = vec![Float::new(prec); n];
    for j in (0..n).rev() {
        let mut v = qtb[j].clone();
        for k in j + 1..n {
            tmp.assign(&rr[j][k] * &x[k]);
            v -= &tmp;
        }
        v /= &rr[j][j];
        x[j] = v;
    }
    let mut res = Float::new(prec);
    for i in 0..m {
        let mut v = Float::with_val(prec, Float::dot(a[i].iter().zip(&x)));
        v -= &b[i];
        res += Float::with_val(prec, &v * &v);
    }
    (x, res.sqrt())
}

pub struct SoftFit {
    /// fitted coefficient vector c'_n (scaled, rows (n - n_lo) E + comp)
    pub csoft: Vec<Float>,
    pub rows: usize,
    pub residual_log10: f64,
    pub worst_rel_log10: f64,
}

/// weighted (relative) least squares of the boundary rows on the zero-family basis
pub fn soft_fit(fam0: &FamilyNum, bnd: &[BoundaryRow], e: usize, prec: u32) -> SoftFit {
    let mut amat = Vec::new();
    let mut bvec = Vec::new();
    let rho = Float::with_val(prec, &fam0.spec.rho);
    for br in bnd {
        let row = ((br.n - fam0.spec.n_lo) as usize) * e + br.comp;
        if row >= fam0.basis.len() {
            continue;
        }
        let sc = Float::with_val(prec, (&rho).pow(br.n as i32));
        let vs = Float::with_val(prec, &br.value * &sc);
        let w = if vs.is_zero() {
            Float::with_val(prec, 1)
        } else {
            Float::with_val(prec, 1 / Float::with_val(prec, vs.clone().abs()))
        };
        amat.push(
            fam0.basis[row]
                .iter()
                .map(|x| Float::with_val(prec, x * &w))
                .collect::<Vec<_>>(),
        );
        bvec.push(Float::with_val(prec, &vs * &w));
    }
    let (x0, res0) = lstsq(&amat, &bvec, prec);
    let mut worst = -1e9f64;
    for (i, row) in amat.iter().enumerate() {
        let mut v = Float::with_val(prec, Float::dot(row.iter().zip(&x0)));
        v -= &bvec[i];
        worst = worst.max(log10abs(&v));
    }
    let csoft = fam0
        .basis
        .iter()
        .map(|row| Float::with_val(prec, Float::dot(row.iter().zip(&x0))))
        .collect();
    SoftFit {
        csoft,
        rows: amat.len(),
        residual_log10: log10abs(&res0),
        worst_rel_log10: worst,
    }
}

/// d(s0) and the head integral int_0^s0 d(s)/(s + 1) ds termwise (d = s^lam sum_n c'_n (-s/rho)^n)
pub fn head(
    fam0: &FamilyNum,
    csoft: &[Float],
    s0q: &Rational,
    e: usize,
    lterms: i64,
    prec: u32,
) -> (Vec<Float>, Vec<Float>) {
    let s0 = Float::with_val(prec, s0q);
    let lam = Float::with_val(prec, &fam0.spec.lam);
    let rho = Float::with_val(prec, &fam0.spec.rho);
    let nrow0 = csoft.len() / e;
    let mut dv = vec![Float::new(prec); e];
    let mut hd = vec![Float::new(prec); e];
    let sig = -Float::with_val(prec, &s0 / &rho);
    let s0lam = Float::with_val(prec, (&s0).pow(&lam));
    for b in 0..nrow0 {
        let n = fam0.spec.n_lo + b as i64;
        let sp = Float::with_val(prec, (&sig).pow(n as i32));
        for c in 0..e {
            dv[c] += Float::with_val(prec, &csoft[b * e + c] * &sp);
        }
    }
    for x in dv.iter_mut() {
        *x *= &s0lam;
    }
    // sum_n c'_n (-1/rho)^n sum_l (-1)^l s0^(lam+n+l+1)/(lam+n+l+1)
    let mrho = Float::with_val(prec, -1 / Float::with_val(prec, &rho));
    for b in 0..nrow0 {
        let n = fam0.spec.n_lo + b as i64;
        let mut inner = Float::new(prec);
        for l in 0..lterms {
            let ex = Float::with_val(prec, &lam + (n + l + 1));
            let mut t = Float::with_val(prec, (&s0).pow(&ex));
            t /= &ex;
            if l % 2 == 1 {
                inner -= &t;
            } else {
                inner += &t;
            }
        }
        let f = Float::with_val(prec, Float::with_val(prec, (&mrho).pow(n as i32)) * &inner);
        for c in 0..e {
            hd[c] += Float::with_val(prec, &csoft[b * e + c] * &f);
        }
    }
    (dv, hd)
}

pub struct Matching {
    pub tail: Vec<Float>,
    pub coef: Vec<Float>,
    pub rows: usize,
    pub cols: usize,
    pub residual_log10: f64,
    pub worst_rel_log10: f64,
    pub max_log10_coef: f64,
}

/// one family's matching columns phi_a(S) = S^(-lam) sum_n B[n, ., a] (-1/(S rho))^n and tail vectors (int_S^inf phi_a/(s+1), analytically
/// continued in the exponents), one pair of E-vectors per basis column
pub fn family_columns(
    f: &FamilyNum,
    sq: &Rational,
    e: usize,
    tterms: i64,
    prec: u32,
) -> (Vec<Vec<Float>>, Vec<Vec<Float>>) {
    let s = Float::with_val(prec, sq);
    let nrow = f.basis.len() / e;
    let lam = Float::with_val(prec, &f.spec.lam);
    let rho = Float::with_val(prec, &f.spec.rho);
    let mlam = Float::with_val(prec, -&lam);
    let spow = Float::with_val(prec, (&s).pow(&mlam));
    let sig = Float::with_val(prec, -1 / Float::with_val(prec, &s * &rho));
    let mrho = Float::with_val(prec, -1 / Float::with_val(prec, &rho));
    let mut tf = Vec::with_capacity(nrow);
    for b in 0..nrow {
        let n = f.spec.n_lo + b as i64;
        let mut inner = Float::new(prec);
        for l in 0..tterms {
            let beta = Float::with_val(prec, &mlam - (n + l));
            let mut t = Float::with_val(prec, (&s).pow(&beta));
            t /= &beta;
            if l % 2 == 1 {
                inner += &t;
            } else {
                inner -= &t;
            }
        }
        tf.push(Float::with_val(
            prec,
            Float::with_val(prec, (&mrho).pow(n as i32)) * &inner,
        ));
    }
    let sp: Vec<Float> = (0..nrow)
        .map(|b| Float::with_val(prec, (&sig).pow((f.spec.n_lo + b as i64) as i32)))
        .collect();
    let mut cols = Vec::with_capacity(f.spec.dim);
    let mut tails = Vec::with_capacity(f.spec.dim);
    for a in 0..f.spec.dim {
        let mut col = vec![Float::new(prec); e];
        let mut tl = vec![Float::new(prec); e];
        for b in 0..nrow {
            for c in 0..e {
                let v = &f.basis[b * e + c][a];
                if v.is_zero() {
                    continue;
                }
                col[c] += Float::with_val(prec, v * &sp[b]);
                tl[c] += Float::with_val(prec, v * &tf[b]);
            }
        }
        for x in col.iter_mut() {
            *x *= &spow;
        }
        cols.push(col);
        tails.push(tl);
    }
    (cols, tails)
}

/// row-scaled least squares of d(S) on all families' columns (rows weighted by 1/max(|d_c(S)|, max_a |phi_a,c|)), then the tail
pub fn match_columns(
    cols: &[Vec<Float>],
    tails: &[Vec<Float>],
    ds: &[Float],
    e: usize,
    prec: u32,
) -> Matching {
    let ncol = cols.len();
    let mut amat = Vec::new();
    let mut bvec = Vec::new();
    for c in 0..e {
        let mut sc = ds[c].clone().abs();
        for col in cols {
            if col[c].cmp_abs(&sc) == Some(std::cmp::Ordering::Greater) {
                sc = col[c].clone().abs();
            }
        }
        if sc.is_zero() {
            continue;
        }
        let w = Float::with_val(prec, 1 / &sc);
        amat.push(
            (0..ncol)
                .map(|a| Float::with_val(prec, &cols[a][c] * &w))
                .collect::<Vec<_>>(),
        );
        bvec.push(Float::with_val(prec, &ds[c] * &w));
    }
    let (coef, resm) = lstsq(&amat, &bvec, prec);
    let mut worst = -1e9f64;
    for (i, row) in amat.iter().enumerate() {
        let mut v = Float::with_val(prec, Float::dot(row.iter().zip(&coef)));
        v -= &bvec[i];
        worst = worst.max(log10abs(&v));
    }
    let amax = coef.iter().map(log10abs).fold(-1e9f64, f64::max);
    let mut tail = vec![Float::new(prec); e];
    for a in 0..ncol {
        for c in 0..e {
            tail[c] += Float::with_val(prec, &tails[a][c] * &coef[a]);
        }
    }
    Matching {
        tail,
        coef,
        rows: amat.len(),
        cols: ncol,
        residual_log10: log10abs(&resm),
        worst_rel_log10: worst,
        max_log10_coef: amax,
    }
}

/// columns of every family + least squares (as the oracle)
pub fn match_at_infinity(
    fams: &[FamilyNum],
    ds: &[Float],
    sq: &Rational,
    e: usize,
    tterms: i64,
    prec: u32,
) -> Matching {
    let mut cols = Vec::new();
    let mut tails = Vec::new();
    for f in fams {
        let (c, t) = family_columns(f, sq, e, tterms, prec);
        cols.extend(c);
        tails.extend(t);
    }
    match_columns(&cols, &tails, ds, e, prec)
}

/// J(z = 1) = -(sin(pi d0/2)/pi) I
pub fn prefactor(d0: &Rational, prec: u32) -> Float {
    let pi = Float::with_val(prec, rug::float::Constant::Pi);
    let x = Float::with_val(
        prec,
        Float::with_val(prec, &pi * Float::with_val(prec, d0)) / 2u32,
    );
    let s = x.sin();
    -Float::with_val(prec, &s / &pi)
}

/// residual of the Frobenius recursion (lam + m) c_m - sum_k B_k c_{m-k} = 0 on the kept rows, for up to `max_cols` basis columns
/// (evenly spread), INCLUDING the leading equations m in [n_lo - r, n_lo) (pole constraints, c_m = 0 below n_lo): log10 of
/// max |residual| / max |terms| over the checked orders (m + pole rank <= k_stable)
pub fn family_residual(
    b: &super::laurent::Lmat<Float>,
    fam: &FamilyNum,
    e: usize,
    max_cols: usize,
    prec: u32,
) -> f64 {
    let r = -*b.keys().next().unwrap();
    let lam = Float::with_val(prec, &fam.spec.lam);
    let n_lo = fam.spec.n_lo;
    let nrow = (fam.basis.len() / e) as i64;
    let dim = fam.spec.dim;
    let ncol = dim.min(max_cols.max(1));
    let cols: Vec<usize> = (0..ncol).map(|k| k * dim / ncol).collect();
    let mut worst = -1e9f64;
    for &a in &cols {
        let coef = |n: i64, j: usize| -> Option<&Float> {
            if n < n_lo || n - n_lo >= nrow {
                None
            } else {
                Some(&fam.basis[((n - n_lo) as usize) * e + j][a])
            }
        };
        for m in (n_lo - r)..(n_lo + nrow - r) {
            let mut res = vec![Float::new(prec); e];
            let mut mag = Float::new(prec);
            for i in 0..e {
                if let Some(c) = coef(m, i) {
                    res[i] = Float::with_val(prec, c * Float::with_val(prec, &lam + m));
                    let t = res[i].clone().abs();
                    if t > mag {
                        mag = t;
                    }
                }
            }
            for (k, bk) in b {
                let jb = m - k;
                if jb < n_lo || jb - n_lo >= nrow {
                    continue;
                }
                for (i, j, v) in bk {
                    if let Some(c) = coef(jb, *j) {
                        let t = Float::with_val(prec, v * c);
                        let ta = t.clone().abs();
                        if ta > mag {
                            mag = ta;
                        }
                        res[*i] -= t;
                    }
                }
            }
            if mag.is_zero() {
                continue;
            }
            for x in &res {
                if !x.is_zero() {
                    worst = worst.max(log10abs(x) - log10abs(&mag));
                }
            }
        }
    }
    worst
}

/// W3b: the recursion residual of `family_residual` with two normalisations and its location, so that a real failure is visible
/// and a normalisation artefact can be told apart:
/// - `per_order`: max over orders m of |res_m| / max|terms of order m| (= `family_residual`): an order whose terms are all tiny
///   compared with the column (deep cancellation among the stable-basis combination) can show a large RELATIVE residual that is
///   harmless;
/// - `per_column`: max over m of |res_m| / max over ALL checked orders of |terms| for that column (the column's own scale in the
///   rho-scaled variable, where the terms of all orders are comparable by construction): a real failure of the recursion shows here.
///
/// Also: the order m and column of each worst value, the worst over the leading (pole) equations m < n_lo, and the log10 of the
/// largest term at the worst per-order m relative to the column scale.
pub fn family_residual_report(
    b: &super::laurent::Lmat<Float>,
    fam: &FamilyNum,
    e: usize,
    max_cols: usize,
    prec: u32,
) -> serde_json::Value {
    let r = -*b.keys().next().unwrap();
    let lam = Float::with_val(prec, &fam.spec.lam);
    let n_lo = fam.spec.n_lo;
    let nrow = (fam.basis.len() / e) as i64;
    let dim = fam.spec.dim;
    let ncol = dim.min(max_cols.max(1));
    let cols: Vec<usize> = (0..ncol).map(|k| k * dim / ncol).collect();
    let (mut w_ord, mut w_col, mut w_lead) = (
        (-1e9f64, 0i64, 0usize, 0.0f64),
        (-1e9f64, 0i64, 0usize),
        -1e9f64,
    );
    // per order m, over the checked columns: max log10 |res| and max log10 |terms| (for evaluation-weighted residuals)
    let mut orders: std::collections::BTreeMap<i64, (f64, f64)> = Default::default();
    for &a in &cols {
        let coef = |n: i64, j: usize| -> Option<&Float> {
            if n < n_lo || n - n_lo >= nrow {
                None
            } else {
                Some(&fam.basis[((n - n_lo) as usize) * e + j][a])
            }
        };
        let mut per_m: Vec<(i64, f64, f64)> = Vec::new(); // (m, log10 max|res|, log10 max|term|)
        for m in (n_lo - r)..(n_lo + nrow - r) {
            let mut res = vec![Float::new(prec); e];
            let mut mag = Float::new(prec);
            for i in 0..e {
                if let Some(c) = coef(m, i) {
                    res[i] = Float::with_val(prec, c * Float::with_val(prec, &lam + m));
                    let t = res[i].clone().abs();
                    if t > mag {
                        mag = t;
                    }
                }
            }
            for (k, bk) in b {
                let jb = m - k;
                if jb < n_lo || jb - n_lo >= nrow {
                    continue;
                }
                for (i, j, v) in bk {
                    if let Some(c) = coef(jb, *j) {
                        let t = Float::with_val(prec, v * c);
                        let ta = t.clone().abs();
                        if ta > mag {
                            mag = ta;
                        }
                        res[*i] -= t;
                    }
                }
            }
            if mag.is_zero() {
                continue;
            }
            let rmax = res.iter().map(log10abs).fold(-1e9f64, f64::max);
            per_m.push((m, rmax, log10abs(&mag)));
        }
        let scale = per_m.iter().map(|x| x.2).fold(-1e9f64, f64::max);
        for &(m, rm, mg) in &per_m {
            let o = orders.entry(m).or_insert((-1e9, -1e9));
            o.0 = o.0.max(rm);
            o.1 = o.1.max(mg);
            if rm > -1e8 {
                let ro = rm - mg;
                if ro > w_ord.0 {
                    w_ord = (ro, m, a, mg - scale);
                }
                let rc = rm - scale;
                if rc > w_col.0 {
                    w_col = (rc, m, a);
                }
                if m < n_lo {
                    w_lead = w_lead.max(ro);
                }
            }
        }
    }
    serde_json::json!({"per_order_log10": w_ord.0, "per_order_worst_m": w_ord.1, "per_order_worst_col": w_ord.2,
        "per_order_worst_m_scale_vs_column_log10": w_ord.3,
        "per_column_log10": w_col.0, "per_column_worst_m": w_col.1, "per_column_worst_col": w_col.2,
        "leading_equations_per_order_log10": w_lead, "cols_checked": cols.len(), "n_lo": n_lo, "pole_rank": r,
        "orders_m_res_terms_log10": orders.iter().map(|(m, (a, b))| serde_json::json!([m, (a * 10.0).round() / 10.0, (b * 10.0).round() / 10.0])).collect::<Vec<_>>()})
}

/// the residual of a report weighted for evaluation at t = x / rho (x = s for the zero family, 1/s at infinity): log10 of
/// max_m |res_m| t^m / max_m |terms_m| t^m, i.e. how badly the truncated series violates the recursion where it is USED
pub fn residual_at(report: &serde_json::Value, log10_t: f64) -> f64 {
    let Some(o) = report["orders_m_res_terms_log10"].as_array() else {
        return f64::NAN;
    };
    let (mut r, mut t) = (-1e9f64, -1e9f64);
    for x in o {
        let m = x[0].as_f64().unwrap_or(0.0);
        let (a, b) = (x[1].as_f64().unwrap_or(-1e9), x[2].as_f64().unwrap_or(-1e9));
        if a > -1e8 {
            r = r.max(a + m * log10_t);
        }
        if b > -1e8 {
            t = t.max(b + m * log10_t);
        }
    }
    r - t
}
