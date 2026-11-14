//! Periodic exact checkpoints INSIDE a Frobenius family run (W3b; the zero family and every infinity family take hours at 5 loops).
//!
//! A shed gives the process 2 s after SIGTERM, so a flush on the signal cannot be relied on: the family state after a completed order m
//! (sparse T, q, top, step, ranks and pivots so far, diagnostics, the active pivoting) is written atomically every `every_seconds`; a
//! rerun with the same fingerprint continues at the next order and produces bit-identical results (the remaining orders are computed
//! from the identical state; thread-count independent like the uninterrupted run).
//!
//! Format (binary, little endian): magic, fingerprint, payload, blake3 of everything before it. MPFR values exactly (precision, sign,
//! exponent, mantissa digits).
use super::field::Scalar;
use super::frobenius::{FrobDiag, Pivoting, SparseRow};
use std::io::Write;
use std::path::PathBuf;

const MAGIC: &[u8] = b"W3B-FROBENIUS-CKPT-1\n";

/// where and how often a family run checkpoints
#[derive(Clone, Debug)]
pub struct FrobCheckpoint {
    pub path: PathBuf,
    /// identity of the run: sample fingerprint + family name + precision + pivoting (a mismatching file is an ERROR, never reused)
    pub fingerprint: String,
    pub every_seconds: f64,
    /// test hook: write the checkpoint at order m and stop with an error (resume tests)
    #[doc(hidden)]
    pub stop_at_m: Option<i64>,
}

/// the recursion state after all orders < m_next
pub struct FrobState<S> {
    /// first and one-past-last order of the whole run (progress = (m_next - m_first) / (m_end - m_first))
    pub m_first: i64,
    pub m_end: i64,
    pub m_next: i64,
    pub t: Vec<SparseRow<S>>,
    pub q: usize,
    pub top: i64,
    pub step: usize,
    pub ranks: Vec<usize>,
    pub pivots: Vec<Vec<(u32, u32)>>,
    pub diag: FrobDiag,
    pub pivoting: Pivoting,
}

/// buffered writer that streams to the file in 64 MB pieces and hashes what it writes (no full in-memory copy of a GB-sized T)
struct W {
    buf: Vec<u8>,
    f: std::fs::File,
    h: blake3::Hasher,
}
impl W {
    fn spill(&mut self, force: bool) -> std::io::Result<()> {
        if force || self.buf.len() > (64 << 20) {
            self.h.update(&self.buf);
            self.f.write_all(&self.buf)?;
            self.buf.clear();
        }
        Ok(())
    }
    fn u8(&mut self, x: u8) {
        self.buf.push(x);
    }
    fn u32(&mut self, x: u32) {
        self.buf.extend_from_slice(&x.to_le_bytes());
    }
    fn u64(&mut self, x: u64) {
        self.buf.extend_from_slice(&x.to_le_bytes());
    }
    fn i64(&mut self, x: i64) {
        self.buf.extend_from_slice(&x.to_le_bytes());
    }
    fn f64(&mut self, x: f64) {
        self.buf.extend_from_slice(&x.to_le_bytes());
    }
    fn bytes(&mut self, b: &[u8]) {
        self.u64(b.len() as u64);
        self.buf.extend_from_slice(b);
    }
}
struct Rd<'a> {
    b: &'a [u8],
    p: usize,
}
impl Rd<'_> {
    fn take(&mut self, n: usize) -> Result<&[u8], String> {
        if self.p + n > self.b.len() {
            return Err("truncated Frobenius checkpoint".into());
        }
        let s = &self.b[self.p..self.p + n];
        self.p += n;
        Ok(s)
    }
    fn u8(&mut self) -> Result<u8, String> {
        Ok(self.take(1)?[0])
    }
    fn u32(&mut self) -> Result<u32, String> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn u64(&mut self) -> Result<u64, String> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
    fn i64(&mut self) -> Result<i64, String> {
        Ok(i64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
    fn f64(&mut self) -> Result<f64, String> {
        Ok(f64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
    fn bytes(&mut self) -> Result<Vec<u8>, String> {
        let n = self.u64()? as usize;
        Ok(self.take(n)?.to_vec())
    }
}

fn piv_code(p: Pivoting) -> (u8, i64) {
    match p {
        Pivoting::Complete => (0, 0),
        Pivoting::Threshold(t) => (1, t as i64),
        Pivoting::Structural => (2, 0),
    }
}

/// write the state atomically (tmp + fsync + rename + dir fsync); false if the scalar type cannot be checkpointed
pub fn save<S: Scalar>(ck: &FrobCheckpoint, st: &FrobState<S>) -> Result<bool, String> {
    let tmp = ck.path.with_extension(format!("tmp.{}", std::process::id()));
    let f = std::fs::File::create(&tmp).map_err(|e| format!("{}: {e}", tmp.display()))?;
    let mut w = W { buf: Vec::with_capacity(1 << 20), f, h: blake3::Hasher::new() };
    w.buf.extend_from_slice(MAGIC);
    w.bytes(ck.fingerprint.as_bytes());
    w.i64(st.m_first);
    w.i64(st.m_end);
    w.i64(st.m_next);
    w.u64(st.q as u64);
    w.i64(st.top);
    w.u64(st.step as u64);
    let (pc, pt) = piv_code(st.pivoting);
    w.u8(pc);
    w.i64(pt);
    w.u64(st.ranks.len() as u64);
    for r in &st.ranks {
        w.u64(*r as u64);
    }
    w.u64(st.pivots.len() as u64);
    for p in &st.pivots {
        w.u64(p.len() as u64);
        for (a, b) in p {
            w.u32(*a);
            w.u32(*b);
        }
    }
    let d = &st.diag;
    w.i64(d.pole_rank);
    for x in [d.min_rank_gap_log10, d.min_pivot_log10, d.min_structural_pivot_rel_log10, d.sum_structural_pivot_rel_log10, d.seconds] {
        w.f64(x);
    }
    w.u64(d.t_nnz_max as u64);
    w.i64(d.structural_fallback_step.map(|s| s as i64).unwrap_or(-1));
    w.u64(st.t.len() as u64);
    for row in &st.t {
        w.u32(row.len() as u32);
        for (c, v) in row {
            w.u32(*c);
            if !v.ckpt_put(&mut w.buf) {
                drop(w);
                let _ = std::fs::remove_file(&tmp);
                return Ok(false);
            }
        }
        w.spill(false).map_err(|e| e.to_string())?;
    }
    w.spill(true).map_err(|e| e.to_string())?;
    let h = w.h.finalize();
    w.f.write_all(h.as_bytes()).map_err(|e| e.to_string())?;
    w.f.sync_all().map_err(|e| e.to_string())?;
    drop(w);
    std::fs::rename(&tmp, &ck.path).map_err(|e| e.to_string())?;
    if let Some(dir) = ck.path.parent() {
        if let Ok(d) = std::fs::File::open(if dir.as_os_str().is_empty() { std::path::Path::new(".") } else { dir }) {
            let _ = d.sync_all();
        }
    }
    Ok(true)
}

/// header only (fingerprint, m_first, m_end, m_next) without reading the body: progress reports
pub fn peek(path: &std::path::Path) -> Option<(String, i64, i64, i64)> {
    use std::io::Read;
    let mut f = std::fs::File::open(path).ok()?;
    let mut head = vec![0u8; MAGIC.len() + 8];
    f.read_exact(&mut head).ok()?;
    if &head[..MAGIC.len()] != MAGIC {
        return None;
    }
    let n = u64::from_le_bytes(head[MAGIC.len()..].try_into().ok()?) as usize;
    let mut rest = vec![0u8; n + 24];
    f.read_exact(&mut rest).ok()?;
    let fp = String::from_utf8(rest[..n].to_vec()).ok()?;
    let g = |k: usize| i64::from_le_bytes(rest[n + 8 * k..n + 8 * k + 8].try_into().unwrap());
    Some((fp, g(0), g(1), g(2)))
}

/// Ok(None) = no file; Err = unreadable, corrupt or of another run (never reused silently)
pub fn load<S: Scalar>(ck: &FrobCheckpoint, prec: u32) -> Result<Option<FrobState<S>>, String> {
    let Ok(b) = std::fs::read(&ck.path) else { return Ok(None) };
    let name = ck.path.display();
    if b.len() < MAGIC.len() + 32 || &b[..MAGIC.len()] != MAGIC {
        return Err(format!("{name}: not a Frobenius checkpoint"));
    }
    let (body, h) = b.split_at(b.len() - 32);
    if blake3::hash(body).as_bytes() != h {
        return Err(format!("{name}: checksum mismatch (corrupt Frobenius checkpoint): remove it"));
    }
    let mut r = Rd { b: body, p: MAGIC.len() };
    let fp = String::from_utf8(r.bytes()?).map_err(|e| e.to_string())?;
    if fp != ck.fingerprint {
        return Err(format!("{name}: Frobenius checkpoint of another run ({fp} != {}); refusing to resume: use a fresh workdir", ck.fingerprint));
    }
    let m_first = r.i64()?;
    let m_end = r.i64()?;
    let m_next = r.i64()?;
    let q = r.u64()? as usize;
    let top = r.i64()?;
    let step = r.u64()? as usize;
    let pivoting = match (r.u8()?, r.i64()?) {
        (0, _) => Pivoting::Complete,
        (1, t) => Pivoting::Threshold(t as i32),
        (2, _) => Pivoting::Structural,
        x => return Err(format!("{name}: bad pivoting code {x:?}")),
    };
    let nr = r.u64()? as usize;
    let ranks = (0..nr).map(|_| r.u64().map(|x| x as usize)).collect::<Result<Vec<_>, _>>()?;
    let np = r.u64()? as usize;
    let mut pivots = Vec::with_capacity(np);
    for _ in 0..np {
        let n = r.u64()? as usize;
        pivots.push((0..n).map(|_| Ok((r.u32()?, r.u32()?))).collect::<Result<Vec<_>, String>>()?);
    }
    let mut diag = FrobDiag { pole_rank: r.i64()?, ..Default::default() };
    diag.min_rank_gap_log10 = r.f64()?;
    diag.min_pivot_log10 = r.f64()?;
    diag.min_structural_pivot_rel_log10 = r.f64()?;
    diag.sum_structural_pivot_rel_log10 = r.f64()?;
    diag.seconds = r.f64()?;
    diag.t_nnz_max = r.u64()? as usize;
    diag.structural_fallback_step = match r.i64()? {
        -1 => None,
        s => Some(s as usize),
    };
    let nt = r.u64()? as usize;
    let mut t = Vec::with_capacity(nt);
    for _ in 0..nt {
        let n = r.u32()? as usize;
        let mut row = Vec::with_capacity(n);
        for _ in 0..n {
            let c = r.u32()?;
            let v = S::ckpt_get(r.b, &mut r.p, prec).ok_or_else(|| format!("{name}: bad scalar"))?;
            row.push((c, v));
        }
        t.push(row);
    }
    if r.p != body.len() {
        return Err(format!("{name}: trailing bytes"));
    }
    Ok(Some(FrobState { m_first, m_end, m_next, t, q, top, step, ranks, pivots, diag, pivoting }))
}
