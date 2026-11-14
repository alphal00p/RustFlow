//! Exact checkpoints of a patch transport.
//!
//! A checkpoint is the [`TransportState`](super::TransportState) after a
//! completed patch: patch index, solution vector, accumulated integral and
//! overlap probe as exact hexadecimal MPFR strings (precision kept), plus the
//! deterministic diagnostics. It carries a [`Fingerprint`]; resuming refuses
//! any mismatch (the binary hash only with an explicit override). The thread
//! count is not fingerprinted: results do not depend on it.
//!
//! Files are written atomically (temporary file, fsync, rename, directory
//! fsync), so a kill at any moment leaves either the previous or the new
//! checkpoint.
use crate::{Error, Result};
use rug::{Complex, Float};
use serde::{Deserialize, Serialize};
use std::io::Write;
use std::path::Path;

pub const CHECKPOINT_SCHEMA: &str = "rustflow.patch-kernel.checkpoint.v2";
/// Exit code of a run that stopped on a signal after flushing its checkpoint
/// (`EX_TEMPFAIL`; the campaign's "retry later").
pub const EXIT_RESUMABLE: i32 = 75;

/// What a checkpoint is valid for.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Fingerprint {
    pub kernel_version: String,
    /// Automatic algorithm identity: digest of the kernel sources and the
    /// resolved dependency graph (always checked).
    pub kernel_algorithm: String,
    pub system: String,
    pub path: String,
    /// Path points without orders (two-precision comparisons need equal points).
    pub path_points: String,
    /// Number of patches of the path (completeness of a run).
    pub path_len: usize,
    pub bits: u32,
    pub mode: String,
    pub cauchy_pole: Option<String>,
    pub project_real_axis: bool,
    pub overlap_check: bool,
    /// blake3 of the exact initial vector (at the working precision).
    pub initial: String,
    /// blake3 of the initial vector rounded to 64 bits: "the same initial
    /// data" across working precisions (two-precision comparisons).
    #[serde(default)]
    pub initial_64: String,
    /// blake3 of the running executable (information; checked only on request).
    pub binary: String,
}

impl Fingerprint {
    /// Fields that differ from `other` (the binary only if `check_binary`).
    pub fn mismatches(&self, other: &Self, check_binary: bool) -> Vec<&'static str> {
        let mut out = Vec::new();
        if self.kernel_version != other.kernel_version {
            out.push("kernel_version");
        }
        if self.kernel_algorithm != other.kernel_algorithm {
            out.push("kernel_algorithm");
        }
        if self.system != other.system {
            out.push("system");
        }
        if self.path != other.path {
            out.push("path");
        }
        if self.bits != other.bits {
            out.push("bits");
        }
        if self.mode != other.mode {
            out.push("mode");
        }
        if self.cauchy_pole != other.cauchy_pole {
            out.push("cauchy_pole");
        }
        if self.project_real_axis != other.project_real_axis {
            out.push("project_real_axis");
        }
        if self.overlap_check != other.overlap_check {
            out.push("overlap_check");
        }
        if self.initial != other.initial {
            out.push("initial");
        }
        if check_binary && self.binary != other.binary {
            out.push("binary");
        }
        out
    }
}

/// blake3 of the current executable ("unknown" if it cannot be read).
pub fn binary_hash() -> String {
    std::env::current_exe()
        .and_then(std::fs::read)
        .map(|b| blake3::hash(&b).to_hex().to_string())
        .unwrap_or_else(|_| "unknown".into())
}

/// Exact text form `[prec_re, prec_im, re_hex, im_hex]`.
pub type ExactComplex = (u32, u32, String, String);

fn float_hex(x: &Float) -> String {
    x.to_string_radix(16, None)
}

fn parse_hex(prec: u32, s: &str) -> Result<Float> {
    let parsed = Float::parse_radix(s, 16)
        .map_err(|e| Error::InvalidInput(format!("hex float {s:?}: {e}")))?;
    let x = Float::with_val(prec, parsed);
    Ok(x)
}

pub fn complex_to_exact(z: &Complex) -> ExactComplex {
    let (pr, pi) = z.prec();
    (pr, pi, float_hex(z.real()), float_hex(z.imag()))
}

pub fn complex_from_exact(e: &ExactComplex) -> Result<Complex> {
    let re = parse_hex(e.0, &e.2)?;
    let im = parse_hex(e.1, &e.3)?;
    let z = Complex::with_val((e.0, e.1), (re, im));
    // Exactness guard: the text must reproduce itself.
    if float_hex(z.real()) != e.2 || float_hex(z.imag()) != e.3 {
        return Err(Error::InvalidInput(
            "inexact hex float in checkpoint".into(),
        ));
    }
    Ok(z)
}

pub fn vector_to_exact(v: &[Complex]) -> Vec<ExactComplex> {
    v.iter().map(complex_to_exact).collect()
}

pub fn vector_from_exact(v: &[ExactComplex]) -> Result<Vec<Complex>> {
    v.iter().map(complex_from_exact).collect()
}

/// blake3 of a vector rounded to 64-bit components (precision independent).
pub fn vector_hash_64(v: &[Complex]) -> String {
    let r: Vec<Complex> = v.iter().map(|z| Complex::with_val(64, z)).collect();
    vector_hash(&r)
}

/// blake3 of an exact vector.
pub fn vector_hash(v: &[Complex]) -> String {
    let mut h = blake3::Hasher::new();
    for (pr, pi, re, im) in vector_to_exact(v) {
        h.update(format!("{pr};{pi};{re};{im}|").as_bytes());
    }
    h.finalize().to_hex().to_string()
}

/// Write `bytes` to `path` atomically.
pub fn write_atomic(path: &Path, bytes: &[u8]) -> Result<()> {
    let name = path
        .file_name()
        .ok_or_else(|| Error::InvalidInput("checkpoint path has no file name".into()))?
        .to_string_lossy()
        .into_owned();
    let tmp = path.with_file_name(format!(".{name}.tmp.{}", std::process::id()));
    {
        let mut f = std::fs::File::create(&tmp)?;
        f.write_all(bytes)?;
        f.sync_all()?;
    }
    std::fs::rename(&tmp, path)?;
    if let Some(dir) = path.parent() {
        let dir = if dir.as_os_str().is_empty() {
            Path::new(".")
        } else {
            dir
        };
        if let Ok(d) = std::fs::File::open(dir) {
            let _ = d.sync_all();
        }
    }
    Ok(())
}

/// When to write checkpoints during a run.
///
/// A run never overwrites a checkpoint that is not its own past: before
/// every write the existing file must carry the same fingerprint (binary
/// excepted) and a `next_patch` not beyond the state being written, unless
/// `fresh` is set (then only the first write of the run may replace it).
#[derive(Clone, Debug)]
pub struct CheckpointPolicy {
    pub file: std::path::PathBuf,
    /// Write after every `every_patches` completed patches (0 = never by count).
    pub every_patches: usize,
    /// ... or when `every_seconds` passed since the last write (0 = never by time).
    pub every_seconds: f64,
    /// Allow the first write to replace an unrelated or more advanced file.
    pub fresh: bool,
}

impl CheckpointPolicy {
    pub fn new(file: impl Into<std::path::PathBuf>) -> Self {
        Self {
            file: file.into(),
            every_patches: 1,
            every_seconds: 0.0,
            fresh: false,
        }
    }
}

/// Advisory exclusive lock on `<checkpoint>.lock`, held for a run's lifetime
/// so that two processes never write the same checkpoint.
pub struct CheckpointLock {
    _file: std::fs::File,
}

impl CheckpointLock {
    pub fn acquire(checkpoint: &Path) -> Result<Self> {
        let mut name = checkpoint.as_os_str().to_owned();
        name.push(".lock");
        let file = std::fs::OpenOptions::new()
            .create(true)
            .truncate(false)
            .write(true)
            .open(std::path::PathBuf::from(name))?;
        #[cfg(unix)]
        {
            use std::os::unix::io::AsRawFd;
            // SAFETY: flock on an open descriptor owned by `file`.
            let rc = unsafe { libc::flock(file.as_raw_fd(), libc::LOCK_EX | libc::LOCK_NB) };
            if rc != 0 {
                return Err(Error::InvalidInput(format!(
                    "checkpoint {} is locked by another run",
                    checkpoint.display()
                )));
            }
        }
        Ok(Self { _file: file })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hex_round_trip_is_exact() {
        for bits in [53u32, 200, 581, 643] {
            let mut x = Float::with_val(bits, 1) / 3u32;
            x *= Float::with_val(bits, 1e-300f64);
            let y = -Float::with_val(bits, 2).sqrt() * Float::with_val(bits, 7.25e250f64);
            let z = Complex::with_val(bits, (x, y));
            let e = complex_to_exact(&z);
            let back = complex_from_exact(&e).unwrap();
            assert_eq!(back, z);
            assert_eq!(back.prec(), z.prec());
        }
        let zero = Complex::new(100);
        assert_eq!(complex_from_exact(&complex_to_exact(&zero)).unwrap(), zero);
    }

    #[test]
    fn checkpoint_lock_is_exclusive() {
        let dir = std::env::temp_dir().join(format!("rf-pk-lock-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let file = dir.join("c.json");
        let a = CheckpointLock::acquire(&file).unwrap();
        assert!(CheckpointLock::acquire(&file).is_err());
        drop(a);
        assert!(CheckpointLock::acquire(&file).is_ok());
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
