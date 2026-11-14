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

pub const CHECKPOINT_SCHEMA: &str = "rustflow.patch-kernel.checkpoint.v1";
/// Exit code of a run that stopped on a signal after flushing its checkpoint
/// (`EX_TEMPFAIL`; the campaign's "retry later").
pub const EXIT_RESUMABLE: i32 = 75;

/// What a checkpoint is valid for.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Fingerprint {
    pub kernel_version: String,
    pub system: String,
    pub path: String,
    pub bits: u32,
    pub mode: String,
    pub cauchy_pole: Option<String>,
    pub project_real_axis: bool,
    pub overlap_check: bool,
    /// blake3 of the exact initial vector.
    pub initial: String,
    /// blake3 of the running executable (information + optional check).
    pub binary: String,
}

impl Fingerprint {
    /// Fields that differ from `other` (the binary only if `check_binary`).
    pub fn mismatches(&self, other: &Self, check_binary: bool) -> Vec<&'static str> {
        let mut out = Vec::new();
        if self.kernel_version != other.kernel_version {
            out.push("kernel_version");
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
#[derive(Clone, Debug)]
pub struct CheckpointPolicy {
    pub file: std::path::PathBuf,
    /// Write after every `every_patches` completed patches (0 = never by count).
    pub every_patches: usize,
    /// ... or when `every_seconds` passed since the last write (0 = never by time).
    pub every_seconds: f64,
}

/// Register SIGTERM and SIGINT to set `flag` (the run then flushes its
/// checkpoint and the binary exits with [`EXIT_RESUMABLE`]).
#[cfg(unix)]
pub fn cancel_on_signals(flag: std::sync::Arc<std::sync::atomic::AtomicBool>) -> Result<()> {
    for sig in [signal_hook::consts::SIGTERM, signal_hook::consts::SIGINT] {
        signal_hook::flag::register(sig, flag.clone())?;
    }
    Ok(())
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
}
