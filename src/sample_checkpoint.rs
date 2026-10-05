//! Completed numerical samples are restart data, never verified boundaries.
use crate::{ComplexFloat, Error, FlowOptions, Integral, IntegralFamily, Precision, Result};
use serde::{Deserialize, Serialize};
use std::io::Write;
use std::sync::atomic::{AtomicU64, Ordering};
use symbolica::prelude::*;

const MAGIC: &[u8] = b"INTEGRAL-SAMPLE\0\x01";
static NEXT_FILE: AtomicU64 = AtomicU64::new(0);

#[derive(Serialize, Deserialize)]
struct Envelope {
    key: String,
    digest: String,
    payload: Vec<u8>,
}

#[derive(Serialize, Deserialize)]
struct CompletedSample {
    working_bits: u32,
    values: Vec<ComplexFloat>,
}

fn decode<T: serde::de::DeserializeOwned>(bytes: &[u8]) -> Result<T> {
    let (value, consumed) = bincode::serde::decode_from_slice(bytes, bincode::config::standard())
        .map_err(|e| Error::Cache(e.to_string()))?;
    if consumed != bytes.len() {
        return Err(Error::Cache("trailing completed-sample data".into()));
    }
    Ok(value)
}

/// Checkpoint a complete ordered vector, before algebraic projection or fitting.
/// The exact sample, precision, truncation settings and actual sources are keyed;
/// increasing precision cannot turn a prior result into independent validation.
pub(crate) fn evaluate(
    family: &IntegralFamily,
    integrals: &[Integral],
    epsilon: &Rational,
    options: &FlowOptions,
    backend: &dyn crate::ReductionBackend,
    source_identity: Option<&str>,
    compute: impl FnOnce() -> Result<Vec<ComplexFloat>>,
) -> Result<Vec<ComplexFloat>> {
    let Some(directory) = &options.sample_cache_directory else {
        return compute();
    };
    let mut settings = options.clone();
    settings.cache_directory = None;
    settings.sample_cache_directory = None;
    settings.reuse_samples = true;
    settings.workers = 1;
    let p = Precision::decimal(options.digits + options.guard_digits)?;
    let identity = serde_json::json!([
        "completed-integral-sample-v1",
        crate::cache::system_key(family, integrals, backend, options)?,
        format!("{settings:?}"),
        epsilon.to_string(),
        p.bits,
        source_identity,
    ]);
    let key =
        blake3::hash(&serde_json::to_vec(&identity).map_err(|e| Error::Cache(e.to_string()))?)
            .to_hex()
            .to_string();
    let path = directory.join(format!("sample-{key}.bin"));
    let validate = |entry: &CompletedSample| -> Result<()> {
        if entry.working_bits != p.bits
            || entry.values.len() != integrals.len()
            || entry
                .values
                .iter()
                .any(|v| !p.finite(v) || v.re.prec() < p.bits || v.im.prec() < p.bits)
        {
            return Err(Error::Cache(
                "completed sample dimensions, precision or values do not match".into(),
            ));
        }
        Ok(())
    };
    if options.reuse_samples {
        match std::fs::read(&path) {
            Ok(bytes) => {
                let envelope: Envelope =
                    decode(bytes.strip_prefix(MAGIC).ok_or_else(|| {
                        Error::Cache("incompatible completed-sample codec".into())
                    })?)?;
                if envelope.key != key
                    || envelope.digest != blake3::hash(&envelope.payload).to_hex().as_str()
                {
                    return Err(Error::Cache(
                        "completed-sample identity or checksum mismatch".into(),
                    ));
                }
                let entry: CompletedSample = decode(&envelope.payload)?;
                validate(&entry)?;
                return Ok(entry.values);
            }
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => return Err(e.into()),
        }
    }
    let entry = CompletedSample {
        working_bits: p.bits,
        values: compute()?,
    };
    validate(&entry)?;
    let payload = bincode::serde::encode_to_vec(&entry, bincode::config::standard())
        .map_err(|e| Error::Cache(e.to_string()))?;
    let envelope = Envelope {
        key: key.clone(),
        digest: blake3::hash(&payload).to_hex().to_string(),
        payload,
    };
    std::fs::create_dir_all(directory)?;
    let temporary = directory.join(format!(
        ".sample-{key}.{}.{}.tmp",
        std::process::id(),
        NEXT_FILE.fetch_add(1, Ordering::Relaxed)
    ));
    let write = (|| -> Result<()> {
        let file = std::fs::OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(&temporary)?;
        let mut writer = std::io::BufWriter::new(file);
        writer.write_all(MAGIC)?;
        bincode::serde::encode_into_std_write(&envelope, &mut writer, bincode::config::standard())
            .map_err(|e| Error::Cache(e.to_string()))?;
        writer.flush()?;
        writer.get_ref().sync_all()?;
        std::fs::rename(&temporary, &path)?;
        Ok(())
    })();
    if write.is_err() {
        let _ = std::fs::remove_file(&temporary);
    }
    write?;
    Ok(entry.values)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn samples_restart_without_promoting_accuracy_and_force_bypasses_only_samples() -> Result<()> {
        let directory = std::env::temp_dir().join(format!(
            "integral-samples-{}-{}",
            std::process::id(),
            NEXT_FILE.fetch_add(1, Ordering::Relaxed)
        ));
        let family = crate::IntegralFamily {
            name: "checkpoint_tadpole".into(),
            loops: vec!["l".into()],
            external: vec![],
            external_gram: vec![],
            propagators: vec![crate::Propagator {
                constant: Atom::num(-1),
                scalar_products: vec![Atom::num(1)],
            }],
            physical_propagators: 1,
            epsilon: symbol!("sample_checkpoint_eps"),
            dimension: 4,
        };
        let integrals = [Integral(vec![1])];
        let epsilon = Rational::from((1, 100));
        let mut options = FlowOptions {
            sample_cache_directory: Some(directory.clone()),
            ..Default::default()
        };
        let backend = crate::RustRedBackend::default();
        let p = Precision::decimal(options.digits + options.guard_digits)?;
        let first = vec![p.rational(&Rational::from((1, 3)))];
        assert_eq!(
            evaluate(
                &family,
                &integrals,
                &epsilon,
                &options,
                &backend,
                None,
                || Ok(first.clone())
            )?,
            first
        );
        assert_eq!(
            evaluate(
                &family,
                &integrals,
                &epsilon,
                &options,
                &backend,
                None,
                || Err(Error::Cancelled)
            )?,
            first
        );
        options.reuse_samples = false;
        assert!(matches!(
            evaluate(
                &family,
                &integrals,
                &epsilon,
                &options,
                &backend,
                None,
                || Err(Error::Cancelled)
            ),
            Err(Error::Cancelled)
        ));
        options.reuse_samples = true;
        assert_eq!(
            evaluate(
                &family,
                &integrals,
                &epsilon,
                &options,
                &backend,
                None,
                || unreachable!()
            )?,
            first
        );
        options.guard_digits += 10;
        assert!(matches!(
            evaluate(
                &family,
                &integrals,
                &epsilon,
                &options,
                &backend,
                None,
                || Err(Error::Cancelled)
            ),
            Err(Error::Cancelled)
        ));
        options.guard_digits -= 10;
        let file = std::fs::read_dir(&directory)?.next().unwrap()?.path();
        // A different supplied connection must never inherit samples from the
        // same physical family evaluated through another declaration.
        assert!(matches!(
            evaluate(
                &family,
                &integrals,
                &epsilon,
                &options,
                &backend,
                Some("connection-a"),
                || Err(Error::Cancelled)
            ),
            Err(Error::Cancelled)
        ));
        assert_eq!(
            evaluate(
                &family,
                &integrals,
                &epsilon,
                &options,
                &backend,
                Some("connection-a"),
                || Ok(first.clone())
            )?,
            first
        );
        assert!(matches!(
            evaluate(
                &family,
                &integrals,
                &epsilon,
                &options,
                &backend,
                Some("connection-b"),
                || Err(Error::Cancelled)
            ),
            Err(Error::Cancelled)
        ));
        assert_eq!(
            evaluate(
                &family,
                &integrals,
                &epsilon,
                &options,
                &backend,
                Some("connection-a"),
                || unreachable!()
            )?,
            first
        );
        let mut bytes = std::fs::read(&file)?;
        *bytes.last_mut().unwrap() ^= 1;
        std::fs::write(&file, bytes)?;
        assert!(matches!(
            evaluate(
                &family,
                &integrals,
                &epsilon,
                &options,
                &backend,
                None,
                || unreachable!()
            ),
            Err(Error::Cache(_))
        ));
        std::fs::remove_dir_all(directory)?;
        Ok(())
    }
}
