//! Versioned reduction caches with native Symbolica atom serialization.
use crate::reduction::{Reduction, ReductionBackend};
use crate::{Error, Integral, IntegralFamily, Result, RunContext};
use serde::{Deserialize, Serialize};
use std::{
    io::Write,
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
};
use symbolica::prelude::*;

const VERSION: u32 = 1;
static NEXT_FILE: AtomicU64 = AtomicU64::new(0);

#[derive(Clone, Debug)]
pub struct CachedBackend<B> {
    pub backend: B,
    pub directory: PathBuf,
}

#[derive(Serialize, Deserialize)]
struct Entry {
    version: u32,
    key: String,
    digest: String,
    #[serde(with = "byte_buffer")]
    payload: Vec<u8>,
}

// Bincode's byte-buffer operation has the same wire representation as Vec<u8>,
// but copies whole buffers instead of routing every byte through a writer.
// Keep sequence decoding too, for the existing JSON cache representation.
mod byte_buffer {
    use serde::{
        Deserializer, Serializer,
        de::{SeqAccess, Visitor},
    };

    pub fn serialize<S: Serializer>(bytes: &[u8], serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_bytes(bytes)
    }
    pub fn deserialize<'de, D: Deserializer<'de>>(deserializer: D) -> Result<Vec<u8>, D::Error> {
        struct Bytes;
        impl<'de> Visitor<'de> for Bytes {
            type Value = Vec<u8>;
            fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                formatter.write_str("a byte buffer or byte sequence")
            }
            fn visit_bytes<E: serde::de::Error>(self, bytes: &[u8]) -> Result<Self::Value, E> {
                Ok(bytes.to_vec())
            }
            fn visit_byte_buf<E: serde::de::Error>(self, bytes: Vec<u8>) -> Result<Self::Value, E> {
                Ok(bytes)
            }
            fn visit_seq<A: SeqAccess<'de>>(
                self,
                mut sequence: A,
            ) -> Result<Self::Value, A::Error> {
                let mut bytes = Vec::new();
                while let Some(byte) = sequence.next_element()? {
                    bytes.push(byte);
                }
                Ok(bytes)
            }
        }
        deserializer.deserialize_byte_buf(Bytes)
    }
}
fn write_entry(file: &mut std::fs::File, entry: &Entry) -> Result<()> {
    // JSON serialization emits many small writes, especially for atom byte
    // arrays. Buffer them before touching the filesystem, then explicitly
    // flush and sync before the caller atomically renames the entry.
    let mut writer = std::io::BufWriter::new(file);
    serde_json::to_writer(&mut writer, entry).map_err(|e| Error::Cache(e.to_string()))?;
    writer.flush()?;
    writer.get_ref().sync_all()?;
    Ok(())
}
type StoredRule = (Vec<i16>, Vec<(Vec<i16>, Vec<u8>)>);

#[derive(Serialize, Deserialize)]
struct StoredReduction {
    rules: Vec<StoredRule>,
    residuals: Vec<Vec<i16>>,
    conditions: Vec<Vec<u8>>,
}
fn atom_bytes(a: &Atom) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    a.export(&mut bytes)?;
    Ok(bytes)
}

#[derive(Default)]
struct AtomDecoder {
    states: Vec<(Vec<u8>, symbolica::state::StateMap)>,
}
impl AtomDecoder {
    fn decode(&mut self, bytes: &[u8]) -> Result<Atom> {
        // Atom::export includes a state header for every coefficient. Large
        // tables repeat a handful of identical headers millions of times.
        // Import each complete header once; its symbol map remains valid as
        // Symbolica appends further symbols to the process-wide state.
        let index = if let Some(index) = self
            .states
            .iter()
            .rposition(|(prefix, _)| bytes.starts_with(prefix))
        {
            index
        } else {
            let mut source = bytes;
            let map = State::import(&mut source, None).map_err(|e| Error::Cache(e.to_string()))?;
            let prefix = bytes[..bytes.len() - source.len()].to_vec();
            if self.states.len() == 16 {
                self.states.remove(0);
            }
            self.states.push((prefix, map));
            self.states.len() - 1
        };
        let (prefix, map) = &self.states[index];
        let mut source = &bytes[prefix.len()..];
        if source.get(..8) != Some(1_u64.to_le_bytes().as_slice()) {
            // Streaming multi-term exports retain Symbolica's general reader.
            let mut complete = bytes;
            return Atom::import(&mut complete, None).map_err(|e| Error::Cache(e.to_string()));
        }
        source = &source[8..];
        let atom =
            Atom::import_with_map(&mut source, map).map_err(|e| Error::Cache(e.to_string()))?;
        if !source.is_empty() {
            return Err(Error::Cache("trailing cached Atom data".into()));
        }
        Ok(atom)
    }
}

impl StoredReduction {
    fn encode(r: &Reduction) -> Result<Self> {
        // All atoms already exist, so a single state snapshot covers every
        // symbol and coefficient-variable list they can reference. Retain the
        // ordinary self-contained Atom format for backward compatibility.
        let mut symbols = ahash::HashSet::default();
        for atom in r
            .rules
            .values()
            .flat_map(|terms| terms.values())
            .chain(&r.nonzero_conditions)
        {
            symbols.extend(atom.as_view().get_all_symbols(true));
        }
        let mut state = Vec::new();
        State::export_partial(&mut state, symbols)?;
        let encode = |atom: &Atom| -> Result<Vec<u8>> {
            let mut bytes = state.clone();
            bytes.extend_from_slice(&1_u64.to_le_bytes());
            atom.as_view().write(&mut bytes)?;
            Ok(bytes)
        };
        Ok(Self {
            rules: r
                .rules
                .iter()
                .map(|(i, terms)| {
                    Ok((
                        i.0.clone(),
                        terms
                            .iter()
                            .map(|(j, c)| Ok((j.0.clone(), encode(c)?)))
                            .collect::<Result<_>>()?,
                    ))
                })
                .collect::<Result<_>>()?,
            residuals: r.residuals.iter().map(|i| i.0.clone()).collect(),
            conditions: r
                .nonzero_conditions
                .iter()
                .map(encode)
                .collect::<Result<_>>()?,
        })
    }
    fn decode(self) -> Result<Reduction> {
        let mut decoder = AtomDecoder::default();
        let mut decode = |bytes: Vec<u8>| decoder.decode(&bytes);
        Ok(Reduction {
            rules: self
                .rules
                .into_iter()
                .map(|(i, terms)| {
                    Ok((
                        Integral(i),
                        terms
                            .into_iter()
                            .map(|(j, c)| Ok((Integral(j), decode(c)?)))
                            .collect::<Result<_>>()?,
                    ))
                })
                .collect::<Result<_>>()?,
            residuals: self.residuals.into_iter().map(Integral).collect(),
            nonzero_conditions: self
                .conditions
                .into_iter()
                .map(decode)
                .collect::<Result<_>>()?,
        })
    }
}

impl<B: ReductionBackend> CachedBackend<B> {
    fn key(&self, family: &IntegralFamily, targets: &[Integral]) -> Result<String> {
        // RustRed's typed family fingerprint includes dimension, all coefficients,
        // external Gram matrix, and scalar-product variable identities.
        let converted = family.convert()?;
        let fingerprint = converted.family.fingerprint();
        let mut targets = targets.iter().map(|i| i.0.clone()).collect::<Vec<_>>();
        targets.sort();
        targets.dedup();
        let mut atoms = Vec::new();
        for p in &family.propagators {
            atoms.push(p.constant.to_string());
            atoms.extend(p.scalar_products.iter().map(ToString::to_string));
        }
        atoms.extend(
            family
                .external_gram
                .iter()
                .flatten()
                .map(ToString::to_string),
        );
        let value = serde_json::json!([
            VERSION,
            env!("PORT_SOURCE_DIGEST"),
            env!("CARGO_PKG_VERSION"),
            env!("DEPENDENCY_SOURCE_DIGEST"),
            self.backend.identity(),
            fingerprint,
            family.physical_propagators,
            Atom::var(family.epsilon).to_string(),
            atoms,
            targets
        ]);
        Ok(
            blake3::hash(&serde_json::to_vec(&value).map_err(|e| Error::Cache(e.to_string()))?)
                .to_hex()
                .to_string(),
        )
    }
}

impl<B: ReductionBackend> ReductionBackend for CachedBackend<B> {
    fn reduce_cut(
        &self,
        family: &crate::cuts::CutFamily,
        targets: &[Integral],
        context: &RunContext,
    ) -> Result<Reduction> {
        CachedBackend {
            backend: crate::cuts::CutBackendView {
                backend: &self.backend,
                cuts: family.cuts(),
            },
            directory: self.directory.clone(),
        }
        .reduce(family.family(), targets, context)
    }
    fn reduce_cut_at_epsilon(
        &self,
        family: &crate::cuts::CutFamily,
        targets: &[Integral],
        epsilon: &Rational,
        context: &RunContext,
    ) -> Result<Reduction> {
        CachedBackend {
            backend: crate::reduction::SampledBackend {
                backend: &self.backend,
                epsilon: epsilon.clone(),
            },
            directory: self.directory.clone(),
        }
        .reduce_cut(family, targets, context)
    }
    fn reduce_at_epsilon(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        epsilon: &Rational,
        context: &RunContext,
    ) -> Result<Reduction> {
        CachedBackend {
            backend: crate::reduction::SampledBackend {
                backend: &self.backend,
                epsilon: epsilon.clone(),
            },
            directory: self.directory.clone(),
        }
        .reduce(family, targets, context)
    }
    fn identity(&self) -> String {
        self.backend.identity()
    }
    fn reduce(
        &self,
        family: &IntegralFamily,
        targets: &[Integral],
        context: &RunContext,
    ) -> Result<Reduction> {
        context.cancellation.check()?;
        let key = self.key(family, targets)?;
        let destination = self.directory.join(format!("{key}.json"));
        match std::fs::read(&destination) {
            Ok(bytes) => {
                let entry: Entry =
                    serde_json::from_slice(&bytes).map_err(|e| Error::Cache(e.to_string()))?;
                if entry.version != VERSION
                    || entry.key != key
                    || entry.digest != blake3::hash(&entry.payload).to_hex().as_str()
                {
                    return Err(Error::Cache("incompatible or corrupt cache entry".into()));
                }
                return serde_json::from_slice::<StoredReduction>(&entry.payload)
                    .map_err(|e| Error::Cache(e.to_string()))?
                    .decode();
            }
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => return Err(e.into()),
        }
        let result = self.backend.reduce(family, targets, context)?;
        let payload = serde_json::to_vec(&StoredReduction::encode(&result)?)
            .map_err(|e| Error::Cache(e.to_string()))?;
        let entry = Entry {
            version: VERSION,
            key: key.clone(),
            digest: blake3::hash(&payload).to_hex().to_string(),
            payload,
        };
        std::fs::create_dir_all(&self.directory)?;
        let temporary = self.directory.join(format!(
            ".{key}.{}.{}.tmp",
            std::process::id(),
            NEXT_FILE.fetch_add(1, Ordering::Relaxed)
        ));
        let write = (|| -> Result<()> {
            let mut file = std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&temporary)?;
            write_entry(&mut file, &entry)?;
            std::fs::rename(&temporary, &destination)?;
            Ok(())
        })();
        if write.is_err() {
            let _ = std::fs::remove_file(&temporary);
        }
        write?;
        Ok(result)
    }
}

#[derive(Serialize, Deserialize)]
struct StoredTransformation {
    previous: Vec<Vec<i16>>,
    current: Vec<Vec<i16>>,
    matrix: Vec<Vec<Vec<u8>>>,
}

#[derive(Serialize, Deserialize)]
struct StoredSystem {
    basis: Vec<Vec<i16>>,
    matrix: Vec<Vec<Vec<u8>>>,
    targets: Vec<StoredRule>,
    reduction: StoredReduction,
    refinement: Option<crate::refine::RefinementReport>,
    transformations: Vec<StoredTransformation>,
}

pub(crate) fn system_key(
    family: &IntegralFamily,
    targets: &[Integral],
    backend: &dyn ReductionBackend,
    options: &crate::FlowOptions,
) -> Result<String> {
    let value = serde_json::json!([
        VERSION,
        env!("PORT_SOURCE_DIGEST"),
        env!("DEPENDENCY_SOURCE_DIGEST"),
        backend.identity(),
        format!("{family:?}"),
        targets.iter().map(|i| &i.0).collect::<Vec<_>>(),
        format!("{:?}", options.mass_mode),
        options.refine_basis,
        options.skip_reduction,
        options.dimension
    ]);
    Ok(
        blake3::hash(&serde_json::to_vec(&value).map_err(|e| Error::Cache(e.to_string()))?)
            .to_hex()
            .to_string(),
    )
}

pub(crate) fn read_system(
    directory: &std::path::Path,
    key: &str,
) -> Result<
    Option<(
        crate::reduction::ReducedSystem,
        Option<crate::refine::RefinementReport>,
    )>,
> {
    let path = directory.join(format!("system-{key}.json"));
    let bytes = match std::fs::read(path) {
        Ok(bytes) => bytes,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(e) => return Err(e.into()),
    };
    let entry: Entry = serde_json::from_slice(&bytes).map_err(|e| Error::Cache(e.to_string()))?;
    if entry.version != VERSION
        || entry.key != key
        || entry.digest != blake3::hash(&entry.payload).to_hex().as_str()
    {
        return Err(Error::Cache(
            "incompatible or corrupt symbolic system".into(),
        ));
    }
    let stored: StoredSystem =
        serde_json::from_slice(&entry.payload).map_err(|e| Error::Cache(e.to_string()))?;
    let decode =
        |bytes: &[u8]| Atom::import(&mut &*bytes, None).map_err(|e| Error::Cache(e.to_string()));
    let reduced = stored.reduction.decode()?;
    let targets = stored
        .targets
        .into_iter()
        .map(|(_, terms)| {
            terms
                .into_iter()
                .map(|(i, c)| Ok((Integral(i), decode(&c)?)))
                .collect::<Result<_>>()
        })
        .collect::<Result<_>>()?;
    let matrix = stored
        .matrix
        .iter()
        .map(|r| {
            r.iter()
                .map(|bytes| decode(bytes))
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    if matrix.len() != stored.basis.len() || matrix.iter().any(|r| r.len() != stored.basis.len()) {
        return Err(Error::Cache("cached matrix dimensions".into()));
    }
    Ok(Some((
        crate::reduction::ReducedSystem {
            basis: stored.basis.into_iter().map(Integral).collect(),
            matrix,
            targets,
            nonzero_conditions: reduced.nonzero_conditions,
            candidates: reduced.rules,
            transformations: stored
                .transformations
                .into_iter()
                .map(|t| {
                    Ok(crate::reduction::BasisTransformation {
                        previous: t.previous.into_iter().map(Integral).collect(),
                        current: t.current.into_iter().map(Integral).collect(),
                        matrix: t
                            .matrix
                            .iter()
                            .map(|r| r.iter().map(|a| decode(a)).collect::<Result<_>>())
                            .collect::<Result<_>>()?,
                    })
                })
                .collect::<Result<_>>()?,
        },
        stored.refinement,
    )))
}

pub(crate) fn write_system(
    directory: &std::path::Path,
    key: &str,
    system: &crate::reduction::ReducedSystem,
    refinement: &Option<crate::refine::RefinementReport>,
) -> Result<()> {
    let reduction = Reduction {
        rules: system.candidates.clone(),
        residuals: vec![],
        nonzero_conditions: system.nonzero_conditions.clone(),
    };
    let stored = StoredSystem {
        basis: system.basis.iter().map(|i| i.0.clone()).collect(),
        matrix: system
            .matrix
            .iter()
            .map(|r| r.iter().map(atom_bytes).collect::<Result<_>>())
            .collect::<Result<_>>()?,
        targets: system
            .targets
            .iter()
            .map(|terms| {
                Ok((
                    vec![],
                    terms
                        .iter()
                        .map(|(i, c)| Ok((i.0.clone(), atom_bytes(c)?)))
                        .collect::<Result<_>>()?,
                ))
            })
            .collect::<Result<_>>()?,
        reduction: StoredReduction::encode(&reduction)?,
        refinement: refinement.clone(),
        transformations: system
            .transformations
            .iter()
            .map(|t| {
                Ok(StoredTransformation {
                    previous: t.previous.iter().map(|i| i.0.clone()).collect(),
                    current: t.current.iter().map(|i| i.0.clone()).collect(),
                    matrix: t
                        .matrix
                        .iter()
                        .map(|r| r.iter().map(atom_bytes).collect::<Result<_>>())
                        .collect::<Result<_>>()?,
                })
            })
            .collect::<Result<_>>()?,
    };
    let payload = serde_json::to_vec(&stored).map_err(|e| Error::Cache(e.to_string()))?;
    let entry = Entry {
        version: VERSION,
        key: key.into(),
        digest: blake3::hash(&payload).to_hex().to_string(),
        payload,
    };
    std::fs::create_dir_all(directory)?;
    let temporary = directory.join(format!(
        ".system-{key}.{}.{}.tmp",
        std::process::id(),
        NEXT_FILE.fetch_add(1, Ordering::Relaxed)
    ));
    let result = (|| -> Result<()> {
        let mut file = std::fs::OpenOptions::new()
            .create_new(true)
            .write(true)
            .open(&temporary)?;
        write_entry(&mut file, &entry)?;
        std::fs::rename(&temporary, directory.join(format!("system-{key}.json")))?;
        Ok(())
    })();
    if result.is_err() {
        let _ = std::fs::remove_file(temporary);
    }
    result
}

#[derive(Serialize, Deserialize)]
struct NativeStage {
    visited: Vec<Vec<i16>>,
    pending: Vec<Vec<i16>>,
    reduction: StoredReduction,
}
// The codec version is independent of the content-key version. Keeping the
// frontier first permits compact inspection without decoding all Atom data.
const NATIVE_MAGIC: &[u8] = b"AMFLOW-NATIVE\0\x01";

fn binary_decode<T: serde::de::DeserializeOwned>(bytes: &[u8]) -> Result<T> {
    let (value, consumed) = bincode::serde::decode_from_slice(bytes, bincode::config::standard())
        .map_err(|e| Error::Cache(e.to_string()))?;
    if consumed != bytes.len() {
        return Err(Error::Cache("trailing binary checkpoint data".into()));
    }
    Ok(value)
}
type NativeRestart = (Reduction, Vec<Integral>, Vec<Integral>);
pub(crate) fn read_native_stage(
    directory: &std::path::Path,
    key: &str,
) -> Result<Option<NativeRestart>> {
    let (bytes, binary) = match std::fs::read(directory.join(format!("native-stage-{key}.bin"))) {
        Ok(b) => (b, true),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            match std::fs::read(directory.join(format!("native-stage-{key}.json"))) {
                Ok(b) => (b, false),
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(None),
                Err(e) => return Err(e.into()),
            }
        }
        Err(e) => return Err(e.into()),
    };
    let entry: Entry = if binary {
        binary_decode(
            bytes
                .strip_prefix(NATIVE_MAGIC)
                .ok_or_else(|| Error::Cache("incompatible native checkpoint codec".into()))?,
        )?
    } else {
        serde_json::from_slice(&bytes).map_err(|e| Error::Cache(e.to_string()))?
    };
    if entry.version != VERSION
        || entry.key != key
        || entry.digest != blake3::hash(&entry.payload).to_hex().as_str()
    {
        return Err(Error::Cache(
            "incompatible native restart checkpoint".into(),
        ));
    }
    let stored: NativeStage = if binary {
        binary_decode(&entry.payload)?
    } else {
        serde_json::from_slice(&entry.payload).map_err(|e| Error::Cache(e.to_string()))?
    };
    Ok(Some((
        stored.reduction.decode()?,
        stored.visited.into_iter().map(Integral).collect(),
        stored.pending.into_iter().map(Integral).collect(),
    )))
}
pub(crate) fn write_native_stage(
    directory: &std::path::Path,
    key: &str,
    reduction: &Reduction,
    visited: &[Integral],
    pending: &[Integral],
) -> Result<()> {
    let payload = bincode::serde::encode_to_vec(
        &NativeStage {
            reduction: StoredReduction::encode(reduction)?,
            visited: visited.iter().map(|i| i.0.clone()).collect(),
            pending: pending.iter().map(|i| i.0.clone()).collect(),
        },
        bincode::config::standard(),
    )
    .map_err(|e| Error::Cache(e.to_string()))?;
    let entry = Entry {
        version: VERSION,
        key: key.into(),
        digest: blake3::hash(&payload).to_hex().to_string(),
        payload,
    };
    std::fs::create_dir_all(directory)?;
    let temporary = directory.join(format!(
        ".native-stage-{key}.{}.{}.tmp",
        std::process::id(),
        NEXT_FILE.fetch_add(1, Ordering::Relaxed)
    ));
    let write = (|| -> Result<()> {
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)?;
        let mut writer = std::io::BufWriter::new(&mut file);
        writer.write_all(NATIVE_MAGIC)?;
        bincode::serde::encode_into_std_write(&entry, &mut writer, bincode::config::standard())
            .map_err(|e| Error::Cache(e.to_string()))?;
        writer.flush()?;
        writer.get_ref().sync_all()?;
        std::fs::rename(
            &temporary,
            directory.join(format!("native-stage-{key}.bin")),
        )?;
        Ok(())
    })();
    if write.is_err() {
        let _ = std::fs::remove_file(temporary);
    }
    write
}

#[cfg(test)]
mod native_tests {
    use super::*;

    #[test]
    fn bulk_payload_codec_preserves_binary_and_json_wire_formats() {
        #[derive(Serialize, Deserialize)]
        struct LegacyEntry {
            version: u32,
            key: String,
            digest: String,
            payload: Vec<u8>,
        }
        let old = LegacyEntry {
            version: VERSION,
            key: "bulk".into(),
            digest: "test".into(),
            payload: (0..=255).cycle().take(65539).collect(),
        };
        let entry = Entry {
            version: old.version,
            key: old.key.clone(),
            digest: old.digest.clone(),
            payload: old.payload.clone(),
        };
        let config = bincode::config::standard();
        let legacy = bincode::serde::encode_to_vec(&old, config).unwrap();
        let bytes = bincode::serde::encode_to_vec(&entry, config).unwrap();
        assert_eq!(bytes, legacy);
        assert_eq!(
            binary_decode::<Entry>(&legacy).unwrap().payload,
            old.payload
        );
        assert_eq!(
            binary_decode::<LegacyEntry>(&bytes).unwrap().payload,
            old.payload
        );
        let json = serde_json::to_vec(&entry).unwrap();
        assert_eq!(json, serde_json::to_vec(&old).unwrap());
        assert_eq!(
            serde_json::from_slice::<Entry>(&json).unwrap().payload,
            old.payload
        );
    }

    #[test]
    fn shared_table_state_remains_compatible_with_symbolicas_atom_reader() {
        let coefficients = [
            parse!("table_state::x/(1+table_state::y)"),
            Atom::num((7, 3)),
        ];
        let reduction = Reduction {
            rules: std::collections::BTreeMap::from([(
                Integral(vec![1]),
                coefficients
                    .iter()
                    .enumerate()
                    .map(|(i, c)| (Integral(vec![i as i16 + 2]), c.clone()))
                    .collect(),
            )]),
            nonzero_conditions: vec![parse!("table_state::y")],
            ..Default::default()
        };
        let stored = StoredReduction::encode(&reduction).unwrap();
        for ((_, bytes), expected) in stored.rules[0].1.iter().zip(&coefficients) {
            assert_eq!(
                &Atom::import(&mut bytes.as_slice(), None).unwrap(),
                expected
            );
        }
        assert_eq!(
            Atom::import(&mut stored.conditions[0].as_slice(), None).unwrap(),
            reduction.nonzero_conditions[0]
        );
        assert_eq!(stored.decode().unwrap().rules, reduction.rules);
    }

    #[test]
    fn repeated_atom_headers_preserve_mixed_symbols_and_reject_trailing_data() {
        let atoms = [
            parse!("cache_codec::a+1/3"),
            parse!("cache_codec::f(cache_codec::b)^2"),
            Atom::num(7),
        ];
        let bytes = atoms
            .iter()
            .map(atom_bytes)
            .collect::<Result<Vec<_>>>()
            .unwrap();
        let mut decoder = AtomDecoder::default();
        for _ in 0..16 {
            for (expected, bytes) in atoms.iter().zip(&bytes) {
                assert_eq!(&decoder.decode(bytes).unwrap(), expected);
            }
        }
        assert!(decoder.states.len() <= atoms.len());
        let mut bad = bytes[0].clone();
        bad.push(0);
        assert!(matches!(decoder.decode(&bad), Err(Error::Cache(_))));
        assert!(matches!(
            decoder.decode(&bytes[0][..bytes[0].len() - 1]),
            Err(Error::Cache(_))
        ));
        let mut bad_header = bytes[0].clone();
        bad_header[0] ^= 1;
        assert!(matches!(decoder.decode(&bad_header), Err(Error::Cache(_))));
    }

    #[test]
    fn cached_atom_state_keeps_streaming_sum_compatibility() {
        let a = parse!("cache_stream::a/7");
        let b = parse!("cache_stream::b^2");
        let mut symbols = a.as_view().get_all_symbols(true);
        symbols.extend(b.as_view().get_all_symbols(true));
        let mut bytes = Vec::new();
        State::export_partial(&mut bytes, symbols).unwrap();
        bytes.extend_from_slice(&2_u64.to_le_bytes());
        a.as_view().write(&mut bytes).unwrap();
        b.as_view().write(&mut bytes).unwrap();
        assert_eq!(AtomDecoder::default().decode(&bytes).unwrap(), a + b);
    }
    use std::collections::BTreeMap;

    #[test]
    fn native_binary_checkpoint_migrates_legacy_and_authenticates_payload() {
        let directory =
            std::env::temp_dir().join(format!("amflow-native-codec-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&directory);
        std::fs::create_dir_all(&directory).unwrap();
        let key = "codec-roundtrip";
        let target = Integral(vec![2, 1, -3]);
        let master = Integral(vec![1, 1, 0]);
        let reduction = Reduction {
            rules: BTreeMap::from([(
                target.clone(),
                BTreeMap::from([(master.clone(), parse!("(1-codec_eps)/(codec_eta+1)"))]),
            )]),
            residuals: vec![master.clone()],
            nonzero_conditions: vec![parse!("codec_eta+1")],
        };
        let stored = NativeStage {
            visited: vec![target.0.clone()],
            pending: vec![master.0.clone()],
            reduction: StoredReduction::encode(&reduction).unwrap(),
        };
        let payload = serde_json::to_vec(&stored).unwrap();
        let legacy = Entry {
            version: VERSION,
            key: key.into(),
            digest: blake3::hash(&payload).to_hex().to_string(),
            payload,
        };
        let legacy_path = directory.join(format!("native-stage-{key}.json"));
        std::fs::write(&legacy_path, serde_json::to_vec(&legacy).unwrap()).unwrap();
        let (old, done, pending) = read_native_stage(&directory, key).unwrap().unwrap();
        assert_eq!(old.rules, reduction.rules);
        assert_eq!(old.nonzero_conditions, reduction.nonzero_conditions);
        write_native_stage(&directory, key, &old, &done, &pending).unwrap();
        let binary_path = directory.join(format!("native-stage-{key}.bin"));
        assert!(
            std::fs::metadata(&binary_path).unwrap().len()
                < std::fs::metadata(&legacy_path).unwrap().len()
        );
        let (new, new_done, new_pending) = read_native_stage(&directory, key).unwrap().unwrap();
        assert_eq!(new.rules, old.rules);
        assert_eq!(new.residuals, old.residuals);
        assert_eq!(new.nonzero_conditions, old.nonzero_conditions);
        assert_eq!(new_done, done);
        assert_eq!(new_pending, pending);
        let mut bytes = std::fs::read(&binary_path).unwrap();
        *bytes.last_mut().unwrap() ^= 1;
        std::fs::write(&binary_path, &bytes).unwrap();
        // A corrupt new entry must not silently fall back to the old entry.
        assert!(matches!(
            read_native_stage(&directory, key),
            Err(Error::Cache(_))
        ));
        write_native_stage(&directory, key, &old, &done, &pending).unwrap();
        let mut bytes = std::fs::read(&binary_path).unwrap();
        bytes.extend_from_slice(b"trailing");
        std::fs::write(&binary_path, &bytes).unwrap();
        assert!(matches!(
            read_native_stage(&directory, key),
            Err(Error::Cache(_))
        ));
        write_native_stage(&directory, key, &old, &done, &pending).unwrap();
        let mut bytes = std::fs::read(&binary_path).unwrap();
        bytes[NATIVE_MAGIC.len() - 1] = 2;
        std::fs::write(&binary_path, &bytes).unwrap();
        assert!(matches!(
            read_native_stage(&directory, key),
            Err(Error::Cache(_))
        ));
        std::fs::remove_dir_all(directory).unwrap();
    }
}
