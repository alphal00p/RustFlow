//! Authenticated external mathematics for the optional Higgs-plus-jet workflows.
use crate::{Error, Result};
use std::{
    collections::BTreeMap,
    sync::{Arc, LazyLock, RwLock},
};

include!(concat!(env!("OUT_DIR"), "/gg_hg_data_manifest.rs"));
/// Data revision is independent of code revisions; payload hashes authenticate
/// the exact inputs used to build the library and its boundary fingerprints.
pub const SOURCE: &str = "https://raw.githubusercontent.com/alphal00p/RustFlow/c04ea606049007efb780683332a3440dc20317e7/fixtures/gg-hg/";
static LOADED: LazyLock<RwLock<BTreeMap<String, Arc<str>>>> =
    LazyLock::new(|| RwLock::new(BTreeMap::new()));

pub fn digest(name: &str) -> Result<&'static str> {
    FILES
        .iter()
        .find(|(file, _)| *file == name)
        .map(|(_, hash)| *hash)
        .ok_or_else(|| Error::InvalidInput(format!("unknown Higgs-jet data file: {name}")))
}

/// Authenticate every payload before publishing any of them. Invalid downloads
/// cannot change a live system's mathematical identity or its cache contents.
pub fn install(documents: BTreeMap<String, String>) -> Result<()> {
    for (name, text) in &documents {
        if blake3::hash(text.as_bytes()).to_hex().as_str() != digest(name)? {
            return Err(Error::InvalidInput(format!(
                "Higgs-jet data checksum mismatch: {name}"
            )));
        }
    }
    LOADED.write().unwrap().extend(
        documents
            .into_iter()
            .map(|(name, text)| (name, Arc::from(text))),
    );
    Ok(())
}

pub fn is_loaded(name: &str) -> bool {
    LOADED.read().unwrap().contains_key(name)
}

/// Rust callers may supply authenticated documents directly or explicitly name
/// a local data directory. Python/browser callers use the asynchronous loader.
pub(crate) fn get(name: &str) -> Result<Arc<str>> {
    digest(name)?;
    if let Some(text) = LOADED.read().unwrap().get(name) {
        return Ok(text.clone());
    }
    let directory = std::env::var_os("RUSTFLOW_HIGGS_JET_DATA_DIR").map(std::path::PathBuf::from);
    #[cfg(test)]
    let directory = directory.or_else(|| {
        Some(std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("fixtures/gg-hg"))
    });
    if let Some(directory) = directory {
        let text = std::fs::read_to_string(directory.join(name))
            .map_err(|error| Error::InvalidInput(format!("Higgs-jet data {name}: {error}")))?;
        install(BTreeMap::from([(name.to_owned(), text)]))?;
        return Ok(LOADED.read().unwrap()[name].clone());
    }
    Err(Error::InvalidInput(format!(
        "Higgs-jet data {name} has not been loaded; call await load_higgs_jet_data(form_factors=True) in Python, or install external data with gg_hg::data::install"
    )))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rejects_untrusted_data_before_installing_valid_entries() {
        let name = "physical-configurations.json";
        let text = std::fs::read_to_string(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("fixtures/gg-hg")
                .join(name),
        )
        .unwrap();
        let documents = BTreeMap::from([
            (name.to_owned(), text),
            ("generic-form-factors.json".to_owned(), "{}".to_owned()),
        ]);
        assert!(
            install(documents)
                .unwrap_err()
                .to_string()
                .contains("checksum mismatch")
        );
        assert!(digest("../../other.json").is_err());
    }
    #[test]
    fn external_inputs_keep_the_existing_mathematical_fingerprints() {
        let before = super::super::PluginBasisMap::input_fingerprint();
        let text = get("plugin-physical-map.json").unwrap();
        assert_eq!(before, blake3::hash(text.as_bytes()).to_hex().to_string());
        let system = super::super::HiggsJetIntegralSystem::load(
            super::super::PluginFamilyKind::Planar,
            "external_higgs_jet_data_test",
        )
        .unwrap();
        assert!(!system.mathematical_fingerprint().is_empty());
    }
}
