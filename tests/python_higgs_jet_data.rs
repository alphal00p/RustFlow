#![cfg(feature = "python_api")]

use pyo3::{prelude::*, types::PyModule};
use std::collections::BTreeMap;

#[test]
fn community_data_loader_exports_authenticate_before_installing() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let module = PyModule::new(py, "symbolica.community.hep_integration_native")?;
        symbolica_amflow::python::register(&module)?;
        let exports: Vec<String> = module.getattr("__all__")?.extract()?;
        for name in ["higgs_jet_data_manifest", "install_higgs_jet_data"] {
            assert!(exports.iter().any(|export| export == name));
        }
        let manifest = || -> PyResult<BTreeMap<String, (String, String, bool)>> {
            module
                .getattr("higgs_jet_data_manifest")?
                .call0()?
                .extract()
        };
        let before = manifest()?;
        assert_eq!(before.len(), 7);
        let fixtures = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/gg-hg");
        for (name, (url, digest, loaded)) in &before {
            let text = std::fs::read(fixtures.join(name)).unwrap();
            assert_eq!(*digest, blake3::hash(&text).to_hex().to_string());
            assert_eq!(
                url,
                &format!("{}{name}", symbolica_amflow::gg_hg::data::SOURCE)
            );
            assert!(!loaded);
        }
        // The valid entry sorts first, so an incremental installer would expose
        // it before discovering the second entry's invalid checksum.
        let name = "integral-systems.json";
        let text = std::fs::read_to_string(fixtures.join(name)).unwrap();
        let install = module.getattr("install_higgs_jet_data")?;
        let error = install
            .call1((BTreeMap::from([
                (name, text.as_str()),
                ("physical-configurations.json", "{}"),
            ]),))
            .unwrap_err();
        assert!(error.is_instance(py, &module.getattr("InvalidInputError")?));
        assert!(error.to_string().contains("checksum mismatch"));
        assert_eq!(
            manifest()?,
            before,
            "invalid batch partially installed data"
        );

        install.call1((BTreeMap::from([(name, text)]),))?;
        let after = manifest()?;
        for (file, (url, digest, loaded)) in &after {
            assert_eq!((url, digest), (&before[file].0, &before[file].1));
            assert_eq!(*loaded, file == name);
        }
        Ok(())
    })
    .unwrap();
}
