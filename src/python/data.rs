use super::*;
#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyfunction(module = "symbolica.community.hep.integration")
)]
#[pyfunction]
pub(super) fn higgs_jet_data_manifest() -> BTreeMap<String, (String, String, bool)> {
    crate::gg_hg::data::FILES
        .iter()
        .map(|&(name, digest)| {
            (
                name.to_owned(),
                (
                    format!("{}{name}", crate::gg_hg::data::SOURCE),
                    digest.to_owned(),
                    crate::gg_hg::data::is_loaded(name),
                ),
            )
        })
        .collect()
}

#[cfg_attr(
    feature = "python_stubgen",
    pyo3_stub_gen::derive::gen_stub_pyfunction(module = "symbolica.community.hep.integration")
)]
#[pyfunction]
pub(super) fn install_higgs_jet_data(documents: BTreeMap<String, String>) -> PyResult<()> {
    crate::gg_hg::data::install(documents).map_err(error)
}
