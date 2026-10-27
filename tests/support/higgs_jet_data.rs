//! Supply repository fixtures through the same authenticated API as external callers.
use std::{collections::BTreeMap, path::Path, sync::Once};
use symbolica_amflow::gg_hg::data;

pub fn install() {
    static INSTALLED: Once = Once::new();
    INSTALLED.call_once(|| {
        let directory = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/gg-hg");
        let documents = data::FILES
            .iter()
            .map(|&(name, _)| {
                let text = std::fs::read_to_string(directory.join(name))
                    .expect("read Higgs-plus-jet test fixture");
                (name.to_owned(), text)
            })
            .collect::<BTreeMap<_, _>>();
        data::install(documents).expect("authenticate Higgs-plus-jet test fixtures");
    });
}
