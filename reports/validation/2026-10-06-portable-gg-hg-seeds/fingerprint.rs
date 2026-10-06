use std::path::Path;
fn field(hash: &mut blake3::Hasher, bytes: &[u8]) {
    hash.update(&(bytes.len() as u64).to_le_bytes());
    hash.update(bytes);
}
fn main() {
    let root = std::env::args().nth(1).expect("fixture directory");
    for topology in ["Planar_EW1", "NP_EW1"] {
        let mut hash = blake3::Hasher::new();
        for text in [
            "symbolica-hep-integration:higgs-jet-mathematics:v1",
            topology,
            "D=4-2*eps;+i0;measure=exp(2*eps*EulerGamma);mV2=mu2=1;basis=ordered-plugin-canonical;root-sheets=principal-or-opposite",
        ] { field(&mut hash, text.as_bytes()); }
        for name in ["integral-systems.json", "plugin-physical-map.json", "physical-configurations.json"] {
            field(&mut hash, name.as_bytes());
            field(&mut hash, &std::fs::read(Path::new(&root).join(name)).unwrap());
        }
        println!("{topology} {}", hash.finalize().to_hex());
    }
}
