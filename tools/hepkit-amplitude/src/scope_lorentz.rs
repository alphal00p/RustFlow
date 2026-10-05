mod scope;
use std::{collections::BTreeMap, path::Path};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    let directory = args.get(1).ok_or("provide a scope directory")?;
    let base = Path::new(directory);
    let inputs: BTreeMap<String, String> =
        serde_json::from_slice(&std::fs::read(base.join("original-lorentz.json"))?)?;
    let mut output = BTreeMap::new();
    let mut counts = BTreeMap::new();
    for (name, expression) in inputs {
        let (converted, count) = scope::preserve_momentum_denominators(&expression)?;
        if count > 0 {
            output.insert(name.clone(), converted);
            counts.insert(name, count);
        }
    }
    std::fs::write(
        base.join("scoped-lorentz.json"),
        serde_json::to_vec_pretty(&output)?,
    )?;
    std::fs::write(
        base.join("scope-counts.json"),
        serde_json::to_vec_pretty(&counts)?,
    )?;
    println!("preserved original momentum denominator scopes: {counts:?}");
    Ok(())
}
