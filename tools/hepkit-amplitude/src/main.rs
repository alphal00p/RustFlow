fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    let input = args
        .get(1)
        .ok_or("provide a native-compatible UFO directory")?;
    let output = args.get(2).ok_or("provide an output model JSON path")?;
    let model = pyo3::Python::attach(|py| {
        feynkit_ufo::UfoLoader::new()
            .restriction_name("full")
            .simplify_model(false)
            .load(py, input)
    })?;
    std::fs::write(output, model.model.to_json_pretty()?)?;
    println!(
        "Native UfoLoader exported {} particles and {} vertices",
        model.diagnostics.particle_count, model.diagnostics.vertex_rule_count
    );
    Ok(())
}
