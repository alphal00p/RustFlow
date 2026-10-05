use feynkit_generator::{GenerationOptions, Process};
use feynkit_model::Model;
use idenso::tensor::{AlgebraContraction, AlgebraSettings, ContractSettings, SymbolicTensor};
use spenso::structure::partial::PartialStructure;
use std::{collections::BTreeMap, sync::Arc};
use symbolica::prelude::*;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    let model_path = args.get(1).ok_or("provide the scoped native model JSON")?;
    let output = args.get(2).ok_or("provide an output directory")?;
    let base = std::path::Path::new(&output);
    std::fs::create_dir_all(base)?;
    let model = Arc::new(Model::from_json(&std::fs::read_to_string(model_path)?)?);
    eprintln!(
        "native model loaded: {} particles, {} vertices",
        model.particles().len(),
        model.vertex_rules().len()
    );
    let process = Process::new(["G", "G"], ["G", "H"]).with_filters(
        vec![],
        Some(vec!["GGGHEWZZ".into(), "GGGHEWWW".into()]),
        vec![],
    );
    let options = GenerationOptions::default()
        .with_loop_count(0, 0)?
        .threads(1);
    let amplitude = process.generate_amplitude(model, &options, Default::default())?;
    for (index, graph) in amplitude.diagrams().iter().enumerate() {
        std::fs::write(base.join(format!("generated-{index}.dot")), graph.to_dot()?)?;
    }
    eprintln!(
        "generated {} effective diagrams",
        amplitude.diagrams().len()
    );
    let expression = amplitude.expression();
    std::fs::write(
        base.join("amplitude.atom"),
        expression.to_canonical_string(),
    )?;
    eprintln!(
        "native amplitude {} bytes, {} external legs",
        expression.to_canonical_string().len(),
        amplitude.legs().len()
    );
    let sums = amplitude
        .squared()?
        .sum_spins(&[0, 1, 2, 3], true, &BTreeMap::new(), &BTreeMap::new())?
        .sum_colors(&[0, 1, 2, 3], true)?;
    eprintln!(
        "squared and state-summed {} bytes",
        sums.expression().to_canonical_string().len()
    );
    let tensor = SymbolicTensor::<PartialStructure>::infer(sums.expression().clone())?;
    let contracted = tensor.contract(ContractSettings::default().without_rank_one_tensors())?;
    std::fs::write(
        base.join("contracted.atom"),
        contracted.expression().to_canonical_string(),
    )?;
    eprintln!(
        "contracted {} bytes",
        contracted.expression().to_canonical_string().len()
    );
    assert!(tensor.is_scalar());
    let expanded = tensor.expanded(None, false)?;
    eprintln!(
        "expanded tensor {} bytes",
        expanded.expression().to_canonical_string().len()
    );
    let simplified = expanded.simplify_algebra(&AlgebraSettings {
        contract: AlgebraContraction::Dots,
        ..AlgebraSettings::hep()
    })?;
    assert!(simplified.is_scalar());
    std::fs::write(
        base.join("scalar.atom"),
        simplified.expression().to_canonical_string(),
    )?;
    eprintln!(
        "native scalar simplification {} bytes, status {:?}",
        simplified.expression().to_canonical_string().len(),
        simplified.reduction_status()
    );
    Ok(())
}
