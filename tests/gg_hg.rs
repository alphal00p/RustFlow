//! Complete supplied mixed QCD–EW Higgs-plus-jet systems; no Mathematica required.
#[path = "common/gg_hg.rs"]
mod gg_hg;
use symbolica_amflow::{Precision, Result};

#[test]
fn complete_planar_gg_hg_system_matches_original_and_restarts() -> Result<()> {
    let fixture = gg_hg::fixture();
    gg_hg::evaluate_case(&fixture, 0, 40, 64)?;
    Ok(())
}

#[test]
#[ignore = "full48/61 systems, four W/Z permutations and precision/order refinement; run in release mode"]
fn complete_crossed_gg_hg_systems_refine_and_match_original() -> Result<()> {
    let fixture = gg_hg::fixture();
    let p = Precision::decimal(120)?;
    assert_eq!(fixture["cases"].as_array().unwrap().len(), 16);
    let mut count = 0;
    let mut boundaries = Vec::new();
    for index in 0..16 {
        let low = gg_hg::evaluate_case(&fixture, index, 40, 64)?;
        let high = gg_hg::evaluate_case(&fixture, index, 60, 96)?;
        for (a, b) in low
            .coefficients
            .iter()
            .flatten()
            .zip(high.coefficients.iter().flatten())
        {
            assert!(
                p.norm(&p.sub(a, b)) < p.tolerance(24),
                "case{index} refinement"
            );
            count += 1;
        }
        boundaries.push(high);
    }
    assert_eq!(count, 4360);
    gg_hg::check_projections(&fixture, &boundaries)?;
    Ok(())
}
