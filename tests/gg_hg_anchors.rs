//! Reference-loading checks only; neither reference is an AMF seed provider.
#[path = "support/gg_hg_anchors.rs"]
mod anchors;
use anchors::{EuclideanAnchor, first_planar_canonical};
use symbolica_amflow::gg_hg::{HiggsJetIntegralSystem, PluginFamilyKind};
use symbolica_amflow::symbolica::prelude::*;
use symbolica_amflow::{Precision, Result};

#[test]
fn euclidean_anchor_conventions_and_unknown_taylor_remainder_are_explicit() -> Result<()> {
    let p = Precision::decimal(70)?;
    for kind in [PluginFamilyKind::Planar, PluginFamilyKind::Nonplanar] {
        let input = HiggsJetIntegralSystem::load(kind, "gg_hg_anchor_test")?;
        let anchor = EuclideanAnchor::load(&input, p)?;
        assert_eq!(anchor.root_germ.sheets.len(), input.basis_map().roots.len());
        assert_eq!(anchor.verified_digits, 40);
        let at_zero = anchor.truncated(&Rational::from(0), p);
        assert_eq!(at_zero.values, anchor.coefficients[0]);
        assert_eq!(at_zero.coefficient_errors, anchor.absolute_errors[0]);
        assert_eq!(at_zero.first_unknown_power, 5);
        assert!(anchor.provenance.contains("sha256="));
    }
    Ok(())
}

#[test]
fn first_planar_all_epsilon_formula_matches_fourth_order_anchor_with_expected_remainder()
-> Result<()> {
    let p = Precision::decimal(70)?;
    let input = HiggsJetIntegralSystem::load(PluginFamilyKind::Planar, "gg_hg_analytic_check")?;
    let anchor = EuclideanAnchor::load(&input, p)?;
    // Independent 80-digit Mathematica evaluation retained in the first-native
    // component report; this checks native MPFR gamma normalization as well.
    let oracle = p.parse(
        "1.03948953194649787895621722388065767777646246176884141120465612461445885426733977500448504",
        "0",
    )?;
    let analytic = first_planar_canonical(&Rational::from((1, 101)), &Rational::from((-1, 50)), p)?;
    assert!(p.close(&analytic, &oracle, 65));
    let small = Rational::from((1, 1001));
    let large = &small * &Rational::from(2);
    let difference = |e: &Rational| -> Result<Float> {
        let exact = first_planar_canonical(e, &Rational::from((-1, 50)), p)?;
        Ok(p.norm(&p.sub(&exact, &anchor.truncated(e, p).values[0])))
    };
    let ratio = difference(&large)? / difference(&small)?;
    // A missing normalization or low-order coefficient would scale as eps^0..4;
    // the first unknown fifth-order term instead scales as 2^5.
    assert!(ratio > p.real(31) && ratio < p.real(33));
    Ok(())
}
