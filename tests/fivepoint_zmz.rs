//! Complete75-master ZMZ benchmark from pinned2005.04195v2 ancillary data.
#[path = "common/fivepoint.rs"]
mod fivepoint;
use symbolica_amflow::{Precision, Result};
#[test]
#[ignore = "full75 two-loop supplied-boundary benchmark; run explicitly in release mode"]
fn full_zmz_refines_and_matches_original() -> Result<()> {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../fixtures/diffexp/fivepoint-zmz-full75.json"
    ))
    .unwrap();
    assert_eq!(
        fixture["original_master_indices_one_based"],
        serde_json::json!((1..=75).collect::<Vec<_>>())
    );
    let roots = &["tr5", "sqrtG3", "sqrtG3nc"];
    let low = fivepoint::evaluate_canonical(&fixture, roots, 60, 56)?;
    let high = fivepoint::evaluate_canonical(&fixture, roots, 70, 64)?;
    let p = Precision::decimal(90)?;
    assert_eq!(low.len(), 5);
    assert_eq!(high.len(), 5);
    for (k, (a, b)) in low.iter().zip(&high).enumerate() {
        assert_eq!(a.len(), 75);
        assert_eq!(b.len(), 75);
        for (j, (a, b)) in a.iter().zip(b).enumerate() {
            assert!(
                p.norm(&p.sub(a, b)) < p.tolerance(20),
                "refinement epsilon{k} master{}",
                j + 1
            );
        }
    }
    Ok(())
}
