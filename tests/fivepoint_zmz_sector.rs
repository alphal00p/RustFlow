//! Genuine connected two-loop sector of the pinned 75-master ZMZ family.
#[path = "common/fivepoint.rs"]
mod fivepoint;
use symbolica_amflow::{Precision, Result};
#[test]
#[ignore = "connected two-loop benchmark; run explicitly in release mode"]
fn zmz_coupled_sector_refines_and_matches_original() -> Result<()> {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../fixtures/diffexp/fivepoint-zmz-sector14-15.json"
    ))
    .unwrap();
    assert_eq!(
        fixture["original_master_indices_one_based"],
        serde_json::json!([14, 15, 49, 50, 51, 52, 55, 56, 59, 60, 65, 68, 69])
    );
    let low = fivepoint::evaluate(&fixture, &["tr5", "sqrtG3"], 40, 48)?;
    let high = fivepoint::evaluate(&fixture, &["tr5", "sqrtG3"], 60, 80)?;
    let p = Precision::decimal(80)?;
    assert_eq!(low.len(), 5);
    assert_eq!(high.len(), 5);
    for (k, (a, b)) in low.iter().zip(&high).enumerate() {
        assert_eq!(a.len(), 13);
        assert_eq!(b.len(), 13);
        for (j, (a, b)) in a.iter().zip(b).enumerate() {
            assert!(
                p.norm(&p.sub(a, b)) < p.tolerance(20),
                "precision/order refinement epsilon{k} original master{}",
                fixture["original_master_indices_one_based"][j]
            );
        }
    }
    Ok(())
}
