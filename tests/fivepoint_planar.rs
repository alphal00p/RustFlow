//! Full planar one-loop PH1-to-PH6 benchmark.
#[path = "common/fivepoint.rs"]
mod fivepoint;
use symbolica_amflow::{Precision, Result};
#[test]
#[ignore = "full13-master physical benchmark; run explicitly in release mode"]
fn full_fivepoint_ph1_to_ph6_refines_and_matches_original() -> Result<()> {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../fixtures/diffexp/fivepoint-planar-1loop.json"
    ))
    .unwrap();
    let economical = fivepoint::evaluate(&fixture, &["tr5", "sqrtG3"], 40, 48)?;
    let low = fivepoint::evaluate(&fixture, &["tr5", "sqrtG3"], 60, 80)?;
    let high = fivepoint::evaluate(&fixture, &["tr5", "sqrtG3"], 80, 112)?;
    let p = Precision::decimal(90)?;
    assert_eq!(economical.len(), 5);
    assert_eq!(low.len(), 5);
    assert_eq!(high.len(), 5);
    for (k, (low, high)) in low.iter().zip(&high).enumerate() {
        assert_eq!(economical[k].len(), 13);
        assert_eq!(low.len(), 13);
        assert_eq!(high.len(), 13);
        for (j, (a, b)) in low.iter().zip(high).enumerate() {
            assert!(
                p.norm(&p.sub(&economical[k][j], b)) < p.tolerance(20),
                "economical-profile refinement epsilon{k} master{}",
                j + 1
            );
            assert!(
                p.norm(&p.sub(a, b)) < p.tolerance(40),
                "refinement epsilon{k} master{}",
                j + 1
            );
        }
    }
    Ok(())
}
