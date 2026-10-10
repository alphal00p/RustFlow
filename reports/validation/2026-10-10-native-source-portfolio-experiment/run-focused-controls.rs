fn main() {
    let passed = rustred::solver::run_guarded_portfolio_controls();
    assert_eq!(passed.len(), 19);
    let output = serde_json::json!({
        "status": "passed", "passed": passed.len(), "failed": 0,
        "controls": passed, "scope": "Isolated point-only native portfolio branch controls; no physical values or full closure.",
        "production_modified": false,
    });
    if let Some(path) = std::env::args().nth(1) {
        std::fs::write(path, serde_json::to_vec_pretty(&output).unwrap()).unwrap();
    }
    println!("{}", output);
}
