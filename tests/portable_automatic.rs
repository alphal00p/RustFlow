#![cfg(all(feature = "wasm", feature = "automatic"))]

use symbolica::prelude::*;
use symbolica_amflow::*;

#[test]
fn portable_reduction_rejects_requested_threads_before_starting_search() {
    let family = IntegralFamily {
        name: "portable_tadpole".into(),
        loops: vec!["k".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![Propagator {
            constant: Atom::num(-1),
            scalar_products: vec![Atom::one()],
        }],
        physical_propagators: 1,
        epsilon: symbol!("portable_automatic::epsilon"),
        dimension: 4,
    };
    let backend = RustRedBackend {
        native_workers: 2,
        ..RustRedBackend::default()
    };
    for factorized in [true, false] {
        let backend = RustRedBackend {
            factorized,
            ..backend.clone()
        };
        let result = backend.reduce(&family, &[Integral(vec![2])], &RunContext::default());
        assert!(
            matches!(result, Err(Error::Unsupported(ref reason)) if reason.contains("one worker")),
            "{result:?}",
        );
    }
}
