use serde_json::Value;
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
    time::Instant,
};
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{CanonicalAlgebraicSystem, SquareRoot};
use symbolica_amflow::contour::PolynomialPrescription;
use symbolica_amflow::physical_transport::PhysicalRoute;
use symbolica_amflow::transport_cache::*;
use symbolica_amflow::{
    Error, FlowOptions, Precision, Prescription, Progress, Result, RunContext, RustFlow,
};
fn atom(s: &str) -> Result<Atom> {
    Atom::parse(s, "fivepoint_prescribed_cache", Default::default())
        .map_err(|e| Error::InvalidInput(e.to_string()))
}
fn symbol(s: &str) -> Symbol {
    if let AtomView::Var(v) = atom(s).unwrap().as_view() {
        v.get_symbol()
    } else {
        panic!("symbol")
    }
}
fn coordinates(f: &Value, key: &str) -> Result<BTreeMap<Symbol, Atom>> {
    f["coordinate_names"]
        .as_array()
        .unwrap()
        .iter()
        .map(|n| {
            let n = n.as_str().unwrap();
            Ok((symbol(n), atom(f[key][n].as_str().unwrap())?))
        })
        .collect()
}
fn values(p: Precision, f: &Value) -> Result<Vec<Vec<symbolica_amflow::ComplexFloat>>> {
    f.as_array()
        .unwrap()
        .iter()
        .map(|r| {
            r.as_array()
                .unwrap()
                .iter()
                .map(|v| {
                    p.parse(
                        v["real"].as_str().unwrap(),
                        v["imaginary"].as_str().unwrap(),
                    )
                })
                .collect()
        })
        .collect()
}
#[test]
#[ignore = "full13 physical prescribed-cache acceptance; bounded explicit run"]
fn full13_ph1_to_ph6_through_the_prescribed_physical_cache() -> Result<()> {
    let clock = Instant::now();
    let fixture: Value = serde_json::from_str(include_str!(
        "../fixtures/diffexp/fivepoint-planar-1loop.json"
    ))
    .unwrap();
    let variables = fixture["coordinate_names"]
        .as_array()
        .unwrap()
        .iter()
        .map(|n| symbol(n.as_str().unwrap()))
        .collect::<Vec<_>>();
    let used = fixture["sparse_rational_coefficients"]
        .as_array()
        .unwrap()
        .iter()
        .map(|e| e["indices"][0].as_u64().unwrap() as usize)
        .collect::<BTreeSet<_>>();
    let mut letters = Vec::new();
    let mut matrices = Vec::new();
    for index in used {
        letters.push(atom(fixture["letters"][index - 1].as_str().unwrap())?);
        let mut matrix = vec![vec![Atom::new(); 13]; 13];
        for entry in fixture["sparse_rational_coefficients"].as_array().unwrap() {
            let ix = entry["indices"].as_array().unwrap();
            if ix[0].as_u64().unwrap() as usize == index {
                matrix[ix[1].as_u64().unwrap() as usize - 1]
                    [ix[2].as_u64().unwrap() as usize - 1] =
                    atom(entry["coefficient"].as_str().unwrap())?;
            }
        }
        matrices.push(matrix);
    }
    let roots = ["tr5", "sqrtG3"]
        .into_iter()
        .map(|name| {
            Ok(SquareRoot {
                symbol: symbol(name),
                radicand: atom(
                    fixture["root_definitions"][name]
                        .as_str()
                        .unwrap()
                        .strip_prefix("sqrt(")
                        .unwrap()
                        .strip_suffix(')')
                        .unwrap(),
                )?,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let system = CanonicalAlgebraicSystem::new(
        symbol("eps"),
        &variables,
        &letters,
        &matrices,
        roots.clone(),
    )?;
    let mut declarations = variables
        .iter()
        .map(|&s| PolynomialPrescription {
            polynomial: Atom::var(s),
            prescription: Prescription::PlusI0,
        })
        .collect::<Vec<_>>();
    declarations.extend(roots.iter().map(|r| PolynomialPrescription {
        polynomial: r.radicand.clone(),
        prescription: Prescription::PlusI0,
    }));
    let basis = (1..=13)
        .map(|i| atom(&format!("I{i}")))
        .collect::<Result<Vec<_>>>()?;
    let flow=RustFlow::new_canonical(system,&basis,&Atom::one(),Prescription::PlusI0,"supplied PH1 physical normalization and named planar basis")?.with_prescribed_continuation(PhysicalContinuation{prescriptions:declarations,unprescribed_side:Prescription::PlusI0,domain:"PH1 to PH6 fixture: physical invariants and tr5/sqrtG3 polynomials +i0, no extra winding; all unprescribed singularities use upper parameter bypass".into()})?;
    eprintln!(
        "prescribed-cache identity prepared at {:.3}s",
        clock.elapsed().as_secs_f64()
    );
    let p = Precision::decimal(160)?;
    let range = EpsilonRange::new(0, 4)?;
    let germ = RootGerm {
        sheets: roots
            .iter()
            .map(|r| (r.symbol, RootSheet::Principal))
            .collect(),
    };
    let errors = fixture["boundary_absolute_errors_by_epsilon"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| {
            row.as_array()
                .unwrap()
                .iter()
                .map(|v| p.parse(v.as_str().unwrap(), "0").map(|z| z.re))
                .collect()
        })
        .collect::<Result<Vec<Vec<_>>>>()?;
    let mut cache = RustFlowCache::default();
    cache.insert(CachedBoundary{identity:flow.identity().clone(),point:CachedPoint::Exact(coordinates(&fixture,"source_point")?).with_root_germ(germ.clone())?,kind:PointKind::Physical,range,coefficients:values(p,&fixture["boundary_by_epsilon"])? ,accuracy:BoundaryAccuracy::supplied(128,p.bits,errors,"supplied original boundary; 132 absolute source digits, conservative 1e-130 inexact errors; no precision inference")?})?;
    let options = FlowOptions {
        digits: 20,
        guard_digits: 30,
        series_order: 64,
        max_steps: 5000,
        ..Default::default()
    };
    let context = RunContext {
        progress: Some(Arc::new(|event| {
            if let Progress::Step { index } = event
                && index % 100 == 0
            {
                eprintln!("prescribed-cache step {index}");
            }
        })),
        ..Default::default()
    };
    let admission = |_: &CachedBoundary, _: &CachedPoint, route: &PhysicalRoute| {
        eprintln!(
            "prescribed-cache route: {} crossings, {} exact vertices",
            route.crossings.len(),
            route.waypoints.len()
        );
        eprintln!(
            "prescribed-cache exact_route_json={}",
            serde_json::json!({
                "parameter": Atom::var(route.path.parameter).to_canonical_string(),
                "coordinates": route.path.coordinates.iter().map(|(symbol, value)| (Atom::var(*symbol).to_canonical_string(), value.to_canonical_string())).collect::<BTreeMap<_,_>>(),
                "vertices": route.waypoints.iter().map(Atom::to_canonical_string).collect::<Vec<_>>(),
            })
        );
        Ok(true)
    };
    let result = flow.evaluate_prescribed_to(
        &mut cache,
        &coordinates(&fixture, "destination_point")?,
        Some(&germ),
        range,
        &options,
        &context,
        &ScaledDistance {
            scales: BTreeMap::new(),
            admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
        },
        &admission,
    )?;
    eprintln!(
        "prescribed-cache accepted_result_json={}",
        serde_json::json!({
            "coefficients_by_epsilon": result.boundary.coefficients.iter().map(|row| row.iter().map(|value| serde_json::json!({"real": value.re.to_string(), "imaginary": value.im.to_string()})).collect::<Vec<_>>()).collect::<Vec<_>>(),
            "propagated_comparison_errors_by_epsilon": result.boundary.accuracy.comparison_errors().iter().map(|row|row.iter().map(ToString::to_string).collect::<Vec<_>>()).collect::<Vec<_>>(),
            "verified_digits": result.boundary.accuracy.verified_digits(),
            "working_bits": result.boundary.accuracy.working_bits(),
            "cached_points": cache.len(),
            "accepted_steps": result.transport.as_ref().unwrap().diagnostics.steps,
            "rejected_steps": result.transport.as_ref().unwrap().diagnostics.rejected_steps,
        })
    );
    for key in ["reference_by_epsilon", "original_by_epsilon"] {
        let expected = values(p, &fixture[key])?;
        let mut maximum = p.real(0);
        for (a, b) in result
            .boundary
            .coefficients
            .iter()
            .flatten()
            .zip(expected.iter().flatten())
        {
            let error = p.norm(&p.sub(a, b));
            if error > maximum {
                maximum = error;
            }
        }
        eprintln!("prescribed-cache {key} absolute error {maximum}");
        assert!(maximum < p.tolerance(20));
    }
    eprintln!(
        "prescribed-cache full13 passed at {:.3}s; {} cached points",
        clock.elapsed().as_secs_f64(),
        cache.len()
    );
    Ok(())
}
