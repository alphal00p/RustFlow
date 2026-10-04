use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use std::time::Instant;
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{
    AnalyticOriginOptions, CanonicalAlgebraicSystem, RootSeed, SquareRoot,
};
use symbolica_amflow::{
    ComplexFloat, Error, FlowOptions, Precision, Prescription, Progress, Result, RunContext,
    kinematics,
};
fn atom(s: &str) -> Result<Atom> {
    Atom::parse(s, "nonplanar108", Default::default())
        .map_err(|e| Error::InvalidInput(e.to_string()))
}
fn symbol(s: &str) -> Symbol {
    if let AtomView::Var(v) = atom(s).unwrap().as_view() {
        v.get_symbol()
    } else {
        panic!("symbol")
    }
}
fn evaluate(digits: u32, order: usize, denominator: i64) -> Result<Vec<Vec<ComplexFloat>>> {
    let fixture: Value =
        serde_json::from_str(include_str!("../fixtures/diffexp/nonplanar108.json")).unwrap();
    let p = Precision::decimal(digits)?;
    let all_clock = Instant::now();
    let x = symbol("x");
    let eps = symbol("eps");
    let vars = fixture["coordinate_names"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| symbol(v.as_str().unwrap()))
        .collect::<Vec<_>>();
    let mut coordinates = BTreeMap::new();
    for name in fixture["coordinate_names"].as_array().unwrap() {
        let name = name.as_str().unwrap();
        let a = atom(fixture["source_point"][name].as_str().unwrap())?;
        let b = atom(fixture["destination_point"][name].as_str().unwrap())?;
        coordinates.insert(symbol(name), &a + (b - &a) * Atom::var(x));
    }
    let path = kinematics::KinematicPath {
        parameter: x,
        coordinates,
    };
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
        let mut matrix = vec![vec![Atom::new(); 108]; 108];
        for entry in fixture["sparse_rational_coefficients"].as_array().unwrap() {
            let indices = entry["indices"].as_array().unwrap();
            if indices[0].as_u64().unwrap() as usize == index {
                matrix[indices[1].as_u64().unwrap() as usize - 1]
                    [indices[2].as_u64().unwrap() as usize - 1] =
                    atom(entry["coefficient"].as_str().unwrap())?;
            }
        }
        matrices.push(matrix);
    }
    let canonical = CanonicalAlgebraicSystem::new(
        eps,
        &vars,
        &letters,
        &matrices,
        vec![SquareRoot {
            symbol: symbol("r5"),
            radicand: atom(fixture["root_definitions"]["r5"].as_str().unwrap())?,
        }],
    )?;
    let system = canonical.pullback(&path, 4)?;
    let preparation_seconds = all_clock.elapsed().as_secs_f64();
    eprintln!(
        "PREPARE {preparation_seconds}s source_guards={} pulled_guards={}",
        canonical.nonzero_conditions().len(),
        system.nonzero_conditions.len()
    );
    let constants = ahash::HashMap::from_iter([
        (
            atom("pi_const")?,
            ComplexFloat::new(p.log(&p.i(-1)).im, p.real(0)),
        ),
        (atom("log3_const")?, p.log(&p.i(3))),
        (atom("i_const")?, p.complex(0, 1)),
    ]);
    let exact_boundary = fixture["boundary_by_epsilon"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| {
            row.as_array()
                .unwrap()
                .iter()
                .map(|v| atom(v.as_str().unwrap()))
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    let seeds = BTreeMap::from([(symbol("r5"), RootSeed::I0(Prescription::PlusI0))]);
    let offset = p.rational(&Rational::from((1, denominator)));
    let initialization_clock = Instant::now();
    let seed = system.analytic_origin(
        &exact_boundary,
        &constants,
        &seeds,
        p,
        &offset,
        &AnalyticOriginOptions {
            order,
            check_digits: 25,
            ..Default::default()
        },
    )?;
    let initialization_seconds = initialization_clock.elapsed().as_secs_f64();
    eprintln!(
        "INITIALIZE {initialization_seconds}s offset={} defect={}",
        offset, seed.maximum_defect
    );
    let compiled = system.compile(p)?;
    let interior = compiled
        .singularities()
        .iter()
        .filter(|z| {
            z.im < p.tolerance(25)
                && z.im > -p.tolerance(25)
                && z.re > offset.re
                && z.re < p.real(1)
        })
        .collect::<Vec<_>>();
    if !interior.is_empty() {
        return Err(Error::Unsupported(format!(
            "interior poles require explicit prescription: {interior:?}"
        )));
    }
    let options = FlowOptions {
        digits: 20,
        guard_digits: digits - 20,
        series_order: order,
        max_steps: 20000,
        ..Default::default()
    };
    let context = RunContext {
        progress: Some(Arc::new(|event| {
            if let Progress::Step { index } = event
                && index % 25 == 0
            {
                eprintln!("STEP {index}");
            }
        })),
        ..Default::default()
    };
    let clock = Instant::now();
    let answer = compiled.transport(
        &seed.boundary,
        &[p.i(1)],
        &seed.seeds,
        &options,
        &context,
        false,
    )?;
    let transport_seconds = clock.elapsed().as_secs_f64();
    let mut errors = Vec::new();
    for field in ["reference_by_epsilon", "original_by_epsilon"] {
        let mut maximum = p.real(0);
        for (k, row) in answer.solution.coefficients.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                let reference =
                    p.eval(&atom(fixture[field][k][j].as_str().unwrap())?, &constants)?;
                let error = p.norm(&p.sub(value, &reference));
                if error > maximum {
                    maximum = error;
                }
            }
        }
        errors.push(maximum);
    }
    let pass = errors.iter().all(|e| *e < p.tolerance(20));
    eprintln!(
        "full108 digits={digits} order={order} offset=1/{denominator}: preparation={preparation_seconds}s initialization={initialization_seconds}s transport={transport_seconds}s steps={} rejected={} source_error={} original_error={}",
        answer.solution.diagnostics.steps,
        answer.solution.diagnostics.rejected_steps,
        errors[0],
        errors[1]
    );
    if !pass {
        return Err(Error::Accuracy(format!(
            "108-master reference disagreement: {errors:?}"
        )));
    }
    Ok(answer.solution.coefficients)
}

#[test]
#[ignore = "full108 scientific benchmark; run explicitly in release mode (about18minutes)"]
fn full108_singular_origin_matches_original_and_independent_refinement() -> Result<()> {
    let low = evaluate(40, 48, 128)?;
    let high = evaluate(60, 64, 256)?;
    let p = Precision::decimal(80)?;
    assert_eq!(low.len(), 5);
    assert_eq!(high.len(), 5);
    for k in 0..5 {
        assert_eq!(low[k].len(), 108);
        assert_eq!(high[k].len(), 108);
        for i in 0..108 {
            assert!(
                p.norm(&p.sub(&low[k][i], &high[k][i])) < p.tolerance(20),
                "epsilon={k} master={}",
                i + 1
            );
        }
    }
    Ok(())
}
