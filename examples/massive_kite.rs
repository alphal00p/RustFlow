//! Cold, single-worker AMFlow seed for the finite two-equal-mass kite.
use std::{collections::BTreeMap, sync::Arc, time::Instant};
use symbolica::prelude::*;
use symbolica_amflow::{
    transport_cache::{RustFlowCache, ScaledDistance},
    *,
};

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let number = |index: usize, fallback: u32| {
        args.get(index).map_or(Ok(fallback), |arg| {
            arg.parse()
                .map_err(|_| Error::InvalidInput("expected integer option".into()))
        })
    };
    let rho = symbol!("massive_kite_benchmark::rho");
    let gram = vec![vec![-Atom::var(rho)]];
    let family = IntegralFamily {
        name: "two_equal_mass_kite".into(),
        loops: vec!["k".into(), "l".into()],
        external: vec!["p".into()],
        propagators: vec![
            Propagator::quadratic(&[1, 0], &[0], Atom::num(1), &gram)?,
            Propagator::quadratic(&[1, 0], &[-1], Atom::num(0), &gram)?,
            Propagator::quadratic(&[0, 1], &[0], Atom::num(number(8, 1)? as i64), &gram)?,
            Propagator::quadratic(&[0, 1], &[-1], Atom::num(0), &gram)?,
            Propagator::quadratic(&[1, -1], &[0], Atom::num(0), &gram)?,
        ],
        external_gram: gram,
        physical_propagators: 5,
        epsilon: symbol!("massive_kite_benchmark::epsilon"),
        dimension: 4,
    };
    let options = FlowOptions {
        mass_mode: MassMode::from_selection(args.get(1).map(String::as_str), None)?,
        digits: number(2, 8)?,
        guard_digits: number(3, 24)?,
        series_order: number(4, 70)? as usize,
        workers: 1,
        sampled_reduction: args.get(6).is_none_or(|value| value != "symbolic"),
        residual_arithmetic: args
            .get(7)
            .map_or(Ok(ResidualArithmetic::Ball), |v| v.parse())?,
        ..Default::default()
    };
    let started = Instant::now();
    let trace = std::env::var_os("AMFLOW_BENCH_TRACE").is_some();
    let context = RunContext {
        progress: Some(Arc::new(move |event| {
            if trace || matches!(event, Progress::Sample { .. } | Progress::Prepared { .. }) {
                eprintln!("{:.3}s {event:?}", started.elapsed().as_secs_f64());
            }
        })),
        ..Default::default()
    };
    let backend = RustRedBackend::default();
    let prepared = PreparedPhysicalFamily::new(
        &family,
        &[Integral(vec![1; 5])],
        &[rho],
        &backend,
        &options,
        "positive Euclidean virtuality; two squared masses equal one",
        &context,
    )?;
    let preparation_seconds = started.elapsed().as_secs_f64();
    let point = BTreeMap::from([(rho, Atom::num(Rational::from((5, 2))))]);
    let range =
        prepared.required_master_range(&point, transport_cache::EpsilonRange::new(0, 0)?)?;
    let mut cache = RustFlowCache::default();
    let boundary = prepared.seed_cache(
        &mut cache,
        &point,
        range.last,
        number(5, 8)?,
        &options,
        &backend,
        &context,
    )?;
    let target = prepared.project_targets(
        &boundary,
        transport_cache::EpsilonRange::new(0, 0)?,
        options.digits,
    )?;
    let value = -target[0].coefficients[&0].re.clone();
    let seed_seconds = started.elapsed().as_secs_f64();
    let mut transported_values = Vec::new();
    if args.get(9).is_some_and(|arg| arg == "transport") {
        let policy = ScaledDistance {
            scales: BTreeMap::new(),
            admissible: |_: &transport_cache::CachedBoundary, _: &transport_cache::CachedPoint| {
                Ok(true)
            },
        };
        for destination in [Rational::from(2), Rational::from((5, 2)), Rational::from(3)] {
            let point = BTreeMap::from([(rho, Atom::num(destination.clone()))]);
            let mut destination_cache = RustFlowCache::default();
            destination_cache.insert(boundary.clone())?;
            let transported = prepared.flow().evaluate_to(
                &mut destination_cache,
                &point,
                range,
                &options,
                &context,
                &policy,
            )?;
            let target = prepared.project_targets(
                &transported.boundary,
                transport_cache::EpsilonRange::new(0, 0)?,
                options.digits,
            )?;
            transported_values.push(serde_json::json!({"rho": destination.to_string(), "value": (-target[0].coefficients[&0].re.clone()).as_raw().to_string(), "absolute_error": target[0].absolute_errors[&0].as_raw().to_string()}));
        }
    }
    println!(
        "{}",
        serde_json::json!({
            "seconds": started.elapsed().as_secs_f64(), "seed_seconds": seed_seconds, "transported": transported_values, "preparation_seconds": preparation_seconds,
            "mass_mode": options.mass_mode.name(), "digits": options.digits,
        "seed_digits": boundary.accuracy.verified_digits(), "masters": prepared.basis().iter().map(|integral| &integral.0).collect::<Vec<_>>(),
            "leading": range.leading, "last": range.last,
            "value": value.as_raw().to_string(), "absolute_error": target[0].absolute_errors[&0].as_raw().to_string(),
            "workers": options.workers,
        })
    );
    Ok(())
}
