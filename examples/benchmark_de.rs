//! Shared-input differential-equation benchmark; see scripts/benchmark_de.py.
use serde::Deserialize;
use std::{collections::BTreeMap, time::Instant};
use symbolica::prelude::*;
use symbolica_amflow::*;

#[derive(Deserialize)]
struct Input {
    matrix: Vec<Vec<String>>,
    epsilon: String,
    start: [String; 2],
    end: [String; 2],
    boundary: Vec<[String; 2]>,
    digits: u32,
    working_decimal_digits: u32,
    working_bits: u32,
    order: usize,
}

fn atom(text: &str) -> Result<Atom> {
    Atom::parse(text, "performance", Default::default())
        .map_err(|error| Error::InvalidInput(error.to_string()))
}

fn main() -> Result<()> {
    let path = std::env::args()
        .nth(1)
        .ok_or_else(|| Error::InvalidInput("usage: benchmark_de INPUT.json".into()))?;
    let input: Input = serde_json::from_reader(std::fs::File::open(path)?)
        .map_err(|error| Error::InvalidInput(error.to_string()))?;
    if input.working_decimal_digits < input.digits || input.working_bits < 64 {
        return Err(Error::InvalidInput("invalid benchmark precision".into()));
    }
    let started = Instant::now();
    // Match the original C++ solver's Boost decimal-to-binary conversion;
    // unlike Precision::decimal, this adds no extra Rust-only guard bits.
    let p = Precision {
        bits: input.working_bits,
    };
    let substitutions = BTreeMap::from([(atom("eps")?, atom(&input.epsilon)?)]);
    let system = DifferentialSystem {
        variable: symbol!("performance::eta"),
        matrix: input
            .matrix
            .iter()
            .map(|row| {
                row.iter()
                    .map(|entry| {
                        Ok(family::substitute(&atom(entry)?, &substitutions)
                            .together()
                            .cancel())
                    })
                    .collect::<Result<_>>()
            })
            .collect::<Result<_>>()?,
    };
    let compiled = system.compile(p, &Default::default())?;
    let boundary = BoundaryData {
        point: p.parse(&input.start[0], &input.start[1])?,
        values: input
            .boundary
            .iter()
            .map(|value| p.parse(&value[0], &value[1]))
            .collect::<Result<_>>()?,
    };
    let prepared_ns = started.elapsed().as_nanos();
    let running = Instant::now();
    let result = compiled.transport(
        &boundary,
        &[p.parse(&input.end[0], &input.end[1])?],
        &FlowOptions {
            digits: input.digits,
            guard_digits: input.working_decimal_digits - input.digits,
            series_order: input.order,
            ..Default::default()
        },
        &RunContext::default(),
    )?;
    let transport_ns = running.elapsed().as_nanos();
    let values: Vec<_> = result
        .values
        .iter()
        .map(|value| [value.re.as_raw().to_string(), value.im.as_raw().to_string()])
        .collect();
    println!(
        "{}",
        serde_json::json!({
            "values": values,
            "working_bits": p.bits,
            "prepare_ns": prepared_ns,
            "transport_ns": transport_ns,
            "steps": result.diagnostics.steps,
            "rejected_steps": result.diagnostics.rejected_steps,
        })
    );
    Ok(())
}
