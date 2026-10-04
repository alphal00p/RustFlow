//! Shared-input singular matching benchmark; see scripts/benchmark_singular.py.
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
    Atom::parse(text, "singular_performance", Default::default())
        .map_err(|error| Error::InvalidInput(error.to_string()))
}

fn main() -> Result<()> {
    let path = std::env::args()
        .nth(1)
        .ok_or_else(|| Error::InvalidInput("usage: benchmark_singular INPUT.json".into()))?;
    let input: Input = serde_json::from_reader(std::fs::File::open(path)?)
        .map_err(|error| Error::InvalidInput(error.to_string()))?;
    if input.digits == 0
        || input.working_decimal_digits < input.digits
        || input.working_bits < 64
        || input.order < 2
    {
        return Err(Error::InvalidInput(
            "invalid benchmark precision/order".into(),
        ));
    }
    let started = Instant::now();
    let p = Precision {
        bits: input.working_bits,
    };
    let epsilon = symbol!("singular_performance::eps");
    let substitutions = BTreeMap::from([(Atom::var(epsilon), atom(&input.epsilon)?)]);
    let system = DifferentialSystem {
        variable: symbol!("singular_performance::eta"),
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
    let start = p.parse(&input.start[0], &input.start[1])?;
    let end = p.parse(&input.end[0], &input.end[1])?;
    let boundary = input
        .boundary
        .iter()
        .map(|v| p.parse(&v[0], &v[1]))
        .collect::<Result<Vec<_>>>()?;
    let prepare_ns = started.elapsed().as_nanos();
    let running = Instant::now();
    let basis = system.frobenius(p, &Default::default(), input.order)?;
    let frobenius_ns = running.elapsed().as_nanos();
    let matching = Instant::now();
    let constants = basis.match_values(&start, &boundary, &Default::default())?;
    let matching_ns = matching.elapsed().as_nanos();
    let singular_matching_ns = running.elapsed().as_nanos();
    let evaluation = Instant::now();
    let matrix = basis.evaluate(&end, &Default::default())?;
    let values = matrix
        .iter()
        .map(|row| {
            row.iter()
                .zip(&constants)
                .fold(p.zero(), |sum, (a, b)| p.add(&sum, &p.mul(a, b)))
        })
        .collect::<Vec<_>>();
    let endpoint = basis.physical_limit(&constants, epsilon)?;
    let evaluation_ns = evaluation.elapsed().as_nanos();
    let strings = |values: &[ComplexFloat]| {
        values
            .iter()
            .map(|v| [v.re.as_raw().to_string(), v.im.as_raw().to_string()])
            .collect::<Vec<_>>()
    };
    println!(
        "{}",
        serde_json::json!({
            "values": strings(&values),
            "endpoint": strings(&endpoint),
            "working_bits": p.bits,
            "prepare_ns": prepare_ns,
            "frobenius_ns": frobenius_ns,
            "matching_ns": matching_ns,
            "singular_matching_ns": singular_matching_ns,
            "evaluation_ns": evaluation_ns,
            "exponents": basis.columns.iter().map(|c| c.exponent.to_string()).collect::<Vec<_>>(),
            "max_log_power": basis.columns.iter().flat_map(|c| &c.coefficients)
                .map(|coefficient| coefficient.len().saturating_sub(1)).max().unwrap_or(0),
        })
    );
    Ok(())
}
