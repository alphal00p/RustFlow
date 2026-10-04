//! Export the integrity-checked paper-family matrix at epsilon=1/2700 for a benchmark.
//! Link with the immutable release dependencies; see docs/performance.md.
use serde::Deserialize;
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};
use symbolica::{
    poly::{PolyVariable, polynomial::MultivariatePolynomial},
    prelude::*,
};
use symbolica_amflow::{DifferentialSystem, Precision, family};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
#[derive(Deserialize)]
struct Entry {
    version: u32,
    key: String,
    digest: String,
    payload: Vec<u8>,
}
#[derive(Deserialize)]
struct Stored {
    basis: Vec<Vec<i16>>,
    matrix: Vec<Vec<Vec<u8>>>,
}
fn polynomial(value: &MultivariatePolynomial<IntegerRing, u16>) -> String {
    if value.is_zero() {
        return "0".into();
    }
    let mut out = String::new();
    for (coefficient, powers) in value.coefficients.iter().zip(value.exponents_iter()) {
        let text = coefficient.to_string();
        if !out.is_empty() && !text.starts_with('-') {
            out.push('+');
        }
        out.push_str(&text);
        if powers[0] > 0 {
            out.push_str("*eta");
            if powers[0] > 1 {
                out.push('^');
                out.push_str(&powers[0].to_string());
            }
        }
    }
    out
}
fn main() -> Result<()> {
    let arguments = std::env::args().collect::<Vec<_>>();
    if arguments.len() != 3 {
        return Err("usage: export_cached_de INPUT_SYSTEM.json OUTPUT_CASE.json".into());
    }
    let raw = std::fs::read(&arguments[1])?;
    let entry: Entry = serde_json::from_slice(&raw)?;
    if entry.version != 1 || blake3::hash(&entry.payload).to_hex().as_str() != entry.digest {
        return Err("incompatible or corrupt cache envelope".into());
    }
    let stored: Stored = serde_json::from_slice(&entry.payload)?;
    let n = stored.basis.len();
    if n == 0 || stored.matrix.len() != n || stored.matrix.iter().any(|row| row.len() != n) {
        return Err("cached matrix dimensions".into());
    }
    let mut symbols = BTreeSet::new();
    let mut matrix = Vec::new();
    for row in stored.matrix {
        let mut decoded = Vec::new();
        for bytes in row {
            let mut source = bytes.as_slice();
            let atom = Atom::import(&mut source, None)?;
            if !source.is_empty() {
                return Err("trailing serialized Atom data".into());
            }
            symbols.extend(atom.as_view().get_all_symbols(true));
            decoded.push(atom);
        }
        matrix.push(decoded);
    }
    let old = symbol!("symbolica_amflow::eta");
    if symbols != BTreeSet::from([old]) {
        return Err(format!(
            "expected eta-only matrix, found {:?}",
            symbols.iter().map(|s| s.get_name()).collect::<Vec<_>>()
        )
        .into());
    }
    let eta = symbol!("performance::eta");
    let substitution = BTreeMap::from([(Atom::var(old), Atom::var(eta))]);
    let variables = Arc::new(vec![PolyVariable::Symbol(eta)]);
    let mut exported = Vec::new();
    let mut exact = Vec::new();
    let mut nonzero = 0;
    for row in matrix {
        let mut strings = Vec::new();
        let mut atoms = Vec::new();
        for atom in row {
            let renamed = family::substitute(&atom, &substitution);
            let rational: RationalPolynomial<IntegerRing, u16> =
                renamed.try_to_rational_polynomial(&Q, &Z, Some(variables.clone()))?;
            let text = if rational.is_zero() {
                "0".into()
            } else {
                nonzero += 1;
                format!(
                    "({})/({})",
                    polynomial(&rational.numerator),
                    polynomial(&rational.denominator)
                )
            };
            let parsed = Atom::parse(&text, "performance", Default::default())?;
            let difference: RationalPolynomial<IntegerRing, u16> =
                (&parsed - &renamed).try_to_rational_polynomial(&Q, &Z, Some(variables.clone()))?;
            if !difference.is_zero() {
                return Err("export roundtrip changed a matrix entry".into());
            }
            atoms.push(renamed);
            strings.push(text);
        }
        exact.push(atoms);
        exported.push(strings);
    }
    let system = DifferentialSystem {
        variable: eta,
        matrix: exact,
    };
    let mut pole_sets = Vec::new();
    let mut radius = 4_i64;
    for bits in [201, 267] {
        let p = Precision { bits };
        let compiled = system.compile(p, &Default::default())?;
        for pole in &compiled.poles {
            if !p.finite(pole) {
                return Err("nonfinite pole".into());
            }
        }
        while compiled
            .poles
            .iter()
            .any(|pole| p.norm(pole) >= p.real(radius / 4))
        {
            radius = radius
                .checked_mul(2)
                .ok_or("pole radius exceeds exact endpoint range")?;
        }
        pole_sets.push(serde_json::json!({"working_bits":bits,"poles":compiled.poles.iter().map(|pole|[pole.re.as_raw().to_string(),pole.im.as_raw().to_string()]).collect::<Vec<_>>()}));
    }
    // All numerically found poles lie inside |eta| < R/4. The segment from
    // -2 R i to -R i therefore stays at least 3 R/4 away from those poles.
    let start = radius.checked_mul(2).ok_or("start endpoint overflow")?;
    let case = serde_json::json!({"name":"paper_27_integral_system","matrix":exported,"epsilon":"1/2700","start":["0".to_owned(),format!("-{start}")],"end":["0".to_owned(),format!("-{radius}")],"boundary":(0..n).map(|i|[(i+1).to_string(),((i%5) as i64-2).to_string()]).collect::<Vec<_>>(),"provenance":{"scope":"Actual sampled paper-family differential matrix; deterministic supplied boundary vector is not a physical integral boundary.","cache_file":arguments[1],"cache_key":entry.key,"cache_payload_blake3_verified":entry.digest,"basis":stored.basis,"nonzero_entries":nonzero,"exact_matrix_roundtrip_verified_entries":n*n,"boundary_rule":"component j=0..n-1 is (j+1)+((j mod5)-2)i, exact integers","pole_checks":pole_sets,"endpoint_radius":radius,"pole_clearance_lower_bound_from_numerically_found_poles":format!("{}/4",3_i128*i128::from(radius)),"endpoint_selection":"Smallest power of two R>=4 enclosing both computed pole sets in |eta|<R/4; straight segment -2Ri to -Ri.","blocks":system.blocks()?}});
    std::fs::write(&arguments[2], serde_json::to_vec_pretty(&case)?)?;
    println!(
        "{}",
        serde_json::json!({"dimension":n,"nonzero_entries":nonzero,"radius":radius,"output":arguments[2]})
    );
    Ok(())
}
