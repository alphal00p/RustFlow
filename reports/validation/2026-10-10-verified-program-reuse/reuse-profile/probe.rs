#![allow(dead_code)]
//! Diagnostic only: reconstruct the complete historical source corpus embedded
//! in a locally generated program, then ask native RustRed to replay and search.
use bincode::Decode;
use rustred::persistence::{
    BinaryIoLimits, CoefficientId, DecodedCoefficientTable, SectionTag, inspect_program,
};
use rustred::solver::guarded::{
    GuardedProgram, GuardedSource, GuardedSourceSystem, IndexBounds, IndexDomain, IndexRole,
};
use rustred::solver::{Integral, Power, SearchOptions, Term};
use std::sync::Arc;
type Domain = Vec<(Option<i64>, Option<i64>)>;
type Label = Vec<(bool, i16)>;
#[derive(Decode)]
struct Tr {
    integral: Label,
    coefficient: usize,
}
#[derive(Decode)]
struct Sr {
    id: String,
    domain: Domain,
    conditions: Vec<usize>,
    terms: Vec<Tr>,
}
#[derive(Decode)]
struct Sd {
    row: usize,
    integral: Label,
    shifts: Vec<i16>,
}
#[derive(Decode)]
struct Rr {
    fixed: Vec<Option<i16>>,
    target: Label,
    rhs: Vec<Tr>,
    sources: Vec<Sd>,
    domain: Domain,
    discovery_domain: Domain,
    conditions: Vec<usize>,
    sector: Vec<bool>,
    permutation: Option<Vec<usize>>,
}
#[derive(Decode)]
struct Rec {
    schema: String,
    measure: String,
    roles: Vec<u8>,
    indices: Vec<usize>,
    sources: Vec<Sr>,
    zero_domains: Vec<Domain>,
    rules: Vec<Rr>,
    terminals: Vec<Vec<i64>>,
}
fn domain<const N: usize>(v: &Domain) -> IndexDomain<N> {
    IndexDomain::new(
        v.iter()
            .map(|&(lo, hi)| IndexBounds::new(lo, hi).unwrap())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap(),
    )
    .unwrap()
}
fn label<const N: usize>(v: &Label) -> Integral<N> {
    Integral::new(
        v.iter()
            .map(|&(s, p)| Power::new(s, p).unwrap())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap(),
    )
}
fn bounds<const N: usize>(v: &IndexDomain<N>) -> Domain {
    v.bounds().iter().map(|b| (b.lower(), b.upper())).collect()
}

use std::time::Instant;
use serde_json::json;
fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(args.len(), 3);
    let out = std::path::PathBuf::from(&args[2]);
    std::fs::create_dir_all(&out).unwrap();
    let mut timings = Vec::new();
    macro_rules! timed { ($name:expr, $body:expr) => {{
        std::fs::write(out.join("progress.json"), serde_json::to_vec_pretty(&json!({"phase":$name,"completed_timings":timings})).unwrap()).unwrap();
        println!("begin {}", $name);
        let t = Instant::now(); let value = $body;
        let seconds = t.elapsed().as_secs_f64();
        timings.push(json!({"phase":$name,"wall_seconds":seconds}));
        println!("end {} {:.9}s", $name, seconds); value
    }}; }
    let bytes = std::fs::read(&args[1]).unwrap();
    let limits = BinaryIoLimits::default();
    let envelope = inspect_program(&bytes, limits).unwrap();
    let (record, used): (Rec, usize) = bincode::decode_from_slice(
        envelope.section(SectionTag::PROGRAM).unwrap(), bincode::config::standard()).unwrap();
    assert_eq!(used, envelope.section(SectionTag::PROGRAM).unwrap().len());
    assert_eq!(record.schema, "rustred.guarded-source-program.v2");
    let table = DecodedCoefficientTable::import_generated_normalized(
        envelope.section(SectionTag::SYMBOLICA_STATE).unwrap(),
        envelope.section(SectionTag::COEFFICIENTS).unwrap(), limits).unwrap();
    let polynomial = |id| { let c = table.coefficient(CoefficientId::try_from_index(id).unwrap()).unwrap().clone(); assert!(c.denominator.is_one()); c.numerator };
    let roles: [IndexRole;16] = record.roles.iter().map(|r| match r {0=>IndexRole::Ordinary,1=>IndexRole::RequiredCut,2=>IndexRole::Occupation,_=>panic!()}).collect::<Vec<_>>().try_into().unwrap();
    let rows = record.sources.iter().map(|s| GuardedSource::new(s.id.clone(),
        s.terms.iter().map(|t| Term { integral:label(&t.integral), coefficient:polynomial(t.coefficient) }).collect(),
        domain(&s.domain)).with_nonzero_conditions(s.conditions.iter().map(|&id| polynomial(id)).collect())).collect();
    let sources = Arc::new(GuardedSourceSystem::new(record.measure.clone(), roles,
        record.indices.clone().try_into().unwrap(), rows).unwrap()
        .with_zero_domains(record.zero_domains.iter().map(domain).collect()).unwrap());
    let program = timed!("decode_generated_and_replay_all_rules", GuardedProgram::decode_generated(&bytes, sources.clone(), limits).unwrap());
    assert!(program.terminals().is_empty());
    assert_eq!(program.rules().len(), record.rules.len());
    let count = program.rules().len();
    let before = timed!("encode_loaded_program", program.encode_native(limits).unwrap());
    let empty = GuardedProgram::new(sources.clone(), vec![], []).unwrap();
    let program = timed!("union_verified_with_empty_same_source", program.union_verified(empty, [], count).unwrap());
    let verified = timed!("encode_verified_union", program.encode_native(limits).unwrap());
    assert_eq!(before, verified, "verified union changed persisted bytes");
    let empty = GuardedProgram::new(sources.clone(), vec![], []).unwrap();
    let program = timed!("union_replayed_with_empty_same_source", program.union_replayed(empty, [], count).unwrap());
    let replayed = timed!("encode_replayed_union", program.encode_native(limits).unwrap());
    assert_eq!(verified, replayed, "replayed and verified union differ");
    assert_eq!(program.rules().len(), count);
    assert!(program.terminals().is_empty());
    std::fs::write(out.join("unchanged-program.bin"), &verified).unwrap();
    let result = json!({
        "scope":"One completed immutable native program per arm. Initial decode independently replays all proofs, followed by verified and replayed empty union on the same unchanged program. No discovery, new closure, period or end-to-end speed claim.",
        "source_schema":record.schema,"source_identity":record.measure,"source_count":record.sources.len(),
        "rule_count":count,"original_input_bytes":bytes.len(),"encoded_bytes":verified.len(),
        "initial_loaded_proof_replay":true,"verified_union_equals_loaded_encoded_bytes":true,
        "verified_union_equals_replayed_encoded_bytes":true,
        "third_decode_omitted":"The verified bytes exactly equal the immediately replayed loaded program; the subsequent replayed union independently revalidates the unchanged proofs. Small unit controls separately decode verified unions.",
        "timings":timings,"timing_scope":"Single samples in one process on a shared host. Empty-program construction and encoding are outside each union timer. Library/harness compilation excluded. Tiny verified duration is subject to clock noise; no system-wide speedup is inferred."});
    std::fs::write(out.join("result.json"), serde_json::to_vec_pretty(&result).unwrap()).unwrap();
    std::fs::write(out.join("progress.json"), serde_json::to_vec_pretty(&json!({"phase":"complete","completed_timings":timings})).unwrap()).unwrap();
    println!("PASS: {count} rules, byte-identical verified/replayed unions");
}
