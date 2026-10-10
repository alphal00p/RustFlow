//! Isolated finite-data candidate fitter. No physical identity is certified.
use serde_json::{Value, json};
use std::{fs, str::FromStr};
use symbolica::{
    domains::{Ring, RingOps},
    prelude::{Integer, Q, Rational},
    tensors::matrix::Matrix,
};
const TRAIN: usize = 48;
const HOLDOUT: usize = 16;
const MAX_TERMS: usize = 200_000;
const MAX_BYTES: usize = 128 * 1024 * 1024;
#[derive(Clone, Debug)]
struct Channel {
    id: String,
    beta: Rational,
    c: Vec<Rational>,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Theta,
    Recurrence,
}
impl Kind {
    fn name(self) -> &'static str {
        match self {
            Self::Theta => "theta",
            Self::Recurrence => "recurrence",
        }
    }
    fn parse(s: &str) -> Result<Self, String> {
        match s {
            "theta" => Ok(Self::Theta),
            "recurrence" => Ok(Self::Recurrence),
            _ => Err("unknown fit kind".into()),
        }
    }
}
#[derive(Clone, Debug)]
struct Candidate {
    kind: Kind,
    order: usize,
    degree: usize,
    a: Vec<Rational>,
}
fn rat(s: &str) -> Result<Rational, String> {
    if s.len() > MAX_BYTES {
        return Err("rational exceeds retained input bound".into());
    }
    let mut parts = s.split('/');
    let n = Integer::from_str(parts.next().ok_or("empty rational")?).map_err(|e| e.to_string())?;
    let d = Integer::from_str(parts.next().unwrap_or("1")).map_err(|e| e.to_string())?;
    if parts.next().is_some() || d == Integer::from(0) {
        return Err("invalid rational denominator".into());
    }
    Ok(Q.to_element(n, d, true))
}
fn channels(v: &Value, expected: usize) -> Result<Vec<Channel>, String> {
    if v["schema"] != "integrated-rational-prefix.v1"
        || v["dimension"] != "13/2"
        || (v["q"] != 2 && v["q"] != 1)
    {
        return Err("unsupported sequence schema/dimension/grading".into());
    }
    if v["series_identity"].as_str().is_none_or(str::is_empty) {
        return Err("missing nonempty source-bound series identity".into());
    }
    if v["first_index"].as_u64() != Some(if expected == TRAIN { 0 } else { TRAIN as u64 }) {
        return Err("incorrect coefficient start index".into());
    }
    if v["q"] == 1
        && v["grading_proof_identity"]
            .as_str()
            .is_none_or(str::is_empty)
    {
        return Err("q=1 requires source-bound parity/offset proof identity".into());
    }
    let ch = v["channels"].as_array().ok_or("missing channels")?;
    if ch.is_empty() || ch.len() > 32 {
        return Err("channel count outside cap".into());
    }
    let mut out = Vec::new();
    for x in ch {
        let id = x["id"].as_str().ok_or("missing channel id")?.to_owned();
        if out.iter().any(|y: &Channel| y.id == id) {
            return Err("duplicate channel id".into());
        }
        if x["log_power"].as_u64() != Some(0) {
            return Err("first fitter supports log-free channels only; no logs dropped".into());
        }
        let c = x["coefficients"]
            .as_array()
            .ok_or("missing exact coefficients")?;
        if c.len() != expected {
            return Err(format!(
                "expected exactly {expected} coefficients; unknown tail is not zero"
            ));
        }
        out.push(Channel {
            id,
            beta: rat(x["beta"].as_str().ok_or("missing t exponent")?)?,
            c: c.iter()
                .map(|v| rat(v.as_str().ok_or("coefficients must be rational strings")?))
                .collect::<Result<_, _>>()?,
        });
    }
    Ok(out)
}
fn check_values<'a>(values: impl Iterator<Item = &'a Rational>) -> Result<(), String> {
    let mut count = 0usize;
    let mut bytes = 0usize;
    for x in values {
        count = count
            .checked_add(usize::from(!Q.is_zero(x)))
            .ok_or("monomial counter overflow")?;
        bytes = bytes
            .checked_add(x.to_string().len())
            .ok_or("payload counter overflow")?;
        if count > MAX_TERMS || bytes > MAX_BYTES {
            return Err("retained canonical algebra payload cap".into());
        }
    }
    Ok(())
}
fn input_values(ch: &[Channel]) -> impl Iterator<Item = &Rational> {
    ch.iter()
        .flat_map(|ch| ch.c.iter().chain(std::iter::once(&ch.beta)))
}

fn row(ch: &Channel, n: usize, kind: Kind, r: usize, d: usize) -> Vec<Rational> {
    let mut v = Vec::with_capacity((r + 1) * (d + 1));
    match kind {
        Kind::Theta => {
            for s in 0..=d {
                for j in 0..=r {
                    let a = if n < s {
                        Q.zero()
                    } else {
                        let lambda = Q.add(&Rational::from((n - s) as i64), &ch.beta);
                        Q.mul(&lambda.pow(j as u64), &ch.c[n - s])
                    };
                    v.push(a);
                }
            }
        }
        Kind::Recurrence => {
            for j in 0..=r {
                for k in 0..=d {
                    v.push(Q.mul(&Rational::from(n as i64).pow(k as u64), &ch.c[n + j]));
                }
            }
        }
    }
    v
}
fn fit(ch: &[Channel], kind: Kind) -> Result<(Option<Candidate>, Value), String> {
    if ch.iter().any(|x| x.c.len() != TRAIN) {
        return Err("main training protocol requires 48 known terms".into());
    }
    if ch.iter().all(|x| x.c.iter().all(|v| Q.is_zero(v))) {
        return Ok((None, json!({"status":"zero_prefix_not_a_zero_certificate"})));
    }
    check_values(input_values(ch))?;
    let mut shapes = (1..=4)
        .flat_map(|r| (0..=4).map(move |d| ((r + 1) * (d + 1), r, d)))
        .collect::<Vec<_>>();
    shapes.sort();
    let mut attempts = Vec::new();
    for (u, r, d) in shapes {
        if u > 25 || attempts.len() >= 24 {
            break;
        }
        let end = if kind == Kind::Theta {
            TRAIN
        } else {
            TRAIN - r
        };
        let mut rows = ch
            .iter()
            .flat_map(|c| (0..end).map(move |n| row(c, n, kind, r, d)))
            .collect::<Vec<_>>();
        let total_rows = rows.len();
        rows.retain(|x| x.iter().any(|v| !Q.is_zero(v)));
        if total_rows > 2048 || rows.len() < u - 1 + 8 {
            attempts.push(json!({"order":r,"degree":d,"unknowns":u,"nonzero_rows":rows.len(),"status":"insufficient_excess_rows_or_row_cap"}));
            continue;
        }
        check_values(input_values(ch).chain(rows.iter().flatten()))?;
        // Native matrix takes ownership of rows; no retained training-matrix clone.
        let mut matrix = Matrix::from_nested_vec(rows, Q)?;
        let rank = matrix.row_reduce(u as u32);
        let rr = matrix.into_vec();
        check_values(input_values(ch).chain(rr.iter()))?;
        let pivots = rr
            .chunks(u)
            .filter_map(|r| r.iter().position(|v| !Q.is_zero(v)))
            .collect::<Vec<_>>();
        let free = (0..u).filter(|j| !pivots.contains(j)).collect::<Vec<_>>();
        attempts.push(json!({"order":r,"degree":d,"unknowns":u,"rows":total_rows,"rank":rank,"nullity":free.len()}));
        if free.is_empty() {
            continue;
        }
        if free.len() != 1 {
            return Ok((
                None,
                json!({"status":"training_ambiguity","attempts":attempts}),
            ));
        }
        let mut a = vec![Q.zero(); u];
        a[free[0]] = Q.one();
        for (i, p) in pivots.iter().enumerate() {
            a[*p] = Q.neg(&rr[i * u + free[0]]);
        }
        let normal = a
            .iter()
            .find(|a| !Q.is_zero(a))
            .ok_or("zero candidate")?
            .clone();
        for x in &mut a {
            *x = &*x / &normal;
        }
        check_values(input_values(ch).chain(rr.iter()).chain(a.iter()))?;
        // Reconstruct one verification equation at a time from training, never holdout.
        for channel in ch {
            for n in 0..end {
                let equation = row(channel, n, kind, r, d);
                if !Q.is_zero(&dot(&equation, &a)) {
                    return Err("native nullspace verification failed".into());
                }
            }
        }
        return Ok((
            Some(Candidate {
                kind,
                order: r,
                degree: d,
                a,
            }),
            json!({"status":"training_candidate_only","attempts":attempts}),
        ));
    }
    Ok((
        None,
        json!({"status":"no_candidate_within_fixed_shapes","attempts":attempts}),
    ))
}
fn dot(a: &[Rational], b: &[Rational]) -> Rational {
    a.iter()
        .zip(b)
        .fold(Q.zero(), |s, (a, b)| Q.add(&s, &Q.mul(a, b)))
}
fn candidate_json(c: &Candidate) -> Value {
    json!({"kind":c.kind.name(),"order":c.order,"degree":c.degree,"coefficients":c.a.iter().map(ToString::to_string).collect::<Vec<_>>(),"coefficient_order":if c.kind==Kind::Theta{"t_degree_major_theta_order_minor"}else{"shift_major_n_degree_minor"},"status":"finite_training_candidate_not_identity"})
}
fn candidate(v: &Value) -> Result<Candidate, String> {
    let kind = Kind::parse(v["kind"].as_str().ok_or("missing kind")?)?;
    let r = v["order"].as_u64().ok_or("missing order")? as usize;
    let d = v["degree"].as_u64().ok_or("missing degree")? as usize;
    if r == 0 || r > 4 || d > 4 {
        return Err("candidate shape exceeds frozen caps".into());
    }
    let a = v["coefficients"]
        .as_array()
        .ok_or("missing operator coefficients")?
        .iter()
        .map(|x| rat(x.as_str().ok_or("non-rational operator coefficient")?))
        .collect::<Result<Vec<_>, _>>()?;
    if a.len() != (r + 1) * (d + 1) || a.iter().all(|x| Q.is_zero(x)) {
        return Err("invalid operator size/zero operator".into());
    }
    Ok(Candidate {
        kind,
        order: r,
        degree: d,
        a,
    })
}
fn active_weights(
    c: &Candidate,
    beta: &Rational,
    n: usize,
) -> Result<std::collections::BTreeMap<usize, Rational>, String> {
    let mut weights = std::collections::BTreeMap::new();
    for (at, a) in c.a.iter().enumerate() {
        if Q.is_zero(a) {
            continue;
        }
        let (index, power) = match c.kind {
            Kind::Theta => {
                let s = at / (c.order + 1);
                let j = at % (c.order + 1);
                if n < s {
                    continue;
                }
                (
                    n - s,
                    Q.add(&Rational::from((n - s) as i64), beta).pow(j as u64),
                )
            }
            Kind::Recurrence => {
                let shift = at / (c.degree + 1);
                let degree = at % (c.degree + 1);
                (n + shift, Rational::from(n as i64).pow(degree as u64))
            }
        };
        let value = Q.mul(a, &power);
        let previous = weights.remove(&index).unwrap_or_else(|| Q.zero());
        let sum = Q.add(&previous, &value);
        if !Q.is_zero(&sum) {
            weights.insert(index, sum);
        }
    }
    Ok(weights)
}
fn validate_holdout(train: &[Channel], tail: &[Channel], c: &Candidate) -> Result<Value, String> {
    if train.len() != tail.len() {
        return Err("channel set changed at holdout".into());
    }
    let mut residuals = Vec::new();
    let mut missing = Vec::new();
    let mut validated = 0usize;
    let support =
        c.a.iter()
            .enumerate()
            .filter(|(_, a)| !Q.is_zero(a))
            .map(|(i, _)| match c.kind {
                Kind::Theta => i / (c.order + 1),
                Kind::Recurrence => i / (c.degree + 1),
            })
            .collect::<Vec<_>>();
    let start = match c.kind {
        Kind::Theta => TRAIN + support.iter().min().ok_or("zero operator")?,
        Kind::Recurrence => TRAIN - support.iter().max().ok_or("zero operator")?,
    };
    for (a, b) in train.iter().zip(tail) {
        if a.id != b.id || a.beta != b.beta || a.c.len() != TRAIN || b.c.len() != HOLDOUT {
            return Err("holdout channel identity/length mismatch".into());
        }
        let mut both = a.c.clone();
        both.extend(b.c.iter().cloned());
        let mut touched = std::collections::BTreeSet::new();
        for n in start..start + HOLDOUT {
            let weights = active_weights(c, &a.beta, n)?;
            if weights.keys().any(|i| *i >= both.len()) {
                return Err("operator needs unknown coefficient".into());
            }
            let held = weights
                .keys()
                .filter(|i| **i >= TRAIN)
                .copied()
                .collect::<Vec<_>>();
            let value = weights
                .iter()
                .fold(Q.zero(), |s, (i, w)| Q.add(&s, &Q.mul(w, &both[*i])));
            let genuine = !held.is_empty();
            if genuine {
                validated += 1;
                touched.extend(held.iter().copied());
            }
            residuals.push(json!({"channel":a.id,"row":n,"heldout_indices_with_nonzero_weight":held,"genuine_heldout_row":genuine,"residual":value.to_string(),"zero":Q.is_zero(&value)}));
        }
        for i in TRAIN..TRAIN + HOLDOUT {
            if !touched.contains(&i) {
                missing.push(json!({"channel":a.id,"coefficient_index":i}));
            }
        }
    }
    let status = if residuals.iter().any(|x| x["zero"] != true) {
        "heldout_rejected_protocol_stops"
    } else if validated != HOLDOUT * train.len() || !missing.is_empty() {
        "heldout_inconclusive_operator_singularity"
    } else {
        "heldout_pass_candidate_only"
    };
    Ok(
        json!({"status":status,"rows":residuals,"genuine_heldout_rows":validated,"unconstrained_heldout_coefficients":missing,"identity_certified":false}),
    )
}

fn read(path: &str) -> Result<(Vec<u8>, Value), String> {
    if fs::metadata(path).map_err(|e| e.to_string())?.len() > MAX_BYTES as u64 {
        return Err("input byte cap before read".into());
    }
    let b = fs::read(path).map_err(|e| e.to_string())?;
    if b.len() > MAX_BYTES {
        return Err("input byte cap".into());
    }
    let v = serde_json::from_slice(&b).map_err(|e| e.to_string())?;
    Ok((b, v))
}
fn write_new(path: &str, v: &Value) -> Result<(), String> {
    use std::io::Write;
    let mut f = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|e| e.to_string())?;
    f.write_all(serde_json::to_string_pretty(v).unwrap().as_bytes())
        .map_err(|e| e.to_string())
}
fn run(args: &[String]) -> Result<(), String> {
    match args.get(1).map(String::as_str){
        Some("digest") if args.len()==3=>{let (bytes,_)=read(&args[2])?;println!("{}",blake3::hash(&bytes).to_hex());Ok(())},
        Some("fit") if args.len()==5=>{let (bytes,v)=read(&args[2])?;let ch=channels(&v,TRAIN)?;let kind=Kind::parse(&args[4])?;let(c,audit)=fit(&ch,kind)?;write_new(&args[3],&json!({"schema":"integrated-fit-candidate.v1","training_blake3":blake3::hash(&bytes).to_hex().to_string(),"series_identity":v["series_identity"],"candidate":c.as_ref().map(candidate_json),"audit":audit,"holdout_read":false,"identity_certified":false}))},
        Some("validate") if args.len()==7=>{let (cb,cv)=read(&args[2])?;if blake3::hash(&cb).to_hex().as_str()!=args[3]{return Err("frozen candidate hash mismatch before holdout access".into());}let(tb,tv)=read(&args[4])?;if cv["training_blake3"]!=blake3::hash(&tb).to_hex().to_string(){return Err("training bytes changed".into());}let c=candidate(&cv["candidate"])?;let(_,hv)=read(&args[5])?;if tv["q"]!=hv["q"] || tv["dimension"]!=hv["dimension"] || tv["grading_proof_identity"]!=hv["grading_proof_identity"] || tv["series_identity"]!=hv["series_identity"] || cv["series_identity"]!=tv["series_identity"]{return Err("holdout series identity mismatch".into());}write_new(&args[6],&validate_holdout(&channels(&tv,TRAIN)?,&channels(&hv,HOLDOUT)?,&c)?)},
        _=>Err("usage: fit TRAIN CANDIDATE theta|recurrence; validate CANDIDATE FROZEN_BLAKE3 TRAIN HOLDOUT REPORT".into())
    }
}
fn main() {
    if let Err(e) = run(&std::env::args().collect::<Vec<_>>()) {
        eprintln!("{e}");
        std::process::exit(2);
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    fn geometric(a: i64, n: usize) -> Vec<Rational> {
        (0..n).map(|n| Rational::from(a).pow(n as u64)).collect()
    }
    fn channel(c: Vec<Rational>) -> Channel {
        Channel {
            id: "synthetic".into(),
            beta: Q.zero(),
            c,
        }
    }
    #[test]
    fn exact_native_theta_fit_and_untouched_holdout() {
        let full = geometric(2, 64);
        let tr = channel(full[..48].to_vec());
        let (c, _) = fit(&[tr.clone()], Kind::Theta).unwrap();
        let c = c.unwrap();
        assert_eq!((c.order, c.degree), (1, 1));
        assert_eq!(
            c.a,
            vec![Q.zero(), Q.one(), Rational::from(-2), Rational::from(-2)]
        );
        let h = validate_holdout(&[tr], &[channel(full[48..].to_vec())], &c).unwrap();
        assert_eq!(h["status"], "heldout_pass_candidate_only");
        assert_eq!(h["identity_certified"], false);
    }
    #[test]
    fn same_training_false_extension_is_rejected_without_refit() {
        let full = geometric(2, 64);
        let tr = channel(full[..48].to_vec());
        let (c, _) = fit(&[tr.clone()], Kind::Theta).unwrap();
        let c = c.unwrap();
        let before = candidate_json(&c);
        let mut tail = full[48..].to_vec();
        tail[0] = Q.add(&tail[0], &Q.one());
        assert_eq!(
            validate_holdout(&[tr], &[channel(tail)], &c).unwrap()["status"],
            "heldout_rejected_protocol_stops"
        );
        assert_eq!(candidate_json(&c), before);
    }
    #[test]
    fn recurrence_and_frobenius_exponent_are_explicit() {
        let full = geometric(3, 64);
        let tr = channel(full[..48].to_vec());
        let (c, _) = fit(&[tr.clone()], Kind::Recurrence).unwrap();
        let c = c.unwrap();
        assert_eq!((c.order, c.degree), (1, 0));
        assert_eq!(
            validate_holdout(&[tr], &[channel(full[48..].to_vec())], &c).unwrap()["status"],
            "heldout_pass_candidate_only"
        );
        let mut x = channel(geometric(2, 48));
        x.beta = rat("1/2").unwrap();
        let (c, _) = fit(&[x.clone()], Kind::Theta).unwrap();
        let c = c.unwrap();
        for n in 0..48 {
            assert!(Q.is_zero(&dot(&row(&x, n, c.kind, c.order, c.degree), &c.a)));
        }
    }
    #[test]
    fn zero_and_missing_prefixes_are_not_certificates() {
        assert!(
            fit(&[channel(vec![Q.zero(); 48])], Kind::Theta)
                .unwrap()
                .0
                .is_none()
        );
        assert!(fit(&[channel(vec![Q.one(); 24])], Kind::Theta).is_err());
        assert!(rat("1/0").is_err());
    }
    #[test]
    fn exact_factorial_sequence_uses_native_rational_matrix() {
        let mut all = vec![Q.one()];
        for n in 1..64 {
            all.push(&all[n - 1] / &Rational::from(n as i64));
        }
        let tr = channel(all[..48].to_vec());
        let (c, _) = fit(&[tr.clone()], Kind::Theta).unwrap();
        let c = c.unwrap();
        assert_eq!(c.a, vec![Q.zero(), Q.one(), Rational::from(-1), Q.zero()]);
        assert_eq!(
            validate_holdout(&[tr], &[channel(all[48..].to_vec())], &c).unwrap()["status"],
            "heldout_pass_candidate_only"
        );
    }
    #[test]
    fn formal_channels_are_fit_jointly() {
        let mut a = channel(geometric(2, 48));
        a.id = "hard-master-a".into();
        let mut b = channel(geometric(3, 48));
        b.id = "hard-master-b".into();
        let (op, _) = fit(&[a.clone(), b.clone()], Kind::Theta).unwrap();
        let op = op.expect("common operator within fixed cap");
        assert!(op.order > 1);
        let mut ah = channel(geometric(2, 64)[48..].to_vec());
        ah.id = a.id.clone();
        let mut bh = channel(geometric(3, 64)[48..].to_vec());
        bh.id = b.id.clone();
        assert_eq!(
            validate_holdout(&[a, b], &[ah, bh], &op).unwrap()["status"],
            "heldout_pass_candidate_only"
        );
    }
    #[test]
    fn schema_preserves_unknowns_logs_and_channel_identity() {
        let mut v = json!({"schema":"integrated-rational-prefix.v1","dimension":"13/2","q":2,"first_index":0,"series_identity":"synthetic-test-only","channels":[{"id":"x","beta":"0","log_power":0,"coefficients":vec!["1";48]}]});
        assert!(channels(&v, 48).is_ok());
        v["channels"][0]["log_power"] = json!(1);
        assert!(channels(&v, 48).is_err());
        v["channels"][0]["log_power"] = json!(0);
        v["first_index"] = json!(48);
        assert!(channels(&v, 48).is_err());
        v["first_index"] = json!(0);
        v["series_identity"] = json!("");
        assert!(channels(&v, 48).is_err());
    }
    #[test]
    fn candidate_hash_is_checked_before_holdout_access_and_no_overwrite() {
        let dir =
            std::env::temp_dir().join(format!("integrated-fit-synthetic-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let candidate_path = dir.join("candidate.json");
        let output_path = dir.join("output.json");
        // No training/holdout paths exist: bad candidate binding must fail first.
        fs::write(&candidate_path, b"{}").unwrap();
        let args = vec![
            "fit".into(),
            "validate".into(),
            candidate_path.to_str().unwrap().into(),
            "wrong-hash".into(),
            dir.join("absent-training.json").to_str().unwrap().into(),
            dir.join("absent-holdout.json").to_str().unwrap().into(),
            output_path.to_str().unwrap().into(),
        ];
        assert_eq!(
            run(&args).unwrap_err(),
            "frozen candidate hash mismatch before holdout access"
        );
        assert!(!output_path.exists());
        assert!(write_new(candidate_path.to_str().unwrap(), &json!({"different":true})).is_err());
        assert_eq!(fs::read(&candidate_path).unwrap(), b"{}");
        fs::remove_dir_all(dir).unwrap();
    }
    #[test]
    fn effective_support_controls_holdout_indices() {
        let all = geometric(2, 64);
        let tr = channel(all[..48].to_vec());
        let tail = channel(all[48..].to_vec());
        // A nominal shift-4 recurrence whose actual support is only shifts0/1.
        let c = Candidate {
            kind: Kind::Recurrence,
            order: 4,
            degree: 0,
            a: vec![Rational::from(-2), Q.one(), Q.zero(), Q.zero(), Q.zero()],
        };
        let v = validate_holdout(&[tr.clone()], &[tail.clone()], &c).unwrap();
        assert_eq!(v["status"], "heldout_pass_candidate_only");
        assert_eq!(v["rows"][0]["row"], 47);
        assert_eq!(v["rows"][15]["row"], 62);
        // t^2 times the same valid theta equation has min shift2, not0.
        let c = Candidate {
            kind: Kind::Theta,
            order: 1,
            degree: 3,
            a: vec![
                Q.zero(),
                Q.zero(),
                Q.zero(),
                Q.zero(),
                Q.zero(),
                Q.one(),
                Rational::from(-2),
                Rational::from(-2),
            ],
        };
        let v = validate_holdout(&[tr], &[tail], &c).unwrap();
        assert_eq!(v["status"], "heldout_pass_candidate_only");
        assert_eq!(v["rows"][0]["row"], 50);
        assert_eq!(v["rows"][15]["row"], 65);
    }
    #[test]
    fn singular_heldout_operator_row_is_not_counted_as_a_test() {
        let all = geometric(2, 64);
        // (n-47)*(c[n+1]-2c[n]) has a zero first nominal heldout row.
        let c = Candidate {
            kind: Kind::Recurrence,
            order: 1,
            degree: 1,
            a: vec![
                Rational::from(94),
                Rational::from(-2),
                Rational::from(-47),
                Q.one(),
            ],
        };
        let v = validate_holdout(
            &[channel(all[..48].to_vec())],
            &[channel(all[48..].to_vec())],
            &c,
        )
        .unwrap();
        assert_eq!(v["status"], "heldout_inconclusive_operator_singularity");
        assert_eq!(v["genuine_heldout_rows"], 15);
        assert_eq!(v["rows"][0]["genuine_heldout_row"], false);
    }
    #[test]
    fn integer_grading_requires_declared_proof_without_losing_fractional_branch() {
        let mut v = json!({"schema":"integrated-rational-prefix.v1","dimension":"13/2","q":1,"first_index":0,"series_identity":"synthetic-q1","channels":[{"id":"branch","beta":"-1/4","log_power":0,"coefficients":vec!["1";48]}]});
        assert!(channels(&v, 48).is_err());
        v["grading_proof_identity"] = json!("synthetic-hard-reversal-offset-proof");
        let ch = channels(&v, 48).unwrap();
        assert_eq!(ch[0].beta, rat("-1/4").unwrap());
        let (op, _) = fit(&ch, Kind::Theta).unwrap();
        let op = op.unwrap();
        for n in 0..48 {
            assert!(Q.is_zero(&dot(&row(&ch[0], n, op.kind, op.order, op.degree), &op.a)));
        }
    }
}
