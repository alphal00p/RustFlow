//! Component validation; recorded form factors are not used by production seed generation.
use std::{collections::BTreeMap, sync::Arc};
use symbolica_amflow::symbolica::prelude::*;
use symbolica_amflow::{
    Precision, RunContext,
    gg_hg::amplitude::{HiggsJetAmplitude, HiggsJetFormFactors},
};

fn atom(s: &str) -> Atom {
    Atom::parse(s, "UFO", Default::default()).unwrap()
}

#[test]
fn native_coherent_amplitude_and_reused_kinematic_kernel() {
    let input: serde_json::Value =
        serde_json::from_str(include_str!("../fixtures/gg-hg/amplitude-validation.json")).unwrap();
    let model = Arc::new(
        feynkit_model::Model::from_json(include_str!("../fixtures/gg-hg/native-model.json"))
            .unwrap(),
    );
    let context = RunContext::default();
    let kernel = HiggsJetAmplitude::new(model.clone(), &context).unwrap();
    assert!(!kernel.diagrams().is_empty());
    assert!(
        kernel
            .diagrams()
            .iter()
            .all(|g| std::ptr::eq(g.model(), model.as_ref()))
    );
    let p = Precision::decimal(110).unwrap();
    let mut values = vec![];
    let mut errors = vec![];
    for tag in ["W", "Z"] {
        let block = input["form_factors"]
            .as_array()
            .unwrap()
            .iter()
            .find(|block| block["mass"] == tag)
            .unwrap();
        for index in 1..=4 {
            let coefficient = block["values"]
                .as_array()
                .unwrap()
                .iter()
                .find(|v| v["index"] == index)
                .unwrap();
            values.push(
                p.parse(
                    coefficient["real"].as_str().unwrap(),
                    coefficient["imaginary"].as_str().unwrap(),
                )
                .unwrap(),
            );
            errors.push(
                p.parse(coefficient["absolute_error"].as_str().unwrap(), "0")
                    .unwrap()
                    .re,
            );
        }
    }
    let factors = HiggsJetFormFactors {
        values: values.try_into().unwrap(),
        absolute_errors: errors.try_into().unwrap(),
        provenance: "archived native coherent transport; amplitude component test only".into(),
    };
    let parameters = [
        ("aEWM1", "128"),
        ("aS", "118/1000"),
        ("MZ", "(7775/14631)^(1/2)"),
        (
            "Gf",
            "𝜋*(1/128)*(7775/14631)/(2^(1/2)*(5399/13074)*((7775/14631)-(5399/13074)))",
        ),
    ]
    .into_iter()
    .map(|(name, value)| (atom(name), atom(value)))
    .collect::<BTreeMap<_, _>>();
    let point = input["physical_s_t_MH_squared"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| atom(v.as_str().unwrap()))
        .collect::<Vec<_>>();
    let result = kernel
        .evaluate(
            &point[0],
            &point[1],
            &point[2],
            &factors,
            &parameters,
            20,
            40,
            &context,
        )
        .unwrap();
    for (label, result) in [
        ("electroweak_squared", &result.electroweak_squared),
        ("interference", &result.interference),
        ("effective_squared", &result.effective_squared),
    ] {
        let expected = p
            .parse(input["expected_observables"][label].as_str().unwrap(), "0")
            .unwrap();
        let difference = p.norm(&p.sub(&Complex::new(result.value.clone(), p.real(0)), &expected));
        assert!(
            difference
                <= p.mul(&expected, &Complex::new(p.tolerance(30), p.real(0)))
                    .re,
            "{label}: {} vs {expected}",
            result.value
        );
    }
    // The archived EW inputs support 19 relative digits, even though arithmetic
    // agrees with the independent oracle at more than 30 digits.
    assert_eq!(
        result.electroweak_squared.verified_relative_digits,
        Some(19)
    );
    assert!(result.interference.verified_relative_digits.unwrap() >= 20);

    // Reuse the native kernel at two more points, including a changed Higgs mass.
    // The archived EW inputs are not physical at these points; only HEFT is checked.
    // The HEFT square
    // has the independent crossing-symmetric shape (s^4+t^4+u^4+MH^8)/(s*t*u).
    let shape = |s: &Atom, t: &Atom, h: &Atom| {
        let u = h - s - t;
        (s.pow(4) + t.pow(4) + u.pow(4) + h.pow(4)) / (s * t * u)
    };
    for (s, t, h) in [("2", "-1/3", "1"), ("8", "-4/3", "4")] {
        let [s, t, h] = [s, t, h].map(atom);
        let additional = kernel
            .evaluate(&s, &t, &h, &factors, &parameters, 20, 40, &context)
            .unwrap();
        let ratio = p
            .eval(
                &(shape(&s, &t, &h) / shape(&point[0], &point[1], &point[2])),
                &Default::default(),
            )
            .unwrap();
        let expected = p.mul(
            &Complex::new(result.effective_squared.value.clone(), p.real(0)),
            &ratio,
        );
        assert!(
            p.norm(&p.sub(
                &Complex::new(additional.effective_squared.value, p.real(0)),
                &expected
            )) <= p
                .mul(&expected, &Complex::new(p.tolerance(35), p.real(0)))
                .re
        );
    }
}
