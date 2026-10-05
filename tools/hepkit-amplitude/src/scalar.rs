use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use idenso::color::ColorSimplifier;
use std::{collections::HashMap, sync::Arc};
use symbolica::prelude::*;

pub fn atom(s: &str) -> Atom {
    Atom::parse(s, "UFO", Default::default()).unwrap()
}
pub fn decimal_atom(s: &str) -> Atom {
    // Preserve supplied decimal mantissas exactly. Symbolica’s ordinary Atom
    // parser intentionally turns decimal literals into rounded Float values.
    let (mantissa, exponent) = s
        .split_once('e')
        .or_else(|| s.split_once('E'))
        .map(|(a, b)| (a, b.parse::<i64>().unwrap()))
        .unwrap_or((s, 0));
    let decimals = mantissa.split_once('.').map_or(0, |(_, b)| b.len() as i64);
    atom(&format!(
        "{}*10^({})",
        mantissa.replace('.', ""),
        exponent - decimals
    ))
}
pub fn abs(x: &Float) -> Float {
    Float::with_val(384, x.as_raw().clone().abs())
}
pub fn validate_input(input: &serde_json::Value) -> Result<(), Box<dyn std::error::Error>> {
    let coordinates = input["physical_coordinates"]["physical_s_t_MH_squared"]
        .as_array()
        .ok_or("missing physical coordinates")?;
    if coordinates.len() != 3
        || coordinates.iter().any(|v| !v.is_string())
        || atom(coordinates[2].as_str().unwrap()) != Atom::one()
    {
        return Err("this benchmark requires exact (s,t,MH_squared=1) coordinates".into());
    }
    let blocks = input["form_factors"]
        .as_array()
        .ok_or("missing form factors")?;
    let mut masses = std::collections::BTreeSet::new();
    for block in blocks {
        let mass = block["mass"].as_str().ok_or("missing W/Z mass label")?;
        if !["W", "Z"].contains(&mass) || !masses.insert(mass) {
            return Err("form factors require one W block and one Z block".into());
        }
        let factors = block["form_factors"]
            .as_array()
            .ok_or("missing form-factor block")?;
        let mut indices = std::collections::BTreeSet::new();
        for factor in factors {
            let index = factor["form_factor"]
                .as_u64()
                .ok_or("missing form-factor index")?;
            if !(1..=4).contains(&index) || !indices.insert(index) {
                return Err("each mass requires distinct form factors 1 through 4".into());
            }
            for field in ["real", "imaginary"] {
                factor["normalized_value"][field]
                    .as_str()
                    .ok_or("missing decimal form-factor component")?;
            }
            factor["normalized_absolute_error"]
                .as_str()
                .ok_or("missing form-factor error allowance")?;
        }
        if indices.len() != 4 {
            return Err("each mass requires four form factors".into());
        }
    }
    if masses.len() != 2 {
        return Err("both W and Z form factors are required".into());
    }
    Ok(())
}

pub struct ScalarPreparation {
    pub expression: Atom,
    pub ff_values: Vec<(Atom, Atom)>,
    pub ff_errors: Vec<serde_json::Value>,
    pub ff_names: Vec<String>,
    pub s: Atom,
    pub t: Atom,
    pub u: Atom,
}
/// Native process/tensor/model preparation shared by both optional observables.
pub fn prepare_scalar(
    model: Arc<Model>,
    input: &serde_json::Value,
    include_heft: bool,
) -> Result<ScalarPreparation, Box<dyn std::error::Error>> {
    validate_input(input)?;
    let process = feynkit_generator::Process::new(["G", "G"], ["G", "H"]).with_filters(
        vec![],
        Some(if include_heft {
            vec!["GGGHEWZZ".into(), "GGGHEWWW".into(), "GGGHHEFT".into()]
        } else {
            vec!["GGGHEWZZ".into(), "GGGHEWWW".into()]
        }),
        vec![],
    );
    let amplitude = process.generate_amplitude(
        model.clone(),
        &feynkit_generator::GenerationOptions::default()
            .with_loop_count(0, 0)?
            .threads(1),
        Default::default(),
    )?;
    let registered = amplitude
        .squared()?
        .sum_spins(
            &[0, 1, 2, 3],
            true,
            &Default::default(),
            &Default::default(),
        )?
        .sum_colors(&[0, 1, 2, 3], true)?;
    let scalar =
        idenso::tensor::SymbolicTensor::<spenso::structure::partial::PartialStructure>::infer(
            registered.expression().clone(),
        )?
        .expanded(None, false)?
        .simplify_algebra(&idenso::tensor::AlgebraSettings {
            contract: idenso::tensor::AlgebraContraction::Dots,
            ..idenso::tensor::AlgebraSettings::hep()
        })?
        .into_expression();
    let s = atom(
        input["physical_coordinates"]["physical_s_t_MH_squared"][0]
            .as_str()
            .unwrap(),
    );
    let t = atom(
        input["physical_coordinates"]["physical_s_t_MH_squared"][1]
            .as_str()
            .unwrap(),
    );
    let u = Atom::one() - &s - &t;
    let p: Vec<_> = (0..4)
        .map(|i| atom(&format!("gammalooprs::P({i})")))
        .collect();
    let products = [
        [Atom::Zero, &s / 2, -&t / 2, (&s + &t) / 2],
        [&s / 2, Atom::Zero, -&u / 2, (&s + &u) / 2],
        [-&t / 2, -&u / 2, Atom::Zero, (&s - Atom::one()) / 2],
        [
            (&s + &t) / 2,
            (&s + &u) / 2,
            (&s - Atom::one()) / 2,
            Atom::one(),
        ],
    ];
    let mut kinematics = Kinematics::new();
    for i in 0..4 {
        for j in i..4 {
            kinematics = kinematics.with_scalar_product(&p[i], &p[j], products[i][j].clone())?;
        }
    }
    let reduced = kinematics.apply(&scalar.to_cof_dimension_invariants());
    eprintln!(
        "native kinematic and SU(3) reduction {} bytes",
        reduced.to_canonical_string().len()
    );
    let mut expression = model.expand_couplings(&reduced);
    let definitions: Vec<_> = model
        .parameters()
        .iter()
        .filter_map(|p| {
            p.expression
                .as_ref()
                .map(|e| Replacement::new(atom(&p.name).to_pattern(), e.clone()))
        })
        .collect();
    for _ in 0..=model.parameters().len() {
        let next = expression.replace_multiple(&definitions);
        if next == expression {
            break;
        }
        expression = next;
    }
    let external = [
        ("aEWM1", "128"),
        ("aS", "118/1000"),
        ("MZ", "(7775/14631)^(1/2)"),
        ("MH", "1"),
        (
            "Gf",
            "𝜋*(1/128)*(7775/14631)/(2^(1/2)*(5399/13074)*((7775/14631)-(5399/13074)))",
        ),
    ];
    let rules: Vec<_> = external
        .iter()
        .map(|(n, e)| Replacement::new(atom(n).to_pattern(), atom(e)))
        .collect();
    expression = expression.replace_multiple(&rules);
    // Use Symbolica's native scalar conjugation after tensor algebra is closed.
    let wildcard = symbol!("scope_eval::value_");
    expression = expression
        .replace(function!(symbol!("spenso::conj"), wildcard))
        .with(Atom::var(wildcard).conj());
    // All FF components are declared real external UFO parameters.
    let mut ff_values = Vec::new();
    let mut ff_errors = Vec::new();
    let mut ff_names = Vec::new();
    for block in input["form_factors"].as_array().unwrap() {
        let tag = if block["mass"] == "W" {
            "GGGHEWWW"
        } else {
            "GGGHEWZZ"
        };
        for f in block["form_factors"].as_array().unwrap() {
            let number = f["form_factor"].as_u64().unwrap();
            for (suffix, component) in [("RE", "real"), ("IM", "imaginary")] {
                let variable = atom(&format!("{tag}_ForFac{number}_{suffix}"));
                expression = expression.replace(variable.conj()).with(variable.clone());
                let value = decimal_atom(f["normalized_value"][component].as_str().unwrap());
                ff_names.push(format!("{tag}_ForFac{number}_{suffix}"));
                ff_values.push((variable, value));
                ff_errors.push(f["normalized_absolute_error"].clone());
            }
        }
    }
    Ok(ScalarPreparation {
        expression,
        ff_values,
        ff_errors,
        ff_names,
        s,
        t,
        u,
    })
}
pub fn evaluate_exact(
    expression: &Atom,
    bits: u32,
) -> Result<Complex<Float>, Box<dyn std::error::Error>> {
    let empty: HashMap<Atom, Complex<Float>> = HashMap::new();
    Ok(expression.evaluate_with_prec(&empty, bits)?)
}
