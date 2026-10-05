use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use idenso::color::ColorSimplifier;
use std::{collections::HashMap, path::Path};
use symbolica::prelude::*;

fn atom(s: &str) -> Atom {
    Atom::parse(s, "UFO", Default::default()).unwrap()
}
fn decimal_atom(s: &str) -> Atom {
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
fn abs(x: &Float) -> Float {
    Float::with_val(384, x.as_raw().clone().abs())
}
fn validate_input(input: &serde_json::Value) -> Result<(), Box<dyn std::error::Error>> {
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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    let directory = args.get(1).ok_or("provide a scope directory")?;
    let base = Path::new(directory);
    let model = std::sync::Arc::new(Model::from_json(&std::fs::read_to_string(
        base.join("native-model.json"),
    )?)?);
    let input: serde_json::Value = serde_json::from_slice(&std::fs::read(
        args.get(2)
            .ok_or("provide the coherent amplitude input JSON")?,
    )?)?;
    validate_input(&input)?;
    let process = feynkit_generator::Process::new(["G", "G"], ["G", "H"]).with_filters(
        vec![],
        Some(vec!["GGGHEWZZ".into(), "GGGHEWWW".into()]),
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
    std::fs::write(
        base.join("me-symbolic.atom"),
        expression.to_canonical_string(),
    )?;
    eprintln!(
        "exact model-expanded ME {} bytes",
        expression.to_canonical_string().len()
    );
    let replacements: Vec<_> = ff_values
        .iter()
        .map(|(v, x)| Replacement::new(v.to_pattern(), x.clone()))
        .collect();
    let value_expression = expression.replace_multiple(&replacements);
    std::fs::write(
        base.join("me-values.atom"),
        value_expression.to_canonical_string(),
    )?;
    let mut results = Vec::new();
    for bits in [256, 384] {
        let empty: HashMap<Atom, Complex<Float>> = HashMap::new();
        let value = value_expression.evaluate_with_prec(&empty, bits)?;
        println!("ME {bits} {value}");
        results.push(serde_json::json!({"bits":bits,"real":value.re.to_string(),"imaginary":value.im.to_string()}));
    }
    let variables: Vec<_> = ff_values.iter().map(|(v, _)| v.clone()).collect();
    let variable_symbols: Vec<_> = variables
        .iter()
        .map(|v| match v.as_view() {
            AtomView::Var(v) => v.get_symbol(),
            _ => unreachable!(),
        })
        .collect();
    let grouped = expression.expand().to_polynomial_in_vars::<i32>(&variables);
    let mut term_count = 0;
    for term in grouped.into_iter() {
        assert!(term.exponents.iter().all(|e| *e >= 0));
        assert_eq!(
            term.exponents.iter().sum::<i32>(),
            2,
            "ME must be exactly homogeneous quadratic in the real FF components"
        );
        let symbols = term.coefficient.get_all_symbols(false);
        assert!(
            variable_symbols.iter().all(|v| !symbols.contains(v)),
            "nonpolynomial FF dependency hidden in coefficient"
        );
        term_count += 1;
    }
    eprintln!(
        "native polynomial certifies {term_count} homogeneous quadratic terms in 16 real FF components"
    );
    let bits = 384;
    let empty: HashMap<Atom, Complex<Float>> = HashMap::new();
    let eval = |a: &Atom| -> Result<Complex<Float>, Box<dyn std::error::Error>> {
        Ok(a.evaluate_with_prec(&empty, bits)?)
    };
    let errors: Vec<_> = ff_errors
        .iter()
        .map(|e| eval(&decimal_atom(e.as_str().unwrap())).map(|v| v.re))
        .collect::<Result<_, _>>()?;
    if errors
        .iter()
        .any(|v| !v.as_raw().is_finite() || v < &Float::with_val(bits, 0))
    {
        return Err("form-factor absolute allowances must be finite and nonnegative".into());
    }
    let mut linear = Float::with_val(bits, 0);
    let mut quadratic = Float::with_val(bits, 0);
    let mut gradients = Vec::new();
    for (i, &v) in variable_symbols.iter().enumerate() {
        let derivative = expression.derivative(v);
        let gradient = eval(&derivative.replace_multiple(&replacements))?;
        assert!(gradient.im.is_zero());
        let weight = abs(&gradient.re);
        linear += weight * &errors[i];
        gradients.push(serde_json::json!({"parameter":ff_names[i],"derivative":gradient.re.to_string(),"component_absolute_allowance":errors[i].to_string()}));
        for (j, &w) in variable_symbols.iter().enumerate() {
            let hessian = eval(&derivative.derivative(w))?;
            assert!(hessian.im.is_zero());
            quadratic += abs(&hessian.re) * &errors[i] * &errors[j] / Float::with_val(bits, 2);
        }
    }
    let rounding = eval(&atom("10^(-65)"))?.re;
    let total = linear.clone() + &quadratic + &rounding;
    let value = eval(&value_expression)?;
    let relative = total.clone() / abs(&value.re);
    let oracle: serde_json::Value = serde_json::from_slice(&std::fs::read(
        args.get(3)
            .ok_or("provide the independent QP oracle result JSON")?,
    )?)?;
    let oracle_string = oracle["binary128"]["SMATRIX"][0].as_str().unwrap();
    let oracle_value = eval(&decimal_atom(oracle_string))?.re;
    let difference = abs(&(value.re.clone() - &oracle_value));
    if !value.re.as_raw().is_finite() || !value.im.is_zero() || value.re.is_zero() {
        return Err("expected a finite, real, nonzero squared matrix element".into());
    }
    let arithmetic_change =
        abs(&(eval(&decimal_atom(results[0]["real"].as_str().unwrap()))?.re - &value.re));
    if arithmetic_change > rounding {
        return Err("256/384-bit arithmetic refinement exceeded its reserve".into());
    }
    let oracle_relative_tolerance = eval(&atom("10^(-30)"))?.re;
    if difference > oracle_relative_tolerance * abs(&oracle_value) {
        return Err(
            "native/oracle matrix elements disagree at the 30-digit arithmetic check".into(),
        );
    }
    let mut verified_relative_digits = 0;
    for digits in 1..=100 {
        if relative <= eval(&atom(&format!("10^(-{digits})")))?.re {
            verified_relative_digits = digits;
        } else {
            break;
        }
    }
    let report = serde_json::json!({
        "status":"passed", "observable":"Pure |A_W+A_Z|^2, incoming two-gluon spin/color average only (IDEN256), final-state sum; no HEFT interference",
        "native_owner_chain":"FeynKit Process -> Amplitude::squared -> sum_spins/sum_colors -> Idenso typed expanded/simplify_algebra -> native SU(3) Casimir formulas -> Kinematics::apply -> Model::expand_couplings and exact model parameter definitions -> Symbolica evaluate_with_prec",
        "evaluations":results,"verified_relative_digits":verified_relative_digits,"arithmetic_refinement_absolute_difference":arithmetic_change.to_string(),"oracle_acceptance_relative_tolerance":"1e-30","oracle_binary128":oracle_string,"oracle_absolute_difference":difference.to_string(),"oracle_relative_difference":(difference.clone()/abs(&value.re)).to_string(),
        "uncertainty":{"metric":"absolute error, conditional on supplied complex FF absolute allowances; each real/imag component conservatively gets the same complex allowance", "formula":"sum_i |df/dx_i| delta_i + 1/2 sum_ij |d2f/dx_i dx_j| delta_i delta_j; exact quadratic dependence certified by native polynomial grouping", "linear":linear.to_string(),"quadratic":quadratic.to_string(),"arithmetic_reserve":rounding.to_string(),"absolute_total":total.to_string(),"relative_total":relative.to_string(),"source_scope":"The model/scientific phase-point inputs are exact definitions; supplied FF numerical evidence remains the limiting uncertainty. No 20-digit claim follows from working precision or oracle agreement alone.","gradients":gradients},
        "homogeneous_quadratic_term_count":term_count,"physical_s":s.to_string(),"physical_t":t.to_string(),"physical_u":u.to_string(),"MH_squared":"1","MW_squared":"5399/13074","MZ_squared":"7775/14631","alpha":"1/128","alpha_s":"118/1000"
    });
    std::fs::write(
        base.join("me-result.json"),
        serde_json::to_vec_pretty(&report)?,
    )?;
    println!(
        "absolute propagated bound {total}; relative {relative}; oracle absolute {difference}"
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn scientific_decimal_inputs_remain_exact() {
        for (input, expected) in [
            ("-1.25e-3", "-1/800"),
            ("3.25E+2", "325"),
            ("-.5", "-1/2"),
            ("0.000", "0"),
        ] {
            assert_eq!(decimal_atom(input), atom(expected));
        }
    }
    #[test]
    fn input_requires_complete_distinct_w_z_form_factors_at_unit_higgs_mass() {
        let factor = |i| serde_json::json!({"form_factor":i,"normalized_value":{"real":"0","imaginary":"0"},"normalized_absolute_error":"1e-30"});
        let block = |mass| serde_json::json!({"mass":mass,"form_factors":(1..=4).map(factor).collect::<Vec<_>>()});
        let good = serde_json::json!({"physical_coordinates":{"physical_s_t_MH_squared":["4","-1","1"]},"form_factors":[block("W"),block("Z")]});
        assert!(validate_input(&good).is_ok());
        let mut bad = good.clone();
        bad["form_factors"][1]["mass"] = "W".into();
        assert!(validate_input(&bad).is_err());
        let mut bad = good.clone();
        bad["form_factors"][1]["form_factors"][3]["form_factor"] = 1.into();
        assert!(validate_input(&bad).is_err());
        let mut bad = good;
        bad["physical_coordinates"]["physical_s_t_MH_squared"][2] = "2".into();
        assert!(validate_input(&bad).is_err());
    }
}
