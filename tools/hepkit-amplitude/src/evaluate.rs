mod scalar;
use feynkit_model::Model;
use scalar::{
    ScalarPreparation, abs, atom, decimal_atom, evaluate_exact, prepare_scalar, validate_input,
};
use std::path::Path;
use symbolica::prelude::*;

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
    let ScalarPreparation {
        expression,
        ff_values,
        ff_errors,
        ff_names,
        s,
        t,
        u,
    } = prepare_scalar(model, &input, false)?;
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
        let value = evaluate_exact(&value_expression, bits)?;
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
    let eval = |a: &Atom| -> Result<Complex<Float>, Box<dyn std::error::Error>> {
        evaluate_exact(a, bits)
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
