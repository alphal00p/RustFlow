mod scalar;
use feynkit_model::Model;
use scalar::{abs, atom, decimal_atom, evaluate_exact, prepare_scalar};
use std::{collections::BTreeSet, path::Path};
use symbolica::prelude::*;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn exact_rational(text: &str) -> Result<Atom> {
    let mut parts = text.split('/');
    let numerator = parts.next().ok_or("missing rational numerator")?;
    let denominator = parts.next();
    let integer = |s: &str| {
        let s = s
            .strip_prefix('-')
            .or_else(|| s.strip_prefix('+'))
            .unwrap_or(s);
        !s.is_empty() && s.bytes().all(|c| c.is_ascii_digit())
    };
    if parts.next().is_some()
        || !integer(numerator)
        || denominator.is_some_and(|d| !integer(d) || atom(d).is_zero())
    {
        return Err(
            "HEFT values and phase coordinates must be finite exact integer/rational strings"
                .into(),
        );
    }
    Ok(atom(text))
}

fn read_heft(input: &serde_json::Value, s: &Atom, t: &Atom) -> Result<Vec<Atom>> {
    for (name, expected) in [
        ("physical_s", s),
        ("physical_t", t),
        ("MH_squared", &Atom::one()),
    ] {
        let value = input[name]
            .as_str()
            .ok_or("missing exact HEFT phase coordinate")?;
        if exact_rational(value)? != *expected {
            return Err("HEFT coefficients belong to a different phase point".into());
        }
    }
    let entries = input["coefficients"]
        .as_array()
        .ok_or("missing exact HEFT coefficients")?;
    if entries.len() != 4 {
        return Err("exact HEFT input requires four coefficients".into());
    }
    let mut seen = BTreeSet::new();
    let mut values = vec![Atom::Zero; 4];
    for entry in entries {
        let index = entry["index"]
            .as_u64()
            .ok_or("missing HEFT coefficient index")?;
        if !(1..=4).contains(&index) || !seen.insert(index) {
            return Err("HEFT indices must be distinct and cover1 through4".into());
        }
        values[index as usize - 1] = exact_rational(
            entry["value"]
                .as_str()
                .ok_or("missing exact HEFT coefficient")?,
        )?;
    }
    Ok(values)
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 6 {
        return Err("usage: evaluate-heft <scope-directory> <coherent-input.json> <exact-heft.json> <oracle.json> <output-directory>".into());
    }
    let scope = Path::new(&args[1]);
    let output = Path::new(&args[5]);
    let model = std::sync::Arc::new(Model::from_json(&std::fs::read_to_string(
        scope.join("native-model.json"),
    )?)?);
    let input: serde_json::Value = serde_json::from_slice(&std::fs::read(&args[2])?)?;
    let heft_input: serde_json::Value = serde_json::from_slice(&std::fs::read(&args[3])?)?;
    let oracle: serde_json::Value = serde_json::from_slice(&std::fs::read(&args[4])?)?;
    // Parse complete numerical input before native tensor work.
    scalar::validate_input(&input)?;
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
    let heft_values = read_heft(&heft_input, &s, &t)?;
    let mut prepared = prepare_scalar(model, &input, true)?;
    let marker = atom("heft_marker");
    let marker_symbol = match marker.as_view() {
        AtomView::Var(v) => v.get_symbol(),
        _ => unreachable!(),
    };
    let mut rules = Vec::new();
    for (i, value) in heft_values.iter().enumerate() {
        let variable = atom(&format!("GGGH_HEFT_ForFac{}", i + 1));
        prepared.expression = prepared
            .expression
            .replace(variable.conj())
            .with(variable.clone());
        rules.push(Replacement::new(variable.to_pattern(), &marker * value));
    }
    let marked = prepared.expression.replace_multiple(&rules).expand();
    let polynomial = marked.to_polynomial_in_vars::<i32>(std::slice::from_ref(&marker));
    let mut projected = [Atom::Zero, Atom::Zero, Atom::Zero];
    for term in polynomial.into_iter() {
        let degree = term.exponents[0];
        if !(0..=2).contains(&degree)
            || term
                .coefficient
                .get_all_symbols(false)
                .contains(&marker_symbol)
        {
            return Err("native scalar expression is not quadratic in the HEFT marker".into());
        }
        projected[degree as usize] += term.coefficient;
    }
    let reconstruction = &projected[0] + &marker * &projected[1] + marker.pow(2) * &projected[2];
    if !(&marked - reconstruction).expand().is_zero() {
        return Err("exact scalar coupling-order reconstruction failed".into());
    }
    let replacements: Vec<_> = prepared
        .ff_values
        .iter()
        .map(|(v, x)| Replacement::new(v.to_pattern(), x.clone()))
        .collect();
    let variables: Vec<_> = prepared.ff_values.iter().map(|(v, _)| v.clone()).collect();
    let symbols: Vec<_> = variables
        .iter()
        .map(|v| match v.as_view() {
            AtomView::Var(v) => v.get_symbol(),
            _ => unreachable!(),
        })
        .collect();
    let reserve = evaluate_exact(&atom("10^(-65)"), 384)?.re;
    let oracle_tolerance = evaluate_exact(&atom("10^(-30)"), 384)?.re;
    let mut results = Vec::new();
    for (index, expression) in projected.iter().enumerate() {
        let label = ["ew_squared", "heft_ew_interference", "heft_squared"][index];
        let value_expression = expression.replace_multiple(&replacements);
        let low = evaluate_exact(&value_expression, 256)?;
        let high = evaluate_exact(&value_expression, 384)?;
        let change = abs(&(high.re.clone() - &low.re)) + abs(&(high.im.clone() - &low.im));
        if !high.re.as_raw().is_finite()
            || !high.im.as_raw().is_finite()
            || high.re.is_zero()
            || abs(&high.im) > reserve
            || abs(&low.im) > reserve
            || change > reserve
        {
            return Err("native256/384-bit arithmetic refinement or reality check failed".into());
        }
        let mut result = serde_json::json!({"observable":label,"evaluations":[{"bits":256,"real":low.re.to_string(),"imaginary":low.im.to_string()},{"bits":384,"real":high.re.to_string(),"imaginary":high.im.to_string()}],"arithmetic_change_l1":change.to_string(),"arithmetic_reserve":reserve.to_string()});
        if index > 0 {
            let source = oracle["binary128"][label]
                .as_str()
                .ok_or("missing independent binary128 oracle observable")?;
            let expected = evaluate_exact(&decimal_atom(source), 384)?;
            let difference = abs(&(high.re.clone() - &expected.re));
            if !expected.re.as_raw().is_finite()
                || !expected.im.is_zero()
                || expected.re.is_zero()
                || difference > oracle_tolerance.clone() * abs(&expected.re)
            {
                return Err(
                    "native/original observable fails30-digit relative arithmetic comparison"
                        .into(),
                );
            }
            result["original_binary128"] = source.into();
            result["oracle_absolute_difference"] = difference.to_string().into();
            result["oracle_relative_difference"] =
                (difference / abs(&expected.re)).to_string().into();
            result["oracle_relative_tolerance"] = "1e-30".into();
        }
        if index == 1 {
            let polynomial = expression.to_polynomial_in_vars::<i32>(&variables);
            let mut term_count = 0;
            for term in polynomial.into_iter() {
                if term.exponents.iter().sum::<i32>() != 1
                    || term.exponents.iter().any(|e| *e < 0)
                    || symbols
                        .iter()
                        .any(|v| term.coefficient.get_all_symbols(false).contains(v))
                {
                    return Err(
                        "interference is not exactly linear in the real EWFF components".into(),
                    );
                }
                term_count += 1;
            }
            let mut bound = reserve.clone();
            let mut gradients = Vec::new();
            for (i, v) in symbols.iter().enumerate() {
                let derivative = expression.derivative(*v);
                for w in &symbols {
                    if !derivative.derivative(*w).expand().is_zero() {
                        return Err("nonzero EW interference Hessian".into());
                    }
                }
                let gradient = evaluate_exact(&derivative.replace_multiple(&replacements), 384)?;
                let error = evaluate_exact(
                    &decimal_atom(
                        prepared.ff_errors[i]
                            .as_str()
                            .ok_or("missing EWFF allowance")?,
                    ),
                    384,
                )?
                .re;
                if !error.as_raw().is_finite() || error < Float::with_val(384, 0) {
                    return Err("EWFF allowances must be finite and nonnegative".into());
                }
                bound += (abs(&gradient.re) + abs(&gradient.im)) * &error;
                gradients.push(serde_json::json!({"parameter":prepared.ff_names[i],"real":gradient.re.to_string(),"imaginary":gradient.im.to_string(),"component_absolute_allowance":error.to_string()}));
            }
            let relative = bound.clone() / abs(&high.re);
            let mut digits = 0;
            for d in 1..=100 {
                if relative <= evaluate_exact(&atom(&format!("10^(-{d})")), 384)?.re {
                    digits = d
                } else {
                    break;
                }
            }
            result["uncertainty"] = serde_json::json!({"metric":"absolute error; relative error divides by this observable's magnitude","formula":"sum_i (abs(Re gradient_i)+abs(Im gradient_i))*delta_i +1e-65; exact linear polynomial and vanishing Hessian certified","conditional_verified_relative_digits":digits,"absolute_bound":bound.to_string(),"relative_bound":relative.to_string(),"linear_terms":term_count,"gradients":gradients,"scope":"Conditional on supplied complex EW form-factor allowances. Exact HEFT rational coefficients/model/phase-point inputs remain exact until MPFR evaluation. Each real and imaginary FF component conservatively receives the full complex allowance."});
        }
        if index == 2
            && symbols
                .iter()
                .any(|v| !expression.derivative(*v).expand().is_zero())
        {
            return Err("HEFT LO unexpectedly depends on EW form factors".into());
        }
        results.push(result);
    }
    std::fs::create_dir_all(output)?;
    let report = serde_json::json!({"status":"passed","observable_scope":"Infinite-top HEFT LO and2Re(A_HEFT conjugate(A_W+A_Z)); incoming two-gluon spin/color averages and final-state sums; no finite-top QCD claim","native_owner_chain":"Shared Process/Amplitude::squared/state sums/Idenso contraction/Kinematics/Model coupling expansion, then native Symbolica scalar polynomial projection and exact-input MPFR evaluation","exact_scalar_order_reconstruction":true,"physical_s":prepared.s.to_string(),"physical_t":prepared.t.to_string(),"physical_u":prepared.u.to_string(),"exact_heft_inputs":heft_input,"results":results,"limits":"HEFT LO precision is demonstrated by256/384-bit refinement and113-bit independent original arithmetic. Interference physical accuracy is limited by explicitly propagated supplied EWFF allowances; oracle agreement alone does not establish it."});
    std::fs::write(
        output.join("heft-result.json"),
        serde_json::to_vec_pretty(&report)?,
    )?;
    println!(
        "HEFT LO and interference arithmetic checks passed; see heft-result.json for propagated EW-input accuracy"
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn heft_inputs_are_exact() {
        assert_eq!(exact_rational("-3/7").unwrap(), atom("-3/7"));
        for s in ["0.5", "1e-20", "x", "1/0", "1/2/3", "NaN", ""] {
            assert!(exact_rational(s).is_err(), "{s}");
        }
    }
    #[test]
    fn heft_requires_matching_point_and_complete_unique_indices() {
        let mut input = serde_json::json!({"physical_s":"4","physical_t":"-1","MH_squared":"1","coefficients":(1..=4).map(|i|serde_json::json!({"index":i,"value":"1/3"})).collect::<Vec<_>>()});
        assert!(read_heft(&input, &atom("4"), &atom("-1")).is_ok());
        assert!(read_heft(&input, &atom("5"), &atom("-1")).is_err());
        input["coefficients"][3]["index"] = 1.into();
        assert!(read_heft(&input, &atom("4"), &atom("-1")).is_err());
    }
}
