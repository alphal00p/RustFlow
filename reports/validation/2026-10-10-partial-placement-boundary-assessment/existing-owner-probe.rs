//! Standalone checks of existing boundary conversion and ordinary zero owners.
//! This is not a partial-placement boundary admission or an occupied zero test.
use symbolica::prelude::*;
use symbolica_amflow::{IntegralFamily, integrand, recursive};
fn template(loops: &[&str], external: &[&str], gram: Vec<Vec<Atom>>) -> IntegralFamily {
    IntegralFamily {
        name: "partial_boundary_owner_probe".into(),
        loops: loops.iter().map(|x| (*x).into()).collect(),
        external: external.iter().map(|x| (*x).into()).collect(),
        external_gram: gram,
        propagators: vec![],
        physical_propagators: 0,
        epsilon: symbol!("partial_boundary_probe::eps"),
        dimension: 4,
    }
}
fn inspect(
    name: &str,
    expression: Atom,
    variables: &[Atom],
    family: &IntegralFamily,
) -> serde_json::Value {
    let terms = integrand::to_integrals(&expression, variables, family, 4096).unwrap();
    let mut rebuilt = Atom::zero();
    let mut output = Vec::new();
    for term in terms {
        term.family.validate().unwrap();
        let mut value = term.coefficient.clone();
        let mut denominators = Vec::new();
        for (d, &n) in term.family.propagators.iter().zip(&term.integral.0) {
            let denominator = d
                .scalar_products
                .iter()
                .zip(variables)
                .fold(d.constant.clone(), |s, (c, v)| s + c * v);
            value *= denominator.clone().pow(-i64::from(n));
            denominators.push(denominator.to_canonical_string());
        }
        rebuilt += value;
        output.push(serde_json::json!({"coefficient":term.coefficient.to_canonical_string(),"powers":term.integral.0,"denominators":denominators,"native_unrestricted_proved_zero":recursive::scaleless(&term.family,&term.integral).unwrap()}));
    }
    assert!((rebuilt - &expression).together().cancel().is_zero());
    serde_json::json!({"name":name,"expression":expression.to_canonical_string(),"exact_reconstruction":true,"terms":output})
}
fn main() {
    let x = Atom::var(symbol!("pb_x"));
    let y = Atom::var(symbol!("pb_y"));
    let z = Atom::var(symbol!("pb_z"));
    let s = Atom::var(symbol!("pb_s"));
    let f2 = template(&["K", "L"], &[], vec![]);
    let vars = vec![x.clone(), z.clone(), y.clone()];
    let mut tests = Vec::new();
    tests.push(inspect(
        "hard_duplicate_momentum_unit_mass_gap",
        parse!("1/((pb_x-1)*pb_x*pb_y*(pb_x+pb_y-2*pb_z))"),
        &vars,
        &f2,
    ));
    tests.push(inspect(
        "hard_raised_duplicate_momentum_with_numerator",
        parse!("(pb_z+pb_y)/((pb_x-1)^2*pb_x^3*pb_y*(pb_x+pb_y-2*pb_z))"),
        &vars,
        &f2,
    ));
    tests.push(inspect(
        "ordinary_massless_two_loop_vacuum",
        parse!("1/(pb_x*pb_y*(pb_x+pb_y-2*pb_z))"),
        &vars,
        &f2,
    ));
    tests.push(inspect(
        "ordinary_free_virtual_direction",
        parse!("pb_z^2/(pb_y-1)"),
        &vars,
        &f2,
    ));
    let f1 = template(&["K"], &["q"], vec![vec![s.clone()]]);
    let v1 = vec![x.clone(), z.clone()];
    tests.push(inspect(
        "translated_massless_tadpole",
        parse!("1/(pb_x-2*pb_z+pb_s)"),
        &v1,
        &f1,
    ));
    tests.push(inspect(
        "translated_massive_tadpole_negative",
        parse!("1/(pb_x-2*pb_z+pb_s-1)"),
        &v1,
        &f1,
    ));
    tests.push(inspect(
        "incompatible_translations_external_scale_negative",
        parse!("1/(pb_x*(pb_x-2*pb_z+pb_s))"),
        &v1,
        &f1,
    ));
    let all_zero = |case: &serde_json::Value| {
        case["terms"]
            .as_array()
            .unwrap()
            .iter()
            .all(|t| t["native_unrestricted_proved_zero"] == true)
    };
    assert!(all_zero(&tests[2]) && all_zero(&tests[3]) && all_zero(&tests[4]));
    assert!(!all_zero(&tests[5]) && !all_zero(&tests[6]));
    let hard = tests[0]["terms"].as_array().unwrap();
    assert_eq!(hard.len(), 2);
    assert_eq!(
        hard.iter()
            .filter(|t| t["native_unrestricted_proved_zero"] == true)
            .count(),
        1
    );
    std::fs::write(std::env::args().nth(1).unwrap(),serde_json::to_vec_pretty(&serde_json::json!({"scope":"Existing public native partial fractions and ordinary unrestricted zero owner only; no occupied measure admission or numerical period","tests":tests,"passed":7})).unwrap()).unwrap();
}
