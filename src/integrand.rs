//! Convert rational boundary integrands to independent denominator families.
use crate::algebra::{inverse, nullspace, rref};
use crate::family::substitute;
use crate::*;
use std::collections::BTreeMap;
use symbolica::prelude::*;

#[derive(Clone, Debug)]
pub struct IntegralTerm {
    pub coefficient: Atom,
    pub family: IntegralFamily,
    pub integral: Integral,
}

pub(crate) fn factors(a: &Atom) -> Vec<(Atom, i64)> {
    let list = if let AtomView::Mul(m) = a.as_view() {
        m.iter().map(|v| v.to_owned()).collect()
    } else {
        vec![a.clone()]
    };
    list.into_iter()
        .map(|v| {
            if let AtomView::Pow(w) = v.as_view() {
                let (base, exponent) = w.get_base_exp();
                if let Ok(power) = exponent.to_string().parse::<i64>() {
                    return (base.to_owned(), power);
                }
            }
            (v, 1)
        })
        .collect()
}

pub(crate) fn powers(monomial: &Atom, variables: &[Atom]) -> Result<Vec<i16>> {
    let mut out = vec![0_i16; variables.len()];
    for (base, power) in factors(monomial) {
        if base.is_one() {
            continue;
        }
        let j = variables
            .iter()
            .position(|v| v == &base)
            .ok_or_else(|| Error::Unsupported("nonpolynomial boundary numerator".into()))?;
        out[j] =
            i16::try_from(power).map_err(|_| Error::Limit("numerator power overflow".into()))?;
    }
    Ok(out)
}

fn affine(a: &Atom, variables: &[Atom]) -> Result<Propagator> {
    let mut out = Propagator {
        constant: Atom::new(),
        scalar_products: vec![Atom::new(); variables.len()],
    };
    for (m, c) in a.coefficient_list::<i32>(variables) {
        if m.is_one() {
            out.constant += c;
        } else if let Some(j) = variables.iter().position(|v| v == &m) {
            out.scalar_products[j] += c;
        } else {
            return Err(Error::Unsupported("nonaffine boundary denominator".into()));
        }
    }
    Ok(out)
}

/// Exact multivariate partial fractions, ISP completion, and numerator conversion.
/// `template` supplies the loop/external space; its propagators are replaced.
pub fn to_integrals(
    expression: &Atom,
    variables: &[Atom],
    template: &IntegralFamily,
    budget: usize,
) -> Result<Vec<IntegralTerm>> {
    if expression.is_zero() {
        return Ok(vec![]);
    }
    let mut numerator = Atom::num(1);
    let mut denominators = Vec::new();
    let mut indices = Vec::new();
    for (base, power) in factors(&expression.together().factor()) {
        let depends = variables.iter().any(|v| {
            !base
                .derivative(match v.as_view() {
                    AtomView::Var(v) => v.get_symbol(),
                    _ => unreachable!(),
                })
                .is_zero()
        });
        if power < 0 && depends {
            let mut propagator = affine(&base, variables)?;
            // Normalize the quadratic coefficient, preserving the physical sign.
            let scale = propagator
                .scalar_products
                .iter()
                .find(|v| !v.is_zero())
                .cloned()
                .ok_or_else(|| Error::InvalidInput("constant denominator".into()))?;
            propagator.constant = (&propagator.constant / &scale).together().cancel();
            for c in &mut propagator.scalar_products {
                *c = (&*c / &scale).together().cancel();
            }
            numerator *= scale.pow(power);
            denominators.push(propagator);
            indices.push(
                i16::try_from(-power)
                    .map_err(|_| Error::Limit("denominator power overflow".into()))?,
            );
        } else {
            numerator *= base.pow(power);
        }
    }
    let mut out = Vec::new();
    let mut remaining = budget;
    partial_fraction(
        denominators,
        indices,
        numerator,
        variables,
        template,
        &mut remaining,
        &mut out,
    )?;
    Ok(out)
}

fn partial_fraction(
    denominators: Vec<Propagator>,
    indices: Vec<i16>,
    numerator: Atom,
    variables: &[Atom],
    template: &IntegralFamily,
    remaining: &mut usize,
    out: &mut Vec<IntegralTerm>,
) -> Result<()> {
    *remaining = remaining
        .checked_sub(1)
        .ok_or_else(|| Error::Limit("boundary partial fraction budget exhausted".into()))?;
    let active = indices
        .iter()
        .enumerate()
        .filter_map(|(i, &n)| (n > 0).then_some(i))
        .collect::<Vec<_>>();
    let transpose = (0..variables.len())
        .map(|j| {
            active
                .iter()
                .map(|&i| denominators[i].scalar_products[j].clone())
                .collect()
        })
        .collect();
    let relations = nullspace(transpose);
    if let Some(relation) = relations.first() {
        let constant = active
            .iter()
            .zip(relation)
            .fold(Atom::new(), |a, (&i, c)| a + c * &denominators[i].constant)
            .together()
            .cancel();
        if !constant.is_zero() {
            for (&i, c) in active.iter().zip(relation) {
                if c.is_zero() {
                    continue;
                }
                let mut next = indices.clone();
                next[i] -= 1;
                partial_fraction(
                    denominators.clone(),
                    next,
                    (&numerator * c / &constant).together().cancel(),
                    variables,
                    template,
                    remaining,
                    out,
                )?;
            }
        } else {
            let pivot = relation.iter().rposition(|c| !c.is_zero()).unwrap();
            for (j, (&i, c)) in active.iter().zip(relation).enumerate() {
                if j == pivot || c.is_zero() {
                    continue;
                }
                let mut next = indices.clone();
                next[i] -= 1;
                next[active[pivot]] = next[active[pivot]]
                    .checked_add(1)
                    .ok_or_else(|| Error::Limit("partial fraction index overflow".into()))?;
                partial_fraction(
                    denominators.clone(),
                    next,
                    (-&numerator * c / &relation[pivot]).together().cancel(),
                    variables,
                    template,
                    remaining,
                    out,
                )?;
            }
        }
        return Ok(());
    }
    let mut propagators = active
        .iter()
        .map(|&i| denominators[i].clone())
        .collect::<Vec<_>>();
    let physical = propagators.len();
    let mut indices = active.iter().map(|&i| indices[i]).collect::<Vec<_>>();
    let n = variables.len();
    for j in 0..n {
        let mut row = vec![Atom::new(); n];
        row[j] = Atom::num(1);
        let mut matrix = propagators
            .iter()
            .map(|p| p.scalar_products.clone())
            .collect::<Vec<_>>();
        matrix.push(row.clone());
        if rref(matrix).1.len() > propagators.len() {
            propagators.push(Propagator {
                constant: Atom::new(),
                scalar_products: row,
            });
            indices.push(0);
        }
    }
    let inverse = inverse(
        &propagators
            .iter()
            .map(|p| p.scalar_products.clone())
            .collect::<Vec<_>>(),
    )?;
    let d = (0..n)
        .map(|j| Atom::var(symbol!(format!("symbolica_amflow::boundary_d_{j}"))))
        .collect::<Vec<_>>();
    let rules = variables
        .iter()
        .enumerate()
        .map(|(i, v)| {
            (
                v.clone(),
                (0..n).fold(Atom::new(), |a, j| {
                    a + &inverse[i][j] * (&d[j] - &propagators[j].constant)
                }),
            )
        })
        .collect::<BTreeMap<_, _>>();
    let numerator = substitute(&numerator, &rules).expand();
    let mut family = template.clone();
    family.propagators = propagators;
    family.physical_propagators = physical;
    for (m, c) in numerator.coefficient_list::<i32>(&d) {
        if c.is_zero() {
            continue;
        }
        let powers = powers(&m, &d)?;
        let integral = Integral(
            indices
                .iter()
                .zip(powers)
                .map(|(&a, b)| {
                    a.checked_sub(b)
                        .ok_or_else(|| Error::Limit("boundary index overflow".into()))
                })
                .collect::<Result<_>>()?,
        );
        out.push(IntegralTerm {
            coefficient: c,
            family: family.clone(),
            integral,
        });
    }
    Ok(())
}

#[derive(Clone, Debug)]
pub struct FactorizedTerm {
    pub coefficient: Atom,
    pub factors: Vec<IntegralTerm>,
}

/// Integrate out mixed hard/soft tensor structure and split the resulting
/// scalar integrands into vacuum and soft families.
pub fn factor_region(
    expression: &Atom,
    coordinates: &[Atom],
    family: &IntegralFamily,
    hard: &[bool],
    budget: usize,
) -> Result<Vec<FactorizedTerm>> {
    let loops = family.loops.len();
    if hard.len() != loops {
        return Err(Error::InvalidInput("hard mask dimensions".into()));
    }
    let mut labels = Vec::new();
    for i in 0..loops {
        for j in i..loops {
            labels.push((i, j));
        }
    }
    for i in 0..loops {
        for e in 0..family.external.len() {
            labels.push((i, loops + e));
        }
    }
    if coordinates.len() != labels.len() {
        return Err(Error::InvalidInput("region coordinate dimensions".into()));
    }
    let is_hard = |i: usize| i < loops && hard[i];
    let mut numerator = Atom::num(1);
    let mut denominator_h = Atom::num(1);
    let mut denominator_s = Atom::num(1);
    let symbols = coordinates
        .iter()
        .map(|a| match a.as_view() {
            AtomView::Var(v) => Ok(v.get_symbol()),
            _ => Err(Error::InvalidInput("coordinate must be a symbol".into())),
        })
        .collect::<Result<Vec<_>>>()?;
    for (base, power) in factors(&expression.together().factor()) {
        if power >= 0 {
            numerator *= base.pow(power);
            continue;
        }
        let present = symbols
            .iter()
            .enumerate()
            .filter_map(|(j, &s)| (!base.derivative(s).is_zero()).then_some(j))
            .collect::<Vec<_>>();
        if present.is_empty() {
            numerator *= base.pow(power);
        } else if present
            .iter()
            .all(|&j| is_hard(labels[j].0) && is_hard(labels[j].1))
        {
            denominator_h *= base.pow(power);
        } else if present
            .iter()
            .all(|&j| !is_hard(labels[j].0) && !is_hard(labels[j].1))
        {
            denominator_s *= base.pow(power);
        } else {
            return Err(Error::Unsupported(
                "mixed hard/soft boundary denominator".into(),
            ));
        }
    }
    let scalar = |i: usize, j: usize| -> Atom {
        let (i, j) = if i <= j { (i, j) } else { (j, i) };
        if i >= loops {
            family.external_gram[i - loops][j - loops].clone()
        } else {
            coordinates[labels.iter().position(|&pair| pair == (i, j)).unwrap()].clone()
        }
    };
    let mut projector = crate::tensor::TensorProjector::new(
        Atom::num(family.dimension) - Atom::num(2) * Atom::var(family.epsilon),
    );
    let mut projected = Atom::new();
    for (m, c) in numerator.expand().coefficient_list::<i32>(coordinates) {
        let mut h = Vec::new();
        let mut s = Vec::new();
        let mut unmixed = c;
        for (j, n) in powers(&m, coordinates)?.into_iter().enumerate() {
            let (a, b) = labels[j];
            if is_hard(a) != is_hard(b) {
                let (a, b) = if is_hard(a) { (a, b) } else { (b, a) };
                if n < 0 {
                    return Err(Error::Unsupported("nonpolynomial mixed numerator".into()));
                }
                h.extend(std::iter::repeat_n(a, n as usize));
                s.extend(std::iter::repeat_n(b, n as usize));
            } else {
                unmixed *= coordinates[j].clone().pow(n as i64);
            }
        }
        let gram = |v: &[usize]| {
            v.iter()
                .map(|&i| v.iter().map(|&j| scalar(i, j)).collect())
                .collect::<Vec<Vec<_>>>()
        };
        projected += unmixed * projector.project(&gram(&h), &gram(&s))?;
    }
    let mut pieces = Vec::new();
    for selected_hard in [true, false] {
        let selected = (0..loops)
            .filter(|&j| hard[j] == selected_hard)
            .collect::<Vec<_>>();
        let vars = labels
            .iter()
            .enumerate()
            .filter_map(|(j, &(a, b))| {
                (selected.contains(&a) && (selected.contains(&b) || (!selected_hard && b >= loops)))
                    .then_some(coordinates[j].clone())
            })
            .collect::<Vec<_>>();
        let template = IntegralFamily {
            name: "amflow_boundary".into(),
            loops: selected.iter().map(|&j| family.loops[j].clone()).collect(),
            external: if selected_hard {
                vec![]
            } else {
                family.external.clone()
            },
            external_gram: if selected_hard {
                vec![]
            } else {
                family.external_gram.clone()
            },
            propagators: vec![],
            physical_propagators: 0,
            epsilon: family.epsilon,
            dimension: family.dimension,
        };
        pieces.push((vars, template));
    }
    let mut out = Vec::new();
    for (m, c) in projected.expand().coefficient_list::<i32>(coordinates) {
        let mut expressions = [denominator_h.clone(), denominator_s.clone()];
        for (j, n) in powers(&m, coordinates)?.into_iter().enumerate() {
            if n == 0 {
                continue;
            }
            let side = usize::from(!is_hard(labels[j].0));
            expressions[side] *= coordinates[j].clone().pow(n as i64);
        }
        let mut products = vec![FactorizedTerm {
            coefficient: c,
            factors: vec![],
        }];
        for (side, (vars, template)) in pieces.iter().enumerate() {
            if template.loops.is_empty() {
                for p in &mut products {
                    p.coefficient *= &expressions[side];
                }
                continue;
            }
            let terms = to_integrals(&expressions[side], vars, template, budget)?;
            let mut next = Vec::new();
            for product in &products {
                for term in &terms {
                    let mut p = product.clone();
                    p.coefficient *= &term.coefficient;
                    let mut term = term.clone();
                    term.coefficient = Atom::num(1);
                    p.factors.push(term);
                    next.push(p);
                }
            }
            products = next;
        }
        out.extend(products);
        if out.len() > budget {
            return Err(Error::Limit(
                "factorized boundary term budget exhausted".into(),
            ));
        }
    }
    Ok(out)
}
