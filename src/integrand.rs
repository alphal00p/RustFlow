//! Convert rational boundary integrands with exact scalar coefficients to independent families.
//! Gaussian-rational coefficients use HEPKit; other admitted expressions use the exact AtomField route.
use crate::algebra::{inverse, nullspace, rref};
use crate::family::substitute;
use crate::*;
use std::collections::BTreeMap;
use symbolica::prelude::*;

mod native;

#[derive(Clone, Debug)]
pub struct IntegralTerm {
    pub coefficient: Atom,
    pub family: IntegralFamily,
    pub integral: Integral,
}

pub(crate) use crate::coefficient::{factors, powers};

fn affine(a: &Atom, variables: &[Atom]) -> Result<Propagator> {
    let mut out = Propagator {
        constant: Atom::new(),
        scalar_products: vec![Atom::new(); variables.len()],
    };
    for (m, c) in crate::coefficient::exact_coefficient_list(a, variables)? {
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
    native::validate_coordinates(variables, template)?;
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
    if native::supports(&denominators) {
        return native::to_integrals(
            &denominators,
            &indices,
            &numerator,
            variables,
            template,
            budget,
        );
    }
    // Keep the existing exact AtomField route for coefficients outside the
    // structurally admitted native Gaussian-rational field. Native errors and
    // limits never enter this route; i is handled by the native exact field.
    let mut out = Vec::new();
    let mut remaining = budget;
    atom_field_partial_fraction(
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

fn atom_field_partial_fraction(
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
                atom_field_partial_fraction(
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
                atom_field_partial_fraction(
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
    let d = denominator_labels(&numerator, &denominators, variables)?;
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
    append_numerator_terms(&numerator, &d, &family, &indices, out)
}

fn denominator_labels(
    numerator: &Atom,
    denominators: &[Propagator],
    variables: &[Atom],
) -> Result<Vec<Atom>> {
    let sources = boundary_sources(numerator, denominators, variables);
    (0..variables.len())
        .map(|j| fresh_boundary_symbol(&format!("symbolica_amflow::boundary_d_{j}"), &sources))
        .collect()
}

fn boundary_sources<'a>(
    numerator: &'a Atom,
    denominators: &'a [Propagator],
    variables: &'a [Atom],
) -> Vec<&'a Atom> {
    std::iter::once(numerator)
        .chain(variables)
        .chain(
            denominators
                .iter()
                .flat_map(|d| std::iter::once(&d.constant).chain(&d.scalar_products)),
        )
        .collect()
}

// Internal notation must not turn a user's scalar parameter into a momentum
// or denominator coordinate. All callers use distinct stems.
fn fresh_boundary_symbol(stem: &str, sources: &[&Atom]) -> Result<Atom> {
    let mut suffix = 0_usize;
    loop {
        let label = Atom::var(symbol!(if suffix == 0 {
            stem.to_owned()
        } else {
            format!("{stem}_fresh_{suffix}")
        }));
        if !sources.iter().any(|a| a.contains(label.as_view())) {
            return Ok(label);
        }
        suffix = suffix
            .checked_add(1)
            .ok_or_else(|| Error::Limit("boundary symbol allocation exhausted".into()))?;
    }
}

fn append_numerator_terms(
    numerator: &Atom,
    d: &[Atom],
    family: &IntegralFamily,
    indices: &[i16],
    out: &mut Vec<IntegralTerm>,
) -> Result<()> {
    for (m, c) in crate::coefficient::exact_coefficient_list(numerator, d)? {
        if c.is_zero() {
            continue;
        }
        let powers = powers(&m, d)?;
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

/// Scalar coordinates for one side of a tensor-projected region.
/// Indices refer to the supplied region coordinates and loop routing, before
/// selecting hard or soft loops. Region Jacobians and measure/distribution
/// metadata remain with the caller; projection does not change either measure.
#[derive(Clone, Debug)]
pub struct RegionFactorSpace {
    pub coordinates: Vec<Atom>,
    pub source_coordinate_indices: Vec<usize>,
    pub source_loop_indices: Vec<usize>,
    /// Loop/external space for conversion, with no denominator slots yet.
    pub template: IntegralFamily,
}

/// A separated scalar term after hard tensor projection and before integral
/// conversion. In particular, soft denominators are preserved as expressions.
#[derive(Clone, Debug)]
pub struct ProjectedRegionTerm {
    pub hard: Atom,
    pub soft: Atom,
}

#[derive(Clone, Debug)]
pub struct ProjectedRegion {
    pub hard: RegionFactorSpace,
    pub soft: RegionFactorSpace,
    pub terms: Vec<ProjectedRegionTerm>,
}

/// Integrate out mixed hard/soft tensor structure and split the resulting
/// scalar expressions without choosing an integration measure for the soft
/// factors. No partial fractions, scaleless-sector tests, or ordinary integral
/// conversion are performed here. A weighted-measure owner must integrate or
/// recursively reduce the retained soft expression, including its denominators.
/// `budget` bounds the number of separated expression products.
pub fn projected_factor_region(
    expression: &Atom,
    coordinates: &[Atom],
    family: &IntegralFamily,
    hard: &[bool],
    budget: usize,
) -> Result<ProjectedRegion> {
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
    for (m, c) in crate::coefficient::exact_coefficient_list(&numerator.expand(), coordinates)? {
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
        let source_coordinate_indices = labels
            .iter()
            .enumerate()
            .filter_map(|(j, &(a, b))| {
                (selected.contains(&a) && (selected.contains(&b) || (!selected_hard && b >= loops)))
                    .then_some(j)
            })
            .collect::<Vec<_>>();
        let vars = source_coordinate_indices
            .iter()
            .map(|&j| coordinates[j].clone())
            .collect();
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
        pieces.push(RegionFactorSpace {
            coordinates: vars,
            source_coordinate_indices,
            source_loop_indices: selected,
            template,
        });
    }
    // Keep the larger polynomial intact while collecting monomials in the
    // smaller coordinate space. Expanding both sides first repeats denominator
    // conversion and misses cancellations in denominator coordinates. All-hard
    // and all-soft regions each need only one conversion of the full polynomial.
    let split_side = usize::from(pieces[0].coordinates.len() > pieces[1].coordinates.len());
    let mut terms = Vec::new();
    for (monomial, polynomial) in crate::coefficient::exact_coefficient_list(
        &projected.expand(),
        &pieces[split_side].coordinates,
    )? {
        let mut expressions = [denominator_h.clone(), denominator_s.clone()];
        expressions[split_side] *= monomial;
        expressions[1 - split_side] *= polynomial;
        if terms.len() == budget {
            return Err(Error::Limit(format!(
                "projected boundary term budget exhausted: more than {budget} expression products"
            )));
        }
        let [hard, soft] = expressions;
        terms.push(ProjectedRegionTerm { hard, soft });
    }
    let mut pieces = pieces.into_iter();
    Ok(ProjectedRegion {
        hard: pieces.next().expect("hard coordinate space"),
        soft: pieces.next().expect("soft coordinate space"),
        terms,
    })
}

/// Integrate out mixed hard/soft tensor structure and convert the projected
/// expressions to ordinary vacuum and soft families. Equivalent integral
/// products are combined exactly before applying `budget` to the output size.
/// The same budget bounds each partial-fraction decomposition.
pub fn factor_region(
    expression: &Atom,
    coordinates: &[Atom],
    family: &IntegralFamily,
    hard: &[bool],
    budget: usize,
) -> Result<Vec<FactorizedTerm>> {
    // Ordinary callers historically bound combined integral products, which
    // need not have the same count as intermediate projected expressions.
    let projected = projected_factor_region(expression, coordinates, family, hard, usize::MAX)?;
    let pieces = [&projected.hard, &projected.soft];
    let mut combined = BTreeMap::<Vec<Vec<(Vec<Atom>, i64)>>, FactorizedTerm>::new();
    for term in &projected.terms {
        let expressions = [&term.hard, &term.soft];
        let mut products = vec![FactorizedTerm {
            coefficient: Atom::num(1),
            factors: vec![],
        }];
        for (side, space) in pieces.iter().enumerate() {
            if space.template.loops.is_empty() {
                for p in &mut products {
                    p.coefficient *= expressions[side];
                }
                continue;
            }
            let terms = to_integrals(
                expressions[side],
                &space.coordinates,
                &space.template,
                budget,
            )?;
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
        for product in products {
            // All families on a side use the same loop/external coordinates.
            // Ignore unused completion slots and denominator ordering so that
            // partial fractions with different active sets share an exact key.
            let key = product
                .factors
                .iter()
                .map(|factor| {
                    let mut powers = BTreeMap::<Vec<Atom>, i64>::new();
                    for (d, &n) in factor.family.propagators.iter().zip(&factor.integral.0) {
                        if n != 0 {
                            let shape = std::iter::once(d.constant.clone())
                                .chain(d.scalar_products.iter().cloned())
                                .collect();
                            *powers.entry(shape).or_default() += i64::from(n);
                        }
                    }
                    powers.into_iter().filter(|(_, n)| *n != 0).collect()
                })
                .collect();
            match combined.entry(key) {
                std::collections::btree_map::Entry::Occupied(mut entry) => {
                    entry.get_mut().coefficient += product.coefficient;
                }
                std::collections::btree_map::Entry::Vacant(entry) => {
                    entry.insert(product);
                }
            }
        }
    }
    let out = combined
        .into_values()
        .filter_map(|mut term| {
            term.coefficient = term.coefficient.together().cancel();
            (!term.coefficient.is_zero()).then_some(term)
        })
        .collect::<Vec<_>>();
    if out.len() > budget {
        return Err(Error::Limit(format!(
            "factorized boundary term budget exhausted: {} distinct terms exceed {budget}",
            out.len()
        )));
    }
    Ok(out)
}
