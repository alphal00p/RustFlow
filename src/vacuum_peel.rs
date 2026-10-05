//! Generic single-mass vacuum reduction by homogeneity and one radial integral.
use crate::{ComplexFloat as C, Error, Integral, IntegralFamily, Precision, Result};
use feynkit_graph::IntegralFamily as NativeFamily;
use feynkit_kinematics::Kinematics;
use symbolica::prelude::*;

#[derive(Clone, Debug)]
pub(crate) struct SingleMassPlan {
    pub child: IntegralFamily,
    pub target: Integral,
    #[cfg(test)]
    pub massive_line: usize,
    #[cfg(test)]
    pub transformation: Vec<Vec<Atom>>,
    pub determinant: Rational,
    pub mass_squared: Rational,
    pub massive_scale: Rational,
    pub massive_power: u32,
    pub remaining_power: i64,
    pub parent_loops: usize,
}

fn rational(a: &Atom) -> Option<Rational> {
    let AtomView::Num(n) = a.as_view() else {
        return None;
    };
    let symbolica::coefficient::Coefficient::Complex(c) = n.get_coeff_view().to_owned() else {
        return None;
    };
    c.im.is_zero().then_some(c.re)
}

impl SingleMassPlan {
    /// A topology-independent radial integral, multiplied by the exact routing
    /// Jacobian. The child includes the signs of every remaining propagator.
    pub(crate) fn prefactor(&self, epsilon: &Rational, p: Precision) -> Result<C> {
        let dimension = Rational::from(self.child.dimension) - epsilon * &Rational::from(2);
        let half = &dimension / &Rational::from(2);
        let m = Rational::from(self.remaining_power)
            - Rational::from((self.parent_loops - 1) as i64) * &half;
        let a = Rational::from(self.massive_power as i64);
        let gamma = |x: &Rational| p.gamma_real(&p.rational(x).re);
        let ratio = p.div(
            &p.mul(&gamma(&(&half - &m))?, &gamma(&(&a + &m - &half))?),
            &p.mul(&gamma(&half)?, &gamma(&a)?),
        );
        let mass = p.pow(
            &p.rational(&self.mass_squared),
            &p.rational(&(&half - &m - &a)),
        );
        let scale = p.powi(
            &p.rational(&self.massive_scale),
            -i64::from(self.massive_power),
        );
        let jacobian = p.pow(
            &p.rational(&self.determinant.clone().abs()),
            &p.rational(&dimension),
        );
        let value = p.mul(&p.mul(&ratio, &mass), &p.mul(&scale, &jacobian));
        let value = if self.massive_power.is_multiple_of(2) {
            value
        } else {
            p.neg(&value)
        };
        if !p.finite(&value) {
            return Err(Error::Numerical(
                "nonfinite single-mass radial prefactor".into(),
            ));
        }
        Ok(value)
    }

    /// None preserves the caller's existing fallback. The caller must continue
    /// checking custom terminal providers and certified scaleless sectors first.
    pub(crate) fn find(family: &IntegralFamily, target: &Integral) -> Result<Option<Self>> {
        family.validate_integral(target)?;
        let loops = family.loops.len();
        if loops < 2 || !family.external.is_empty() || target.0.iter().any(|&n| n < 0) {
            return Ok(None);
        }
        family.validate()?;
        let active = target
            .0
            .iter()
            .enumerate()
            .filter_map(|(i, &n)| (n > 0).then_some(i))
            .collect::<Vec<_>>();
        let massive = active
            .iter()
            .copied()
            .filter(|&i| !family.propagators[i].constant.is_zero())
            .collect::<Vec<_>>();
        if massive.len() != 1 {
            return Ok(None);
        }
        let massive_line = massive[0];
        let mut branches = Vec::new();
        for &i in &active {
            let d = &family.propagators[i];
            if std::iter::once(&d.constant)
                .chain(&d.scalar_products)
                .any(|a| rational(a).is_none())
            {
                return Ok(None);
            }
            let branch = match crate::regions::branch(d, loops) {
                Ok(b) => b,
                Err(Error::Unsupported(_)) => return Ok(None),
                Err(e) => return Err(e),
            };
            // A positive quadratic normalization preserves the declared +i0
            // when its scalar factor is extracted or Wick rotated.
            let mut offset = 0;
            let positive = (0..loops).find_map(|row| {
                let value = rational(&d.scalar_products[offset]).unwrap();
                offset += loops - row;
                (!value.is_zero()).then_some(value > Rational::zero())
            });
            if positive != Some(true) {
                return Ok(None);
            }
            branches.push((i, branch));
        }
        let mut basis = vec![
            branches
                .iter()
                .find(|(i, _)| *i == massive_line)
                .unwrap()
                .1
                .clone(),
        ];
        for (_, branch) in &branches {
            if basis.len() == loops {
                break;
            }
            let mut candidate = basis.clone();
            candidate.push(branch.clone());
            if crate::algebra::rref(candidate.clone()).1.len() > basis.len() {
                basis = candidate;
            }
        }
        if basis.len() != loops {
            return Ok(None);
        }
        let transformation = crate::algebra::inverse(&basis)?;
        let determinant = rational(&crate::algebra::determinant(transformation.clone()))
            .ok_or_else(|| Error::Numerical("nonrational exact single-mass routing".into()))?;
        let transformed = family.transform_loops(&transformation)?;
        let massive_scale =
            rational(&transformed.propagators[massive_line].scalar_products[0]).unwrap();
        let mass_squared =
            -rational(&transformed.propagators[massive_line].constant).unwrap() / &massive_scale;
        if mass_squared <= Rational::zero() {
            return Ok(None);
        }

        let dimension_symbol = (0..)
            .map(|i| symbol!(format!("symbolica_amflow::vacuum_peel_dimension_{i}")))
            .find(|&d| d != family.epsilon)
            .unwrap();
        let external = symbol!("symbolica_amflow::vacuum_peel_external").call(0);
        let child_loops = (0..loops - 1)
            .map(|i| symbol!("symbolica_amflow::vacuum_peel_loop").call(i as i64))
            .collect::<Vec<_>>();
        let kinematics = Kinematics::in_dimension(&Atom::var(dimension_symbol))
            .map_err(|e| Error::InvalidInput(e.to_string()))?
            .with_momenta(
                child_loops
                    .iter()
                    .cloned()
                    .chain(std::iter::once(external.clone())),
            )
            .and_then(|k| k.with_mass_squared(&external, Atom::num(-1)))
            .map_err(|e| Error::InvalidInput(e.to_string()))?;
        // q is exactly the first transformed loop. Native Kinematics owns
        // q²=-1 and scalar-product Atom identities in the child family.
        let momenta = std::iter::once(external.clone())
            .chain(child_loops.iter().cloned())
            .collect::<Vec<_>>();
        let mut products = Vec::new();
        for i in 0..loops {
            for j in i..loops {
                products.push(
                    kinematics
                        .scalar_product(&momenta[i], &momenta[j])
                        .map_err(|e| Error::InvalidInput(e.to_string()))?,
                );
            }
        }
        let child_indices = active
            .iter()
            .copied()
            .filter(|&i| i != massive_line)
            .collect::<Vec<_>>();
        let denominators = child_indices
            .iter()
            .map(|&i| {
                let d = &transformed.propagators[i];
                d.scalar_products
                    .iter()
                    .zip(&products)
                    .fold(d.constant.clone(), |sum, (c, product)| sum + c * product)
                    .expand()
            })
            .collect::<Vec<_>>();
        let native = NativeFamily::new(child_loops, vec![external], denominators, &kinematics)
            .map_err(|e| Error::InvalidInput(e.to_string()))?;
        if !native.is_independent() {
            return Ok(None);
        }
        let native = native
            .complete(&[])
            .map_err(|e| Error::InvalidInput(e.to_string()))?;
        let mut child = IntegralFamily::from_hepkit(
            &native,
            child_indices.len(),
            family.epsilon,
            family.dimension,
        )?;
        child.name = format!("{}_single_mass_child", family.name);
        let mut child_powers = child_indices
            .iter()
            .map(|&i| target.0[i])
            .collect::<Vec<_>>();
        let remaining_power = child_powers.iter().map(|&n| i64::from(n)).sum();
        child_powers.resize(child.propagators.len(), 0);
        let target_child = Integral(child_powers);
        child.validate_integral(&target_child)?;
        Ok(Some(Self {
            child,
            target: target_child,
            #[cfg(test)]
            massive_line,
            #[cfg(test)]
            transformation,
            determinant,
            mass_squared,
            massive_scale,
            massive_power: target.0[massive_line] as u32,
            remaining_power,
            parent_loops: loops,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::SingleMassPlan;
    use crate::{ComplexFloat, Integral, IntegralFamily, Precision, Propagator, vacuum};
    use symbolica::prelude::*;
    fn six_line(mass: i64) -> IntegralFamily {
        IntegralFamily {
            name: "connected_six_line_vacuum".into(),
            loops: vec!["q".into(), "l".into(), "k".into()],
            external: vec![],
            external_gram: vec![],
            propagators: [
                ([1, 0, 0], mass),
                ([0, 1, 0], 0),
                ([-1, 1, 0], 0),
                ([0, 0, 1], 0),
                ([-1, 0, 1], 0),
                ([0, 1, -1], 0),
            ]
            .into_iter()
            .map(|(v, m)| Propagator::quadratic(&v, &[], Atom::num(m), &[]).unwrap())
            .collect(),
            physical_propagators: 6,
            epsilon: symbol!("single_mass_test::epsilon"),
            dimension: 4,
        }
    }
    fn sunset(mass: Rational) -> IntegralFamily {
        IntegralFamily {
            name: "sunset".into(),
            loops: vec!["q".into(), "l".into()],
            external: vec![],
            external_gram: vec![],
            propagators: [[1, 0], [0, 1], [-1, 1]]
                .into_iter()
                .enumerate()
                .map(|(i, v)| {
                    Propagator::quadratic(
                        &v,
                        &[],
                        if i == 0 {
                            Atom::num(mass.clone())
                        } else {
                            Atom::new()
                        },
                        &[],
                    )
                    .unwrap()
                })
                .collect(),
            physical_propagators: 3,
            epsilon: symbol!("single_mass_test::epsilon"),
            dimension: 4,
        }
    }
    fn exact(a: &Atom, b: &Atom) {
        assert!((a - b).together().cancel().is_zero(), "{a} != {b}");
    }
    fn bubble_value(b: u32, c: u32, dimension: &Rational, p: Precision) -> ComplexFloat {
        let h = dimension / &Rational::from(2);
        let b = Rational::from(b as i64);
        let c = Rational::from(c as i64);
        let gamma = |x: &Rational| p.gamma_real(&p.rational(x).re).unwrap();
        let numerator = p.mul(
            &gamma(&(&b + &c - &h)),
            &p.mul(&gamma(&(&h - &b)), &gamma(&(&h - &c))),
        );
        let denominator = p.mul(
            &gamma(&b),
            &p.mul(&gamma(&c), &gamma(&(dimension - &b - &c))),
        );
        let sign = if ((&b + &c).to_string().parse::<i64>().unwrap()) % 2 == 0 {
            1
        } else {
            -1
        };
        p.mul(&p.i(sign), &p.div(&numerator, &denominator))
    }

    #[test]
    fn connected_six_line_peels_to_complete_two_loop_five_line_family() {
        let family = six_line(3);
        let plan = SingleMassPlan::find(&family, &Integral(vec![1; 6]))
            .unwrap()
            .unwrap();
        assert_eq!(plan.parent_loops, 3);
        assert_eq!(plan.child.loops.len(), 2);
        assert_eq!(plan.child.external.len(), 1);
        assert_eq!(plan.child.propagators.len(), 5);
        assert_eq!(plan.child.physical_propagators, 5);
        assert_eq!(plan.target, Integral(vec![1; 5]));
        assert_eq!(plan.remaining_power, 5);
        assert_eq!(plan.mass_squared, Rational::from(3));
        exact(&plan.child.external_gram[0][0], &Atom::num(-1));
        let expected = [
            (0, [1, 0, 0, 0, 0]),
            (-1, [1, 0, 0, -2, 0]),
            (0, [0, 0, 1, 0, 0]),
            (-1, [0, 0, 1, 0, -2]),
            (0, [1, -2, 1, 0, 0]),
        ];
        for (d, (constant, coefficients)) in plan.child.propagators.iter().zip(expected) {
            exact(&d.constant, &Atom::num(constant));
            for (a, b) in d.scalar_products.iter().zip(coefficients) {
                exact(a, &Atom::num(b));
            }
        }
    }

    #[test]
    fn native_completion_fills_only_zero_power_child_isp_slots() {
        let mut family = six_line(3);
        family.propagators = [
            ([1, 0, 0], 3),
            ([0, 1, 0], 0),
            ([0, 0, 1], 0),
            ([1, 1, 1], 0),
        ]
        .into_iter()
        .map(|(v, m)| Propagator::quadratic(&v, &[], Atom::num(m), &[]).unwrap())
        .collect();
        for coordinate in [1, 2] {
            let mut scalar_products = vec![Atom::new(); 6];
            scalar_products[coordinate] = Atom::num(1);
            family.propagators.push(Propagator {
                constant: Atom::new(),
                scalar_products,
            });
        }
        family.physical_propagators = 4;
        let plan = SingleMassPlan::find(&family, &Integral(vec![1, 1, 1, 1, 0, 0]))
            .unwrap()
            .unwrap();
        assert_eq!(plan.child.loops.len(), 2);
        assert_eq!(plan.child.physical_propagators, 3);
        assert_eq!(plan.child.propagators.len(), 5);
        assert_eq!(plan.target, Integral(vec![1, 1, 1, 0, 0]));
        plan.child.validate().unwrap();
    }

    #[test]
    fn radial_factor_matches_existing_sunset_with_raised_powers_and_general_dimension() {
        let p = Precision::decimal(70).unwrap();
        let eps = Rational::from((1, 10));
        for d0 in [4, 6] {
            let dimension = Rational::from(d0) - &eps * &Rational::from(2);
            for powers in [[1, 1, 1], [2, 1, 1], [1, 2, 1], [2, 3, 1]] {
                let mut family = sunset(Rational::from((3, 2)));
                family.dimension = d0;
                let plan =
                    SingleMassPlan::find(&family, &Integral(powers.map(|x| x as i16).to_vec()))
                        .unwrap()
                        .unwrap();
                let value = p.mul(
                    &plan.prefactor(&eps, p).unwrap(),
                    &bubble_value(powers[1], powers[2], &dimension, p),
                );
                let expected = vacuum::single_mass_sunset(
                    powers,
                    &p.rational(&Rational::from((3, 2))),
                    &dimension,
                    p,
                )
                .unwrap();
                assert!(
                    p.close(&value, &expected, 50),
                    "{powers:?} D0={d0}: {value} != {expected}"
                );
            }
        }
    }

    #[test]
    fn six_line_mass_scaling_and_raised_massive_index_follow_exact_homogeneity() {
        let p = Precision::decimal(70).unwrap();
        let base = SingleMassPlan::find(&six_line(1), &Integral(vec![1; 6]))
            .unwrap()
            .unwrap();
        let mass3 = SingleMassPlan::find(&six_line(3), &Integral(vec![1; 6]))
            .unwrap()
            .unwrap();
        let raised = SingleMassPlan::find(&six_line(3), &Integral(vec![2, 1, 1, 1, 1, 1]))
            .unwrap()
            .unwrap();
        assert_eq!(
            mass3.child.integral_key(&mass3.target).unwrap(),
            raised.child.integral_key(&raised.target).unwrap()
        );
        for eps in [Rational::from((1, 10)), Rational::from((1, 11))] {
            let a = base.prefactor(&eps, p).unwrap();
            let b = mass3.prefactor(&eps, p).unwrap();
            let expected = p.mul(
                &a,
                &p.pow(&p.i(3), &p.rational(&(-Rational::from(3) * &eps))),
            );
            assert!(p.close(&b, &expected, 50));
            assert!(p.close(
                &p.div(&raised.prefactor(&eps, p).unwrap(), &b),
                &p.rational(&(-eps)),
                50
            ));
        }
    }

    #[test]
    fn denominator_order_scaling_and_nonunit_routing_preserve_the_integral() {
        let p = Precision::decimal(70).unwrap();
        let eps = Rational::from((1, 10));
        let dimension = Rational::from(4) - &eps * &Rational::from(2);
        let original = sunset(Rational::from(3));
        let transformation = vec![
            vec![Atom::num(2), Atom::num(1)],
            vec![Atom::num(1), Atom::num(3)],
        ];
        let transformed = original.transform_loops(&transformation).unwrap();
        let weights = [2, 3, 5];
        let mut modified = transformed.clone();
        for (d, c) in modified.propagators.iter_mut().zip(weights) {
            d.constant = &d.constant * Atom::num(c);
            for a in &mut d.scalar_products {
                *a = &*a * Atom::num(c);
            }
        }
        let order = [2, 0, 1];
        modified.propagators = order.map(|i| modified.propagators[i].clone()).to_vec();
        let plan = SingleMassPlan::find(&modified, &Integral(vec![1; 3]))
            .unwrap()
            .unwrap();
        assert_eq!(plan.massive_line, 1);
        assert_ne!(plan.determinant.clone().abs(), Rational::one());
        let routed = modified.transform_loops(&plan.transformation).unwrap();
        let massive = &routed.propagators[plan.massive_line];
        exact(
            &massive.scalar_products[0],
            &Atom::num(plan.massive_scale.clone()),
        );
        assert!(massive.scalar_products[1..].iter().all(Atom::is_zero));
        exact(
            &massive.constant,
            &Atom::num(-plan.massive_scale.clone() * &plan.mass_squared),
        );
        // Reconstruct every active denominator in native child coordinates. The
        // remaining two are a scaled ordinary bubble after this chosen routing.
        let child = &plan.child;
        let mut scales = Vec::new();
        for d in &child.propagators {
            scales.push(match d.scalar_products[0].as_view() {
                AtomView::Num(n) => match n.get_coeff_view().to_owned() {
                    symbolica::coefficient::Coefficient::Complex(c) => c.re,
                    _ => panic!(),
                },
                _ => panic!(),
            });
        }
        let child_value = p.div(
            &bubble_value(1, 1, &dimension, p),
            &p.rational(&(&scales[0] * &scales[1])),
        );
        // This child has q²=-1 but its two branch momenta can differ by a multiple
        // of q. Derive that exact transfer square from each affine quadratic.
        let shifts = child
            .propagators
            .iter()
            .zip(&scales)
            .map(|(d, c)| {
                (&d.scalar_products[1] / Atom::num(c.clone()) / Atom::num(2))
                    .together()
                    .cancel()
            })
            .collect::<Vec<_>>();
        for ((d, scale), shift) in child.propagators.iter().zip(&scales).zip(&shifts) {
            exact(&d.constant, &(-Atom::num(scale.clone()) * shift.pow(2)));
        }
        let transfer = (&shifts[0] - &shifts[1]).pow(2).together().cancel();
        let transfer = p.eval(&transfer, &Default::default()).unwrap();
        let child_value = p.mul(
            &child_value,
            &p.pow(
                &transfer,
                &p.rational(&(&dimension / &Rational::from(2) - Rational::from(2))),
            ),
        );
        let value = p.mul(&plan.prefactor(&eps, p).unwrap(), &child_value);
        let expected = vacuum::single_mass_sunset([1, 1, 1], &p.i(3), &dimension, p).unwrap();
        let expected = p.div(
            &expected,
            &p.mul(&p.i(30), &p.pow(&p.i(5), &p.rational(&dimension))),
        );
        assert!(p.close(&value, &expected, 45), "{value} != {expected}");
    }

    #[test]
    fn unsupported_domains_leave_the_existing_fallback_available() {
        let mut family = six_line(3);
        assert!(
            SingleMassPlan::find(&family, &Integral(vec![1, 1, 1, 1, 1, -1]))
                .unwrap()
                .is_none()
        );
        family.propagators[1].constant = Atom::num(-2);
        assert!(
            SingleMassPlan::find(&family, &Integral(vec![1; 6]))
                .unwrap()
                .is_none()
        );
        family = six_line(-3);
        assert!(
            SingleMassPlan::find(&family, &Integral(vec![1; 6]))
                .unwrap()
                .is_none()
        );
        family = six_line(3);
        family.propagators[1].scalar_products[0] = Atom::num(1);
        assert!(
            SingleMassPlan::find(&family, &Integral(vec![1; 6]))
                .unwrap()
                .is_none()
        );
        family = six_line(3);
        family.propagators[1]
            .scalar_products
            .iter_mut()
            .for_each(|a| *a = -a.clone());
        assert!(
            SingleMassPlan::find(&family, &Integral(vec![1; 6]))
                .unwrap()
                .is_none()
        );
    }
}
