//! Uniform polynomial seeds for one occupied shell, with Euclidean spatial
//! measure d^d q/(2*pi)^d. This owner never integrates a surviving virtual pole.
//! Incomplete beta series are summed with a geometric tail bound and an explicit
//! work limit. The caller still verifies working-precision stability separately.
use super::measure::EnergyResidue;
use crate::coefficient::{exact_coefficient_list, powers};
use crate::{ComplexFloat as C, Error, Precision, Result};
use ahash::HashMap;
use symbolica::prelude::*;

#[derive(Clone, Debug)]
pub struct CompactShell {
    pub mass_squared: Rational,
    /// Nonnegative occupied energy bound after fixing charge orientation.
    pub chemical_potential: Rational,
}

impl CompactShell {
    fn validate(&self) -> Result<()> {
        if self.mass_squared < 0 || self.chemical_potential < 0 {
            return Err(Error::InvalidInput(
                "compact seed requires nonnegative mass squared and oriented chemical potential"
                    .into(),
            ));
        }
        Ok(())
    }

    fn angular(dimension: &Rational, p: Precision) -> Result<C> {
        let pi = C::new(p.real(1).pi(), p.real(0));
        let half_d = dimension.clone() / Rational::from(2);
        Ok(p.div(
            &p.scale(&p.pow(&pi, &p.rational(&half_d)), 2, 1),
            &p.mul(
                &p.pow(&p.scale(&pi, 2, 1), &p.rational(dimension)),
                &p.gamma_real(&p.rational(&half_d).re)?,
            ),
        ))
    }

    /// Positive shell moment integral theta(mu-E) E^r (q^2)^k/(2E).
    /// `spatial_dimension` is kept exact; r may be negative for raised-shell
    /// bulk terms. This is a convergent compact integral, not dimensional
    /// continuation through an infrared divergence.
    pub fn moment(
        &self,
        spatial_dimension: &Rational,
        energy_power: i32,
        radial_power: u16,
        digits: u32,
        max_terms: usize,
        p: Precision,
    ) -> Result<C> {
        self.validate()?;
        let a = spatial_dimension.clone() / Rational::from(2) + Rational::from(radial_power);
        if a <= 0 {
            return Err(Error::Unsupported(
                "compact radial integral requires Re(d/2+k)>0".into(),
            ));
        }
        let gap = &self.chemical_potential * &self.chemical_potential - &self.mass_squared;
        if gap < 0 || gap.is_zero() && !self.mass_squared.is_zero() {
            return Ok(p.zero());
        }
        let b = Rational::from((i64::from(energy_power) - 1, 2));
        let angular = Self::angular(spatial_dimension, p)?;
        if self.mass_squared.is_zero() {
            let exponent = &a + &b;
            if exponent <= 0 {
                return Err(Error::Unsupported(format!(
                    "massless compact seed has an infrared divergent radial exponent {}; requires dimensional continuation before mass removal",
                    exponent * Rational::from(2) - Rational::from(1)
                )));
            }
            if gap.is_zero() {
                return Ok(p.zero());
            }
            return Ok(p.scale(
                &p.mul(
                    &angular,
                    &p.div(
                        &p.pow(&p.rational(&gap), &p.rational(&exponent)),
                        &p.rational(&exponent),
                    ),
                ),
                1,
                4,
            ));
        }
        // y=q^2 and t=y/(m^2+y): integral is
        // (m^2)^(a+b) z^a/a 2F1(a,a+b+1;a+1;z).
        let z = gap.clone() / (&gap + &self.mass_squared);
        let second = &a + &b + Rational::from(1);
        let deviation = (&a + &b).abs();
        let mut term = p.i(1);
        let mut sum = term.clone();
        let tolerance = p.tolerance(digits);
        for n in 0..max_terms {
            let n = i64::try_from(n)
                .map_err(|_| Error::Limit("compact summation index overflow".into()))?;
            let ratio =
                z.clone() * (a.clone() + Rational::from(n)) * (second.clone() + Rational::from(n))
                    / ((a.clone() + Rational::from(n + 1)) * Rational::from(n + 1));
            let next = p.mul(&term, &p.rational(&ratio));
            // For every later j>=n, |ratio_j| <= z*(1+|a+b|/(n+1)).
            let bound = z.clone() * (Rational::from(1) + deviation.clone() / Rational::from(n + 1));
            if bound < 1 {
                let tail = p.div(&next, &p.rational(&(Rational::from(1) - bound)));
                if p.norm(&tail) <= tolerance.clone() * p.norm(&sum) {
                    let prefactor = p.div(
                        &p.mul(
                            &p.pow(&p.rational(&self.mass_squared), &p.rational(&(&a + &b))),
                            &p.pow(&p.rational(&z), &p.rational(&a)),
                        ),
                        &p.rational(&a),
                    );
                    return Ok(p.scale(&p.mul(&angular, &p.mul(&prefactor, &sum)), 1, 4));
                }
            }
            sum = p.add(&sum, &next);
            term = next;
        }
        Err(Error::Limit(format!(
            "compact beta series did not meet {digits} digits within {max_terms} terms"
        )))
    }

    /// Occupied correction for an independently raised Euclidean line with
    /// original numerator (q0/i)^r (spatial q^2)^k on the continued Euclidean
    /// shell q0=iE. A Euclidean medium numerator q0^r therefore multiplies this
    /// seed by i^r; it must not be silently identified with E^r. The upper
    /// surface derivatives are integrated exactly, without finite differences.
    pub fn raised_moment(
        &self,
        spatial_dimension: &Rational,
        power: u16,
        energy_power: u16,
        radial_power: u16,
        digits: u32,
        max_terms: usize,
        p: Precision,
    ) -> Result<C> {
        self.validate()?;
        if power == 0 {
            return Err(Error::InvalidInput(
                "occupied line power must be positive".into(),
            ));
        }
        let gap = &self.chemical_potential * &self.chemical_potential - &self.mass_squared;
        if gap < 0 {
            return Ok(p.zero());
        }
        if gap.is_zero() {
            if self.mass_squared.is_zero() {
                let radial_degree = spatial_dimension.clone()
                    + Rational::from(energy_power)
                    + Rational::from(2 * u32::from(radial_power))
                    + Rational::from(1)
                    - Rational::from(2 * u32::from(power));
                if radial_degree > 0 {
                    return Ok(p.zero());
                }
                return Err(Error::Unsupported(format!(
                    "joint massless zero-density endpoint with radial degree {radial_degree} requires a distributional dimensional limit"
                )));
            }
            let threshold_exponent = spatial_dimension.clone() / Rational::from(2)
                + Rational::from(radial_power)
                - Rational::from(power - 1);
            if !self.mass_squared.is_zero() && threshold_exponent > 0 {
                // Both one-sided limits vanish; no arbitrary theta(0) value.
                return Ok(p.zero());
            }
            if threshold_exponent.is_zero() && power > 1 {
                // Positive mass, alpha=d/2+k_radial=power-1>=1. Regulate the
                // occupation by n_F((E-mu)/T) before taking mass derivatives.
                // The leading radial integral is proportional to
                // T^alpha int x^(alpha-1)n_F(x-gap/(2mT)) dx. Its alpha-th
                // gap derivative gives Gamma(alpha)*n_F(0), with n_F(0)=1/2;
                // derivatives of smooth prefactors vanish as T->0+. This is
                // the regulated finite jump, not an assignment of theta(0).
                let angular = Self::angular(spatial_dimension, p)?;
                let mass_power = Rational::from((i64::from(energy_power) - 1, 2));
                return Ok(p.scale(
                    &p.mul(
                        &angular,
                        &p.pow(&p.rational(&self.mass_squared), &p.rational(&mass_power)),
                    ),
                    -1,
                    8 * i64::from(power - 1),
                ));
            }
            return Err(Error::Numerical(format!(
                "raised occupied threshold has gap exponent {threshold_exponent}; no finite continuous pointwise value (thermal distributional limit required)"
            )));
        }
        let energy = symbol!("rustflow_density::seed_energy");
        let e = Atom::var(energy);
        let mut residue = EnergyResidue {
            bulk: -e.clone().pow(i64::from(energy_power) - 1) / Atom::num(2),
            upper_surface: vec![],
        };
        for n in 1..power {
            residue = residue.raise(energy, n)?;
        }
        let mut result = p.zero();
        for (monomial, coefficient) in
            exact_coefficient_list(&residue.bulk, std::slice::from_ref(&e))?
        {
            let degree = powers(&monomial, std::slice::from_ref(&e))?[0];
            let moment = self.moment(
                spatial_dimension,
                i32::from(degree) + 1,
                radial_power,
                digits,
                max_terms,
                p,
            )?;
            result = p.add(
                &result,
                &p.mul(
                    &p.scale(&p.eval(&coefficient, &HashMap::default())?, 2, 1),
                    &moment,
                ),
            );
        }
        let radial_exponent = (spatial_dimension.clone() - Rational::from(2)) / Rational::from(2)
            + Rational::from(radial_power);
        let jacobian = &e
            * (e.clone().pow(2) - Atom::num(self.mass_squared.clone()))
                .pow(Atom::num(radial_exponent));
        let params = HashMap::from_iter([(e.clone(), p.rational(&self.chemical_potential))]);
        let angular = Self::angular(spatial_dimension, p)?;
        for (order, coefficient) in residue.upper_surface.iter().enumerate() {
            let mut integrand = coefficient * &jacobian;
            for _ in 0..order {
                integrand = integrand.derivative(energy);
            }
            result = p.add(&result, &p.mul(&angular, &p.eval(&integrand, &params)?));
        }
        if !p.finite(&result) {
            return Err(Error::Numerical("nonfinite occupied seed".into()));
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finite_jump_threshold_matches_independent_thermal_radial_integrals() {
        // Differentiate the original fixed-spatial-momentum thermal kernels
        // before integration. With f(E)=n_F((E-mu)/T), d=2 gives
        // I2=1/(8pi) int_m^infty [f'/E-f/E²] dE,
        // I3[q²]=-1/(32pi) int_m^infty (E²-m²)
        //                       [f''/E²-3f'/E³+3f/E⁴] dE.
        // Independent integration by parts gives -f(m)/(8pi*m) and
        // -f(m)/(16pi*m), respectively. Integrate both kernels numerically,
        // including smooth bulk and endpoint contributions, at finite T.
        let mass = 0.5_f64;
        let pi = std::f64::consts::PI;
        for temperature in [0.125_f64, 0.03125, 0.0078125] {
            for offset in [-2.0_f64, 0.0, 2.0] {
                let mu = mass + offset * temperature;
                let intervals = 32768;
                let step = 80.0 / intervals as f64;
                let mut integrals = [0.0_f64; 2];
                for i in 0..=intervals {
                    let x = i as f64 * step;
                    let energy = mass + temperature * x;
                    let f = 1.0 / (1.0 + (x + (mass - mu) / temperature).exp());
                    let first = -f * (1.0 - f) / temperature;
                    let second = f * (1.0 - f) * (1.0 - 2.0 * f) / (temperature * temperature);
                    let kernels = [
                        temperature / (8.0 * pi) * (first / energy - f / energy.powi(2)),
                        -temperature / (32.0 * pi)
                            * (energy * energy - mass * mass)
                            * (second / energy.powi(2) - 3.0 * first / energy.powi(3)
                                + 3.0 * f / energy.powi(4)),
                    ];
                    let weight = if i == 0 || i == intervals {
                        1.0
                    } else if i % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    };
                    for (sum, kernel) in integrals.iter_mut().zip(kernels) {
                        *sum += weight * kernel;
                    }
                }
                let f_at_mass = 1.0 / (1.0 + (-offset).exp());
                for (integral, expected) in integrals.into_iter().zip([
                    -f_at_mass / (8.0 * pi * mass),
                    -f_at_mass / (16.0 * pi * mass),
                ]) {
                    assert!((step * integral / 3.0 - expected).abs() < 1e-10);
                }
            }
        }
        let shell = CompactShell {
            mass_squared: Rational::from((1, 4)),
            chemical_potential: Rational::from((1, 2)),
        };
        let p = Precision::decimal(50).unwrap();
        for (power, radial_power, expected) in [(2, 0, "-1/(8*pi)"), (3, 1, "-1/(16*pi)")] {
            let actual = shell
                .raised_moment(&Rational::from(2), power, 0, radial_power, 40, 10000, p)
                .unwrap();
            let expected = p.eval(&parse!(expected), &HashMap::default()).unwrap();
            assert!(p.close(&actual, &expected, 40));
        }
        assert!(
            matches!(shell.raised_moment(&Rational::from(2), 3, 0, 0, 40, 10000, p),
            Err(Error::Numerical(message)) if message.contains("gap exponent -1"))
        );
    }

    #[test]
    fn support_raised_surface_and_threshold_limits() {
        let p = Precision::decimal(45).unwrap();
        let d = Rational::from(3);
        let mut shell = CompactShell {
            mass_squared: Rational::from(1),
            chemical_potential: Rational::from(2),
        };
        // Direct radial integration at d=3, including both upper endpoint and
        // bulk derivative: I_2^occ=-acosh(mu/m)/(8*pi^2).
        let value = shell.raised_moment(&d, 2, 0, 0, 35, 10000, p).unwrap();
        let reference = p
            .eval(&parse!("-log(2+3^(1/2))/(8*pi^2)"), &HashMap::default())
            .unwrap();
        assert!(p.norm(&p.sub(&value, &reference)) < p.tolerance(32));
        // The simple moment independently checks the spatial measure factor.
        let value = shell.moment(&d, 0, 0, 35, 10000, p).unwrap();
        let reference = p
            .eval(
                &parse!("(2*3^(1/2)-log(2+3^(1/2)))/(8*pi^2)"),
                &HashMap::default(),
            )
            .unwrap();
        assert!(p.norm(&p.sub(&value, &reference)) < p.tolerance(32));
        shell.chemical_potential = Rational::from((1, 2));
        assert!(
            p.norm(&shell.raised_moment(&d, 2, 0, 0, 35, 10000, p).unwrap())
                .is_zero()
        );
        shell.chemical_potential = Rational::from(1);
        assert!(
            p.norm(&shell.raised_moment(&d, 2, 0, 0, 35, 10000, p).unwrap())
                .is_zero()
        );
        assert!(
            matches!(shell.raised_moment(&d, 3, 0, 0, 35, 10000, p), Err(Error::Numerical(message)) if message.contains("gap exponent -1/2"))
        );
    }
    #[test]
    fn massless_and_medium_seeds_have_defined_domains() {
        let p = Precision::decimal(40).unwrap();
        let shell = CompactShell {
            mass_squared: Rational::from(0),
            chemical_potential: Rational::from(2),
        };
        let value = shell.moment(&Rational::from(3), 2, 0, 30, 100, p).unwrap();
        let reference = p.eval(&parse!("1/pi^2"), &HashMap::default()).unwrap();
        assert!(p.norm(&p.sub(&value, &reference)) < p.tolerance(32));
        assert!(
            shell
                .raised_moment(&Rational::from(3), 2, 0, 0, 30, 100, p)
                .is_err()
        );
        let zero = CompactShell {
            chemical_potential: Rational::from(0),
            ..shell
        };
        assert!(
            p.norm(
                &zero
                    .raised_moment(&Rational::from(3), 1, 0, 0, 30, 100, p)
                    .unwrap()
            )
            .is_zero()
        );
        assert!(zero.moment(&Rational::from(3), -2, 0, 30, 100, p).is_err());
    }
}
