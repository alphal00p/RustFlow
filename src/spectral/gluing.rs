//! Assembly of glued integrals from one spectral sample.
//!
//! A [`GluingMap`] writes each glued integral as a sum of terms over the
//! targets of one auxiliary-mass family whose deformed line has mass^2
//! `z = 1 + eta`:
//!
//! * `spectral`: `coeff * (1/2 pi i) int_Gamma W(z) (z - 1)^p V_T(z) dz` with
//!   `W` a product of bubbles in series (`B_{1,1} = B_E`, `B_{2,1} = B21_E`)
//!   and `p = eta_power - chain_n`;
//! * `chain_const`: `coeff * beta_{m-1} * V_T(eta = 0)` with `beta` the Taylor
//!   coefficients of `B_E` at `z = 1` (a bubble chain contracted to a point);
//! * `product_T1`: `coeff * (-Gamma(eps - 1)) * V_T(eta = 0)` (a factorized
//!   one-loop tadpole).
//!
//! The contour integrals use the lower arm `z = c - i t` (real `eps`):
//! nodes on `t <= T0`, and the analytic tail of the large-`z` expansion built
//! from the target's `eta -> infinity` classes.
use std::collections::{BTreeMap, BTreeSet};

use serde::Deserialize;
use symbolica::domains::float::{Complex, Float};

use super::*;
use crate::engine::SpectralSample;
use crate::{ComplexFloat, Integral};

/// Gluing map (schema `sgw4-map-v1`); unknown fields are ignored.
#[derive(Clone, Debug, Deserialize)]
pub struct GluingMap {
    pub masters: BTreeMap<String, GluedIntegral>,
}

#[derive(Clone, Debug, Deserialize)]
pub struct GluedIntegral {
    pub terms: Vec<GluingTerm>,
    /// Loops of the glued integral for the `exp(loops * EulerGamma * eps)`
    /// normalization; the caller's default when absent.
    #[serde(default)]
    pub loops_norm: Option<i64>,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum GluingTerm {
    Spectral {
        coeff: String,
        /// Bubble line powers in pairs; each pair sums to 2 (`B_{1,1}`) or 3
        /// (`B_{2,1}`); several pairs are bubbles in series.
        weight: Vec<i64>,
        chain_n: i64,
        eta_power: i64,
        target: Vec<i16>,
    },
    ChainConst {
        coeff: String,
        weight: Vec<i64>,
        m: usize,
        target: Vec<i16>,
    },
    #[serde(rename = "product_T1")]
    ProductT1 { coeff: String, target: Vec<i16> },
}

impl GluingMap {
    /// Every target any term needs, in a canonical order.
    pub fn targets(&self) -> Vec<Integral> {
        let mut out = BTreeSet::new();
        for master in self.masters.values() {
            for term in &master.terms {
                out.insert(match term {
                    GluingTerm::Spectral { target, .. }
                    | GluingTerm::ChainConst { target, .. }
                    | GluingTerm::ProductT1 { target, .. } => target.clone(),
                });
            }
        }
        out.into_iter().map(Integral).collect()
    }
}

/// The lower contour arm and its quadrature.
#[derive(Clone, Debug)]
pub struct ArmSpec {
    /// Crossing point `c` of the contour with the real `z` axis (rational).
    pub c: String,
    /// End `T0` of the numerical part of the arm.
    pub t0: i64,
    /// Second split point for the tail-split self check (a panel edge).
    pub t_check: i64,
    /// Quadrature digits of the Gauss-Legendre panels.
    pub digits: u32,
}

impl Default for ArmSpec {
    fn default() -> Self {
        Self {
            c: "5/2".into(),
            t0: 120,
            t_check: 64,
            digits: 0,
        }
    }
}

/// Arm nodes plus the `eta` path and nodes for
/// [`crate::PreparedFlow::evaluate_spectral`].
pub struct SpectralArm {
    pub spec: ArmSpec,
    pub c: R,
    pub arm: Arm,
}

impl SpectralArm {
    /// Panels sized for the bubble singularities at `z = 0` and `z = 4`.
    pub fn new(ctx: &Ctx, spec: ArmSpec) -> Result<Self> {
        // Arm::new needs T0 > 1 for its geometric panel edges [0, 1, 2, ...],
        // and quadrature digits size the panels (0 would give 12-point panels).
        if spec.t0 <= 1 || spec.t_check <= 0 || spec.t_check > spec.t0 || spec.digits == 0 {
            return Err(Error::InvalidInput(format!(
                "spectral arm needs T0 > 1, 0 < t_check <= T0 and positive quadrature digits \
                 (T0 {}, t_check {}, digits {})",
                spec.t0, spec.t_check, spec.digits
            )));
        }
        let c = ctx.rational(&spec.c)?;
        let arm = Arm::new(ctx, &c, &ctx.r(spec.t0), &[0.0, 4.0], spec.digits);
        Ok(Self { spec, c, arm })
    }
    /// The nodes in `eta = z - 1`.
    pub fn eta_nodes(&self) -> Vec<ComplexFloat> {
        self.arm
            .z
            .iter()
            .map(|z| from_rug(&(z.clone() - 1i32)))
            .collect()
    }
    /// The transport waypoints in `eta`: `c - 1 - i T0`, then `c - 1`.
    pub fn eta_path(&self, ctx: &Ctx) -> Vec<ComplexFloat> {
        let bottom = ctx.cc(self.c.clone() - 1i32, -ctx.r(self.spec.t0));
        let top = ctx.cc(self.c.clone() - 1i32, ctx.r(0));
        vec![from_rug(&bottom), from_rug(&top)]
    }
}

pub fn to_rug(c: &ComplexFloat, bits: u32) -> C {
    C::with_val(bits, (c.re.as_raw(), c.im.as_raw()))
}
pub fn from_rug(c: &C) -> ComplexFloat {
    Complex::new(Float::from(c.real().clone()), Float::from(c.imag().clone()))
}

#[derive(Clone, Debug)]
pub struct TargetDiagnostic {
    pub target: Vec<i16>,
    pub classes: usize,
    /// Expansion versus transported value at the last arm node.
    pub asym_vs_node_at_t0_digits: f64,
    /// For targets without the deformed line: node value versus endpoint.
    pub eta_independent_node_vs_endpoint_digits: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct ContourIntegral {
    /// `"{target:?}|w{weight}|p{power}"`.
    pub key: String,
    pub x: R,
    pub tail_split_digits: f64,
    /// Coefficient of `z^-1` of the integrand expansion (its integer part).
    pub integer_z_inverse_coefficient: Option<C>,
}

#[derive(Clone, Debug)]
pub struct GluedValue {
    pub raw: R,
    pub normalized: R,
    pub tail_split_digits: f64,
}

#[derive(Clone, Debug)]
pub struct GluedSample {
    pub masters: BTreeMap<String, GluedValue>,
    pub contour_integrals: Vec<ContourIntegral>,
    pub diagnostics: Vec<TargetDiagnostic>,
    pub bubble_self_check_digits: f64,
    pub quality: GluingQuality,
}

/// The internal consistency of one glued sample, summarized for acceptance.
/// Digits are significant digits of agreement (999 for exact agreement).
#[derive(Clone, Debug)]
pub struct GluingQuality {
    /// Smallest tail-split agreement of a glued master.
    pub min_master_tail_split_digits: f64,
    /// Smallest tail-split agreement of a contour integral.
    pub min_contour_tail_split_digits: f64,
    /// Smallest agreement of a target's large-z expansion with its value
    /// transported to the last arm node.
    pub min_asym_vs_node_digits: f64,
    /// Smallest node-versus-endpoint agreement of the eta-independent
    /// targets, if there are any.
    pub min_eta_independent_digits: Option<f64>,
    pub bubble_self_check_digits: f64,
    /// [`SpectralSample::log_terms_residual`]: dropped logarithmic part of the
    /// large-z expansions relative to the kept part.
    pub log_terms_residual: f64,
    /// Largest |Im| of an endpoint value that `chain_const` / `product_T1`
    /// terms use through its real part, relative to the largest |Re| of all
    /// endpoint values (0 when no such term occurs).
    pub endpoint_imaginary_relative: f64,
}

/// Acceptance thresholds for [`GluingQuality`].
#[derive(Clone, Copy, Debug)]
pub struct GluingAcceptance {
    /// Every agreement must reach `min_digits`; the relative log residual
    /// and the relative imaginary part must stay below `10^-min_digits`.
    pub min_digits: f64,
}

/// Whether [`glue`] refuses a sample whose quality misses the thresholds.
#[derive(Clone, Copy, Debug)]
pub enum Acceptance {
    Require(GluingAcceptance),
    /// Explicit opt-out: return the sample with its measured quality only.
    Skip,
}

impl GluingAcceptance {
    /// The violated thresholds, empty when the sample is acceptable.
    pub fn violations(&self, quality: &GluingQuality) -> Vec<String> {
        let mut out = Vec::new();
        let mut digits = |name: &str, value: f64| {
            // NaN (an undefined agreement) is a violation.
            if !(value >= self.min_digits) {
                out.push(format!("{name} {value:.1} < {:.1} digits", self.min_digits));
            }
        };
        digits("master tail split", quality.min_master_tail_split_digits);
        digits("contour tail split", quality.min_contour_tail_split_digits);
        digits("expansion vs arm node", quality.min_asym_vs_node_digits);
        if let Some(value) = quality.min_eta_independent_digits {
            digits("eta-independent node vs endpoint", value);
        }
        digits("bubble self check", quality.bubble_self_check_digits);
        let bound = 10f64.powf(-self.min_digits);
        for (name, value) in [
            ("log terms residual", quality.log_terms_residual),
            (
                "endpoint imaginary part",
                quality.endpoint_imaginary_relative,
            ),
        ] {
            if !(value <= bound) {
                out.push(format!("{name} {value:.3e} > {bound:.1e}"));
            }
        }
        out
    }

    pub fn check(&self, quality: &GluingQuality) -> Result<()> {
        let violations = self.violations(quality);
        if violations.is_empty() {
            Ok(())
        } else {
            Err(Error::Accuracy(format!(
                "glued sample fails acceptance: {}",
                violations.join("; ")
            )))
        }
    }
}

/// Bubble-weight key: pair sums joined by `x` (`"2"` = `B_{1,1}`, `"3"` =
/// `B_{2,1}`, `"2x2"` = chain `B_{1,1} B_{1,1}`).
fn weight_key(w: &[i64]) -> Result<String> {
    if w.is_empty() || w.len() % 2 != 0 {
        return Err(Error::InvalidInput(
            "gluing weight must list bubble power pairs".into(),
        ));
    }
    // Only B_{1,1} and B_{2,1} = B_{1,2} are implemented; a pair such as
    // (0, 2) has the same sum but is not a bubble.
    if let Some(pair) = w
        .chunks(2)
        .find(|c| !matches!((c[0], c[1]), (1, 1) | (1, 2) | (2, 1)))
    {
        return Err(Error::Unsupported(format!(
            "gluing bubble line powers {pair:?}: only (1, 1), (1, 2) and (2, 1)"
        )));
    }
    Ok(w.chunks(2)
        .map(|c| (c[0] + c[1]).to_string())
        .collect::<Vec<_>>()
        .join("x"))
}

fn integral_key(target: &[i16], weight: &str, power: i64) -> String {
    format!("{target:?}|w{weight}|p{power}")
}

/// Assemble every glued integral of `map` from one spectral sample.
///
/// `targets` are the flow targets in the order of `sample`; `deformed_slot` is
/// the slot carrying `eta`; `order` truncates the large-`z` expansions at
/// `z^-order`; `loops_norm` is the default normalization loop count.
#[allow(clippy::too_many_arguments)]
pub fn glue(
    ctx: &Ctx,
    eps: &str,
    arm: &SpectralArm,
    order: usize,
    targets: &[Integral],
    deformed_slot: usize,
    sample: &SpectralSample,
    map: &GluingMap,
    loops_norm: i64,
    acceptance: Acceptance,
) -> Result<GluedSample> {
    if exact_rational(eps)? != exact_rational(&sample.epsilon.to_string())? {
        return Err(Error::InvalidInput(format!(
            "gluing at epsilon {eps} but the spectral sample is at {}",
            sample.epsilon
        )));
    }
    let e = ctx.rational(eps)?;
    let nmin = -(order as i64);
    let bubble = Bubble::new(*ctx, &e, &arm.c)?;
    let b_asym = bubble.asym(order as i64);
    let b21_asym = bubble.asym21(order as i64);
    let index: BTreeMap<&[i16], usize> = targets
        .iter()
        .enumerate()
        .map(|(i, t)| (t.0.as_slice(), i))
        .collect();
    let lookup = |target: &[i16]| {
        index.get(target).copied().ok_or_else(|| {
            Error::InvalidInput(format!("gluing target {target:?} is not a flow target"))
        })
    };
    let last_node = arm
        .arm
        .z
        .last()
        .ok_or_else(|| Error::InvalidInput("empty spectral arm".into()))?;
    if sample.node_values.len() != arm.arm.z.len()
        || sample.expansions.len() != targets.len()
        || sample.endpoint.len() != targets.len()
        || sample
            .node_values
            .iter()
            .any(|row| row.len() != targets.len())
    {
        return Err(Error::InvalidInput(
            "spectral sample does not match the arm or the targets".into(),
        ));
    }
    // Per-target large-z expansions and diagnostics.
    let mut v_asym = Vec::with_capacity(targets.len());
    let mut diagnostics = Vec::with_capacity(targets.len());
    for (ti, classes) in sample.expansions.iter().enumerate() {
        let mut integer_classes = vec![];
        for class in classes {
            let b = class.base.round();
            let k = class.slope.round();
            if (b - class.base).abs() > 1e-9 || (k - class.slope).abs() > 1e-9 {
                return Err(Error::Numerical(format!(
                    "non-integer exponent class base {} slope {}",
                    class.base, class.slope
                )));
            }
            integer_classes.push((
                b as i64,
                k as i64,
                class
                    .coefficients
                    .iter()
                    .map(|x| to_rug(x, ctx.bits))
                    .collect::<Vec<_>>(),
            ));
        }
        let va = eta_classes_to_z(ctx, &e, &integer_classes, nmin);
        let vlast = to_rug(&sample.node_values.last().unwrap()[ti], ctx.bits);
        let at_t0 = agree(&va.eval(ctx, &e, last_node), &vlast);
        let vend = to_rug(&sample.endpoint[ti], ctx.bits);
        // A negative power on the deformed slot is a numerator D + eta, so
        // only targets without the slot are eta independent.
        let eta_free = targets[ti]
            .0
            .get(deformed_slot)
            .is_none_or(|power| *power == 0);
        diagnostics.push(TargetDiagnostic {
            target: targets[ti].0.clone(),
            classes: classes.len(),
            asym_vs_node_at_t0_digits: at_t0,
            eta_independent_node_vs_endpoint_digits: eta_free
                .then(|| agree(&to_rug(&sample.node_values[0][ti], ctx.bits), &vend)),
        });
        v_asym.push(va);
    }
    // Contour integrals for every (target, weight, power) the map uses.
    let mut wanted = BTreeSet::new();
    for master in map.masters.values() {
        for term in &master.terms {
            if let GluingTerm::Spectral {
                weight,
                chain_n,
                eta_power,
                target,
                ..
            } = term
            {
                wanted.insert((target.clone(), weight_key(weight)?, eta_power - chain_n));
            }
        }
    }
    let t_check = ctx.r(arm.spec.t_check);
    let mut integrals: BTreeMap<String, (R, R)> = BTreeMap::new();
    let mut contour_integrals = vec![];
    for (target, weight, power) in &wanted {
        let ti = lookup(target)?;
        let sums = weight
            .split('x')
            .map(|x| {
                x.parse::<i64>()
                    .map_err(|e| Error::InvalidInput(e.to_string()))
            })
            .collect::<Result<Vec<_>>>()?;
        let mut wa = zm1_pow(ctx, 0, nmin);
        for sum in &sums {
            let one = match sum {
                2 => &b_asym,
                3 => &b21_asym,
                _ => {
                    return Err(Error::Unsupported(
                        "gluing bubble weight beyond B_{2,1}".into(),
                    ));
                }
            };
            wa = wa.mul(one, nmin);
        }
        let asym = wa
            .mul(&zm1_pow(ctx, *power, nmin), nmin)
            .mul(&v_asym[ti], nmin);
        let f = arm
            .arm
            .z
            .iter()
            .zip(&sample.node_values)
            .map(|(z, values)| {
                let mut wz = ctx.ci(1, 0);
                for sum in &sums {
                    wz *= if *sum == 2 {
                        bubble.b(z)?
                    } else {
                        bubble.b21(z)?
                    };
                }
                let zp = if *power == 0 {
                    ctx.ci(1, 0)
                } else {
                    (z.clone() - 1i32).pow(*power as i32)
                };
                Ok(wz * zp * to_rug(&values[ti], ctx.bits))
            })
            .collect::<Result<Vec<_>>>()?;
        let (_, x) = arm.arm.integrate(ctx, &e, &f, &asym);
        let (_, x2) = arm.arm.integrate_upto(ctx, &e, &f, &asym, &t_check);
        let key = integral_key(target, weight, *power);
        contour_integrals.push(ContourIntegral {
            key: key.clone(),
            x: x.clone(),
            tail_split_digits: agree_r(&x2, &x),
            integer_z_inverse_coefficient: asym.terms.get(&(0, -1)).cloned(),
        });
        integrals.insert(key, (x, x2));
    }
    // Assembly.
    let t1_raw = -ctx.gamma(&(e.clone() - 1i32));
    let taylor = bubble.taylor_at_one(8)?;
    let mut masters = BTreeMap::new();
    for (name, master) in &map.masters {
        let norm = (ctx.euler() * &e * master.loops_norm.unwrap_or(loops_norm)).exp();
        let mut raw = ctx.r(0);
        let mut raw2 = ctx.r(0);
        for term in &master.terms {
            match term {
                GluingTerm::Spectral {
                    coeff,
                    weight,
                    chain_n,
                    eta_power,
                    target,
                } => {
                    lookup(target)?;
                    let coeff = ctx.rational(coeff)?;
                    let key = integral_key(target, &weight_key(weight)?, eta_power - chain_n);
                    let (x, x2) = &integrals[&key];
                    raw += coeff.clone() * x;
                    raw2 += coeff * x2;
                }
                GluingTerm::ChainConst {
                    coeff,
                    weight,
                    m,
                    target,
                } => {
                    if weight_key(weight)? != "2" {
                        return Err(Error::Unsupported(
                            "gluing chain with a dotted bubble".into(),
                        ));
                    }
                    if *m == 0 || *m > taylor.len() {
                        return Err(Error::InvalidInput(format!("gluing chain order {m}")));
                    }
                    let coeff = ctx.rational(coeff)?;
                    let vend = to_rug(&sample.endpoint[lookup(target)?], ctx.bits)
                        .real()
                        .clone();
                    let cm = taylor[m - 1].real().clone();
                    raw += coeff.clone() * &cm * &vend;
                    raw2 += coeff * cm * &vend;
                }
                GluingTerm::ProductT1 { coeff, target } => {
                    let coeff = ctx.rational(coeff)?;
                    let vend = to_rug(&sample.endpoint[lookup(target)?], ctx.bits)
                        .real()
                        .clone();
                    raw += coeff.clone() * &t1_raw * &vend;
                    raw2 += coeff * &t1_raw * &vend;
                }
            }
        }
        let normalized = raw.clone() * &norm;
        masters.insert(
            name.clone(),
            GluedValue {
                tail_split_digits: agree_r(&raw2, &raw),
                raw,
                normalized,
            },
        );
    }
    let minimum = |values: &mut dyn Iterator<Item = f64>| values.fold(f64::INFINITY, f64::min);
    let largest_real = sample
        .endpoint
        .iter()
        .map(|v| to_rug(v, ctx.bits).real().clone().abs())
        .fold(ctx.r(0), |m, v| if v > m { v } else { m });
    let mut endpoint_imaginary = ctx.r(0);
    for master in map.masters.values() {
        for term in &master.terms {
            if let GluingTerm::ChainConst { target, .. } | GluingTerm::ProductT1 { target, .. } =
                term
            {
                let value = to_rug(&sample.endpoint[lookup(target)?], ctx.bits);
                let size = value.imag().clone().abs();
                if size > endpoint_imaginary {
                    endpoint_imaginary = size;
                }
            }
        }
    }
    let quality = GluingQuality {
        min_master_tail_split_digits: minimum(&mut masters.values().map(|m| m.tail_split_digits)),
        min_contour_tail_split_digits: minimum(
            &mut contour_integrals.iter().map(|x| x.tail_split_digits),
        ),
        min_asym_vs_node_digits: minimum(
            &mut diagnostics.iter().map(|d| d.asym_vs_node_at_t0_digits),
        ),
        min_eta_independent_digits: diagnostics
            .iter()
            .filter_map(|d| d.eta_independent_node_vs_endpoint_digits)
            .reduce(f64::min),
        bubble_self_check_digits: bubble.self_check_digits,
        log_terms_residual: sample.log_terms_residual,
        endpoint_imaginary_relative: if endpoint_imaginary.is_zero() {
            0.0
        } else if largest_real.is_zero() {
            f64::INFINITY
        } else {
            (endpoint_imaginary / largest_real).to_f64()
        },
    };
    if let Acceptance::Require(thresholds) = acceptance {
        thresholds.check(&quality)?;
    }
    Ok(GluedSample {
        masters,
        contour_integrals,
        diagnostics,
        bubble_self_check_digits: bubble.self_check_digits,
        quality,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bubble_weights_and_arm_specifications_are_validated() {
        assert_eq!(weight_key(&[1, 1]).unwrap(), "2");
        assert_eq!(weight_key(&[1, 2, 2, 1]).unwrap(), "3x3");
        for invalid in [&[0, 2][..], &[3, -1], &[0, 3], &[1], &[]] {
            assert!(weight_key(invalid).is_err(), "{invalid:?}");
        }
        let ctx = Ctx::digits(30);
        assert!(SpectralArm::new(&ctx, ArmSpec::default()).is_err());
        let spec = |t0, t_check| ArmSpec {
            t0,
            t_check,
            digits: 30,
            ..ArmSpec::default()
        };
        assert!(SpectralArm::new(&ctx, spec(120, 64)).is_ok());
        assert!(SpectralArm::new(&ctx, spec(1, 1)).is_err());
        assert!(SpectralArm::new(&ctx, spec(120, 121)).is_err());
    }

    #[test]
    fn acceptance_names_every_violated_threshold() {
        // The worst production values (sgw4 prod a0-000 and its reproduction).
        let production = GluingQuality {
            min_master_tail_split_digits: 103.0,
            min_contour_tail_split_digits: 102.7,
            min_asym_vs_node_digits: 96.1,
            min_eta_independent_digits: Some(155.4),
            bubble_self_check_digits: 139.4,
            log_terms_residual: 1.5e-194,
            endpoint_imaginary_relative: 1e-120,
        };
        let acceptance = GluingAcceptance { min_digits: 52.0 };
        assert!(acceptance.check(&production).is_ok());
        let broken = GluingQuality {
            min_asym_vs_node_digits: 12.0,
            min_eta_independent_digits: Some(f64::NAN),
            log_terms_residual: 1e-20,
            endpoint_imaginary_relative: f64::INFINITY,
            ..production
        };
        let violations = acceptance.violations(&broken);
        assert_eq!(violations.len(), 4, "{violations:?}");
        assert!(matches!(acceptance.check(&broken), Err(Error::Accuracy(_))));
        assert!(
            GluingAcceptance { min_digits: 104.0 }
                .violations(&production)
                .iter()
                .any(|v| v.starts_with("master tail split"))
        );
    }
}
