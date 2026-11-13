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
}

/// Bubble-weight key: pair sums joined by `x` (`"2"` = `B_{1,1}`, `"3"` =
/// `B_{2,1}`, `"2x2"` = chain `B_{1,1} B_{1,1}`).
fn weight_key(w: &[i64]) -> Result<String> {
    if w.is_empty() || w.len() % 2 != 0 {
        return Err(Error::InvalidInput(
            "gluing weight must list bubble power pairs".into(),
        ));
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
) -> Result<GluedSample> {
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
    if sample.node_values.len() != arm.arm.z.len() || sample.expansions.len() != targets.len() {
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
                    if weight.iter().sum::<i64>() != 2 || weight.len() != 2 {
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
    Ok(GluedSample {
        masters,
        contour_integrals,
        diagnostics,
        bubble_self_check_digits: bubble.self_check_digits,
    })
}
