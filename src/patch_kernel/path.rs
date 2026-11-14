//! Planned paths of Taylor patches.
//!
//! A path lives in the path variable `z = x + offset` (`x` = system variable;
//! the dispersive mode uses `z = 1 + eta`, offset 1). Patch `k` expands at
//! `centre_k` and steps by `step_k`; consecutive patches must be contiguous
//! exactly (`centre_{k+1} = centre_k + step_k`).
//!
//! [`plan_axis_detours`] is a port of the campaign planner `zplan.py`
//! (values/d5l/v9c/tools/zplan.py): from `-s0` to `-S` along the negative
//! real axis, a 12-vertex polygon above every singular point on the axis,
//! steps `R / ratio` rounded to 24 significant bits (exact dyadics), order
//! `ceil((W + 8) / log10(R / |h|)) + 4`. The same points with the orders of
//! another `W` give the second run of a two-precision comparison.
use crate::{Error, Result};
use rug::{Float, Rational};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PatchSpec {
    pub centre: (Rational, Rational),
    pub step: (Rational, Rational),
    pub order: usize,
}

impl PatchSpec {
    pub fn is_real(&self) -> bool {
        self.centre.1 == 0 && self.step.1 == 0
    }
    pub fn end(&self) -> (Rational, Rational) {
        (
            Rational::from(&self.centre.0 + &self.step.0),
            Rational::from(&self.centre.1 + &self.step.1),
        )
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PlannedPath {
    /// `z = x + offset`.
    pub offset: Rational,
    pub patches: Vec<PatchSpec>,
}

#[derive(Serialize, Deserialize)]
struct PathFile {
    #[serde(default)]
    offset: Option<String>,
    patches: Vec<(String, String, String, String, usize)>,
    #[serde(flatten)]
    extra: serde_json::Map<String, serde_json::Value>,
}

fn q(text: &str) -> Result<Rational> {
    text.parse::<Rational>()
        .map_err(|e| Error::InvalidInput(format!("rational {text:?}: {e}")))
}

impl PlannedPath {
    /// Read a zplan-style file `{"patches": [[zr, zi, hr, hi, order], ..]}`;
    /// `offset` from the file if present, else the given default.
    pub fn from_json(value: serde_json::Value, default_offset: &Rational) -> Result<Self> {
        let f: PathFile = serde_json::from_value(value)
            .map_err(|e| Error::InvalidInput(format!("path file: {e}")))?;
        let offset = match &f.offset {
            Some(o) => q(o)?,
            None => default_offset.clone(),
        };
        let patches = f
            .patches
            .iter()
            .map(|(a, b, c, d, n)| {
                Ok(PatchSpec {
                    centre: (q(a)?, q(b)?),
                    step: (q(c)?, q(d)?),
                    order: *n,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let path = Self { offset, patches };
        path.check_contiguous()?;
        Ok(path)
    }

    pub fn to_json(&self) -> serde_json::Value {
        let patches: Vec<serde_json::Value> = self
            .patches
            .iter()
            .map(|p| {
                serde_json::json!([
                    p.centre.0.to_string(),
                    p.centre.1.to_string(),
                    p.step.0.to_string(),
                    p.step.1.to_string(),
                    p.order
                ])
            })
            .collect();
        serde_json::json!({
            "offset": self.offset.to_string(),
            "patches": patches,
            "n": self.patches.len(),
            "n_complex": self.patches.iter().filter(|p| !p.is_real()).count(),
            "sum_nord": self.patches.iter().map(|p| p.order).sum::<usize>(),
        })
    }

    pub fn check_contiguous(&self) -> Result<()> {
        for (k, w) in self.patches.windows(2).enumerate() {
            if w[0].end() != w[1].centre {
                return Err(Error::InvalidInput(format!(
                    "path not contiguous between patches {k} and {}",
                    k + 1
                )));
            }
        }
        Ok(())
    }

    /// Sub-path of patches `range` (keeps the offset).
    pub fn segment(&self, range: std::ops::Range<usize>) -> Result<Self> {
        if range.end > self.patches.len() || range.start >= range.end {
            return Err(Error::InvalidInput("bad path segment".into()));
        }
        Ok(Self {
            offset: self.offset.clone(),
            patches: self.patches[range].to_vec(),
        })
    }

    /// blake3 of offset, points, steps and orders.
    pub fn fingerprint(&self) -> String {
        let mut h = blake3::Hasher::new();
        h.update(b"rustflow.patch-kernel.path.v1;");
        h.update(self.offset.to_string().as_bytes());
        for p in &self.patches {
            for v in [&p.centre.0, &p.centre.1, &p.step.0, &p.step.1] {
                h.update(b";");
                h.update(v.to_string().as_bytes());
            }
            h.update(format!(";{}", p.order).as_bytes());
        }
        h.finalize().to_hex().to_string()
    }

    /// The same points with orders for `digits` working digits:
    /// `ceil((W + 8) / log10(R / |h|)) + 4` with `R` the distance from the
    /// centre to the nearest singular point (path variable).
    pub fn with_orders(&self, singular_points: &[(f64, f64)], digits: f64) -> Self {
        let mut out = self.clone();
        for p in &mut out.patches {
            let c = (to_f64(&p.centre.0), to_f64(&p.centre.1));
            let r = radius(singular_points, c);
            let h = to_f64(&p.step.0).hypot(to_f64(&p.step.1));
            p.order = order_rule(digits, r, h);
        }
        out
    }
}

/// Correctly rounded conversion (as Python's `float(Fraction)`).
pub fn to_f64(x: &Rational) -> f64 {
    Float::with_val(53, x).to_f64()
}

fn radius(pts: &[(f64, f64)], z: (f64, f64)) -> f64 {
    pts.iter()
        .map(|s| (z.0 - s.0).hypot(z.1 - s.1))
        .fold(f64::INFINITY, f64::min)
}

fn order_rule(digits: f64, r: f64, h: f64) -> usize {
    ((digits + 8.0) / (r / h).log10()).ceil() as usize + 4
}

/// `x` rounded to 24 significant bits, as an exact dyadic (zplan `rr`).
fn rr(x: f64) -> Rational {
    if x == 0.0 {
        return Rational::new();
    }
    // math.frexp: x = m 2^e with 0.5 <= |m| < 1
    let e = frexp_exp(x);
    let k = 24 - e;
    let scaled = (x * 2f64.powi(k)).round_ties_even();
    let mut r = Rational::from_f64(scaled).expect("finite");
    if k >= 0 {
        r >>= k as u32;
    } else {
        r <<= (-k) as u32;
    }
    r
}

fn frexp_exp(x: f64) -> i32 {
    let bits = x.abs().to_bits();
    let exp = ((bits >> 52) & 0x7ff) as i32;
    if exp == 0 {
        // subnormal: x = m 2^-1074 = (m / 2^L) 2^(L - 1074), L = bit length of m
        let m = bits & ((1u64 << 52) - 1);
        (64 - m.leading_zeros() as i32) - 1074
    } else {
        exp - 1022
    }
}

/// Inputs of the axis planner (the `plan`/`singular_points`/`s0`/`S` data
/// of the dispersive preparation).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PlanInput {
    /// Start `z = -s0` (exact rational or float string).
    pub s0: String,
    /// End `z = -S`.
    pub s_end: String,
    /// All finite singular points of the system in `z`.
    pub singular_points: Vec<(f64, f64)>,
    /// Singular points on the negative real axis, in path order.
    pub axis_singular: Vec<f64>,
    #[serde(default = "default_ratio")]
    pub ratio: f64,
}

fn default_ratio() -> f64 {
    3.0
}

fn parse_point(text: &str) -> Result<Rational> {
    if let Ok(r) = text.parse::<Rational>() {
        return Ok(r);
    }
    let f: f64 = text
        .parse()
        .map_err(|e| Error::InvalidInput(format!("point {text:?}: {e}")))?;
    Rational::from_f64(f).ok_or_else(|| Error::InvalidInput(format!("point {text:?}")))
}

/// Port of zplan.py (see module docs); `digits` = working digits W.
pub fn plan_axis_detours(input: &PlanInput, digits: f64, offset: Rational) -> Result<PlannedPath> {
    let pts = &input.singular_points;
    let ratio = input.ratio;
    let s0 = parse_point(&input.s0)?;
    let s_end = parse_point(&input.s_end)?;
    let mut patches = Vec::new();
    let zero = Rational::new();
    let push = |patches: &mut Vec<PatchSpec>,
                c: (Rational, Rational),
                h: (Rational, Rational),
                n: usize| {
        patches.push(PatchSpec {
            centre: c,
            step: h,
            order: n,
        });
    };
    let step_real =
        |patches: &mut Vec<PatchSpec>, mut zq: Rational, target: &Rational| -> Rational {
            while zq > *target {
                let zc = to_f64(&zq);
                let r = radius(pts, (zc, 0.0));
                let h = r / ratio;
                let hq = if zc - h * 1.5 <= to_f64(target) {
                    Rational::from(target - &zq)
                } else {
                    -rr(h)
                };
                let n = order_rule(digits, r, to_f64(&hq).abs());
                push(
                    patches,
                    (zq.clone(), Rational::new()),
                    (hq.clone(), Rational::new()),
                    n,
                );
                zq += &hq;
            }
            zq
        };
    let mut zq = -s0;
    let neg_s = to_f64(&Rational::from(-&s_end));
    for &xp in &input.axis_singular {
        if !(neg_s < xp && xp < to_f64(&zq)) {
            continue;
        }
        let near = pts
            .iter()
            .filter(|s| (s.0 - xp).hypot(s.1) > 1e-9 * (1.0 + xp.abs()))
            .map(|s| (xp - s.0).hypot(-s.1))
            .chain(std::iter::once(xp.abs()))
            .fold(f64::INFINITY, f64::min);
        let rd = 0.4 * near;
        let start = rr(xp + rd);
        let end = rr(xp - rd);
        if to_f64(&start) < to_f64(&zq) || start < zq {
            zq = step_real(&mut patches, zq, &start);
        }
        let nv = 12;
        let mut verts: Vec<(Rational, Rational)> = (1..nv)
            .map(|k| {
                let a = std::f64::consts::PI * k as f64 / nv as f64;
                (rr(xp + rd * a.cos()), rr(rd * a.sin()))
            })
            .collect();
        verts.push((end, Rational::new()));
        let (mut cr, mut ci) = (zq.clone(), Rational::new());
        for (vr, vi) in verts {
            while (&cr, &ci) != (&vr, &vi) {
                let cur = (to_f64(&cr), to_f64(&ci));
                let r = radius(pts, cur);
                let dv = (
                    to_f64(&Rational::from(&vr - &cr)),
                    to_f64(&Rational::from(&vi - &ci)),
                );
                let adv = dv.0.hypot(dv.1);
                let hq = if adv <= 1.5 * r / ratio {
                    (Rational::from(&vr - &cr), Rational::from(&vi - &ci))
                } else {
                    let hh = ((dv.0 / adv) * r / ratio, (dv.1 / adv) * r / ratio);
                    (rr(hh.0), rr(hh.1))
                };
                let hn = to_f64(&hq.0).hypot(to_f64(&hq.1));
                let n = order_rule(digits, r, hn);
                push(&mut patches, (cr.clone(), ci.clone()), hq.clone(), n);
                cr += &hq.0;
                ci += &hq.1;
            }
        }
        if ci != zero {
            return Err(Error::Numerical(
                "detour does not return to the axis".into(),
            ));
        }
        zq = cr;
    }
    let target = Rational::from(-&s_end);
    zq = step_real(&mut patches, zq, &target);
    if zq != target {
        return Err(Error::Numerical("planner overshoot".into()));
    }
    let path = PlannedPath { offset, patches };
    path.check_contiguous()?;
    Ok(path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rr_rounds_to_24_bits() {
        assert_eq!(rr(1.0), Rational::from(1));
        assert_eq!(rr(-0.75), Rational::from((-3, 4)));
        let x = 1.0 / 3.0;
        let r = rr(x);
        assert_eq!(r.denom().significant_bits(), 26); // 2^25
        assert!((to_f64(&r) - x).abs() < 2f64.powi(-25));
    }

    #[test]
    fn frexp_matches_python() {
        assert_eq!(frexp_exp(1.0), 1);
        assert_eq!(frexp_exp(0.5), 0);
        assert_eq!(frexp_exp(0.75), 0);
        assert_eq!(frexp_exp(3.0), 2);
        assert_eq!(frexp_exp(1.4155343563970746e-15), -49);
    }
}
