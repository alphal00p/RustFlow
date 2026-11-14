//! Exact univariate polynomial helpers over Q (ascending coefficients, rug numbers).
//!
//! Polynomials are kept in two forms: `QPoly` (rational coefficients) and
//! `ZPoly` (content x primitive integer polynomial). Shifts by -1 (z = 0 and
//! z = infinity expansions) and by dyadic complex centres (Taylor patches) are
//! done exactly in integer arithmetic before any rounding.
use rug::{Complex, Float, Integer, Rational};

pub type QPoly = Vec<Rational>;

pub fn trim(p: &mut QPoly) {
    while p.len() > 1 && *p.last().unwrap() == 0 {
        p.pop();
    }
}

pub fn is_zero(p: &[Rational]) -> bool {
    p.iter().all(|c| *c == 0)
}

pub fn add(a: &[Rational], b: &[Rational]) -> QPoly {
    let n = a.len().max(b.len());
    let mut out = vec![Rational::new(); n];
    for (i, x) in a.iter().enumerate() {
        out[i] += x;
    }
    for (i, x) in b.iter().enumerate() {
        out[i] += x;
    }
    trim(&mut out);
    out
}

pub fn mul(a: &[Rational], b: &[Rational]) -> QPoly {
    let mut out = vec![Rational::new(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if *x == 0 {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            if *y != 0 {
                out[i + j] += Rational::from(x * y);
            }
        }
    }
    trim(&mut out);
    out
}

pub fn scale(a: &[Rational], c: &Rational) -> QPoly {
    let mut out: QPoly = a.iter().map(|x| Rational::from(x * c)).collect();
    trim(&mut out);
    out
}

/// index of the first nonzero coefficient
pub fn valuation(p: &[Rational]) -> Option<usize> {
    p.iter().position(|c| *c != 0)
}

/// content x primitive integer polynomial (exact): p = content * sum prim_k x^k, gcd(prim) = 1
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ZPoly {
    pub content: Rational,
    pub prim: Vec<Integer>,
}

impl ZPoly {
    pub fn from_q(p: &[Rational]) -> ZPoly {
        let mut den = Integer::from(1);
        for c in p {
            den.lcm_mut(c.denom());
        }
        let mut ints: Vec<Integer> = p
            .iter()
            .map(|c| c.numer() * Integer::from(&den / c.denom()))
            .collect();
        let mut g = Integer::new();
        for x in &ints {
            g.gcd_mut(x);
        }
        if g == 0 {
            return ZPoly {
                content: Rational::new(),
                prim: vec![Integer::new()],
            };
        }
        for x in ints.iter_mut() {
            x.div_exact_mut(&g);
        }
        while ints.len() > 1 && *ints.last().unwrap() == 0 {
            ints.pop();
        }
        ZPoly {
            content: Rational::from((g, den)),
            prim: ints,
        }
    }

    pub fn to_q(&self) -> QPoly {
        self.prim
            .iter()
            .map(|x| Rational::from(x * &self.content.clone()))
            .collect()
    }

    pub fn degree(&self) -> usize {
        self.prim.len() - 1
    }
}

/// p(x + m) for an integer m, exact (Horner Taylor shift: additions and small multiplications only)
pub fn taylor_shift_int(p: &[Integer], m: &Integer) -> Vec<Integer> {
    let mut b = p.to_vec();
    let n = b.len();
    if n <= 1 || *m == 0 {
        return b;
    }
    let mut t = Integer::new();
    for i in 0..n - 1 {
        for j in (i..n - 1).rev() {
            // b[j] += m * b[j + 1]
            t.assign_mul(m, &b[j + 1]);
            b[j] += &t;
        }
    }
    b
}

trait AssignMul {
    fn assign_mul(&mut self, a: &Integer, b: &Integer);
}
impl AssignMul for Integer {
    #[inline]
    fn assign_mul(&mut self, a: &Integer, b: &Integer) {
        use rug::Assign;
        self.assign(a * b);
    }
}

/// p(x + m) for a Gaussian integer m = mr + i mi, exact; returns (re, im) coefficient vectors
pub fn taylor_shift_gauss(
    p: &[Integer],
    mr: &Integer,
    mi: &Integer,
) -> (Vec<Integer>, Vec<Integer>) {
    if *mi == 0 {
        let re = taylor_shift_int(p, mr);
        let im = vec![Integer::new(); re.len()];
        return (re, im);
    }
    let n = p.len();
    let mut re = p.to_vec();
    let mut im = vec![Integer::new(); n];
    let mut t = Integer::new();
    for i in 0..n.saturating_sub(1) {
        for j in (i..n - 1).rev() {
            // (re + i im)[j] += (mr + i mi) (re + i im)[j + 1]
            let (r1, i1) = (re[j + 1].clone(), im[j + 1].clone());
            t.assign_mul(mr, &r1);
            re[j] += &t;
            t.assign_mul(mi, &i1);
            re[j] -= &t;
            t.assign_mul(mr, &i1);
            im[j] += &t;
            t.assign_mul(mi, &r1);
            im[j] += &t;
        }
    }
    (re, im)
}

/// dyadic decomposition of a rational with a power-of-two denominator: x = m / 2^e (e >= 0)
pub fn dyadic(x: &Rational) -> Option<(Integer, u32)> {
    let d = x.denom();
    if !d.is_power_of_two() {
        return None;
    }
    Some((x.numer().clone(), d.significant_bits() - 1))
}

/// Taylor coefficients of p at the complex rational centre c = (cr, ci) (both dyadic or both arbitrary):
/// returns the coefficients of p(c + t) as (re, im) rationals. Dyadic centres use the exact integer shift.
pub fn shift_complex(p: &ZPoly, cr: &Rational, ci: &Rational) -> (QPoly, QPoly) {
    let deg = p.degree();
    if let (Some((mr, er)), Some((mi, ei))) = (dyadic(cr), dyadic(ci)) {
        let e = er.max(ei);
        let mr = mr << (e - er);
        let mi = mi << (e - ei);
        // b_k = a_k 2^(e (deg - k));  sum b_k (m + w)^k = 2^(e deg) p(c + w / 2^e)
        let b: Vec<Integer> = p
            .prim
            .iter()
            .enumerate()
            .map(|(k, a)| Integer::from(a << (e as usize * (deg - k)) as u32))
            .collect();
        let (re, im) = taylor_shift_gauss(&b, &mr, &mi);
        // coefficient of t^n: 2^(-e deg) * 2^(e n) * shifted[n] = shifted[n] / 2^(e (deg - n))
        let conv = |v: Vec<Integer>| -> QPoly {
            v.into_iter()
                .enumerate()
                .map(|(n, x)| {
                    let mut r =
                        Rational::from((x, Integer::from(1) << (e as usize * (deg - n)) as u32));
                    r *= &p.content;
                    r
                })
                .collect()
        };
        return (conv(re), conv(im));
    }
    // generic exact rational Horner (rarely used)
    let q = p.to_q();
    let mut re: QPoly = Vec::new();
    let mut im: QPoly = Vec::new();
    for a in q.iter().rev() {
        let k = re.len();
        let mut nr = vec![Rational::new(); k + 1];
        let mut ni = vec![Rational::new(); k + 1];
        for i in 0..k {
            nr[i + 1] += &re[i];
            ni[i + 1] += &im[i];
            nr[i] += Rational::from(&re[i] * cr) - Rational::from(&im[i] * ci);
            ni[i] += Rational::from(&re[i] * ci) + Rational::from(&im[i] * cr);
        }
        nr[0] += a;
        re = nr;
        im = ni;
    }
    (re, im)
}

/// p(x - 1) exact, as rationals
pub fn shift_minus_one(p: &ZPoly) -> QPoly {
    let s = taylor_shift_int(&p.prim, &Integer::from(-1));
    s.into_iter()
        .map(|x| Rational::from(x) * &p.content)
        .collect()
}

/// exact remainder test / division over Q of primitive integer polynomials: Some(a / b) primitive when b | a
pub fn exact_div(a: &[Integer], b: &[Integer]) -> Option<Vec<Integer>> {
    let db = b.len() - 1;
    if a.len() < b.len() {
        return None;
    }
    let mut r: Vec<Rational> = a.iter().map(Rational::from).collect();
    let mut q = vec![Rational::new(); a.len() - db];
    let lead = Rational::from(b[db].clone());
    for i in (db..a.len()).rev() {
        let c = Rational::from(&r[i] / &lead);
        if c != 0 {
            for j in 0..=db {
                r[i - db + j] -= Rational::from(&c * &b[j]);
            }
        }
        q[i - db] = c;
    }
    if r[..db].iter().any(|x| *x != 0) {
        return None;
    }
    Some(ZPoly::from_q(&q).prim_signed())
}

impl ZPoly {
    /// primitive part with a positive leading coefficient
    pub fn prim_signed(self) -> Vec<Integer> {
        let mut p = self.prim;
        if *p.last().unwrap() < 0 {
            for x in p.iter_mut() {
                *x = Integer::from(-&*x);
            }
        }
        p
    }
}

/// remainder of a by b mod p (p prime, b's leading coefficient invertible mod p) is zero
pub fn divides_mod_p(a: &[Integer], b: &[Integer], p: u64) -> bool {
    let db = b.len() - 1;
    if a.len() < b.len() {
        return false;
    }
    let md = |x: &Integer| -> u64 { x.mod_u(p as u32) as u64 };
    let bl = md(&b[db]);
    if bl == 0 {
        return false;
    }
    let inv = pow_mod(bl, p - 2, p);
    let bm: Vec<u64> = b.iter().map(md).collect();
    let mut am: Vec<u64> = a.iter().map(md).collect();
    for i in (db..am.len()).rev() {
        let c = am[i] * inv % p;
        if c != 0 {
            for j in 0..=db {
                let k = i - db + j;
                am[k] = (am[k] + p - c * bm[j] % p) % p;
            }
        }
    }
    am[..db].iter().all(|&x| x == 0)
}

pub fn pow_mod(mut b: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1u64;
    b %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % p;
        }
        b = b * b % p;
        e >>= 1;
    }
    r
}

/// rational -> residue mod p (denominator must be invertible)
pub fn rat_mod_p(x: &Rational, p: u64) -> u64 {
    let n = x.numer().mod_u(p as u32) as u64;
    let d = x.denom().mod_u(p as u32) as u64;
    assert!(d != 0, "denominator divisible by the prime");
    n * pow_mod(d, p - 2, p) % p
}

pub fn fl(prec: u32, x: &Rational) -> Float {
    Float::with_val(prec, x)
}

pub fn cx(prec: u32, re: &Rational, im: &Rational) -> Complex {
    Complex::with_val(prec, (re, im))
}

/// parse "a/b" or "a" (also accepts a JSON number)
pub fn parse_rational(v: &serde_json::Value) -> Result<Rational, String> {
    match v {
        serde_json::Value::String(s) => s
            .trim()
            .parse::<Rational>()
            .map_err(|e| format!("rational {s:?}: {e}")),
        serde_json::Value::Number(n) => n.as_i64().map(Rational::from).ok_or_else(|| {
            format!("non-integer JSON number {n} where a rational string was expected")
        }),
        _ => Err(format!("not a rational: {v}")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn q(v: &[i64]) -> QPoly {
        v.iter().map(|&x| Rational::from(x)).collect()
    }
    #[test]
    fn shifts_are_exact() {
        // p = 3 - 2x + x^3 / 4
        let p: QPoly = vec![
            Rational::from(3),
            Rational::from(-2),
            Rational::new(),
            Rational::from((1, 4)),
        ];
        let z = ZPoly::from_q(&p);
        assert_eq!(z.to_q(), p);
        // p(x - 1)
        let s = shift_minus_one(&z);
        let x = Rational::from((5, 7));
        let ev = |c: &[Rational], x: &Rational| {
            c.iter().rev().fold(Rational::new(), |acc, a| acc * x + a)
        };
        assert_eq!(ev(&s, &x), ev(&p, &(x.clone() - 1)));
        // complex dyadic centre c = 3/8 + i 5/16: p(c + t) at t = 1/3 against direct complex evaluation
        let (cr, ci) = (Rational::from((3, 8)), Rational::from((5, 16)));
        let (re, im) = shift_complex(&z, &cr, &ci);
        let t = Rational::from((1, 3));
        let (vr, vi) = (ev(&re, &t), ev(&im, &t));
        // direct: (a + ib)^k
        let (ar, ai) = (cr.clone() + &t, ci.clone());
        let (mut pr, mut pi) = (Rational::from(1), Rational::new());
        let (mut sr, mut si) = (Rational::new(), Rational::new());
        for c in &p {
            sr += Rational::from(c * &pr);
            si += Rational::from(c * &pi);
            let nr = Rational::from(&pr * &ar) - Rational::from(&pi * &ai);
            let ni = Rational::from(&pr * &ai) + Rational::from(&pi * &ar);
            pr = nr;
            pi = ni;
        }
        assert_eq!((vr, vi), (sr, si));
        let _ = q(&[1]);
    }
    #[test]
    fn factor_division() {
        let a = vec![Integer::from(-2), Integer::from(-1), Integer::from(1)]; // (x - 2)(x + 1)
        let b = vec![Integer::from(1), Integer::from(1)];
        assert!(divides_mod_p(&a, &b, 1048573));
        assert_eq!(
            exact_div(&a, &b).unwrap(),
            vec![Integer::from(-2), Integer::from(1)]
        );
        assert!(exact_div(&a, &[Integer::from(3), Integer::from(1)]).is_none());
    }
}
