//! Exact Taylor shift of integer polynomials to a Gaussian-rational centre.
//!
//! For `c = (gr + i gi) / q` and `P(x) = sum_k A_k x^k`,
//! `P(c + x) = q^-n sum_m R_m q^m x^m` where `R(u) = B(u + g)`,
//! `B_k = A_k q^(n - k)`; `R` is obtained by integer synthetic division (only
//! multiply-adds by the small Gaussian integer `g`, no gcds). The coefficient
//! `R_m / q^(n - m)` is then rounded once to a guarded precision.
use rug::{Assign, Complex, Float, Integer, Rational};

/// A Gaussian-rational expansion centre.
#[derive(Clone, Debug)]
pub struct Centre {
    pub gr: Integer,
    pub gi: Integer,
    /// Common positive denominator.
    pub q: Integer,
    /// `q = 2^s` (exact binary scaling).
    q_log2: Option<u32>,
}

impl Centre {
    pub fn new(re: &Rational, im: &Rational) -> Self {
        let q = Integer::from(re.denom().lcm_ref(im.denom()));
        let gr = re.numer() * Integer::from(&q / re.denom());
        let gi = im.numer() * Integer::from(&q / im.denom());
        let q_log2 = q.is_power_of_two().then(|| q.significant_bits() - 1);
        Self { gr, gi, q, q_log2 }
    }
    pub fn is_real(&self) -> bool {
        self.gi == 0
    }
}

/// `p(c + x)` coefficients, each rounded once to `gbits` (complex; imaginary
/// parts are exactly zero for a real centre).
pub fn shift_guarded(p: &[Integer], c: &Centre, gbits: u32) -> Vec<Complex> {
    let n = p.len() - 1;
    // B_k = A_k q^(n-k)
    let mut br: Vec<Integer> = Vec::with_capacity(n + 1);
    if let Some(s) = c.q_log2 {
        for (k, a) in p.iter().enumerate() {
            br.push(Integer::from(a << (s * (n - k) as u32)));
        }
    } else {
        let mut pw = Integer::from(1);
        let mut rev = Vec::with_capacity(n + 1);
        for a in p.iter().rev() {
            rev.push(Integer::from(a * &pw));
            pw *= &c.q;
        }
        rev.reverse();
        br = rev;
    }
    let real = c.is_real();
    let mut bi: Vec<Integer> = if real {
        Vec::new()
    } else {
        vec![Integer::new(); n + 1]
    };
    if c.gr != 0 || c.gi != 0 {
        let mut t = Integer::new();
        for i in 0..n {
            for k in (i..n).rev() {
                let (lo, hi) = br.split_at_mut(k + 1);
                if real {
                    lo[k] += &hi[0] * &c.gr;
                } else {
                    let (li, hii) = bi.split_at_mut(k + 1);
                    // (r + i m)(gr + i gi) = (r gr - m gi) + i (r gi + m gr)
                    t.assign(&hi[0] * &c.gi);
                    li[k] += &t;
                    li[k] += &hii[0] * &c.gr;
                    lo[k] += &hi[0] * &c.gr;
                    lo[k] -= &hii[0] * &c.gi;
                }
            }
        }
    }
    let mut out = Vec::with_capacity(n + 1);
    for m in 0..=n {
        let e = n - m;
        // Correctly rounded R_m / q^e: exact binary scaling, or the quotient
        // of two exactly represented integers.
        let conv = |v: &Integer| -> Float {
            if e == 0 || c.q_log2.is_some() {
                let mut f = Float::with_val(gbits, v);
                if let Some(s) = c.q_log2 {
                    f >>= s * e as u32;
                }
                f
            } else {
                let d = Integer::from(rug::ops::Pow::pow(&c.q, e as u32));
                let exact = |z: &Integer| Float::with_val(z.significant_bits().max(1), z);
                Float::with_val(gbits, &exact(v) / &exact(&d))
            }
        };
        let re = conv(&br[m]);
        let im = if real {
            Float::new(gbits)
        } else {
            conv(&bi[m])
        };
        out.push(Complex::with_val(gbits, (re, im)));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn eval_exact(p: &[Integer], re: &Rational, im: &Rational) -> (Rational, Rational) {
        let (mut ar, mut ai) = (Rational::new(), Rational::new());
        for a in p.iter().rev() {
            let nr = Rational::from(&ar * re) - Rational::from(&ai * im) + a;
            let ni = Rational::from(&ar * im) + Rational::from(&ai * re);
            ar = nr;
            ai = ni;
        }
        (ar, ai)
    }

    #[test]
    fn shift_matches_exact_taylor_coefficients() {
        let p: Vec<Integer> = [3, -7, 0, 11, 5, -2]
            .iter()
            .map(|&v| Integer::from(v))
            .collect();
        for (re, im) in [
            ("-51/36028797018963968", "0"),
            ("-6468430547/8589934592", "3/7"),
            ("5/3", "-1/9"),
            ("0", "0"),
        ] {
            let re: Rational = re.parse().unwrap();
            let im: Rational = im.parse().unwrap();
            let c = Centre::new(&re, &im);
            let s = shift_guarded(&p, &c, 300);
            // value p(c) and derivative p'(c)
            let (v_re, v_im) = eval_exact(&p, &re, &im);
            assert_eq!(Float::with_val(300, &v_re), *s[0].real());
            assert_eq!(Float::with_val(300, &v_im), *s[0].imag());
            let dp: Vec<Integer> = p
                .iter()
                .enumerate()
                .skip(1)
                .map(|(k, a)| Integer::from(a * k as u32))
                .collect();
            let (d_re, d_im) = eval_exact(&dp, &re, &im);
            assert_eq!(Float::with_val(300, &d_re), *s[1].real());
            assert_eq!(Float::with_val(300, &d_im), *s[1].imag());
            // leading coefficient unchanged
            assert_eq!(Float::with_val(300, -2), *s[5].real());
        }
    }
}
