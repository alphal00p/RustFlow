//! Exact rational normal form for dense identity hashing, not a system rewrite.
use super::*;
use symbolica::domains::rational_polynomial::FromNumeratorAndDenominator;

// The pinned converter casts powers through i32 even when E=u32. Bound
// numerator/denominator degree conservatively before entering that converter;
// sums/products add bounds and integer powers multiply them. Larger supported
// Atom expressions keep their structural identity instead of risking a cast.
fn conversion_degree(value: AtomView<'_>) -> Option<u64> {
    const MAX: u64 = i32::MAX as u64;
    let degree = match value {
        AtomView::Num(_) => 0,
        AtomView::Var(_) => 1,
        AtomView::Add(a) => a
            .iter()
            .try_fold(0u64, |sum, a| sum.checked_add(conversion_degree(a)?))?,
        AtomView::Mul(a) => a
            .iter()
            .try_fold(0u64, |sum, a| sum.checked_add(conversion_degree(a)?))?,
        AtomView::Pow(p) => {
            let (base, exponent) = p.get_base_exp();
            let exponent = Rational::try_from(exponent).ok()?;
            if !exponent.is_integer() {
                return None;
            }
            let magnitude = exponent.numerator().to_i64()?.unsigned_abs();
            if magnitude > MAX {
                return None;
            }
            conversion_degree(base)?.checked_mul(magnitude)?
        }
        _ => return None,
    };
    (degree <= MAX).then_some(degree)
}

/// Normalize Q(symbols) independently of polynomial variable registration order.
/// Preserve structural identity for supported nonrational/complex Atom entries;
/// treating algebraic powers as independent polynomial symbols would change it.
pub(super) fn coefficient(value: &Atom) -> Result<String> {
    let structural = || {
        serde_json::to_string(&("atom-v1", value.to_canonical_string()))
            .map_err(|e| Error::Cache(e.to_string()))
    };
    if conversion_degree(value.as_view()).is_none() {
        return structural();
    }
    let Ok(fraction): std::result::Result<RationalPolynomial<IntegerRing, u32>, _> =
        value.try_to_rational_polynomial(&Q, &Z, None)
    else {
        return structural();
    };
    if fraction
        .numerator
        .variables()
        .iter()
        .any(|v| !matches!(v, PolyVariable::Symbol(_)))
    {
        return structural();
    }
    let variables = fraction.numerator.variables();
    let mut order = (0..variables.len()).collect::<Vec<_>>();
    order.sort_by_cached_key(|&i| match &variables[i] {
        PolyVariable::Symbol(s) => Atom::var(*s).to_canonical_string(),
        _ => unreachable!("non-symbol variables retained structurally above"),
    });
    let (numerator, denominator) = if order.is_empty() {
        (fraction.numerator, fraction.denominator)
    } else {
        (
            fraction.numerator.rearrange(&order, false),
            fraction.denominator.rearrange(&order, false),
        )
    };
    // This public constructor cancels the polynomial/content gcd and fixes the
    // denominator's leading sign in the selected qualified-name variable order.
    let fraction = RationalPolynomial::from_num_den(numerator, denominator, &Z, true);
    serde_json::to_string(&(
        "rational-v1",
        fraction.numerator.to_expression().to_canonical_string(),
        fraction.denominator.to_expression().to_canonical_string(),
    ))
    .map_err(|e| Error::Cache(e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::process::Command;

    fn parse(s: &str) -> Atom {
        Atom::parse(
            s,
            "dense_hash_regression",
            symbolica::parser::ParseSettings::default(),
        )
        .unwrap()
    }
    fn identity(
        value: Atom,
        domain: &str,
        prescription: Prescription,
        guards: &[Atom],
    ) -> BoundaryIdentity {
        let x = symbol!("dense_hash_regression::x");
        BoundaryIdentity::with_conditions(
            &KinematicSystem {
                epsilon: symbol!("dense_hash_regression::epsilon"),
                derivatives: BTreeMap::from([(x, vec![vec![value]])]),
            },
            &[parse("I(1)")],
            &Atom::one(),
            prescription,
            domain,
            guards,
        )
        .unwrap()
    }

    #[test]
    fn actual_gl080_sign_unit_pair() {
        let old = parse(ACTUAL_OLD);
        let new = parse(ACTUAL_NEW);
        assert!((&old - &new).together().cancel().is_zero());
        assert_eq!(coefficient(&old).unwrap(), coefficient(&new).unwrap());
        assert_ne!(
            coefficient(&old).unwrap(),
            coefficient(&(&new + Atom::one())).unwrap()
        );
        let mut derivatives = BTreeMap::new();
        for name in ["A", "C", "s"] {
            derivatives.insert(
                symbol!(&format!("closed_top_open_box_unit::{name}")),
                vec![vec![Atom::zero()]],
            );
        }
        let variable = symbol!("closed_top_open_box_unit::s");
        derivatives.get_mut(&variable).unwrap()[0][0] = old;
        let a = IdentitySystem::Dense(KinematicSystem {
            epsilon: symbol!("closed_top_open_box::epsilon"),
            derivatives: derivatives.clone(),
        })
        .canonical_key()
        .unwrap();
        derivatives.get_mut(&variable).unwrap()[0][0] = new;
        let b = IdentitySystem::Dense(KinematicSystem {
            epsilon: symbol!("closed_top_open_box::epsilon"),
            derivatives,
        })
        .canonical_key()
        .unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn expanded_factored_sign_and_content_forms() {
        let forms = [
            "(x+y)*(x-y)/((x+1)*(y+2))",
            "(x^2-y^2)/(x*y+2*x+y+2)",
            "-(-x^2+y^2)/(x*y+2*x+y+2)",
            "(7*x^2-7*y^2)/(7*x*y+14*x+7*y+14)",
            "(x^2-y^2)*(x+3)/((x*y+2*x+y+2)*(x+3))",
        ];
        let expected = coefficient(&parse(forms[0])).unwrap();
        for form in forms {
            assert_eq!(expected, coefficient(&parse(form)).unwrap(), "{form}");
        }
        for constant in ["0", "1", "-1", "2/3", "-4/6"] {
            assert!(coefficient(&parse(constant)).is_ok());
        }
        assert_eq!(
            coefficient(&parse("2/3")).unwrap(),
            coefficient(&parse("4/6")).unwrap()
        );
    }

    #[test]
    fn genuine_change_does_not_merge() {
        let value = parse("(x+1)/(x+2)");
        assert_ne!(
            coefficient(&value).unwrap(),
            coefficient(&(&value + Atom::one())).unwrap()
        );
        assert_ne!(
            identity(value.clone(), "domain", Prescription::PlusI0, &[]).key(),
            identity(value + Atom::one(), "domain", Prescription::PlusI0, &[]).key()
        );
    }

    #[test]
    fn identity_metadata_remains_distinct() {
        let value = parse("x+1");
        let key = identity(value.clone(), "domain", Prescription::PlusI0, &[])
            .key()
            .to_owned();
        assert_ne!(
            key,
            identity(value.clone(), "other-domain", Prescription::PlusI0, &[]).key()
        );
        assert_ne!(
            key,
            identity(value.clone(), "domain", Prescription::MinusI0, &[]).key()
        );
        assert_ne!(
            key,
            identity(
                value.clone(),
                "domain",
                Prescription::PlusI0,
                &[parse("x+2")]
            )
            .key()
        );
        let x = symbol!("dense_hash_regression::x");
        let epsilon = symbol!("dense_hash_regression::epsilon");
        let system = KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(x, vec![vec![value]])]),
        };
        for (basis, normalization) in [(parse("I(2)"), Atom::one()), (parse("I(1)"), Atom::num(2))]
        {
            assert_ne!(
                key,
                BoundaryIdentity::new(
                    &system,
                    &[basis],
                    &normalization,
                    Prescription::PlusI0,
                    "domain"
                )
                .unwrap()
                .key()
            );
        }
        let changed = KinematicSystem {
            epsilon: symbol!("dense_hash_regression::other_epsilon"),
            derivatives: system.derivatives,
        };
        assert_ne!(
            key,
            BoundaryIdentity::new(
                &changed,
                &[parse("I(1)")],
                &Atom::one(),
                Prescription::PlusI0,
                "domain"
            )
            .unwrap()
            .key()
        );
    }

    #[test]
    fn supported_algebraic_and_complex_entries_keep_structural_identity() {
        for value in [parse("x^(1/2)"), Atom::i() * parse("x")] {
            let serialized = coefficient(&value).unwrap();
            let fields: serde_json::Value = serde_json::from_str(&serialized).unwrap();
            assert_eq!(fields[0], "atom-v1");
            assert_eq!(fields[1], value.to_canonical_string());
            assert!(identity(value, "domain", Prescription::PlusI0, &[]).dimension() == 1);
        }
    }

    #[test]
    fn registered_roots_remain_distinct() {
        let x = symbol!("dense_hash_regression::x");
        let root = symbol!("dense_hash_regression::root");
        let system = |radicand| AlgebraicKinematicSystem {
            epsilon: symbol!("dense_hash_regression::epsilon"),
            derivatives: BTreeMap::from([(x, vec![vec![Atom::var(root)]])]),
            roots: vec![SquareRoot {
                symbol: root,
                radicand,
            }],
        };
        let a = BoundaryIdentity::with_algebraic_conditions(
            &system(Atom::var(x)),
            &[parse("I(1)")],
            &Atom::one(),
            Prescription::PlusI0,
            "domain",
            &[],
        )
        .unwrap();
        let b = BoundaryIdentity::with_algebraic_conditions(
            &system(Atom::var(x) + Atom::one()),
            &[parse("I(1)")],
            &Atom::one(),
            Prescription::PlusI0,
            "domain",
            &[],
        )
        .unwrap();
        assert_ne!(a.key(), b.key());
    }

    #[test]
    fn integer_exponent_conversion_boundaries() {
        for expression in ["x^2147483647", "x^(-2147483647)"] {
            let value = parse(expression);
            assert_eq!(conversion_degree(value.as_view()), Some(i32::MAX as u64));
            let fields: serde_json::Value =
                serde_json::from_str(&coefficient(&value).unwrap()).unwrap();
            assert_eq!(fields[0], "rational-v1");
        }
        for expression in [
            "x^2147483648",
            "x^(-2147483648)",
            "x^4294967296",
            "x^(-4294967296)",
            "x^9223372036854775808",
            "x^(-9223372036854775808)",
            "(x^1073741824+1)^2",
        ] {
            let value = parse(expression);
            assert_eq!(conversion_degree(value.as_view()), None, "{expression}");
            let fields: serde_json::Value =
                serde_json::from_str(&coefficient(&value).unwrap()).unwrap();
            assert_eq!(fields[0], "atom-v1");
            assert_eq!(fields[1], value.to_canonical_string());
        }
    }

    #[test]
    fn previously_supported_large_power_identity_stays_supported() {
        let value = parse("x^2147483648");
        let original = identity(value.clone(), "large power", Prescription::PlusI0, &[]);
        let changed = identity(
            value + Atom::one(),
            "large power",
            Prescription::PlusI0,
            &[],
        );
        assert_ne!(original.key(), changed.key());
        assert_eq!(original.dimension(), 1);
    }

    #[test]
    fn variable_registration_order() {
        const ENV: &str = "RUSTFLOW_DENSE_HASH_REGISTRATION_TEST";
        if let Ok(order) = std::env::var(ENV) {
            let names = if order == "az" {
                ["a", "z"]
            } else {
                ["z", "a"]
            };
            for name in names {
                symbol!(&format!("dense_hash_registration::{name}"));
            }
            let value = Atom::parse(
                "(a+z)/(a-z)",
                "dense_hash_registration",
                symbolica::parser::ParseSettings::default(),
            )
            .unwrap();
            println!("DENSE_HASH_KEY={}", coefficient(&value).unwrap());
            return;
        }
        let mut keys = Vec::new();
        for order in ["az", "za"] {
            let output = Command::new(std::env::current_exe().unwrap())
                .args([
                    "--exact",
                    "transport_cache::dense_key::tests::variable_registration_order",
                    "--nocapture",
                    "--test-threads=1",
                ])
                .env(ENV, order)
                .env("SYMBOLICA_HIDE_BANNER", "1")
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "isolated registration check failed: {}",
                output.status
            );
            let stdout = String::from_utf8(output.stdout).unwrap();
            // libtest can print its test-name prefix before the child's marker.
            const MARKER: &str = "DENSE_HASH_KEY=";
            assert_eq!(
                stdout.matches(MARKER).count(),
                1,
                "expected exactly one child key marker"
            );
            let key = stdout.split_once(MARKER).unwrap().1.lines().next().unwrap();
            println!("DENSE_HASH_REGISTRATION order={order} key={key}");
            keys.push(key.to_owned());
        }
        assert_eq!(keys[0], keys[1]);
    }

    // Exact strings from the two saved GL080 contracts; no runtime input files.
    const ACTUAL_OLD: &str = "(-1*closed_top_open_box_unit::{}::s^2+-2*closed_top_open_box::{}::epsilon+-2*closed_top_open_box::{}::epsilon*closed_top_open_box_unit::{}::A^2+2*closed_top_open_box::{}::epsilon*closed_top_open_box_unit::{}::A*closed_top_open_box_unit::{}::s+2*closed_top_open_box::{}::epsilon*closed_top_open_box_unit::{}::s+4*closed_top_open_box::{}::epsilon*closed_top_open_box_unit::{}::A+closed_top_open_box_unit::{}::A*closed_top_open_box_unit::{}::s+closed_top_open_box_unit::{}::s)*(-2*closed_top_open_box_unit::{}::A+-2*closed_top_open_box_unit::{}::A*closed_top_open_box_unit::{}::s+-2*closed_top_open_box_unit::{}::s+1+closed_top_open_box_unit::{}::A^2+closed_top_open_box_unit::{}::s^2)^(-1)*closed_top_open_box_unit::{}::s^(-1)";
    const ACTUAL_NEW: &str = "(-1*closed_top_open_box_unit::{}::A*closed_top_open_box_unit::{}::s+-1*closed_top_open_box_unit::{}::s+-2*closed_top_open_box::{}::epsilon*closed_top_open_box_unit::{}::A*closed_top_open_box_unit::{}::s+-2*closed_top_open_box::{}::epsilon*closed_top_open_box_unit::{}::s+-4*closed_top_open_box::{}::epsilon*closed_top_open_box_unit::{}::A+2*closed_top_open_box::{}::epsilon+2*closed_top_open_box::{}::epsilon*closed_top_open_box_unit::{}::A^2+closed_top_open_box_unit::{}::s^2)*(-2*closed_top_open_box_unit::{}::A+-2*closed_top_open_box_unit::{}::A*closed_top_open_box_unit::{}::s+-2*closed_top_open_box_unit::{}::s+1+closed_top_open_box_unit::{}::A^2+closed_top_open_box_unit::{}::s^2)^(-1)*-1*closed_top_open_box_unit::{}::s^(-1)";
}
