//! Private importer experiment: preserve original UFO momentum-denominator scopes.
//! Token owns parsing; Spenso/Idenso own tensor inference and contractions.
use symbolica::{
    parser::{Operator, ParseSettings, Token},
    prelude::*,
};

fn contains_momentum(token: &Token) -> bool {
    match token {
        Token::Fn(_, _, args) => {
            matches!(args.first(), Some(Token::ID(s)) if s.as_str().rsplit("::").next()==Some("P"))
                || args.iter().skip(1).any(contains_momentum)
        }
        Token::Op(_, _, _, args) => args.iter().any(contains_momentum),
        _ => false,
    }
}

fn scoped(token: &Token) -> bool {
    matches!(token, Token::Fn(_, _, args) if matches!(args.first(), Some(Token::ID(s)) if s.as_str()=="spenso::bracket"))
}

fn visit(token: &mut Token) -> Result<usize, String> {
    let mut changed = 0;
    match token {
        Token::Fn(_, _, args) => {
            // An existing native scope is deliberate and must remain intact.
            if matches!(args.first(), Some(Token::ID(s)) if s.as_str()=="spenso::bracket") {
                return Ok(0);
            }
            for arg in args.iter_mut().skip(1) {
                changed += visit(arg)?;
            }
        }
        Token::Op(_, _, operator, args) => {
            for arg in args.iter_mut() {
                changed += visit(arg)?;
            }
            let reciprocal = match operator {
                Operator::Inv => args.len() == 1,
                Operator::Pow if args.len() == 2 => {
                    let power = Atom::parse(args[1].to_string(), "UFO", ParseSettings::default())?;
                    Rational::try_from(power.as_view()).is_ok_and(|r| r.is_negative())
                }
                _ => false,
            };
            if reciprocal && contains_momentum(&args[0]) && !scoped(&args[0]) {
                args[0] = Token::Fn(
                    false,
                    false,
                    vec![Token::ID("spenso::bracket".into()), args[0].clone()],
                );
                changed += 1;
            }
        }
        _ => {}
    }
    Ok(changed)
}

pub fn preserve_momentum_denominators(source: &str) -> Result<(String, usize), String> {
    let mut token = Token::parse(
        &source.replace("**", "^"),
        ParseSettings::default().convert_mul_to_atom(false),
    )?;
    let count = visit(&mut token)?;
    Ok((token.to_string(), count))
}

#[cfg(test)]
mod tests {
    use super::*;
    use idenso::{
        IndexTooling,
        tensor::{ContractSettings, SymbolicTensor},
    };
    use spenso::{
        network::tags::SPENSO_TAG,
        structure::{abstract_index::AbstractIndex, partial::PartialStructure},
    };

    fn lower(source: &str) -> Atom {
        let p = SPENSO_TAG.rank_one_tensor_symbol("scope_tests::P");
        let raw = Atom::parse(source, "scope_tests", ParseSettings::default()).unwrap();
        let _ = p;
        raw
    }

    #[test]
    fn closed_denominator_retains_external_port_and_native_dot() {
        let source = "P(3,spenso::mink(4,nu))/(P(1,spenso::mink(4,mu))*P(2,spenso::mink(4,mu)))";
        assert!(lower(source).list_dangling::<AbstractIndex>().is_err());
        let (converted, n) = preserve_momentum_denominators(source).unwrap();
        assert_eq!(n, 1);
        let tensor = SymbolicTensor::<PartialStructure>::infer(lower(&converted)).unwrap();
        let compact = tensor
            .contract(ContractSettings::default())
            .unwrap()
            .to_dots()
            .unwrap();
        assert_eq!(
            compact
                .expression()
                .list_dangling::<AbstractIndex>()
                .unwrap()
                .len(),
            1
        );
        let expected =
            lower("P(3,spenso::mink(4,nu))/spenso::dot(P(1,spenso::mink(4)),P(2,spenso::mink(4)))");
        assert_eq!(compact.expression(), &expected);
    }

    #[test]
    fn negative_powers_nested_division_and_existing_scope() {
        for source in [
            "(P(1,spenso::mink(4,mu))*P(2,spenso::mink(4,mu)))^(-2)",
            "1/(2+1/(P(1,spenso::mink(4,mu))*P(2,spenso::mink(4,mu))))",
        ] {
            let (converted, n) = preserve_momentum_denominators(source).unwrap();
            assert!(n > 0);
            assert!(
                SymbolicTensor::<PartialStructure>::infer(lower(&converted))
                    .unwrap()
                    .is_scalar()
            );
            assert_eq!(preserve_momentum_denominators(&converted).unwrap().1, 0);
        }
    }

    #[test]
    fn open_denominator_stays_invalid_and_scalar_algebra_unchanged() {
        let (converted, n) = preserve_momentum_denominators("1/P(1,spenso::mink(4,mu))").unwrap();
        assert_eq!(n, 1);
        assert!(SymbolicTensor::<PartialStructure>::infer(lower(&converted)).is_err());
        let source = "(a+b)/(m^2+s)^2";
        let (converted, n) = preserve_momentum_denominators(source).unwrap();
        assert_eq!(n, 0);
        assert_eq!(lower(&converted), lower(source));
    }
}
