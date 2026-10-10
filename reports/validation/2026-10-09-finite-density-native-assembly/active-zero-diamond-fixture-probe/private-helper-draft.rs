use std::sync::Arc;
use rustred::algebra::CoefficientContext;
use rustred::solver::{Integral,IntegralOrder,Seed,Case,RuleCandidate,SearchStats,SeedSource,Term};
use rustred::solver::guarded::{GuardedSourceSystem,GuardedSource,GuardedProgram,IndexRole,IndexDomain,IndexBounds,GuardedReductionLimits};
    fn guarded_pending_diamond(cancel: bool, canonical_rhs: bool) -> GuardedProgram<1> {
        let context = CoefficientContext::new(["pending_n", "pending_x", "pending_y"]);
        let x = context.parameter("pending_x").unwrap();
        let y = context.parameter("pending_y").unwrap();
        let target = Integral::symbolic([0]).unwrap();
        let branches = [
            (
                6,
                context.one(),
                if canonical_rhs { vec![(-1, context.one()), (-4, context.one())] } else { vec![(-4, context.one()), (-1, context.one())] },
            ),
            (5, context.one(), vec![(-1, context.one())]),
            (
                4,
                context.one(),
                vec![(
                    -2,
                    if cancel {
                        -context.one()
                    } else {
                        context.one()
                    },
                )],
            ),
            (2, y.clone(), vec![(-1, &context.one() / &y)]),
        ];
        let sources = Arc::new(
            GuardedSourceSystem::new(
                "pending-native-order-guarded-diamond",
                [IndexRole::Occupation],
                [0],
                branches
                    .iter()
                    .enumerate()
                    .map(|(ordinal, (point, pivot, rhs))| {
                        let mut row = vec![Term {
                            integral: target,
                            coefficient: pivot.numerator.clone(),
                        }];
                        row.extend(rhs.iter().map(|(shift, coefficient)| Term {
                            integral: Integral::symbolic([*shift]).unwrap(),
                            coefficient: (-(coefficient * pivot)).numerator,
                        }));
                        let source = GuardedSource::new(
                            format!("pending-{ordinal}"),
                            row,
                            IndexDomain::new([IndexBounds::fixed(*point)]).unwrap(),
                        );
                        if *point == 5 {
                            source.with_nonzero_conditions(vec![x.numerator.clone()])
                        } else {
                            source
                        }
                    })
                    .collect(),
            )
            .unwrap(),
        );
        let rules = branches
            .into_iter()
            .enumerate()
            .map(|(ordinal, (point, _, rhs))| {
                sources
                    .seal_candidate(
                        RuleCandidate {
                            case: Case::generic(),
                            target,
                            rhs: rhs
                                .into_iter()
                                .map(|(shift, coefficient)| Term {
                                    integral: Integral::symbolic([shift]).unwrap(),
                                    coefficient,
                                })
                                .collect(),
                            sources: vec![SeedSource {
                                basis_row: ordinal,
                                seed: Seed {
                                    integral: target,
                                    shifts: [0],
                                },
                            }],
                            stats: SearchStats::default(),
                        },
                        IntegralOrder::new([true], [false])
                            .with_roles([IndexRole::Occupation])
                            .unwrap(),
                        IndexDomain::new([IndexBounds::fixed(point)]).unwrap(),
                    )
                    .unwrap()
            })
            .collect();
        GuardedProgram::new(sources, rules, [[1]]).unwrap()
    }

fn main(){
for cancel in [true,false] {
 let original=std::panic::catch_unwind(||guarded_pending_diamond(cancel,false));assert!(original.is_err());
 let p=guarded_pending_diamond(cancel,true);let r=p.reduce([6],GuardedReductionLimits{max_rule_applications:4,max_pending_integrals:10}).unwrap();
 assert!(r.unresolved.is_empty());assert_eq!(r.rule_applications,if cancel{3}else{4});
 let c=CoefficientContext::new(["pending_n","pending_x","pending_y"]);let x=c.parameter("pending_x").unwrap();let y=c.parameter("pending_y").unwrap();assert!(r.nonzero_conditions.contains(&x.numerator));
 if cancel{assert!(r.terms.is_empty());assert!(!r.nonzero_conditions.contains(&y.numerator));}else{assert_eq!(r.terms[&[1]],&(&c.one()+&c.one())/&y);assert!(r.nonzero_conditions.contains(&y.numerator));}
 println!("cancel={cancel} original_rhs_rejected=true canonical_rhs_replayed=true applications={} conditions={:?}",r.rule_applications,r.nonzero_conditions.iter().map(|v|v.to_expression().to_string()).collect::<Vec<_>>());
}
}