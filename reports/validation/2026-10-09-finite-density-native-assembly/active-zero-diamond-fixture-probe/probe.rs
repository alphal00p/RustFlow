use std::sync::Arc;
use rustred::algebra::CoefficientContext;
use rustred::solver::{Integral,IntegralOrder,Seed,Case,RuleCandidate,SearchStats,SeedSource,Term};
use rustred::solver::guarded::{GuardedSourceSystem,GuardedSource,GuardedProgram,IndexRole,IndexDomain,IndexBounds,GuardedReductionLimits};
    fn guarded_pending_diamond(cancel: bool, canonical_rhs: bool) -> GuardedProgram<1> {
        let context = CoefficientContext::try_new(["pending_n", "pending_x", "pending_y"]).unwrap();
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
        let domains=branches.iter().map(|(point,_,_)|IndexDomain::new([IndexBounds::fixed(*point)]).unwrap()).collect::<Vec<_>>();
        let solution=sources.solve_domains(domains,rustred::solver::SearchOptions{max_depth:Some(2),..Default::default()},32).unwrap();
        for r in &solution.rules { println!("discovered target={:?} rhs={:?}",r.candidate().target,r.candidate().rhs.iter().map(|t|t.integral).collect::<Vec<_>>()); }
        let rules=solution.rules;
        GuardedProgram::new(sources, rules, [[1]]).unwrap()
    }

fn main(){
for cancel in [true,false] {

 let p=guarded_pending_diamond(cancel,true);let r=p.reduce([6],GuardedReductionLimits{max_rule_applications:4,max_pending_integrals:10}).unwrap();
 assert!(r.unresolved.is_empty());assert_eq!(r.rule_applications,if cancel{3}else{4});
 let c=CoefficientContext::try_new(["pending_n","pending_x","pending_y"]).unwrap();let x=c.parameter("pending_x").unwrap();let y=c.parameter("pending_y").unwrap();assert!(r.nonzero_conditions.contains(&x.numerator));
 if cancel{assert!(r.terms.is_empty());assert!(!r.nonzero_conditions.contains(&y.numerator));}else{assert_eq!(r.terms[&[1]],&(&c.one()+&c.one())/&y);assert!(r.nonzero_conditions.contains(&y.numerator));}
 println!("cancel={cancel} public_native_discovery_replayed=true applications={} conditions={:?}",r.rule_applications,r.nonzero_conditions.iter().map(|v|v.to_expression().to_string()).collect::<Vec<_>>());
}
}