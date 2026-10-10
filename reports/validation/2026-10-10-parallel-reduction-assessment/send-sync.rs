//! Compile-time trait gate only. No reductions or worker threads execute.
fn shared<T: Send + Sync>() {}
fn evidence() {
 shared::<rustred::solver::guarded::GuardedProgram<16>>();
 shared::<rustred::solver::guarded::GuardedReduction<16>>();
 shared::<rustred::algebra::Coefficient>();
 shared::<rustred::algebra::CoefficientPolynomial>();
 shared::<symbolica_amflow::finite_density::guarded::GuardedReductionProgram<16>>();
 shared::<symbolica_amflow::finite_density::guarded::GuardedContext<16>>();
 shared::<symbolica_amflow::RunContext>();
 shared::<symbolica::prelude::Atom>();
}
