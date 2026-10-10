use rustred::solver::guarded::{GuardedProgram,GuardedSourceSystem};
use std::sync::Arc;
fn forge(s: Arc<GuardedSourceSystem<1>>) -> GuardedProgram<1> { GuardedProgram { sources:s, rules:vec![], terminals:Default::default() } }
