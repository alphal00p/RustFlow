use rustred::solver::guarded::{GuardedProgram,GuardedSourceSystem};
use std::sync::Arc;
fn mutate(p: &mut GuardedProgram<1>) { p.terminals().clear(); }
