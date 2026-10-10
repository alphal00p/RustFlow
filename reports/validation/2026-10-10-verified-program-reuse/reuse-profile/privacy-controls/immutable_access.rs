use rustred::solver::guarded::{GuardedProgram,GuardedSourceSystem};
use std::sync::Arc;
fn inspect(p: &GuardedProgram<1>) { let _ = (p.sources(),p.rules(),p.terminals()); }
