use symbolica_amflow::*;
use symbolica_amflow::algebraic::{AlgebraicSystem,SquareRoot};
use symbolica_amflow::symbolica::prelude::*;
use symbolica_amflow::symbolica::domains::float::{FloatField,ComplexBall,RealBall};
use symbolica_amflow::symbolica::tensors::matrix::Matrix;
fn main() {
 let p=Precision::decimal(55).unwrap();
 let system=AlgebraicSystem::ordinary(DifferentialSystem { variable:symbol!("z"),matrix:vec![vec![parse!("1/r")]] },vec![SquareRoot { symbol:symbol!("r"),radicand:parse!("z") }]);
 let prepared=system.prepare_frobenius(8,&RunContext::default()).unwrap();
 let basis=prepared.rational_preparation().evaluate(p,&Default::default(),32,&RunContext::default()).unwrap();
 for winding in [0,1,-1] {
  let matrix=Matrix::from_nested_vec(basis.evaluate_with_winding(&p.rational(&Rational::from((1,16))),&Default::default(),winding).unwrap(),FloatField::from_rep(p.zero())).unwrap();
  let inverse=matrix.inv().unwrap();
  println!("winding {winding} matrix {matrix:?} inverse {inverse:?}");
  let ball=|a:&ComplexFloat| { let mut value=ComplexBall::new(RealBall::exact(a.re.clone()),RealBall::exact(a.im.clone())); value.set_precision(p.bits); value };
  let field=FloatField::from_rep(ball(&p.zero()));
  let residual=&Matrix::identity(2,field.clone())-&(&matrix.map(ball,field.clone())*&inverse.map(ball,field));
  println!("residual {residual:?}");
 }
}
