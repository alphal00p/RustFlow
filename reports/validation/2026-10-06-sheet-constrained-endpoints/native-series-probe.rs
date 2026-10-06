use symbolica::prelude::*;
fn main() {
 let z=symbol!("sheet_probe::z");
 for expr in ["(1-z)^(-1/2)","(z-z^2)*(1-z)^(-1/2)"] {
  let a=Atom::parse(expr,"sheet_probe",Default::default()).unwrap();
  let series=a.series(z,0,17).unwrap();
  println!("expression={a}");
  for (power, coefficient) in series.terms() { println!("power={power} coefficient={coefficient}"); }
 }
}
