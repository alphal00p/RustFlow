//! Symbolic counterpart of boundary.rs compact polynomial moments.
//! Returns coefficient after removing one A(D)*mu^(D-2) per shell,
//! A(D)=1/(sqrt(pi)*Gamma((D-1)/2)); no period is supplied here.
use super::*;
use symbolica_amflow::integrand::RegionFactorSpace;
use symbolica_amflow::finite_density::boundary::OccupiedBoundaryDistribution;
use symbolica_amflow::family::substitute;
const BUDGET:usize=200000;
fn poly(a:&Atom,vars:&[Atom])->Result<Vec<(Vec<i16>,Atom)>>{
 let terms=coefficient::exact_coefficient_list(a,vars)?;
 if terms.len()>BUDGET{return Err(Error::Limit("moment monomial budget".into()))}
 terms.into_iter().map(|(m,c)|{let p=coefficient::powers(&m,vars)?;if p.iter().any(|&p|p<0)||vars.iter().any(|v|{let AtomView::Var(s)=v.as_view() else{return true};!c.derivative(s.get_symbol()).is_zero()}){return Err(Error::Unsupported("residual compact or virtual-soft denominator".into()))}Ok((p,c))}).collect()
}
pub fn radial(d:&Atom,n:i16,s:i16,l:i16,r:i16,j:i16,mu:&Rational)->Result<Atom>{
 if n<=0{return Ok(Atom::zero())}if s<0||l<0||r<0||j<0||mu<=&Rational::zero(){return Err(Error::InvalidInput("massless polynomial moment roles/support".into()))}
 if n>32||s>32||l>32{return Err(Error::Limit("moment jet cap32".into()))}if l>0{return Ok(Atom::zero())}
 let mut choose=Atom::one();for k in 0..n-1{choose*= (d/Atom::num(2)+Atom::num(j-1-k))/Atom::num(k+1)}
 let beta=d+Atom::num(r+2*j+1-2*n);let offset=i64::from(r+2*j+2-2*n);
 let mut result=Atom::num(if (n-1)%2==0{1}else{-1})*choose;
 if s==0{result=result/beta*Atom::num(mu.clone()).pow(offset)}else{
  for k in 0..s-1{result*=(-&beta+Atom::num(1+k))/Atom::num(k+1)}
  result*=Atom::num(mu.clone()).pow(offset-i64::from(s));
 }
 Ok(result.together().cancel())
}
pub fn integrate(expression:&Atom,space:&RegionFactorSpace,distributions:&[OccupiedBoundaryDistribution])->Result<Atom>{
 let loops=space.template.loops.len();if loops!=distributions.len()||space.template.external.len()!=1||space.template.external_gram!=vec![vec![Atom::one()]]{return Err(Error::Unsupported("compact factor geometry".into()))}
 let pairs=(0..loops).flat_map(|i|(i..loops).map(move|j|(i,j))).collect::<Vec<_>>();
 if space.coordinates.len()!=pairs.len()+loops{return Err(Error::InvalidInput("moment coordinate count".into()))}poly(expression,&space.coordinates)?;
 let es=(0..loops).map(|i|Atom::var(symbol!(format!("integrated_moment::E_{i}")))).collect::<Vec<_>>();let gs=pairs.iter().map(|&(i,j)|Atom::var(symbol!(format!("integrated_moment::R_{i}_{j}")))).collect::<Vec<_>>();
 let mut rules=BTreeMap::new();for (k,&(i,j)) in pairs.iter().enumerate(){rules.insert(space.coordinates[k].clone(),&es[i]*&es[j]-&gs[k]);}for i in 0..loops{rules.insert(space.coordinates[pairs.len()+i].clone(),es[i].clone());}
 let mut p=substitute(expression,&rules).expand();let d=Atom::num(space.template.dimension-1)-Atom::num(2)*Atom::var(space.template.epsilon);let mut projector=tensor::TensorProjector::new(d.clone());
 let scalar=|i:usize,j:usize|gs[pairs.iter().position(|&x|x==(i.min(j),i.max(j))).unwrap()].clone();
 for selected in 0..loops{let mut next=Atom::zero();for (powers,c) in poly(&p,&gs)?{let mut partners=Vec::new();let mut unmixed=c;for (k,(&power,&(a,b))) in powers.iter().zip(&pairs).enumerate(){if a==selected&&b>selected{if partners.len()+power as usize>32{return Err(Error::Limit("compact tensor rank32".into()))}partners.extend(std::iter::repeat_n(b,power as usize));}else{unmixed*=gs[k].pow(i64::from(power));}}
 let h=vec![vec![scalar(selected,selected);partners.len()];partners.len()];let s=partners.iter().map(|&i|partners.iter().map(|&j|scalar(i,j)).collect()).collect::<Vec<_>>();next+=unmixed*projector.project(&h,&s)?;}p=next.expand();}
 let vars=es.iter().chain(&gs).cloned().collect::<Vec<_>>();let mut result=Atom::zero();for(powers,c)in poly(&p,&vars)?{let mut v=c;for (i,dist)in distributions.iter().enumerate(){if !dist.shell.mass_squared.is_zero(){return Err(Error::Unsupported("symbolic adapter massless scope".into()))}let j=pairs.iter().position(|&p|p==(i,i)).unwrap();v*=radial(&d,dist.cut_index,dist.upper_index,dist.lower_index,powers[i],powers[loops+j],&dist.shell.chemical_potential)?;}for(k,&(i,j))in pairs.iter().enumerate(){if i!=j&&powers[loops+k]!=0{return Err(Error::InvalidInput("unprojected spatial contraction".into()))}}result+=v;}
 Ok(result.together().cancel())
}
pub fn validate_polynomial(expression:&Atom,space:&RegionFactorSpace)->Result<usize>{Ok(poly(expression,&space.coordinates)?.len())}
