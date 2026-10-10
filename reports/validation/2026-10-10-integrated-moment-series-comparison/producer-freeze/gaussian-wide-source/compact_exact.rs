//! Generic independent compact polynomial functional. No selected invariant.
use super::*;
use std::collections::BTreeMap;
pub struct CompactMoments {
 k:usize,d:Atom,dist:Vec<(i16,i16,i16,Rational)>,pairs:Vec<(usize,usize)>,angular:angular::AngularAverager,cache:BTreeMap<Vec<i16>,Atom>,pub evaluated:usize,budget:usize,
}
impl CompactMoments{
 pub fn new(k:usize,d:Atom,dist:Vec<(i16,i16,i16,Rational)>,budget:usize)->Result<Self>{if dist.len()!=k{return Err(Error::InvalidInput("compact distribution arity".into()))}Ok(Self{k,d:d.clone(),dist,pairs:(0..k).flat_map(|i|(i..k).map(move|j|(i,j))).collect(),angular:angular::AngularAverager::new(k,d,budget).map_err(Error::Limit)?,cache:BTreeMap::new(),evaluated:0,budget})}
 pub fn integrate(&mut self,p:&Atom,coords:&[Atom])->Result<Atom>{let terms=poly(p,coords,self.budget)?;let mut out=Atom::zero();for(degrees,c)in terms{out+=c*self.monomial(&degrees)?;}Ok(out.together().cancel())}
 pub fn monomial(&mut self,p:&[i16])->Result<Atom>{
 if p.len()!=self.pairs.len()+self.k||p.iter().any(|&v|v<0){return Err(Error::InvalidInput("compact polynomial degree".into()))}if let Some(v)=self.cache.get(p){return Ok(v.clone())}
 let mut ds=self.dist.clone();for (j,&(a,b))in self.pairs.iter().enumerate(){if a==b{ds[a].0-=p[j];if ds[a].0<=0{return Ok(Atom::zero())}}}
 let mut terms=BTreeMap::from([((p[self.pairs.len()..].to_vec(),vec![0u16;self.pairs.len()]),Atom::one())]);
 for(j,&(a,b))in self.pairs.iter().enumerate(){if a==b||p[j]==0{continue}let mut next=BTreeMap::new();let mut choose=Atom::one();for s in 0..=p[j]{for((energy,spatial),c)in &terms{let mut e=energy.clone();e[a]=e[a].checked_add(p[j]-s).ok_or_else(||Error::Limit("compact energy degree".into()))?;e[b]=e[b].checked_add(p[j]-s).ok_or_else(||Error::Limit("compact energy degree".into()))?;let mut g=spatial.clone();g[j]=s as u16;*next.entry((e,g)).or_insert_with(Atom::zero)+=c*&choose*Atom::num(if s%2==0{1}else{-1});}if s<p[j]{choose=choose*Atom::num(p[j]-s)/Atom::num(s+1)}}if next.len()>self.budget{return Err(Error::Limit("compact expansion term cap".into()))}terms=next;}
 let mut out=Atom::zero();for((e,g),c)in terms{let(radial,a)=self.angular.monomial(&g).map_err(Error::Limit)?;if a.is_zero(){continue}let mut v=c*a;for i in 0..self.k{v*=moments::radial(&self.d,ds[i].0,ds[i].1,ds[i].2,e[i],i16::try_from(radial[i]).map_err(|_|Error::Limit("compact radial degree".into()))?,&ds[i].3)?;}out+=v;}
 let out=out.together().cancel();self.evaluated+=1;if self.cache.len()>=self.budget{return Err(Error::Limit("compact memo entries".into()))}self.cache.insert(p.to_vec(),out.clone());Ok(out)
 }
}
