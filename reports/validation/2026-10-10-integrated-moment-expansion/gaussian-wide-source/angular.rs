//! Independent unit-sphere Gram moments, with equivalent Wick pairings combined.
//! No preferred invariant or topology; explicit degree and retained-term limits.
use std::collections::BTreeMap;
use symbolica::prelude::*;

pub struct AngularAverager {
    loops: usize,
    d: Atom,
    slots: Vec<(usize, usize)>,
    pair_cache: BTreeMap<Vec<u16>, BTreeMap<Vec<u16>, Atom>>,
    average_cache: BTreeMap<Vec<u16>, Atom>,
    used: usize,
    budget: usize,
    max_total_degree: u32,
}

impl AngularAverager {
    pub fn new(loops: usize, d: Atom, budget: usize) -> Result<Self, String> {
        if loops == 0 || loops > 16 || budget == 0 {
            return Err("angular geometry or retained-term budget invalid".into());
        }
        let slots = (0..loops).flat_map(|i| (i..loops).map(move |j| (i,j))).collect();
        Ok(Self { loops, d, slots, pair_cache: BTreeMap::new(), average_cache: BTreeMap::new(), used: 0, budget, max_total_degree: 512 })
    }

    /// Set a tighter preflight for the host's polynomial degree budget.
    pub fn with_max_total_degree(mut self, degree: u32) -> Result<Self, String> {
        if degree == 0 || degree > 512 { return Err("angular total-degree limit must be in 1..=512".into()); }
        self.max_total_degree = degree; Ok(self)
    }

    fn retain(&mut self, count: usize) -> Result<(), String> {
        self.used = self.used.checked_add(count).ok_or("angular budget overflow")?;
        if self.used > self.budget { return Err("aggregated angular retained-term budget exhausted".into()); }
        Ok(())
    }

    pub fn retained_entries(&self) -> usize { self.used }

    /// Upper-triangular Gram powers. Returns powers of r_i^2 and an exact Q(d)
    /// coefficient. An odd degree in any independent vector gives zero.
    pub fn monomial(&mut self, powers: &[u16]) -> Result<(Vec<u16>, Atom), String> {
        if powers.len() != self.slots.len() { return Err("angular Gram arity mismatch".into()); }
        let total_degree = powers.iter().try_fold(0u32, |s,&n| s.checked_add(2*u32::from(n))).ok_or("angular total degree overflow")?;
        if total_degree > self.max_total_degree { return Err("angular total-degree preflight exceeded".into()); }
        let mut degree = vec![0u32; self.loops];
        let mut graph = powers.to_vec();
        for (s, &(i,j)) in self.slots.iter().enumerate() {
            degree[i] += u32::from(powers[s]); degree[j] += u32::from(powers[s]);
            if i == j { graph[s] = 0; }
        }
        if degree.iter().any(|x| x % 2 != 0) { return Ok((vec![0; self.loops], Atom::zero())); }
        let radial = degree.into_iter().map(|x| u16::try_from(x/2).map_err(|_| "radial exponent overflow".to_string())).collect::<Result<Vec<_>,_>>()?;
        Ok((radial, self.average_graph(&graph)?))
    }

    fn pairings(&mut self, counts: &[u16]) -> Result<BTreeMap<Vec<u16>, Atom>, String> {
        if let Some(value) = self.pair_cache.get(counts) { return Ok(value.clone()); }
        let mut result = BTreeMap::new();
        if let Some(i) = counts.iter().position(|&n| n != 0) {
            let mut rest = counts.to_vec(); rest[i] -= 1;
            for j in i..self.loops {
                let multiplicity = rest[j];
                if multiplicity == 0 { continue; }
                rest[j] -= 1;
                for (mut graph, weight) in self.pairings(&rest)? {
                    if i != j {
                        let s = self.slots.iter().position(|&pair| pair == (i,j)).unwrap();
                        graph[s] = graph[s].checked_add(1).ok_or("pairing exponent overflow")?;
                    }
                    if !result.contains_key(&graph) && result.len() >= self.budget.saturating_sub(self.used) {
                        return Err("aggregated angular intermediate-term budget exhausted".into());
                    }
                    let value = result.entry(graph).or_insert_with(Atom::zero);
                    *value += Atom::num(i64::from(multiplicity)) * weight;
                }
                rest[j] += 1;
            }
        } else { result.insert(vec![0; self.slots.len()], Atom::one()); }
        self.retain(result.len())?;
        self.pair_cache.insert(counts.to_vec(), result.clone());
        Ok(result)
    }

    fn average_graph(&mut self, graph: &[u16]) -> Result<Atom, String> {
        if let Some(value) = self.average_cache.get(graph) { return Ok(value.clone()); }
        let Some(i) = self.slots.iter().enumerate().find_map(|(s,&(i,j))| (i != j && graph[s] != 0).then_some(i)) else { return Ok(Atom::one()); };
        let mut counts = vec![0u16; self.loops];
        let mut rest = graph.to_vec();
        for (s,&(a,b)) in self.slots.iter().enumerate() {
            if a == i && b != i { counts[b] = graph[s]; rest[s] = 0; }
        }
        let rank: u32 = counts.iter().map(|&x| u32::from(x)).sum();
        if rank % 2 != 0 { return Ok(Atom::zero()); }
        let mut value = Atom::zero();
        for (added, weight) in self.pairings(&counts)? {
            let next = rest.iter().zip(added).map(|(&a,b)| a.checked_add(b).ok_or("angular graph exponent overflow".to_string())).collect::<Result<Vec<_>,_>>()?;
            value += weight * self.average_graph(&next)?;
        }
        let denominator = (0..rank/2).fold(Atom::one(), |a,j| a * (&self.d + Atom::num(i64::from(2*j))));
        value = (value / denominator).together().cancel();
        self.retain(1)?;
        self.average_cache.insert(graph.to_vec(), value.clone());
        Ok(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn equal(a: Atom,b: Atom) { assert!((a-b).together().cancel().is_zero()); }
    #[test]
    fn independent_parity_and_radial_factors() {
        let mut a=AngularAverager::new(3,parse!("d"),200000).unwrap();
        assert!(a.monomial(&[2,1,0,0,0,0]).unwrap().1.is_zero());
        let (radial,c)=a.monomial(&[3,0,0,4,0,5]).unwrap();
        assert_eq!(radial,vec![3,4,5]); equal(c,Atom::one());
    }
    #[test]
    fn coupled_triangles_and_four_vector_cycle() {
        let d=parse!("d");
        let mut a=AngularAverager::new(3,d.clone(),200000).unwrap();
        equal(a.monomial(&[0,1,1,0,1,0]).unwrap().1,Atom::one()/d.pow(2));
        equal(a.monomial(&[0,2,2,0,2,0]).unwrap().1,parse!("(d+8)/(d^2*(d+2)^2)"));
        let mut b=AngularAverager::new(4,parse!("d"),200000).unwrap();
        // g01*g12*g23*g03, with four independent compact directions.
        equal(b.monomial(&[0,1,0,1,0,1,0,0,1,0]).unwrap().1,parse!("1/d^3"));
    }
    #[test]
    fn high_rank_without_pairing_enumeration() {
        let mut a=AngularAverager::new(2,parse!("11/2"),200000).unwrap();
        let (radial,c)=a.monomial(&[0,128,0]).unwrap();
        assert_eq!(radial,vec![64,64]);
        let expected=(0..64).fold(Atom::one(),|a,j|a*Atom::num(2*j+1)/(parse!("11/2")+Atom::num(2*j)));
        equal(c,expected); assert!(a.retained_entries()<200);
    }
    #[test]
    fn budget_failure_is_explicit() {
        let mut a=AngularAverager::new(3,parse!("d"),1).unwrap();
        assert!(a.monomial(&[0,2,2,0,2,0]).is_err());
    }
    #[test]
    fn degree_preflight_precedes_recursive_allocation() {
        let mut a=AngularAverager::new(2,parse!("d"),200000).unwrap();
        assert!(a.monomial(&[0,65534,0]).unwrap_err().contains("preflight"));
        assert_eq!(a.retained_entries(),0);
        let mut b=AngularAverager::new(2,parse!("d"),200000).unwrap().with_max_total_degree(16).unwrap();
        assert!(b.monomial(&[0,10,0]).is_err());
    }
}
