//! Cache infrastructure benchmark on 33,000 exact analytic constant seeds.
//! Run: cargo run --release --example cache_bank_benchmark -- 33000 target/cache-bank-benchmark
//! This measures cache selection, not integral evaluation or Monte Carlo throughput.
use std::{
    collections::BTreeMap,
    sync::atomic::{AtomicUsize, Ordering},
    time::Instant,
};
use symbolica::prelude::*;
use symbolica_amflow::{kinematics::KinematicSystem, transport_cache::*, *};

fn memory() -> serde_json::Value {
    let status = std::fs::read_to_string("/proc/self/status").unwrap();
    let get = |key: &str| {
        status
            .lines()
            .find_map(|line| line.strip_prefix(key))
            .unwrap()
            .split_whitespace()
            .next()
            .unwrap()
            .parse::<u64>()
            .unwrap()
    };
    serde_json::json!({"resident_kib":get("VmRSS:"),"high_water_kib":get("VmHWM:")})
}
fn emit(value: serde_json::Value) {
    println!("{value}");
}
struct Guarded<'a, F> {
    identity: &'a BoundaryIdentity,
    distance: &'a ScaledDistance<F>,
    checks: AtomicUsize,
    lazy: bool,
}
impl<F: Fn(&CachedBoundary, &CachedPoint) -> Result<bool>> TransportCost for Guarded<'_, F> {
    fn cost(
        &self,
        source: &CachedBoundary,
        target: &CachedPoint,
        p: Precision,
    ) -> Result<Option<Float>> {
        let count = self.checks.fetch_add(1, Ordering::Relaxed) + 1;
        if count.is_multiple_of(5000) {
            eprintln!("checked {count} guarded candidates");
        }
        if !self
            .identity
            .conditions_admit_straight_path(&source.point, target, p, 20)?
        {
            return Ok(None);
        }
        self.distance.cost(source, target, p)
    }
    fn lower_bound(
        &self,
        source: &CachedBoundary,
        target: &CachedPoint,
        p: Precision,
    ) -> Result<Option<Float>> {
        if self.lazy {
            self.distance.lower_bound(source, target, p)
        } else {
            Ok(None)
        }
    }
    fn compare_tied_costs(
        &self,
        left: &CachedBoundary,
        right: &CachedBoundary,
        target: &CachedPoint,
        p: Precision,
    ) -> Result<std::cmp::Ordering> {
        self.distance.compare_tied_costs(left, right, target, p)
    }
}
fn main() -> Result<()> {
    let n = std::env::args()
        .nth(1)
        .map(|s| s.parse::<usize>().unwrap())
        .unwrap_or(33000);
    assert!(n > 0 && n <= 33000);
    let (s, t, epsilon) = symbol!(
        "cache_benchmark::s",
        "cache_benchmark::t",
        "cache_benchmark::epsilon"
    );
    let flow = RustFlow::with_conditions(
        KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([
                (s, vec![vec![Atom::new()]]),
                (t, vec![vec![Atom::new()]]),
            ]),
        },
        &[parse!("cache_benchmark::one")],
        &Atom::one(),
        Prescription::PlusI0,
        "s>0,t>0; exact constant solution",
        &[Atom::var(s) * Atom::var(t)],
    )?;
    let p = Precision::decimal(90)?;
    let query_p = Precision::decimal(60)?;
    let range = EpsilonRange::new(0, 2)?;
    emit(
        serde_json::json!({"phase":"start","points":n,"coordinates":2,"components":1,"epsilon_coefficients":3,"working_bits":p.bits,"query_bits":query_p.bits,"requested_digits":20,"memory":memory(),"scope":"synthetic exact analytic Y(s,t,epsilon)=1, not an integral benchmark; no f64 coordinates or values"}),
    );
    let start = Instant::now();
    let mut cache = RustFlowCache::default();
    for index in 0..n {
        let coordinates = BTreeMap::from([
            (
                s,
                Atom::num(Rational::from((((index % 220) + 1) as i64, 10))),
            ),
            (
                t,
                Atom::num(Rational::from((((index / 220) + 1) as i64, 10))),
            ),
        ]);
        cache.insert(CachedBoundary {identity:flow.identity().clone(),point:CachedPoint::Exact(coordinates),kind:PointKind::Physical,range,
            coefficients:vec![vec![p.i(1)],vec![p.zero()],vec![p.zero()]],accuracy:BoundaryAccuracy::supplied(80,p.bits,vec![vec![p.real(0)];3],"Analytic solution identically one: epsilon^0=1 and every other coefficient=0 exactly. Evidence comes from the formula, not mantissa length.")?})?;
    }
    emit(
        serde_json::json!({"phase":"cold_indexed_insert","elapsed_ns":start.elapsed().as_nanos(),"entries":cache.len(),"memory":memory()}),
    );
    let directory = std::env::args()
        .nth(2)
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::path::PathBuf::from("target/cache-bank-benchmark"));
    std::fs::create_dir_all(&directory)?;
    let start = Instant::now();
    cache.save(&directory)?;
    emit(
        serde_json::json!({"phase":"binary_save","elapsed_ns":start.elapsed().as_nanos(),"bytes":std::fs::metadata(directory.join("physical-boundaries.bin"))?.len(),"memory":memory()}),
    );
    drop(cache);
    let start = Instant::now();
    let cache = RustFlowCache::load(&directory)?;
    assert_eq!(cache.len(), n);
    emit(
        serde_json::json!({"phase":"binary_load","elapsed_ns":start.elapsed().as_nanos(),"entries":cache.len(),"memory":memory()}),
    );
    let distance = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    let last = cache
        .entries()
        .last()
        .unwrap()
        .point
        .rounded_coordinates_as_exact()?;
    for exact in [false, true] {
        let target = CachedPoint::Exact(
            last.iter()
                .map(|(&key, value)| {
                    (
                        key,
                        if exact {
                            value.clone()
                        } else {
                            value + Atom::num((1, 100))
                        },
                    )
                })
                .collect(),
        );
        let query = BoundaryQuery::new(flow.identity(), &target, range, 20)?;
        let guarded = Guarded {
            identity: flow.identity(),
            distance: &distance,
            checks: AtomicUsize::new(0),
            lazy: false,
        };
        let start = Instant::now();
        let chosen = cache.best(&query, &guarded, query_p)?.unwrap();
        let duration = start.elapsed().as_nanos();
        assert_eq!(chosen.boundary.point.rounded_coordinates_as_exact()?, last);
        let expected_cost = chosen.cost.clone();
        emit(
            serde_json::json!({"phase":if exact {"baseline_exact_hit_selection"} else {"baseline_nearest_selection"},"elapsed_ns":duration,"guard_checks":guarded.checks.load(Ordering::Relaxed),"cost":expected_cost.to_string(),"memory":memory()}),
        );
        let lazy = Guarded {
            identity: flow.identity(),
            distance: &distance,
            checks: AtomicUsize::new(0),
            lazy: true,
        };
        let start = Instant::now();
        let chosen = cache.best(&query, &lazy, query_p)?.unwrap();
        let duration = start.elapsed().as_nanos();
        assert_eq!(chosen.cost, expected_cost);
        assert_eq!(chosen.boundary.point.rounded_coordinates_as_exact()?, last);
        emit(
            serde_json::json!({"phase":if exact {"lazy_exact_hit_selection"} else {"lazy_nearest_selection"},"elapsed_ns":duration,"guard_checks":lazy.checks.load(Ordering::Relaxed),"compatible_candidates":n,"cost":chosen.cost.to_string(),"memory":memory()}),
        );
    }
    emit(serde_json::json!({"phase":"complete","memory":memory()}));
    Ok(())
}
