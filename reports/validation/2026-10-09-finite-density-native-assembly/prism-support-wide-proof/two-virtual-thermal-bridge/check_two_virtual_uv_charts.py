"""Exact audit of rank-two UV charts against already archived native evidence."""
from fractions import Fraction as Q
from pathlib import Path
import gzip
import hashlib
import json

root = Path('/common/dev/rustflow_fermi/reports/validation/2026-10-09-finite-density-native-assembly/prism-support-wide-proof/native-full-kinematic')
runs = []
for a, b in [(1, 5), (1, 7), (5, 7)]:
    path = root / f'prism-cuts-{a}-{b}-supports.json.gz'
    raw = path.read_bytes()
    data = json.loads(gzip.decompress(raw))
    counts = dict(supports=0, full_rank=0, rank_deficient=0, primary_charts=0,
                  transverse_subcharts=0, uv_operation_bound=0, native_operations=0)
    for proof in data['supports']:
        counts['supports'] += 1
        counts['native_operations'] += proof['resolution_operations']
        full = proof['full_kinematics']
        if full is None:
            counts['rank_deficient'] += 1
            continue
        counts['full_rank'] += 1
        slots = full['parameter_slots']
        n = len(slots)
        rows = [[Q(x) for x in proof['routed_rows'][s][2:]] for s in slots]
        u = {tuple(term['exponents']): Q(term['coefficient']) for term in full['u']}
        assert all(sum(e) == 2 and c > 0 for e, c in u.items())
        all_terms = list(full['u'])
        all_terms += [t for p in full['pairs'] for t in p['terms']]
        all_terms += [t for p in full['diagonal_masses'] + full['virtual_masses'] for t in p]
        counts['uv_operation_bound'] += 16*n*n*len(all_terms)
        for primary, r in enumerate(rows):
            transverse = []
            for pivot, s in enumerate(rows):
                determinant = r[0]*s[1]-r[1]*s[0]
                exponent = tuple(int(j == primary)+int(j == pivot) for j in range(n))
                assert u.get(exponent, Q(0)) == determinant**2
                if determinant:
                    transverse.append(pivot)
            assert transverse
            assert all(sum(term['exponents'][j] for j in transverse) >= 1 for term in all_terms)
            counts['primary_charts'] += 1
            counts['transverse_subcharts'] += len(transverse)
    reservation = 16*8*8*20000
    phase = (200000000-reservation)//3
    assert counts['native_operations'] <= phase
    assert counts['uv_operation_bound'] <= reservation
    runs.append(dict(file=str(path), sha256=hashlib.sha256(raw).hexdigest(),
                     construction_phase_budget=phase, uv_operation_reservation=reservation,
                     **counts))
result = dict(status='pass', scope='exact UV chart audit only; no new physical admission or period',
              theorem='generic two-virtual-loop maximum-parameter UV radial charts', runs=runs)
Path('/tmp/two-virtual-uv-chart-audit.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
