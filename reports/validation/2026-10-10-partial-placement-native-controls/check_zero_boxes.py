#!/usr/bin/env python3
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
checks = []
for mode in ['single-all', 'single-slot0', 'double-all', 'double-slots12', 'double-slots24']:
    r = json.loads((HERE / mode / 'result.json').read_text())
    physical_arity = r['physical_arity']
    cuts = r['cuts']
    def contains(box, point):
        return all((lo is None or n >= lo) and (hi is None or n <= hi)
                   for n, (lo, hi) in zip(point, box))
    for box in r['zero_domains']:
        assert box[:5] == [[1, None]] * 5
        assert box[5:9] == [[None, 0]] * 4
        assert box[physical_arity:] == [[0, 0]] * (16 - physical_arity)
        point = [1]*5 + [0]*11
        for i, (lo, hi) in enumerate(box):
            if lo is not None:
                point[i] = max(point[i], lo)
        assert contains(box, point)
        n = 1
        for slot in [i for i in range(5) if i not in cuts]:
            deleted = point.copy(); deleted[slot] = 0
            assert not any(contains(b, deleted) for b in r['zero_domains'])
            n += 1
        for slot in range(5, 9):
            invalid = point.copy(); invalid[slot] = 1
            assert not any(contains(b, invalid) for b in r['zero_domains'])
            n += 1
        for slot in range(9, physical_arity):
            invalid = point.copy(); invalid[slot] = -1
            assert not any(contains(b, invalid) for b in r['zero_domains'])
            n += 1
        for slot in range(physical_arity, 16):
            invalid = point.copy(); invalid[slot] = 1
            assert not any(contains(b, invalid) for b in r['zero_domains'])
            n += 1
        checks.append({'mode': mode, 'box': box, 'passed': n})
out = HERE / 'zero-box-checks.json'
assert not out.exists()
report = {'scope': 'Exact saved-domain audit only, including exclusion of deleted uncut supports; native application replay is separately recorded in result.json.',
          'checks': checks, 'total_passed': sum(x['passed'] for x in checks),
          'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
out.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({'zero_domain_checks': report['total_passed'], 'result': 'PASS'}))
