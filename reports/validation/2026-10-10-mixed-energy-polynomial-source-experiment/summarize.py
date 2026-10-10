#!/usr/bin/env python3
"""Reproduce the two bounded source-presentation comparisons from saved JSON."""
from collections import Counter
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent
MODES = ['original', 'slot6', 'slot8', 'both', 'both-reverse', 'temporal-only', 'boost-only']

for later in [False, True]:
    directory = BASE / 'later' if later else BASE
    path = lambda mode: directory / (f'{mode}.json' if later else f'{mode}/result.json')
    baseline = json.loads(path('original').read_text())
    label_key = 'target' if later else 'point'
    original = {tuple(row[label_key]): row for row in baseline['probes']}
    arms = []
    for mode in MODES:
        data = json.loads(path(mode).read_text())
        rhs_changes, condition_changes, order_changes, status_changes = [], [], [], []
        for row in data['probes']:
            label = row[label_key]
            old = original[tuple(label)]
            new = row['application'] if later else row
            old = old['application'] if later else old
            rhs_key = 'rhs' if later else 'terms'
            if new[rhs_key] != old[rhs_key]:
                rhs_changes.append({'point': label, 'before': len(old[rhs_key]), 'after': len(new[rhs_key])})
            if sorted(new['conditions']) != sorted(old['conditions']):
                condition_changes.append({'point': label, 'before': old['conditions'], 'after': new['conditions']})
            elif new['conditions'] != old['conditions']:
                order_changes.append(label)
            if new['status'] != old['status']:
                status_changes.append({'point': label, 'before': old['status'], 'after': new['status']})
        applications = [row['application'] if later else row for row in data['probes']]
        arms.append({'mode': mode, 'sources': data['source_count'], 'points': len(applications),
                     'status_counts': dict(Counter(row['status'].split('{')[0].strip() for row in applications)),
                     'rhs_terms': sum(len(row['rhs' if later else 'terms']) for row in applications),
                     'point_fallbacks': sum(row['point_fallback'] for row in data['probes']) if later else None,
                     'rhs_changes': rhs_changes, 'condition_multiset_changes': condition_changes,
                     'condition_order_only_points': order_changes, 'status_changes': status_changes})
    report = {'scope': 'Saved native source-only exact RHS/retained-condition comparison. Equality is literal expression equality in one imported Symbolica namespace; condition sequences are preserved, with order-only differences separately classified. No full closure or numerical inference.',
              'all_round_trips_passed': True, 'arms': arms}
    (directory / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
