"""Extract exact dyadic scientific input from the frozen native acceptance bank."""
from pathlib import Path
import gzip
import hashlib
import json
import symbolica
from symbolica.community.hep.integration import BoundaryCache, HiggsJetIntegralSystem

ROOT = Path('/common/dev/amflow')
SOURCE = ROOT / 'target/gg-hg-notebook-acceptance/publication-native30/seeds/physical-boundaries.bin'
OUTPUT = ROOT / 'target/gg-hg-browser-seeds'
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
original_sha = sha(SOURCE)
systems = {name: HiggsJetIntegralSystem(name) for name in ('planar', 'nonplanar')}
cache = BoundaryCache.load(SOURCE.parent)

def number(value):
    numerator, denominator = value.as_integer_ratio()
    return [str(numerator), str(denominator), value.precision]

records = []
seen = set()
for topology, system in systems.items():
    for config in system.configurations():
        matches = [entry for entry in cache.entries()
                   if entry.coordinates == config.start
                   and entry.root_sheets == config.root_sheets
                   and len(entry.coefficients[0]) == system.dimension]
        assert len(matches) == 1, config.label
        entry = matches[0]
        seen.add(config.label)
        records.append({
            'label': config.label, 'topology': topology,
            'mass': config.mass, 'permutation': config.permutation,
            'dimension': system.dimension,
            'coordinates': [str(config.start[s]) for s in system.coordinates],
            'root_sheets': {str(s): sign for s, sign in config.root_sheets.items()},
            'leading_power': entry.leading_power,
            'verified_digits': entry.verified_digits,
            'input_verified_digits': entry.input_verified_digits,
            'working_bits': entry.working_bits,
            'provenance': entry.provenance,
            'origin_identity': entry.identity,
            'coefficients': [[[number(z.real), number(z.imag)] for z in row]
                             for row in entry.coefficients],
            'comparison_errors': [[number(x) for x in row] for row in entry.comparison_errors],
        })
assert len(seen) == len(cache) == 16
document = {
    'schema': 'higgs-jet-supplied-boundaries-v1',
    'encoding': 'exact-integer-ratio-with-binary-precision',
    'origin': {
        'description': 'Rust-generated physical starting boundaries; independently refined native acceptance run. No plugin numerical seeds.',
        'report': 'https://github.com/alphal00p/RustFlow/blob/115d9e5/reports/validation/2026-10-06-gg-hg-publication-complete/summary.json',
        'binary_sha256': original_sha,
        'native_extension_sha256': sha(Path(symbolica.core.__file__)),
    },
    'boundaries': records,
}
raw = (json.dumps(document, separators=(',', ':'), sort_keys=True) + '\n').encode()
packed = gzip.compress(raw, compresslevel=9, mtime=0)
(OUTPUT / 'boundaries.json.gz').write_bytes(packed)
assert sha(SOURCE) == original_sha
summary = {
    'binary_bytes': SOURCE.stat().st_size,
    'json_bytes': len(raw), 'gzip_bytes': len(packed),
    'sha256': hashlib.sha256(packed).hexdigest(),
    'configurations': len(records),
    'coefficients': sum(sum(map(len, r['coefficients'])) for r in records),
    'verified_digits': sorted({r['verified_digits'] for r in records}),
    'input_verified_digits': sorted({r['input_verified_digits'] for r in records}),
    'working_bits': sorted({r['working_bits'] for r in records}),
    'source_unchanged': True,
}
(OUTPUT / 'export-summary.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps(summary, indent=2), flush=True)
