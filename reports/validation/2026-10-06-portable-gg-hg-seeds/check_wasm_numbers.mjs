import { readFile, writeFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import { basename, join } from 'node:path';
import { pathToFileURL } from 'node:url';
const [runtime, wheel, bundle, loader, reportPath] = process.argv.slice(2);
const { loadPyodide } = await import(pathToFileURL(join(runtime, 'pyodide.mjs')));
const py = await loadPyodide({indexURL: runtime, env: {
  SYMBOLICA_LICENSE: process.env.SYMBOLICA_LICENSE || '',
}});
await py.loadPackage('micropip');
const wheelBytes = await readFile(wheel);
py.FS.writeFile('/' + basename(wheel), wheelBytes);
py.FS.writeFile('/boundaries.json.gz', await readFile(bundle));
py.FS.writeFile('/loader.py', await readFile(loader));
py.globals.set('wheel_uri', 'emfs:/' + basename(wheel));
const result = await py.runPythonAsync(`
import micropip
await micropip.install(wheel_uri)
import json, gzip, runpy
from time import perf_counter_ns
from symbolica import ComplexFloat
number = runpy.run_path('/loader.py')['_number']
started = perf_counter_ns()
document = json.loads(gzip.decompress(open('/boundaries.json.gz','rb').read()))
decoded = perf_counter_ns()
values = errors = 0
for record in document['boundaries']:
    for row in record['coefficients']:
        for z in row:
            value = ComplexFloat(number(z[0]), number(z[1]))
            assert value.real.as_integer_ratio() == tuple(map(int, z[0][:2]))
            assert value.imag.as_integer_ratio() == tuple(map(int, z[1][:2]))
            assert value.real.precision == value.imag.precision == record['working_bits']
            values += 1
    for row in record['comparison_errors']:
        for x in row:
            value = number(x)
            assert value.as_integer_ratio() == tuple(map(int, x[:2]))
            assert value.precision == x[2]
            errors += 1
ended = perf_counter_ns()
json.dumps(dict(status='passed', coefficients=values, absolute_errors=errors,
    decode_ns=decoded-started, number_import_and_exact_verification_ns=ended-decoded,
    total_ns=ended-started, working_bits=415,
    scope='Actual single-core Pyodide arbitrary-precision data import; no RustFlow transport or amplitude calculation.'))
`);
const report = JSON.parse(result);
report.wheel_sha256 = createHash('sha256').update(wheelBytes).digest('hex');
report.wheel_path = wheel;
report.runtime_path = runtime;
report.node_version = process.version;
await writeFile(reportPath, JSON.stringify(report, null, 2)+'\n');
console.log(JSON.stringify(report, null, 2));
