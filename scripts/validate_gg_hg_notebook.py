"""Execute the unchanged notebook from a copied, committed cold-cache archive.

Local acceptance helper. It never reads numerical values from reference seeds,
changes the live installation, or writes inside the live run directory.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
from time import perf_counter_ns


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--notebook', type=Path, required=True)
    parser.add_argument('--run-directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--profile', action='store_true')
    args = parser.parse_args()
    notebook, run, output = (p.resolve() for p in (args.notebook, args.run_directory, args.output))
    if output == run or run in output.parents or output in run.parents:
        raise ValueError('Output must be separate from the live run directory.')
    if output.exists():
        raise FileExistsError('Use a new output directory.')
    if (run / 'force-preparation-pending.json').exists():
        raise RuntimeError('Cold numerical archive is not yet committed.')
    generation_bytes = (run / 'numerical-generation.json').read_bytes()
    generation = json.loads(generation_bytes)
    archive = (run / generation['archive_directory']).resolve()
    if archive.parent != run / 'numerical-history' or archive.name != generation['generation']:
        raise RuntimeError('Unexpected numerical archive identity.')
    acceptance_path = run / 'acceptance.json'
    acceptance_bytes = acceptance_path.read_bytes()
    acceptance = json.loads(acceptance_bytes)
    required = {'boundaries', 'transport', 'amplitude'}
    completed = {s['stage'] for s in acceptance['stages'] if s.get('mode') == 'cold'}
    if not required.issubset(completed):
        raise RuntimeError('Cold native physics stages have not all completed.')
    if acceptance['reference_comparison']['coefficients'] != 4360:
        raise RuntimeError('Expected full cold coefficient comparison.')
    import symbolica.core as core
    attestation = acceptance['runtime_attestation']
    if digest(core.__file__) != attestation['loaded_extension']['sha256']:
        raise RuntimeError('Loaded native runtime does not match the archived calculation.')
    for key, source in [('notebook', notebook), ('controller', notebook.with_name('gg_hg_support.py'))]:
        if digest(source) != attestation['execution_sources'][key]['sha256']:
            raise RuntimeError(f'{key} does not match the cold calculation.')
    original = {}
    for name in ('seeds', 'transport'):
        path = archive / name / 'physical-boundaries.bin'
        original[name] = {'path': str(path), 'sha256': digest(path), 'bytes': path.stat().st_size}
    output.mkdir(parents=True)
    bank = output / 'cache'
    bank.mkdir()
    for name in original:
        shutil.copytree(archive / name, bank / name)
        if digest(bank / name / 'physical-boundaries.bin') != original[name]['sha256']:
            raise RuntimeError('Copied cache does not match source.')
    if (run / 'force-preparation-pending.json').exists() or (run / 'numerical-generation.json').read_bytes() != generation_bytes:
        raise RuntimeError('Generation changed while copying; retry with a new directory.')
    support = module(notebook.with_name('gg_hg_support.py'), 'gg_hg_populated_support')
    report = {
        'schema': 'rustflow-populated-notebook-copy-v1',
        'recorded_utc': datetime.now(timezone.utc).isoformat(),
        'scope': 'Actual original notebook execution with a prepared native session from a copy of its completed cold bank; no browser pixel or interactive-click claim.',
        'status': 'running', 'archive': str(archive), 'source_banks': original,
        'acceptance_observation_sha256': hashlib.sha256(acceptance_bytes).hexdigest(),
        'runtime_attestation': attestation, 'profile': {},
    }
    (output / 'acceptance-observation.json').write_bytes(acceptance_bytes)
    (output / 'generation-observation.json').write_bytes(generation_bytes)

    def save_report():
        (output / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')

    def timed(label, fn):
        started = perf_counter_ns()
        result = fn()
        report['profile'][label] = {'elapsed_ns': perf_counter_ns() - started}
        save_report()
        print(label, report['profile'][label], flush=True)
        return result

    session = None
    try:
        session = timed('construct_and_load', lambda: support.CalculationSession(
            notebook.parent / 'data/gg_hg/native-model.json', bank,
            seed_digits=30, digits=20, workers=1, boundary_workers=1))
        before = support.boundary_evidence(session.cache)
        if len(before) != acceptance['binary_restart']['entries']:
            raise RuntimeError('Archive has unexpected cache entry count.')
        # A direct API pass separates cache lookup/evidence work from persistence.
        if args.profile:
            direct = {}
            def evaluate_only():
                for topology, configuration in session.configurations:
                    system = session.systems[topology]
                    result = system.evaluate(session.cache, session._destination(system, configuration),
                                             configuration.root_sheets, options=session._options())
                    if not result.cache_hit or result.steps or result.inserted_points:
                        raise RuntimeError('Profile was not an unchanged exact hit.')
                    direct[configuration.label] = result
                return direct
            timed('sixteen_exact_evaluations_without_save', evaluate_only)
            if support.boundary_evidence(session.cache) != before:
                raise RuntimeError('Exact evaluations changed cached values/evidence.')
            timed('sixteen_unchanged_saves', lambda: [session.cache.save(output / 'save-profile') for _ in session.configurations])
        results = timed('original_controller_warm_transport', lambda: session.submit('transport').result())
        if len(results) != 16 or any(not r.cache_hit or r.steps or r.inserted_points for r in results.values()):
            raise RuntimeError('Notebook preparation did not yield sixteen exact hits.')
        if support.boundary_evidence(session.cache) != before:
            raise RuntimeError('Controller exact hits changed cached values/evidence.')
        if args.profile:
            for key, result in results.items():
                other = direct[key]
                if (result.coefficients, result.comparison_errors, result.provenance) != (other.coefficients, other.comparison_errors, other.provenance):
                    raise RuntimeError('Direct and controller results differ.')
        timed('original_controller_amplitude', lambda: session.submit('amplitude').result())
        app = module(notebook, 'gg_hg_populated_notebook').app
        outputs, definitions = timed('original_notebook_execution', lambda: app.run(defs={'session': session}))
        if definitions['session'] is not session or definitions['state']['done'] is not True:
            raise RuntimeError('Notebook did not use the completed native session.')
        fragments = []
        for index, value in enumerate(outputs):
            if value is None:
                continue
            if not hasattr(value, '_repr_html_'):
                raise TypeError(f'Unexpected notebook output {index}: {type(value)}')
            html = value._repr_html_()
            fragments.append(html)
            (output / f'cell-{index:02}.html').write_text(html)
        combined = '\n'.join(fragments)
        for required_text in ['Native HEPKit diagrams', 'Native W/Z form factors', 'Coherent observables', 'Transport and accumulated boundary reuse']:
            if required_text not in combined:
                raise RuntimeError(f'Missing populated output: {required_text}')
        diagrams = session.amplitude.diagrams
        for index, diagram in enumerate(diagrams):
            svg = diagram._repr_svg_()
            if '<svg' not in svg:
                raise RuntimeError('Native diagram did not render SVG.')
            (output / f'diagram-{index:02}.svg').write_text(svg)
        if not diagrams:
            raise RuntimeError('Native amplitude has no diagrams.')
        report['notebook'] = {
            'output_count': len(fragments), 'native_diagram_count': len(diagrams),
            'transport_configurations': len(results), 'form_factor_count': sum(len(f.values) for f in session.form_factors.values()),
            'observable_digits': session.observables.verified_relative_digits,
            'rendered_populated_tables': True,
            'serialized_output_bytes': len(combined.encode()),
        }
        report['status'] = 'passed'
    except BaseException as error:
        report['status'] = 'failed'
        report['error'] = f'{type(error).__name__}: {error}'
        raise
    finally:
        if session is not None:
            session._pool.shutdown(wait=True)
        report['source_banks_unchanged'] = all(digest(v['path']) == v['sha256'] for v in original.values())
        save_report()
        if not report['source_banks_unchanged']:
            raise RuntimeError('Immutable source cache changed during validation.')


if __name__ == '__main__':
    main()
