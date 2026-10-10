#!/usr/bin/env python3
"""Replay saved exact checks in scratch storage and bind an independent review."""
import contextlib
import gzip
import hashlib
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SINGLE = HERE.parent/'2026-10-10-singleton-energy-ward-assessment'
K = HERE.parent/'2026-10-10-requested-ray-point-integration'
read_bytes = Path.read_bytes
read_text = Path.read_text
sha = lambda p: hashlib.sha256(read_bytes(p)).hexdigest()
read = lambda p: json.loads(read_text(p))


def main():
    assert not (HERE/'artifact-manifest.json').exists()
    assert not (SINGLE/'artifact-manifest.json').exists()
    archive_path = K/'requested-ray-point-native-archives.json'
    archive = read(archive_path)
    mapping = {str(K/row['original_path']): row for row in archive['entries']}
    def bytes_with_archive(path):
        if path.exists(): return read_bytes(path)
        row = mapping[str(path)]
        encoded = read_bytes(K/row['archive_path'])
        assert hashlib.sha256(encoded).hexdigest() == row['archive_sha256']
        data = gzip.decompress(encoded)
        assert hashlib.sha256(data).hexdigest() == row['original_sha256']
        return data
    def text_with_archive(path, *args, **kwargs):
        if path.exists(): return read_text(path,*args,**kwargs)
        assert not args and not kwargs
        return bytes_with_archive(path).decode()
    def load(name,path):
        spec=importlib.util.spec_from_file_location(name,path)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        return module
    saved_inputs = {str(p.relative_to(ROOT)):sha(p) for directory in (HERE,SINGLE)
        for p in sorted(directory.glob('*')) if p.is_file() and p.name != Path(__file__).name}
    runs=[]
    sys.path.insert(0,str(HERE))
    for i,(folder,script,result) in enumerate([
        (HERE,'check.py','exact-checks.json'),
        (HERE,'redundancy_check.py','redundancy-checks.json'),
        (SINGLE,'check.py','checks.json')]):
        module=load('ward_review_'+str(i),folder/script)
        with tempfile.TemporaryDirectory(prefix='ward-check-review-') as scratch:
            module.HERE=Path(scratch)
            try:
                Path.read_bytes=bytes_with_archive
                Path.read_text=text_with_archive
                module.main()
            finally:
                Path.read_bytes=read_bytes
                Path.read_text=read_text
            reproduced=read(Path(scratch)/result)
            expected=read(folder/result)
            assert reproduced == expected, 'saved exact evidence differs'
        # Verify every reported input against raw or restored bytes.
        for path,expected in read(folder/result).get('inputs_sha256',{}).items():
            p=Path(path)
            if not p.is_absolute():p=ROOT/p
            assert hashlib.sha256(bytes_with_archive(p)).hexdigest()==expected
        runs.append({'script':str((folder/script).relative_to(ROOT)),
            'script_sha256':sha(folder/script),'result':str((folder/result).relative_to(ROOT)),
            'result_sha256':sha(folder/result),'exact_saved_json_reproduced':True})
    first=read(HERE/'exact-checks.json');second=read(HERE/'redundancy-checks.json');third=read(SINGLE/'checks.json')
    assert first['checks']==216
    assert second['generator_row_checks']==second['raised_cut_moment_checks']==64
    assert third['checks']==210
    for path,expected in saved_inputs.items():assert sha(ROOT/path)==expected
    common={'status':'passed','kind':'read-only mathematical/source-presentation audit plus exact saved-check replay',
        'production_modified':False,'new_native_replay':False,'new_numerical_predictions':0,
        'saved_checks_reproduced':runs,'input_files_sha256':saved_inputs,
        'archive_map_sha256':sha(archive_path),'audit_script_sha256':sha(Path(__file__)),
        'findings':[
            'At a fully factorized selected leg, local dilation gives D-2n; the exact multiplication identity is q^2 C_(n+1)=C_n. The earlier tentative C1-only new-source suggestion is explicitly superseded.',
            'The singleton C1 normalized global boost annihilates all Gram factors. Pure cE completions of total polynomial degree k add +k, giving D-2-sum(a_energy). Nonzero rational scaling cancels from the identity.',
            'Use the existing regulated joint high-D origin prescription before meromorphic continuation. These are not pointwise inverse-energy products, reciprocal source admission, or an unchanged Cn>1 theorem.',
            'For global dilation, the displayed minus 2 eta derivative and minus sum(mu derivative) signs agree with D_j=q_j^2-eta and normalized upper distributions. Exact homogeneous completion degree and all physical shift/mass assumptions are necessary.',
            'Rebasing mixed completions can change the finite-eta family. A complete target map and original-mass coefficient derivatives must be bound before any new numerical comparison.',
            'Moment/generator arithmetic does not establish native closure benefit. The redundancy audit reads native source descriptors and re-evaluates the documented generator, not native coefficient tables or generated proof rules.'
        ],'remaining_requirements':[
            'Any singleton-energy extension requires a separately reviewed source experiment and full native replay under exact pure-energy/base guards.',
            'General completion rebasing remains a design, not tested production admission.',
            'No three-loop or four-loop numerical acceptance follows from these checks.'
        ]}
    for folder in (HERE,SINGLE):
        report=dict(common)
        report['report_scope']='factorized compact-leg, existing dilation and rebase assessment' if folder==HERE else 'pure occupied-energy polynomial extension of SingletonGerm C1 Ward'
        (folder/'independent-review.json').write_text(json.dumps(report,indent=2)+'\n')
        files={str(p.relative_to(folder)):{'sha256':sha(p),'bytes':p.stat().st_size}
            for p in sorted(folder.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}
        manifest={'status':'frozen_assessment_only','production_admission_changed':False,
            'native_closure_claim':False,'files':files,'file_count':len(files)}
        (folder/'artifact-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        print(json.dumps({'report':str(folder.relative_to(ROOT)),'files':len(files),
            'manifest_sha256':sha(folder/'artifact-manifest.json')}))

if __name__=='__main__':main()
