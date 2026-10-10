#!/usr/bin/env python3
"""Bind one successful selected-sector pipeline and archive completed diagnostics."""
import gzip
import hashlib
import json
import re
import subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
K=HERE.parent/'2026-10-10-requested-ray-point-integration'
read=lambda p:json.loads(p.read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    assert not (HERE/'artifact-manifest.json').exists()
    binding=read(HERE/'run-binding.json');resource=read(HERE/'resources.json')
    folder=HERE/'selected-cut';prediction=folder/'prediction-18-60-8.json';data=read(prediction)
    assert binding['exit_code']==resource['exit_code']==0 and binding['post_run_source_and_executable_unchanged']
    assert not (folder/'failure.json').exists()
    assert data['cut_slots']==[0] and data['full_amplitude'] is False
    assert (data['digits'],data['series_order'],data['occupied_start_scale'],data['guard_digits'],data['epsilon'])==(18,60,8,40,'-5/4')
    assert data['independent_reference_comparisons']==0
    assert data['report']['values']==['(0+0i)','(0+0i)']
    assert data['report']['construction']=='weighted_amf' and data['report']['basis_size']==7
    assert data['report']['source_options']=={'free_virtual_zero_sectors':True,'policy':'polynomial-closure-v1','positive_compact_energy_powers':False}
    assert 'integrated_coefficients: 7, integrated_products: 7' in data['report']['boundary']
    assert read(folder/'input.json')==read(ROOT/binding['settings']['RUSTFLOW_WEIGHTED_INPUT'])
    log=HERE/'resources.log';text=log.read_text()
    times=re.findall(r'test result: ok\..*?finished in ([0-9.]+)s',text);assert len(times)==1
    stages=['compiling the auxiliary-mass differential equation (7 basis integrals, 3 loops)',
        'generating and matching native asymptotic boundary regions (7 basis integrals, 3 loops)',
        'transporting the auxiliary-mass solution (7 basis integrals, 3 loops)',
        'extracting the dimensionally regulated physical limit (7 basis integrals, 3 loops)']
    assert all(s in text for s in stages)
    corpus=folder/'native-closure';closed=corpus/'round-015-closed.json';row=read(closed)
    assert row['status']=='closed' and len(row['frontier'])==7 and row['round']==15
    toolpath=K/'native-digest-build.json';tool=read(toolpath);exe=Path(tool['executable']['path'])
    assert sha(exe)==tool['executable']['sha256']
    def verify_proof(path,meta):
        digest,n=subprocess.check_output([str(exe),str(path)],text=True).split()
        assert digest==meta['program_blake3'] and int(n)==meta['program_bytes']==path.stat().st_size
    proof=closed.with_suffix('.bin');verify_proof(proof,row)
    link=row['active_state']['requested_discovery'];assert re.fullmatch(r'requested-discovery-[0-9]+[.]json',link['path'])
    transaction=corpus/link['path'];tx=read(transaction)
    assert subprocess.check_output([str(exe),str(transaction)],text=True).split()[0]==link['blake3']
    provisional=corpus/'round-015-provisional.json';pr=read(provisional);pp=provisional.with_suffix('.bin');verify_proof(pp,pr)
    assert pr['active_state']==row['active_state'] and tx['status']=='complete'
    assert tx['native_program']['blake3']==pr['program_blake3'] and tx['native_program']['bytes']==pr['program_bytes']
    assert tx['source_measure_id']==row['measure_id'] and tx['completed_search_implies_coverage'] is False
    bound=[HERE/'run.py',HERE/'run-binding.json',HERE/'resources.json',HERE/'resources.provenance.json',log,
        prediction,folder/'input.json',folder/'configuration.json',folder/'sector.json',closed,proof,provisional,pp,transaction,toolpath,
        K/'requested-ray-point-source-hashes.json',K/'requested-ray-point-build-provenance.json']
    summary={'status':'selected_sector_pipeline_passed','test_count':1,'cut_slots':[0],
        'full_amplitude':False,'profile':{'digits':18,'series_order':60,'start_scale':8,'epsilon':'-5/4','guard_digits':40},
        'native_closed_basis':7,'native_closure_rounds':16,'native_rule_count':row['native_rule_count'],
        'boundary_integrated_coefficients':7,'boundary_integrated_products':7,'verified_log_stages':stages,
        'construction':'weighted_amf','observed_values':data['report']['values'],
        'harness_seconds':float(times[0]),'wall_seconds':resource['wall_seconds'],'peak_child_rss_kib':resource['peak_child_rss_kib'],
        'independent_reference_comparisons':0,'profile_refinements':0,'supplied_oracle_comparisons':0,
        'scope':'Fresh selected original cut-0 pipeline at one profile, including recursive ordinary hard children. Exact-zero outputs alone do not validate nonzero hard-boundary accuracy, a full three-loop amplitude, precision stability, Laurent coefficients or speed. No prior rules were imported.',
        'bindings_sha256':{str(p.relative_to(ROOT)):sha(p) for p in bound}}
    (HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    preserve={closed,proof,provisional,pp,transaction,corpus/'round-000-provisional.json',corpus/'round-000-provisional.bin'}
    entries=[]
    for path in sorted([p for p in corpus.rglob('*') if p.is_file() and p not in preserve]+[log]):
        raw=path.read_bytes();target=path.with_name(path.name+'.gz');assert not target.exists()
        target.write_bytes(gzip.compress(raw,compresslevel=9,mtime=0));assert gzip.decompress(target.read_bytes())==raw
        entries.append({'original_path':str(path.relative_to(HERE)),'archive_path':str(target.relative_to(HERE)),
            'original_sha256':hashlib.sha256(raw).hexdigest(),'archive_sha256':sha(target),'original_bytes':len(raw),'archive_bytes':target.stat().st_size})
        path.unlink()
    archive={'scope':'Completed intermediate proofs/log only; predictions, input/configuration, final closed and same-round provisional/discovery proof chain retained raw.',
        'gzip_mtime':0,'restore_checks_passed':True,'entries':entries,
        'preserved_raw':{str(p.relative_to(HERE)):{'sha256':sha(p),'bytes':p.stat().st_size}for p in sorted(preserve)},
        'file_count':len(entries),'original_bytes':sum(r['original_bytes']for r in entries),'archive_bytes':sum(r['archive_bytes']for r in entries)}
    (HERE/'archive-map.json').write_text(json.dumps(archive,indent=2)+'\n')
    mapping={r['original_path']:r for r in entries}
    for name,digest in summary['bindings_sha256'].items():
        p=ROOT/name;raw=p.read_bytes()if p.exists()else gzip.decompress((HERE/mapping[str(p.relative_to(HERE))]['archive_path']).read_bytes())
        assert hashlib.sha256(raw).hexdigest()==digest
    files={str(p.relative_to(HERE)):{'sha256':sha(p),'bytes':p.stat().st_size}for p in sorted(HERE.rglob('*'))if p.is_file()and'__pycache__'not in p.parts}
    (HERE/'artifact-manifest.json').write_text(json.dumps({'scope':summary['scope'],'files':files,'file_count':len(files),'archive_and_binding_checks':'passed'},indent=2)+'\n')
    print(json.dumps({'status':'frozen','files':len(files),'archives':len(entries),'manifest_sha256':sha(HERE/'artifact-manifest.json')}))

if __name__=='__main__':main()
