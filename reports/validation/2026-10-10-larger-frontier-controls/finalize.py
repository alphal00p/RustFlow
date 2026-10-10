#!/usr/bin/env python3
"""Verify completed bounded runs, preserve unpublished work and seal evidence."""
import gzip
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
K=HERE.parent/'2026-10-10-requested-ray-point-integration'
sys.path.insert(0,str(K))
from frozen_sources import verify_snapshot
read=lambda p:json.loads(p.read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    assert not (HERE/'artifact-manifest.json').exists()
    snapshot=K/'requested-ray-point-source-hashes.json';buildpath=K/'requested-ray-point-build-provenance.json'
    verify_snapshot(snapshot)
    build=read(buildpath);runtime,=[name for name in build['files']if Path(name).name.startswith('finite_density_runtime_flow-')]
    assert sha(ROOT/runtime)==build['files'][runtime]['sha256']
    toolpath=K/'native-digest-build.json';tool=read(toolpath);exe=Path(tool['executable']['path']);assert sha(exe)==tool['executable']['sha256']
    def b3(path):
        h,n=subprocess.check_output([str(exe),str(path)],text=True).split();assert int(n)==path.stat().st_size;return h
    bound=[snapshot,buildpath,toolpath,HERE/'run.py',HERE/'resume-assessment.md',Path(__file__)]
    preserve=set();cases=[];candidates=[]
    for arm,expected in [('single3',[58,458,1469,2023,2510,2357,1851,1301,872]),('double',[27,188,822,1820,2758])]:
        binding_path=HERE/(arm+'-binding.json');resources_path=HERE/(arm+'-resources.json')
        binding=read(binding_path);resources=read(resources_path);folder=HERE/arm;corpus=folder/'native-closure'
        assert binding['exit_code']==resources['exit_code']==124 and binding['post_run_source_and_executable_unchanged']
        assert binding['build_provenance_sha256']==sha(buildpath) and binding['source_snapshot_sha256']==sha(snapshot)
        assert binding['executable_sha256']==sha(ROOT/runtime) and binding['launcher_sha256']==sha(HERE/'run.py')
        assert binding['input_sha256']==sha(ROOT/binding['settings']['RUSTFLOW_WEIGHTED_INPUT'])
        assert not list(folder.glob('prediction-*.json')) and not (folder/'failure.json').exists()
        assert not list(corpus.glob('*-closed.*'))
        rows=sorted(corpus.glob('round-*-provisional.json'));front=[len(read(p)['frontier'])for p in rows];assert front==expected
        linked=set();rounds=[]
        for path in rows:
            row=read(path);proof=path.with_suffix('.bin');assert proof.stat().st_size==row['program_bytes']and b3(proof)==row['program_blake3']
            link=row['active_state']['requested_discovery'];transaction=corpus/link['path'];linked.add(transaction)
            assert b3(transaction)==link['blake3'];tx=read(transaction)
            assert tx['status']=='complete' and tx['native_program']['blake3']==row['program_blake3']and tx['native_program']['bytes']==row['program_bytes']
            assert tx['source_measure_id']==row['measure_id'] and tx['completed_search_implies_coverage']is False
            rounds.append({'round':row['round'],'frontier':len(row['frontier']),'native_rules':row['native_rule_count'],
                'program_sha256':sha(proof),'metadata_sha256':sha(path),'transaction_sha256':sha(transaction),
                'program_and_transaction_internal_digests_verified':True})
            bound.extend([path,proof,transaction])
        last=rows[-1];latest=read(last);lasttx=corpus/latest['active_state']['requested_discovery']['path']
        alltx=set(corpus.glob('requested-discovery-*.json'));unlinked=sorted(alltx-linked);assert len(unlinked)==1
        unpublished=[]
        for path in unlinked:
            tx=read(path);assert tx['status']=='complete'
            assert tx['native_program']['blake3'] not in {read(p)['program_blake3']for p in rows}
            unpublished.append({'path':str(path.relative_to(HERE)),'status_in_record':tx['status'],
                'native_program_descriptor':tx['native_program'],
                'scope':'Completed discovery record, but no matching published active-round program or metadata. It is not a resumable or completed closure-round state.'})
            bound.append(path)
        # Sidecars are deliberately not promoted to a committed memo prefix.
        direct=sorted(corpus.glob('direct-zero-*.json'))
        for p in direct:bound.append(p)
        parts=sorted(corpus.glob('*.part'))
        preserve.update([last,last.with_suffix('.bin'),lasttx,*unlinked,*direct,*parts,corpus/'active-request-history.json'])
        log=HERE/(arm+'-resources.log');candidates.append(log)
        bound.extend([binding_path,resources_path,HERE/(arm+'-resources.provenance.json'),log,
            folder/'input.json',folder/'configuration.json',folder/'sector.json'])
        cases.append({'arm':arm,'cut_slots':read(folder/'sector.json')['cut_slots'],'exit_code':124,
            'status':'timeout_without_prediction','wall_seconds':resources['wall_seconds'],
            'peak_child_rss_kib':resources['peak_child_rss_kib'],'timeout_seconds':1200,
            'published_provisional_rounds':len(rows),'last_published_round':latest['round'],
            'frontier_progression':front,'last_published_native_rules':latest['native_rule_count'],
            'published_round_integrity':rounds,'unpublished_discovery_records':unpublished,
            'direct_zero_records':len(direct),'partial_files':[p.name for p in parts],
            'complete_native_connection':False,'resumable_state_proved':False,'predictions_saved':0})
    for arm in ['single3','double']:
        candidates += [p for p in (HERE/arm/'native-closure').rglob('*')if p.is_file()and p not in preserve]
    summary={'status':'both_bounded_controls_timed_out','source_or_build_changes':False,
        'changes_from_K':{'frontier':[1024,4096],'historical_requests':[4096,16384],'per_round_domains':[8192,32768]},
        'other_control_changes':'Selected occupied-flow tests, one18:60:8 profile, separate report paths and1200s caps instead of the full-amplitude request. No isolated single-budget or speed inference.',
        'cases':cases,'numerical_predictions':0,'full_amplitude_acceptance':False,
        'resume_scope':'Diagnostic schema1 evidence only. See resume-assessment.md; no existing checkpoint is relabeled resumable.',
        'bindings_sha256':{str(p.relative_to(ROOT)):sha(p)for p in sorted(set(bound))}}
    (HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    entries=[]
    for p in sorted(candidates):
        raw=p.read_bytes();target=p.with_name(p.name+'.gz');assert not target.exists()
        target.write_bytes(gzip.compress(raw,compresslevel=9,mtime=0));assert gzip.decompress(target.read_bytes())==raw
        entries.append({'original_path':str(p.relative_to(HERE)),'archive_path':str(target.relative_to(HERE)),
            'original_sha256':hashlib.sha256(raw).hexdigest(),'archive_sha256':sha(target),'original_bytes':len(raw),'archive_bytes':target.stat().st_size})
        p.unlink()
    archive={'scope':'Completed historical proof payloads and process logs. Latest published program/metadata pairs, their transactions, all direct-zero sidecars, unlinked later transactions and request histories remain raw. No publication or resume status is inferred.',
        'gzip_mtime':0,'restore_checks_passed':True,'entries':entries,
        'preserved_raw':{str(p.relative_to(HERE)):{'sha256':sha(p),'bytes':p.stat().st_size}for p in sorted(preserve)},
        'file_count':len(entries),'original_bytes':sum(r['original_bytes']for r in entries),'archive_bytes':sum(r['archive_bytes']for r in entries)}
    (HERE/'archive-map.json').write_text(json.dumps(archive,indent=2)+'\n')
    mapping={r['original_path']:r for r in entries}
    for name,digest in summary['bindings_sha256'].items():
        p=ROOT/name;raw=p.read_bytes()if p.exists()else gzip.decompress((HERE/mapping[str(p.relative_to(HERE))]['archive_path']).read_bytes())
        assert hashlib.sha256(raw).hexdigest()==digest
    (HERE/'final-verification.json').write_text(json.dumps({'status':'passed','frozen_source_assets':len(read(snapshot)['files']),
        'source_snapshot_sha256':sha(snapshot),'runtime_executable_sha256':sha(ROOT/runtime),
        'historical_bound_files_verified':len(summary['bindings_sha256']),'archives_restored':len(entries),
        'scope':'Integrity and source/build verification after both timeout exits; no additional native or numerical computation.'},indent=2)+'\n')
    files={str(p.relative_to(HERE)):{'sha256':sha(p),'bytes':p.stat().st_size}for p in sorted(HERE.rglob('*'))if p.is_file()and'__pycache__'not in p.parts}
    manifest={'status':'frozen_timeouts','scope':summary['resume_scope'],'file_count':len(files),'files':files}
    (HERE/'artifact-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'status':'frozen','files':len(files),'archives':len(entries),'raw_bytes':archive['original_bytes'],
        'gzip_bytes':archive['archive_bytes'],'manifest_sha256':sha(HERE/'artifact-manifest.json')}))

if __name__=='__main__':main()
