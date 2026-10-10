"""Capture only after the source owner explicitly announces freeze."""
import argparse,json,re,subprocess
from datetime import datetime,timezone
from pathlib import Path
from frozen_sources import ROOT,BASE,PREDECESSOR,digest,inventory,verify_snapshot

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix',required=True)
    parser.add_argument('--confirm-frozen',action='store_true',help='Operator confirms the source owner has announced freeze')
    parser.add_argument('--scope',required=True)
    parser.add_argument('--include-source',action='append',default=[])
    args=parser.parse_args()
    if not args.confirm_frozen:parser.error('source freeze has not been confirmed')
    if not re.fullmatch(r'[a-z0-9-]+',args.prefix):parser.error('unsafe prefix')
    output=BASE/(args.prefix+'-source-hashes.json')
    if output.exists():raise ValueError('preserve the existing source snapshot')
    names=inventory(args.include_source)
    previous=json.loads(PREDECESSOR.read_text())['files']
    report={'schema':1,'base_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'captured_utc':datetime.now(timezone.utc).isoformat(),'scope':args.scope,
        'inventory_policy':'predecessor inventory plus current Rust/Cargo sources, test Python and finite-density definitions; explicit additional source assets if supplied',
        'predecessor':{'path':str(PREDECESSOR.relative_to(ROOT)),'sha256':digest(PREDECESSOR)},
        'additional_sources':args.include_source,'added_since_predecessor':sorted(set(names)-set(previous)),
        'removed_since_predecessor':sorted(set(previous)-set(names)),
        'files':{name:digest(ROOT/name) for name in names},'capture_script_sha256':digest(Path(__file__))}
    temporary=output.with_suffix('.json.tmp');temporary.write_text(json.dumps(report,indent=2)+'\n')
    verify_snapshot(temporary);temporary.replace(output)
    print(json.dumps({'snapshot':str(output.relative_to(ROOT)),'files':len(names)}))

if __name__=='__main__':main()
