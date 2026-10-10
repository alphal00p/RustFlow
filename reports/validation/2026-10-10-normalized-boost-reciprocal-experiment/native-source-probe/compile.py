import pathlib,subprocess,json,hashlib,datetime,os,sys,time
root=pathlib.Path.cwd();p=root/'reports/validation/2026-10-10-normalized-boost-reciprocal-experiment/native-source-probe';name=sys.argv[1]
deps={'symbolica':root/'target/release/deps/libsymbolica-02ed301b2d4bd7f2.rlib','serde_json':root/'target/release/deps/libserde_json-84255895f7dd0f90.rlib'}
deps['symbolica_amflow']=root/'target/release/deps/libsymbolica_amflow-bc386786e4706ace.rlib'
if True:
 deps['rustred']=pathlib.Path('/tmp/rustflow-zero-projection-native-library/librustred-392075d0c6995fd1.rlib');deps['bincode']=root/'target/release/deps/libbincode-30fae4fd0bc9d1b5.rlib'
meta=lambda f:{'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
source=p/(name+'.rs');exe=pathlib.Path('/tmp/rustflow-reciprocal-global-boost-'+name)
cmd=['rustc','--edition=2024','--crate-name','reciprocal_'+name,'-C','opt-level=0','-C','debuginfo=0']
for k,v in deps.items():cmd+=['--extern',k+'='+str(v)]
cmd+=['-L','dependency='+str(root/'target/release/deps')]
for d in ['gmp-mpfr-sys-2e6d668657be102d','gmp-mpfr-sys-c789de89ed67d3db']:cmd+=['-L','native='+str(root/'target/release/build'/d/'out/lib')]
cmd += [str(source),'-o',str(exe)]
bind={'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':meta(source),'dependencies':{k:meta(v) for k,v in deps.items()},'command':cmd,'cargo_invoked':False};(p/(name+'-build-binding.json')).write_text(json.dumps(bind,indent=2)+'\n');start=time.monotonic()
with (p/(name+'-build.log')).open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env={**os.environ,'CARGO_CRATE_NAME':'reciprocal_'+name})
out={'exit_code':r.returncode,'wall_seconds':time.monotonic()-start}
if r.returncode==0:out['executable']=meta(exe)
(p/(name+'-build.json')).write_text(json.dumps(out,indent=2)+'\n');print(out);sys.exit(r.returncode)
