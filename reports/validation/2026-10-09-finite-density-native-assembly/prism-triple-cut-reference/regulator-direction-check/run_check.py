from pathlib import Path
import shutil,subprocess,sys,json,hashlib,gzip
root=Path('reports/validation/2026-10-09-finite-density-native-assembly/prism-triple-cut-reference').resolve();out=root/'regulator-direction-check';out.mkdir(exist_ok=True);gen=out/'generator';gen.mkdir(exist_ok=True)
for f in (root/'generator').glob('*.py'):shutil.copyfile(f,gen/f.name)
shutil.copyfile(root/'monomials.json',out/'monomials.json')
p=gen/'prism_barnes_symbolic_regulator.py';s=p.read_text().replace('[7,11,13]','[41,43,47]').replace('[17,19,23]','[53,59,61]').replace('[29,31,37]','[67,71,73]');p.write_text(s)
initial=out/'initial.json';cont=out/'contours.json'
with (out/'construction.log').open('w')as log:
 for cmd in [[sys.executable,str(gen/'prism_barnes_symbolic_regulator.py'),'1/101','1/100003',str(initial)],[sys.executable,str(gen/'prism_barnes_symbolic_optimize.py'),str(initial),str(cont)],[sys.executable,str(gen/'prism_barnes_regulator_analysis.py'),str(cont),str(out/'regulator.json')],[sys.executable,str(gen/'prism_barnes_tail_check.py'),str(cont),str(out/'tails.json')]]:subprocess.run(cmd,check=True,stdout=log,stderr=subprocess.STDOUT,timeout=120)
wrapper=root.parent/'run-resource-command.py';subprocess.run([sys.executable,str(wrapper),str(out/'resources.json'),'timeout','300','/tmp/prism_barnes_full_reference_v2',str(cont),'55','1/30','14',str(out/'reference.json')],check=True)
archives=[]
for p in [initial,cont]:
 blob=p.read_bytes();sha=hashlib.sha256(blob).hexdigest();dest=p.with_suffix(p.suffix+'.gz');gzip.open(dest,'wb').write(blob);archives.append({'file':dest.name,'uncompressed_sha256':sha,'gzip_sha256':hashlib.sha256(dest.read_bytes()).hexdigest()});p.unlink()
(out/'archives.json').write_text(json.dumps(archives,indent=2)+'\n')
