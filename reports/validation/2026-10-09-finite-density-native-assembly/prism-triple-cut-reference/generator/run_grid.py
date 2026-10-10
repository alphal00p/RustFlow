"""Independent reference grid; native arithmetic executable, no AMF predictions."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import argparse,gzip,hashlib,json,subprocess,sys,time
p=argparse.ArgumentParser();p.add_argument('executable');p.add_argument('output');p.add_argument('--denominator',type=int,default=10000);p.add_argument('--count',type=int,default=12);p.add_argument('--digits',type=int,default=70);p.add_argument('--step',default='1/30');p.add_argument('--cutoff',default='14');p.add_argument('--sample-timeout',type=int,default=300);p.add_argument('--start',type=int,default=1);args=p.parse_args()
here=Path(__file__).resolve().parent;repo=next(x for x in here.parents if (x/'Cargo.toml').is_file());out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=True)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
exe=Path(args.executable).resolve();fixture=repo/'examples/finite_density/triangular_prism.json';definition=json.loads(fixture.read_text());assert definition['targets'][1]=={'powers':[1,2,1,1,1,1,1,1,1],'numerator':'g1_2^2+g1_3*g2_4'};assert definition['chemical_potentials']==['1'];assert definition['numerator_convention']=='shifted_euclidean';assert all(e['mass_squared']=='0' for e in definition['edges']);assert [e['routing']for e in definition['edges']]==[[str(x)for x in v]for v in [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,0,0,-1],[0,1,0,-1],[0,0,1,-1],[1,-1,0,0],[1,0,-1,0]]]
meta={'status':'running independent reference grid','configuration':vars(args),'executable_sha256':digest(exe),'source_sha256':{str(x.relative_to(repo)):digest(x)for x in [*here.glob('*.py'),here/'prism_barnes_full_reference.rs',here.parent/'monomials.json',fixture,repo/'fixtures/finite_density/prism_virtual_candidate_relations.json']},'native_amf_predictions_read':0,'oracle_answers_read':0,'scope':'supplemental original raised numerator prism target only; empirical reference refinement pending','completed_samples':[]};save(out/'provenance.json',meta)
started=time.monotonic()
for i in range(args.start,args.count+1):
 stem=f'sample-{i:02d}';eps=str(F(i,args.denominator));raw=out/(stem+'-initial.json');cont=out/(stem+'-contours.json');res=out/(stem+'.json');t=time.monotonic()
 if res.exists():print('preserving existing',res,flush=True);continue
 with (out/(stem+'-construction.log')).open('w')as log:
  for cmd in [[sys.executable,str(here/'prism_barnes_symbolic_regulator.py'),eps,'1/100003',str(raw)],[sys.executable,str(here/'prism_barnes_symbolic_optimize.py'),str(raw),str(cont)],[sys.executable,str(here/'prism_barnes_regulator_analysis.py'),str(cont),str(out/(stem+'-regulator.json'))],[sys.executable,str(here/'prism_barnes_tail_check.py'),str(cont),str(out/(stem+'-tails.json'))]]:
   subprocess.run(cmd,check=True,cwd=repo,stdout=log,stderr=subprocess.STDOUT,timeout=args.sample_timeout)
 # Resource wrapper hashes the actual executable and exact continued input before launch.
 wrapper=repo/'reports/validation/2026-10-09-finite-density-native-assembly/run-resource-command.py'
 subprocess.run([sys.executable,str(wrapper),str(out/(stem+'-resources.json')),'timeout',str(args.sample_timeout),str(exe),str(cont),str(args.digits),args.step,args.cutoff,str(res)],check=True,cwd=repo,timeout=args.sample_timeout+60,stdout=subprocess.DEVNULL)
 archives=[]
 for x in [raw,cont]:
  sha=digest(x);dest=x.with_suffix(x.suffix+'.gz');gzip.open(dest,'wb').write(x.read_bytes());archives.append({'file':dest.name,'sha256_uncompressed':sha,'sha256_gzip':digest(dest)});x.unlink()
 save(out/(stem+'-archives.json'),archives);meta['completed_samples'].append({'epsilon':eps,'file':res.name,'sha256':digest(res),'elapsed_seconds':time.monotonic()-t});save(out/'provenance.json',meta);print('completed',i,eps,'seconds',time.monotonic()-t,flush=True)
meta['status']='sample grid complete; refinements and Laurent acceptance separate';meta['wall_seconds']=time.monotonic()-started;save(out/'provenance.json',meta)
# Exact rational Lagrange weights; only the already-saved numerical values enter the fit.
with localcontext()as c:
 c.prec=args.digits+40;xs=[];ys=[]
 for i in range(1,args.count+1):
  d=json.loads((out/f'sample-{i:02d}.json').read_text());assert F(d['epsilon'])==F(i,args.denominator)
  value=d['full_msbar_reference_normalization'][1:-1];real,imag=value.split('+',1);assert imag.endswith('i');real=Decimal(real);imag=Decimal(imag[:-1]);xs.append(i);ys.append((real,imag))
 def decimal(q):return Decimal(q.numerator)/Decimal(q.denominator)
 coefficients=[(Decimal(0),Decimal(0))for _ in xs]
 for i,(vr,vi)in zip(xs,ys):
  poly=[F(1)];den=F(1)
  for j in xs:
   if i==j:continue
   nxt=[F(0)]*(len(poly)+1)
   for k,v in enumerate(poly):nxt[k]-=j*v;nxt[k+1]+=v
   poly=nxt;den*=i-j
  scale=decimal(F(i,args.denominator)**6)
  for k,w in enumerate(poly):
   wr=decimal(w/den*args.denominator**k)*scale;ar,ai=coefficients[k];coefficients[k]=(ar+wr*vr,ai+wr*vi)
 report={'status':'Laurent interpolation only; no acceptance until independent refinements','working_pole_floor':-6,'requested_orders':[-4,0],'normalization':d['normalization'],'coefficients':{str(k-6):{'re':str(a),'im':str(b)}for k,(a,b)in enumerate(coefficients)},'configuration':vars(args),'sample_hashes':{f'sample-{i:02d}.json':digest(out/f'sample-{i:02d}.json')for i in xs},'arithmetic':'exact rational Lagrange weights and Decimal precision '+str(c.prec)};save(out/'laurent-fit.json',report);print(json.dumps(report['coefficients'],indent=2),flush=True)
