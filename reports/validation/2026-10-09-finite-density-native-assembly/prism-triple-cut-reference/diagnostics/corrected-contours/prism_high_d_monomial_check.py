exec(open('/tmp/prism_triple_high_d_qmc.py').read().split('start=time.monotonic()')[0])
rows=json.load(open('/tmp/prism-triple-cut-monomials.json'))['rows'];sums=[0. for _ in rows];primes=[2,3,5,7,11,13,17,19]
for i in range(1,count+1):
 us=[radical(i,b)for b in primes]
 for surface in[False,True]:
  r1=1. if surface else us[0]**(1/p);r2=us[1]**(1/p);r3=us[2]**(1/p)
  x=2*us[3]-1;y=2*us[4]-1;z=2*us[5]-1;cos23=x*y+math.sqrt((1-x*x)*(1-y*y))*z
  hh=[2*r1*r2*(1-x),2*r1*r3*(1-y),2*r2*r3*(1-cos23)];xx=[us[6],(1-us[6])*us[7],(1-us[6])*(1-us[7])];rr=[r1,r2,r3];F=xx[0]*xx[1]*hh[0]+xx[0]*xx[2]*hh[1]+xx[1]*xx[2]*hh[2]
  weight=(1-x*x)**((p-2)/2)*(1-y*y)**((p-2)/2)*(1-z*z)**((p-3)/2)*8*normalxy**2*normalz
  weight*=math.prod(x**(b-1)for x,b in zip(xx,beta))*(1-xx[0])*math.prod(h**d for h,d in zip(hh,delta))*math.prod(r**k for r,k in zip(rr,kappa))/p**(2 if surface else 3)/math.prod(math.gamma(b)for b in beta)
  for j,row in enumerate(rows):
   if (row['sector']=='surface')!=surface:continue
   val=weight*float(Q(row['coefficient']))*math.gamma(nu+row['nu_shift'])*F**(-nu-row['nu_shift'])*math.prod(h**a for h,a in zip(hh,row['h_powers']))*math.prod(r**k for r,k in zip(rr,row['energy_powers']))*math.prod(x**a for x,a in zip(xx,row['simplex_powers']));sums[j]+=val
vals=[x/count for x in sums];old=json.load(open('/tmp/prism-barnes-high-d-numeric.json'))['rows']
for j,(actual,expected)in enumerate(zip(vals,old)):print(j,actual,expected[0],actual/expected[0]if expected[0] else 0,flush=True)
print('TOTAL',sum(vals),sum(x[0]for x in old));Path('/tmp/prism-high-d-monomials-qmc.json').write_text(json.dumps({'values':vals,'points':count},indent=2)+'\n')
