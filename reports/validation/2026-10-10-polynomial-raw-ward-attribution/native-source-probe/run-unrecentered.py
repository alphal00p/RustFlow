import pathlib,json,subprocess,sys
root=pathlib.Path.cwd();p=root/'reports/validation/2026-10-10-polynomial-raw-ward-attribution/native-source-probe';B=root/'reports/validation/2026-10-09-finite-density-native-assembly'
commands={}
for mode in ['original','ward-only','boost-ward','ward-boost']:
 commands[mode]=[sys.executable,str(B/'run-resource-command.py'),str(p/('unrecentered-'+mode+'-resources.json')),'timeout','90s','env','RAW_WARD_RECENTER=0','/tmp/rustflow-raw-ward-probe',str(B/'three-loop-singleton-source-diagnostic/round-006-provisional.bin'),str(B/'three-loop-singleton-source-diagnostic/points.json'),str(root/'examples/finite_density/massless_three_loop_chain.json'),str(p/('unrecentered-'+mode)),mode]
(p/'unrecentered-planned-commands.json').write_text(json.dumps(commands,indent=2)+'\n')
for mode,cmd in commands.items():
 r=subprocess.run(cmd);print(mode,r.returncode,flush=True)
 if r.returncode:sys.exit(r.returncode)
