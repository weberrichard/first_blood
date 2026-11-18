"""Full-network smoke tests, each at most 5 ms of simulated time."""
from pathlib import Path
import csv,json,math,shutil,subprocess,tempfile

ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT/'projects/custom_output'
BASE=ROOT/'models/Bathsheba_CO2_pul_2026_03_11'
REPORT=Path('/home/wr/first_blood-runs/dev_gas_conservative_20261118/smoke.json')
created=[];report=[]
def fixture(duration):
    p=Path(tempfile.mkdtemp(prefix='_gas_implicit_ms_',dir=ROOT/'models'));created.append(p)
    shutil.copytree(BASE,p,dirs_exist_ok=True)
    main=p/'main.csv'
    main.write_text('\n'.join('time,'+str(duration) if l.startswith('time,') else l for l in main.read_text().splitlines())+'\n')
    (p/'custom_output.txt').write_text('save_dt,0.00025\np8,lum,transport\nheart_kim_lit,lum,transport\n')
    return p

def run(p,duration,success=True,label='implicit'):
    result=subprocess.run(['./custom_output.out',p.name],cwd=WORK,capture_output=True,text=True,timeout=60)
    if not success:
        assert result.returncode==1,(label,result.returncode,result.stdout[-500:],result.stderr)
        assert 'custom_output:' in result.stderr
        assert '+ run: OK' not in result.stdout
        report.append(dict(case=label,rejected=True,message=result.stderr.strip()))
        return
    assert result.returncode==0,(label,result.returncode,result.stdout[-1000:],result.stderr)
    assert '+ run: OK' in result.stdout
    files=list((WORK/'results'/p.name).rglob('*.txt'))
    assert len(files)==136,(label,len(files))
    signature={};last_time=0.
    for f in files:
        rows=[[float(v) for v in r] for r in csv.reader(f.open()) if r]
        assert len(rows)>=2 and rows[0][0]==0. and rows[-1][0]>=duration-.00051,(str(f),rows[-1])
        assert all(len(r)==3 and all(math.isfinite(x) for x in r) for r in rows),str(f)
        assert all(rows[i][0]>rows[i-1][0] for i in range(1,len(rows))),str(f)
        if f.parent.name=='HB_O2':assert all(0.<=x<=1. for r in rows for x in r[1:]),str(f)
        else:assert all(x>=0. for r in rows for x in r[1:]),str(f)
        if f.parent.name=='RBC':assert all(abs(x/4.9e15-1.)<1.e-6 for r in rows for x in r[1:]),str(f)
        signature[str(f.relative_to(WORK/'results'/p.name))]=rows[-1][1:]
        last_time=max(last_time,rows[-1][0])
    report.append(dict(case=label,duration_s=duration,histories=len(files),last_saved_s=last_time,finite_nonnegative=True))
    print(label+': '+str(duration*1000)+' ms; 136 histories finite, nonnegative and bounded',flush=True)
    return signature

try:
    run(fixture(.002),.002,label='implicit_2ms')
    baseline=run(fixture(.005),.005,label='implicit_5ms')
    p=fixture(.005)
    o=p/'O2_parameters.txt';c=p/'CO2_parameters.txt'
    lines=o.read_text().splitlines();gas_rows=[l for l in lines if l.startswith('gas_')]
    o.write_text('\n'.join(l for l in lines if not l.startswith('gas_'))+'\n')
    c.write_text(c.read_text()+'\n'+'\n'.join(gas_rows)+'\n')
    moved=run(p,.005,label='parameters_from_CO2_file')
    assert baseline.keys()==moved.keys()
    assert all(abs(a-b)<1.e-10*max(1.,abs(a)) for key in baseline for a,b in zip(baseline[key],moved[key])),'Parameter file location changed results'
    for mode in ['legacy','albanese2016']:
        p=fixture(.002);o=p/'O2_parameters.txt'
        o.write_text(o.read_text()+'\ngas_solver,legacy\ngas_model,'+mode+'\n')
        run(p,.002,label='regression_'+mode)
    for key,replacement,label in [
        ('gas_relaxation,',None,'missing_relaxation'),
        ('gas_partition,','gas_partition,1','malformed_partition'),
        ('gas_relaxation,','gas_relaxation,.7,.7,0,.3,.001,.07','zero_exchange_time'),
        ('gas_solver,','gas_solver,unsupported','unknown_solver')]:
        p=fixture(.002);o=p/'O2_parameters.txt'
        o.write_text('\n'.join((replacement or '') if l.startswith(key) else l for l in o.read_text().splitlines())+'\n')
        run(p,.002,False,label)
    REPORT.write_text(json.dumps(report,indent=2))
    print('All millisecond integration, regression and input failure checks passed',flush=True)
finally:
    for p in created:
        assert p.parent.resolve()==(ROOT/'models').resolve() and p.name.startswith('_gas_implicit_ms_')
        shutil.rmtree(p)
        results=WORK/'results'/p.name
        if results.exists():shutil.rmtree(results)
