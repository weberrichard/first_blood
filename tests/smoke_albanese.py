"""Short integration and error-path checks. No long simulation is started."""
from pathlib import Path
import csv
import math
import shutil
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT/'projects/custom_output'
BASE=ROOT/'models/Bathsheba_CO2_pul_2026_03_11'
SOURCE=ROOT/'models/Bathsheba_CO2_pul'
CREATED=[]

def fixture(mode):
    path=Path(tempfile.mkdtemp(prefix='_gas_smoke_',dir=ROOT/'models'))
    CREATED.append(path)
    shutil.copytree(BASE,path,dirs_exist_ok=True)
    main=path/'main.csv'
    lines=main.read_text().splitlines()
    main.write_text('\n'.join('time,0.05' if l.startswith('time,') else l for l in lines)+'\n')
    for name in ['O2_parameters.txt','CO2_parameters.txt']:
        if not (path/name).exists(): shutil.copyfile(SOURCE/name,path/name)
    params=path/'O2_parameters.txt'
    params.write_text(params.read_text()+'\ngas_model,'+mode+'\n')
    (path/'custom_output.txt').write_text('save_dt,0.005\np8,lum,transport\nheart_kim_lit,lum,transport\n')
    return path

def run(path,success=True):
    result=subprocess.run(['./custom_output.out',path.name],cwd=WORK,
                          capture_output=True,text=True,timeout=45)
    if success:
        assert result.returncode==0,(result.returncode,result.stdout,result.stderr)
        assert '+ run: OK' in result.stdout
    else:
        assert result.returncode==1,(result.returncode,result.stdout,result.stderr)
        assert 'custom_output:' in result.stderr
        assert '+ run: OK' not in result.stdout
    return result

try:
    fingerprints={}
    for mode in ['legacy','albanese2016']:
        path=fixture(mode)
        run(path)
        files=list((WORK/'results'/path.name).rglob('*.txt'))
        assert len(files)==136,(mode,len(files))
        for file in files:
            rows=[[float(v) for v in row] for row in csv.reader(file.open()) if row]
            assert len(rows)>=2 and rows[0][0]==0. and rows[-1][0]>.035,str(file)
            assert all(len(row)==3 for row in rows),str(file)
            assert all(math.isfinite(v) for row in rows for v in row),str(file)
            assert all(rows[i][0]>rows[i-1][0] for i in range(1,len(rows))),str(file)
            if file.parent.name=='RBC':
                assert all(abs(v/4.9e15-1)<1e-6 for row in rows for v in row[1:]),str(file)
            elif file.parent.name=='HB_O2':
                assert all(0<=v<=1 for row in rows for v in row[1:]),str(file)
            else:
                assert all(v>=0 for row in rows for v in row[1:]),str(file)
        hb=WORK/'results'/path.name/'p8/HB_O2/capillary.txt'
        fingerprints[mode]=[float(v) for v in list(csv.reader(hb.open()))[-1]][1:]
        print(mode+': 0.05 s integration and 136 output histories passed',flush=True)
    assert max(abs(a-b) for a,b in zip(fingerprints['legacy'],fingerprints['albanese2016']))>1.e-6, 'Albanese selector did not change sources'
    path=fixture('legacy')
    params=path/'CO2_parameters.txt'
    params.write_text(params.read_text()+'\ngas_model,albanese2016\n')
    run(path)
    hb=WORK/'results'/path.name/'p8/HB_O2/capillary.txt'
    carbon_file_result=[float(v) for v in list(csv.reader(hb.open()))[-1]][1:]
    assert all(abs(a-b)<1.e-10 for a,b in zip(carbon_file_result,fingerprints['albanese2016'])), 'CO2 file model selection differs from O2 file'
    print('Both parameter files and changed source response verified',flush=True)
    path=fixture('albanese2016')
    params=path/'O2_parameters.txt'
    params.write_text(params.read_text()+'\nalbanese_o2,9,0,0.03198,0.008275,14.99\n')
    assert 'invalid dissociation' in run(path,False).stderr
    path=fixture('albanese2016')
    (path/'CO2_parameters.txt').unlink()
    assert 'perifco2' in run(path,False).stderr
    print('Malformed and missing parameter failures passed',flush=True)
finally:
    for path in CREATED:
        assert path.parent.resolve()==(ROOT/'models').resolve() and path.name.startswith('_gas_smoke_')
        shutil.rmtree(path)
        results=WORK/'results'/path.name
        if results.exists(): shutil.rmtree(results)
