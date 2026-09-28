from pathlib import Path
import shutil, subprocess, tempfile, sys, re, json, csv
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'src'/'pic_2010_recovered.cpp'
CFG0=ROOT/'config'/'validation_quick_collisionless.txt'
CFG1=ROOT/'config'/'validation_quick_collisional.txt'

def parse_h4(path):
    lines=path.read_text(errors='replace').strip().splitlines()
    if len(lines)<2: raise RuntimeError('missing fitting result')
    return float(lines[-1].split()[0])

def run_case(exe,cfg,td,name,passes=3):
    d=td/name; d.mkdir()
    txt=cfg.read_text(encoding='utf-8')
    txt=re.sub(r'pass_loop\s*=\s*3', f'pass_loop= {passes}', txt)
    (d/'input.txt').write_text(txt,encoding='utf-8')
    cp=subprocess.run([str(exe)],cwd=d,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE,text=True,timeout=120)
    if cp.returncode: raise RuntimeError(f'{name} exit={cp.returncode}: {cp.stderr[-1000:]}')
    out=(d/'output.txt').read_text(errors='replace')
    for needle in ['collision_alpha=3','P_weighting=1',f'pass_loop={passes}']:
        if needle not in out: raise RuntimeError(f'{name}: output missing {needle}')
    return parse_h4(d/'fitting_number(omega_p*tau_R).txt')

def main():
    compiler=shutil.which('g++')
    if not compiler: raise SystemExit('THESIS_VALIDATION_NEEDS_GXX')
    with tempfile.TemporaryDirectory(prefix='pic_thesis_validate_') as x:
        td=Path(x); exe=td/('pic_validate.exe' if sys.platform.startswith('win') else 'pic_validate')
        subprocess.run([compiler,'-std=c++17','-O2',str(SRC),'-o',str(exe)],check=True,stdout=subprocess.DEVNULL)
        t0=run_case(exe,CFG0,td,'collisionless',3)
        t1=run_case(exe,CFG1,td,'collisional',3)
        passes=3
        # User-approved protocol: default 3; only escalate to 5 if stochastic direction is unclear.
        if not (t1<t0):
            t0=run_case(exe,CFG0,td,'collisionless_5',5)
            t1=run_case(exe,CFG1,td,'collisional_5',5)
            passes=5
        result={'status':'PASS' if t1<t0 else 'FAIL','alpha':3,'weighting':1,'passes':passes,'collisionless_H4_tau':t0,'collisional_H4_tau':t1,'ratio_collisional_to_collisionless':t1/t0}
        out=ROOT/'validation'/'latest_one_click_validation.json'; out.write_text(json.dumps(result,indent=2),encoding='utf-8')
        if result['status']!='PASS': raise SystemExit('THESIS_VALIDATION_FAIL: representative collision direction not reproduced')
        print('THESIS_VALIDATION_PASS')
        print(f"alpha=3 weighting=1 passes={passes} collisionless_H4_tau={t0:.6g} collisional_H4_tau={t1:.6g} ratio={t1/t0:.4f}")
if __name__=='__main__': main()
