"""Fresh single-case DEVSIM check; numerical results/logs kept separate."""
from pathlib import Path
import sys,json,contextlib,math
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'support/scripts'))
OUT=ROOT/'generated/smoke';OUT.mkdir(parents=True,exist_ok=True)
with (OUT/'solver.log').open('w') as log,contextlib.redirect_stdout(log):
    from pin_device import PIN,Config
    p=PIN(Config(thickness_um=35.6));p.bias(7.86)
    dark=p.current();p.light(1);photo=p.current()-dark;ac=p.ac(1e5)
old=json.loads((ROOT/'pin/runs/d35.6_v7.86_time0.5/summary.json').read_text())['metrics']
result={'case':'d35.6_v7.86','scope':'Fresh DEVSIM DC and 100 kHz AC only; not a transient or full receiver reproduction','photocurrent_a':photo,'capacitance_f':ac['capacitance_f'],'reference_photocurrent_a':old['photocurrent_a'],'reference_capacitance_f':old['capacitance_light_100khz_f']}
result['pass']=math.isclose(photo,old['photocurrent_a'],rel_tol=1e-6) and math.isclose(ac['capacitance_f'],old['capacitance_light_100khz_f'],rel_tol=1e-6)
(OUT/'result.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2));assert result['pass']
