"""Fresh receiver AC check using preserved device operating-point parameters."""
from pathlib import Path
import sys,json,math,os
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'support/scripts'))
import tia_interface as ti
ti.SHARED_LOCK=Path(os.environ.get('LTSPICE_LOCK',str(ROOT/'.ltspice_serial.lock')))
source=json.loads((ROOT/'pin/receiver_runs/d35.6_v7.86/ac/summary.json').read_text())
case={k:v for k,v in source['conditions'].items() if k in ['case_id','capacitance_f','photocurrent_a','dark_current_a','responsivity_a_w','reverse_bias_v','absorber_thickness_um']}
settings={**ti.DEFAULTS,**{k:source['conditions'][k] for k in ti.DEFAULTS if k in source['conditions']}}
out=ROOT/'generated/smoke_receiver'
result=ti.run_case(case,settings,'ac',(ROOT/'tia/original_tia.lib').read_text(),out,force=True)
check={'scope':'Fresh LTspice AC using preserved device operating point, not full device/transient reproduction','upper_3db_hz':result['upper_3db_hz'],'reference_upper_3db_hz':source['upper_3db_hz'],'pass':math.isclose(result['upper_3db_hz'],source['upper_3db_hz'],rel_tol=1e-5)}
(out/'verification.json').write_text(json.dumps(check,indent=2));print(json.dumps(check,indent=2));assert check['pass']
