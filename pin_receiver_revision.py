"""Simulate the TIA readout and analyse its coupling to generated PIN responses."""
from pathlib import Path
import sys,json,csv,hashlib,argparse,contextlib
import numpy as np
ROOT=Path(__file__).resolve().parent;WORK=ROOT.parent
sys.dont_write_bytecode=True
HELPERS=ROOT/'support/scripts'
sys.path.insert(0,str(HELPERS))
import tia_interface as ti
ti.SHARED_LOCK=ROOT/'.ltspice_serial.lock'
from pin_device import PIN,Config,write_csv
from run_pin_revision import latest,DS,VS,OUT

def read(p):return json.loads(Path(p).read_text())
def cases():
    selected=read(ROOT/'tia/selection.json')
    settings={**ti.DEFAULTS,'feedback_resistance_ohm':selected['feedback_resistance_ohm'],
              'feedback_capacitance_f':selected['feedback_capacitance_f'],
              'temperature_c':26.85,'transient_max_step_s':.25e-9}
    rows=[]
    for d,v in [(DS[0],VS[1]),(DS[1],VS[1])]+[(DS[2],v) for v in VS]:
        _,folder,detail=latest(d,v);m=detail['metrics']
        rows.append(dict(case_id=f'd{d:g}_v{v:g}',capacitance_f=m['capacitance_light_100khz_f'],
            photocurrent_a=m['photocurrent_a'],dark_current_a=m['dark_current_a'],
            responsivity_a_w=m['responsivity_a_w'],reverse_bias_v=v,absorber_thickness_um=d,
            device_summary_path=str(folder/'summary.json'),
            device_summary_sha256=hashlib.sha256((folder/'summary.json').read_bytes()).hexdigest(),
            optical_pulse_csv=str(folder/'photocurrent_replay.csv')))
    return rows,settings

def simulate(case_ids=None):
    rows,settings=cases(); model=(ROOT/'tia/original_tia.lib').read_text()
    if case_ids:
        rows=[c for c in rows if c['case_id'] in case_ids]
    (OUT/'receiver_inputs.json').write_text(json.dumps({'settings':settings,'cases':rows},indent=2))
    for c in rows:
        for analysis in ['op','ac','noise']:
            ti.run_case(c,settings,analysis,model,OUT/'receiver_runs')
        pulse=Path(c['optical_pulse_csv'])
        with pulse.open() as handle:first=next(csv.DictReader(handle))
        tc={**c,'photocurrent_a':float(first['photocurrent_a'])}
        ti.run_case(tc,settings,'tran',model,OUT/'receiver_runs',pulse)
    # Independent circuit time discretisation check on slow detector case.
    slow=next((c for c in rows if c['absorber_thickness_um']==DS[2] and c['reverse_bias_v']==VS[0]),None)
    if slow:
        with Path(slow['optical_pulse_csv']).open() as h:first=next(csv.DictReader(h))
        tc={**slow,'photocurrent_a':float(first['photocurrent_a'])}
        ti.run_case(tc,{**settings,'transient_max_step_s':.125e-9},'tran',model,OUT/'receiver_halfstep',Path(slow['optical_pulse_csv']))
    ti.export_metrics_csv(OUT/'receiver_runs')

def admittance():
    rows,_=cases(); dest=OUT/'dense_admittance';dest.mkdir(exist_ok=True)
    for c in rows[1:]:
        path=dest/(c['case_id']+'.csv')
        if path.exists():continue
        with (dest/(c['case_id']+'.log')).open('w') as log,contextlib.redirect_stdout(log):
            p=PIN(Config(thickness_um=c['absorber_thickness_um']));p.bias(c['reverse_bias_v']);p.light(1)
            freq=np.unique(np.r_[np.geomspace(1e3,2e7,31),2e6,4e6,7e6,8e6])
            records=[{'case_id':c['case_id'],'illumination_scale':1,**p.ac(float(f))} for f in freq]
            write_csv(path,records)
        print('ADMITTANCE_DONE',c['case_id'],flush=True)
    if all((OUT/'receiver_runs'/c['case_id']/'ac/test.raw').exists() for c in rows[1:]):
        ti.audit_dense_admittance(dest,OUT/'receiver_runs',OUT/'admittance_audit.csv')
        # The helper reports the sampled band through 4 MHz. Also summarise
        # directly simulated Y points through 8 MHz, including the explicit
        # 8 MHz sample.
        with (OUT/'admittance_audit.csv').open() as handle:
            full=list(csv.DictReader(handle))
        audit=[]
        for c in rows[1:]:
            sample=[r for r in full if r['case_id']==c['case_id'] and float(r['frequency_hz'])<=8e6]
            worst=max(sample,key=lambda r:float(r['corrected_transfer_complex_change_pct']))
            audit.append({'case_id':c['case_id'],'sampled_low_hz':min(float(r['frequency_hz']) for r in sample),
                'sampled_high_hz':max(float(r['frequency_hz']) for r in sample),
                'max_complex_transfer_change_pct':float(worst['corrected_transfer_complex_change_pct']),
                'worst_frequency_hz':float(worst['frequency_hz']),
                'max_absolute_magnitude_change_pct':max(abs(float(r['corrected_transfer_magnitude_change_pct'])) for r in sample),
                'scope':'CALCULATED linear loading correction, maximum among directly simulated device Y points <=8MHz; not continuous-band proof or nonlinear cosimulation'})
        (OUT/'admittance_audit_8mhz.json').write_text(json.dumps(audit,indent=2))
        print('ADMITTANCE_AUDIT_8MHZ',json.dumps(audit),flush=True)
    else:
        print('DEVICE_ADMITTANCE_SAVED; transfer audit awaits actual receiver AC runs.',flush=True)

def upper(f,a):
    level=a[0]/np.sqrt(2);hit=np.flatnonzero(a<=level)
    if not len(hit):return None
    i=hit[0]
    if i==0:return float(f[0])
    frac=(level-a[i-1])/(a[i]-a[i-1]);return float(np.exp(np.log(f[i-1])+frac*np.log(f[i]/f[i-1])))

def summarise():
    rows,settings=cases();result=[]
    for c in rows:
        d,v=c['absorber_thickness_um'],c['reverse_bias_v'];_,folder,detail=latest(d,v)
        t=detail['transient'];s=np.genfromtxt(folder/'optical_step.csv',delimiter=',',names=True)
        norm=(s['total_current_a']-t['initial_current_a'])/(t['final_current_a']-t['initial_current_a'])
        ac=ti.parse_raw(OUT/'receiver_runs'/c['case_id']/'ac/test.raw')
        idx=(ac['frequency'].real>=1e3)&(ac['frequency'].real<=1e8)
        f=ac['frequency'].real[idx];z=ac['v(out)'][idx]
        mid=(s['time_s'][1:]+s['time_s'][:-1])/2
        h=np.exp(-2j*np.pi*f[:,None]*mid[None,:])@np.diff(norm)
        tr=read(OUT/'receiver_runs'/c['case_id']/'tran/summary.json')
        rr={'case_id':c['case_id'],'device_rise_ns':t['rise_10_90_s']*1e9,
            'device_settling_ns':t['settling_2pct_s']*1e9,
            'tia_electrical_bw_mhz':upper(f,np.abs(z))/1e6,
            'receiver_optical_bw_mhz':upper(f,np.abs(h*z))/1e6,
            'optical_bandwidth_method':'CALCULATED linear cascade of simulated optical step increments and SPICE AC; not nonlinear TCAD/SPICE cosimulation',
            'transient_summary':tr}
        result.append(rr)
        write_csv(OUT/f"optical_transfer_{c['case_id']}.csv",[{'frequency_hz':float(x),'device_abs':float(abs(y)),'receiver_abs_ohm':float(abs(zz*y))} for x,y,zz in zip(f,h,z)])
    (OUT/'receiver_summary.json').write_text(json.dumps(result,indent=2))
    flat=[{k:v for k,v in r.items() if not isinstance(v,dict)} for r in result]
    write_csv(OUT/'receiver_summary.csv',flat)
    print(json.dumps(flat,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=['simulate','admittance','summary'],required=True)
    ap.add_argument('--case-ids',nargs='*');a=ap.parse_args()
    if a.mode=='simulate':simulate(a.case_ids)
    else:{'admittance':admittance,'summary':summarise}[a.mode]()
