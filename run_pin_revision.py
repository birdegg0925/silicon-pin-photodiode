"""Run the PIN grid, time-step refinement and validation in generated/pin/."""
from pathlib import Path
import sys, json, traceback, csv, hashlib, argparse, shutil, contextlib
import numpy as np

ROOT=Path(__file__).resolve().parent
WORK=ROOT.parent
sys.dont_write_bytecode=True
HELPERS=ROOT/'support/scripts'
sys.path.insert(0,str(HELPERS))
from run_study import case
from pin_device import Config, PIN, Q, write_csv

OUT=ROOT/'generated/pin'
PROTOCOL=json.loads((ROOT/'protocol.json').read_text())['pin']
DS=PROTOCOL['thickness_um']; VS=PROTOCOL['reverse_bias_v']

def ident(d,v): return f'd{d:g}_v{v:g}'
def read(p):return json.loads(p.read_text())

def execute(d,v,name,transient=False,time_scale=1,mesh_scale=1):
    folder=OUT/'runs'/name
    if (folder/'summary.json').exists():
        print('PRESERVED',name,flush=True);return read(folder/'summary.json')
    if (folder/'failure.txt').exists():
        print('PRESERVED_FAILURE',name,flush=True);return None
    folder.mkdir(parents=True,exist_ok=True)
    try:
        case(Config(thickness_um=d,mesh_scale=mesh_scale),v,folder,
             transient=transient,time_scale=time_scale)
        return read(folder/'summary.json')
    except Exception:
        (folder/'failure.txt').write_text(traceback.format_exc())
        print('FAILED',name,flush=True);return None

def latest(d,v):
    found=[]
    for p in (OUT/'runs').glob(ident(d,v)+'*/summary.json'):
        r=read(p)
        if r['config']['mesh_scale']==1 and 'transient' in r:
            found.append((r['transient']['time_scale'],p.parent,r))
    return min(found,key=lambda r:r[0]) if found else None

def initial():
    for d in DS:
        for v in VS:execute(d,v,ident(d,v),True,1)

def refine():
    # Only middle/thick cases qualify for responsivity target or headline timing.
    for d in DS[1:]:
        for v in VS:
            old=read(OUT/'runs'/ident(d,v)/'summary.json')
            for scale in [.5,.25,.125]:
                new=execute(d,v,ident(d,v)+f'_time{scale:g}',True,scale)
                if not new:break
                keys=['rise_10_90_s','settling_2pct_s']
                delta={k:100*(new['transient'][k]-old['transient'][k])/old['transient'][k]
                       if old['transient'][k] and new['transient'][k] else None for k in keys}
                verdict={'case_id':ident(d,v),'time_scale':scale,'comparison_time_scale':old['transient']['time_scale'],
                         'differences_pct':delta,'pass':all(x is not None and abs(x)<=2 for x in delta.values()),
                         'definition':'Both10-90rise and sampled+/-2%-of-optical-increment settling must change<=2%.'}
                (OUT/'runs'/(ident(d,v)+f'_time{scale:g}')/'refinement_check.json').write_text(json.dumps(verdict,indent=2))
                print('CONVERGENCE',json.dumps(verdict),flush=True)
                if verdict['pass']:break
                old=new

def summarise():
    rows=[]
    for d in DS:
        for v in VS:
            item=latest(d,v)
            if not item:
                rows.append({'case_id':ident(d,v),'status':'failed'});continue
            scale,folder,r=item;m=r['metrics'];t=r['transient']
            audit_path=folder/'refinement_check.json'
            converged=read(audit_path)['pass'] if audit_path.exists() else False
            row={'case_id':ident(d,v),'status':'simulated','thickness_um':d,'reverse_bias_v':v,
                 'responsivity_a_w':m['responsivity_a_w'],'capacitance_light_pf':m['capacitance_light_100khz_f']*1e12,
                 'photocurrent_ua':m['photocurrent_a']*1e6,'dark_current_na':m['dark_current_a']*1e9,
                 'rise_ns':t['rise_10_90_s']*1e9 if t['rise_10_90_s'] else None,
                 'settling_ns':t['settling_2pct_s']*1e9 if t['settling_2pct_s'] else None,
                 'time_scale':scale,'timing_convergence_pass':converged,
                 'collected_fraction_of_absorbed_photons':m['collected_fraction_of_absorbed_photons'],
                 'max_field_v_cm':m['max_field_v_cm'],'terminal_current_imbalance_a':m['current_imbalance_a'],
                 'transient_final_relative_error':t['final_relative_error'],
                 'source_folder':str(folder.relative_to(OUT))}
            row['physics_checks_pass']=(0<row['collected_fraction_of_absorbed_photons']<=1.002 and
                row['max_field_v_cm']<5e4 and abs(m['current_imbalance_a'])<1e-5*(abs(m['photocurrent_a'])+abs(m['dark_current_a'])) and
                t['final_relative_error']<.001)
            rows.append(row)
    OUT.mkdir(exist_ok=True)
    (OUT/'summary.json').write_text(json.dumps(rows,indent=2))
    keys=list(dict.fromkeys(k for r in rows for k in r))
    write_csv(OUT/'summary.csv',[{k:r.get(k) for k in keys} for r in rows])
    passing=[r for r in rows if r.get('physics_checks_pass') and r['timing_convergence_pass'] and r['responsivity_a_w']>=.28 and r['settling_ns']<=10]
    choice=min(passing,key=lambda r:(r['reverse_bias_v'],-r['responsivity_a_w'])) if passing else None
    selection={'criterion':PROTOCOL['selection'],'chosen':choice,'eligible_case_ids':[r['case_id'] for r in passing],
               'qualification':'Provisional until mesh and capacitance crosschecks; not a fabricated or calibrated device.'}
    if choice and (OUT/'mesh_check.json').exists() and (OUT/'capacitance_crosscheck.json').exists():
        mc=read(OUT/'mesh_check.json');cc=read(OUT/'capacitance_crosscheck.json')
        if mc['case_id']==choice['case_id'] and mc['pass'] and cc['case_id']==choice['case_id']:
            selection['qualification']='Selected reference passed the predeclared DC/capacitance mesh check; independent carrier-inventory capacitance crosscheck completed. See both audit files. Dynamic spatial-mesh convergence, fabrication and calibrated device accuracy are not claimed.'
    (OUT/'selection.json').write_text(json.dumps(selection,indent=2))
    print('SUMMARY',json.dumps(selection),flush=True)
    return rows,choice

def audit():
    rows,choice=summarise()
    if not choice:return
    d,v=choice['thickness_um'],choice['reverse_bias_v']
    r=execute(d,v,ident(d,v)+'_meshhalf',False,1,.5)
    base=read(OUT/choice['source_folder']/'summary.json')
    check={'case_id':choice['case_id'],'mesh_scale':.5,'relative_differences_pct':{}}
    for key in ['responsivity_a_w','capacitance_light_100khz_f']:
        check['relative_differences_pct'][key]=100*(r['metrics'][key]-base['metrics'][key])/base['metrics'][key]
    check['pass']=all(abs(x)<2 for x in check['relative_differences_pct'].values())
    # Independent carrier inventory derivative, not local contact-field charge.
    with (OUT/'capacitance_crosscheck.log').open('w') as log,contextlib.redirect_stdout(log):
        import devsim as ds
        p=PIN(Config(thickness_um=d));p.bias(v);p.light(1)
        measured=p.ac(1e5)['capacitance_f']
        inventories=[]
        for voltage in [v-.001,v+.001]:
            p.bias(voltage)
            fields=p.fields()
            n=np.array(ds.get_node_model_values(device=p.d,region=p.r,name='Electrons'))
            h=np.array(ds.get_node_model_values(device=p.d,region=p.r,name='Holes'))
            volumes=np.array(ds.get_node_model_values(device=p.d,region=p.r,name='NodeVolume'))
            inventories.append((float(np.sum(n*volumes)),float(np.sum(h*volumes))))
        electron=-Q*p.c.area_cm2*(inventories[1][0]-inventories[0][0])/.002
        hole=-Q*p.c.area_cm2*(inventories[1][1]-inventories[0][1])/.002
    cap={'case_id':choice['case_id'],'ac_light_100khz_f':measured,
         'electron_inventory_capacitance_f':electron,'hole_inventory_capacitance_f':hole,
         'electron_relative_difference_pct':100*(electron-measured)/measured,
         'hole_relative_difference_pct':100*(hole-measured)/measured,
         'note':'Each carrier-inventory derivative is compared separately; do not sum. Finite-frequency/transport effects may cause discrepancy, which must be reported.'}
    (OUT/'mesh_check.json').write_text(json.dumps(check,indent=2))
    (OUT/'capacitance_crosscheck.json').write_text(json.dumps(cap,indent=2))
    print('AUDIT',json.dumps({'mesh':check,'capacitance':cap}),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=['initial','refine','summary','audit','all'],default='all');a=ap.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    if a.mode in ['initial','all']:initial()
    if a.mode in ['refine','all']:refine()
    if a.mode in ['summary','all']:summarise()
    if a.mode in ['audit','all']:audit()
