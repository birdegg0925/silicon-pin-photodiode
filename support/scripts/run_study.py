"""Device-simulation helpers; the repository README distinguishes the older standalone grid."""
from pathlib import Path
from dataclasses import asdict
import argparse, contextlib, hashlib, json, math, traceback
import numpy as np
from pin_device import PIN, Config, Q, H, C0, EPS, write_csv

ROOT=Path(__file__).resolve().parents[1]


def case(c,bias,out,transient=False,time_scale=1):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    with (out/'solver.log').open('w') as log,contextlib.redirect_stdout(log):
        p=PIN(c); p.bias(bias)
        dark=p.current(); dark_check=p.terminal_checks()
        dark_ac=[p.ac(f) for f in (1e3,1e5,1e6,1e7)]
        # Local electrode-field charge probe, NOT total diode capacitance:
        # screened ohmic contacts leave this derivative near zero. Preserved
        # as a rejected initial validation idea. Correct independent carrier-
        # inventory capacitance crosscheck is in run_pin_revision.py (audit).
        p.bias(bias-1e-3); qm=p.charge()
        p.bias(bias+1e-3); qp=p.charge()
        p.bias(bias)
        cfd=(qp-qm)/.002
        p.light(1)
        current=p.current(); ac=[p.ac(f) for f in (1e3,1e5,1e6,1e7)]
        fields=p.fields(); check=p.terminal_checks()
        write_csv(out/'fields.csv',[dict(zip(fields,vals)) for vals in zip(*fields.values())])
        write_csv(out/'admittance.csv',[dict(illumination_scale=s,**r) for s,vals in [(0,dark_ac),(1,ac)] for r in vals])
        photon_limit=Q*c.wavelength_nm*1e-9/(H*C0)*(1-c.reflection)*(1-math.exp(-c.alpha_cm*(c.thickness_um+2*c.contact_um)*1e-4))
        row=dict(thickness_um=c.thickness_um,reverse_bias_v=bias,mesh_scale=c.mesh_scale,
                 lifetime_s=c.lifetime_s,incident_power_w=c.optical_power_w,
                 dark_current_a=dark,photocurrent_a=current-dark,responsivity_a_w=(current-dark)/c.optical_power_w,
                 capacitance_dark_100khz_f=dark_ac[1]['capacitance_f'],capacitance_light_100khz_f=ac[1]['capacitance_f'],
                 dark_charge_derivative_f=cfd,geometric_absorber_capacitance_f=EPS*c.area_cm2/(c.thickness_um*1e-4),
                 absorbed_photon_responsivity_ceiling_a_w=photon_limit,
                 collected_fraction_of_absorbed_photons=(current-dark)/(c.optical_power_w*photon_limit),
                 max_field_v_cm=check['max_field_v_cm'],current_imbalance_a=check['current_imbalance_a'],
                 nodes=len(fields['x']))
        detail=dict(config=asdict(c),metrics=row,dc_dark_checks=dark_check,dc_light_checks=check,
                    source_sha256=hashlib.sha256((ROOT/'scripts/pin_device.py').read_bytes()).hexdigest(),
                    engine='DEVSIM 2.11.0 / UMFPACK',status='SIMULATED except analytic checks explicitly named')
        if transient:
            step,metrics=p.optical_step(time_scale=time_scale)
            write_csv(out/'optical_step.csv',step)
            (out/'transient_summary.json').write_text(json.dumps(metrics,indent=2))
            # TIA replay uses dark-subtracted EXTERNAL current at fixed bias. Cd
            # describes electrical terminal response, not a second optical pulse.
            write_csv(out/'photocurrent_replay.csv',[dict(time_s=r['time_s'],photocurrent_a=r['total_current_a']-dark) for r in step])
            row.update({k:metrics[k] for k in ('rise_10_90_s','settling_2pct_s','final_relative_error','source_displacement_check_max_a')})
            detail['transient']=metrics
        (out/'summary.json').write_text(json.dumps(detail,indent=2))
    print('DONE',out.name,json.dumps(row),flush=True)
    return row


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--mode',choices=['grid','transient'],default='grid')
    args=ap.parse_args(); rows=[]
    cases=[(d,v) for d in (5,10,20,40) for v in (.5,1,2,5,10)] if args.mode=='grid' else [(5,5),(10,5),(20,5),(40,5),(40,.5),(40,10)]
    for d,v in cases:
        name=f'd{d:g}_v{v:g}'
        out=ROOT/'runs'/args.mode/name
        try: rows.append(case(Config(thickness_um=d),v,out,transient=args.mode=='transient'))
        except Exception as e:
            out.mkdir(parents=True,exist_ok=True)
            (out/'failure.txt').write_text(traceback.format_exc())
            print('FAILED',name,repr(e),flush=True)
    if rows: write_csv(ROOT/'results'/f'{args.mode}_summary.csv',rows)


if __name__=='__main__': main()
