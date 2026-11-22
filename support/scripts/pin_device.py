"""One-dimensional planar silicon PIN: Poisson + drift diffusion in DEVSIM.

All internal lengths are cm, densities cm^-3 and current densities A/cm^2.
Terminal area is applied ONCE at export. This is a device model, not a process
simulation or a calibrated commercial photodiode. See references/physical_parameters.md.
"""
from dataclasses import dataclass, asdict
from pathlib import Path
import argparse
import csv
import json
import math
import numpy as np
import devsim as ds
from devsim.python_packages import simple_physics as sp
from devsim.python_packages.model_create import (
    CreateSolution, CreateNodeModel, CreateNodeModelDerivative,
    CreateEdgeModel, CreateEdgeModelDerivatives,
)

Q = 1.602176634e-19
K = 1.380649e-23
H = 6.62607015e-34
C0 = 299792458.0
EPS = 11.7 * 8.8541878128e-14


@dataclass
class Config:
    thickness_um: float = 20.0
    area_cm2: float = 0.0025
    bulk_doping_cm3: float = 1e13
    contact_doping_cm3: float = 1e17
    contact_um: float = 0.3
    junction_smoothing_um: float = 0.02
    mesh_scale: float = 1.0
    lifetime_s: float = 1e-6
    wavelength_nm: float = 850.0
    alpha_cm: float = 534.0
    refractive_index: float = 3.641
    optical_power_w: float = 1e-6
    temperature_k: float = 300.0
    intrinsic_cm3: float = 9.65e9

    @property
    def reflection(self):
        extinction = self.alpha_cm * self.wavelength_nm * 1e-7 / (4*math.pi)
        return ((self.refractive_index-1)**2+extinction**2)/((self.refractive_index+1)**2+extinction**2)


class PIN:
    def __init__(self, config: Config):
        self.c = c = config
        if c.temperature_k != 300.0:
            raise ValueError('This optical/transport parameter set is documented only at 300 K.')
        self.d = 'pin'
        self.r = 'silicon'
        if self.d in ds.get_device_list():
            ds.delete_device(device=self.d)
        ds.reset_devsim()
        from devsim.umfpack.umfshim import local_solver_callback
        ds.set_parameter(name='direct_solver',value='custom')
        ds.set_parameter(name='solver_callback',value=local_solver_callback)
        # Extended precision reduces cancellation when extracting tiny dark currents.
        for key in ('extended_solver', 'extended_model', 'extended_equation'):
            ds.set_parameter(name=key, value=True)
        ds.create_1d_mesh(mesh='mesh')
        a, b = c.contact_um*1e-4, (c.contact_um+c.thickness_um)*1e-4
        total = b+a
        fine = 0.005e-4*c.mesh_scale
        coarse = 0.10e-4*c.mesh_scale
        anchors = [(0, fine, 'front'), (a-0.1e-4,fine,''),
                   (a,fine,''),(a+0.1e-4,fine,''),
                   ((a+b)/2,coarse,''), (b-0.1e-4,fine,''),
                   (b,fine,''),(b+0.1e-4,fine,''),(total,fine,'back')]
        for pos, spacing, tag in anchors:
            ds.add_1d_mesh_line(mesh='mesh', pos=pos, ps=spacing, tag=tag)
        ds.add_1d_contact(mesh='mesh', name='anode', tag='front', material='metal')
        ds.add_1d_contact(mesh='mesh', name='cathode', tag='back', material='metal')
        ds.add_1d_region(mesh='mesh', material='Si', region=self.r, tag1='front', tag2='back')
        ds.finalize_mesh(mesh='mesh')
        ds.create_device(mesh='mesh', device=self.d)
        sp.SetSiliconParameters(self.d,self.r,c.temperature_k)
        pars = dict(Permittivity=EPS, ElectronCharge=Q, n_i=c.intrinsic_cm3,
                    V_t=K*c.temperature_k/Q,kT=K*c.temperature_k,
                    n1=c.intrinsic_cm3,p1=c.intrinsic_cm3,
                    taun=c.lifetime_s,taup=c.lifetime_s, optical_scale=0.0)
        for name,value in pars.items():
            ds.set_parameter(device=self.d,region=self.r,name=name,value=value)
        width=c.junction_smoothing_um*1e-4
        doping=f'{c.bulk_doping_cm3} - {c.contact_doping_cm3}*0.5*erfc((x-{a})/{width}) + {c.contact_doping_cm3}*0.5*erfc(({b}-x)/{width})'
        CreateNodeModel(self.d,self.r,'NetDoping',doping)
        sp.CreateSiliconPotentialOnly(self.d,self.r)
        ds.set_parameter(device=self.d,name='anode_bias',value=0.0)
        ds.circuit_element(name='Vbias',n1='cathode_bias',n2='0',value=0.0,acreal=1.0,acimag=0.0)
        sp.CreateSiliconPotentialOnlyContact(self.d,self.r,'anode')
        sp.CreateSiliconPotentialOnlyContact(self.d,self.r,'cathode',is_circuit=True)
        ds.solve(type='dc',absolute_error=1.0,relative_error=1e-10,maximum_iterations=100)
        for solution,initial in [('Electrons','IntrinsicElectrons'),('Holes','IntrinsicHoles')]:
            CreateSolution(self.d,self.r,solution)
            ds.set_node_values(device=self.d,region=self.r,name=solution,init_from=initial)
        # High-resistivity <100> bulk, Scharf & Klanner 2015 Table 4 (300 K).
        # This empirical field law is also used in the thin contacts; no doping-
        # dependent mobility/BGN/contact resistance are claimed.
        fields = {'AbsField':'(ElectricField^2+1)^0.5',
                  'HighHoleField':'max(AbsField-2640,0)',
                  'ElectronMobility':'1530/(1+1530*AbsField/1.03e7)',
                  'HoleMobility':'1/(1/464+9.57e-8*HighHoleField-3.31e-13*HighHoleField^2)'}
        for name,expr in fields.items():
            CreateEdgeModel(self.d,self.r,name,expr)
            CreateEdgeModelDerivatives(self.d,self.r,name,expr,'Potential')
        sp.CreateSiliconDriftDiffusion(self.d,self.r,'ElectronMobility','HoleMobility')
        sp.CreateSiliconDriftDiffusionAtContact(self.d,self.r,'anode')
        sp.CreateSiliconDriftDiffusionAtContact(self.d,self.r,'cathode',is_circuit=True)
        photon_energy=H*C0/(c.wavelength_nm*1e-9)
        g0=c.alpha_cm*(1-c.reflection)*c.optical_power_w/(c.area_cm2*photon_energy)
        CreateNodeModel(self.d,self.r,'OpticalGeneration',f'optical_scale*{g0}*exp(-{c.alpha_cm}*x)')
        for name,expr in [('ElectronGeneration','ElectronCharge*(OpticalGeneration-USRH)'),
                          ('HoleGeneration','ElectronCharge*(USRH-OpticalGeneration)')]:
            CreateNodeModel(self.d,self.r,name,expr)
            for variable in ('Electrons','Holes'):
                CreateNodeModelDerivative(self.d,self.r,name,expr,variable)
        self.solve()
        self.bias_v=0.0
        self.light_scale=0.0

    def solve(self, **kwargs):
        ds.solve(type='dc',absolute_error=1e-8,relative_error=1e-9,maximum_iterations=80,**kwargs)

    def bias(self, voltage):
        count=max(1,int(math.ceil(abs(voltage-self.bias_v)/0.2)))
        for val in np.linspace(self.bias_v,voltage,count+1)[1:]:
            ds.circuit_alter(name='Vbias',value=float(val))
            self.solve()
        self.bias_v=voltage

    def light(self, scale):
        ds.set_parameter(device=self.d,region=self.r,name='optical_scale',value=float(scale))
        self.solve()
        self.light_scale=scale

    def current(self,contact='cathode'):
        return self.c.area_cm2*sum(ds.get_contact_current(device=self.d,contact=contact,equation=e)
            for e in ('ElectronContinuityEquation','HoleContinuityEquation'))

    def charge(self):
        return self.c.area_cm2*ds.get_contact_charge(device=self.d,contact='cathode',equation='PotentialEquation')

    def ac(self,frequency):
        ds.solve(type='ac',frequency=frequency)
        # Current through the exciting voltage source has opposite terminal sign.
        y=-self.c.area_cm2*complex(ds.get_circuit_node_value(node='Vbias.I',solution='ssac_real'),
                                 ds.get_circuit_node_value(node='Vbias.I',solution='ssac_imag'))
        return {'frequency_hz':frequency,'admittance_real_s':y.real,'admittance_imag_s':y.imag,
                'capacitance_f':y.imag/(2*math.pi*frequency)}

    def fields(self):
        names=['x','Potential','NetDoping','Electrons','Holes','OpticalGeneration','USRH']
        return {name:np.array(ds.get_node_model_values(device=self.d,region=self.r,name=name)) for name in names}

    def terminal_checks(self):
        ef=np.array(ds.get_edge_model_values(device=self.d,region=self.r,name='ElectricField'))
        jn=np.array(ds.get_edge_model_values(device=self.d,region=self.r,name='ElectronCurrent'))
        jp=np.array(ds.get_edge_model_values(device=self.d,region=self.r,name='HoleCurrent'))
        return {'anode_current_a':self.current('anode'),'cathode_current_a':self.current(),
                'current_imbalance_a':self.current()+self.current('anode'),
                'max_field_v_cm':float(np.max(np.abs(ef))),
                'current_density_spread_a_cm2':float(np.ptp(jn+jp))}

    def optical_step(self, initial=0.9, final=1.0, time_scale=1.0, time_end_s=2e-6):
        """Small optical step at fixed voltage; BDF1 external total current.

        The DC optical operating point is 1 uW. A 0.9->1.0 uW step is used to
        avoid treating a large carrier-injection event as a linear impulse.
        Nonuniform steps resolve sub-ns drift and long neutral-region tails.
        """
        self.light(initial)
        i0=self.current(); self.light(final); i1=self.current()
        self.light(initial)
        ds.solve(type='transient_dc',absolute_error=1e-8,relative_error=1e-9,maximum_iterations=80)
        times=np.unique(np.r_[0, np.arange(0.025e-9*time_scale,5e-9,0.025e-9*time_scale),
                      np.arange(5e-9,30e-9,0.125e-9*time_scale),
                      np.arange(30e-9,200e-9,1e-9*time_scale),
                      np.arange(200e-9,2e-6,10e-9*time_scale),2e-6])
        times=np.unique(np.r_[times[times<time_end_s],time_end_s])
        ds.set_parameter(device=self.d,region=self.r,name='optical_scale',value=final)
        prev_charge=self.charge()
        rows=[dict(time_s=0,total_current_a=i0,conduction_current_a=i0,
                   displacement_current_a=0,source_current_check_a=i0,
                   optical_power_w=initial*self.c.optical_power_w)]
        for t0,t1 in zip(times[:-1],times[1:]):
            dt=float(t1-t0)
            # DEVSIM charge_error tests time-step charge change, not a DC Newton
            # residual. A discontinuous optical input fails a 1e-6 threshold.
            # Accept up to 50% per step; independently halve the time grid and
            # check total-current/charge consistency before accepting metrics.
            ds.solve(type='transient_bdf1',tdelta=dt,charge_error=0.5,
                     absolute_error=1e-8,relative_error=1e-8,maximum_iterations=80)
            charge=self.charge()
            conduction=self.current()
            displacement=(charge-prev_charge)/dt
            source=-self.c.area_cm2*ds.get_circuit_node_value(node='Vbias.I',solution='dcop')
            rows.append(dict(time_s=float(t1),total_current_a=source,
                             conduction_current_a=conduction,displacement_current_a=displacement,
                             source_current_check_a=conduction+displacement,
                             optical_power_w=final*self.c.optical_power_w))
            prev_charge=charge
        self.light_scale=final
        y=np.array([r['total_current_a'] for r in rows]); norm=(y-i0)/(i1-i0)
        # First threshold crossing, valid only after checking monotonicity/overshoot.
        def cross(level):
            inds=np.where(norm>=level)[0]
            if not len(inds): return None
            i=inds[0]
            return float(np.interp(level,norm[max(0,i-1):i+1],times[max(0,i-1):i+1]))
        t10,t90=cross(.1),cross(.9)
        bad=np.where(np.abs(norm-1)>.02)[0]
        settle=float(times[bad[-1]+1]) if len(bad) and bad[-1]+1<len(times) else None
        metrics=dict(step_initial_scale=initial,step_final_scale=final,time_scale=time_scale,time_end_s=time_end_s,
                     initial_current_a=i0,final_current_a=i1,t10_s=t10,t90_s=t90,
                     rise_10_90_s=(t90-t10 if t90 is not None and t10 is not None else None),
                     settling_2pct_s=settle,final_relative_error=float(abs(norm[-1]-1)),
                     max_normalised=float(norm.max()),
                     source_displacement_check_max_a=float(max(abs(r['total_current_a']-r['source_current_check_a']) for r in rows)))
        return rows,metrics


def write_csv(path,rows):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--thickness',type=float,default=20)
    p.add_argument('--bias',type=float,default=5)
    p.add_argument('--mesh-scale',type=float,default=1)
    p.add_argument('--output',default='generated/device_smoke')
    p.add_argument('--transient',action='store_true')
    a=p.parse_args(); out=Path(a.output); out.mkdir(parents=True,exist_ok=True)
    c=Config(thickness_um=a.thickness,mesh_scale=a.mesh_scale)
    pin=PIN(c); pin.bias(a.bias)
    dark=pin.current(); pin.light(1)
    result=dict(config=asdict(c),reverse_bias_v=a.bias,dark_current_a=dark,
                photocurrent_a=pin.current()-dark,responsivity_a_w=(pin.current()-dark)/c.optical_power_w,
                ac=[pin.ac(f) for f in (1e3,1e5,1e6,1e7)],checks=pin.terminal_checks())
    fields=pin.fields()
    write_csv(out/'fields.csv',[dict(zip(fields,vals)) for vals in zip(*fields.values())])
    (out/'summary.json').write_text(json.dumps(result,indent=2))
    if a.transient:
        rows,metrics=pin.optical_step()
        write_csv(out/'optical_step.csv',rows)
        (out/'transient_summary.json').write_text(json.dumps(metrics,indent=2))
        print('TRANSIENT',json.dumps(metrics))
    print('SUMMARY',json.dumps(result))


if __name__=='__main__': main()
