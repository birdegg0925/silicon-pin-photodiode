"""Check saved scalar evidence without requiring DEVSIM or LTspice."""
from pathlib import Path
import json,csv,ast,math,hashlib
ROOT=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((ROOT/'pin/receiver_summary.csv').open()))
by={r['case_id']:r for r in rows}
assert math.isclose(float(by['d35.6_v1.55']['device_rise_ns']),78.78287789118956,rel_tol=1e-12)
assert math.isclose(float(by['d35.6_v7.86']['device_rise_ns']),1.699899756277028,rel_tol=1e-12)
assert math.isclose(float(by['d35.6_v1.55']['receiver_optical_bw_mhz']),3.728334289642624,rel_tol=1e-12)
assert math.isclose(float(by['d35.6_v7.86']['receiver_optical_bw_mhz']),5.89934255698705,rel_tol=1e-12)
assert json.loads((ROOT/'pin/mesh_check.json').read_text())['pass']
json_rows={r['case_id']:r for r in json.loads((ROOT/'pin/receiver_summary.json').read_text())}
device_rows={r['case_id']:r for r in json.loads((ROOT/'pin/summary.json').read_text())}
for r in rows:
    cid=r['case_id']
    for key in ['device_rise_ns','device_settling_ns','tia_electrical_bw_mhz','receiver_optical_bw_mhz']:
        assert math.isclose(float(r[key]),json_rows[cid][key],rel_tol=1e-12),(cid,key)
    detail=json.loads((ROOT/'pin'/device_rows[cid]['source_folder']/'summary.json').read_text())
    assert math.isclose(detail['transient']['rise_10_90_s']*1e9,float(r['device_rise_ns']),rel_tol=1e-12),cid
    transient=json.loads((ROOT/'pin'/device_rows[cid]['source_folder']/'transient_summary.json').read_text())
    assert transient==detail['transient'],cid
models=[ROOT/'tia/original_tia.lib',ROOT/'support/models/original_tia.lib']
assert models[0].read_bytes()==models[1].read_bytes()
expected=json.loads((ROOT/'tia/selection.json').read_text())['original_model_sha256']
assert hashlib.sha256(models[0].read_bytes()).hexdigest()==expected
for p in ROOT.rglob('*.py'):
    if '.venv' not in p.parts:ast.parse(p.read_text(),filename=str(p))
for p in ['support/scripts/pin_device.py','support/scripts/run_study.py','support/scripts/tia_interface.py','tia/original_tia.lib','tia/selection.json','protocol.json']:
    assert (ROOT/p).is_file(),p
print('PASS: CSV/JSON/transient cross-consistency, model hashes, headline values, mesh-check record, Python syntax and bundled files')
