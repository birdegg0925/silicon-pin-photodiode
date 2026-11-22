#!/usr/bin/env python3
"""One-way PIN reduced-model interface to the seven-MOS LTspice TIA.

Input JSON::

  {"settings": {"temperature_c": 27, "feedback_capacitance_f": 5e-13},
   "cases": [{"case_id": "example", "capacitance_f": 5e-12,
              "photocurrent_a": 2e-6, "responsivity_a_w": 0.5,
              "reverse_bias_v": 2.0, "dark_current_a": 0.0}]}

Optional per-case optical_pulse_csv has columns time_s,photocurrent_a. It is
the TOTAL photocurrent (not an increment), with strictly increasing times
starting at zero; dark_current_a is added separately. This is fixed-C,
one-way waveform replay around an extracted device operating condition,
NOT nonlinear device/circuit cosimulation. A frequency-independent
responsivity is only reported as a DC conversion, not as an optical BW.

Use --prepare-only to inspect decks without simulation. Default analyses
are op/ac/noise. All numerical results are extracted from actual LTspice
ASCII raw files. Detector shot noise is a separately labelled calculation;
ideal photocurrent and dark-current sources contribute no SPICE noise.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
SOURCE_MODEL = ROOT / "models/original_tia.lib"
SHARED_LOCK = ROOT.parent / ".ltspice_serial.lock"
RUNNER = os.environ.get("LTSPICE_RUNNER", "/Applications/LTspice.app/Contents/SharedSupport/ltspice/LTspice/run_ltspice")
ELECTRON_CHARGE = 1.602176634e-19
DEFAULTS = dict(supply_v=1.8, reference_v=0.9, feedback_resistance_ohm=1e5,
                feedback_capacitance_f=0.5e-12, output_load_f=1e-12,
                internal_compensation_f=1e-12, temperature_c=27.0,
                ac_start_hz=10.0, ac_stop_hz=1e9, ac_points_per_decade=160,
                noise_low_hz=1e3, noise_high_hz=10e6,
                transient_max_step_s=1e-9, timeout_s=180)


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def clean(value):
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, np.generic):
        value = value.item()
    return None if isinstance(value, float) and not math.isfinite(value) else value


def write_json(path, data):
    path.write_text(json.dumps(clean(data), indent=2) + "\n")


def validate_case(case):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", case["case_id"]):
        raise ValueError("case_id must be a safe alphanumeric filename")
    for key in ["capacitance_f", "photocurrent_a", "responsivity_a_w", "reverse_bias_v"]:
        value = float(case[key])
        if not math.isfinite(value) or value < 0 or (key == "capacitance_f" and value == 0):
            raise ValueError(f"Invalid {key}: {value}")
    dark = float(case.get("dark_current_a", 0.0))
    if not math.isfinite(dark) or dark < 0:
        raise ValueError("dark_current_a must be finite/nonnegative")


def export_grid_cases(grid_root, destination):
    """Derive interface inputs from completed device runs."""
    cases, temperatures = [], set()
    for path in sorted(Path(grid_root).glob("d*/summary.json")):
        source = json.loads(path.read_text())
        metrics, config = source["metrics"], source["config"]
        temperatures.add(float(config["temperature_k"]) - 273.15)
        case = dict(case_id=path.parent.name,
                    capacitance_f=metrics["capacitance_light_100khz_f"],
                    photocurrent_a=metrics["photocurrent_a"],
                    responsivity_a_w=metrics["responsivity_a_w"],
                    reverse_bias_v=metrics["reverse_bias_v"],
                    dark_current_a=metrics["dark_current_a"],
                    absorber_thickness_um=metrics["thickness_um"],
                    incident_optical_power_w=metrics["incident_power_w"],
                    capacitance_extraction_frequency_hz=1e5,
                    device_summary_path=str(path.resolve()),
                    device_summary_sha256=sha256_bytes(path.read_bytes()),
                    dark_current_scope="Uncalibrated assumed-lifetime SRH model prediction; not a measured or commercial-device dark-current specification",
                    reduction_scope="Fixed illuminated differential capacitance at100kHz plus extracted DC photocurrent and dark current; electrical conductance/frequency dependence neglected; validity must be checked against device admittance")
        validate_case(case)
        cases.append(case)
    if not cases or len(temperatures) != 1:
        raise ValueError("Expected completed cases at one common temperature")
    cases.sort(key=lambda c: (c["absorber_thickness_um"], c["reverse_bias_v"]))
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    write_json(destination, dict(settings={"temperature_c": temperatures.pop()}, cases=cases,
                                 provenance="Interface parameters extracted without rounding from completed DEVSIM summaries"))
    return cases


def export_step_cases(destination):
    """Extract the older six-case layout described in the repository README."""
    grid = json.loads((ROOT / "results/tia_grid_cases.json").read_text())
    by_id = {case["case_id"]: case for case in grid["cases"]}
    cases = []
    for thickness, bias in [(5, 5), (10, 5), (20, 5), (40, 5), (40, .5), (40, 10)]:
        case_id = f"d{thickness:g}_v{bias:g}"
        candidates = [ROOT / "runs/transient" / case_id] + list((ROOT / "runs/validation").glob("*"))
        available = []
        for folder in candidates:
            summary_path = folder / "summary.json"
            if not summary_path.exists() or not (folder / "photocurrent_replay.csv").exists():
                continue
            source = json.loads(summary_path.read_text())
            config, metrics = source["config"], source["metrics"]
            if (config["thickness_um"] == thickness and metrics["reverse_bias_v"] == bias
                    and config["mesh_scale"] == 1 and "transient" in source):
                with (folder / "photocurrent_replay.csv").open(newline="") as handle:
                    pulse_rows = list(csv.DictReader(handle))
                # Short sub-ns timing-refinement probes do not cover the TIA
                # settling interval; require the full2µs device waveform.
                if pulse_rows and float(pulse_rows[-1]["time_s"]) >= 2e-6*(1-1e-9):
                    available.append((source["transient"]["time_scale"], folder))
        if not available:
            raise ValueError(f"No completed device pulse for {case_id}")
        _, selected = min(available)
        pulse_path = selected / "photocurrent_replay.csv"
        with pulse_path.open(newline="") as handle:
            initial = float(next(csv.DictReader(handle))["photocurrent_a"])
        source = json.loads((selected / "summary.json").read_text())
        case = {**by_id[case_id], "photocurrent_a": initial,
                "optical_pulse_csv": str(pulse_path.resolve()),
                "optical_pulse_sha256": sha256_bytes(pulse_path.read_bytes()),
                "transient_source_summary_sha256": sha256_bytes((selected / "summary.json").read_bytes()),
                "device_transient_time_scale": source["transient"]["time_scale"],
                "initial_incident_optical_power_w": .9 * by_id[case_id]["incident_optical_power_w"],
                "final_incident_optical_power_w": by_id[case_id]["incident_optical_power_w"],
                "pulse_model_scope": "One-way DEVSIM total-photocurrent replay excluding dark (added separately); capacitance held at1µW illuminated100kHz extraction across0.9–1µW optical step; final-DC responsivity not an assumed instantaneous relation"}
        validate_case(case)
        read_pulse(pulse_path, initial)
        cases.append(case)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    write_json(destination, dict(settings={**grid["settings"], "transient_max_step_s": .5e-9}, cases=cases))
    return cases


def read_pulse(path, dc_current):
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    times = np.array([float(r["time_s"]) for r in rows])
    currents = np.array([float(r["photocurrent_a"]) for r in rows])
    if len(times) < 2 or times[0] != 0 or np.any(np.diff(times) <= 0):
        raise ValueError("Pulse must start at zero and contain strictly increasing times")
    if not np.all(np.isfinite(times)) or not np.all(np.isfinite(currents)) or np.any(currents < 0):
        raise ValueError("Pulse times/current must be finite and photocurrent nonnegative")
    if not np.isclose(currents[0], dc_current, rtol=1e-5, atol=1e-18):
        raise ValueError("Pulse initial photocurrent must equal case.photocurrent_a for the DC operating point")
    return times, currents


def deck(case, settings, analysis, model_text, pulse=None):
    """Cathode=DET; anode=ANODE. Reverse current goes DET -> ANODE.

    ANODE is held at VREF - Vreverse by an ideal noiseless EXTERNAL rail.
    Bias > VREF consequently requires a negative external rail. Its
    implementation/noise/power conversion efficiency are not modelled.
    Extracted capacitance replaces, rather than parallels, the old 5 pF.
    """
    p = {**DEFAULTS, **settings}
    c = case
    source = f"DC {c['photocurrent_a']:.15g} AC 1"
    if analysis == "tran":
        if pulse is None:
            raise ValueError("Transient requested without optical_pulse_csv")
        source = "PWL FILE=\"photocurrent.pwl\" AC 1"
    text = f"""* Extracted PIN model -> original seven-MOS TIA; no nonlinear TCAD cosimulation.
* CASE {c['case_id']}; Cd replaces original5p; no extra original detector capacitance.
* Positive photogenerated conventional current: n/cathode DET -> p/anode ANODE.
* VBIAS is external ideal quiet bias; its implementation is NOT in1.8Vcore power.
* Photocurrent/dark current sources are noiseless; no implicit shot noise.
.options plotwinsize=0 numdgt=15 reltol=1e-5 abstol=1e-13 vntol=1e-8 chgtol=1e-20 trtol=1 method=trap
.temp {p['temperature_c']:.15g}
VDDTIA VDD 0 {p['supply_v']:.15g}
VREFTIA REF 0 {p['reference_v']:.15g}
VBIAS ANODE 0 {p['reference_v'] - c['reverse_bias_v']:.15g}
IPHOTO DET ANODE {source}
IDARK DET ANODE {c.get('dark_current_a', 0):.15g}
CDET DET ANODE {c['capacitance_f']:.15g}
RF OUT DET {p['feedback_resistance_ohm']:.15g}
CF OUT DET {p['feedback_capacitance_f']:.15g}
XTIA REF DET OUT VDD 0 CMOS_OPAMP CCOMP={p['internal_compensation_f']:.15g} BIAS_SCALE=1
CL OUT 0 {p['output_load_f']:.15g}
* Verbatim original educational amplifier model follows.
{model_text}
"""
    grid = f"dec {int(p['ac_points_per_decade'])} {p['ac_start_hz']:.15g} {p['ac_stop_hz']:.15g}"
    if analysis == "op":
        text += ".options logopinfo\n.op\n"
    elif analysis == "ac":
        text += f".ac {grid}\n.save V(OUT) V(DET) V(ANODE) I(IPHOTO) I(VDDTIA) I(VBIAS)\n"
    elif analysis == "noise":
        text += f".noise V(OUT) IPHOTO {grid}\n"
    elif analysis == "tran":
        stop = float(pulse[0][-1])
        text += f".tran 0 {stop:.15g} 0 {p['transient_max_step_s']:.15g}\n.save V(OUT) V(DET) V(ANODE) I(IPHOTO) I(VDDTIA) I(VBIAS)\n"
    else:
        raise ValueError(analysis)
    return text + ".end\n"


def parse_raw(path):
    """Parse LTspice ASCII raw files."""
    data = Path(path).read_bytes()
    text = data.decode("utf-16") if data[:2] in [b"\xff\xfe", b"\xfe\xff"] else data.decode("latin1")
    lines = text.splitlines()
    a, b = lines.index("Variables:"), lines.index("Values:")
    names = [line.split()[1].lower() for line in lines[a + 1:b] if line.strip()]
    def number(value):
        pieces = value.split(",")
        return complex(*map(float, pieces)) if len(pieces) == 2 else float(value)
    rows = []
    for i in range(b + 1, len(lines), len(names)):
        block = lines[i:i + len(names)]
        if len(block) == len(names):
            rows.append([number(block[0].split()[-1])] + [number(v.strip()) for v in block[1:]])
    array = np.asarray(rows)
    if not len(array):
        raise ValueError(f"Empty raw file {path}")
    return {name: array[:, i] for i, name in enumerate(names)}


def save_raw_csv(path, record):
    fields, arrays = [], []
    for key, values in record.items():
        if np.iscomplexobj(values) and np.any(values.imag != 0):
            fields += [key + "_real", key + "_imag"]
            arrays += [values.real, values.imag]
        else:
            fields.append(key)
            arrays.append(values.real)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        writer.writerows(zip(*arrays))


def log_interp(f, y, frequency):
    if not f[0] <= frequency <= f[-1]:
        raise ValueError("Interpolation requested outside sampled frequency range")
    return float(np.interp(np.log10(frequency), np.log10(f), y))


def integrate_asd(f, asd, lo, hi):
    selected = (f > lo) & (f < hi)
    ff = np.r_[lo, f[selected], hi]
    yy = np.r_[log_interp(f, asd, lo), asd[selected], log_interp(f, asd, hi)]
    return float(np.sqrt(np.trapezoid(yy * yy, ff)))


def summarise(record, case, settings, analysis, folder):
    p = {**DEFAULTS, **settings}
    base = dict(case_id=case["case_id"], analysis=analysis,
                provenance="CALCULATED metrics from SIMULATED LTspice outputs", conditions={**p, **case})
    if analysis == "op":
        val = lambda key: float(record[key][0].real)
        base.update(output_dc_v=val("v(out)"), detector_cathode_dc_v=val("v(det)"),
                    detector_anode_dc_v=val("v(anode)"),
                    actual_reverse_bias_v=val("v(det)")-val("v(anode)"),
                    delivered_core_current_a=-val("i(vddtia)"),
                    core_supply_power_w=-p["supply_v"]*val("i(vddtia)"),
                    bias_rail_signed_delivered_power_w=-(p["reference_v"]-case["reverse_bias_v"])*val("i(vbias)"),
                    power_scope="Seven-MOS core only; ideal bias/reference generators and external detector-bias implementation excluded")
    elif analysis == "ac":
        f, z = record["frequency"].real, record["v(out)"]
        gain = abs(z)
        reference = gain[0]
        crossing = np.flatnonzero((gain[:-1] >= reference / np.sqrt(2)) & (gain[1:] < reference / np.sqrt(2)))
        bw = None
        if len(crossing):
            i = int(crossing[0])
            fraction = (reference/np.sqrt(2)-gain[i])/(gain[i+1]-gain[i])
            bw = float(10**(np.log10(f[i])+fraction*np.log10(f[i+1]/f[i])))
        base.update(low_frequency_transimpedance_ohm=float(reference),
                    transimpedance_1mhz_ohm=log_interp(f, gain, 1e6), upper_3db_hz=bw,
                    peaking_db=float(20*np.log10(gain.max()/reference)),
                    optical_dc_gain_v_w=float(reference*case["responsivity_a_w"]),
                    optical_gain_scope="DC only; device optical-current frequency response must be combined separately")
        with (folder / "transimpedance.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["frequency_hz", "z_real_ohm", "z_imag_ohm", "z_magnitude_ohm", "z_phase_rad"])
            writer.writerows(zip(f, z.real, z.imag, gain, np.angle(z)))
    elif analysis == "noise":
        f = record["frequency"].real
        output = record["v(onoise)"].real
        input_noise = record["i(inoise)"].real
        z = abs(record["gain"])
        shot = math.sqrt(2*ELECTRON_CHARGE*(case["photocurrent_a"]+case.get("dark_current_a", 0)))
        total = np.sqrt(output**2 + (shot*z)**2)
        lo, hi = p["noise_low_hz"], p["noise_high_hz"]
        base.update(output_circuit_noise_100khz_v_sqrt_hz=log_interp(f, output, 1e5),
                    input_circuit_noise_100khz_a_sqrt_hz=log_interp(f, input_noise, 1e5),
                    calculated_detector_shot_asd_a_sqrt_hz=shot,
                    circuit_output_rms_v=integrate_asd(f, output, lo, hi),
                    circuit_plus_assumed_shot_output_rms_v=integrate_asd(f, total, lo, hi),
                    integration_low_hz=lo, integration_high_hz=hi,
                    shot_scope="CALCULATED sqrt(2*q*(Iphoto+Idark)); independent Poisson currents assumed; not SPICE generated; no generation/recombination or bias-rail noise")
        with (folder / "noise.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["frequency_hz", "circuit_output_v_sqrt_hz", "circuit_input_a_sqrt_hz", "calculated_shot_output_v_sqrt_hz", "calculated_quadrature_total_v_sqrt_hz"])
            writer.writerows(zip(f, output, input_noise, shot*z, total))
    elif analysis == "tran":
        t, v = record["time"].real, record["v(out)"].real
        base.update(output_start_v=float(v[0]), output_min_v=float(v.min()), output_max_v=float(v.max()),
                    duration_s=float(t[-1]), max_recorded_step_s=float(np.diff(t).max()),
                    terminal_reverse_bias_min_v=float((record["v(det)"]-record["v(anode)"]).real.min()),
                    terminal_reverse_bias_max_v=float((record["v(det)"]-record["v(anode)"]).real.max()),
                    scope="One-way total photocurrent replay with fixed extracted capacitance; no device feedback cosimulation")
        amplitude = float(v[-1]-v[0])
        if amplitude:
            normalised = (v-v[0])/amplitude
            def crossing(level):
                indices = np.flatnonzero((normalised[:-1] < level) & (normalised[1:] >= level))
                if not len(indices):
                    return None
                i = int(indices[0])
                return float(t[i]+(level-normalised[i])/(normalised[i+1]-normalised[i])*(t[i+1]-t[i]))
            t10, t90 = crossing(.1), crossing(.9)
            bad = np.flatnonzero(abs(normalised-1) > .02)
            settling = float(t[bad[-1]+1]) if len(bad) and bad[-1]+1 < len(t) else (0.0 if not len(bad) else None)
            tail = v[t >= .9*t[-1]]
            base.update(output_final_v=float(v[-1]), output_step_amplitude_v=amplitude,
                        output_t10_s=t10, output_t90_s=t90,
                        output_rise_10_90_s=t90-t10 if t10 is not None and t90 is not None else None,
                        output_settling_2pct_s=settling,
                        final_10pct_duration_output_range_fraction_of_step=float(np.ptp(tail)/abs(amplitude)),
                        output_overshoot_pct=float(max(0, normalised.max()-1)*100),
                        timing_scope="10–90% and±2% relative to first and last recorded output values; final10%ofduration flatness separately checked; finite observation window")
    return base


def run_case(case, settings, analysis, model, output, pulse_path=None, prepare_only=False, force=False):
    folder = output / case["case_id"] / analysis
    pulse = read_pulse(pulse_path, case["photocurrent_a"]) if analysis == "tran" else None
    content = deck(case, settings, analysis, model, pulse)
    fingerprint = sha256_bytes((content + (sha256_bytes(Path(pulse_path).read_bytes()) if pulse else "")).encode())
    meta_path = folder / "run.json"
    if meta_path.exists() and not force and not prepare_only:
        old = json.loads(meta_path.read_text())
        if old.get("success") and old.get("sha256") == fingerprint and (folder / "test.raw").exists():
            result = summarise(parse_raw(folder / "test.raw"), case, settings, analysis, folder)
            write_json(folder / "summary.json", result)
            print("CACHED", case["case_id"], analysis, flush=True)
            return result
    if meta_path.exists():
        archive = folder.with_name(analysis + "_archived_" + str(time.time_ns()))
        folder.rename(archive)
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "test.cir").write_text(content)
    if pulse:
        np.savetxt(folder / "photocurrent.pwl", np.column_stack(pulse), fmt="%.15g")
    info = dict(sha256=fingerprint, model_sha256=sha256_bytes(model.encode()), source_model=str(SOURCE_MODEL),
                case=case, settings=settings, analysis=analysis, success=False)
    if prepare_only:
        info["prepared_only"] = True
        write_json(meta_path, info)
        print("PREPARED", case["case_id"], analysis, flush=True)
        return info
    start = time.time()
    command = [RUNNER, "-b", "-ascii", "test.cir"]
    try:
        with SHARED_LOCK.open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            result = subprocess.run(command, cwd=folder, env={**os.environ, "MVK_CONFIG_LOG_LEVEL": "0", "WINEDEBUG": "-all"},
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=settings["timeout_s"])
        output_text, returncode = result.stdout, result.returncode
    except subprocess.TimeoutExpired as error:
        output_text, returncode = error.stdout or b"", 124
    (folder / "runner.log").write_bytes(output_text)
    log = (folder / "test.log").read_text(errors="replace") if (folder / "test.log").exists() else ""
    info.update(command=command, returncode=returncode, elapsed_s=time.time()-start,
                success=returncode == 0 and (folder / "test.raw").exists() and "Total elapsed time" in log)
    write_json(meta_path, info)
    if not info["success"]:
        raise RuntimeError(f"LTspice failed: {folder}\n{log[-1800:]}\n{output_text[-1000:].decode(errors='replace')}")
    record = parse_raw(folder / "test.raw")
    save_raw_csv(folder / "raw.csv", record)
    summary = summarise(record, case, settings, analysis, folder)
    write_json(folder / "summary.json", summary)
    print("FINISHED", case["case_id"], analysis, round(info["elapsed_s"], 2), flush=True)
    return summary


def export_metrics_csv(output):
    """Combine only successful, saved analyses without starting the simulator."""
    output = Path(output)
    rows = []
    for case_dir in sorted(p for p in output.iterdir() if p.is_dir()):
        combined = {"case_id": case_dir.name}
        completed = []
        for analysis in ["op", "ac", "noise", "tran"]:
            folder = case_dir / analysis
            if not (folder / "run.json").exists() or not (folder / "summary.json").exists():
                continue
            if not json.loads((folder / "run.json").read_text()).get("success"):
                continue
            result = json.loads((folder / "summary.json").read_text())
            for key in ["absorber_thickness_um", "reverse_bias_v", "capacitance_f", "photocurrent_a", "dark_current_a", "responsivity_a_w", "incident_optical_power_w", "temperature_c"]:
                if key in result["conditions"]:
                    combined[key] = result["conditions"][key]
            for key, value in result.items():
                if key not in ["case_id", "analysis", "conditions"]:
                    combined[key] = value
            completed.append(analysis)
        if completed:
            combined["completed_analyses"] = ";".join(completed)
            rows.append(combined)
    if rows:
        fields = list(dict.fromkeys(key for row in rows for key in row))
        with (output / "metrics.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fields)
            writer.writeheader()
            writer.writerows(clean(rows))
    return rows


def audit_reduced_admittance(grid_root, run_root, destination):
    """Exact linear port-loading correction at simulated admittance frequencies.

    SPICE IPHOTO is positive LEAVING DET. Thus VDET/IPHOTO is the
    negative of the usual impedance for current INJECTED into DET.
    With Zin=-VDET/IPHOTO, Zcorrected=Zoriginal/(1+Zin*(Ydevice-j*w*C)).
    This post-processing is CALCULATED, not an additional nonlinear TCAD
    cosimulation. No electrical-admittance interpolation between the four
    supplied device frequencies is used. The optical-current transfer must
    still be included separately to obtain an optical system response.
    """
    rows = []
    for source_path in sorted(Path(grid_root).glob("d*/admittance.csv")):
        folder = Path(run_root) / source_path.parent.name / "ac"
        if not (folder / "run.json").exists() or not (folder / "summary.json").exists():
            continue
        if not json.loads((folder / "run.json").read_text()).get("success"):
            continue
        summary = json.loads((folder / "summary.json").read_text())
        source_summary = source_path.with_name("summary.json")
        expected_hash = summary["conditions"].get("device_summary_sha256")
        if expected_hash and sha256_bytes(source_summary.read_bytes()) != expected_hash:
            raise ValueError(f"Device summary changed after TIA extraction: {source_summary}")
        cap = summary["conditions"]["capacitance_f"]
        raw = parse_raw(folder / "test.raw")
        admittance_hash = sha256_bytes(source_path.read_bytes())
        raw_hash = sha256_bytes((folder / "test.raw").read_bytes())
        frequencies = raw["frequency"].real
        with source_path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if float(row["illumination_scale"]) != 1:
                    continue
                frequency = float(row["frequency_hz"])
                index = int(np.argmin(abs(np.log(frequencies/frequency))))
                if not np.isclose(frequencies[index], frequency, rtol=1e-8, atol=0):
                    raise ValueError("Device frequency absent in SPICE grid; refusing an unreported interpolation")
                current = raw["i(iphoto)"][index]
                z_original = raw["v(out)"][index] / current
                z_in = -raw["v(det)"][index] / current
                y_device = complex(float(row["admittance_real_s"]), float(row["admittance_imag_s"]))
                y_fixed = 2j*np.pi*frequency*cap
                correction = 1/(1+z_in*(y_device-y_fixed))
                z_corrected = z_original*correction
                rows.append(dict(case_id=source_path.parent.name, frequency_hz=frequency,
                                 capacitance_fixed_f=cap, device_y_real_s=y_device.real, device_y_imag_s=y_device.imag,
                                 fixed_y_imag_s=y_fixed.imag, usual_zin_real_ohm=z_in.real, usual_zin_imag_ohm=z_in.imag,
                                 fixed_model_z_magnitude_ohm=abs(z_original), corrected_z_magnitude_ohm=abs(z_corrected),
                                 full_y_relative_difference_pct=100*abs(y_device-y_fixed)/abs(y_device),
                                 corrected_transfer_complex_change_pct=100*abs(correction-1),
                                 corrected_transfer_magnitude_change_pct=100*(abs(correction)-1),
                                 corrected_transfer_phase_change_deg=float(np.angle(correction, deg=True)),
                                 device_admittance_sha256=admittance_hash, tia_raw_sha256=raw_hash,
                                 provenance="CALCULATED linear port correction from SIMULATED TCAD admittance and LTspice port response; not nonlinear cosimulation"))
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if rows:
        with destination.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, list(rows[0]))
            writer.writeheader()
            writer.writerows(clean(rows))
    return rows


def audit_dense_admittance(device_directory, run_root, destination):
    """Keep every directly simulated device Y point; interpolate only dense TIA AC.

    LTspice real and imaginary response components are linearly interpolated
    on its logarithmic frequency axis (160points/decade by default). This
    permits checking additional device frequencies without interpolating the
    simulated device admittance.
    """
    rows, summaries = [], []
    for path in sorted(Path(device_directory).glob("*.csv")):
        case_id = path.stem
        folder = Path(run_root) / case_id / "ac"
        if not json.loads((folder / "run.json").read_text()).get("success"):
            raise ValueError(f"No successful TIA AC run for {case_id}")
        summary = json.loads((folder / "summary.json").read_text())
        cap = summary["conditions"]["capacitance_f"]
        raw = parse_raw(folder / "test.raw")
        f_spice = raw["frequency"].real
        log_axis = np.log10(f_spice)
        def at(values, frequency):
            if not f_spice[0] <= frequency <= f_spice[-1]:
                raise ValueError("Dense admittance frequency outside SPICE range")
            return complex(np.interp(np.log10(frequency), log_axis, values.real),
                           np.interp(np.log10(frequency), log_axis, values.imag))
        device_hash = sha256_bytes(path.read_bytes())
        raw_hash = sha256_bytes((folder / "test.raw").read_bytes())
        case_rows = []
        with path.open(newline="") as handle:
            for item in csv.DictReader(handle):
                if float(item["illumination_scale"]) != 1:
                    continue
                frequency = float(item["frequency_hz"])
                current = at(raw["i(iphoto)"], frequency)
                z_in = -at(raw["v(det)"], frequency)/current
                z_original = at(raw["v(out)"], frequency)/current
                y_device = complex(float(item["admittance_real_s"]), float(item["admittance_imag_s"]))
                y_fixed = 2j*np.pi*frequency*cap
                correction = 1/(1+z_in*(y_device-y_fixed))
                row = dict(case_id=case_id, frequency_hz=frequency,
                           capacitance_fixed_f=cap, device_y_real_s=y_device.real, device_y_imag_s=y_device.imag,
                           fixed_y_imag_s=y_fixed.imag, usual_zin_real_ohm=z_in.real, usual_zin_imag_ohm=z_in.imag,
                           fixed_model_z_magnitude_ohm=abs(z_original), corrected_z_magnitude_ohm=abs(z_original*correction),
                           full_y_relative_difference_pct=100*abs(y_device-y_fixed)/abs(y_device),
                           corrected_transfer_complex_change_pct=100*abs(correction-1),
                           corrected_transfer_magnitude_change_pct=100*(abs(correction)-1),
                           corrected_transfer_phase_change_deg=float(np.angle(correction, deg=True)),
                           device_admittance_sha256=device_hash, tia_raw_sha256=raw_hash,
                           provenance="CALCULATED linear port correction; direct SIMULATED dense deviceY, real/imag LTspice response interpolated onlog-frequency; not nonlinear cosimulation")
                rows.append(row)
                case_rows.append(row)
        within = [r for r in case_rows if r["frequency_hz"] <= 4e6]
        worst = max(within, key=lambda r: r["corrected_transfer_complex_change_pct"])
        summaries.append(dict(case_id=case_id, simulated_device_frequency_points=len(case_rows),
                              sampled_band_low_hz=min(r["frequency_hz"] for r in within), sampled_band_high_hz=max(r["frequency_hz"] for r in within),
                              max_complex_transfer_change_pct_below4mhz=worst["corrected_transfer_complex_change_pct"],
                              worst_frequency_hz=worst["frequency_hz"],
                              max_absolute_magnitude_change_pct_below4mhz=max(abs(r["corrected_transfer_magnitude_change_pct"]) for r in within),
                              max_absolute_phase_change_deg_below4mhz=max(abs(r["corrected_transfer_phase_change_deg"]) for r in within),
                              scope="Maximum among sampled device frequencies≤4MHz; not proof of a continuous-band supremum"))
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if rows:
        with destination.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, list(rows[0]))
            writer.writeheader()
            writer.writerows(clean(rows))
    write_json(destination.with_name(destination.stem + "_summary.json"), summaries)
    return rows, summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("cases_json", type=Path)
    parser.add_argument("--analyses", default="op,ac,noise")
    parser.add_argument("--output", type=Path, default=ROOT / "runs/tia")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--extract-grid", type=Path, help="Derive cases_json from completed d*/summary.json device runs before executing")
    parser.add_argument("--extract-steps", action="store_true", help="Derive six step cases from latest completed device time refinements; use --analyses tran")
    parser.add_argument("--case-ids", help="Comma-separated subset of input case IDs")
    parser.add_argument("--transient-max-step", type=float, help="Explicit transient maximum step in seconds for a convergence check")
    args = parser.parse_args()
    if args.extract_grid:
        export_grid_cases(args.extract_grid, args.cases_json)
    if args.extract_steps:
        export_step_cases(args.cases_json)
    data = json.loads(args.cases_json.read_text())
    settings = {**DEFAULTS, **data.get("settings", {})}
    if args.transient_max_step is not None:
        if args.transient_max_step <= 0:
            raise ValueError("Transient maximum step must be positive")
        settings["transient_max_step_s"] = args.transient_max_step
    model = SOURCE_MODEL.read_text()
    summaries = []
    for case in data["cases"]:
        if args.case_ids and case["case_id"] not in args.case_ids.split(","):
            continue
        validate_case(case)
        path = case.get("optical_pulse_csv")
        if path and not Path(path).is_absolute():
            path = args.cases_json.resolve().parent / path
        for analysis in args.analyses.split(","):
            summaries.append(run_case(case, settings, analysis, model, args.output, path, args.prepare_only, args.force))
    args.output.mkdir(parents=True, exist_ok=True)
    write_json(args.output / "summary.json", dict(input_file=str(args.cases_json.resolve()),
                input_sha256=sha256_bytes(args.cases_json.read_bytes()), summaries=summaries))
    if not args.prepare_only:
        export_metrics_csv(args.output)


if __name__ == "__main__":
    main()
