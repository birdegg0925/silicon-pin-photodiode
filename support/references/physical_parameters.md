# PIN model: implemented inputs, sources and limits

This note describes the model used by the current PIN study. The device equations and fixed defaults are in [pin_device.py](../scripts/pin_device.py); the nine thickness/bias combinations are in [protocol.json](../../protocol.json). Saved run configurations in `pin/runs/*/summary.json` record the inputs used for each result. All device and circuit responses are simulations. Analytical checks and the combined receiver response are calculations.

## Geometry and operating conditions

The simulated structure is planar p+–n−–n+ silicon. The lightly doped middle region is the absorber; it is not assigned zero residual doping. The one-dimensional domain includes both contact layers.

| Input | Implemented value |
|---|---|
| Absorber thickness | 7.9, 20.8 and 35.6 µm |
| Reverse bias | 1.55, 4.15 and 7.86 V, at each thickness |
| Electrical and illuminated area | 0.0025 cm² = 0.25 mm² |
| Bulk donor density | 1e13 cm^-3 |
| Contact doping scale | 1e17 cm^-3 |
| Contact-layer thickness | 0.3 µm on each side |
| Junction smoothing length | 0.02 µm |
| Temperature | 300 K |
| SRH electron and hole lifetimes | 1 µs each |
| Incident optical power | 1 µW at the illuminated surface |
| Wavelength | 850 nm |

These geometry, doping, lifetime and operating values are model assumptions. The contact transitions use complementary-error-function profiles rather than a simulated fabrication process. Contacts are ideal ohmic boundaries. Surface/perimeter leakage, process variation, contact resistance, band-gap narrowing and avalanche multiplication are absent. The implementation rejects temperatures other than 300 K because the optical and transport inputs are fixed at that temperature.

The mesh uses nominal spacings of 0.005 µm near contacts and junctions and 0.10 µm in the absorber interior. `mesh_scale` multiplies those spacings. The archived refinement at `mesh_scale=0.5` checks DC responsivity and capacitance for the selected reference; it is not a transient spatial-mesh convergence study.

## Optical inputs and generation

The absorption coefficient and refractive index are taken from M. A. Green, *Self-consistent optical parameters of intrinsic silicon at 300 K including temperature coefficients*, Solar Energy Materials and Solar Cells 92 (2008), 1305–1310: [DOI](https://doi.org/10.1016/j.solmat.2008.06.009), [full paper](https://subversion.xray.aps.anl.gov/AGIPD_TCAD/PyTCADTools/Attenuation/Silicon_absorption_coeff_visible_photon.pdf), Table 1, printed page 1307.

| Quantity | Value and use |
|---|---|
| Absorption coefficient alpha | 534 cm^-1; used directly in generation |
| Refractive index n | 3.641; used in front-surface reflectance |
| Tabulated extinction coefficient | 0.004, rounded in the source |
| Extinction coefficient used in reflectance | alpha × wavelength / (4π), approximately 0.003612 |
| Absorption length | 1/alpha, approximately 18.73 µm |
| Bare-surface normal-incidence reflectance | Approximately 0.324, calculated from n and the reconstructed extinction coefficient |

The reflectance calculation is `R = ((n-1)^2+k^2)/((n+1)^2+k^2)`. Generation is

    G(x) = alpha * (1-R) * P_inc / (A * h*c/lambda) * exp(-alpha*x)

Here x starts at the outer p+ illuminated surface, and generation is evaluated throughout the silicon, including both contact layers. There is no separate excluded dead layer. The front-transmission factor is applied once. Lengths in the solver are centimetres, area is cm², and G is pairs cm^-3 s^-1.

The model assumes one electron–hole pair per absorbed photon. The analytical responsivity ceiling is q*lambda/(h*c), approximately 0.686 A/W at 850 nm; including bare-surface reflection and infinite absorption gives approximately 0.464 A/W. The finite-thickness photon check includes the absorber and both contact layers. These values are analytical bounds, not measured detector performance.

The optical model has no antireflection coating, texture, optical cavity, interference, back reflector or guided mode. Green reports absorption-coefficient uncertainty of approximately 4% near this wavelength; the numerical results do not calibrate that material uncertainty.

## Electrical material properties

The implemented relative permittivity is 11.7. With epsilon_0 = 8.8541878128e-14 F/cm, this gives epsilon_Si approximately 1.036e-12 F/cm. The literature source is W. C. Dunlap Jr. and R. L. Watters, *Direct Measurement of the Dielectric Constants of Silicon and Germanium*, Physical Review 92 (1953), 1396: [DOI](https://journals.aps.org/pr/abstract/10.1103/PhysRev.92.1396). This is a bulk-material input, not a contact-layer calibration.

The implemented intrinsic carrier density is 9.65e9 cm^-3 at 300 K, from P. P. Altermatt et al., *Reassessment of the intrinsic carrier density in crystalline silicon in view of band-gap narrowing*, Journal of Applied Physics 93 (2003), 1598–1604, DOI 10.1063/1.1529297: [full paper](https://iis-people.ee.ethz.ch/~schenk/nijap.pdf). The model uses this fixed value for n_i, n1 and p1; it does not implement a temperature-dependent material model.

Field-dependent mobilities use the 300 K Table 4 fit from C. Scharf and R. Klanner, *Measurement of the drift velocities of electrons and holes in high-ohmic <100> silicon*: [paper](https://arxiv.org/pdf/1503.08656). With E in V/cm and mobility in cm²/(V s), the code evaluates

    Eabs = sqrt(E^2 + 1)
    mu_e = 1530 / (1 + 1530*Eabs/1.03e7)
    z = max(Eabs - 2640, 0)
    mu_h = 1 / (1/464 + 9.57e-8*z - 3.31e-13*z*z)

The regularised field magnitude keeps the expressions smooth near zero. The cited fit extends to 50 kV/cm; the study's physics check requires the maximum simulated field to remain below that value. The same empirical mobility laws are used in the thin contact regions. Their use there is an approximation: doping-dependent mobility and contact transport are not separately modelled.

DEVSIM's [simple_physics.py](https://github.com/devsim/devsim/blob/438ccb280aec7308522cc0e4b01de3c3033bd78b/lib/devsim/python_packages/simple_physics.py) supplies the Poisson, drift–diffusion, SRH and ohmic-contact helpers. The PIN implementation overrides the helper defaults for permittivity, intrinsic density and lifetime, and supplies the field-dependent mobility expressions above. Saved records identify DEVSIM 2.11.0 with UMFPACK. The solver uses centimetre units; terminal current, charge and admittance are multiplied by area once.

## Transient response and numerical checks

The device transient is a 0.9 to 1.0 µW optical step at fixed reverse bias. It uses DEVSIM BDF1 integration over a 2 µs observation window, with finer time intervals at the start. Reported response uses the external total current, including displacement current. Rise time is the first 10–90% crossing interval; settling is the first sample after the last departure from ±2% of the optical-current increment.

[run_pin_revision.py](../../run_pin_revision.py) evaluates the full 3×3 grid. Time-step refinement applies to the 20.8 and 35.6 µm cases, halving the time scale until both rise and settling times change by at most 2% between consecutive runs, or the available refinement levels are exhausted. The 7.9 µm timing records do not pass this consecutive-refinement qualification and are excluded from the main timing comparison.

The saved checks cover:

- Terminal-current balance, collected-photon fraction, maximum field and transient endpoint error.
- A DC responsivity/capacitance mesh refinement for the selected reference.
- A separate electron- and hole-inventory capacitance derivative compared with the 100 kHz terminal-admittance capacitance. These derivatives are compared separately, not added.
- Sampled-frequency loading corrections for the detector-to-TIA interface.

The local contact-field charge derivative also remains in the archived outputs. Screening at the ohmic contacts makes it unsuitable as the total-diode-capacitance estimate; the carrier-inventory crosscheck is the relevant separate check. Analytical depletion-width and epsilon*A/W estimates provide context but do not replace those simulated checks.

The included runs use one SRH lifetime and one mobility law. No 0.1/10 µs lifetime-sensitivity sweep or alternative-mobility comparison is included. Dynamic spatial-mesh convergence, measured-device calibration and fabrication validation are outside the recorded results.

## Detector-to-receiver model

[pin_receiver_revision.py](../../pin_receiver_revision.py) uses the seven-MOS TIA model and fixed selected settings in `tia/`: 81.6 kΩ feedback resistance, 0.36 pF feedback capacitance, 1.8 V core supply, 0.9 V reference, 1 pF internal compensation and 1 pF output load. The detector reverse-bias rail is separate from the amplifier supply. Its physical implementation, noise and power conversion are not modelled.

The interface replays the simulated external photocurrent after subtracting dark current, then adds a separate dark-current source. Detector capacitance is fixed at its illuminated 100 kHz operating-point value. Circuit voltage changes do not feed back into the nonlinear device solution.

The optical receiver transfer function is calculated from the simulated optical-step increments and the TIA's complex AC response. A DC responsivity or terminal-admittance capacitance alone is not an optical bandwidth. The admittance checks apply a linear port-loading correction at sampled device frequencies; the saved 4 MHz and 8 MHz summaries describe different sampled bands, not a continuous-band guarantee. This remains a reduced model, not a self-consistent nonlinear TCAD/SPICE co-simulation.

Large raw waveforms and some original plots are absent from the repository; [the manifest](../../packaging_manifest.json) lists them. Scalar summaries and validation records are included. The checks in [verification_status.json](../../verification_status.json) cover the stated single-case DC/AC reproductions and consistency checks; the full device/transient and receiver sweeps were not rerun in that verification.
