# Conservative coupled gas sources

`gas_solver,implicit` selects advection followed by a local, coupled backward-Euler
source solve. `gas_model,albanese2016` must also be enabled. The existing transported
variables are retained. The previous explicit path is available with
`gas_solver,legacy` or when the new selector is absent.

## Concentrations and equations

Oxygen `O` is gas volume per blood volume, `S` is Hb saturation, and `T` is oxygen
per tissue volume. Carbon species `Cp,Bp` are per plasma volume and `Cr,Br,H` are
per RBC volume. Gas volumes use the Albanese model's STPD molar-volume convention.

Let `phi_p=1-phi_r`, `qB=beta_r/beta_p`, `qC=W_r/W_p`. Internal exchange fluxes,
per blood volume per second, are

    JC = phi_r/tau_C_rp * (qC*Cp-Cr)
    JB = phi_r/tau_B_rp * (qB*Bp-Br).

The opposite-direction times are derived, not independently fitted:

    tau_C_pr = tau_C_rp*phi_p/(phi_r*qC)
    tau_B_pr = tau_B_rp*phi_p/(phi_r*qB).

Albanese's coupled curves supply total equilibrium contents. Dissolved CO2 is
subtracted and the remaining content is partitioned using `albanese_options`:

    etaB = B_eq-(phi_p*Bp+phi_r*Br)
    etaH = H_eq-phi_r*H
    RB = phi_r*etaB/tau_B_reaction
    RH = phi_r*etaH/tau_H
    Cp' = -JC/phi_p + RQ*phi_t/(phi_c*phi_p)*M(T)   [systemic]
    Cr' = (JC-RB-RH)/phi_r
    Bp' = -JB/phi_p
    Br' = (JB+RB)/phi_r
    H'  = RH/phi_r.

Thus `phi_p*(Cp+Bp)+phi_r*(Cr+Br+H)` is conserved by internal exchange/reaction.
One Michaelis-Menten rate `M(T)=Mmax*T/(T+C50)` at the implicit iterate supplies
both the oxygen sink and the carbon production. No separately evolved tissue CO2
pool is used in this path.

With `A_H=C_RBC*Z`, `dS=S-S_eq` and `JO=A_H*dS/tau_O`:

    O' = JO-F
    S' = -dS/tau_O
    T' = phi_c/phi_t*F-M(T),

where `F=kc/hc*(S_cap/V_cap)*(PO2-PtO2)` with the same dilation multiplier in
both exchange partners. Internal transfer conserves `O+A_H*S+phi_t/phi_c*T`.

In the lung, tissue transfer/metabolism are replaced by
`K_O*sin²(pi*x/L)*(PO2_alv-PO2)` and
`K_C*sin²(pi*x/L)*(PCO2_alv-PCO2)` in the receiving O2 and plasma-CO2 equations.
`K_O,K_C` and pressure differences use Pa. Both equilibrium curves are recomputed
inside each Newton iteration, preserving their O2-CO2 coupling.

## Input rows

Rows may appear in either `O2_parameters.txt` or `CO2_parameters.txt`. The selected
model supplies all required new groups. New groups are applied after both files
have been read, so their position relative to the legacy parameter rows does not
change precedence. Existing per-region `Mmax`, `C50`, `Z`, wall thickness and
permeability remain configurable by the existing records. In particular `Mmax`
continues to come from each region's model CSV.

| Record | Ordered values |
|---|---|
| `gas_solver` | `implicit` or `legacy` |
| `gas_fractions` | capillary fraction, tissue fraction, RBC fraction |
| `gas_partition` | beta_plasma, beta_RBC |
| `gas_solubility` | alpha_O2_blood/mmHg, alpha_O2_tissue/mmHg, alpha_CO2_water mol/(L mmHg), W_plasma, W_RBC, Pa/mmHg |
| `gas_relaxation` | tau_O2_systemic, tau_O2_pulmonary, tau_CO2_RBC_plasma, tau_HCO3_RBC_plasma, tau_HCO3_reaction, tau_HbCO2 [s] |
| `gas_pulmonary` | alveolar PO2/mmHg, alveolar PCO2/mmHg, Kmax_O2/(Pa s), Kmax_CO2/(Pa s) |
| `gas_numerics` | scaled Newton tolerance, iteration limit, maximum binary step subdivisions |

Missing required groups, non-finite/negative parameters, inconsistent fractions
and invalid solver controls fail before integration. CO2 solubilities are derived
as `alpha_CO2_water*Vm*W/Pa_per_mmHg`, consistent with the exchange ratio. The
normalization `beta_plasma=1` is arbitrary; only its ratio to `beta_RBC` is used.

The adopted numbers are the 2026-11-18 literature-derived starting estimates,
not a calibrated physiological parameter set. In particular the O2 time uses
Richardson 2020's 23 C mean with an approximate 1.3 rate correction, the reaction
times use a fixed-target local interpretation of measured rates, and the lung
coefficients normalize chosen whole-lung conductances by the model's effective
storage volume. Their assumptions and sources are recorded in the parameter matrix.

## Time stepping and scheduling

All gas inlet nodes are updated first. `advance_coupled_gas` then updates each
connected lumped model once per hydraulic edge step. The new path bypasses the
two old separate capillary updates. It uses one shared blood/RBC/plasma velocity,
corresponding to xi=phi; it does not apply those volume corrections twice.

Upwind advection is subdivided to CFL<=0.8. RBC number and Hb-bound oxygen are
advected together before recovering saturation. The stationary tissue field is
not advected. Backward Euler uses a damped Newton solve with physical bounds;
failed source steps are retried with binary subdivision, without committing a
partial failed step. An unresolved nonlinear failure raises an error instead of
clipping concentrations or suppressing the failure.

The optimized implementation uses analytic source derivatives and an exact
algebraic reduction to three local Newton unknowns. All eight original equations
are checked before acceptance; the full solver remains as a fallback. GCC builds
of `custom_output` enable optional cell-level OpenMP parallelism. See
[gas_optimization.md](gas_optimization.md) for conservation checks, staged timings
and the recommended six-thread configuration on whserver. No physical parameters
or solver tolerances were changed by these optimizations.

## Short validation

`tests/test_conservative_gas.cpp` checks isolated conservation, beta normalization,
equilibrium stationarity, stiff positive relaxation, metabolic stoichiometry,
pulmonary external-flux balance and input rejection. `tests/smoke_conservative_gas.py`
runs the full network for only 2 or 5 ms, including both parameter-file locations,
legacy regressions, all 136 saved histories and malformed-input cases.

Validated on whserver on 2026-11-18:

- `projects/custom_output`: full Clang and GCC builds passed.
- Core tests passed both normally and with AddressSanitizer/UndefinedBehaviorSanitizer.
- Full-network implicit tests requested 2 and 5 ms; all 136 saved histories were
  finite and nonnegative, and Hb saturation remained in [0,1]. The existing
  scheduler saves its final sample 0.5 ms beyond the requested stop time.
- Moving the new records from the O2 file into the CO2 file produced matching
  terminal states for all histories.
- Both legacy paths passed 2 ms regression smoke tests; four malformed/missing
  input cases were rejected with the expected diagnostic.

Logs and machine-readable results are on whserver under
`/home/wr/first_blood-runs/dev_gas_conservative_20261118`. These short checks do
not establish physiological calibration or long-time accuracy/stability.
