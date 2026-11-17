# Albanese equilibrium sources

The existing eight transport fields (RBC, Hb oxygen saturation, dissolved O2,
and the five CO2 components) are retained. Hydraulics, advection, membrane
exchange, tissue metabolism and alveolar exchange retain their existing form.
Only local chemical-equilibrium targets use the coupled curves in Albanese
et al. (2016), DOI 10.1152/ajpheart.00230.2014, equations A46-A49.

## Parameter files

The following optional rows can be placed in `O2_parameters.txt` or
`CO2_parameters.txt` in the model directory, alongside the existing rows.
The usual whitespace removal and comma-separated format apply. Put each
setting in one file; later rows replace earlier settings. Legacy rows such as
`peripheries`, `haemoglobin`, `pulmonary` and `perifco2` remain supported.

```text
gas_model, albanese2016
albanese_o2, 9, 0.3836, 0.03198, 0.008275, 14.99
albanese_co2, 86.11, 1.819, 0.05591, 0.03255, 194.4
albanese_options, 0.05, 22.414, 0.201
```

- `gas_model`: `legacy` (default when absent) or `albanese2016`.
- Curve rows: saturation content [mmol/l blood], h [-], alpha [1/mmHg],
  beta [1/mmHg], K [mmHg]. The exponent is **1/h**, not h.
- Options: carbamino fraction of total CO2, molar gas volume [l/mol, STPD],
  reference Hb-bound oxygen capacity [m3 O2/m3 blood]. Defaults are shown.
- The existing CO2 Henry coefficients remain in 1/Pa; their partial pressures
  are divided by the solver's Pa/mmHg constant before curve evaluation.
- `mmol/l * 0.001 * molar_volume_l_mol` converts to m3 gas/m3 blood.

All eight gas/RBC transport switches must be on in a participating lumped
model. Existing `perifco2` parameters are still required for CO2 transport;
missing files now produce a controlled error instead of uninitialized values.
Nonpositive time constants/solubilities, invalid fractions and malformed
Albanese parameter rows are rejected with an error message.

## Source construction and additional assumptions

The oxygen reference curve is scaled by `(RBC * Z) / reference_hb_capacity`
to match the current field's Hb oxygen capacity. The default reference is
15 g Hb/dl times 1.34 ml O2/g divided by 100 = 0.201. The equilibrium
saturation is `(scaled_total_O2 - dissolved_O2)/(RBC*Z)`. Its projection to
[0,1] bounds extrapolation of the empirical curve; the curves and the capacity
must still be calibrated consistently. Carbon uses the reference empirical
curve, whose parameters correspond to the paper's reference blood composition.

The empirical CO2 curve does not identify the five individual components.
The additional partition assumption is:

```text
HbCO2_target_blood = f_Hb * total_CO2
HCO3_target_blood = (1-f_Hb)*total_CO2
                    - (phi_pla*alpha_CO2_pla + phi_rbc*alpha_CO2_rbc)*PCO2_Pa
```

Negative target bicarbonate is rejected. The previous bicarbonate partition
between plasma and RBC and the existing relaxation times are retained.
The Hb oxygen reaction transfers exactly opposite amounts to saturation and
dissolved O2. The RBC carbon sink is paired with equal bicarbonate and
carbamino sources on the same blood-volume basis. No external gas is created
by these internal reactions.

`first_blood::run` captures both gases' local reaction departures before either
gas is advanced at an attached lumped compartment. These cached departures are
used for the interior and both outgoing endpoint branches of systemic and
pulmonary updates. This is an explicit local source evaluation on the common
pre-transport state, not a new transported total-gas field or a full reproduction
of Albanese's ventilation, shunt and neural-control model.

## Build and short checks

```sh
mkdir -p bin
cd projects/custom_output
make -j2 -f custom_output.mk CXX=clang++
cd ../..
g++ -std=c++17 -O2 tests/test_albanese_gas.cpp -o bin/test_albanese_gas.out
./bin/test_albanese_gas.out
python3 tests/smoke_albanese.py
```

The smoke runner copies the selected network to disposable test directories,
sets simulation time to 0.05 seconds, checks legacy and Albanese modes, output
values and invalid-input failures, and never alters the original 130-second
model duration. No long integration is part of these checks.

The transport output history now includes its initial t=0 state, matching the
hydraulic time vector. This fixes the previously observed final-interval
out-of-bounds interpolation. `custom_output` reports model/run errors with a
nonzero exit status instead of claiming a successful run.

Sanitizer checks also exposed two pre-existing startup problems on this model:
optional node CSV columns were indexed without size checks, and the cardiac
period/phase was initialized only after its first use by the elastance law.
Optional columns are now guarded, and the period, heart rate and zero initial
phase are set before the first elastance evaluation.
