# Gas solver optimization, 2026-11-19

Base: `7a2621dc7e709a599d6ec2e8d784fb19de40d2d3`, branch `dev_gas`.
All physical parameter files, including every relaxation time, are unchanged.
The transported fields, source equations, backward-Euler discretization,
Newton tolerance and step-subdivision limits are retained.

## Changes, implemented and validated in order

1. **Analytic derivatives.** `AlbaneseGas` optionally supplies derivatives of
   the same equilibrium values used previously, including the active saturation
   bounds. `ConservativeGas::rates` assembles the complete analytic source
   Jacobian. The full eight-variable Newton solver uses it instead of eight
   finite-difference perturbations. At zero pressure the existing finite-difference
   path remains available because the empirical carbon curve has a singular
   pressure derivative there.
2. **Algebraic reduction.** For trial free O2, plasma CO2 and tissue O2, the Hb
   saturation, carbamino pool, two bicarbonate pools and RBC free CO2 are obtained
   from their exact linear backward-Euler equations. Only three nonlinear
   unknowns remain; pulmonary tissue O2 is a passive identity equation, so only
   two pulmonary unknowns are nonlinear. The reduced Jacobian is analytic.
   Every accepted solution must satisfy **all eight original equations** at
   the original scaled tolerance and have physical bounds. Failed reduced solves
   retry the full analytic solver, then the existing time-subdivision mechanism.
   Eliminated linear pools may be outside physical bounds at intermediate Newton
   trials; only positive free pressures/nonnegative tissue O2 are used for curve
   evaluation, and the full converged state is checked before acceptance.
3. **Parallel cells.** After the old advection state is prepared, each capillary
   cell writes its own tentative state. OpenMP uses a static schedule for meshes
   with at least 64 points. Exceptions are captured separately per cell and the
   lowest failing index is reported after the workers finish. A failing substep
   commits none of its tentative gas, RBC or tissue fields. The network event
   ordering and inter-model coupling remain serial.

No hydraulic QR replacement, transport-vector-copy optimization, buffer reuse or
constant-precomputation optimization was included.

## Measured wall time

whserver, Ryzen 5 3600 (6 physical cores / 12 hardware threads), GCC `-O3` solver,
`custom_output`, model `Bathsheba_CO2_pul_2026_03_11`. Original output configuration:
p8 and heart_kim_lit transport, save interval 0.01 s. Runs used temporary model
copies. Timing includes program startup and output, but excludes model copying,
compilation and validation. Stages were measured sequentially, with three repeats
each; medians are below. `OMP_DYNAMIC=FALSE`, `OMP_WAIT_POLICY=ACTIVE`.

| Validated stage | Threads | Wall time for 0.2 s simulation | Speedup vs base | Wall / simulation time |
|---|---:|---:|---:|---:|
| Original | 1 | 14.8300 s | 1.00x | 74.15 |
| Analytic 8-variable Jacobian | 1 | 9.7304 s | 1.52x | 48.65 |
| Reduced nonlinear system | 1 | 7.6496 s | 1.94x | 38.25 |
| Reduced + OpenMP, serial setting | 1 | 7.6753 s | 1.93x | 38.38 |
| Reduced + OpenMP | 2 | 6.1132 s | 2.43x | 30.57 |
| Reduced + OpenMP | 4 | 5.3059 s | 2.79x | 26.53 |
| Reduced + OpenMP | 6 | **5.0838 s** | **2.92x** | **25.42** |

A separate 1 s simulation, covering more than the 0.9242 s cardiac period, gave:

| Stage | Wall time | Speedup vs base |
|---|---:|---:|
| Original | 75.7511 s | 1.00x |
| Analytic Jacobian | 50.1094 s | 1.51x |
| Reduced system | 39.5417 s | 1.92x |
| Parallel reduced system, 6 threads | **26.4951 s** | **2.86x** |

The 1 s results are single verification runs, not three-repeat medians. The
existing scheduler may slightly overrun the requested simulation endpoint.

Six threads are recommended on whserver. Using the runtime's default thread
count (12 on this machine) gave a median 6.3253 s for the 0.2 s test. Explicit
12-thread ACTIVE runs gave 6.3342 s. Six threads without an explicit wait policy
gave 5.1028 s, close to 5.0838 s with ACTIVE. More hardware threads did not help
this workload. The measured improvements are lower than the earlier optimistic
estimates; these tables report actual results, not those estimates.

Machine-readable repeats, compiler checks, parameter hashes and validation logs
are in [benchmarks/gas_optimization_20261119.json](benchmarks/gas_optimization_20261119.json).
Full server-side stage snapshots, executables and logs are in
`/home/wr/first_blood-runs/gas_optimization_20261119`.

## Conservation and equivalence checks

Every source-stage check used 160 seeded systemic/pulmonary states and four
time steps (1 us, 0.1 ms, 1 ms, 5 ms): **640 source solves per stage**. The tests
compare with the original finite-difference implementation and independently
evaluate all original backward-Euler residuals. The conservation quantities are

    oxygen = O + capacity*S + (phi_t/phi_c)*T
    carbon = phi_p*(Cp+Bp) + phi_r*(Cr+Br+H).

Systemic balance includes the same Michaelis-Menten sink and RQ-weighted source;
pulmonary balance includes the backward-Euler external alveolar fluxes.

| Check | Analytic stage | Reduced / parallel stage |
|---|---:|---:|
| Maximum absolute source O2/CO2 balance error | 1.64e-16 | 2.22e-12 |
| Maximum scaled original equation residual | 9.97e-11 | 9.97e-11 |
| Maximum scaled source-state difference from original | 9.93e-11 | 9.92e-11 |
| Maximum analytic Jacobian comparison error | 3.86e-8 | 3.86e-8 |
| Maximum reduced Jacobian comparison error | n/a | 2.60e-9 |

All 640 reduced cases converged without the full-solver fallback. The source
balance error after reduction is larger than the analytic full-system result,
but remains within the unchanged residual tolerance; no tolerance was relaxed.

`test_gas_capillary.cpp` additionally tests systemic/pulmonary capillaries,
forward/reverse advection and spatially varying RBC concentration. With the
external sources closed, the change in the complete discrete capillary totals
matches the upwind inlet-minus-outlet flux. Maximum absolute O2/CO2 sum-balance
error was **3.93e-12**, and relative RBC sum-balance error **8.50e-18**. Serial and
six-thread fields were exactly identical. Injected worker failures identified
the failing advected cell and left the failing substep's fields uncommitted.

All 136 saved histories were finite, nonnegative and Hb saturation remained in
[0,1]. Analytic-stage saved values were identical to baseline at output precision.
Reduced/parallel maximum scaled saved-value difference was **4.73e-8** for the
0.2 s runs and **7.38e-8** for the 1 s runs (scaling floor 1e-4). Saved arrays were
exactly identical between the 1/2/4/6-thread configurations and the extra thread-
policy checks. These comparisons include every saved row, not only endpoints.

Also passed:

- zero-pressure stationary states and the full-solver fallback;
- existing Albanese curves, coupling, units and legacy source tests;
- both legacy integration paths and parameter-file-location/input-error smoke tests;
- GCC AddressSanitizer + UndefinedBehaviorSanitizer source and six-thread capillary
  tests; leak checking was disabled for consistency with existing solver ownership;
- a separate Clang serial build and 1 s result comparison;
- unchanged O2/CO2 parameter-file bytes and `git diff --check`.

## Building and running

From `projects/custom_output`, GCC enables OpenMP by default:

```sh
make -j2 -f custom_output.mk CXX=g++
OMP_NUM_THREADS=6 OMP_DYNAMIC=FALSE OMP_WAIT_POLICY=ACTIVE ./custom_output.out Bathsheba_CO2_pul_2026_03_11
```

The command above uses the original model's configured duration (130 s); the
validation runs changed duration only in temporary model copies. `OMP_NUM_THREADS=1`
selects serial execution. `OPENMP_FLAGS=` builds without OpenMP. Clang defaults to
serial so it remains usable without an installed OpenMP runtime; with a suitable
runtime installed, `OPENMP_FLAGS=-fopenmp` may be supplied explicitly. Other project
makefiles are unchanged and use the serial fallback unless they enable OpenMP.

Standalone checks from the repository root (Bash):

```sh
g++ -std=c++17 -O2 tests/test_conservative_gas.cpp -o /tmp/test_conservative_gas
/tmp/test_conservative_gas
g++ -std=c++17 -O2 tests/test_gas_optimization.cpp -o /tmp/test_gas_optimization
/tmp/test_gas_optimization
g++ -std=c++17 -O2 -fopenmp tests/test_gas_capillary.cpp \
  bin/{file_io,first_blood,moc_edge,moc_node,solver_lumped,solver_lumped_io,solver_moc,solver_moc_io,statistics}.o \
  -o /tmp/test_gas_capillary
/tmp/test_gas_capillary
```

For the optional independent old-Newton comparison:

```sh
git show 7a2621d:source/conservative_gas.h > /tmp/reference_conservative_gas.h
g++ -std=c++17 -O2 -Isource \
  -DGAS_REFERENCE_HEADER='"/tmp/reference_conservative_gas.h"' \
  tests/test_gas_optimization.cpp -o /tmp/test_gas_reference
/tmp/test_gas_reference
```
