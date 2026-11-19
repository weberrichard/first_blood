# Junction optimization and current runtime profile — 2026-11-19

Baseline: `0944bbe422e3b209f6cfb1e183b514a1aa96efba`, branch `dev_gas`.
Project `custom_output`, model `Bathsheba_CO2_pul_2026_03_11`, whserver AMD Ryzen 5 3600.
GCC `-O3`, six OpenMP threads, `OMP_DYNAMIC=FALSE`, `OMP_WAIT_POLICY=ACTIVE`.
All model files, gas time constants, Newton iteration limits and stopping tolerances are unchanged.

## Changes, applied and measured separately

1. The accepted scalar junction flow is evaluated from the matching incident edge context already prepared for the current event. The old final boundary call redundantly reconstructed characteristics and evaluated an unused finite-difference derivative. A flow-only method now avoids those operations. The final flow is evaluated at the accepted, updated pressure, rather than reused from the preceding Newton iterate.
2. In the prepared flow and analytic derivative, fractional area powers are replaced by equivalent square-root expressions. For Olufsen, with `d=1-(P-P0)/F`, `A=An/d²` gives `A^(-1/4)=An^(-1/4)*sqrt(abs(d))`. For the linear wall law `A=r²`, `A^(1/4)=sqrt(abs(r))`. The nominal power is computed once per prepared context. Exceptional area underflow/overflow retains the original pow expression. Absolute values preserve the existing squared-law branches even for nonphysical Newton trial pressures.

No solver replacement, approximate constitutive law, cross-event characteristic cache or tolerance relaxation was introduced.

## Native timings

Five 1 s simulations per version, with cyclically rotated execution order: baseline/flow/sqrt, flow/sqrt/baseline, sqrt/baseline/flow. Timings include normal model loading and output. Compilation, instrumentation and validation are excluded.

| Version | Median wall time for 1 s | Time reduction vs previous | Total time reduction | Maximum field NRMS difference vs baseline |
|---|---:|---:|---:|---:|
| baseline | 12.302073 s | 0.000% | 0.000% | 0% |
| flow | 11.607104 s | 5.649% | 5.649% | 3.997999931e-08% |
| sqrt | 11.523180 s | 0.723% | 6.331% | 4.122874739e-08% |

Total speedup: **1.0676x**. The wall/simulated-time ratio is **11.5232:1**.

The independent five-repeat 0.2 s medians were baseline: 2.472006 s, flow: 2.333345 s, sqrt: 2.310968 s.

Second-step gain is small. The paired 1 s time reductions for sqrt vs flow were 1.673%, -0.131%, 0.188%, -0.279%, 1.370%. These are observations on one machine, not a universal speedup guarantee.

## Solution and conservation checks

Complete spatial states were compared at 0.05, 0.1, 0.2, 0.5 and 1 s: pressure, flow, area, velocity, oxygen, carbon components, RBC and saturation. NRMS% = `100*sqrt(sum((new-reference)^2)/sum(reference^2))`. Pointwise relative errors use a denominator floor of `1e-6*max(abs(reference_field))` to handle zero crossings. All checked values were finite; concentrations/areas were nonnegative and saturation bounded.

| Check, through 1 s | Baseline | Flow reuse | Plus square roots |
|---|---:|---:|---:|
| Maximum junction continuity residual, m³/s | 1.69042365e-17 | 1.72557552e-17 | 1.95300387e-17 |
| Maximum normalized junction imbalance, % | 6.53093557e-10 | 6.72199655e-10 | 1.0167949e-09 |
| Maximum concentrated-node continuity residual, ml/s | 2.27373675e-13 | 2.27373675e-13 | 2.27373675e-13 |
| Maximum mixed-unit hydraulic equation norm | 8.46722416e-09 | 8.4672238e-09 | 8.46722492e-09 |

The external audit uses the original boundary functions and recomputes all concentrated hydraulic equations, independently of the optimized prepared flow law. Each version checked 969537 scalar-junction balances and 1613678 concentrated-node balances. No Newton nonconvergence was observed.

Final maximum field NRMS difference: **4.122874739e-08%**. Maximum pointwise relative difference: **5.632656855e-05%**. Maximum absolute pressure difference: 0.00011758288 Pa. O2/carbon NRMS differences: 1.56642284e-08% / 3.22386026e-09%.

The 136 complete native saved histories through 1 s were checked on every repeat. baseline: maximum per-file NRMS 0%; flow: maximum per-file NRMS 2.76921105e-07%; sqrt: maximum per-file NRMS 0%. Exact agreement at output precision does not imply bitwise agreement of the internal floating-point states.

Additional checks passed:

- Linear and Olufsen boundary laws, both ends, 524 physiological pressure cases and 28 signed wall-law trials: maximum flow difference 8.68e-19 m³/s and analytic derivative error 6.74e-10 relative to central differences. The accepted flow matches the residual flow exactly in these cases.
- 640 stiff gas-source steps: original backward-Euler residual at most 9.97e-11 and O2/CO2 balance error 2.22e-12.
- Capillary discrete O2/CO2 balance error at most 3.93e-12; RBC relative balance error 8.51e-18; serial/parallel capillary states identical.
- GCC AddressSanitizer and UndefinedBehaviorSanitizer unit/integration tests; 2/5 ms full-network, legacy-path and malformed-input smoke tests. Existing ownership/leak behavior was outside scope; leak checking was disabled.
- Independent Clang serial build and complete 1 s saved-history comparison.

## Fresh full-program profile

Instrumentation is confined to an external source copy; the production executable has no timers. Native and instrumented saved outputs agree exactly at output precision for both 0.2 and 1 s. Three repeats per duration and executable were measured.

Measurement overhead: 5.06% at 0.2 s and 6.61% at 1 s. Percentages below are approximate fractions of the instrumented `run()` wall time, not summed CPU time across workers. Internal subregions are included in their parent and are not added again.

| Nonoverlapping phase | Median seconds during 1 simulated second | Fraction of run |
|---|---:|---:|
| Network junctions and hydraulic boundaries | 1.518565 | 12.427% |
| Concentrated hydraulic Newton solves | 3.348775 | 27.404% |
| Coupled capillary gases, including parallel synchronization | 4.372447 | 35.781% |
| Network O2/CO2/RBC transport | 1.483331 | 12.139% |
| MacCormack interior points | 0.878727 | 7.191% |
| Next-event selection | 0.319931 | 2.618% |
| Transport history storage | 0.069980 | 0.573% |
| Edge state and CFL timestep updates | 0.050649 | 0.414% |
| Other control and initialization within run | 0.177527 | 1.453% |

The checked LU/QR linear solve is included in concentrated hydraulics: 2.421207 s, 19.814% of run. It is not an additional top-level category.
Scalar junction calls: 969537; mean Newton iterations 2.591008, maximum 4. General junction QR calls: 0.

Reproduction data: [benchmarks/junction_optimization_20261119.json](benchmarks/junction_optimization_20261119.json). Full traces, build logs, native binaries and external audit/profile copies are retained on whserver under `/home/wr/first_blood-runs/boundary_optimization_20261119`. Validation covers up to one simulated second per run; long-term drift was not assessed in this task.
