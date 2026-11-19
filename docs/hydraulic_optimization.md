# Hydraulic and transport optimization, 2026-11-19

Baseline: `584396b0e62eaa5f96fa6cbd86c2eafb27d5fc12`, branch `dev_gas`.
This is the already optimized gas solver, before the hydraulic/transport work.
All 55 model files are byte-for-byte unchanged, including every gas relaxation
time. The gas equations, hydraulic constitutive laws, CFL control and nonlinear
iteration limits are retained.

## Accepted changes, implemented and tested separately

1. **Junction warm start.** `moc_node` stores current pressure and flow separately
   from optional output history. Initialization, solved junctions and prescribed
   boundaries update the current state. Scalar and general junction Newton solves
   start from that state even when hydraulic history saving is disabled.
2. **Prepared scalar boundary laws.** During one scalar junction Newton solve,
   old-state characteristic/source data are prepared once for each incident edge.
   A trial pressure evaluates the same linear/Olufsen flow law and its analytic
   pressure derivative. The original boundary functions remain available and are
   used independently by the conservation audit. Preparation is local to a single
   event; no state-dependent characteristics are cached across time steps.
3. **Checked partial-pivot LU.** Smooth concentrated models use Eigen partial-pivot
   LU. The linear residual must be finite and no larger than
   `1e-10*max(1, norm(f))`, otherwise the original column-pivot QR is used. Models
   containing diode, quadratic-valve or pump elements retain QR. The original
   nonlinear stopping procedure, including its final correction, is retained.
4. **Transport references.** Velocity, area and new-velocity getters return const
   references; the 1D transport function takes velocity by const reference. No
   concentration arithmetic is reordered. References must not be retained across
   modifications of the edge's vectors.
5. **Typed mixing access.** Transport configuration stores typed member pointers
   for each species' incoming/outgoing edges and node concentration. Runtime node
   mixing, edge connection and virtual-edge boundary access no longer dispatch
   through chains of field-name string comparisons. Incoming/outgoing summation
   order, closed-edge rules, master-node rules and zero-inflow behavior are unchanged.

## Runtime and percentage differences after every accepted step

whserver, AMD Ryzen 5 3600, GCC `-O3`, `custom_output`, model
`Bathsheba_CO2_pul_2026_03_11`. Six threads, `OMP_DYNAMIC=FALSE`,
`OMP_WAIT_POLICY=ACTIVE`. Native timing uses the original output configuration
(p8 and heart_kim_lit transport, output interval 0.01 s). Each row is the median
of three sequential 0.2 s simulations. Copying, compiling and audits are excluded.

| Step | Wall time | Speedup vs baseline | Maximum field NRMS difference vs baseline | Maximum field NRMS difference vs previous retained step |
|---|---:|---:|---:|---:|
| Baseline | 5.0943 s | 1.000x | 0% | 0% |
| Junction warm start | 4.3948 s | 1.159x | 0.000003920% | 0.000003920% |
| Prepared scalar boundary laws | 3.3253 s | 1.532x | 0.000002388% | 0.000002620% |
| Checked LU | 2.9661 s | 1.717x | 0.000002388% | 0.0000000000258% |
| Transport references | 2.8970 s | 1.758x | 0.000002388% | **0%** |
| Typed mixing access | **2.4550 s** | **2.075x** | **0.000002388%** | **0%** |

These are measured medians, not predictions. The reference-getter improvement is
modest; all repeats and their ranges are included in the JSON log. The production
executable after final cleanup was byte-identical to the final benchmark executable.

Independent native 1 s runs gave 26.9962 / 26.6474 s for baseline and
12.0536 / 12.0413 s for final: medians **26.8218 vs 12.0474 s**, approximately
**2.23x**. The different ratios for 0.2 s and 1 s reflect workload variation within
the cardiac cycle and normal timing variation.

## What the percentages mean

The full spatial states were captured at 0.05, 0.1, 0.2, 0.5 and 1 s without
enabling extra hydraulic output history. Groups include all 1D pressure, area,
velocity and derived flow fields, concentrated pressures/flows, all gas-component
fields, stationary tissue O2, RBC concentration and Hb saturation. Gauge pressure
is also compared separately by subtracting the unchanged 100000 Pa reference.

For a field group:

    NRMS_percent = 100 * sqrt(sum((new-reference)^2) / sum(reference^2)).

Tables report the largest NRMS percentage among the field groups. The JSON also
includes per-field values, absolute maxima and pointwise relative percentages.
To handle zero crossings, the pointwise relative denominator has a floor of
`1e-6*max(abs(reference_field))`. Maximum-difference / field-peak percentages are
provided separately. None of these error norms is a physical conservation norm.

For the final 1 s spatial comparison, the largest pointwise relative percentage
was **0.003176%** (velocity). The maximum O2/carbon state differences were
1.257e-8 / 1.364e-6 in their respective concentration conventions. The maximum
absolute pressure difference was 0.00657 Pa. Checkpoint global-time mismatch was
at most 1.53e-14 s. Each field retains its asynchronous solver time semantics.

The full set of **136 saved histories**, including every saved point through 1 s,
was also compared using the uninstrumented native executables. Maximum per-file
NRMS difference was **0.0000198062%**, and maximum pointwise relative difference
was **0.000113063%**. The worst file was p8/HCO3_rbc/arteriole.txt. These history
comparisons use the original output precision, unlike the 17-digit spatial traces.

## Conservation and equation checks

An external audit copy independently recomputes the original scalar-junction
flow laws at the converged pressure and all concentrated hydraulic equations
at the accepted state, before time/state substitution. Any diagnostic changes
to diode flags, interpolation indices, residuals and Jacobians are restored.
The original nonlinear/linear equations are used rather than the prepared scalar
derivative. No extra node history is enabled, since that would alter the baseline
Newton starting guesses. Audited timings are not used as native benchmark timings.

| Check, across the 1 s audit | Baseline | Final |
|---|---:|---:|
| Maximum junction continuity residual | 3.91e-15 m3/s | **1.69e-17 m3/s** |
| Maximum normalized junction imbalance | 0.0000003377% | **0.0000000006531%** |
| Maximum concentrated-node continuity residual | 2.27e-13 ml/s | **2.27e-13 ml/s** |
| Maximum normalized concentrated-node imbalance | 3.27e-14% | **3.27e-14%** |
| Maximum recomputed mixed-unit hydraulic equation norm | 8.47e-9 | **8.47e-9** |

The 1 s audit checks 969537 scalar-node balances and 1613678 concentrated-node
balances. Normalized junction throughput has a 1e-6 m3/s floor; concentrated
throughput has a 1 ml/s floor. Ground pressure constraints are excluded from the
flow-continuity maximum, but included in the complete equation norm.

Every accepted stage also passed:

- 640 source steps with analytic/reduced derivative checks, physical bounds,
  original backward-Euler equation residuals and O2/CO2 stoichiometric balances;
- full capillary forward/reverse-flow balance tests with varying RBC concentration:
  maximum O2/CO2 discrete sum-balance error **3.93e-12**, relative RBC balance error
  **8.50e-18**; serial and parallel capillary states exactly match;
- independent comparison of prepared scalar flow laws with the original linear
  and Olufsen laws: maximum flow difference **6.44e-19 m3/s**; analytic derivative
  relative error against central differences **4.38e-10**;
- current node state vs optional history, all eight typed species, closed edges,
  reversed/stagnant inflow and master-node preservation checks.

## Three-second validation and compiler/memory checks

The final and baseline versions were additionally compared at 0.2, 0.5, 1, 2 and
3 s, covering several cardiac cycles. Largest field NRMS difference was
**0.000002069%**; O2 and carbon NRMS were **0.00000004756%** and
**0.00000001651%**. Maximum pointwise relative difference was **0.001228%**.
Across that run, maximum final junction imbalance was **1.73e-17 m3/s** and
concentrated-node imbalance remained **2.27e-13 ml/s**.

GCC AddressSanitizer + UndefinedBehaviorSanitizer builds passed the hydraulic,
typed-mixing, gas-source, capillary and existing Albanese tests, plus 2/5 ms
full-network smoke tests, both legacy paths and malformed-input cases. Leak
checking was disabled because existing solver ownership remains unchanged.
A separate Clang serial build passed a complete 1 s saved-history comparison;
its maximum history differences matched the GCC result at output precision.
All saved values were finite/nonnegative and Hb saturation remained in [0,1].

## Attempts deliberately not retained

| Attempt | Measured 0.2 s wall time | Outcome |
|---|---:|---|
| Skip final QR at original 1e-5 mixed-unit residual tolerance | 2.9197 s | **Rejected:** O2 NRMS 2.344%, carbon 0.764%, saturation 1.109%. |
| Skip only in smooth models at 1e-10 residual and 1e-12 ml/s continuity | 3.3171 s | Accurate, but only about 0.25% median gain over 3.3253 s; not significant relative to timing variation. |
| Additional exact characteristic reuse in concentrated Newton solves | 2.9641 s | Exactly identical spatial results to LU stage, but only 0.07% median gain over 2.9661 s; not retained. |

Consequently the final implementation retains the original nonlinear stopping
procedure and final correction. It does not contain either QR early-exit
attempt or the extra concentrated-characteristic cache. The scalar junction
preparation, which provided a clear improvement, is retained.

## Logs and reproduction

[benchmarks/hydraulic_optimization_20261119.json](benchmarks/hydraulic_optimization_20261119.json)
contains all timing repeats, baseline and previous-step percentage differences,
gauge-pressure comparisons, failed-attempt records, conservation summaries,
parameter hashes, the 136-history comparison and sanitizer logs.

Server-side source snapshots, native/audited executables, full spatial traces
and build logs are under
`/home/wr/first_blood-runs/hydro_optimization_20261119`. Audit/benchmark orchestration
scripts are also retained in the shared analysis workspace under `implementations/`.

Build from `projects/custom_output`:

```sh
make -j2 -f custom_output.mk CXX=g++
OMP_NUM_THREADS=6 OMP_DYNAMIC=FALSE OMP_WAIT_POLICY=ACTIVE ./custom_output.out Bathsheba_CO2_pul_2026_03_11
```

This command uses the original model duration; all validation runs used temporary
copies with 0.2, 1 or 3 s durations. Six threads remain recommended on whserver.
Clang can be built in serial mode with the existing makefile behavior.
New regression tests are `tests/test_hydraulic_optimization.cpp` and
`tests/test_transport_mixing.cpp`; the latter links the same core objects as
`tests/test_gas_capillary.cpp`.
