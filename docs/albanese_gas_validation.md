# Short validation of dev_gas, 2026-11-17

Base revision: 4916dd19cac70134bbd9ffa64d62989d852874b2.
Checks ran on whserver (Ubuntu, Eigen 3), using the `custom_output` project.

- Optimized Clang and GCC builds.
- Standalone gas-curve tests: units/reference value, Bohr and Haldane response
  directions, compatibility with the existing RBC oxygen capacity, reconstruction
  of total carbon content, cancellation of internal reaction sources, and
  rejection of malformed/nonfinite model parameters.
- Full network smoke integrations in legacy and Albanese modes, each **0.05 s**
  of simulated time. The network is copied to disposable model directories;
  the original model's 130 s duration is not changed or executed.
- Each mode produces 136 transport histories (p8 and pulmonary circulation).
  Checks require finite, nonnegative gas values, saturation in [0,1], constant
  RBC concentration, initial t=0 records and increasing timestamps.
- A further 0.05 s integration confirms that selecting the model through
  `CO2_parameters.txt` gives the same response as through `O2_parameters.txt`.
  The new mode must change Hb reaction response relative to legacy mode.
- Malformed dissociation parameters and a missing `perifco2` configuration must
  fail with an explanatory message and nonzero exit status.
- GCC AddressSanitizer and UndefinedBehaviorSanitizer builds run the same short
  checks with `ASAN_OPTIONS=detect_leaks=0` and
  `UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1`. Leak detection is excluded
  because the existing solver has unmanaged allocations outside this change.
  Clang's sanitizer runtime is not installed on this server, so GCC provides
  the sanitizer check. The final executable is restored to optimized Clang.

The checks exposed and fixed out-of-bounds optional CSV column reads and
use of uninitialized cardiac phase/rate during the first elastance evaluation.
The previously observed transport history/time mismatch is also fixed.

Existing compiler warnings about incomplete transport-enum switches remain.
These checks establish short-run numerical/code health, not long-time stability,
physiological calibration, convergence, or full reproduction of the Albanese
cardiopulmonary model. No long simulation was performed.
