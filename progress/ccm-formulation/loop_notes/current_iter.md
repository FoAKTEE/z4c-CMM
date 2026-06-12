# Iteration 19 — N10 final: 4-GPU MPI verification; 64-core builds

1. **Paper anchor** — N10 (AthenaK); multi-rank/multi-device execution path.
2. **What shipped** — MPI+CUDA build (Athena_ENABLE_MPI, cuda 12.8 toolkit, -j64); launcher support in the battery script; 4-GPU battery ALL PASS (UCX PML; ob1 segfault root-caused and ledgered); 1-GPU vs 4-GPU time series identical to full .hst precision; device snapshot proves one rank per A100 (gpu4_devices.txt); stray-keystroke fix in z4c_ccm.hpp line 1 (entered the previous vendor commit; non-compiling there, fixed here).
3. **Next-3 roadmap** — CCE-coupled datum provider; exact U⁻(E,B) in-code; puncture-background CCM run.
4. **Simplification flag** — n/a.
5. **Verifier outputs** (certified, results/numerical/athenak_ccm/):
   gpu4: T1 flatness 8.36e-33 | T2 peak RMS 2.477e-2, pre-pulse = datum tail | T3 ratio 2.005536 — OVERALL PASS
   gpu4_agreement: 35/35 signal samples identical to printed precision vs 1-GPU
   gpu4_devices: 4 distinct GPU UUIDs, one rank PID each (454 MiB)
