# Iteration 18 — N10 admitted: Z4c-CCM implemented in AthenaK, verified on CPU and GPU

1. **Paper anchor** — DAG node N10; arXiv:2409.10383 (medium-effort decomposition, reformulate/z4c-CMM/paper_2409.10383/); formulation rows 1–2 realized in production code.
2. **What shipped** — vendored AthenaK at pinned SHA (athenak/PROVENANCE.md; tracked in-repo by user directive); athenak/src/z4c/z4c_ccm.hpp (device-inline Kokkos CCM injection + Type-II boost + analytic psi0 provider); wiring in z4c.hpp/z4c.cpp/z4c_Sbc.cpp (12 boundary sites); test input + battery (scripts/test_athenak_ccm.sh, check_athenak_ccm.py); CPU(OpenMP) and GPU(CUDA/A100) builds + batteries ALL PASS; CPU=GPU to full hst precision on all signal samples; ledger rows (3 root-caused fails: user-BC enrollment, nvcc_wrapper sm_70 vs CUDA 13, datum-tail threshold; 2 passes + agreement).
3. **Next-3 roadmap** — CCE-coupled datum provider (replace Z4cCCMDatumPsi0; partner: src/z4c/cce/ worldtube infra); exact U⁻(E,B) Bjørhus correction in-code (N9-adjacent); puncture-background CCM run (one_puncture + ccm).
4. **Simplification flag** — header-only injection avoided CMake churn; single-function datum replacement point by design.
5. **Verifier outputs** (certified, results/numerical/athenak_ccm/):
   CPU: T1 flatness 6.78e-21 | T2 peak RMS 2.477e-2 (pre-pulse 1.1e-8 of peak = datum tail) | T3 linearity 2.005536 — OVERALL PASS
   GPU: T1 flatness 2.58e-33 | T2/T3 identical to CPU — OVERALL PASS
   Agreement: 35/35 signal samples equal to full printed precision (small-value diffs = 6-digit hst floor)
