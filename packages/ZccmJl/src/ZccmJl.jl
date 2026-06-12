"""
ZccmJl — Julia migration of the z4c-CCM verification toolchain (mission-3
user directive, iter 29): one language for GPU (CUDA.jl, added when the
solver lands), symbolic/exact derivations (Rational{BigInt} series algebra;
Symbolics.jl when needed), and arbitrary precision (BigFloat).

Migration contract: each ported component carries unit tests that reproduce
the admitted Python-verified numbers EXACTLY (provenance baseline:
scripts/*.py + results/numerical/ transcripts). Further development is
Julia-only once the component's tests pass.

Ported components (v0.1.0):
  teukolsky.jl    — F^(n) Hermite derivatives + the exact finite-radius
                    rpsi4(2,0) waveform (N12 table, AthenaK convention).
  linseries.jl    — exact retarded-series algebra over Rational{BigInt}
                    in the basis F^(n)(u)/r^k (the A1/A2 machinery).
  psi0chain.jl    — the stage-A1 linear psi0 chain (J_full -> psi0_lin)
                    with the C1-C4 gates as unit tests.
"""
module ZccmJl

include("teukolsky.jl")
include("linseries.jl")
include("psi0chain.jl")
include("bonditransform.jl")
include("hierarchy.jl")
include("bondisolver.jl")
include("bondisolver_reg.jl")

end # module
