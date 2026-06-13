#!/usr/bin/env bash
# S0 (mission 6): SpECTRE SolveXcts runtime-MKL fix.
# The build's DT_NEEDED pulls the unversioned system libmkl_gf_lp64.so
# (/lib/x86_64-linux-gnu) which lacks mkl_lapack_zhesvxx. The consistent
# MKL set lives in the shim /data/jiaxiwu/spectre_mkl_shim (symlinks to
# /home/jiaxiwu/miniconda3/lib MKL 2.x). Fix: resolve the NEEDED interface
# lib from the shim (LD_LIBRARY_PATH) AND put the full layered set in the
# global symbol scope (LD_PRELOAD) so MKL's lazy cross-references resolve.
# Source this, then run: SolveXcts --input-file <yaml>
#
# S1 addendum (mission 6): libmkl_gnu_thread.so calls omp_get_num_procs from the
# GNU OpenMP runtime, which is NOT pulled by the MKL libs themselves. Under a
# real solve (threaded MKL BLAS in the Gmres/direct linear solve) this aborts
# with "symbol lookup error: ... undefined symbol: omp_get_num_procs" (exit 127)
# even though --help/--check-options succeed (they never call threaded MKL).
# Fix: preload libgomp.so.1 FIRST so the OpenMP symbols are in the global scope.
# Pin OMP/MKL threads to 1 so the N MPI ranks (charmrun +pN) do not each spawn a
# full OpenMP team and oversubscribe the cores.
export SPECTRE_MKL_SHIM=/data/jiaxiwu/spectre_mkl_shim
export LD_LIBRARY_PATH="$SPECTRE_MKL_SHIM:$LD_LIBRARY_PATH"
export LD_PRELOAD="$SPECTRE_MKL_SHIM/libgomp.so.1 $SPECTRE_MKL_SHIM/libmkl_gf_lp64.so $SPECTRE_MKL_SHIM/libmkl_gnu_thread.so $SPECTRE_MKL_SHIM/libmkl_core.so $SPECTRE_MKL_SHIM/libmkl_def.so"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export SPECTRE_BUILD=/data/haiyangw/nr/spectre/build
