#!/usr/bin/env bash
# S0 (mission 6): SpECTRE SolveXcts runtime-MKL fix.
# The build's DT_NEEDED pulls the unversioned system libmkl_gf_lp64.so
# (/lib/x86_64-linux-gnu) which lacks mkl_lapack_zhesvxx. The consistent
# MKL set lives in the shim /data/jiaxiwu/spectre_mkl_shim (symlinks to
# /home/jiaxiwu/miniconda3/lib MKL 2.x). Fix: resolve the NEEDED interface
# lib from the shim (LD_LIBRARY_PATH) AND put the full layered set in the
# global symbol scope (LD_PRELOAD) so MKL's lazy cross-references resolve.
# Source this, then run: SolveXcts --input-file <yaml>
export SPECTRE_MKL_SHIM=/data/jiaxiwu/spectre_mkl_shim
export LD_LIBRARY_PATH="$SPECTRE_MKL_SHIM:$LD_LIBRARY_PATH"
export LD_PRELOAD="$SPECTRE_MKL_SHIM/libmkl_gf_lp64.so $SPECTRE_MKL_SHIM/libmkl_gnu_thread.so $SPECTRE_MKL_SHIM/libmkl_core.so $SPECTRE_MKL_SHIM/libmkl_def.so"
export SPECTRE_BUILD=/data/haiyangw/nr/spectre/build
