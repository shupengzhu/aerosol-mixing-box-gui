#!/bin/bash
# Build the CESM-SCRAM hetfix aerosol kernel as a standalone static library (libscramcore.a).
#
#   SCRAM_SRC   CESM-SCRAM source tree (default: hetfix release on ts-106)
#   ISO_VARIANT v1 (SCRAM-revised ISORROPIA v1, default) | v23 (adds -DUSE_ISORROPIA2)
#
# The kernel sources are snapshot-copied into build/<variant>/src so the box app always
# records exactly which CESM-SCRAM revision it was built from (build/<variant>/KERNEL_VERSION).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
SCRAM_SRC=${SCRAM_SRC:-/data/shupengz/CESM_SCRAM/releases/cesm-scram-hetfix-20261005-r3/source}
# CAM_SRC: the CAM component checkout (defaults to the one inside SCRAM_SRC); point it at a
# CAM worktree to build the box kernel from a branch.
CAM_SRC=${CAM_SRC:-$SCRAM_SRC/components/cam}
ISO_VARIANT=${ISO_VARIANT:-v1}
CARMA=$CAM_SRC/src/physics/carma
CORE=$CARMA/models/scram/scram_core
B=$HERE/build/$ISO_VARIANT
FC=${FC:-gfortran}

mkdir -p $B/src $B/mod $B/obj
rm -f $B/obj/*.o $B/mod/*.mod $B/libscramcore.a

# ---- snapshot of the kernel sources; like CAM's bld/configure, exactly one ISORROPIA
# directory is compiled: isorropia (v1) or isorropia2 (v2.3, -DUSE_ISORROPIA2)
if [ "$ISO_VARIANT" = "v23" ]; then ISODIR=isorropia2; else ISODIR=isorropia; fi
if [ "$ISODIR" = "isorropia2" ] && [ ! -f $CORE/isorropia2/isocom2.f ]; then
  echo "isorropia2 sources missing: run python3 $CORE/isorropia2/prepare_isorropia2.py <iso_2_3Src>"; exit 1
fi
for d in . scram-1.1 AtmoData RDB VSRM SIMPLE_AQUEOUS INC $ISODIR; do
  mkdir -p $B/src/$d
  find $CORE/$d -maxdepth 1 -type f \( -name '*.F90' -o -name '*.f90' -o -name '*.f' -o -name '*.inc' -o -name '*.INC' \) \
       -exec cp {} $B/src/$d/ \;
done
# test_ppm_dualmoment.f90 is a standalone program; test_mass_nb.f90 is a library routine and stays
rm -f $B/src/RDB/test_ppm_dualmoment.f90
cp $CARMA/models/scram/scram_cfg_io_mod.F90 $B/src/
cp $CARMA/base/carma_precision_mod.F90 $B/src/
cp $HERE/stubs/*.F90 $B/src/
( cd $CAM_SRC && echo "cam $(git rev-parse --abbrev-ref HEAD) $(git log -1 --format='%H %ad %s' --date=short)"; \
  [ -n "$(git status --porcelain --untracked-files=no)" ] && echo "WARNING: uncommitted changes in $CAM_SRC" ) > $B/KERNEL_VERSION
echo "source $CAM_SRC" >> $B/KERNEL_VERSION
echo "isorropia $ISO_VARIANT ($ISODIR)" >> $B/KERNEL_VERSION
[ -f $B/src/$ISODIR/ISORROPIA2_MANIFEST.txt ] || cp $CORE/$ISODIR/ISORROPIA2_MANIFEST.txt $B/src/$ISODIR/ 2>/dev/null
[ -f $B/src/$ISODIR/ISORROPIA2_MANIFEST.txt ] && sed 's/^/  /' $B/src/$ISODIR/ISORROPIA2_MANIFEST.txt >> $B/KERNEL_VERSION

DEFS=""
[ "$ISO_VARIANT" = "v23" ] && DEFS="-DUSE_ISORROPIA2"

INCS="-I$B/src/INC -I$B/src/scram-1.1 -I$B/src/$ISODIR -I$B/src/RDB -I$B/src/VSRM -I$B/src/SIMPLE_AQUEOUS -I$B/src"
# -fopenmp also for F77: SCRAM declares ISORROPIA common blocks THREADPRIVATE, so every
# object must agree on thread-local storage (CESM compiles everything with OpenMP).
FFLAGS_FIX="-O2 -fPIC -w -std=legacy -fopenmp -fallow-argument-mismatch -ffixed-line-length-none"
FFLAGS_FREE="-O2 -fPIC -cpp -fopenmp -ffree-line-length-none -fallow-argument-mismatch $DEFS"

# ---- fixed-form Fortran 77 (no modules)
nfail=0
for f in $(find $B/src -name '*.f'); do
  o=$B/obj/$(basename ${f%.*})_f77.o
  $FC $FFLAGS_FIX $INCS -J $B/mod -c $f -o $o 2>> $B/compile_f77.log || { echo "F77 FAIL $f"; nfail=$((nfail+1)); }
done

# ---- free-form sources: iterate until every module dependency is satisfied
pending=$(find $B/src -name '*.F90' -o -name '*.f90')
for pass in 1 2 3 4 5 6 7 8 9 10 11 12; do
  next=""
  for f in $pending; do
    o=$B/obj/$(basename ${f%.*}).o
    if $FC $FFLAGS_FREE $INCS -J $B/mod -I$B/mod -c $f -o $o > $B/last_err.log 2>&1; then :; else next="$next $f"; cp $B/last_err.log $B/err_$(basename $f).log; fi
  done
  [ -z "$next" ] && break
  [ "$next" = " $pending" ] && break
  pending=$next
done

if [ -n "$next" ]; then
  echo "== FREE-FORM FILES THAT STILL FAIL:"
  for f in $next; do echo "--- $f"; grep -m5 -E 'Error|error' $B/err_$(basename $f).log; done
fi
ar rcs $B/libscramcore.a $B/obj/*.o
echo "== objects: $(ls $B/obj | wc -l)   f77 failures: $nfail   free-form failures: $(echo $next | wc -w)"

# ---- whole-archive link check: every symbol must resolve with BLAS + OpenMP + libm
printf 'program link_check\n  print *, "kernel links"\nend program link_check\n' > $B/link_check.f90
if $FC -fopenmp $B/link_check.f90 -Wl,--whole-archive $B/libscramcore.a -Wl,--no-whole-archive \
      -lblas -lm -o $B/link_check > $B/link_check.log 2>&1; then
  echo "== link check: OK ($($B/link_check))"
else
  echo "== link check: FAILED"; grep -E 'undefined reference' $B/link_check.log | sed 's/.*undefined reference to //' | sort -u | head -20
fi

# ---- box driver executable
if $FC $FFLAGS_FREE -I$B/mod -J $B/obj $HERE/driver/scram_box.F90 $B/libscramcore.a \
      -lblas -lm -o $B/scram_box > $B/driver_build.log 2>&1; then
  echo "== driver: $B/scram_box"
else
  echo "== driver build FAILED"; grep -m10 -E 'Error|error' $B/driver_build.log
fi
cat $B/KERNEL_VERSION
