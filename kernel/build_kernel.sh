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
ISO_VARIANT=${ISO_VARIANT:-v1}
CARMA=$SCRAM_SRC/components/cam/src/physics/carma
CORE=$CARMA/models/scram/scram_core
B=$HERE/build/$ISO_VARIANT
FC=${FC:-gfortran}

mkdir -p $B/src $B/mod $B/obj
rm -f $B/obj/*.o $B/mod/*.mod $B/libscramcore.a

# ---- snapshot of the kernel sources
for d in . scram-1.1 AtmoData RDB VSRM SIMPLE_AQUEOUS INC isorropia; do
  mkdir -p $B/src/$d
  find $CORE/$d -maxdepth 1 -type f \( -name '*.F90' -o -name '*.f90' -o -name '*.f' -o -name '*.inc' -o -name '*.INC' \) \
       -exec cp {} $B/src/$d/ \;
done
# test_ppm_dualmoment.f90 is a standalone program; test_mass_nb.f90 is a library routine and stays
rm -f $B/src/RDB/test_ppm_dualmoment.f90
cp $CARMA/models/scram/scram_cfg_io_mod.F90 $B/src/
cp $CARMA/base/carma_precision_mod.F90 $B/src/
cp $HERE/stubs/*.F90 $B/src/
( cd $SCRAM_SRC/components/cam && echo "cam $(git log -1 --format='%H %ad %s' --date=short)" ) > $B/KERNEL_VERSION
echo "source $SCRAM_SRC" >> $B/KERNEL_VERSION
echo "isorropia $ISO_VARIANT" >> $B/KERNEL_VERSION

DEFS=""
if [ "$ISO_VARIANT" = "v23" ]; then
  DEFS="-DUSE_ISORROPIA2"
  mkdir -p $B/src/isorropia2
  cp $HERE/isorropia2/*.f $HERE/isorropia2/*.inc $HERE/isorropia2/*.F90 $B/src/isorropia2/ 2>/dev/null
fi

INCS="-I$B/src/INC -I$B/src/scram-1.1 -I$B/src/isorropia -I$B/src/RDB -I$B/src/VSRM -I$B/src/SIMPLE_AQUEOUS -I$B/src"
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
pending=$(find $B/src -name '*.F90' -o -name '*.f90' | grep -v '/isorropia/' )
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
cat $B/KERNEL_VERSION
