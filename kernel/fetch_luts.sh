#!/bin/bash
# Copy the SCRAM optics look-up tables used by the post-processor (scrambox/post/optics.py)
# into kernel/data/lut. The tables are not part of this repository:
#   reduced-shell (Optics_mixing 2): scram_reducedshell_optics_lut_cdf5.nc
#   core-shell    (Optics_mixing 3): scram_coreshell_optics_lut.nc
# Usage: LUT_SRC=<directory holding the two tables> bash kernel/fetch_luts.sh
set -eu
HERE=$(cd "$(dirname "$0")" && pwd)
LUT_SRC=${LUT_SRC:?set LUT_SRC to the directory holding the optics tables}
mkdir -p "$HERE/data/lut"
cp -v "$LUT_SRC/scram_reducedshell_optics_lut_cdf5.nc" "$LUT_SRC/scram_coreshell_optics_lut.nc" "$HERE/data/lut/"
( cd "$HERE/data/lut" && sha256sum *.nc > SHA256SUMS && echo "source: $LUT_SRC" > SOURCE.txt )
