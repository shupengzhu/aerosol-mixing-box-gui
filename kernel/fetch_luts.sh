#!/bin/bash
# Copy the CESM-SCRAM optics look-up tables used by the post-processor (scrambox/post/optics.py)
# into kernel/data/lut. Same files as the CESM runtime:
#   reduced-shell (Optics_mixing 2, production): scram/install_package/runtime/scram_reducedshell_optics_lut_cdf5.nc
#   core-shell    (Optics_mixing 3):             components/cam/tools/lut_output/scram_coreshell_optics_lut.nc
# Usage: SCRAM_SRC=<CESM-SCRAM source tree> bash kernel/fetch_luts.sh
set -eu
HERE=$(cd "$(dirname "$0")" && pwd)
SCRAM_SRC=${SCRAM_SRC:-/data/shupengz/CESM_SCRAM/releases/cesm-scram-hetfix-20261005-r3/source}
mkdir -p "$HERE/data/lut"
cp -v "$SCRAM_SRC/scram/install_package/runtime/scram_reducedshell_optics_lut_cdf5.nc" "$HERE/data/lut/"
cp -v "$SCRAM_SRC/components/cam/tools/lut_output/scram_coreshell_optics_lut.nc" "$HERE/data/lut/"
( cd "$HERE/data/lut" && sha256sum *.nc > SHA256SUMS && echo "source: $SCRAM_SRC" > SOURCE.txt )
