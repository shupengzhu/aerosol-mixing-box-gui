#!/bin/bash
# Build a self-contained SCRAM Box package for a single Linux x86-64 machine.
#
#   bash tools/package_standalone.sh            -> dist/scram-box-<date>.tar.gz (v1 / no-ISO kernel)
#   INCLUDE_V23=1 bash tools/package_standalone.sh   also ships the ISORROPIA II v2.3 kernel
#
# The kernels are relinked fully static (gfortran runtime, OpenMP, BLAS), so the target
# needs only Python 3.10+ with numpy and netCDF4 (pip install -r requirements.txt).
# ISORROPIA II is licensed for research use and may not be redistributed: build with
# INCLUDE_V23=1 only for machines covered by that licence (the user's own group).
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
FC=${FC:-gfortran}
STAMP=$(date +%Y%m%d)
NAME=scram-box-$STAMP
OUT=$ROOT/dist/$NAME
rm -rf "$OUT"; mkdir -p "$OUT/kernel/build" "$OUT/kernel/data"

variants="v1"
[ "${INCLUDE_V23:-0}" = "1" ] && variants="v1 v23"
for v in $variants; do
  B=$ROOT/kernel/build/$v
  [ -f "$B/libscramcore.a" ] || { echo "missing $B/libscramcore.a - run kernel/build_kernel.sh"; exit 1; }
  mkdir -p "$OUT/kernel/build/$v" "$OUT/tmp_obj_$v"
  $FC -O2 -fPIC -cpp -fopenmp -ffree-line-length-none -fallow-argument-mismatch \
      $([ "$v" = v23 ] && echo -DUSE_ISORROPIA2) -I"$B/mod" -J "$OUT/tmp_obj_$v" \
      "$ROOT/kernel/driver/scram_box.F90" "$B/libscramcore.a" -static -lblas -lm \
      -o "$OUT/kernel/build/$v/scram_box" 2> >(grep -v "dlopen\|gomp_target_init" >&2)
  strip "$OUT/kernel/build/$v/scram_box"
  rm -rf "$OUT/tmp_obj_$v"
  cp "$B/KERNEL_VERSION" "$OUT/kernel/build/$v/"
  echo "== static kernel $v: $(ls -la "$OUT/kernel/build/$v/scram_box" | awk '{print $5}') bytes"
done

cp -r "$ROOT/scrambox" "$OUT/"
find "$OUT/scrambox" -name __pycache__ -prune -exec rm -rf {} +
cp -r "$ROOT/scenarios" "$ROOT/docs" "$OUT/"
cp "$ROOT"/kernel/data/*.dat "$ROOT"/kernel/data/*.json "$OUT/kernel/data/"
mkdir -p "$OUT/kernel/data/lut"
cp "$ROOT"/kernel/data/lut/*.nc "$OUT/kernel/data/lut/" 2>/dev/null || echo "warning: no LUTs (run kernel/fetch_luts.sh); optics will be unavailable"
cp "$ROOT/README.md" "$OUT/"
printf 'numpy>=1.24\nnetCDF4>=1.6\n' > "$OUT/requirements.txt"
cat > "$OUT/scram-box" <<'EOF'
#!/bin/bash
# Start the SCRAM Box web interface:  ./scram-box [--port 8765] [--host 127.0.0.1]
HERE=$(cd "$(dirname "$0")" && pwd)
export SCRAMBOX_USER_DIR=${SCRAMBOX_USER_DIR:-$HOME/.scrambox}
cd "$HERE" && exec python3 -m scrambox serve "$@"
EOF
chmod +x "$OUT/scram-box"

# smoke test of the package itself (no repository paths)
( cd "$OUT" && SCRAMBOX_USER_DIR="$OUT/.selftest" python3 -m scrambox run scenarios/tutorial_minimal.json -o "$OUT/.selftest/run" > "$OUT/.selftest.log" 2>&1 ) \
  && echo "== package self-test (v1/no-ISO): OK" || { echo "== package self-test FAILED"; cat "$OUT/.selftest.log"; exit 1; }
if [ "${INCLUDE_V23:-0}" = "1" ]; then
  ( cd "$OUT" && SCRAMBOX_USER_DIR="$OUT/.selftest" python3 -m scrambox run scenarios/seasalt_nitrate.json -o "$OUT/.selftest/run23" > "$OUT/.selftest.log" 2>&1 ) \
    && echo "== package self-test (v2.3): OK" || { echo "== v2.3 self-test FAILED"; cat "$OUT/.selftest.log"; exit 1; }
fi
rm -rf "$OUT/.selftest" "$OUT/.selftest.log"
tar -C "$ROOT/dist" -czf "$ROOT/dist/$NAME.tar.gz" "$NAME"
echo "== $(ls -la "$ROOT/dist/$NAME.tar.gz" | awk '{print $5}') bytes: dist/$NAME.tar.gz (kernels: $variants)"
