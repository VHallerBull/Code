#!/usr/bin/env bash
# build_model_image.sh MODEL -- build MODEL_sa.sif from the folder stage_model.sh assembled (run in WSL).
#
#   bash ~/sa-build/coralpop/build_model_image.sh coralpop      # after stage_model.sh coralpop
#
# Unprivileged (--fakeroot), no sudo. Mostly Julia precompilation. Runs the image's %test, writes a
# sha256 and prints the upload commands.
set -euo pipefail

MODEL="${1:?usage: build_model_image.sh MODEL}"
BUILD_DIR="${BUILD_DIR:-$HOME/sa-build/$MODEL}"
SIF="$BUILD_DIR/${MODEL}_sa.sif"
HPC="${HPC:-vhaller@hpc-l001.aims.gov.au}"
SA_BASE_HPC="${SA_BASE_HPC:-CoralModelReview/sa}"   # relative to your HPC home

command -v apptainer >/dev/null || { echo "Apptainer is not installed in this WSL distribution." >&2; exit 1; }
[ -f "$BUILD_DIR/${MODEL}_sa.def" ] && [ -d "$BUILD_DIR/sa/Code/sa" ] || {
    echo "Nothing staged in $BUILD_DIR -- run stage_model.sh $MODEL first." >&2; exit 1; }

cd "$BUILD_DIR"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$HOME/sa-build/.apptainer-cache}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-$HOME/sa-build/.apptainer-tmp}"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

echo "Building $SIF (log: $BUILD_DIR/build.log) -- started $(date '+%H:%M')"
start=$(date +%s)
apptainer build --fakeroot --force "$SIF" "${MODEL}_sa.def" 2>&1 | tee "$BUILD_DIR/build.log"
echo "Build took $(( ($(date +%s) - start) / 60 )) min"

echo "Testing the image"
apptainer test "$SIF"
(cd "$(dirname "$SIF")" && sha256sum "$(basename "$SIF")") | tee "$SIF.sha256"
cp -f "$BUILD_DIR/sa/STAGED_SOURCES.md" "$BUILD_DIR/upload/STAGED_SOURCES_${MODEL}.md"
ls -lh "$SIF"

cat <<EOF

Image ready. Upload to the HPC (Code/sa/hpc/RUNNING_SA_ON_HPC.md, step 4):

  ssh $HPC "mkdir -p ~/$SA_BASE_HPC/{images,jobs,domains,logs}"
  rsync -avP $SIF $SIF.sha256 $HPC:~/$SA_BASE_HPC/images/
  rsync -avP --chmod=go-w $BUILD_DIR/upload/jobs/ $HPC:~/$SA_BASE_HPC/jobs/
  rsync -avP $BUILD_DIR/upload/STAGED_SOURCES_${MODEL}.md $HPC:~/$SA_BASE_HPC/images/
EOF
