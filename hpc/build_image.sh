#!/usr/bin/env bash
# build_image.sh -- build coral_engines.sif from the folder stage_sources.sh assembled (run in WSL).
#
#   bash ".../Code/hpc/build_image.sh"            # after stage_sources.sh
#
# Builds unprivileged (--fakeroot; needs an /etc/subuid entry for you, which Ubuntu creates by default),
# so no sudo. Takes ~50 min (47 min on this laptop), mostly Julia precompilation; the image is ~1.8 GB. Then
# runs the image's %test (all six engine packages load, Python found) and prints the upload commands.
set -euo pipefail

BUILD_DIR="${BUILD_DIR:-$HOME/coral-build}"
SIF="$BUILD_DIR/coral_engines.sif"
HPC="${HPC:-vhaller@hpc-l001.aims.gov.au}"  # ssh login of the AIMS cluster; override as user@host
HPC_BASE="${HPC_BASE:-CoralModelReview}"    # relative to your HPC home

command -v apptainer >/dev/null || {
    echo "Apptainer is not installed in this WSL distribution. Install it with:"
    echo "  sudo add-apt-repository -y ppa:apptainer/ppa && sudo apt-get update && sudo apt-get install -y apptainer"
    exit 1
}
[ -f "$BUILD_DIR/coral_engines.def" ] && [ -d "$BUILD_DIR/coral/Code" ] || {
    echo "Nothing staged in $BUILD_DIR -- run stage_sources.sh first." >&2; exit 1; }

cd "$BUILD_DIR"
# Apptainer's own cache/tmp default to ~/.apptainer and /tmp; keep both on the big ext4 disk.
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$BUILD_DIR/.apptainer-cache}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-$BUILD_DIR/.apptainer-tmp}"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

echo "Building $SIF (log: $BUILD_DIR/build.log) -- started $(date '+%H:%M')"
start=$(date +%s)
apptainer build --fakeroot --force "$SIF" coral_engines.def 2>&1 | tee "$BUILD_DIR/build.log"
echo "Build took $(( ($(date +%s) - start) / 60 )) min"

echo "Testing the image"
apptainer test "$SIF"

# File name only (no directory), so `sha256sum -c` also works next to the uploaded copy on the HPC.
(cd "$(dirname "$SIF")" && sha256sum "$(basename "$SIF")") | tee "$SIF.sha256"
cp -f "$BUILD_DIR/coral/STAGED_SOURCES.md" "$BUILD_DIR/upload/STAGED_SOURCES.md"
ls -lh "$SIF"

cat <<EOF

Image ready. Upload to the HPC (see RUNNING_ON_HPC.md, section 5):

  ssh $HPC "mkdir -p ~/$HPC_BASE/{runs,logs,jobs}"
  rsync -avP $SIF $SIF.sha256 $HPC:~/$HPC_BASE/
  rsync -avP --chmod=go-w $BUILD_DIR/upload/jobs/ $HPC:~/$HPC_BASE/jobs/
  rsync -avP $BUILD_DIR/upload/STAGED_SOURCES.md $HPC:~/$HPC_BASE/

Only for the fallback ReefMod route (host MATLAB outside the container):
  rsync -avP --chmod=go-w $BUILD_DIR/upload/reefmod_host/ $HPC:~/$HPC_BASE/reefmod_host/
EOF
