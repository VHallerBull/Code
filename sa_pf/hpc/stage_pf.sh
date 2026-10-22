#!/usr/bin/env bash
# stage_pf.sh MODEL -- assemble the build folder of the model's pick-and-freeze image (run in WSL).
#
#   bash /mnt/c/Users/vhaller/DecisionSupportProjects/CoralModelReview/Code/sa_pf/hpc/stage_pf.sh coralpop
#
# The same image as Code/sa/hpc/stage_model.sh builds (same definition, packages and repos; that script
# does the staging), plus Code/sa_pf in /opt/sa/Code/sa_pf and the pick-and-freeze job scripts. Staged
# into its own folder (default ~/sa-build/MODEL_pf) so the Sobol' image's build folder stays as it is;
# upload the result as images/MODEL_pf.sif. Then:
#
#   BUILD_DIR=~/sa-build/coralpop_pf bash ~/sa-build/coralpop_pf/build_model_image.sh coralpop
set -euo pipefail

MODEL="${1:?usage: stage_pf.sh MODEL (e.g. coralpop)}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PF_SRC="$(cd "$HERE/.." && pwd)"
CODE_SRC="$(cd "$PF_SRC/.." && pwd)"
export BUILD_DIR="${BUILD_DIR:-$HOME/sa-build/${MODEL}_pf}"

bash "$CODE_SRC/sa/hpc/stage_model.sh" "$MODEL"     # repos, Code/sa, definition, Sobol' job scripts

echo "  Code/sa_pf (part of the CoralModelReview-Code repo recorded above)"
rsync -a --delete --exclude=.claude/ --exclude=__pycache__/ "$PF_SRC/" "$BUILD_DIR/sa/Code/sa_pf/"
rsync -a "$HERE/jobs/" "$BUILD_DIR/upload/jobs/"
find "$BUILD_DIR/upload/jobs" "$BUILD_DIR/sa/Code/sa_pf/hpc" -type f \( -name '*.sh' -o -name '*.slurm' \) \
    -exec sed -i 's/\r$//' {} +
chmod +x "$BUILD_DIR/upload/jobs/"*.sh "$BUILD_DIR/upload/jobs/"*.slurm
echo
echo "Next: BUILD_DIR=$BUILD_DIR bash $BUILD_DIR/build_model_image.sh $MODEL"
echo "      then upload $BUILD_DIR/${MODEL}_sa.sif as images/${MODEL}_pf.sif (Code/sa_pf/hpc/RUNNING_PF_ON_HPC.md)"
