#!/bin/bash
# reefmod_A.slurm -- FALLBACK ReefMod route, stage A (container, Julia): sample + write MATLAB inputs.
# Only used with `submit_all.sh PROJECT --reefmod-fallback` (MATLAB cannot start inside the container).
# Stage B (reefmod_B.slurm, host MATLAB, one array task per scenario) and C (reefmod_C.slurm) follow.
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.out

set -euo pipefail
source "${CMR_BASE:-$HOME/CoralModelReview}/jobs/common.sh"
load_singularity
mkdir -p "$RUNS/$CMR_PROJECT_NAME"
container julia --startup-file=no /opt/coral/Code/hpc/reefmod_fallback.jl prepare
