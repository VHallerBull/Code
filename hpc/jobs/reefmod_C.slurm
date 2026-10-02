#!/bin/bash
# reefmod_C.slurm -- FALLBACK ReefMod route, stage C (container, Julia): read the MATLAB outputs into a
# ReefModResultSet and write <project>/reefmod/ in the same layout 08_RUN_ENGINE.jl writes, so the
# analysis job treats it like any other engine.
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.out

set -euo pipefail
source "${CMR_BASE:-$HOME/CoralModelReview}/jobs/common.sh"
load_singularity
container julia --startup-file=no /opt/coral/Code/hpc/reefmod_fallback.jl finish
