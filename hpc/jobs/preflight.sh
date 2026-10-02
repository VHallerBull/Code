#!/bin/bash
# preflight.sh -- one-time checks on the HPC before the first run (RUNNING_ON_HPC.md section 2).
# Run on the LOGIN node from ~/CoralModelReview:   bash jobs/preflight.sh
# Nothing here changes anything; it prints what it finds and flags what needs attention.
set -uo pipefail
source "${CMR_BASE:-$HOME/CoralModelReview}/jobs/common.sh"

ok()   { echo "  [ok]   $*"; }
warn() { echo "  [WARN] $*"; }
bad()  { echo "  [FAIL] $*"; }

echo "== Layout ($CMR_BASE)"
for d in "$CMR_BASE" "$RUNS" "$CMR_BASE/logs" "$JOBS_DIR"; do
    [ -d "$d" ] && ok "$d" || bad "missing folder $d (mkdir -p it)"
done
[ -f "$SIF" ] && ok "image $(ls -lh "$SIF" | awk '{print $5}') $SIF" || bad "no image at $SIF (upload coral_engines.sif)"
if [ -f "$SIF.sha256" ] && [ -f "$SIF" ]; then
    (cd "$(dirname "$SIF")" && sha256sum -c --status "$(basename "$SIF").sha256") && ok "image checksum matches" || bad "image checksum MISMATCH (re-upload)"
fi

echo "== Singularity"
module load "$SINGULARITY_MODULE" 2>/dev/null || true
command -v singularity >/dev/null && ok "$(singularity --version)" || bad "no singularity (module $SINGULARITY_MODULE)"
if [ -f "$SIF" ] && command -v singularity >/dev/null; then
    singularity exec "$SIF" julia --version >/dev/null 2>&1 && ok "image runs: $(singularity exec "$SIF" julia --version)" || bad "image does not start"
fi

echo "== Slurm"
command -v sbatch >/dev/null && ok "sbatch available" || bad "no sbatch"
command -v sacct >/dev/null && sacct -n -X -S "$(date -d '-1 day' +%F)" -u "$USER" >/dev/null 2>&1 && ok "sacct works (benchmarks can use Slurm accounting)" || warn "sacct not usable: benchmarks fall back to in-process numbers"
sinfo -p defq -h -o "%P %a %l %D %c %m" 2>/dev/null | sed 's/^/  defq: /' || warn "no partition defq"
[ -x /usr/bin/time ] && ok "/usr/bin/time present" || warn "/usr/bin/time missing: per-job time logs are skipped"

echo "== MATLAB (ReefMod)"
if module avail "$MATLAB_MODULE" 2>&1 | grep -qi matlab; then
    module avail "$MATLAB_MODULE" 2>&1 | grep -i matlab | sed 's/^/  available: /'
    module load "$MATLAB_MODULE"
    exe="$(readlink -f "$(command -v matlab)")"
    ok "matlab: $exe"
    echo "  licence variables: MLM_LICENSE_FILE=${MLM_LICENSE_FILE:-<unset>} LM_LICENSE_FILE=${LM_LICENSE_FILE:-<unset>}"
    echo "  asking MATLAB for its version and licences (can take a minute)..."
    matlab -nodisplay -nosplash -batch "fprintf('version %s\n', version); fprintf('Statistics toolbox licensed: %d\n', license('test','Statistics_Toolbox')); fprintf('MATLAB Compiler licensed: %d\n', license('test','Compiler'))" 2>&1 | sed 's/^/  matlab> /'
    echo "  ReefMod-GBR's README asks for MATLAB R2023b or earlier; set MATLAB_MODULE=matlab/<version> to pick one."
    command -v lmutil >/dev/null && lmutil lmstat -a 2>/dev/null | grep -iE "Users of MATLAB|Statistics" | sed 's/^/  licences: /'
else
    bad "no '$MATLAB_MODULE' module: ReefMod cannot run (set MATLAB_MODULE to the right name)"
fi

echo
echo "Next: the ReefMod check on a compute node:"
echo "  srun -p defq -c 1 --mem=16G -t 01:00:00 --pty bash $JOBS_DIR/reefmod_check.sh"
