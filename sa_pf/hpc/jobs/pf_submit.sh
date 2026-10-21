#!/bin/bash
# pf_submit.sh -- submit a model's pick-and-freeze walks for one band as a Slurm array (both design sets
# and all chosen domains in each task), then 04_collect.jl once every task has succeeded.
#
#   bash ~/CoralModelReview/sa/jobs/pf_submit.sh --model coralpop --band 0.5 --domains SingleReef_SobolDHW
#
# Options:
#   --model M          coralpop (needs images/M_pf.sif and the band's designs, 02_design.jl)
#   --band B           input band of the design, e.g. 0.5 (folder band50)
#   --domains "A B"    default: SingleReef_SobolDHW
#   --sets "S ..."     main (varied DHW) and/or zero (no heat stress); default "main zero"
#   --walks K          base samples per array task (default 50: about 20 min per task on SingleReef with
#                      both sets; use about 15 on NineReefs, which is 3-4x slower per run)
#   --first F          first base sample (default 1); --last L (default: the design's n_base). For an
#                      extended design run only the new range, e.g. --first 2501
#   --time HH:MM:SS    wall time per task (default 01:30:00)
#   --no-collect       submit the walks only
#   --dry-run          print the sbatch commands only
#
# Writes <band dir>/jobs.txt (name, job id). Slurm logs go to ~/CoralModelReview/sa/logs.
set -euo pipefail
export SA_BASE="${SA_BASE:-$HOME/CoralModelReview/sa}"

MODEL="" BAND="" DOMAINS="SingleReef_SobolDHW" SETS="main zero" WALKS=50 FIRST=1 LAST="" TIME=01:30:00 COLLECT=1 DRY=0
while [ $# -gt 0 ]; do
    case "$1" in
        --model) MODEL="$2"; shift 2 ;;
        --band) BAND="$2"; shift 2 ;;
        --domains) DOMAINS="$2"; shift 2 ;;
        --sets) SETS="$2"; shift 2 ;;
        --walks) WALKS="$2"; shift 2 ;;
        --first) FIRST="$2"; shift 2 ;;
        --last) LAST="$2"; shift 2 ;;
        --time) TIME="$2"; shift 2 ;;
        --no-collect) COLLECT=0; shift ;;
        --dry-run) DRY=1; shift ;;
        *) echo "unknown option $1" >&2; exit 1 ;;
    esac
done
if [ -z "$MODEL" ] || [ -z "$BAND" ]; then
    echo "usage: pf_submit.sh --model M --band B [--domains \"...\"] [--sets \"main zero\"] [--walks K] [--first F] [--last L] [--dry-run]" >&2
    exit 1
fi
export MODEL
source "$SA_BASE/jobs/pf_common.sh"
[ -f "$SIF" ] || { echo "ERROR: image not found: $SIF" >&2; exit 1; }
bdir="$(pf_band_dir "$BAND")"

N=""
JOBS=""
for s in $SETS; do
    sfx=""
    [ "$s" = main ] || sfx="_$s"
    toml="$bdir/design_pf$sfx.toml"
    if [ ! -f "$toml" ] || [ ! -f "$bdir/design_pf$sfx.jld2" ]; then
        echo "ERROR: $bdir/design_pf$sfx.jld2 / .toml not found -- upload the designs (02_design.jl)" >&2
        exit 1
    fi
    n=$(sed -n 's/^n_base = //p' "$toml")
    if [ -n "$N" ] && [ "$n" != "$N" ]; then
        echo "ERROR: the designs of the sets have different n_base ($N and $n)" >&2
        exit 1
    fi
    N=$n
    for d in $DOMAINS; do
        [ -d "$SA_BASE/domains/$d" ] || { echo "ERROR: domain $SA_BASE/domains/$d not found" >&2; exit 1; }
        if [ "$s" = main ]; then
            [ -f "$CMR_PF_INPUTS_DIR/$d/dhw_draws.jld2" ] || { echo "ERROR: $CMR_PF_INPUTS_DIR/$d/dhw_draws.jld2 not found" >&2; exit 1; }
            JOBS="$JOBS $BAND:$d"
        else
            JOBS="$JOBS $BAND:$d:$s"
        fi
    done
done
JOBS="${JOBS# }"
LAST="${LAST:-$N}"
[ "$LAST" -le "$N" ] || { echo "ERROR: --last $LAST is beyond the design's n_base $N" >&2; exit 1; }
[ "$FIRST" -le "$LAST" ] || { echo "ERROR: --first $FIRST is after --last $LAST" >&2; exit 1; }
# The last task may reach past LAST; 03_run_chunk.jl only caps at n_base.
if [ "$LAST" -lt "$N" ] && [ $(( (LAST - FIRST + 1) % WALKS )) -ne 0 ]; then
    echo "ERROR: with --last below n_base, (--last - --first + 1) must be a multiple of --walks" >&2
    exit 1
fi
TASKS=$(( (LAST - FIRST + WALKS) / WALKS ))

mkdir -p "$SA_BASE/logs"
cd "$SA_BASE/logs"
sub() {   # sub <sbatch args...> -> job id (or a placeholder on --dry-run)
    if [ "$DRY" = "1" ]; then echo "sbatch $*" >&2; echo "DRY$RANDOM"; else sbatch --parsable "$@"; fi
}
tag="${MODEL}_band${bdir##*band}"
j=$(sub --job-name="pf_${tag}" --array=1-"$TASKS" --time="$TIME" \
    --export=ALL,WALKS="$WALKS",FIRST="$FIRST",JOBS="$JOBS" "$SA_BASE/jobs/pf_chunk.slurm")
echo "walks $FIRST-$LAST of n_base $N, $WALKS per task -> $TASKS tasks, job $j"
echo "  jobs: $JOBS"
c=""
if [ "$COLLECT" = "1" ]; then
    c=$(sub --job-name="pf_${tag}_collect" --dependency=afterok:"$j" --export=ALL,BAND="$BAND" "$SA_BASE/jobs/pf_collect.slurm")
    echo "collect: job $c (after $j)"
fi
if [ "$DRY" = "0" ]; then
    mkdir -p "$bdir"
    echo "walks_${FIRST}-${LAST} $j" >> "$bdir/jobs.txt"
    [ -z "$c" ] || echo "collect $c" >> "$bdir/jobs.txt"
fi
echo
echo "Follow with: squeue -u ${USER:-$(whoami)}    (logs: $SA_BASE/logs; per-task time and memory: $bdir/logs)"
