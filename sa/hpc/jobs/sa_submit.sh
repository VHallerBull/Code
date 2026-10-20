#!/bin/bash
# sa_submit.sh -- submit one model's Sobol' runs: per domain, the main design as a Slurm array, the
# zero-DHW control as a second array, then collect + figures once both have finished.
#
#   bash ~/CoralModelReview/sa/jobs/sa_submit.sh --model coralpop --domains "SingleReef_SobolDHW NineReefs_SobolDHW"
#
# Options:
#   --model M          coralpop | coralblox | ... (needs images/M_sa.sif and M's designs, see 04_design.jl)
#   --domains "A B"    default: both synthetic domains
#   --chunk K          runs per array task (default 1024)
#   --n N              run only the first N rows of the design (test runs; default: all rows)
#   --n-control C      rows of the zero-DHW control (default 1024; capped at N); 0 = no control
#   --set S            main (default) or an extension block b2, b3, ... (no control is run for blocks)
#   --dry-run          print the sbatch commands only
#
# Writes <run dir>/jobs.txt (name, job id). Slurm logs go to ~/CoralModelReview/sa/logs.
set -euo pipefail
export SA_BASE="${SA_BASE:-$HOME/CoralModelReview/sa}"

MODEL="" DOMAINS="SingleReef_SobolDHW NineReefs_SobolDHW" CHUNK=1024 N="" N_CONTROL=1024 SET=main DRY=0
while [ $# -gt 0 ]; do
    case "$1" in
        --model) MODEL="$2"; shift 2 ;;
        --domains) DOMAINS="$2"; shift 2 ;;
        --chunk) CHUNK="$2"; shift 2 ;;
        --n) N="$2"; shift 2 ;;
        --n-control) N_CONTROL="$2"; shift 2 ;;
        --set) SET="$2"; shift 2 ;;
        --dry-run) DRY=1; shift ;;
        *) echo "unknown option $1" >&2; exit 1 ;;
    esac
done
[ -n "$MODEL" ] || { echo "usage: sa_submit.sh --model M [--domains \"...\"] [--chunk K] [--n N] [--set S] [--dry-run]" >&2; exit 1; }
export MODEL
source "$SA_BASE/jobs/sa_common.sh"
[ -f "$SIF" ] || { echo "ERROR: image not found: $SIF" >&2; exit 1; }
[ "$SET" = "main" ] || N_CONTROL=0

mkdir -p "$SA_BASE/logs"
cd "$SA_BASE/logs"
sub() {   # sub <sbatch args...> -> job id (or a placeholder on --dry-run)
    if [ "$DRY" = "1" ]; then echo "sbatch $*" >&2; echo "DRY$RANDOM"; else sbatch --parsable "$@"; fi
}

for d in $DOMAINS; do
    dir="$(run_dir "$d")"
    design="$dir/design.csv"; [ "$SET" = "main" ] || design="$dir/design_$SET.csv"
    [ -f "$design" ] || { echo "ERROR: $design not found -- upload the designs (04_design.jl)" >&2; exit 1; }
    [ -d "$SA_BASE/domains/$d" ] || { echo "ERROR: domain $SA_BASE/domains/$d not found" >&2; exit 1; }
    rows=$(( $(wc -l < "$design") - 1 ))
    n=${N:-$rows}; [ "$n" -le "$rows" ] || n=$rows
    nc=$N_CONTROL; [ "$nc" -le "$n" ] || nc=$n
    partial=false; [ "$n" -lt "$rows" ] && partial=true
    tag="${MODEL}_${d%%_*}"

    ids=()
    j=$(sub --job-name="sa_${tag}_$SET" --array=1-$(( (n + CHUNK - 1) / CHUNK )) \
        --export=ALL,DOMAIN="$d",SET="$SET",CHUNK_SIZE="$CHUNK" "$SA_BASE/jobs/sa_chunk.slurm")
    ids+=("$j"); echo "$d $SET: $n runs in $(( (n + CHUNK - 1) / CHUNK )) tasks -> job $j"
    if [ "$nc" -gt 0 ]; then
        c=$(sub --job-name="sa_${tag}_control" --array=1-$(( (nc + CHUNK - 1) / CHUNK )) \
            --export=ALL,DOMAIN="$d",SET=control,CHUNK_SIZE="$CHUNK",CMR_SA_N_CONTROL="$nc" "$SA_BASE/jobs/sa_chunk.slurm")
        ids+=("$c"); echo "$d control: $nc runs -> job $c"
    fi
    dep=$(IFS=:; echo "${ids[*]}")
    p=$(sub --job-name="sa_${tag}_post" --dependency=afterok:"$dep" \
        --export=ALL,DOMAIN="$d",CMR_SA_N_CONTROL="$nc",CMR_SA_ALLOW_PARTIAL="$partial" "$SA_BASE/jobs/sa_post.slurm")
    echo "$d collect + figures: job $p (after $dep)"
    if [ "$DRY" = "0" ]; then
        mkdir -p "$dir"
        { echo "${SET} $j"; [ "$nc" -gt 0 ] && echo "control $c"; echo "post $p"; } >> "$dir/jobs.txt"
    fi
done
echo
echo "Follow with: squeue -u $USER    (logs: $SA_BASE/logs)"
