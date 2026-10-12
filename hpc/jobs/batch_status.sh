#!/bin/bash
# batch_status.sh -- how every project under runs/ is doing: each job's Slurm state (from jobs.txt and
# sacct), then, once the analysis has written FIGURES.md, the scenarios each engine completed.
#
#   cd ~/CoralModelReview
#   bash jobs/batch_status.sh                          # every project
#   bash jobs/batch_status.sh SingleReef NineReefsLine # only these domains' (or projects') folders
#
# States: ok = COMPLETED, run = RUNNING, wait = PENDING, FAIL / OOM (out of memory) / TIMEOUT /
# CANCEL = needs a re-run (6d), ? = not in sacct (too old, or not submitted). An engine job can be "ok"
# with fewer scenarios than requested (ReefMod skips failed scenarios), hence the second line.
set -euo pipefail
source "${CMR_BASE:-$HOME/CoralModelReview}/jobs/common.sh"
cd "$RUNS"
roots=("$@")
[ ${#roots[@]} -gt 0 ] || roots=(.)

short() {
    case "$1" in
        COMPLETED) echo ok ;; RUNNING) echo run ;; PENDING) echo wait ;; OUT_OF_MEMORY) echo OOM ;;
        CANCELLED) echo CANCEL ;; "") echo "?" ;; *) echo "$1" ;;
    esac
}

n=0; n_done=0; n_bad=0
while read -r f; do
    p=$(dirname "${f#./}")
    ids=$(awk '{print $2}' "$f" | paste -sd, -)
    [ -n "$ids" ] || continue
    n=$((n + 1))
    unset st; declare -A st=()
    # one state per job; array jobs (reefmod_B) report per task: keep any task that is not COMPLETED
    while IFS='|' read -r id state; do
        k=${id%%_*}; state=${state%% *}
        [ -z "${st[$k]:-}" ] || [ "${st[$k]}" = COMPLETED ] && st[$k]=$state
    done < <(sacct -X -n -P -j "$ids" --format=JobID,State 2>/dev/null || true)
    line=""; bad=0; busy=0
    while read -r name id; do
        s=$(short "${st[$id]:-}")
        line+="$name $s  "
        case "$s" in ok|"?") ;; run|wait) busy=1 ;; *) bad=1 ;; esac
    done < "$f"
    if [ -f "$p/FIGURES.md" ] && [ "$busy" = 0 ]; then
        status=DONE; n_done=$((n_done + 1))
    elif [ "$busy" = 1 ]; then status=busy
    else status="NO FIGURES"; fi
    [ "$bad" = 1 ] && { status="$status, CHECK"; n_bad=$((n_bad + 1)); }
    printf '%-26s %-20s %s\n' "$p" "[$status]" "$line"
    if [ -f "$p/FIGURES.md" ]; then
        sc=$(grep -m1 '^- Scenarios:' "$p/FIGURES.md" | sed 's/^- Scenarios: //') || true
        miss=$(grep -m1 'Not included' "$p/FIGURES.md" | sed 's/^- //; s/\*\*//g') || true
        nfail=$(awk '/^## Figures that failed/ { f = 1; next } f && /^- / { n++ } END { print n + 0 }' "$p/FIGURES.md")
        printf '%-47s %s%s%s\n' "" "scenarios: ${sc:-?}" "${miss:+; $miss}" \
            "$( [ "${nfail:-0}" -gt 0 ] && echo "; $nfail failed figure(s)" )"
    fi
done < <(find "${roots[@]}" -name jobs.txt | sort)

echo "$n projects: $n_done done, $n_bad with a failed/cancelled job (re-run with: bash jobs/submit_all.sh <project> --engines <engine>, same CMR_* settings)"
