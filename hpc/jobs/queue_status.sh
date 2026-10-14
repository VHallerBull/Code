#!/bin/bash
# queue_status.sh -- how busy the cluster is, where your jobs stand, and how long they should take.
#
#   cd ~/CoralModelReview
#   bash jobs/queue_status.sh              # partition defq, run-time history of the last 60 days
#   bash jobs/queue_status.sh defq 120     # another partition / a longer history
#
#   1. Cluster load (sinfo): nodes by state, CPUs allocated / idle.
#   2. The queue (squeue): everybody's running and pending jobs, and why the pending ones wait.
#   3. Typical run time per job name (cmr_coralblox, cmr_cscape, ...) from YOUR finished jobs (sacct).
#   4. Your jobs now, and an estimate of when each finishes: Slurm's expected start (squeue --start)
#      for a pending job, or the time already run for a running one, plus the typical run time.
#
# It is a guide, not a promise. Slurm's start times move as other people's jobs finish or are submitted,
# and run time depends on the domain (a nine-reef domain runs longer than a single reef), which the
# per-name history does not know about -- compare the min / max columns with the median.
set -uo pipefail

PART="${1:-defq}"
DAYS="${2:-60}"
ME="${USER:-$(id -un)}"
command -v squeue >/dev/null 2>&1 || { echo "squeue not found: run this on the HPC login node." >&2; exit 1; }

fmt() {   # seconds -> "2d 03h 05m" / "3h 05m" / "4m"
    local s=$1 d h m
    (( s < 0 )) && s=0
    s=$((s + 30))   # round to the nearest minute
    d=$((s / 86400)); h=$((s % 86400 / 3600)); m=$((s % 3600 / 60))
    if (( d > 0 )); then printf '%dd %02dh %02dm' "$d" "$h" "$m"
    elif (( h > 0 )); then printf '%dh %02dm' "$h" "$m"
    else printf '%dm' "$m"; fi
}
hms() {   # Slurm [D-][HH:]MM:SS -> seconds
    local t=$1 d=0 s=0 p
    [[ $t == *-* ]] && { d=${t%%-*}; t=${t#*-}; }
    IFS=: read -ra p <<< "$t"
    case ${#p[@]} in
        3) s=$((10#${p[0]} * 3600 + 10#${p[1]} * 60 + 10#${p[2]})) ;;
        2) s=$((10#${p[0]} * 60 + 10#${p[1]})) ;;
        1) s=$((10#${p[0]})) ;;
    esac
    echo $((d * 86400 + s))
}
median() { printf '%s\n' "$@" | sort -n | awk '{a[NR] = $1} END {print a[int((NR + 1) / 2)]}'; }
when() { date -d "@$1" '+%a %d %b %H:%M'; }

echo "== 1. Cluster, partition $PART"
sinfo -h -p "$PART" -o "%D %T" 2>/dev/null | awk '{n[$2] += $1} END {for (s in n) printf "   %4d nodes %s\n", n[s], s}' | sort -k3
sinfo -h -p "$PART" -o "%C" 2>/dev/null | awk -F/ '{a += $1; i += $2; o += $3; t += $4}
    END {printf "   CPUs: %d allocated, %d idle, %d unavailable, %d total -> %d%% of the usable CPUs busy\n", a, i, o, t, (a + i > 0 ? 100 * a / (a + i) : 0)}'

echo
echo "== 2. Queue, partition $PART (all users)"
squeue -h -p "$PART" -o "%T|%u|%C|%r" 2>/dev/null | awk -F'|' '
    $1 == "RUNNING" {r++; rc += $3; if (!($2 in ru)) {ru[$2] = 1; nru++}}
    $1 == "PENDING" {p++; pc += $3; why[$4]++; if (!($2 in pu)) {pu[$2] = 1; npu++}}
    END {
        printf "   running: %d jobs, %d CPUs, %d users\n", r, rc, nru
        printf "   pending: %d jobs, %d CPUs asked for, %d users\n", p, pc, npu
        for (w in why) printf "      pending for %-24s %d jobs\n", w, why[w]
    }'

echo
echo "== 3. Typical run time of your finished jobs (last $DAYS days)"
declare -A HIST MED
while IFS='|' read -r name sec; do
    [[ $name == cmr_* && $sec =~ ^[0-9]+$ ]] && HIST[$name]+="$sec "
done < <(sacct -u "$ME" -S "now-${DAYS}days" -X -n -P -s CD --format=JobName,ElapsedRaw 2>/dev/null)
if (( ${#HIST[@]} == 0 )); then
    echo "   no finished cmr_* jobs in the accounting history: run times below cannot be estimated."
else
    printf '   %-18s %5s %9s %9s %9s\n' job runs min median max
    for name in $(printf '%s\n' "${!HIST[@]}" | sort); do
        read -ra v <<< "${HIST[$name]}"
        MED[$name]=$(median "${v[@]}")
        sorted=($(printf '%s\n' "${v[@]}" | sort -n))
        printf '   %-18s %5d %9s %9s %9s\n' "$name" "${#v[@]}" "$(fmt "${sorted[0]}")" "$(fmt "${MED[$name]}")" "$(fmt "${sorted[-1]}")"
    done
fi

echo
echo "== 4. Your jobs ($ME)"
mapfile -t MINE < <(squeue -h -u "$ME" -o "%i|%j|%T|%M" 2>/dev/null | grep -E "^[^|]*\|cmr_")   # this project's jobs only
if (( ${#MINE[@]} == 0 )); then
    echo "   no cmr_* jobs in the queue."
    exit 0
fi
declare -A START
while IFS='|' read -r id st; do START[$id]=$st; done < <(squeue -h -u "$ME" --start -o "%i|%S" 2>/dev/null)

NOW=$(date +%s)
LAST=$NOW          # latest estimated finish so far
UNKNOWN=0          # jobs whose finish cannot be estimated
declare -A N_STATE
printf '   %-14s %-18s %-9s %-34s %s\n' jobid name state detail "expected finish"

show() { printf '   %-14s %-18s %-9s %-34s %s\n' "$1" "$2" "$3" "$4" "$5"; }
ANALYSE=()
for l in "${MINE[@]}"; do
    IFS='|' read -r id name state el <<< "$l"
    N_STATE[$state]=$(( ${N_STATE[$state]:-0} + 1 ))
    med=${MED[$name]:-}
    if [[ $name == cmr_analyse && $state == PENDING ]]; then ANALYSE+=("$id"); continue; fi
    case $state in
        RUNNING)
            e=$(hms "$el")
            if [[ -n $med ]]; then
                rem=$((med - e)); (( rem < 0 )) && rem=0
                fin=$((NOW + rem)); (( fin > LAST )) && LAST=$fin
                show "$id" "$name" RUNNING "run $(fmt "$e"), ~$(fmt "$rem") left" "$(when "$fin")"
            else
                UNKNOWN=$((UNKNOWN + 1)); show "$id" "$name" RUNNING "run $(fmt "$e"), no history" "?"
            fi ;;
        PENDING)
            st=${START[$id]:-N/A}
            if [[ $st =~ ^[0-9]{4}- ]] && s=$(date -d "$st" +%s 2>/dev/null); then
                if [[ -n $med ]]; then
                    fin=$((s + med)); (( fin > LAST )) && LAST=$fin
                    show "$id" "$name" PENDING "starts $(when "$s"), runs ~$(fmt "$med")" "$(when "$fin")"
                else
                    UNKNOWN=$((UNKNOWN + 1)); show "$id" "$name" PENDING "starts $(when "$s"), no history" "?"
                fi
            else
                UNKNOWN=$((UNKNOWN + 1)); show "$id" "$name" PENDING "no start estimate yet${med:+, runs ~$(fmt "$med")}" "?"
            fi ;;
        *) show "$id" "$name" "$state" "" "" ;;
    esac
done
for id in "${ANALYSE[@]}"; do   # the analysis waits for the engine jobs it depends on
    med=${MED[cmr_analyse]:-}
    if [[ -n $med ]] && (( UNKNOWN == 0 )); then
        fin=$((LAST + med)); LAST=$fin
        show "$id" cmr_analyse PENDING "after the others, runs ~$(fmt "$med")" "$(when "$fin")"
    else
        UNKNOWN=$((UNKNOWN + 1)); show "$id" cmr_analyse PENDING "after the others" "?"
    fi
done

echo
summary=""; for s in "${!N_STATE[@]}"; do summary+="${N_STATE[$s]} $s, "; done
echo "   ${summary%, }."
if (( UNKNOWN == 0 )); then
    echo "   Everything should be finished around $(when "$LAST") (in about $(fmt $((LAST - NOW))))."
else
    echo "   Estimated so far: finished by $(when "$LAST") (in about $(fmt $((LAST - NOW)))), but $UNKNOWN job(s) have no estimate"
    echo "   (no start time from Slurm yet, or no run-time history), so the real finish is later."
fi
