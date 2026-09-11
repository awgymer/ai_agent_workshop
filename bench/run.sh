#!/usr/bin/env bash
# Benchmark mytools against bedtools on the same arguments (SPEC §6, §8 item 6).
# Median of RUNS runs (default 5) of wall-clock time and maximum resident set size, as
# `/usr/bin/time -v` reports it. Not a CI check: it reports, it doesn't fail.
#
# Usage: bench/run.sh <subcommand> <args...>
#   e.g. bench/run.sh merge -i sorted.gtf
#   MYTOOLS overrides the binary (default: target/release/mytools if built, else PATH).
# Targets: peak memory <= bedtools (a miss is a bug); time <= bedtools + 5% (investigate).
set -uo pipefail

if [[ $# -lt 1 ]]; then
  echo "usage: bench/run.sh <subcommand> <args...>" >&2
  exit 2
fi
if [[ ! -x /usr/bin/time ]]; then
  echo "bench/run.sh: needs GNU time at /usr/bin/time" >&2
  exit 2
fi

ROOT=$(cd "$(dirname "$0")/.." && pwd)
RUNS=${RUNS:-5}
if [[ -z ${MYTOOLS:-} ]]; then
  if [[ -x $ROOT/target/release/mytools ]]; then MYTOOLS=$ROOT/target/release/mytools; else MYTOOLS=mytools; fi
fi
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

median() { sort -n | awk '{v[NR] = $1} END {print v[int((NR + 1) / 2)]}'; }

# measure <tool> <args...>: prints "<median seconds> <median max RSS in KB> <last exit code>"
measure() {
  local times=() mems=() rc=0 start end
  for ((i = 0; i < RUNS; i++)); do
    start=$(date +%s%N)
    /usr/bin/time -v -o "$tmp/time" "$@" > /dev/null 2> "$tmp/err"
    rc=$?
    end=$(date +%s%N)
    times+=("$(awk -v ns=$((end - start)) 'BEGIN {printf "%.4f", ns / 1e9}')")
    mems+=("$(awk -F': ' '/Maximum resident set size/ {print $2}' "$tmp/time")")
  done
  echo "$(printf '%s\n' "${times[@]}" | median) $(printf '%s\n' "${mems[@]}" | median) $rc"
}

read -r bt_time bt_mem bt_rc < <(measure bedtools "$@")
read -r my_time my_mem my_rc < <(measure "$MYTOOLS" "$@")

echo "bench: $* (median of $RUNS)"
printf '%-10s %14s %14s %6s\n' tool "wall time (s)" "max RSS (KB)" exit
printf '%-10s %14s %14s %6s\n' bedtools "$bt_time" "$bt_mem" "$bt_rc"
printf '%-10s %14s %14s %6s\n' mytools "$my_time" "$my_mem" "$my_rc"

awk -v bt="$bt_time" -v my="$my_time" -v btm="$bt_mem" -v mym="$my_mem" 'BEGIN {
  if (mym <= btm) print "memory: ok (<= bedtools)";
  else printf "memory: ABOVE bedtools by %d KB (a bug, SPEC §6)\n", mym - btm;
  if (bt == 0) print "time: bedtools took 0 s; no ratio";
  else if (my <= bt * 1.05) printf "time: ok (%+.1f%% vs bedtools)\n", (my / bt - 1) * 100;
  else printf "time: %+.1f%% vs bedtools, over the +5%% target (investigate)\n", (my / bt - 1) * 100;
}'
