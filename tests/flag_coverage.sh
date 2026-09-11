#!/usr/bin/env bash
# Flag coverage (SPEC §8): every flag listed by `bedtools <subcommand> -h` must be used by
# at least one case in tests/golden/<subcommand>.sh. Comment lines don't count.
#
# Usage: ./tests/flag_coverage.sh <subcommand>
# Exits 0 if every flag is covered, 1 listing the uncovered flags, 2 on bad usage.
set -uo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: tests/flag_coverage.sh <subcommand>" >&2
  exit 2
fi
cmd=$1
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cases=$ROOT/tests/golden/$cmd.sh
if [[ ! -f $cases ]]; then
  echo "flag_coverage.sh: no tests/golden/$cmd.sh" >&2
  exit 2
fi

# Option lines in `-h` start with a tab and a dash: "\t-sizeA\t\tSort by ...".
mapfile -t flags < <(bedtools "$cmd" -h 2>&1 | grep -oE $'^\t-[A-Za-z]+' | tr -d '\t' | sort -u)
if [[ ${#flags[@]} -eq 0 ]]; then
  echo "flag_coverage.sh: no flags found in 'bedtools $cmd -h'" >&2
  exit 2
fi

code=$(grep -vE '^[[:space:]]*#' "$cases")
missing=()
for flag in "${flags[@]}"; do
  if ! grep -qE -- "(^|[[:space:]\"'])${flag}([[:space:]\"']|$)" <<< "$code"; then
    missing+=("$flag")
  fi
done

if [[ ${#missing[@]} -eq 0 ]]; then
  echo "$cmd: all ${#flags[@]} flags covered"
  exit 0
fi
echo "$cmd: ${#missing[@]} of ${#flags[@]} flags have no golden case in tests/golden/$cmd.sh:"
printf '  %s\n' "${missing[@]}"
exit 1
