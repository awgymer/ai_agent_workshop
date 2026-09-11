#!/usr/bin/env bash
# Golden tests: diff mytools against real bedtools v2.31.1 (SPEC §8).
#
# Usage: ./tests/run_golden.sh [subcommand...]
#   Runs tests/golden/<subcommand>.sh for each subcommand given, or for all of them.
#   MYTOOLS overrides the binary under test (default: target/release/mytools if built,
#   else mytools on PATH). MYTOOLS=bedtools compares bedtools with itself, which checks
#   the harness.
#
# Helpers available to tests/golden/*.sh:
#   check       <name> -- <args...>          stdout and exit code must match bedtools
#   check_stdin <name> <file> -- <args...>   as check, with <file> on stdin for both
#   check_bam   <name> -- <args...>          BAM output compared after `samtools view -h`
#   check_exit  <name> <code> -- <args...>   mytools only: fixed exit code, empty stdout
#                                            (the SPEC §7 deviations; skipped for
#                                            MYTOOLS=bedtools)
#   sorted      <file>                       prints a path to <file> run through
#                                            `bedtools sort` (merge and -sorted need it)
#   $DATA, $FIXTURES                         data/ and tests/fixtures/
set -uo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
DATA=$ROOT/data
FIXTURES=$ROOT/tests/fixtures
GOLDEN=$ROOT/tests/golden
SUBCOMMANDS=(sort merge intersect window subtract)

if [[ -z ${MYTOOLS:-} ]]; then
  if [[ -x $ROOT/target/release/mytools ]]; then MYTOOLS=$ROOT/target/release/mytools; else MYTOOLS=mytools; fi
fi

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
pass=0; fail=0; skip=0

ok()      { echo "ok   $1"; pass=$((pass + 1)); }
failed()  { echo "FAIL $1"; fail=$((fail + 1)); }
show()    { sed 's/^/      /' "$1" | head -"${2:-20}"; }

# run <tool> <stdin-file|""> <args...>: stdout to $tmp/<tool>.out, stderr to .err, sets rc
run() {
  local tool=$1 label=$2 input=$3; shift 3
  if [[ -n $input ]]; then
    "$tool" "$@" < "$input" > "$tmp/$label.out" 2> "$tmp/$label.err"
  else
    "$tool" "$@" < /dev/null > "$tmp/$label.out" 2> "$tmp/$label.err"
  fi
  rc=$?
}

compare_text() {
  local name=$1 got_rc=$2 want_rc=$3
  if [[ $got_rc -ne $want_rc ]]; then
    failed "$name (exit $got_rc, bedtools gave $want_rc)"; show "$tmp/got.err" 3; return
  fi
  if cmp -s "$tmp/want.out" "$tmp/got.out"; then
    ok "$name"
  else
    failed "$name"; diff -u "$tmp/want.out" "$tmp/got.out" > "$tmp/diff"; show "$tmp/diff"
  fi
}

check() {
  local name=$1; shift 2           # drop the literal --
  run "$MYTOOLS" got "" "$@"; local got_rc=$rc
  run bedtools want "" "$@"; local want_rc=$rc
  compare_text "$name" "$got_rc" "$want_rc"
}

check_stdin() {
  local name=$1 input=$2; shift 3  # drop the literal --
  run "$MYTOOLS" got "$input" "$@"; local got_rc=$rc
  run bedtools want "$input" "$@"; local want_rc=$rc
  compare_text "$name" "$got_rc" "$want_rc"
}

check_bam() {
  local name=$1; shift 2
  run "$MYTOOLS" got "" "$@"; local got_rc=$rc
  run bedtools want "" "$@"; local want_rc=$rc
  if [[ $got_rc -ne $want_rc ]]; then
    failed "$name (exit $got_rc, bedtools gave $want_rc)"; show "$tmp/got.err" 3; return
  fi
  if [[ $want_rc -ne 0 ]]; then ok "$name"; return; fi
  if ! samtools view -h "$tmp/want.out" > "$tmp/want.sam" 2> "$tmp/want.samerr"; then
    failed "$name (bedtools output is not BAM)"; show "$tmp/want.samerr" 3; return
  fi
  if ! samtools view -h "$tmp/got.out" > "$tmp/got.sam" 2> "$tmp/got.samerr"; then
    failed "$name (output is not BAM)"; show "$tmp/got.samerr" 3; return
  fi
  if cmp -s "$tmp/want.sam" "$tmp/got.sam"; then
    ok "$name"
  else
    failed "$name"; diff -u "$tmp/want.sam" "$tmp/got.sam" > "$tmp/diff"; show "$tmp/diff"
  fi
}

check_exit() {
  local name=$1 code=$2; shift 3
  if [[ $(basename "$MYTOOLS") == bedtools ]]; then
    echo "skip $name (deviation from bedtools)"; skip=$((skip + 1)); return
  fi
  run "$MYTOOLS" got "" "$@"
  if [[ $rc -ne $code ]]; then
    failed "$name (exit $rc, expected $code)"; show "$tmp/got.err" 3
  elif [[ -s $tmp/got.out ]]; then
    failed "$name (stdout not empty)"; show "$tmp/got.out" 5
  else
    ok "$name"
  fi
}

sorted() {
  local base; base=$(basename "$1"); base=${base%.gz}
  local out=$tmp/sorted.$base
  [[ -e $out ]] || bedtools sort -i "$1" > "$out"
  echo "$out"
}

if [[ $# -eq 0 ]]; then set -- "${SUBCOMMANDS[@]}"; fi
for cmd in "$@"; do
  if [[ ! -f $GOLDEN/$cmd.sh ]]; then
    echo "run_golden.sh: unknown subcommand $cmd (no tests/golden/$cmd.sh)" >&2
    exit 2
  fi
done

for cmd in "$@"; do
  echo "== $cmd"
  # shellcheck source=/dev/null
  source "$GOLDEN/$cmd.sh"
done

echo "---"
echo "$pass passed, $fail failed, $skip skipped"
[[ $fail -eq 0 ]]
