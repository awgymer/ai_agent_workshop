# Golden cases for `sort`, sourced by tests/run_golden.sh (helpers are documented there).
# Owned by the sort branch (#5). BED is 0-based, half-open.

check "sort a.bed" -- sort -i "$DATA/a.bed"
