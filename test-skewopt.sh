#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
TEST_DIR=$(mktemp -d /tmp/skewopt-test.XXXXXX)
trap 'rm -rf -- "$TEST_DIR"' EXIT
CXX=${CXX:-g++}

"$CXX" -O2 -c "$SCRIPT_DIR/msieve_poly.cpp" -o "$TEST_DIR/msieve_poly.o"
"$CXX" -O2 "$SCRIPT_DIR/skewopt.cpp" "$TEST_DIR/msieve_poly.o" \
    -lgmp -o "$TEST_DIR/skewopt"

degree9_output=$(
    "$TEST_DIR/skewopt" -deg9 \
        17449402268886407318558803753801 \
        -42391158275216203514294433201 \
        -11 0 0 0 0 0 0 0 0 9
)
degree9_skew=$(awk '$1 == "Best" && $2 == "Skew:" { print $3 }' <<<"$degree9_output")
degree9_e=$(awk '$1 == "MurphyE:" { print $2 }' <<<"$degree9_output")
awk -v skew="$degree9_skew" -v murphy="$degree9_e" 'BEGIN {
    exit !(skew > 3.6e-6 && skew < 3.8e-6 &&
           murphy > 4.11e-16 && murphy < 4.12e-16)
}'

large_c0=$(printf '1%0117d' 0)
optimized_output=$(
    "$TEST_DIR/skewopt" -deg9 1 1 "$large_c0" 0 0 0 0 0 0 0 0 1
)
nominal_output=$(
    "$TEST_DIR/skewopt" -deg9 1 1 "$large_c0" 0 0 0 0 0 0 0 0 1 1e13
)
optimized_e=$(awk '$1 == "MurphyE:" { print $2 }' <<<"$optimized_output")
nominal_e=$(awk '$1 == "MurphyE:" { print $2 }' <<<"$nominal_output")
awk -v optimized="$optimized_e" -v nominal="$nominal_e" 'BEGIN {
    exit !(optimized >= nominal)
}'

if "$TEST_DIR/skewopt" 123456789 1 0 0 0 0 0 0 0 0 0 \
    >"$TEST_DIR/zero.out" 2>"$TEST_DIR/zero.err"; then
    echo "all-zero algebraic polynomial was unexpectedly accepted" >&2
    exit 1
fi
grep -Fq "must have degree at least 1" "$TEST_DIR/zero.err"

echo "skewopt: adaptive degree-9, nominal floor, and invalid input checks passed"
