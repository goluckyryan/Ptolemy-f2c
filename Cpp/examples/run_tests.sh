#!/bin/bash
# Comprehensive test suite: compare C++ ptolemy against Cleopatra reference
# Usage: cd Claude_attack/Cpp && bash examples/run_tests.sh

set -e

CPP="./ptolemy"
CLEOPATRA="/home/ryan/ptolemy_2019/digios/analysis/Cleopatra/ptolemy"
TESTDIR="examples"

if [ ! -x "$CPP" ]; then
    echo "ERROR: C++ binary not found. Run 'make' first."
    exit 1
fi
if [ ! -x "$CLEOPATRA" ]; then
    echo "ERROR: Cleopatra binary not found at $CLEOPATRA"
    exit 1
fi

PASS=0
FAIL=0
ERROR=0
TOTAL_TESTS=0

extract_value() {
    local output="$1"
    local val=""
    # Try 0TOTAL: line (transfer/inelastic)
    val=$(echo "$output" | grep "^0TOTAL:" | tail -1 | awk '{print $2}')
    if [ -n "$val" ]; then echo "$val"; return; fi
    # Try TOTAL REACTION CROSS SECTION (elastic)
    val=$(echo "$output" | grep "TOTAL REACTION CROSS SECTION" | grep -o '[0-9][0-9]*\.[0-9]*')
    if [ -n "$val" ]; then echo "$val"; return; fi
    echo ""
}

run_test() {
    local infile="$1"
    local name=$(basename "$infile" .in)
    TOTAL_TESTS=$((TOTAL_TESTS + 1))

    rm -f fort.*
    local cpp_out
    cpp_out=$(timeout 120 $CPP < "$infile" 2>/dev/null) || true
    rm -f fort.*

    local ftn_out
    ftn_out=$(timeout 120 $CLEOPATRA < "$infile" 2>/dev/null) || true
    rm -f fort.*

    local cpp_val=$(extract_value "$cpp_out")
    local ftn_val=$(extract_value "$ftn_out")

    if [ -z "$cpp_val" ] || [ -z "$ftn_val" ]; then
        if echo "$cpp_val$ftn_val" | grep -qi "nan"; then
            printf "  %-40s FAIL  (NaN: C++=%s Ftn=%s)\n" "$name" "$cpp_val" "$ftn_val"
            FAIL=$((FAIL + 1))
        else
            printf "  %-40s ERROR (C++=%s Ftn=%s)\n" "$name" "${cpp_val:-empty}" "${ftn_val:-empty}"
            ERROR=$((ERROR + 1))
        fi
        return
    fi

    if echo "$cpp_val" | grep -qi "nan" || echo "$ftn_val" | grep -qi "nan"; then
        printf "  %-40s FAIL  (NaN: C++=%s Ftn=%s)\n" "$name" "$cpp_val" "$ftn_val"
        FAIL=$((FAIL + 1))
        return
    fi

    if [ "$cpp_val" = "$ftn_val" ]; then
        printf "  %-40s PASS  (%s)\n" "$name" "$cpp_val"
        PASS=$((PASS + 1))
    else
        local pct=$(python3 -c "
a, b = float('$cpp_val'), float('$ftn_val')
if abs(b) > 1e-30:
    print(f'{abs(a-b)/abs(b)*100:.4f}')
elif abs(a) > 1e-30:
    print(f'{abs(a-b)/abs(a)*100:.4f}')
else:
    print('0.0000')
" 2>/dev/null || echo "???")
        if python3 -c "exit(0 if float('$pct') < 0.01 else 1)" 2>/dev/null; then
            printf "  %-40s PASS  (%s vs %s, %s%%)\n" "$name" "$cpp_val" "$ftn_val" "$pct"
            PASS=$((PASS + 1))
        else
            printf "  %-40s FAIL  (%s vs %s, %s%%)\n" "$name" "$cpp_val" "$ftn_val" "$pct"
            FAIL=$((FAIL + 1))
        fi
    fi
}

echo "=== Ptolemy C++ vs Cleopatra Comparison Test Suite ==="
echo ""

echo "--- Elastic ---"
for f in "$TESTDIR"/test_elastic_*.in; do
    [ -f "$f" ] && run_test "$f"
done

echo ""
echo "--- Transfer DWBA ---"
for f in "$TESTDIR"/test_transfer_*.in; do
    [ -f "$f" ] && run_test "$f"
done

echo ""
echo "--- Inelastic DWBA ---"
for f in "$TESTDIR"/test_inel_*.in; do
    [ -f "$f" ] && run_test "$f"
done

echo ""
echo "--- Coupled Channels ---"
for f in "$TESTDIR"/cc_*.in; do
    [ -f "$f" ] && run_test "$f"
done

echo ""
echo "=== Summary ==="
echo "  PASS:  $PASS"
echo "  FAIL:  $FAIL"
echo "  ERROR: $ERROR"
echo "  TOTAL: $TOTAL_TESTS"

if [ $FAIL -eq 0 ] && [ $ERROR -eq 0 ]; then
    echo "  All tests passed!"
    exit 0
else
    echo "  Some tests failed."
    exit 1
fi
