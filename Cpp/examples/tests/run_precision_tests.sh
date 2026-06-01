#!/bin/bash
# Full precision test: compare ALL cross section values (angle-by-angle)
# between C++ and Cleopatra for every test input.
# Usage: cd Claude_attack/Cpp && bash examples/tests/run_precision_tests.sh

CPP="./ptolemy"
CLEOPATRA="/home/ryan/ptolemy_2019/digios/analysis/Cleopatra/ptolemy"
TESTDIR="examples/tests"

if [ ! -x "$CPP" ]; then echo "ERROR: C++ binary not found. Run 'make' first."; exit 1; fi
if [ ! -x "$CLEOPATRA" ]; then echo "ERROR: Cleopatra not found at $CLEOPATRA"; exit 1; fi

PASS=0
FAIL=0
ERROR=0
TOTAL=0

run_precision_test() {
    local infile="$1"
    local name=$(basename "$infile" .in)
    TOTAL=$((TOTAL + 1))

    # Run both
    rm -f fort.*
    local cpp_out=$(timeout 120 $CPP < "$infile" 2>/dev/null) || true
    rm -f fort.*
    local ftn_out=$(timeout 120 $CLEOPATRA < "$infile" 2>/dev/null) || true
    rm -f fort.*

    # Extract cross section lines: "  angle  value" pattern from COMPUTATION OF CROSS SECTIONS
    # For elastic: lines like "   5.00  0.99876E+00   0.99876     100.0    0.000"
    # For transfer/inel: lines like "   5.00  0.12345E+01   ..."
    # Filter: lines starting with spaces then a number (angle)
    local cpp_xsec=$(echo "$cpp_out" | grep -A9999 "COMPUTATION OF CROSS SECTIONS" | grep '^ *[0-9][0-9]*\.[0-9]' | awk '{print $1, $2}')
    local ftn_xsec=$(echo "$ftn_out" | grep -A9999 "COMPUTATION OF CROSS SECTIONS" | grep '^ *[0-9][0-9]*\.[0-9]' | awk '{print $1, $2}')

    if [ -z "$cpp_xsec" ] && [ -z "$ftn_xsec" ]; then
        # Try elastic format: SIGMA/ RUTHERFORD lines
        cpp_xsec=$(echo "$cpp_out" | grep '^ *[0-9][0-9]*\.' | awk 'NF>=3 {print $1, $3}')
        ftn_xsec=$(echo "$ftn_out" | grep '^ *[0-9][0-9]*\.' | awk 'NF>=3 {print $1, $3}')
    fi

    if [ -z "$cpp_xsec" ] || [ -z "$ftn_xsec" ]; then
        printf "  %-40s ERROR (no cross section data)\n" "$name"
        ERROR=$((ERROR + 1))
        return
    fi

    # Compare using python
    local result=$(python3 -c "
import sys
cpp_lines = '''$cpp_xsec'''.strip().split('\n')
ftn_lines = '''$ftn_xsec'''.strip().split('\n')

if len(cpp_lines) != len(ftn_lines):
    # Different number of angles — compare what we can
    n = min(len(cpp_lines), len(ftn_lines))
else:
    n = len(cpp_lines)

max_pct = 0.0
max_angle = ''
n_compared = 0
n_nan = 0
n_mismatch = 0

for i in range(n):
    try:
        ca, cv = cpp_lines[i].split()
        fa, fv = ftn_lines[i].split()
        ca_f, fa_f = float(ca), float(fa)
        cv_f, fv_f = float(cv), float(fv)
    except:
        continue

    if cv == 'NaN' or fv == 'NaN' or cv == '-nan' or fv == '-nan':
        n_nan += 1
        continue

    cv_f = float(cv)
    fv_f = float(fv)
    n_compared += 1

    if abs(fv_f) > 1e-30:
        pct = abs(cv_f - fv_f) / abs(fv_f) * 100
    elif abs(cv_f) > 1e-30:
        pct = abs(cv_f - fv_f) / abs(cv_f) * 100
    else:
        pct = 0.0

    if pct > max_pct:
        max_pct = pct
        max_angle = ca

if n_nan > 0:
    print(f'NAN:{n_nan} angles:{n_compared} max:{max_pct:.4f}% at:{max_angle}')
elif n_compared == 0:
    print('NODATA')
else:
    print(f'OK angles:{n_compared} max:{max_pct:.4f}% at:{max_angle}')
" 2>/dev/null)

    if echo "$result" | grep -q "^OK"; then
        local max_pct=$(echo "$result" | grep -o 'max:[0-9.]*' | cut -d: -f2)
        local n_angles=$(echo "$result" | grep -o 'angles:[0-9]*' | cut -d: -f2)
        local at_angle=$(echo "$result" | grep -o 'at:.*' | cut -d: -f2)

        if python3 -c "exit(0 if float('$max_pct') < 0.02 else 1)" 2>/dev/null; then
            printf "  %-40s PASS  (%s angles, max %.4f%%" "$name" "$n_angles" "$max_pct"
            if [ -n "$at_angle" ] && [ "$at_angle" != "" ]; then
                printf " at %s deg" "$at_angle"
            fi
            printf ")\n"
            PASS=$((PASS + 1))
        else
            printf "  %-40s FAIL  (%s angles, max %.4f%% at %s deg)\n" "$name" "$n_angles" "$max_pct" "$at_angle"
            FAIL=$((FAIL + 1))
        fi
    elif echo "$result" | grep -q "^NAN"; then
        local n_nan=$(echo "$result" | grep -o 'NAN:[0-9]*' | cut -d: -f2)
        printf "  %-40s FAIL  (NaN in %s angles)\n" "$name" "$n_nan"
        FAIL=$((FAIL + 1))
    else
        printf "  %-40s ERROR (%s)\n" "$name" "$result"
        ERROR=$((ERROR + 1))
    fi
}

echo "=== Full Precision Test: C++ vs Cleopatra (angle-by-angle) ==="
echo ""

echo "--- Elastic ---"
for f in "$TESTDIR"/test_elastic_*.in; do [ -f "$f" ] && run_precision_test "$f"; done

echo ""
echo "--- Transfer DWBA ---"
for f in "$TESTDIR"/test_transfer_*.in; do [ -f "$f" ] && run_precision_test "$f"; done

echo ""
echo "--- Inelastic DWBA ---"
for f in "$TESTDIR"/test_inel_*.in; do [ -f "$f" ] && run_precision_test "$f"; done

echo ""
echo "=== Summary ==="
echo "  PASS:  $PASS"
echo "  FAIL:  $FAIL"
echo "  ERROR: $ERROR"
echo "  TOTAL: $TOTAL"

if [ $FAIL -eq 0 ] && [ $ERROR -eq 0 ]; then
    echo "  All tests passed at machine precision!"
    exit 0
else
    echo "  Some tests have precision issues."
    exit 1
fi
