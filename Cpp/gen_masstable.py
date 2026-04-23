#!/usr/bin/env python3
"""Generate masstable.cpp from masstable.f - extracts DATA arrays precisely."""
import re
import sys

FORTRAN_FILE = '../../src/masstable.f'
OUTPUT_FILE = 'src/masstable.cpp'

def read_fortran():
    with open(FORTRAN_FILE) as f:
        return f.readlines()

def extract_data_values(lines, start_pattern, end_patterns):
    """Extract numeric values from Fortran DATA statements."""
    values = []
    in_block = False
    for line in lines:
        stripped = line.rstrip()
        # Check for start
        if start_pattern in stripped:
            in_block = True
        # Check for end
        if in_block:
            for ep in end_patterns:
                if ep in stripped and start_pattern not in stripped:
                    in_block = False
                    break
        if not in_block:
            continue
        # Extract numbers from this line
        # Remove Fortran comment/continuation markers
        text = stripped
        # Remove column 1-6 (label/continuation area) for continuation lines
        if len(text) > 5 and text[5] in ('&', '*', '+'):
            text = text[6:]
        elif text.startswith('     &'):
            text = text[6:]
        # Remove DATA statement header
        text = re.sub(r'DATA\s*\(.*?\)\s*/', '', text)
        text = re.sub(r'DATA\s+LZLO\s*/', '', text)
        text = re.sub(r'DATA\s+LZHI\s*/', '', text)
        # Remove trailing /
        text = text.rstrip()
        if text.endswith('/'):
            text = text[:-1]
        # Remove D0 suffix (Fortran double precision)
        text = re.sub(r'D0\b', '', text)
        text = re.sub(r'D\+?0+\b', '', text)
        # Remove column numbers in column 73+
        # Fortran fixed-form: columns 73+ are sequence numbers
        # We need to be careful - only remove if there's a clear column boundary
        # The format seems to have numbers like "10   " or "100  " at end
        # Let me handle by looking for the pattern
        text = re.sub(r'\s+\d{1,5}\s*$', '', text)
        # Extract all numeric values (integers and floats, possibly negative)
        nums = re.findall(r'[+-]?\d+\.?\d*', text)
        for n in nums:
            values.append(n)
    return values

def extract_lzlo(lines):
    """Extract LZLO array values."""
    values = []
    in_block = False
    for line in lines:
        stripped = line.rstrip()
        if 'DATA LZLO' in stripped:
            in_block = True
        elif 'DATA LZHI' in stripped:
            in_block = False
            continue
        if not in_block:
            continue
        # Get the data part
        text = stripped
        # Remove DATA LZLO /
        text = re.sub(r'DATA\s+LZLO\s*/', '', text)
        # Remove continuation char
        text = re.sub(r'^\s{5}&\s*', '', text)
        text = re.sub(r'^\s+&\s*', '', text)
        # Remove trailing / and sequence numbers
        text = re.sub(r'/\s*$', '', text)
        text = re.sub(r'\s+\d{1,5}\s*$', '', text)
        # Extract integers
        nums = re.findall(r'\d+', text)
        values.extend([int(n) for n in nums])
    return values

def extract_lzhi(lines):
    """Extract LZHI array values."""
    values = []
    in_block = False
    for line in lines:
        stripped = line.rstrip()
        if 'DATA LZHI' in stripped:
            in_block = True
        elif re.match(r'\s+DATA\s*\(\s*XBUF', stripped):
            in_block = False
            continue
        if not in_block:
            continue
        text = stripped
        text = re.sub(r'DATA\s+LZHI\s*/', '', text)
        text = re.sub(r'^\s{5}&\s*', '', text)
        text = re.sub(r'^\s+&\s*', '', text)
        text = re.sub(r'/\s*$', '', text)
        text = re.sub(r'\s+\d{1,5}\s*$', '', text)
        nums = re.findall(r'\d+', text)
        values.extend([int(n) for n in nums])
    return values

def extract_xbuf(lines):
    """Extract XBUF data values."""
    values = []
    in_block = False
    for i, line in enumerate(lines):
        stripped = line.rstrip()
        # Start of XBUF DATA block
        if re.match(r'\s+DATA\s*\(\s*XBUF', stripped):
            in_block = True
            continue  # skip the header line itself
        # End detection: blank line or comment after data, or the lookup code
        if in_block:
            # Check if we've hit the end of all DATA statements
            # The code section starts around line 1770 with comments
            if stripped.startswith('C') or stripped.startswith('c'):
                # Comment lines after data block
                if '******' in stripped:
                    in_block = False
                    continue
                continue  # skip comment lines within
            if stripped == '' or stripped.isspace():
                continue  # skip blank lines
            # continuation of data or new DATA block
            if re.match(r'\s+DATA\s*\(\s*XBUF', stripped):
                continue  # next XBUF DATA segment header
        if not in_block:
            continue
        text = stripped
        # Remove continuation marker
        if len(text) > 5 and text[5] == '&':
            text = text[6:]
        elif text.startswith('     &'):
            text = text[6:]
        elif text.startswith('     *'):
            text = text[6:]
        elif text.startswith('     +'):
            text = text[6:]
        # Remove trailing / (end of DATA block)
        text = text.rstrip()
        if text.endswith('/'):
            text = text[:-1]
        # Remove sequence numbers (columns 73+)
        # These appear as whitespace + number at end of line
        text = re.sub(r'\s{2,}\d{1,5}\s*$', '', text)
        # Remove D0 suffix
        text = text.replace('D0', '')
        # Extract numeric values (including negative, with decimal points)
        nums = re.findall(r'[+-]?\d+\.?\d*', text)
        for n in nums:
            values.append(n)
    return values

def format_array_int(name, values, per_line=10):
    """Format an integer array for C++."""
    lines = []
    lines.append(f'static const int {name}[] = {{')
    lines.append(f'    0, // [0] unused (1-based indexing)')
    for i in range(0, len(values), per_line):
        chunk = values[i:i+per_line]
        line = '    ' + ', '.join(f'{v:4d}' for v in chunk)
        if i + per_line < len(values):
            line += ','
        lines.append(line)
    lines.append('};')
    return '\n'.join(lines)

def format_xbuf(values, per_line=6):
    """Format XBUF double array for C++."""
    lines = []
    lines.append('static double XBUF[9538] = {')
    lines.append('    0.0, // [0] unused (1-based indexing)')
    # Format in groups of 6 (2 triplets per line for readability)
    for i in range(0, len(values), per_line):
        chunk = values[i:i+per_line]
        formatted = []
        for v in chunk:
            # Keep the original precision
            if '.' in v:
                formatted.append(f'{v}')
            else:
                formatted.append(f'{v}.0')
        line = '    ' + ', '.join(formatted)
        if i + per_line < len(values):
            line += ','
        lines.append(line)
    lines.append('};')
    return '\n'.join(lines)

def generate_cpp(lzlo, lzhi, xbuf):
    """Generate the complete masstable.cpp file."""

    cpp = '''// masstable.cpp - Mass table routines from masstable.f
// Verbatim transliteration of Fortran to C++
//
// Contains:
//   EXCESS  - mass excess function
//   MASEXX  - mass excess lookup (2003 Mass Adjustment)
//   GETAS   - get mass numbers for element Z (ENTRY in MASEXX)
//   AZCODE  - decode element string (e.g., "CA40") to Z, A
//   CODEZ   - encode Z to element symbol (ENTRY in AZCODE)

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_forward.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_io.h"
#include <cstring>
#include <cctype>

// ====================================================================
// EXCESS function (lines 18-43 of masstable.f)
// Returns mass excess in MeV for given Z, A
// ====================================================================

double EXCESS(int KZ, int KA, int& NOTTAB) {
    double EX, ER;

    // check for xn or gamma (z=0 n>=0 .. these are ok)
    if (KZ == 0 && KA >= 0) goto L10;
    // check for yp (z=a n=0 ... these are also ok)
    if (KZ == KA) goto L20;

    MASEXX(KZ, KA, EX, ER, NOTTAB);
    return 1.0e-3 * EX;

    // return for neutral particles (xn and gamma)
L10:
    MASEXX(0, 1, EX, ER, NOTTAB);
    return 1.0e-3 * KA * EX;

    // return for y-proton cluster energies
L20:
    MASEXX(1, 1, EX, ER, NOTTAB);
    return 1.0e-3 * KA * EX;
}

// ====================================================================
// MASEXX / GETAS shared data (lines 54-1766 of masstable.f)
// ====================================================================

static const int IZMAX = 119;

'''
    # Add LZLO array
    cpp += format_array_int('LZLO', lzlo) + '\n\n'

    # Add LZHI array
    cpp += format_array_int('LZHI', lzhi) + '\n\n'

    # Add XBUF array
    cpp += format_xbuf(xbuf) + '\n\n'

    cpp += '''// ABUF(i,j) = XBUF[(j-1)*3 + i]  (Fortran EQUIVALENCE mapping)
// Fortran: REAL*8 ABUF(3, 3179), EQUIVALENCE (ABUF(1,1),XBUF(1))
#define ABUF(i,j) XBUF[((j)-1)*3 + (i)]

// ====================================================================
// MASEXX subroutine (lines 1770-1813 of masstable.f)
// Returns mass excess (keV) and uncertainty from 2003 Mass Adjustment
// ====================================================================

void MASEXX(int IZ, int IA, double& EX, double& ER, int& IERR) {
    int IAA, IBA, NDX, IAT;

    // RESET ERROR FLAG
    IERR = 0;
    // INDEX IN LZLO-LZHI PNTR LIST
    NDX = IZ + 1;
    // TST FOR OUT OF RANGE
    if (NDX > IZMAX) goto L100;

    // COMPUTE THE ABUF LO-LIMIT
    IAA = LZLO[NDX];
    // COMPUTE THE ABUF HI-LIMIT
    IBA = LZHI[NDX];

    // LOOK FOR REQUIRED A-NUMBER
    for (int J = IAA; J <= IBA; J++) {
        // CONV TO INTEGER
        IAT = (int)(ABUF(1, J) + 0.5);
        // TST AGAINST REQUEST
        if (IAT == IA) {
            // SET MASS EXCESS
            EX = ABUF(2, J);
            // SET ASSOCIATED UNCERTAINTY
            ER = ABUF(3, J);
            return;
        }
    }

L100:
    IERR = 1;
    return;
}

// ====================================================================
// GETAS - ENTRY point in MASEXX (lines 1815-1825)
// Returns list of mass numbers available for element Z
// ====================================================================

void GETAS(int IZ, int& NUMAS, int* IAS) {
    NUMAS = 0;
    int IZZ = IZ + 1;
    if (IZZ < 1 || IZZ > IZMAX) return;
    int LZ = LZLO[IZZ] - 1;
    NUMAS = LZHI[IZZ] - LZ;
    for (int I = 1; I <= NUMAS; I++) {
        IAS[I] = (int)(ABUF(1, LZ + I) + 0.5);
    }
    return;
}

#undef ABUF

// ====================================================================
// AZCODE subroutine (lines 1828-2026 of masstable.f)
// Decodes string like "CA40" or "40CA" to Z, A, symbol
// ====================================================================

// Element symbols (Z=1..109)
static const int NUMELS = 109;
static const char* ELEMNT[NUMELS + 1] = {
    "",   // [0] unused
    "H ", "HE", "LI", "BE", "B ", "C ", "N ", "O ", "F ",
    "NE", "NA", "MG", "AL", "SI", "P ", "S ", "CL", "AR", "K ", "CA",
    "SC", "TI", "V ", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA",
    "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y ", "ZR", "NB", "MO",
    "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I ",
    "XE", "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD",
    "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W ", "RE",
    "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN",
    "FR", "RA", "AC", "TH", "PA", "U ", "NP", "PU", "AM", "CM", "BK",
    "CF", "ES", "FM", "MD", "NO", "LR", "RF", "HA", "NH", "NS", "HS",
    "MT"
};

// Special one-character cases
static const int NCASE = 7;
static const char* SCASE[NCASE + 1] = {
    "",    // [0] unused
    "G ", "N ", "P ", "D ", "T ", "H ", "A "
};
static const int IZSC[NCASE + 1] = { 0, 0, 0, 1, 1, 1, 2, 2 };
static const int IASC[NCASE + 1] = { 0, 0, 1, 1, 2, 3, 3, 4 };

// Internal implementation for AZCODE + CODEZ entry
static void AZCODE_impl(bool is_codez, char* SIN, int& IZ, int& IA,
                         char* SYMBOL, int& IRTN, char* ZSYM) {
    char S[9];  // 1-based: S[1]..S[8]
    int SVAL[6]; // 1-based: SVAL[1]..SVAL[5]
    char SYM1[3]; // SYM1[1], SYM1[2] as 1-based
    char SYM[3];  // 2-char symbol + null
    int IDCODE, J, IZZ;

    // ---- ENTRY CODEZ ----
    if (is_codez) {
        SYM[0] = ' '; SYM[1] = ' '; SYM[2] = '\\0';
        if (IZ < 0 || IZ > NUMELS) goto L490;
        SYM[0] = 'N'; SYM[1] = ' ';
        if (IZ == 0) goto L490;
        SYM[0] = ELEMNT[IZ][0];
        SYM[1] = ELEMNT[IZ][1];
    L490:
        ZSYM[0] = SYM[0];
        ZSYM[1] = SYM[1];
        return;
    }

    // ---- AZCODE main logic ----
    IDCODE = 0;
    IRTN = 0;
    SYM[0] = ' '; SYM[1] = ' '; SYM[2] = '\\0';
    SYMBOL[0] = ' '; SYMBOL[1] = ' ';
    J = 1;

    // Copy input and convert to upper case
    for (int I = 1; I <= 8; I++) {
        S[I] = SIN[I - 1];  // SIN is 0-based C array
        if (S[I] >= 'a' && S[I] <= 'z')
            S[I] = (char)(S[I] + 'A' - 'a');
    }

    // Compute SVAL for first 5 characters
    for (int I = 1; I <= 5; I++) {
        SVAL[I] = (int)S[I] - (int)'0';
    }

    // LOOK AT FIRST BYTE...IF NUMERIC..STRING IS OF FORM '40CA'
    if (SVAL[1] >= 0 && SVAL[1] <= 9) goto L200;

    // J POINTS TO A NON-NUMERIC SYMBOL...
L100:
    SYM[0] = S[J];
    J = J + 1;

    // CHECK SPECIAL CASES ON G,N,P,D,T,H,A
    if (IDCODE == 0 && IZ != 999 && S[J] == ' ') goto L902;

    // CHECK NEXT CHARACTER..IF NUMERIC OR BLANK..END OF SYMBOL FIELD
    if ((SVAL[J] >= 0 && SVAL[J] <= 9) || S[J] == ' ') goto L110;
    // ELSE THIS IS A TWO CHARACTER SYMBOL
    SYM[1] = S[J];
    J = J + 1;

L110:
    if (IDCODE == 0 && (SVAL[J] < 0 || SVAL[J] > 9)) goto L901;

    // IF ELEMENT SYMBOL IS AT END..CHECK THAT NEXT CHAR IS BLANK
    if (IDCODE == 1 && S[J] != ' ') goto L9003;

    // NOW THE ELEMENT SYMBOL IS IN SYM
    IDCODE = IDCODE + 1;
    // IF BOTH FIELDS DECODED...GO FIND 'IZ'
    if (IDCODE >= 2) goto L300;

    // DECODE ATOMIC MASS...J-TH CHARACTER IS GUARANTEED NUMERIC
L200:
    IA = SVAL[J];
L210:
    J = J + 1;
    // IF NEXT CHARACTER IS NUMERIC...PROCEED WITH 'IA' DECODE
    if (SVAL[J] >= 0 && SVAL[J] <= 9) goto L220;
    // IF 2-ND FIELD AND NEXT CHAR IS BLANK...GO WRAP UP
    if (IDCODE >= 1 && S[J] == ' ') goto L300;
    // IF 2-ND FIELD AND NEXT CHAR IS NON-NUMERIC...BUGGERED STRING
    if (IDCODE >= 1 && (SVAL[J] < 0 || SVAL[J] > 9)) goto L9003;
    // DONE WITH FIRST FIELD OF THE FORM '40CA', NOW DECODE SYM
    IDCODE = IDCODE + 1;
    goto L100;

L220:
    IA = 10 * IA + SVAL[J];
    goto L210;

    // HAVE SYM....GET ITS 'IZ'
L300:
    IZZ = 0;
    SYMBOL[0] = SYM[0];
    SYMBOL[1] = SYM[1];

    // IF IT IS NP WITH N <= 9 THEN IT IS N PROTONS
    if (IA <= 9 && SYM[0] == 'P' && SYM[1] == ' ') goto L380;
    // IF SYM='N'...AND IA<=9...IT IS NEUTRON...ELSE NITROGEN
    if (IA <= 9 && SYM[0] == 'N' && SYM[1] == ' ') goto L390;

    // LOOP ON ELEMNT TABLE TO GET Z OF SYM
    for (IZZ = 1; IZZ <= NUMELS; IZZ++) {
        if (SYM[0] == ELEMNT[IZZ][0] && SYM[1] == ELEMNT[IZZ][1])
            goto L390;
    }
    // FALL THRU MEANS SOME SNAFU
    goto L9002;

L380:
    IZ = IA;
    return;

L390:
    IZ = IZZ;
    return;

    // HANDLE SPECIAL CASES HERE
L902:
    SYMBOL[0] = SYM[0];
    SYMBOL[1] = SYM[1];
    for (int N = 1; N <= NCASE; N++) {
        if (SYM[0] == SCASE[N][0] && SYM[1] == SCASE[N][1]) {
            IA = IASC[N];
            IZ = IZSC[N];
            return;
        }
    }
    // IF NONE OF THESE...GO TO DECODE SYM WITH IA=-1 ON RETURN
    goto L901;

L901:
    IRTN = -1;
    // ONLY A SYMBOL WAS CODED...IZ IS RETURNED BUT IA=-1
    goto L300;

L9002:
    IRTN = -2;
    return;

L9003:
    IRTN = -3;
    return;
}

void AZCODE(char* SIN, int& IZ, int& IA, char* SYMBOL, int& IRTN) {
    AZCODE_impl(false, SIN, IZ, IA, SYMBOL, IRTN, nullptr);
}

void CODEZ(int IZ, char* ZSYM) {
    int dummy_ia = 0, dummy_irtn = 0;
    char dummy_sin[8] = {' ',' ',' ',' ',' ',' ',' ',' '};
    char dummy_symbol[2] = {' ',' '};
    int iz = IZ;
    AZCODE_impl(true, dummy_sin, iz, dummy_ia, dummy_symbol, dummy_irtn, ZSYM);
}
'''
    return cpp

def main():
    lines = read_fortran()

    # Extract LZLO
    lzlo = extract_lzlo(lines)
    print(f"LZLO: {len(lzlo)} values (expected 119)", file=sys.stderr)
    assert len(lzlo) == 119, f"Expected 119 LZLO values, got {len(lzlo)}"

    # Extract LZHI
    lzhi = extract_lzhi(lines)
    print(f"LZHI: {len(lzhi)} values (expected 119)", file=sys.stderr)
    assert len(lzhi) == 119, f"Expected 119 LZHI values, got {len(lzhi)}"

    # Extract XBUF
    xbuf = extract_xbuf(lines)
    print(f"XBUF: {len(xbuf)} values (expected 9537)", file=sys.stderr)
    if len(xbuf) != 9537:
        print(f"WARNING: Expected 9537 XBUF values, got {len(xbuf)}", file=sys.stderr)
        # Print first/last few values for debugging
        print(f"First 12: {xbuf[:12]}", file=sys.stderr)
        print(f"Last 12: {xbuf[-12:]}", file=sys.stderr)

    # Generate C++
    cpp = generate_cpp(lzlo, lzhi, xbuf)

    with open(OUTPUT_FILE, 'w') as f:
        f.write(cpp)

    print(f"Generated {OUTPUT_FILE}", file=sys.stderr)
    print(f"  LZLO: {len(lzlo)} values", file=sys.stderr)
    print(f"  LZHI: {len(lzhi)} values", file=sys.stderr)
    print(f"  XBUF: {len(xbuf)} values", file=sys.stderr)

if __name__ == '__main__':
    main()
