# Ptolemy Output Format

This document describes the output format of Ptolemy cross-section
calculations and how to interpret the diagnostic columns.

---

## 1. Elastic Scattering Output

```
  ANGLE  C.M. MB        /RUTH     %/L   % ERROR  L>LMAX
```

### Columns

| Column | Description |
|---|---|
| ANGLE (C.M.) | Center-of-mass scattering angle (degrees) |
| ANGLE (LAB) | Lab scattering angle (degrees) |
| SIGMA/RUTH | Ratio of cross section to Rutherford |
| C.M. MB | Differential cross section $d\sigma/d\Omega$ in CM frame (mb/sr) |
| LAB | Differential cross section in lab frame (mb/sr) |
| RUTHERFORD | Rutherford cross section (mb/sr) |
| LOW L | % contribution from partial waves $\ell < \ell_{\min}$ |
| HIGH L | % contribution from partial waves $\ell > \ell_{\max}$ |

### Interpreting LOW L and HIGH L

These columns quantify how well the partial wave sum is converged:

$$\text{LOW L} = 50 \times \left(1 - \frac{\sigma_{\text{without } \ell < \ell_{\min}}}{\sigma_{\text{total}}}\right)$$

$$\text{HIGH L} = 50 \times \left(1 - \frac{\sigma_{\text{without } \ell > \ell_{\max}}}{\sigma_{\text{total}}}\right)$$

- **|LOW L| < 1%** and **|HIGH L| < 1%**: well converged
- **|HIGH L| > 1%**: increase `LMAX` or `LMAXADD`
- **|LOW L| > 1%**: decrease `LMIN` (rare)

---

## 2. Transfer (DWBA) Cross-Section Output

```
  ANGLE  C.M.  dsigma/dOmega    LAB     /RUTH
```

### Columns

| Column | Description |
|---|---|
| ANGLE (C.M.) | Center-of-mass angle (degrees) |
| ANGLE (LAB) | Lab angle (degrees) |
| dsigma/dOmega (C.M.) | DWBA differential cross section (mb/sr) |
| dsigma/dOmega (LAB) | Lab-frame cross section (mb/sr) |
| sigma/RUTH | Ratio to Rutherford (if applicable) |

### Note on normalization

The output is the **single-particle** DWBA cross section. To compare with
experiment:

$$\left(\frac{d\sigma}{d\Omega}\right)_{\text{exp}} = \frac{2J_f + 1}{2J_i + 1}\;C^2 S\;\left(\frac{d\sigma}{d\Omega}\right)_{\text{DWBA}}$$

where $C^2 S$ is the spectroscopic factor and $C$ is an isospin
Clebsch-Gordan coefficient.

---

## 3. Inelastic DWBA Output

### 3.1 S-Matrix Elements

```
  L   L  LX    INTEGRAL(0, SUMMAX)      COULOMB    NUCLEAR PART      TOTAL AMPLITUDE
  IN OUT    MAG.    PHASE  CANCELLATION  AMPLITUDE  MAG.  PHASE       MAG.  PHASE
```

| Column | Symbol | Description |
|---|---|---|
| L IN, L OUT | $\ell_i$, $\ell_f$ | Incoming and outgoing partial wave |
| LX | $L_X$ | Multipolarity of excitation |
| MAG (INTEGRAL) | $\|I_{\text{num}}\|$ | Magnitude of numerical integral (0 to SUMMAX) |
| PHASE (INTEGRAL) | $\arg(I_{\text{num}})$ | Phase of numerical integral (radians) |
| CANCELLATION | dB | Loss of significance: $20\log_{10}(\sum\|H\| / \|I\|)$ |
| COULOMB AMPLITUDE | $I_{\text{Coul}}$ | Pure Coulomb excitation amplitude (real) |
| NUCLEAR MAG | $\|I_{\text{nuc}}\|$ | $\|I_{\text{total}} - I_{\text{Coul}}\|$ |
| NUCLEAR PHASE | $\arg(I_{\text{nuc}})$ | Phase of nuclear amplitude |
| TOTAL MAG | $\|I_{\text{total}}\|$ | $\|I_{\text{num}} + I_{\text{tail}}\|$ |
| TOTAL PHASE | $\arg(I_{\text{total}})$ | Phase of total amplitude |

### Interpreting Cancellation

The CANCELLATION column (in dB) shows how much precision is lost in the
Coulomb correction:

- **< 50 dB**: minor cancellation, results reliable
- **50–100 dB**: moderate, use `PARAMETERSET INELOCA2` or higher
- **100–150 dB**: severe, use `INELOCA3` with large `NPCOULOMB`
- **> 150 dB**: extreme (heavy ions), verify with higher `ACCURACY`

### 3.2 Elastic Unitarity Check

```
  L   L' LX  MAGNITUDE   PHASE  COULOMB    ELASTIC  REACTION  RESIDUAL
```

| Column | Description |
|---|---|
| ELASTIC | $\|S_{\text{el}}\|^2$ |
| REACTION | $\sum\|S_{\text{inel}}\|^2$ |
| RESIDUAL | $1 - \text{ELASTIC} - \text{REACTION}$ (should be $\approx 0$) |

Unitarity requires ELASTIC + REACTION $\leq 1$. RESIDUAL should be
small and positive. A large RESIDUAL indicates numerical problems.

### 3.3 Cross-Section Table

```
  C.M.  REACTION     REACTION   LOW L  HIGH L   % FROM
  ANGLE  C.M. MB        /RUTH     %/L   % ERROR  L>LMAX
```

| Column | Description |
|---|---|
| ANGLE (C.M.) | Center-of-mass angle |
| C.M. MB | Inelastic differential cross section (mb/sr) |
| /RUTH | Ratio to Rutherford |
| %/L | Fractional contribution from each L |
| **% ERROR** | **Estimated error from partial wave truncation** |
| L>LMAX | Contribution from extrapolation beyond LMAX |

---

## 4. The % ERROR Column — Reliability Indicator

The `% ERROR` column is the most important diagnostic for judging whether
the cross section is reliable.

### How it's computed

Ptolemy estimates the truncation error by comparing the amplitude sum
including vs excluding the last few partial waves:

$$\text{\% ERROR} \approx 50 \times \left|1 - \frac{\sigma_{\text{truncated}}}{\sigma_{\text{full}}}\right|$$

The Wynn epsilon algorithm may be used to accelerate the sum, and the
error estimate includes the difference between the accelerated and
unaccelerated values.

### Interpretation

| % ERROR | Meaning | Action |
|---|---|---|
| < 0.1% | Excellent convergence | Results are reliable |
| 0.1–1% | Good convergence | Results usable |
| 1–10% | Marginal convergence | Increase `LMAX` or `LMAXADD` |
| 10–100% | Poor convergence | Results may be wrong; increase `LMAX`, check potentials |
| **100%** | **No convergence** | **Results are unreliable — do not use** |

### When % ERROR = 100%

This occurs when:
- Many partial waves failed to extrapolate ("CANNOT EXTRAPOLATE" errors)
- The partial wave sum did not converge
- The outgoing channel energy is very low (high Ex, low beam energy)

The cross sections are formally computed but **meaningless**. The
`**********` in the `/RUTH` column further indicates overflow/unreliable
values.

**Fix:** Use higher beam energy, lower excitation energy, or check that
the optical potentials are physically reasonable for the kinematics.

---

## 5. Common Error Messages

### "CANNOT EXTRAPOLATE FOR CHANNEL"

```
**** ERROR: CANNOT EXTRAPOLATE FOR CHANNEL 2 (JP,JT,LX,LO-LI) = (0/2 4/2 2 0)
```

The L-extrapolation failed for this specific spin channel. The S-matrix
element at high L could not be fit to the extrapolation function. This
channel's contribution beyond LMAX is set to zero.

**Common with:** deuteron (d,d') inelastic (spin-1 tensor channels),
high excitation energies, mismatched optical potentials.

### "PHASE OF S(L) HAS FLUCTUATING SIGN"

```
**** WARNING: PHASE OF S(L) FOR CHANNEL 2 HAS FLUCTUATING SIGN
```

The S-matrix phase changes sign between adjacent L values, preventing
smooth extrapolation. Usually accompanied by "CANNOT EXTRAPOLATE".

### "LXTRP1" or "MINIMIZATION FAILED"

```
*** MINIMIZATION W.R.T. WIDTH HAS FAILED
    EXTRAPOLATION MUST BE SUPPRESSED
```

The least-squares fit of the extrapolation function did not converge.
The L-extrapolation is disabled for this channel.

### "HEADER OF OBJECT N HAS BEEN DESTROYED"

An allocator memory diagnostic (C++ only). A freed object's bookkeeping
tag was overwritten. This is a cosmetic warning — it does not affect
computed cross sections. See [README.md](../README.md) for details.

---

## 6. Named Storage Output

### NSCATALOG

Lists all named objects in the allocator with their sizes and locations.
Useful for debugging memory usage.

### WRITENS4 / WRITENS8

Writes named arrays to `fort.N` files in binary format:
- `WRITENS4 SMAG` — writes S-matrix magnitudes as REAL*4
- `WRITENS8 SMAG` — writes as REAL*8

### DUMPALL / DUMP name

Prints array contents to stdout. `DUMPALL` prints everything;
`DUMP name` prints a specific array.

---

## 7. Cross-Section Extraction

The `extract_xsec.py` script extracts angle-vs-cross-section data from
Ptolemy output:

```bash
python3 extract_xsec.py output.out prefix
# Creates: prefix_xsec_1.dat, prefix_xsec_2.dat, ...
```

Each `.dat` file has two columns: angle (CM) and $d\sigma/d\Omega$ (mb/sr).

**Note:** The extraction script finds lines matching the Ptolemy cross-section
table format. It works for transfer and inelastic output. For elastic
scattering, read the sigma/Rutherford column directly from the output.

---

## 8. Recommended Practices

1. **Always check `% ERROR`** — if it's > 10%, the cross sections at that
   angle are unreliable.

2. **Check `HIGH L`** — if it's > 1%, increase `LMAXADD` or set explicit
   `LMAX`.

3. **Check `CANCELLATION`** (inelastic) — if it's > 100 dB, use
   `PARAMETERSET INELOCA3` with higher `ACCURACY` and `NPCOULOMB`.

4. **Watch for "CANNOT EXTRAPOLATE" errors** — a few are normal for
   deuteron reactions. Many (>20) indicate the calculation is unreliable.

5. **Use `PRINT=2`** to see S-matrix elements and diagnose convergence
   issues. Use `PRINT=999999` only for debugging.
