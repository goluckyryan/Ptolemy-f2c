# Inelastic Scattering (DWBA) — Theory and Calculation Steps

This document describes the collective-model DWBA for inelastic excitation
as implemented in Ptolemy. The target nucleus is excited from its ground
state to a rotational or vibrational state by the deformed nuclear and
Coulomb potentials.

---

## 1. The Inelastic Excitation Process

In inelastic scattering $A(a,a')A^*$, the projectile $a$ excites the
target $A$ from its ground state ($J_i^\pi$) to an excited state
($J_f^\pi$, excitation energy $E_x$) without transferring nucleons:

$$A(a,a')A^*(J_f^\pi,\;E_x)$$

The projectile loses kinetic energy equal to $E_x$ (Q-value = $-E_x$).

---

## 2. Collective Model — Deformed Potential

### 2.1 Nuclear Deformation

In the collective (vibrational/rotational) model, the nuclear surface is
deformed:

$$R(\theta,\phi) = R_0\left[1 + \sum_{\lambda\mu}\beta_\lambda\,Y_{\lambda\mu}(\theta,\phi)\right]$$

where $\beta_\lambda$ is the deformation parameter for multipolarity
$\lambda$ (set by the `BELX` keyword).

### 2.2 Deformed Optical Potential

The optical potential follows the nuclear shape. Expanding to first order
in $\beta_\lambda$:

$$U(r, \theta) = U(r) + \delta U(r, \theta)$$

$$\delta U(r, \theta) = \sum_\lambda \beta_\lambda R_0 \frac{dU}{dr}\,Y_{\lambda 0}(\theta)$$

The coupling interaction that drives the inelastic transition is:

$$V_{\text{coup}}(r) = \beta_\lambda\,R_0\,\frac{dV}{dr}$$

where $dV/dr$ is the derivative of the Woods-Saxon potential.

### 2.3 Form Factor

The nuclear form factor for the transition is:

$$H_{\text{nuc}}(r) = \frac{12E}{h^2}\,R_0\,\frac{dV}{dr}\bigg|_r \cdot w(r)$$

where $w(r)$ is the Gauss quadrature weight and $h$ is the step size
(the factor $12E/h^2$ converts from Ptolemy's internal
$[1 - (h^2/12)(V/E)]$ storage).

In the code (INGRST), this is computed via cubic spline interpolation of
the derivative:

$$H_{\text{nuc}}(r_i) = \frac{12E}{h^2}\,R_{\text{real}}\,w_i \cdot \left[B_n + x(2C_n + 3xD_n)\right]$$

where $B$, $C$, $D$ are spline coefficients and $x = r_i - r_n$ is the
offset from the nearest grid point.

### 2.4 Coulomb Form Factor

The Coulomb excitation form factor is:

$$H_{\text{coul}}(r) = -w(r) \cdot V_C \cdot \begin{cases}
\displaystyle\frac{1}{r^{\lambda+1}} & r \geq R_C \\[8pt]
\displaystyle\frac{r^\lambda}{R_C^{2\lambda+1}} & r < R_C
\end{cases}$$

where $V_C = -3 Z_P Z_T e^2 / (\lambda + 1)$ is the Coulomb coupling
constant (stored as `R2S(4)` in the code).

---

## 3. DWBA Transition Amplitude

### 3.1 The Inelastic T-Matrix

The DWBA transition amplitude for inelastic excitation is:

$$T_{fi}(\theta) = \left\langle \chi_f^{(-)}(k_f, r)\,J_f \left| V_{\text{coup}}(r) \right| J_i\,\chi_i^{(+)}(k_i, r) \right\rangle$$

where $\chi_i^{(+)}$ and $\chi_f^{(-)}$ are distorted waves in the entrance
and exit channels (same optical potential for both since the same
projectile scatters in both channels).

### 3.2 Partial Wave Expansion

Expanding in partial waves:

$$T_{fi}(\theta) = \sum_{\ell_i,\ell_f,L_X} C_{\ell_i \ell_f}^{L_X}\;I(\ell_i, \ell_f, L_X)\;P_{L_X}(\cos\theta)$$

where:
- $L_X = \lambda$ is the multipolarity of the excitation
- $I(\ell_i, \ell_f, L_X)$ is the radial integral (the "I-matrix")
- $C_{\ell_i \ell_f}^{L_X}$ includes Clebsch-Gordan coefficients

### 3.3 Selection Rules

$$|\ell_i - L_X| \leq \ell_f \leq \ell_i + L_X$$

$$\pi_i \cdot \pi_f = (-1)^{L_X}$$

For even-even target nuclei ($J_i = 0^+$), excitation of a $2^+$ state
requires $L_X = 2$ and $\Delta\ell = 0, \pm 2$.

---

## 4. The Radial I-Matrix

### 4.1 Numerical Integration

The I-matrix element is:

$$I(\ell_i, \ell_f, L_X) = \int_0^{R_{\max}} \chi_{\ell_f}^*(k_f, r)\;\left[\beta_N\,H_{\text{nuc}}(r) + \beta_C\,H_{\text{coul}}(r)\right]\;\chi_{\ell_i}(k_i, r)\;dr$$

This integral is split into two regions:

1. **Numerical region** ($0$ to $R_{\text{sum}}$): computed by direct
   summation over the radial grid
2. **Coulomb tail** ($R_{\text{sum}}$ to $\infty$): computed analytically
   using Coulomb function integrals (COULST/COULIN)

### 4.2 The Nuclear Contribution

$$I_{\text{nuc}} = \int_0^{R_{\max}} \chi_{\ell_f}^*(r)\;H_{\text{nuc}}(r)\;\chi_{\ell_i}(r)\;dr$$

Computed by summation over Gauss quadrature points:

$$I_{\text{nuc}} \approx \sum_{i=1}^{N_{\text{pts}}} \chi_{\ell_f}^*(r_i)\;H_{\text{nuc}}(r_i)\;\chi_{\ell_i}(r_i)$$

The quadrature points and weights come from CUBMAP (cubic-sinh mapped
Gauss-Legendre).

### 4.3 The Coulomb Tail Correction

Beyond $R_{\text{sum}}$ (= SUMMAX), the nuclear potential vanishes and
only the Coulomb excitation contributes. The distorted waves become
pure Coulomb functions:

$$\chi_\ell(r) \xrightarrow{r > R_{\text{sum}}} \frac{1}{2}\left[(1+S_\ell)\,F_\ell(\eta,kr) + i(1-S_\ell)\,G_\ell(\eta,kr)\right]$$

The tail integral becomes:

$$I_{\text{tail}}(\ell_i, \ell_f) = \frac{1}{4}\,C\,\beta_C \left[(1+S_f)(1+S_i)\,\mathcal{F}\mathcal{F} - (1-S_f)(1-S_i)\,\mathcal{G}\mathcal{G} + \cdots\right]$$

where the Coulomb function integrals are:

$$\mathcal{F}\mathcal{F}(\ell_i, \ell_f) = \int_{R_{\text{sum}}}^{\infty} \frac{F_{\ell_f}(\eta_f, k_f r)\;F_{\ell_i}(\eta_i, k_i r)}{r^{L_X+1}}\;dr$$

and similarly for $\mathcal{F}\mathcal{G}$, $\mathcal{G}\mathcal{F}$,
$\mathcal{G}\mathcal{G}$.

### 4.4 Pure Coulomb Subtraction

The total I-matrix with Coulomb correction is:

$$I_{\text{total}} = I_{\text{numerical}}(0 \to R_{\text{sum}}) + I_{\text{tail}}(R_{\text{sum}} \to \infty)$$

The **nuclear-only** amplitude (what Ptolemy reports as "NUCLEAR PART") is
obtained by subtracting the pure Coulomb excitation amplitude:

$$I_{\text{nuclear}} = I_{\text{total}} - I_{\text{pure Coulomb}}$$

where the pure Coulomb integral runs from $0$ to $\infty$:

$$I_{\text{pure Coulomb}}(\ell_i, \ell_f) = C\,\beta_C \int_0^{\infty} \frac{F_{\ell_f}(\eta_f, k_f r)\;F_{\ell_i}(\eta_i, k_i r)}{r^{L_X+1}}\;dr$$

---

## 5. Coulomb Function Integrals (COULIN / CLINTS)

### 5.1 The Recursion Method

The integrals $\mathcal{F}\mathcal{F}$, $\mathcal{F}\mathcal{G}$, etc.
are computed by COULIN using an upward recursion in $\ell$. The recursion
relates integrals at $(\ell_i, \ell_f)$ to those at lower $\ell$ values
using the Coulomb function recurrence relations.

The starting values at the lowest $\ell$ are computed by CLINTS using
the Belling asymptotic expansion method.

### 5.2 Belling's Method (CLINTS)

Each Coulomb function can be written as:

$$F_\ell(\eta, \rho) = \frac{\sin\Phi(\rho)}{\sqrt{Z(\rho)}}, \qquad G_\ell(\eta, \rho) = \frac{\cos\Phi(\rho)}{\sqrt{Z(\rho)}}$$

where $Z(\rho)$ is an amplitude function and $\Phi(\rho)$ is a phase
function, both computed by the RCASYM asymptotic expansion.

The product integrals then decompose into slowly and rapidly oscillating
parts:

$$\int \frac{F_f \cdot F_i}{r^N}\,dr = \frac{1}{2}\int\frac{\cos(\Phi_i - \Phi_f) - \cos(\Phi_i + \Phi_f)}{\sqrt{Z_i Z_f}\;r^N}\,dr$$

- **Difference piece** ($\Phi_i - \Phi_f$): slowly oscillating, computed
  numerically by Gauss-Legendre quadrature
- **Sum piece** ($\Phi_i + \Phi_f$): rapidly oscillating, computed by
  Belling's asymptotic expansion in powers of $1/r$

### 5.3 Integration Regions

The integration from $R$ to $\infty$ is divided into three regions:

1. **$R$ to $R_1$**: fully numerical (Gauss quadrature, NPCOULOMB points
   per Coulomb oscillation cycle)
2. **$R_1$ to $R_2$**: sum piece by Belling, difference piece numerical
3. **$R_2$ to $\infty$**: both pieces by Belling's expansion

The boundaries $R_1$, $R_2$ are chosen automatically based on the
convergence of the Belling series.

### 5.4 RCASYM — Asymptotic Coulomb Expansion

The amplitude $Z(\rho)$ and phase $\Phi(\rho)$ are expanded as power
series in $1/\rho$:

$$Z(\rho) = 1 + \sum_{n=2}^{N} z_n\left(\frac{\rho_0}{\rho}\right)^{n-1}$$

$$\Phi(\rho) = \rho - \eta\ln(2\rho) + \sigma_\ell - \frac{\pi\ell}{2} - \rho\sum_{n=2}^{N}\frac{z_n}{n-2}\left(\frac{\rho_0}{\rho}\right)^{n-1}$$

The coefficients $z_n$ are determined recursively from the differential
equation for $Z$. From $Z$ and $\Phi$, the Coulomb functions and their
derivatives are:

$$F = \frac{\sin\Phi}{\sqrt{Z}}, \qquad G = \frac{\cos\Phi}{\sqrt{Z}}$$

$$F' = Z\left(G - \frac{Z'}{2Z^2}F\right), \qquad G' = -Z\left(F + \frac{Z'}{2Z^2}G\right)$$

---

## 6. Beta Ratios and Coupling Strengths

### 6.1 Nuclear and Coulomb Deformation

In general, the nuclear and Coulomb deformation parameters can differ:

$$\beta_N = \texttt{BELX} \quad(\text{or } \texttt{BETA})$$

$$\beta_C = \texttt{BETACOUL} \quad(\text{defaults to } \beta_N \text{ if not set})$$

### 6.2 Beta Ratios

Ptolemy internally works with two ratios:
- **BETANRAT** ($\beta_N / \beta_{\text{ref}}$): nuclear beta ratio
- **BETARATS** ($\beta_C / \beta_{\text{ref}}$): Coulomb beta ratio

These multiply the nuclear and Coulomb form factors respectively in the
integration:

$$H_{\text{total}}(r) = \beta_N\,H_{\text{nuc}}(r) + \beta_C\,H_{\text{coul}}(r)$$

### 6.3 Radius Scaling

The nuclear coupling radius is:

$$R_{\text{real}} = r_0 \cdot A^{1/3} \quad \text{(real potential radius)}$$

$$R_{\text{imag}} = r_{i0} \cdot A^{1/3} \quad \text{(imaginary potential radius)}$$

$$R_C = r_{C0} \cdot A^{1/3} \quad \text{(Coulomb radius)}$$

The Coulomb coupling constant is:

$$V_C = \frac{3\,Z_P\,Z_T\,e^2}{\lambda + 1} = \frac{3\,Z_P\,Z_T\,\hbar c\,\alpha}{\lambda + 1}$$

stored internally as `R2S(4) = -3*ZP*ZT*HBARC/AFINE`.

---

## 7. Cross Section

### 7.1 Inelastic DWBA Cross Section

$$\frac{d\sigma}{d\Omega}(\theta) = \frac{\mu_i\mu_f}{(2\pi\hbar^2)^2}\frac{k_f}{k_i}\frac{1}{(2J_i+1)(2s+1)}\sum_{L_X}\left|\sum_{\ell_i,\ell_f}C_{\ell_i\ell_f}^{L_X}\;I(\ell_i,\ell_f,L_X)\;P_{L_X}(\cos\theta)\right|^2$$

where $s$ is the projectile spin.

### 7.2 Cancellation and Convergence

For Coulomb-dominated excitation (heavy ions, large $Z$), there is
massive cancellation between the nuclear and Coulomb contributions:

$$I_{\text{total}} \approx I_{\text{Coulomb}} + \delta I_{\text{nuclear}}$$

where $|\delta I_{\text{nuclear}}| \ll |I_{\text{Coulomb}}|$. This
cancellation can exceed 100 dB (factor $> 10^{10}$), requiring high
numerical precision. Ptolemy reports the cancellation in the S-matrix
output as a diagnostic.

---

## 8. Calculation Flow

```
Input: OMP, BELX, Ex, Jf, PARAMETERSET INELOCA
          |
          v
+------------------------+
| 1. WAVSET + MAKPOT     |  Build optical potentials for both channels
|    Incoming = Outgoing  |  (same projectile, different k due to Ex)
+----------+-------------+
           |
           v
+------------------------+
| 2. GETSCT              |  Compute distorted waves at Lmin and Lcrit
|    Scattering setup    |  via Numerov (WAVELJ) for grid determination
+----------+-------------+
           |
           v
+------------------------+
| 3. INGRST              |  Set up 1D integration grid (CUBMAP)
|    Form factor setup   |  Compute H_nuc(r): spline derivative of V(r)
|                        |  Compute H_coul(r): r^L / R_C^(2L+1) or 1/r^(L+1)
+----------+-------------+
           |
           v
+------------------------+
| 4. GRDSET              |  Full 2D integration grid (sum + diff + phi)
|    Grid + distorted    |  Compute all distorted waves chi_l(r)
|    waves               |  by Numerov for l = lmin ... lmax
+----------+-------------+
           |
           v
+------------------------+
| 5. COULST              |  Coulomb function integrals:
|    Coulomb integrals   |  FF, FG, GF, GG from R_sum to infinity
|    (COULIN + CLINTS)   |  via Belling expansion + recursion
|                        |  Pure Coulomb FF from 0 to infinity
+----------+-------------+
           |
           v
+------------------------+
| 6. INRDIN              |  For each (l_i, l_f, L_X):
|    Radial integral     |    I = sum over grid of chi_f* H chi_i
|    + Coulomb tail      |    + Coulomb tail (from COULST)
|    + pure Coulomb sub  |    - pure Coulomb (from COULST)
|    -> I-matrix         |    -> I_nuclear (phase * amplitude)
+----------+-------------+
           |
           v
+------------------------+
| 7. INELDC              |  Elastic S-matrix with unitarity:
|    Elastic + unitarity |  |S_el|^2 + sum|S_inel|^2 = 1
+----------+-------------+
           |
           v
+------------------------+
| 8. AMPCAL + XSECTN     |  Partial wave sum -> f(theta)
|    Cross section       |  dsigma/dOmega = |f|^2 [mb/sr]
+------------------------+
```

---

## 9. S-Matrix Output Format

Ptolemy prints for each $(\ell_i, \ell_f, L_X)$:

| Column | Symbol | Description |
|---|---|---|
| INTEGRAL(0, SUMMAX) | $\|I_{\text{num}}\|$ | Magnitude of numerical integral |
| Phase | $\arg(I_{\text{num}})$ | Phase of numerical integral |
| Cancellation | dB | $20\log_{10}(\sum\|H\|/\|I\|)$ — loss of significance |
| COULOMB | $I_{\text{pure Coul}}$ | Pure Coulomb amplitude (real) |
| NUCLEAR PART | $\|I_{\text{nuc}}\|$ | $\|I_{\text{total}} - I_{\text{pure Coul}}\|$ |
| TOTAL AMPLITUDE | $\|I_{\text{total}}\|$ | $\|I_{\text{num}} + I_{\text{tail}}\|$ |

The elastic channel also reports unitarity:

$$|S_{\text{elastic}}|^2 + \sum_{L_X}|S_{\text{inelastic}}|^2 = 1$$

---

## 10. Practical Considerations

### 10.1 PARAMETERSET Choice

| Case | Recommendation |
|---|---|
| Light ions on light targets | INELOCA1 (fast, sufficient) |
| Protons on medium/heavy targets | INELOCA2 |
| Alpha or heavy ions on heavy targets | INELOCA3 (high Coulomb cancellation) |

### 10.2 Additional Keywords for Convergence

For heavy systems with large $\eta$ (Coulomb parameter):

```
ASYMPTOPIA = 50      $ extend matching radius
LMAXADD = 60         $ more partial waves
NPCOULOMB = 12       $ more Coulomb integration points
```

### 10.3 Deformation Parameter

Typical $\beta_2$ values for the first $2^+$ state:

| Nucleus | $\beta_2$ | $E(2^+)$ [MeV] |
|---|---|---|
| $^{12}$C | 0.60 | 4.44 |
| $^{16}$O | 0.36 | 6.92 |
| $^{40}$Ca | 0.12 | 3.35 |
| $^{48}$Ca | 0.10 | 3.83 |
| $^{90}$Zr | 0.09 | 2.19 |
| $^{120}$Sn | 0.11 | 1.17 |
| $^{208}$Pb | 0.06 | 4.09 |

These can be extracted from measured $B(E2)$ values:

$$\beta_\lambda = \frac{4\pi}{3ZR_0^\lambda}\sqrt{\frac{B(E\lambda)\uparrow}{e^2}}$$

---

## References

- G.R. Satchler, *Direct Nuclear Reactions*, Oxford University Press (1983) — Ch. 5-6
- K. Alder et al., Rev. Mod. Phys. **28**, 432 (1956) — Coulomb excitation theory
- J.M. Belling, J. Phys. B **1**, 136 (1968) — asymptotic expansion method for Coulomb integrals
- M.H. Macfarlane and S.C. Pieper, Argonne report ANL-76-11 Rev.1 (1978) — Ptolemy manual
- S.C. Pieper (unpublished notes) — COULIN recursion and CLINTS implementation
- A. Bohr and B.R. Mottelson, *Nuclear Structure*, Vol. II (1975) — collective model
