# Elastic Scattering — Theory and Calculation Steps

This document describes the complete mathematical framework and computational
steps for elastic scattering in Ptolemy, from the optical model potential to
the differential cross section.

---

## 1. The Optical Model Potential

The nuclear scattering is described by a complex optical model potential:

$$U(r) = V_C(r) + V_N(r) + i\,W(r) + V_{SO}(r)\,\boldsymbol{\ell}\cdot\boldsymbol{s} + i\,W_{SO}(r)\,\boldsymbol{\ell}\cdot\boldsymbol{s} + \frac{\hbar^2}{2\mu}\frac{\ell(\ell+1)}{r^2}$$

### 1.1 Woods-Saxon Real Central Potential

$$V_N(r) = -\frac{V}{1 + \exp\!\left(\frac{r - R}{a}\right)}$$

where $R = r_0 \cdot A_T^{1/3}$ (or $r_0(A_P^{1/3}+A_T^{1/3})$ depending on R0 convention).

A generalized power form is also available:

$$V_N(r) = -\frac{V}{\left[1 + \exp\!\left(\frac{r-R}{a}\right)\right]^p}$$

where $p$ is set by the `REALPOWER` keyword (default $p=1$).

### 1.2 Imaginary Volume Potential

$$W_V(r) = -\frac{W}{1 + \exp\!\left(\frac{r - R_I}{a_I}\right)}$$

### 1.3 Imaginary Surface (Derivative Woods-Saxon)

$$W_S(r) = -4\,W_S\,\frac{\exp\!\left(\frac{r-R_{SI}}{a_{SI}}\right)}{\left[1 + \exp\!\left(\frac{r-R_{SI}}{a_{SI}}\right)\right]^2}$$

The total imaginary potential is $W(r) = W_V(r) + W_S(r)$.

### 1.4 Spin-Orbit Potential (Thomas Form)

$$V_{SO}(r) = \left(\frac{\hbar}{m_\pi c}\right)^2 \frac{V_{SO}}{r}\,\frac{d}{dr}\left[\frac{1}{1+\exp\!\left(\frac{r-R_{SO}}{a_{SO}}\right)}\right]$$

The spin-orbit expectation value is:

$$\langle\boldsymbol{\ell}\cdot\boldsymbol{s}\rangle = \frac{1}{2}\left[j(j+1) - \ell(\ell+1) - s(s+1)\right]$$

In the code (WAVELJ), this is computed as:

$$\sigma\cdot\ell = \frac{\frac{1}{4}\left[J_P(J_P+2) - J_{SP}(J_{SP}+2)\right] - \ell(\ell+1)}{J_{SP}}$$

where $J_P = 2j$ and $J_{SP} = 2s$ are doubled quantum numbers.

### 1.5 Coulomb Potential

For two uniform charge spheres of radii $R_1$ and $R_2$ (Poling et al., 
Phys. Rev. C13, 648, 1976):

$$V_C(r) = \begin{cases}
\displaystyle\frac{Z_P Z_T e^2}{r} & r > R_1 + R_2 \\[8pt]
\displaystyle\frac{A + r(B + r(C + r(D + r(E + r^2 F))))}{r} & R_1 - R_2 < r < R_1 + R_2 \\[8pt]
X + Y r^2 & r < R_1 - R_2
\end{cases}$$

where coefficients $A$–$F$, $X$, $Y$ are computed from $R_1$, $R_2$, and
$V_{C0} = Z_P Z_T e^2$. For point-like projectiles ($R_2 = 0$), this
reduces to the standard uniform-sphere formula:

$$V_C(r) = \begin{cases}
\displaystyle\frac{Z_P Z_T e^2}{r} & r > R_C \\[8pt]
\displaystyle\frac{Z_P Z_T e^2}{2R_C}\left(3 - \frac{r^2}{R_C^2}\right) & r < R_C
\end{cases}$$

### 1.6 Centrifugal Term

$$V_{\text{cent}}(r) = \frac{\hbar^2}{2\mu}\frac{\ell(\ell+1)}{r^2}$$

---

## 2. Radial Schrödinger Equation

The radial wave function $u_\ell(r)$ satisfies:

$$\frac{d^2 u_\ell}{dr^2} + \left[k^2 - \frac{\ell(\ell+1)}{r^2} - \frac{2\mu}{\hbar^2}\left(V_C + V_N + iW + (V_{SO}+iW_{SO})\,\boldsymbol{\ell}\cdot\boldsymbol{s}\right)\right] u_\ell = 0$$

where $k = \sqrt{2\mu E}/\hbar$ is the wave number.

---

## 3. Numerov Integration

Ptolemy solves the radial equation using the **Numerov method** (also called
Störmer-Verlet). The equation is rewritten as:

$$u''(r) = f(r)\,u(r)$$

where $f(r) = \frac{2\mu}{\hbar^2}[U(r) - E] + \frac{\ell(\ell+1)}{r^2}$.

### 3.1 Modified Potential Storage

The code stores the potential in a modified form that absorbs the Numerov
correction factor:

$$W_i = 1 - \frac{h^2}{12}\left(\frac{U(r_i)}{E} - 1\right)$$

where $h$ is the step size. For the centrifugal term:

$$C_i = -\frac{h^2}{12}\cdot\frac{1}{r_i^2}$$

The total effective potential at grid point $i$ is:

$$W_i^{\text{eff}} = W_i^{\text{real}} + \ell(\ell+1)\cdot C_i$$

### 3.2 The Three-Point Recursion

Define the auxiliary variable $\xi_i = \frac{1}{12}W_i\,u_i$. The Numerov
recursion is:

$$\xi_{i+1} = -10\,\xi_i + \xi_{i-1} + u_{i-1}$$

$$u_{i+1} = \frac{W_i\,\xi_{i+1}}{|W_i|^2} \cdot 12$$

For complex potentials, $W_i$ and $u_i$ are both complex, and the division
is:

$$u_{i+1} = 12\,\frac{\text{Re}(W_i)\,\text{Re}(\xi_{i+1}) + \text{Im}(W_i)\,\text{Im}(\xi_{i+1})}{|W_i|^2} + i\,\frac{\cdots}{|W_i|^2}$$

### 3.3 Starting Conditions

At the origin, the regular solution behaves as:

$$u_\ell(r) \propto r^{\ell+1} \quad (r \to 0)$$

The code starts at $r = r_{\text{start}} \cdot h$ with:

$$u(r_1) = C \cdot \left(\frac{r_1}{r_1 + h}\right)^{\ell+1}, \qquad u(r_2) = C$$

where $C$ is an arbitrary normalization constant (the absolute normalization
cancels in the S-matrix).

---

## 4. Asymptotic Matching and S-Matrix

### 4.1 Coulomb Wave Functions

At the matching radius $r = R_{\text{match}}$ (= ASYMPTOPIA), the nuclear
potential is negligible and the wave function is a linear combination of
regular and irregular Coulomb functions:

$$u_\ell(r) = A\,F_\ell(\eta, kr) + B\,G_\ell(\eta, kr)$$

where $\eta = Z_P Z_T e^2 \mu / (\hbar^2 k)$ is the Sommerfeld parameter.

Ptolemy computes $F_\ell$ and $G_\ell$ using the RCWFN subroutine
(continued-fraction method of Barnett, Feng, Steed, and Goldfarb,
Comp. Phys. Comm. **8**, 377, 1974).

### 4.2 S-Matrix Extraction

The S-matrix is extracted by matching the numerical wave function at two
points: $r_1 = R_{\text{match}}$ and $r_2 = R_{\text{match}} - N_{\text{back}}\cdot h$:

$$\begin{pmatrix} u(r_1) \\ u(r_2) \end{pmatrix} = \begin{pmatrix} F_1 & G_1 \\ F_2 & G_2 \end{pmatrix} \begin{pmatrix} A \\ B \end{pmatrix}$$

where $F_1 = F_\ell(\eta, kr_1)$, $G_1 = G_\ell(\eta, kr_1)$, etc.

Solving for $A$ and $B$:

$$A_1 = \text{Re}[u(r_1)]\,F_2 + \text{Im}[u(r_1)]\,G_2 - \text{Re}[u(r_2)]\,F_1 - \text{Im}[u(r_2)]\,G_1$$

$$A_2 = -\text{Re}[u(r_1)]\,G_2 + \text{Im}[u(r_1)]\,F_2 + \text{Re}[u(r_2)]\,G_1 - \text{Im}[u(r_2)]\,F_1$$

The S-matrix element is:

$$S_\ell = \frac{C_R + i\,C_I}{A_1 + i\,A_2}$$

where $C_R$ and $C_I$ are constructed similarly from the outgoing-wave
boundary condition.

In terms of the phase shift $\delta_\ell$:

$$S_\ell = e^{2i\delta_\ell}$$

For complex potentials, $\delta_\ell$ is complex and $|S_\ell| < 1$.

### 4.3 Spin-Orbit Coupling

For particles with spin $s > 0$, the S-matrix $S_{lj}$ depends on both
$\ell$ and $j = \ell \pm s$. The physical S-matrix in the $(\ell, L_X)$
basis is obtained via Racah recoupling:

$$S_{L_X}(\ell) = \sum_{j} \sqrt{2L_X+1}\; W(\ell\,\ell'\,s\,s;\,L_X\,j)\; S_{lj}$$

where $W$ is a Racah coefficient (6j symbol).

---

## 5. Scattering Amplitude

### 5.1 Nuclear Scattering Amplitude

The nuclear scattering amplitude for unpolarized particles is:

$$f_N(\theta) = \frac{1}{2ik} \sum_\ell (2\ell+1)\left(S_\ell - 1\right) P_\ell(\cos\theta)$$

For particles with spin, the amplitude generalizes to include spin-flip
and non-spin-flip components summed over $j$.

### 5.2 Coulomb Scattering Amplitude

The point-Coulomb scattering amplitude is:

$$f_C(\theta) = -\frac{\eta}{2k\sin^2(\theta/2)} \exp\!\left[-i\eta\ln\sin^2(\theta/2) + 2i\sigma_0\right]$$

where $\sigma_\ell = \arg\Gamma(\ell+1+i\eta)$ is the Coulomb phase shift.

### 5.3 Total Amplitude

The total elastic amplitude is:

$$f(\theta) = f_C(\theta) + \frac{1}{2ik}\sum_\ell (2\ell+1)\left(S_\ell\,e^{2i\sigma_\ell} - 1\right) P_\ell(\cos\theta)$$

Note that the Coulomb phase shifts $e^{2i\sigma_\ell}$ multiply the nuclear
S-matrix elements.

---

## 6. Differential Cross Section

### 6.1 Elastic Cross Section

$$\frac{d\sigma}{d\Omega}(\theta) = |f(\theta)|^2$$

In millibarns per steradian (Ptolemy's default unit):

$$\frac{d\sigma}{d\Omega} = 10 \times |f(\theta)|^2 \quad [\text{mb/sr}]$$

### 6.2 Rutherford Cross Section

$$\sigma_R(\theta) = \left(\frac{\eta}{2k\sin^2(\theta/2)}\right)^2 = \left(\frac{a}{2\sin^2(\theta/2)}\right)^2$$

where $a = \eta/k = Z_P Z_T e^2/(2E_{\text{CM}})$ is half the distance of
closest approach.

### 6.3 Ratio to Rutherford

$$\frac{\sigma}{\sigma_R}(\theta) = \frac{|f(\theta)|^2}{|f_C(\theta)|^2}$$

This ratio is the primary output of elastic scattering calculations.

### 6.4 Identical Particles

For identical bosons (e.g., $^{12}$C + $^{12}$C):

$$\frac{d\sigma}{d\Omega} = |f(\theta) + f(\pi-\theta)|^2$$

For identical fermions:

$$\frac{d\sigma}{d\Omega} = |f(\theta)|^2 + |f(\pi-\theta)|^2 - \text{Re}[f(\theta)f^*(\pi-\theta)]$$

with the interference term depending on the spin statistics.

---

## 7. WKB Approximation for $\ell_{\min}$

### 7.1 Classical Turning Point

The classical turning point $r_t$ is determined by:

$$E_{\text{CM}} = V_C(r_t) + \frac{\hbar^2\ell(\ell+1)}{2\mu r_t^2}$$

For a pure Coulomb potential:

$$\rho_t = \eta + \sqrt{\eta^2 + \ell(\ell+1) + \tfrac{1}{4}}$$

where $\rho = kr$.

### 7.2 Automatic $\ell$ Range Selection

Ptolemy uses the WKB phase shift to determine the partial wave range:

- **$\ell_{\min}$**: chosen where the absorption factor
  $|S_\ell|^2 \approx e^{-4\,\text{Im}(\delta_\ell)} \lesssim 0.05$ 
  (strong absorption)
- **$\ell_{\max}$**: chosen where $|1 - S_\ell| < \epsilon$ 
  (scattering negligible)

The default automatic selection uses:

$$\ell_{\min} \approx 0.2\,\rho_t - 4, \qquad \ell_{\max} \approx 1.5\,\rho_t$$

---

## 8. Calculation Flow Summary

```
Input: V, W, VSO, R0, a, RC0, ELAB, Z, A
          │
          ▼
┌─────────────────────┐
│ 1. MAKPOT            │  Build W_i = 1 - (h²/12)(U_i/E - 1)
│    Woods-Saxon +     │  for real, imaginary, spin-orbit,
│    Coulomb + cent.   │  centrifugal, and Coulomb potentials
└─────────┬───────────┘
          │
          ▼
┌─────────────────────┐
│ 2. RCWFN             │  Compute F_ℓ(η,kR), G_ℓ(η,kR)
│    Coulomb functions │  at matching radius R = asymptopia
│    at matching point │  for ℓ = ℓ_min ... ℓ_max
└─────────┬───────────┘
          │
          ▼
┌─────────────────────┐
│ 3. WAVELJ            │  For each ℓ (and j if spin-orbit):
│    Numerov from r≈0  │    • Start: u ~ r^(ℓ+1)
│    to R_match        │    • Step: 3-point Numerov recursion
│                      │    • Match: extract S_ℓ from u(R), u(R-h)
└─────────┬───────────┘
          │
          ▼
┌─────────────────────┐
│ 4. JPTOLX            │  If spin-orbit: recouple S_{ℓj}
│    Spin recoupling   │  to S_{L_X}(ℓ) using Racah coefficients
└─────────┬───────────┘
          │
          ▼
┌─────────────────────┐
│ 5. AMPCAL / XSECTN  │  f(θ) = f_C(θ) + Σ_ℓ (S_ℓ e^{2iσ_ℓ} - 1) P_ℓ
│    Partial wave sum  │  σ(θ) = 10 × |f(θ)|²  [mb/sr]
│    → dσ/dΩ(θ)       │  σ/σ_R = |f|² / |f_C|²
└─────────────────────┘
```

---

## 9. Key Physical Constants

| Constant | Symbol | Value |
|---|---|---|
| $\hbar c$ | HBARC | 197.3286 MeV·fm |
| Fine structure | $\alpha$ | 1/137.036 |
| Pion Compton wavelength | $\hbar/m_\pi c$ | 1.414 fm |
| AMU to MeV/c² | | 931.494 |

---

## References

- G.R. Satchler, *Direct Nuclear Reactions*, Oxford University Press (1983)
- M.H. Macfarlane and S.C. Pieper, Argonne report ANL-76-11 Rev.1 (1978)
- A.R. Barnett, D.H. Feng, J.W. Steed, L.J.B. Goldfarb, Comp. Phys. Comm. **8**, 377 (1974) — Coulomb functions
- J.R. Poling, E. Norbeck, R.R. Carlson, Phys. Rev. C **13**, 648 (1976) — folded Coulomb potential
