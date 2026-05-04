# DWBA Transfer Reactions — Theory and Calculation Steps

This document describes the Distorted Wave Born Approximation (DWBA) for
nucleon transfer reactions as implemented in Ptolemy.

---

## 1. The Transfer Reaction

A transfer reaction $A(a,b)B$ involves:

| Particle | Description |
|---|---|
| $A$ | Target nucleus |
| $a$ | Projectile |
| $b$ | Ejectile |
| $B$ | Residual nucleus |
| $x$ | Transferred cluster ($a = b + x$ for stripping, $B = A + x$) |

**Stripping** ($a \to b + x$): the projectile loses nucleons to the target.
Example: $(d,p)$ — deuteron loses a neutron.

**Pickup** ($A \to B + x$): the projectile picks up nucleons from the target.
Example: $(p,d)$ — proton picks up a neutron.

Mass and charge conservation:

$$A_a + A_A = A_b + A_B, \qquad Z_a + Z_A = Z_b + Z_B$$

$$A_x = |A_a - A_b|, \qquad Z_x = |Z_a - Z_b|$$

---

## 2. DWBA Transition Amplitude

The DWBA transition amplitude for a stripping reaction $A(a,b)B$ is:

$$T_{DWBA} = \left\langle \chi_b^{(-)}(\boldsymbol{r}_\beta)\;\phi_B\;\phi_b \;\middle|\; V_{bx} \;\middle|\; \phi_A\;\phi_a\;\chi_a^{(+)}(\boldsymbol{r}_\alpha) \right\rangle$$

where:
- $\chi_a^{(+)}$ is the distorted wave in the entrance channel (outgoing wave boundary condition)
- $\chi_b^{(-)}$ is the distorted wave in the exit channel (incoming wave boundary condition)
- $\phi_a$, $\phi_b$, $\phi_A$, $\phi_B$ are internal wave functions
- $V_{bx}$ is the interaction that causes the transfer
- $\boldsymbol{r}_\alpha$, $\boldsymbol{r}_\beta$ are the relative coordinates in entrance and exit channels

### 2.1 Zero-Range Approximation

In the zero-range approximation, the interaction $V_{bx}$ and the projectile
internal wave function $\phi_a$ are replaced by:

$$\langle \phi_b | V_{bx} | \phi_a \rangle \approx D_0\,\delta(\boldsymbol{r}_{bx})$$

where $D_0$ is the zero-range normalization constant. This reduces the
six-dimensional integral to a three-dimensional one over a single radial
coordinate.

### 2.2 Finite-Range Correction

For finite-range calculations (non-zero range), Ptolemy uses the local
energy approximation (LEA). The finite-range correction factor modifies
the form factor:

$$D(r) = D_0 \cdot \exp\left(-\frac{\mu_{bx}}{2\hbar^2}\,R_{bx}^2\,\left[E_\alpha - U_\alpha(r) - E_\beta + U_\beta(r)\right]\right)$$

where $R_{bx}$ is the range parameter of the $b$-$x$ interaction.

---

## 3. Partial Wave Decomposition

### 3.1 Angular Momentum Coupling

The transition amplitude is expanded in partial waves. For a transferred
particle with orbital angular momentum $\ell_x$ and total angular momentum
$j_x$, the differential cross section involves:

$$\frac{d\sigma}{d\Omega}(\theta) = \frac{\mu_\alpha \mu_\beta}{(2\pi\hbar^2)^2}\frac{k_\beta}{k_\alpha}\frac{1}{(2J_A+1)(2J_a+1)} \sum_{M_A M_a M_B M_b} |T_{fi}|^2$$

### 3.2 Radial Integral

After angular momentum algebra, the transition amplitude reduces to a
one-dimensional radial integral:

$$\beta_{\ell_\alpha \ell_\beta}^{J} = \int_0^\infty \chi_{\ell_\beta}^{(-)*}(k_\beta, r)\;F_{\ell_x j_x}(r)\;\chi_{\ell_\alpha}^{(+)}(k_\alpha, r)\;dr$$

where $F_{\ell_x j_x}(r)$ is the radial form factor.

### 3.3 The Form Factor

The form factor is the overlap of the bound-state wave function with the
target/residual wave functions:

$$F_{\ell_x j_x}(r) = D_0\;\phi_{n\ell_x j_x}(r)$$

where $\phi_{n\ell_x j_x}(r)$ is the bound-state radial wave function of the
transferred particle $x$ in the target (for stripping) or projectile
(for pickup), computed in a Woods-Saxon potential well.

---

## 4. Bound-State Wave Function

### 4.1 Binding Potential

The transferred particle is bound in a Woods-Saxon potential:

$$V_B(r) = -V_0\,f(r) + V_{SO}\,\frac{1}{r}\frac{df}{dr}\,\boldsymbol{\ell}\cdot\boldsymbol{s} + V_C(r)$$

where $f(r) = \left[1 + \exp\left(\frac{r-R}{a}\right)\right]^{-1}$.

The depth $V_0$ is **automatically adjusted** by Ptolemy to reproduce the
correct binding energy:

$$E_B = \begin{cases}
S_n = M(A) - M(B) - M(x) & \text{(separation energy)} \\
Q + E_{\text{ex}}^B - E_{\text{ex}}^A & \text{(from Q-value)}
\end{cases}$$

### 4.2 Quantum Numbers

The bound-state wave function is specified by:
- $n$ — number of radial nodes (`NODES` keyword, $n_r = n - \ell - 1$ in shell-model notation)
- $\ell$ — orbital angular momentum (`L` keyword)
- $j$ — total angular momentum (`JP` keyword)

The node count follows the convention:

| Shell | $n_r$ (NODES) | $\ell$ | Example |
|---|---|---|---|
| 1s | 0 | 0 | $^{16}$O core |
| 1p | 0 | 1 | $^{16}$O valence |
| 2s | 1 | 0 | sd shell |
| 1d | 0 | 2 | sd shell |
| 1f | 0 | 3 | fp shell |
| 2p | 1 | 1 | fp shell |

### 4.3 Asymptotic Form

At large $r$, the bound-state wave function has the Whittaker form:

$$\phi(r) \xrightarrow{r\to\infty} C\;W_{-\eta_B,\ell+1/2}(2\kappa r)$$

where:
- $\kappa = \sqrt{2\mu_x E_B}/\hbar$ is the bound-state wave number
- $\eta_B = Z_x Z_B e^2 \mu_x / (\hbar^2 \kappa)$ is the bound-state Sommerfeld parameter
- $C$ is the asymptotic normalization coefficient (ANC)

### 4.4 Spectroscopic Factor

The physical cross section includes the spectroscopic factor $S$:

$$\left(\frac{d\sigma}{d\Omega}\right)_{\text{exp}} = S \cdot \left(\frac{d\sigma}{d\Omega}\right)_{\text{DWBA}}$$

Ptolemy computes the single-particle (DWBA) cross section; the user
extracts $S$ by comparison with data.

---

## 5. Distorted Waves

### 5.1 Entrance Channel

The distorted wave $\chi_{\ell_\alpha}^{(+)}(r)$ is the solution of:

$$\left[\frac{d^2}{dr^2} + k_\alpha^2 - \frac{\ell_\alpha(\ell_\alpha+1)}{r^2} - \frac{2\mu_\alpha}{\hbar^2}U_\alpha(r)\right]\chi_{\ell_\alpha}^{(+)}(r) = 0$$

with the boundary condition:

$$\chi_{\ell_\alpha}^{(+)}(r) \xrightarrow{r\to\infty} \frac{i}{2}\left[H_{\ell_\alpha}^{(-)}(\eta_\alpha, k_\alpha r) - S_{\ell_\alpha}\,H_{\ell_\alpha}^{(+)}(\eta_\alpha, k_\alpha r)\right]$$

where $H^{(\pm)} = G \pm iF$ are Coulomb-Hankel functions.

### 5.2 Exit Channel

Similarly, $\chi_{\ell_\beta}^{(-)}(r)$ satisfies the same equation with
exit-channel parameters $(k_\beta, \mu_\beta, U_\beta, \eta_\beta)$ and
incoming-wave boundary conditions.

### 5.3 Numerical Solution

Both distorted waves are computed by Numerov integration (same method as
elastic scattering, see [Elastic.md](Elastic.md)), using the optical
potentials specified in the `INCOMING` and `OUTGOING` blocks.

---

## 6. The DWBA Cross Section

### 6.1 General Formula

$$\frac{d\sigma}{d\Omega}(\theta) = \frac{\mu_\alpha \mu_\beta}{(2\pi\hbar^2)^2}\frac{k_\beta}{k_\alpha} \sum_{L_X} \hat{L}_X^2 \left|\sum_{\ell_\alpha \ell_\beta} C_{\ell_\alpha \ell_\beta}^{L_X}\;\beta_{\ell_\alpha \ell_\beta}^{L_X}\;P_{L_X}(\cos\theta)\right|^2$$

where:
- $\hat{L}_X = \sqrt{2L_X+1}$
- $C_{\ell_\alpha \ell_\beta}^{L_X}$ are angular momentum coupling coefficients (Clebsch-Gordan, Racah, 9j symbols)
- $\beta_{\ell_\alpha \ell_\beta}^{L_X}$ is the radial integral
- $P_{L_X}(\cos\theta)$ are Legendre polynomials

### 6.2 Selection Rules

Angular momentum and parity conservation require:

$$|\ell_\alpha - \ell_x| \leq \ell_\beta \leq \ell_\alpha + \ell_x$$

$$\pi_A \cdot \pi_a = \pi_B \cdot \pi_b \cdot (-1)^{\ell_x}$$

$$|J_A - j_x| \leq J_B \leq J_A + j_x$$

### 6.3 Angular Momentum Coupling Coefficients

The coupling coefficient involves:
- **Clebsch-Gordan coefficients**: $\langle \ell_\alpha\,0\;\ell_x\,0 | L_X\,0 \rangle$ (parity selection)
- **Racah W-coefficients**: $W(\ell_\alpha\,j_\alpha\,\ell_\beta\,j_\beta;\,s\,L_X)$ (spin recoupling)
- **9j symbols**: for tensor coupling with spin-orbit

---

## 7. Numerical Integration

### 7.1 Radial Integral Evaluation

The radial integral is computed on a one-dimensional grid using
Gauss-Legendre quadrature:

$$\beta = \sum_{i=1}^{N} w_i\;\chi_\beta^*(r_i)\;F(r_i)\;\chi_\alpha(r_i)$$

The grid uses a cubic-sinh mapping (CUBMAP) to concentrate points in
the nuclear surface region:

$$r = C + B\,\xi + A\,(\xi^2-1)(\xi+1)$$

where $\xi = \sinh(\tau\,x)/\gamma$ and $x \in [-1,1]$ is the Gauss
quadrature variable.

### 7.2 Integration Parameters

| Parameter | Keyword | Typical value |
|---|---|---|
| Grid points | `NPSUM` | 100–300 |
| Lower limit | `SUMMIN` | 0 fm |
| Upper limit | `SUMMAX` | 30–50 fm |
| Midpoint | `SUMMID` | $\sim R_{\text{nuclear}}$ |
| Mapping parameter | `GAMMASUM` | 5 |

### 7.3 Distorted Wave Interpolation

The distorted waves are computed on a uniform radial grid (step size $h$)
by Numerov integration. At the quadrature points $r_i$, the wave function
values are obtained by cubic spline interpolation (SPLNCB).

---

## 8. Two-Body Kinematics

### 8.1 Center-of-Mass Energy

$$E_{\text{CM}} = E_{\text{lab}}\,\frac{A_T}{A_P + A_T}$$

### 8.2 Q-Value

$$Q = (M_a + M_A - M_b - M_B)\,c^2$$

### 8.3 Wave Numbers

$$k_\alpha = \frac{\sqrt{2\mu_\alpha E_{\text{CM}}}}{\hbar}, \qquad k_\beta = \frac{\sqrt{2\mu_\beta (E_{\text{CM}} + Q)}}{\hbar}$$

### 8.4 Sommerfeld Parameters

$$\eta_\alpha = \frac{Z_a Z_A\,e^2\,\mu_\alpha}{\hbar^2\,k_\alpha}, \qquad \eta_\beta = \frac{Z_b Z_B\,e^2\,\mu_\beta}{\hbar^2\,k_\beta}$$

---

## 9. Special Cases

### 9.1 Neutron Transfer: $(d,p)$, $(p,d)$

The transferred particle $x = n$ has $Z_x = 0$, so there is no Coulomb
interaction in the bound state ($\eta_B = 0$). The projectile wave
function in the deuteron is typically computed using the AV18 potential
(`wavefunction av18`).

### 9.2 Proton Transfer: $(d,n)$, $(^3\text{He},d)$

The transferred particle $x = p$ has $Z_x = 1$, contributing a Coulomb
term to the bound-state potential. The outgoing neutron channel has
$Z_b = 0$, so there is no Coulomb distortion in the exit channel
(no `RC0` in outgoing OMP).

### 9.3 Two-Nucleon Transfer: $(p,t)$, $(t,p)$, $(^3\text{He},n)$

The transferred cluster $x$ has $A_x = 2$. In the DWBA, this is treated
as a single entity with appropriate quantum numbers. The form factor
$F(r)$ describes the two-nucleon cluster bound in the target/residual.

### 9.4 Heavy-Ion Transfer: $(^{12}\text{C},^{13}\text{C})$

One-nucleon transfer between heavy ions. The formalism is identical, but:
- Larger $\eta$ → more partial waves needed
- Larger $\ell_{\max}$ (use `LMAXADD`)
- Forward-angle peaking (use small `ANGLESTEP`, restrict `ANGLEMAX`)
- Both channels have strong absorption

### 9.5 Multi-Nucleon Transfer: $(^8\text{Li},\alpha)$

Transfer of a multi-nucleon cluster ($A_x \geq 2$). The form factor
describes the cluster bound state. The exchanged particle mass excess
can be specified with `MXCX`.

---

## 10. Calculation Flow

```
Input: Reaction, Bound-state params, OMPs
          │
          ▼
┌──────────────────────┐
│ 1. BOUND              │  Solve bound-state W.S. equation
│    (PROJECTILE block)  │  for φ_x in projectile: nodes, ℓ, j
│    (TARGET block)      │  for φ_x in target: nodes, ℓ, j
│                        │  Adjust V₀ to get correct E_B
└──────────┬─────────────┘
           │
           ▼
┌──────────────────────┐
│ 2. WAVSET + MAKPOT    │  Build optical potentials on radial grid
│    (INCOMING block)    │  Entrance channel: U_α(r)
│    (OUTGOING block)    │  Exit channel: U_β(r)
└──────────┬─────────────┘
           │
           ▼
┌──────────────────────┐
│ 3. GRDSET             │  Set up integration grids:
│    Integration grid    │  • 1D radial (CUBMAP, Gauss-Legendre)
│                        │  • Angular (Legendre polynomials)
│                        │  Compute distorted waves by Numerov
└──────────┬─────────────┘
           │
           ▼
┌──────────────────────┐
│ 4. WAVEF / WAVELJ     │  For each ℓ_α, j_α:
│    Distorted waves     │  • Numerov integration r=0 → R_match
│    + S-matrix          │  • Match to Coulomb functions
│                        │  • Extract S_ℓ (for elastic part)
│                        │  • Store χ(r) on grid
└──────────┬─────────────┘
           │
           ▼
┌──────────────────────┐
│ 5. LINTRP             │  Radial integral:
│    Form factor ×       │  β = ∫ χ_β*(r) F(r) χ_α(r) dr
│    distorted waves     │  via Gauss quadrature on CUBMAP grid
│    → β_{ℓα,ℓβ}        │  Angular momentum coupling coefficients
└──────────┬─────────────┘
           │
           ▼
┌──────────────────────┐
│ 6. AMPCAL             │  Scattering amplitude:
│    Partial wave sum    │  f(θ) = Σ C·β·P_L(cosθ)
│    + Coulomb           │  Include Coulomb phase shifts σ_ℓ
└──────────┬─────────────┘
           │
           ▼
┌──────────────────────┐
│ 7. XSECTN             │  dσ/dΩ(θ) = (μ_α μ_β)/(2πℏ²)² × (k_β/k_α)
│    Cross section       │              × Σ |f(θ)|²
│    + σ/σ_Ruth          │  Output: angles, σ_CM, σ_lab, σ/σ_R
└──────────────────────┘
```

---

## 11. Output

Ptolemy outputs for each angle:

| Column | Description |
|---|---|
| Angle (CM) | Center-of-mass scattering angle |
| Angle (Lab) | Lab scattering angle |
| $d\sigma/d\Omega$ (CM) | Differential cross section in CM (mb/sr) |
| $d\sigma/d\Omega$ (Lab) | Differential cross section in Lab (mb/sr) |

The cross section is the **single-particle** DWBA value. To compare with
experiment, multiply by the spectroscopic factor:

$$\left(\frac{d\sigma}{d\Omega}\right)_{\text{exp}} = \frac{(2J_B+1)}{(2J_A+1)}\;C^2 S\;\left(\frac{d\sigma}{d\Omega}\right)_{\text{DWBA}}$$

where $C$ is an isospin Clebsch-Gordan coefficient (often $C^2 = 1$ for
simple cases).

---

## References

- G.R. Satchler, *Direct Nuclear Reactions*, Oxford University Press (1983)
- M.H. Macfarlane and S.C. Pieper, Argonne report ANL-76-11 Rev.1 (1978) — Ptolemy manual
- N. Austern, *Direct Nuclear Reaction Theories*, Wiley (1970)
- P.D. Kunz, University of Colorado (unpublished) — DWUCK code documentation
- I.J. Thompson, *Nuclear Reactions for Astrophysics*, Cambridge (2009)
