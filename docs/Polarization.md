# Polarization and Analyzing Power — Theory and Calculation

This document describes the mathematical framework for polarization observables
(analyzing powers) in Ptolemy, covering elastic scattering and DWBA reactions.
The primary reference is B.A. Robson, *The Theory of Polarization Phenomena*
(Clarendon Press, Oxford, 1974).

---

## 1. Scattering Amplitudes with Spin

### 1.1 Spin-Orbit Coupling and the S-Matrix

When the optical potential includes a spin-orbit term, the radial Schrodinger
equation must be solved for each value of the total angular momentum
$j = \ell \pm s_P$, where $s_P$ is the projectile spin and $\ell$ is the orbital
angular momentum. The S-matrix elements $S_{\ell,j}$ then carry both indices.

For a spin-1/2 projectile on a spinless target, the two channels per $\ell$ are
$j = \ell + 1/2$ and $j = \ell - 1/2$. Without spin-orbit, the S-matrix is
$j$-independent:

$$\frac{1}{2s_P + 1}\sum_j (2j+1) S_{\ell,j} = (2\ell+1) S_\ell$$

### 1.2 Elastic Scattering Amplitudes (Spin-1/2 + Spinless Target)

The scattering amplitude decomposes into a spin-independent part $F(\theta)$
and a spin-flip part $B(\theta)$:

$$\frac{d\sigma}{d\Omega} = |F(\theta)|^2 + |B(\theta)|^2$$

where

$$F(\theta) = F_R(\theta) + F_N(\theta)$$

$$F_R(\theta) = -\frac{\eta}{2k[\sin(\theta/2)]^2}\,
e^{2i[\sigma_0 - \eta\ln\sin(\theta/2)]}$$

$$F_N(\theta) = \frac{1}{2ik(2s_P+1)}
\sum_{\ell,j}(2j+1)(S_{\ell,j}-1)\,e^{2i\sigma_\ell}\,P_{\ell,0}(\cos\theta)$$

$$B(\theta) = \frac{1}{2ik}\sum_\ell
\bigl[S_{\ell,\ell+1/2} - S_{\ell,\ell-1/2}\bigr]\,
e^{2i\sigma_\ell}\,P_{\ell,1}(\cos\theta)$$

Here $\sigma_\ell = \arg\,\Gamma(\ell+1+i\eta)$ is the Coulomb phase shift and
$P_{\ell,m}$ are associated Legendre functions (Abramowitz and Stegun convention).

### 1.3 General Case: Multiple Amplitudes

For reactions $A(a,b)B$ with arbitrary spins, the scattering amplitude
is a matrix in spin space. Ptolemy labels each amplitude by a set of
quantum numbers $(L_x, M_x, J_P, J_T)$ through a table-of-contents array:

| Index | Quantum Number | Meaning |
|-------|---------------|---------|
| ITOC(1,k) | $L_{\text{out}} - L_{\text{in}}$ | Orbital angular momentum transfer |
| ITOC(2,k) | $L_x$ | Recoil/multipolarity rank |
| ITOC(3,k) | $J_P$ (×2) | Projectile total angular momentum |
| ITOC(4,k) | $J_T$ (×2) | Target total angular momentum |

The number of independent amplitudes at each angle is NSPL, the number of
valid $(L_x, M_x, J_P, J_T)$ combinations.

Each amplitude is computed by partial-wave summation (subroutine `AMPCAL`):

$$F_k(\theta) = \sum_{L=L_\min}^{L_\max}
P_{L,M_x}(\cos\theta) \cdot \beta_k(L)$$

where $\beta_k(L)$ encodes the S-matrix elements weighted by
Clebsch-Gordan coefficients and Coulomb phases (subroutine `BETCAL`).

---

## 2. The Mueller Matrix Formalism

### 2.1 Definition

All polarization observables for a nuclear reaction are encoded in the
**Mueller matrix** $Z$, which is a generalization of the density matrix
formalism. The Mueller matrix element
$Z(K_a, Q_a, K_b, Q_b, K_A, Q_A, K_B, Q_B)$ describes the correlation
between incoming polarization (ranks $K_a$, $K_A$ for projectile and target)
and outgoing polarization (ranks $K_b$, $K_B$ for ejectile and residual).

The eight indices represent spherical tensor ranks and components for the
four particles. The unpolarized cross section is the special case where
all ranks are zero: $Z(0,0,0,0,0,0,0,0) \propto d\sigma/d\Omega$.

### 2.2 Mueller Matrix in Terms of Amplitudes

The Mueller matrix is bilinear in the scattering amplitudes:

$$Z(K_a,Q_a,K_b,Q_b,K_A,Q_A,K_B,Q_B;\theta) =
\sum_{k,k'} C_{kk'}^{(K_a Q_a K_b Q_b K_A Q_A K_B Q_B)}\;
F_k(\theta)\,F_{k'}^*(\theta)$$

where the coefficients $C_{kk'}$ are angle-independent and computed from
angular momentum recoupling (subroutine `MUELCO`). The coefficients involve
products of three 9-j symbols and Clebsch-Gordan coefficients.

---

## 3. Analyzing Powers

### 3.1 Definition

The **analyzing power** $T_{KQ}$ of rank $K$ and component $Q$ describes the
sensitivity of the cross section to the polarization state of the
incoming beam. For a beam described by a density matrix $\rho$ expanded
in spherical tensor operators $\tau_{KQ}$:

$$\rho = \frac{1}{2s_a+1}\sum_{KQ} t_{KQ}\,\tau_{KQ}^\dagger$$

the cross section becomes:

$$\frac{d\sigma}{d\Omega}(\theta) =
\left(\frac{d\sigma}{d\Omega}\right)_0
\sum_{KQ} t_{KQ}\,T_{KQ}(\theta)$$

where $t_{KQ}$ are the beam polarization parameters, $T_{KQ}$ are the
analyzing powers, and $T_{00} = 1$ by normalization.

### 3.2 Relation to the Mueller Matrix

The analyzing powers are a special case of the Mueller matrix with all
outgoing and target indices set to zero:

$$T_{KQ}(\theta) = (-1)^Q\,
\frac{Z(K, -Q, 0, 0, 0, 0, 0, 0;\theta)}{Z(0, 0, 0, 0, 0, 0, 0, 0;\theta)}$$

Equivalently, using the symmetry $Z(K,Q) = (-1)^{K+Q} Z(K,-Q)^*$:

$$T_{KQ}(\theta) = (-1)^{K}\,
\frac{Z(K, Q, 0, 0, 0, 0, 0, 0;\theta)^*}{Z(0, 0, 0, 0, 0, 0, 0, 0;\theta)}$$

### 3.3 Naming Convention

The rank $K$ runs from 1 to $2s_a$ (twice the projectile spin).
For each $K$, the component $Q$ runs from $Q_{\min}$ to $K$, where:

- $Q_{\min} = 0$ for even $K$
- $Q_{\min} = 1$ for odd $K$

For odd $K$, the analyzing power $Z(K,Q,\ldots)$ is purely imaginary,
so Ptolemy outputs $i\cdot T_{KQ}$ (a real quantity).
The printed labels follow the convention:

| $K$ | $Q$ range | Labels |
|-----|-----------|--------|
| 1 | 1 | $iT_{11}$ |
| 2 | 0, 1, 2 | $T_{20}$, $T_{21}$, $T_{22}$ |
| 3 | 1, 2, 3 | $iT_{31}$, $iT_{32}$, $iT_{33}$ |
| 4 | 0, 1, 2, 3, 4 | $T_{40}$, ..., $T_{44}$ |

The total number of independent analyzing powers is:

$$N_{AP} = \frac{(2s_a + 2)(2s_a)}{2} = \frac{J_a(J_a+2)}{2}$$

where $J_a = 2s_a$ is the doubled spin.

### 3.4 Common Special Cases

**Spin-1/2 projectile** ($s_a = 1/2$, e.g., proton, neutron):

Only $K = 1$ exists, with $Q = 1$. The single observable is the **vector
analyzing power** $iT_{11}$, related to the conventional analyzing power
$A_y$ by:

$$A_y(\theta) = -\sqrt{2}\;iT_{11}(\theta)$$

In terms of the scattering amplitudes $F$ and $B$:

$$A_y(\theta) = \frac{2\,\mathrm{Im}[F(\theta)\,B^*(\theta)]}
{|F(\theta)|^2 + |B(\theta)|^2}$$

**Spin-1 projectile** ($s_a = 1$, e.g., deuteron):

Both $K = 1$ (vector) and $K = 2$ (tensor) exist. The observables are
$iT_{11}$, $T_{20}$, $T_{21}$, $T_{22}$.

The vector analyzing power $iT_{11}$ is related to $A_y$ by:

$$A_y(\theta) = -\sqrt{\frac{3}{2}}\;iT_{11}(\theta)$$

---

## 4. Computation of the Mueller Matrix Coefficients

### 4.1 Structure of the Calculation (MUELCO)

The coefficient of $F_k(\theta)\,F_{k'}^*(\theta)$ in the Mueller matrix is
computed by the subroutine `MUELCO`. For a given set of tensor ranks
$(K_a, K_b, K_A, K_B)$, all $Q$ values are computed simultaneously.

The calculation factorizes into three 9-j symbols:

**Target 9-j symbol** (coupling target spins to tensor rank $K_T$):

$$\mathcal{N}_T \cdot
\begin{Bmatrix}
J_A & J_T & J_B \\
K_A & K_T & K_B \\
J_A & J_T' & J_B
\end{Bmatrix}$$

where

$$\mathcal{N}_T = \sqrt{(2J_A+1)(2J_B+1)(2K_A+1)(2J_T+1)(2J_T'+1)(2K_T+1)}$$

**Projectile 9-j symbol** (coupling projectile spins to tensor rank $K_P$):

$$\mathcal{N}_P \cdot
\begin{Bmatrix}
J_b & J_P & J_a \\
K_b & K_P & K_a \\
J_b & J_P' & J_a
\end{Bmatrix}$$

where

$$\mathcal{N}_P = (2K_P+1)\sqrt{(2J_a+1)(2J_b+1)(2K_b+1)(2J_P+1)(2J_P'+1)}$$

**Orbital 9-j symbol** (coupling orbital angular momenta):

$$\mathcal{N}_L \cdot
\begin{Bmatrix}
L_x & J_P & J_T \\
K & K_P & K_T \\
L_x' & J_P' & J_T'
\end{Bmatrix}$$

where

$$\mathcal{N}_L = (2K+1)\sqrt{2L_x+1}$$

### 4.2 Assembly of the Coefficient

The full coefficient involves:

1. Looping over intermediate ranks $K_T$ and $K_P$, constrained by
   triangle rules from the 9-j symbols
2. For each $(K_T, K_P)$, multiplying the three 9-j symbols together
3. Summing over an intermediate rank $K$ with Clebsch-Gordan coefficients:

$$\text{COEF}(Q_a, k, k') = \sum_{K_T, K_P, K} \mathcal{N}_T \cdot
\text{9j}_T \cdot \mathcal{N}_P \cdot \text{9j}_P \cdot \mathcal{N}_L
\cdot \text{9j}_L \cdot C_{K_T K_P}^{K} \cdot C_{K_a K_b}^{K_P}$$

The Clebsch-Gordan coefficients enforce the coupling of $Q$ values:

$$C^{L_x, \; K, \; L_x'}_{M_x, \; Q_K, \; M_x'} \qquad
C^{K, \; K_P, \; K_T}_{Q_K, \; Q_P, \; Q_T} \qquad
C^{K_a, \; K_P, \; K_b}_{Q_a, \; Q_P-Q_b, \; Q_b} \qquad
C^{K_A, \; K_T, \; K_B}_{Q_A, \; Q_T-Q_B, \; Q_B}$$

### 4.3 Parity Phases

For amplitudes with negative $M_x$ or $M_x'$ (which are not stored
independently), the parity of the reaction enters:

$$F_{L_x, -M_x} = \pi \cdot (-1)^{L_x + M_x} \cdot F_{L_x, M_x}$$

where $\pi = \pm 1$ is the parity change of the reaction.

---

## 5. The ANAPOW Algorithm

The subroutine `ANAPOW` orchestrates the full analyzing power computation.

### 5.1 Overview

```
For K = 1 to 2*s_a:
    Q_min = K mod 2    (0 for even K, 1 for odd K)

    Call MUELCO to get coefficients C(Q, k, k') for all Q, k, k'

    For each amplitude pair (k, k'):
        If diagonal (k = k') and K = K_cross:
            Accumulate cross section: σ(θ) += T · |F_k(θ)|²

        Compute F_k · F_{k'}*:
            Even K: use Re[F_k · F_{k'}*]
            Odd K:  use Im[F_k · F_{k'}*]

        For each Q from Q_min to K:
            T_{KQ}(θ) += C(Q, k, k') · (Re or Im)[F_k · F_{k'}*]

Normalize: T_{KQ}(θ) = T_{KQ}(θ) / σ(θ)
```

### 5.2 Real vs. Imaginary Parts

The analyzing power involves products $F_k(\theta)\,F_{k'}^*(\theta)$.
Whether the real or imaginary part is needed depends on the parity of $K$:

- **Even K** (tensor): The Mueller matrix element is real.
  Use $\mathrm{Re}[F_k\,F_{k'}^*] = \mathrm{Re}(F_k)\,\mathrm{Re}(F_{k'})
  + \mathrm{Im}(F_k)\,\mathrm{Im}(F_{k'})$

- **Odd K** (vector): The Mueller matrix element is imaginary.
  Ptolemy outputs $i \cdot T_{KQ}$ which requires
  $\mathrm{Im}[F_k\,F_{k'}^*] = \mathrm{Re}(F_{k'})\,\mathrm{Im}(F_k)
  - \mathrm{Re}(F_k)\,\mathrm{Im}(F_{k'})$

This follows from the relation $i \cdot T_{KQ} = \mathrm{Im}[Z(K,Q,\ldots)]$
for odd $K$, because the factor $(-1)^K$ cancels the $i$.

### 5.3 Cross Section for Normalization

The cross section used for normalization is accumulated from the diagonal
terms ($k = k'$) during the $K = K_{\text{cross}}$ loop, where
$K_{\text{cross}} = \min(2s_a, 2)$:

$$\sigma(\theta) = \sum_k T_k \cdot |F_k(\theta)|^2$$

where $T_k = 2$ if $M_x \neq 0$ (since $\pm M_x$ amplitudes contribute
equally), and $T_k = 1$ if $M_x = 0$.

---

## 6. Spin-1/2 Elastic Scattering: Explicit Formulas

For the important special case of a spin-1/2 projectile (proton, neutron)
on a spin-0 target, there are exactly two amplitudes: the non-spin-flip
$F(\theta)$ and the spin-flip $B(\theta)$. The observables reduce to
closed-form expressions.

### 6.1 Differential Cross Section

$$\frac{d\sigma}{d\Omega} = |F(\theta)|^2 + |B(\theta)|^2$$

### 6.2 Vector Analyzing Power

$$iT_{11}(\theta) =
-\frac{\sqrt{2}\;\mathrm{Im}[F(\theta)\,B^*(\theta)]}
{|F(\theta)|^2 + |B(\theta)|^2}$$

Or equivalently:

$$A_y(\theta) = \frac{2\,\mathrm{Im}[F(\theta)\,B^*(\theta)]}
{|F(\theta)|^2 + |B(\theta)|^2}$$

### 6.3 Polarization Transfer (Outgoing Polarization)

For an initially unpolarized beam, the outgoing polarization is:

$$P(\theta) = A_y(\theta)$$

This is a consequence of time-reversal invariance for elastic scattering.

---

## 7. Spin-1 Elastic Scattering: Tensor Observables

For a spin-1 projectile (deuteron) on a spin-0 target, there are three
amplitudes corresponding to the three possible total angular momenta
$j = \ell - 1, \ell, \ell + 1$ for each orbital angular momentum $\ell$.

### 7.1 Tensor Analyzing Powers

The tensor analyzing powers $T_{20}$, $T_{21}$, $T_{22}$ describe
the sensitivity to the alignment (as opposed to orientation) of the beam.

For a pure tensor-polarized beam with alignment axis along $\hat{z}$
(beam direction), the cross section is:

$$\frac{d\sigma}{d\Omega} = \left(\frac{d\sigma}{d\Omega}\right)_0
\left[1 + \frac{1}{\sqrt{5}}\,t_{20}\,T_{20}(\theta)\right]$$

where $t_{20}$ characterizes the degree of tensor polarization.

### 7.2 Relation to Cartesian Analyzing Powers

The spherical tensor analyzing powers are related to the Cartesian
frame (Madison convention, $\hat{y}$ normal to scattering plane):

$$A_y = -i\sqrt{3}\;iT_{11}$$

$$A_{yy} = \sqrt{\frac{2}{3}}\;T_{20} + \sqrt{2}\;T_{22}$$

$$A_{xx} = \sqrt{\frac{2}{3}}\;T_{20} - \sqrt{2}\;T_{22}$$

$$A_{xz} = -\sqrt{2}\;T_{21}$$

---

## 8. Analyzing Powers in DWBA Reactions

### 8.1 Inelastic Excitation

For inelastic scattering $A(a,b)B$ where $b = a$ and $B = A^*$, the
analyzing powers are computed in the same Mueller matrix framework.
The amplitudes $G^{L_x}_{M_x}(\theta)$ (Eq. 2.38 of the manual)
carry the multipolarity index $L_x$ and its projection $M_x$, plus
channel-spin quantum numbers $(J_P, J_T)$.

The number of independent amplitudes NSPL grows with the spins of
the particles. For a spin-1/2 projectile exciting a spin-0 target
to a spin-$L_x$ state, the amplitudes are indexed by
$M_x = -L_x, \ldots, L_x$ and $J_P = L_x \pm 1/2$.

### 8.2 Transfer Reactions

For transfer reactions, the DWBA amplitude (Eq. 2.81 of the manual):

$$T(\vec{k}_i \to \vec{k}_{out}) = J \iint d^3r_i\,d^3r_{out}\,
[\chi^-(\vec{k}_{out}, \vec{r}_{out})]^*\,
\langle B,b | V_{eff} | A,a \rangle\,
\chi^+(\vec{k}_i, \vec{r}_i)$$

generates amplitudes labeled by the same quantum numbers. The analyzing
powers are computed identically through the Mueller matrix, with
the coefficients from `MUELCO` accounting for the spin structure of all
four particles.

### 8.3 Projectile vs. Target Analyzing Powers

Ptolemy can compute analyzing powers for the projectile beam (default,
array name `ANPOW`) and separately for target polarization. The
distinction enters through which particle's spin indices carry the
nonzero tensor rank in the Mueller matrix:

- **Projectile analyzing power**: $K_a \neq 0$, $K_b = K_A = K_B = 0$
- **Target analyzing power**: $K_A \neq 0$, $K_a = K_b = K_B = 0$

---

## 9. Summary of Key Relations

| Observable | Formula |
|-----------|---------|
| Unpolarized cross section | $\displaystyle\frac{d\sigma}{d\Omega} = \sum_k T_k\,\|F_k\|^2$ |
| Analyzing power $T_{KQ}$ | $\displaystyle T_{KQ} = (-1)^Q \frac{Z(K,-Q,0,0,0,0,0,0)}{Z(0,0,0,0,0,0,0,0)}$ |
| Mueller matrix element | $\displaystyle Z(\{K,Q\}) = \sum_{k,k'} C_{kk'}\,F_k\,F_{k'}^*$ |
| $A_y$ (spin-1/2) | $\displaystyle \frac{2\,\mathrm{Im}(F B^*)}{|F|^2+|B|^2}$ |
| $iT_{11}$ (spin-1/2) | $\displaystyle -\frac{A_y}{\sqrt{2}}$ |
| $iT_{11}$ (spin-1) | $\displaystyle -\frac{A_y}{i\sqrt{3}}$ |

---

## 10. Source Code Map

| Subroutine | Location | Role |
|-----------|----------|------|
| `ANAPOW` | `source.f:582–914` | Orchestrate analyzing power calculation and printing |
| `MUELCO` | `source.f:27278–27637` | Compute angle-independent Mueller matrix coefficients |
| `AMPCAL` | `source.f:222–580` | Compute scattering amplitudes $F_k(\theta)$ at one angle |
| `BETCAL` | `source.f:3796–4080` | Compute partial-wave coefficients $\beta_k(L)$ from S-matrix |
| `XSECTN` | `source.f:38258+` | Cross section calculation, calls AMPCAL and ANAPOW |
| `CUPSPN` | `source.f:13179–13292` | Channel-spin coupling coefficients (6-j symbols) |
| `WIG9J` | `source.f` | 9-j symbol evaluation |
| `CLEBSH` | `source.f` | Clebsch-Gordan coefficient evaluation |
