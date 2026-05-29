# Coupled-Channels Calculations

This document describes the coupled-channels (CC) method as implemented
in Ptolemy. CC goes beyond the DWBA by solving coupled differential
equations for multiple nuclear states simultaneously, allowing for
multi-step excitation and back-coupling effects.

---

## 1. When to Use Coupled Channels

DWBA treats the inelastic process as a first-order perturbation — the
incoming and outgoing distorted waves are generated independently. This
breaks down when:

- The coupling is strong ($\beta \gtrsim 0.15$)
- Multi-step paths contribute (e.g., $0^+ \to 2^+ \to 4^+$)
- Back-coupling from excited states modifies the elastic channel
- Coulomb excitation dominates (heavy ions, large $Z_P Z_T$)

CC solves these cases by coupling all channels together.

---

## 2. The Coupled Equations

### 2.1 Multi-Channel Schrödinger Equation

For $N$ coupled channels, the radial wave functions $\chi_c(r)$
satisfy:

$$\left[-\frac{\hbar^2}{2\mu}\left(\frac{d^2}{dr^2} - \frac{\ell_c(\ell_c+1)}{r^2}\right) + U_c(r) + V_C(r) - E_c\right]\chi_c(r) = -\sum_{c' \neq c} V_{cc'}(r)\,\chi_{c'}(r)$$

where:
- $c$ labels the channel (combination of nuclear state and partial wave)
- $U_c(r)$ is the diagonal optical potential
- $V_{cc'}(r)$ is the coupling interaction between channels $c$ and $c'$
- $E_c = E_{\text{c.m.}} - \varepsilon_c$ is the channel energy ($\varepsilon_c$ = excitation energy)

### 2.2 Coupling Interaction

For collective excitation, the coupling potential is:

$$V_{cc'}(r) = \beta_\lambda\,R_0\,\frac{dU}{dr}\,\langle J_f \| Y_\lambda \| J_i \rangle$$

where $\langle J_f \| Y_\lambda \| J_i \rangle$ is the reduced matrix
element coupling the nuclear states.

### 2.3 Boundary Conditions

- **Origin:** $\chi_c(0) = 0$ (regularity)
- **Asymptotic:** outgoing waves only in inelastic channels, incoming + outgoing in the elastic channel, defined by the S-matrix $S_{cc'}$

---

## 3. Solution Method

### 3.1 Iterative Coupled-Channels

Ptolemy solves the coupled equations iteratively:

1. Start with DWBA distorted waves (uncoupled solution)
2. Compute coupling terms $V_{cc'} \chi_{c'}$ as source terms
3. Solve the inhomogeneous equation for each channel
4. Iterate until the S-matrix converges to the specified accuracy

The convergence is controlled by `ACCURACY` and `MAXITER`.

### 3.2 S-Matrix Extraction

At each iteration, the asymptotic form of $\chi_c(r)$ is matched to
Coulomb functions $F_\ell$ and $G_\ell$ to extract the S-matrix:

$$\chi_c(r) \xrightarrow{r \to \infty} \frac{i}{2}\left[H_\ell^{(-)}(r)\,\delta_{c,\text{elastic}} - S_{cc'}\,H_\ell^{(+)}(r)\right]$$

### 3.3 Coulomb Born Correction

For heavy-ion systems with strong Coulomb coupling, the S-matrix
converges slowly with $\ell$. Ptolemy accelerates convergence using the
Coulomb Born approximation:

1. **FFACST:** Compute Born S-matrix elements from Coulomb FF2 integrals
   for all coupled channels
2. **SETBRN:** Subtract the Born S-matrix from the CC S-matrix for
   $\ell \leq \ell_{\max}$ (the CC values already include the Born
   contribution)
3. **L-extrapolation (LINTRP):** For $\ell > \ell_{\max}$, add back the
   Born S-matrix values directly (these partial waves are too high for
   the full CC solve, but the Born approximation is accurate there)

This effectively replaces the CC S-matrix at high $\ell$ with the Born
approximation while keeping the exact CC values at low $\ell$.

---

## 4. Input File Format

### 4.1 Basic Structure

```
REACTION: target(proj, proj)target  ELAB = energy
CHANNELS(NAME, ET, JT, EP, JP,  1, E1, J1pi,  E2, J2pi,  2, E3, J3pi, E4, J4pi)
COUPLING(COUPLED, channel-pairs)
TRANSITION(nucleus, MULTIPOLE, BETAN, BETAC, lambda, beta_n, beta_c)
PARAMETERSET  INELOCAn
```

### 4.2 CHANNELS

Defines the coupled states. Each channel specifies projectile and
target states:

```
CHANNELS(NAME, ET, JT, EP, JP,
    channel_number, excitation_target, Jpi_target, excitation_proj, Jpi_proj,
    ...)
```

- Channel numbers: 1 = elastic, 2 = first inelastic, etc.
- `ET, JT` = target excitation and spin-parity
- `EP, JP` = projectile excitation and spin-parity

**Example — two channels (g.s. + first 2+):**
```
CHANNELS(NAME, EP, JP, ET, JT,
    0, 0.0, 0+, 0.0, 0+,
    2, 0.0, 0+, 1.37, 2+)
```

**Example — target excitation only:**
```
CHANNELS(NAME, ET, JT, EP, JP,
    1, 0.0, 0+, 0.0, 0+,
    2, 2.615, 3-, 0.0, 0+)
```

### 4.3 COUPLING

Specifies which channels are coupled and the coupling model:

```
COUPLING(COUPLED, 1-2)
```

With explicit model:
```
COUPLING(COUPLED, ORDER, BETAN, BETAC, NUCLEUS, MULTIPOLE = lambda,
    MODEL = ROTATIONAL,
    0-2, order, beta_n, beta_c ;)
```

**Coupling models:**
- `ROTATIONAL` — rotational nucleus (rigid rotor $\beta$)
- `DEFORMED` — deformed optical potential

### 4.4 TRANSITION

Defines the coupling strength:

```
TRANSITION(TARGET, MULTIPOLE, BETAN, BETAC, lambda, beta_nuclear, beta_coulomb)
```

or for projectile excitation:

```
TRANSITION(PROJECTILE, MULTIPOLE, BELX, lambda, matrix_element)
```

- `TARGET` or `PROJECTILE` — which nucleus is excited
- `lambda` — multipolarity (1 = dipole, 2 = quadrupole, 3 = octupole)
- `BETAN, BETAC` — nuclear and Coulomb deformation parameters
- `BELX` — reduced matrix element (for projectile excitation)

---

## 5. Key Parameters

### 5.1 PARAMETERSET

| Setting | Use case |
|---------|----------|
| `INELOCA1` | Light ions, moderate coupling |
| `INELOCA2` | Medium-mass systems |
| `INELOCA3` | Heavy ions, strong Coulomb (highest precision) |

### 5.2 Convergence Control

| Keyword | Default | Description |
|---------|---------|-------------|
| `MAXITER` | 10 | Maximum CC iterations |
| `ACCURACY` | 0.01 | S-matrix convergence criterion |
| `ASYMPTOPIA` | 30 | Matching radius (fm); increase for heavy ions |

### 5.3 Partial Wave Range

| Keyword | Description |
|---------|-------------|
| `LMIN` | Minimum partial wave (use large values for distant collisions) |
| `LMAX` | Maximum partial wave for full CC solve |
| `LSTEP` | Partial wave step size |

For heavy-ion Coulomb excitation at energies below the barrier, most of
the cross section comes from large $\ell$ (distant collisions). Use
`LMIN` to skip the unimportant low partial waves.

---

## 6. Output

### 6.1 S-Matrix Table

With `PRINT = 2`, Ptolemy prints the CC S-matrix for each partial wave
and channel, including:

- Real and imaginary parts of S-matrix elements
- Unitarity check: $\sum_c |S_{c,\text{elastic}}|^2 = 1$
- Born S-matrix comparison (when Coulomb Born correction is active)

### 6.2 Cross Sections

The output contains angular distributions for each channel:
- **Elastic:** $d\sigma/d\sigma_{\text{Ruth}}$ (ratio to Rutherford)
- **Inelastic:** $d\sigma/d\Omega$ in mb/sr
- **TOTAL:** integrated inelastic cross section in mb

The first `TOTAL: 0.00` in the output is the elastic channel (Ptolemy
does not integrate the elastic cross section). The subsequent `TOTAL:`
values correspond to inelastic channels.

---

## 7. Examples

### 7.1 Heavy-Ion Coulomb Excitation

**208Pb(90Zr,90Zr')208Pb* 3- at 1350 MeV:**

```
REACTION:   208Pb(90Zr, 90Zr)208Pb  ELAB = 1350
CHANNELS(NAME, ET, JT, EP, JP, 1, 0.0, 0+, 0.0, 0+, 2, 2.615, 3-, 0.0, 0+)
COUPLING(COUPLED, 1-2)
TRANSITION(TARGET, MULTIPOLE, BETAN, BETAC, 3, .00321, .00321)
PARAMETERSET  INELOCA3   ASYMPTOPIA = 60
LABANGLES
ANGLESTEP = 0.2  ANGLEMAX = 30
LMIN = 400  LMAX = 1000  LSTEP = 20
INCOMING
V = 50  R0 = 1.133  A = .736
VI = 50  RI0 = 1.135  AI = .736
RC0 = 1.2
;
OUTGOING
;
```

Key points:
- `INELOCA3` for heavy-ion precision
- `LMIN = 400` skips low partial waves (sub-barrier, distant collisions)
- `ASYMPTOPIA = 60` extends the matching radius for the large system
- `LABANGLES` converts angles to lab frame

### 7.2 Light Heavy-Ion Rotational Coupling

**24Mg(12C,12C')24Mg* 2+ at Ecm = 20 MeV:**

```
REACTION   24Mg(12C, 12C)24Mg   ECM = 20
CHANNELS(NAME, EP, JP, ET, JT, 0, 0.0, 0+, 0.0, 0+, 2, 0.0, 0+, 1.37, 2+)
COUPLING(COUPLED, ORDER, BETAN, BETAC, NUCLEUS, MULTIPOLE = 2,
    MODEL = ROTATIONAL,
    0-2, 1, 0.52, 0.498 ;)
PARAMETERSET  INELOCA1   MAXITER = 15
LMAX = 55
SCATASYMP = 40  STEPSPER = 12
ANGLEMAX = 180  ANGLESTEP = 2
INCOMING
V = 36.0  R = 6.207  A = .52  RC0 = 1.24
W = 4.0
;
OUTGOING
R = 6.207  A = .52
;
```

Key points:
- `MODEL = ROTATIONAL` for the 24Mg rotational band
- `ECM` specifies center-of-mass energy (vs `ELAB` for lab energy)
- `ORDER` with explicit `BETAN, BETAC` in the COUPLING block
- Full angular range `ANGLEMAX = 180` for light-ion backscattering

---

## 8. CC vs DWBA Comparison

| Aspect | DWBA | CC |
|--------|------|-----|
| Coupling | First order | All orders |
| Multi-step | No | Yes |
| Back-coupling | No | Yes |
| Speed | Fast | Slower (iterative) |
| Accuracy | Good for weak coupling | Required for strong coupling |
| $\beta$ range | $\beta \lesssim 0.15$ | Any $\beta$ |

For a given system, DWBA and CC should agree when the coupling is weak.
The CC cross section will differ when multi-step processes or
back-coupling are significant.

---

## 9. Troubleshooting

### Convergence Issues

- **"CANNOT CONVERGE"**: Increase `MAXITER` or loosen `ACCURACY`
- **Oscillating S-matrix**: Try `PARAM20 = 0.7` (damping factor)
- **Wrong cross section scale**: Check that `BETAN` and `BETAC` have
  correct signs and magnitudes

### Performance

- Heavy-ion CC with $\ell_{\max} > 1000$ can be slow. Use `LSTEP > 1`
  to sample partial waves and interpolate.
- The Coulomb Born correction (FFACST) adds computation but greatly
  improves convergence for high-$Z$ systems.

### Common Mistakes

- Forgetting `LABANGLES` when comparing to lab-frame data
- Using `ELAB` when the potentials were fit in the c.m. frame (use `ECM`)
- Not setting `LMIN` high enough for sub-barrier heavy-ion reactions
  (wastes time on irrelevant partial waves)

---

## References

- T. Tamura, Rev. Mod. Phys. **37**, 679 (1965) — coupled-channels formalism
- G.R. Satchler, *Direct Nuclear Reactions*, Oxford University Press (1983) — Ch. 10
- M.H. Macfarlane and S.C. Pieper, ANL-76-11 Rev.1 (1978) — Ptolemy manual
- K. Alder and A. Winther, *Electromagnetic Excitation* (1975) — Coulomb excitation theory
