# Ptolemy Usage Guide

Ptolemy is a DWBA (Distorted Wave Born Approximation) nuclear reaction code.
Input is a plain-text file read from stdin. Keywords are case-insensitive
(first 8 characters significant). Lines starting with `$` are comments;
`$+` comments are echoed to output.

```bash
./ptolemy < input.in > output.out
```

**Important:** Always `rm -f fort.*` before running — leftover `fort.15`
causes an infinite loop.

---

## Table of Contents

1. [Input File Structure](#1-input-file-structure)
2. [Reaction Specification](#2-reaction-specification)
3. [Parameter Sets](#3-parameter-sets)
4. [Optical Potential Keywords](#4-optical-potential-keywords)
5. [Bound State Keywords](#5-bound-state-keywords)
6. [Integration & Grid Keywords](#6-integration--grid-keywords)
7. [Angular Momentum & Spin](#7-angular-momentum--spin)
8. [Control Switches](#8-control-switches)
9. [Output Control](#9-output-control)
10. [Linkules (Global OMP)](#10-linkules-global-optical-model-parametrizations)
11. [Fitting Keywords](#11-fitting-keywords)
12. [Example Input Files](#12-example-input-files)

---

## 1. Input File Structure

A typical input file has this structure:

```
REACTION: Target(Proj, Eject)Residual(Jpi Ex) ELAB=energy
PARAMETERSET name
keyword = value
keyword = value

PROJECTILE
  bound-state keywords
;
TARGET
  bound-state keywords
;
INCOMING
  optical potential keywords
;
OUTGOING
  optical potential keywords
;
angular keywords
;
end
```

- Semicolons `;` terminate a block and trigger computation.
- `RETURN` or `end` terminates the input.
- Multiple reactions can be chained in one file using `reset`.

---

## 2. Reaction Specification

### REACTION keyword

```
REACTION: 16O(D,P)17O(5/2+ 0.000) ELAB=20
REACTION: 40CA(P, P)40CA(2+ 4.000) ELAB=25
REACTION: 12C(8LI(2+), A)16N(3- 0.5) ELAB=28
```

Format: `Target(Projectile, Ejectile)Residual(Jpi Excitation) ELAB=energy`

- Particle names: standard chemical symbols (P, D, T, 3HE, A or 4HE, 16O, etc.)
- Jpi: spin and parity (e.g., `5/2+`, `0+`, `3-`, `1/2-`)
- Excitation energy in MeV
- ELAB in MeV (lab frame)

### CHANNEL keyword (alternative)

```
CHANNEL: 16O + 48CA
CHANNEL: P + 48CA
CHANNEL: 12C + 12C
CHANNEL: 12C(1/2 1.2) + 12C(1/2 1.2)   $ identical fermions
```

Used for elastic scattering. Spins and excitation can be specified in
parentheses.

### Supported reaction types

| Type | Example | PROBLM |
|---|---|---|
| Elastic scattering | `ELASTIC SCATTERING` after CHANNEL | 6 or 21 |
| Transfer (stripping) | `16O(D,P)17O` | 20 |
| Transfer (pickup) | `16O(P,D)15O` | 20 |
| 2-nucleon transfer | `18O(P,T)16O` | 20 |
| Heavy-ion transfer | `16O(12C,13C)15O` | 20 |
| Inelastic (DWBA) | `40CA(P,P)40CA(2+ 4.0)` | 22 |
| Coupled channels | `COUPLING (COUPLED, ...)` | 23/24 |

### Deformation parameter (inelastic)

```
BELX = 0.10          $ beta_L for multipolarity LX (usually L=2)
BETA = 0.10           $ same
BETACOUL = 0.08       $ Coulomb deformation (if different from nuclear)
```

### Ground state spin

```
JBIGA = 0+            $ target ground state spin-parity
JBIGB = 7/2-          $ residual ground state spin-parity
```

---

## 3. Parameter Sets

Parameter sets pre-configure integration grids and numerical parameters.

### Transfer DWBA

```
PARAMETERSET dpsb              $ standard DWBA
PARAMETERSET dpsb r0target     $ DWBA with r0 scaled to target mass
PARAMETERSET dpsa              $ alternate DWBA
```

| Name | Description |
|---|---|
| DPSA, DPSB | Standard DWBA parameter sets |
| DPDA, DPDB | Dense DWBA parameter sets |

### Elastic scattering

```
PARAMETERSET EL1     $ moderate accuracy
PARAMETERSET EL2     $ higher accuracy
PARAMETERSET EL3     $ highest accuracy
```

### Inelastic scattering

```
PARAMETERSET INELOCA1     $ moderate (ACCINE=1e-3, NPCOUL=6)
PARAMETERSET INELOCA2     $ good     (ACCINE=1e-4, NPCOUL=8)
PARAMETERSET INELOCA3     $ best     (ACCINE=1e-6, NPCOUL=12)
```

| Parameter | INELOCA1 | INELOCA2 | INELOCA3 |
|---|---|---|---|
| INELASACC | 1e-3 | 1e-4 | 1e-6 |
| DWCUTOFF | 1e-3 | 1e-4 | 1e-5 |
| STEPSPER | 15 | 20 | 25 |
| SUMPOINT | 6 | 8 | 12 |
| ASYMPTOPIA | 20 | 25 | 30 |
| LMAXADD | 20 | 30 | (large) |
| NPCOULOMB | 6 | 8 | 12 |
| ACCURACY | 1e-2 | 5e-3 | 1e-3 |

### Other sets

| Name | Description |
|---|---|
| PB100A-C | Lead at 100 MeV |
| CA60A-B | Calcium at 60 MeV |
| ALPHA1-3 | Alpha particle sets |

---

## 4. Optical Potential Keywords

Used inside `INCOMING` and `OUTGOING` blocks.

### Real central (Woods-Saxon)

```
V = 50.0          $ depth (MeV)
R0 = 1.17         $ radius parameter r0 (fm), R = r0 * A^(1/3)
A = 0.75          $ diffuseness (fm)
```

### Imaginary volume

```
VI = 2.0          $ depth (MeV)    (alias: W)
RI0 = 1.32        $ radius parameter (fm)
AI = 0.53         $ diffuseness (fm)
```

### Imaginary surface (derivative Woods-Saxon)

```
VSI = 8.0         $ surface depth (MeV)
RSI0 = 1.32       $ radius parameter (fm)
ASI = 0.53        $ diffuseness (fm)
```

### Real spin-orbit (Thomas form)

```
VSO = 6.0         $ spin-orbit depth (MeV)
RSO0 = 1.01       $ radius parameter (fm)
ASO = 0.75        $ diffuseness (fm)
TAU = 1.0         $ spin-orbit strength multiplier
```

### Imaginary spin-orbit

```
VSOI = 0.5        $ imaginary spin-orbit depth (MeV)
RSOI0 = 1.01      $ radius parameter (fm)
ASOI = 0.75       $ diffuseness (fm)
TAUI = 1.0        $ imaginary spin-orbit multiplier
```

### Coulomb

```
RC0 = 1.25        $ Coulomb radius parameter (fm), RC = rc0 * A^(1/3)
```

### Tensor potentials (for spin-1 projectiles)

```
VTR = 0.4         $ real TR tensor depth
RTR0 = 1.70       $ TR radius
ATR = 0.55        $ TR diffuseness
VTRI = 0.4        $ imaginary TR tensor
RTRI0 = 1.58      ATRI = 0.21
VTL = 1.0         $ TL tensor
VTPR = 5.0        $ approximate TP tensor
```

### Energy-dependent potentials

```
VE = 0.5          $ dV/dE coefficient
VIE = 0.3         $ dVI/dE coefficient
R0E = 0.01        $ dR0/dE
RI0E = 0.01       $ dRI0/dE
REALPOWER = 1.0   $ power of Woods-Saxon (default 1)
IMAGPOWER = 1.0   $ power of imaginary W.S.
```

---

## 5. Bound State Keywords

Used inside `PROJECTILE` and `TARGET` blocks.

```
NODES = 0         $ number of radial nodes (n_r = n - l - 1)
L = 2             $ orbital angular momentum
JP = 5/2          $ j of transferred particle
R0 = 1.25         $ bound-state W.S. radius parameter (fm)
A = 0.65          $ bound-state W.S. diffuseness (fm)
V = 60            $ W.S. depth (auto-adjusted to get correct binding)
VSO = 6           $ spin-orbit depth (MeV)
RSO0 = 1.10       $ spin-orbit radius parameter
ASO = 0.65        $ spin-orbit diffuseness
RC0 = 1.30        $ Coulomb radius parameter
```

### Wavefunction linkule

```
PROJECTILE
wavefunction av18       $ use AV18 deuteron wave function
r0=1 a=0.5 l=0 rc0=1.2
;
```

---

## 6. Integration & Grid Keywords

```
STEPSIZE = 0.1        $ radial step size (fm)
ASYMPTOPIA = 50       $ matching radius (fm)
SUMMAX = 30           $ upper integration limit
SUMMIN = 0            $ lower integration limit
SUMMID = 15           $ midpoint for grid mapping
LMIN = 0              $ minimum L
LMAX = 30             $ maximum L
LSTEP = 1             $ L step (1 or 2)
LMAXADD = 20          $ extra L values beyond automatic
MAXLEXTRAP = 5000     $ max L for extrapolation
LBACK = 0             $ backward L stepping
DWCUTOFF = 1e-4       $ distorted wave cutoff
ACCURACY = 1e-3       $ overall accuracy
INELASACC = 1e-4      $ inelastic accuracy
STEPSPER = 20         $ steps per oscillation
SUMPOINT = 8          $ Gauss points per Coulomb cycle
NPCOULOMB = 8         $ Coulomb integration points
MAXCOULOMB = 80       $ max Coulomb terms
NPSUM = 200           $ sum integration points
NPDIF = 20            $ difference integration points
NPPHI = 20            $ phi integration points
```

---

## 7. Angular Momentum & Spin

```
J = 3/2               $ total angular momentum
JA = 1                 $ spin of particle a
JB = 0                 $ spin of particle b
JBIGA = 0+             $ spin-parity of big A (target)
JBIGB = 7/2-           $ spin-parity of big B (residual)
JX = 2-                $ spin-parity of exchanged particle
JP = 5/2               $ j of transferred particle
SP = 1                 $ projectile spin (alias for JSP)
ST = 0                 $ target spin (alias for JST)
```

---

## 8. Control Switches

### Reaction type

```
ELASTIC SCATTERING     $ elastic scattering calculation
NOELASTIC              $ suppress elastic
DWBA                   $ DWBA transfer
INELASTIC              $ inelastic excitation
```

### Wave function type

```
REALWAVE               $ real wave functions only (no absorption)
COMPLEXWAVE            $ complex wave functions (default)
```

### Nuclear/Coulomb

```
NUCONLY                $ nuclear interaction only (no Coulomb correction)
USECOULOMB             $ include Coulomb (default)
USECORE                $ nuclear + Coulomb core corrections
```

### Asymptotic behavior

```
SKIPASYMPTOPIA         $ skip Whittaker matching
DOASYMPTOPIA           $ use Whittaker matching (default)
CHECKASYMPT            $ check asymptotic convergence
NOCHECKASYMP           $ don't check (default)
```

### Vertex convention

```
USEPROJECTILE          $ vertex at projectile
USETARGET              $ vertex at target
```

### Mass convention

```
TRUEMASS               $ use real masses
INTEGERMASS            $ use integer masses
```

### R0 convention

```
R0DEFAULT              $ R = r0 * (Aproj^1/3 + Atgt^1/3)
R0TARGET               $ R = r0 * Atgt^1/3
R0SUM                  $ R = r0 * (Ap + At)^1/3
```

### Form factor

```
FORMDERIV              $ derivative form factor (default for inelastic)
FORMASIS               $ direct form factor
```

### Miscellaneous

```
RESET                  $ clear all definitions
CLEAR                  $ clear channel definitions
UNDEFINE keyword       $ undefine a specific keyword
NEWPAGE                $ page break in output
HEADER text            $ title text for output
```

---

## 9. Output Control

```
PRINT = 2              $ print level (0=minimal, 999999=max debug)
ANGLESTEP = 1.0        $ angular step (degrees)
ANGLEMIN = 0.0         $ minimum angle
ANGLEMAX = 180.0       $ maximum angle
LABANGLE               $ output in lab frame
CMANGLES               $ output in CM frame (default)
WRITESTEP = 0.1        $ write wave function every 0.1 fm
PRINTWAVE              $ print bound-state wave function
```

### Named storage (NSCATALOG)

```
NSCATALOG              $ print allocator catalog
NSSTATUS               $ print allocator status
DUMPALL                $ dump all stored arrays
DUMP name              $ dump specific array
WRITENS4 name          $ write array as REAL*4 to fort.* file
WRITENS8 name          $ write array as REAL*8
```

---

## 10. Linkules (Global Optical Model Parametrizations)

Linkules are built-in optical potential parametrizations invoked by name.

```
PROJECTILE
wavefunction av18       $ use AV18 deuteron wavefunction
;
```

### Available linkules

| Name | Description |
|---|---|
| AV18 | Argonne V18 nucleon-nucleon potential |
| BKGPTELP | Becchetti-Greenlees proton/neutron global OMP |
| REID | Reid soft-core potential |
| GAUSSIAN | Gaussian potential form |
| SPLINE | Spline-interpolated potential |
| DEFORMED | Deformed Woods-Saxon |
| JDEPEN | L- and J-dependent potential |
| JDEPENWS | L- and J-dependent Woods-Saxon |
| PHIFFER | Phiffer potential |
| OHTA | Ohta parametrization |
| PARITWOO | Parity-dependent Woods-Saxon |
| LAGRANGE | Lagrange interpolation |
| SHAPE | Shape-dependent potential |
| TWOSHAPE | Two-body shape potential |
| RAWITSCH | Rawitsch parametrization |

### Linkule assignment keywords

```
REALPOTENTIAL linkule_name     $ real central potential
IMAGPOTENTIAL linkule_name     $ imaginary potential
REALSO linkule_name            $ real spin-orbit
IMAGSO linkule_name            $ imaginary spin-orbit
COULOMBPOT linkule_name        $ Coulomb potential
WAVEFUNCTION linkule_name      $ bound-state wavefunction
```

---

## 11. Fitting Keywords

```
FIT ( V  R0=RI0  AI )          $ fit these parameters
FIT ( R0=RI0=RC0  AI  V )      $ tied parameters (=)
FITACCURACY = 1e-4             $ fitting convergence
FITMULTIPLIER = 1.0            $ cross-section multiplier
FITRATIO = 1.0                 $ sigma_theory / sigma_data
MAXITER = 20                   $ maximum iterations
RENORM                         $ include normalization in fit

DATA ( ELAB=104 ANGLE SIGMATORU PERCENTER
  26.87  1.120  15.
  32.21  0.994   5.
  ...
)
```

---

## 12. Example Input Files

### Transfer: 16O(d,p)17O — neutron stripping

```
reset
REACTION: 16O(D,P)17O(5/2+ 0.000) ELAB=20
PARAMETERSET dpsb r0target
lstep=1 lmin=0 lmax=30 maxlextrap=0 asymptopia=50

PROJECTILE
wavefunction av18
r0=1 a=0.5 l=0 rc0=1.2
;
TARGET
JBIGA=0+
nodes=0 l=2 jp=5/2
r0=1.25 a=.65
vso=6 rso0=1.10 aso=.65
rc0=1.3
;
INCOMING
v=89.0 r0=1.15 a=0.75 vi=2.3 ri0=1.34 ai=0.60
vsi=10.2 rsi0=1.39 asi=0.69 vso=3.6 rso0=0.97 aso=1.01 rc0=1.30
;
OUTGOING
v=50.0 r0=1.17 a=0.75 vi=2.0 ri0=1.32 ai=0.53
vsi=8.0 rsi0=1.32 asi=0.53 vso=6.0 rso0=1.01 aso=0.75 rc0=1.25
;
anglemin=0 anglemax=180 anglestep=1
;
end
```

### Transfer: 16O(p,d)15O — neutron pickup

```
reset
REACTION: 16O(P,D)15O(1/2- 0.000) ELAB=25
PARAMETERSET dpsb r0target
lstep=1 lmin=0 lmax=30 asymptopia=50

PROJECTILE
wavefunction av18
r0=1 a=0.5 l=0 rc0=1.2
;
TARGET
JBIGA=0+
nodes=0 l=1 jp=1/2
r0=1.25 a=.65
vso=6 rso0=1.10 aso=.65
rc0=1.3
;
INCOMING
v=50.0 r0=1.17 a=0.75 vi=2.0 ri0=1.32 ai=0.53
vsi=8.0 rsi0=1.32 asi=0.53 vso=6.0 rso0=1.01 aso=0.75 rc0=1.25
;
OUTGOING
v=89.0 r0=1.15 a=0.75 vi=2.3 ri0=1.34 ai=0.60
vsi=10.2 rsi0=1.39 asi=0.69 vso=3.6 rso0=0.97 aso=1.01 rc0=1.30
;
anglemin=0 anglemax=180 anglestep=1
;
end
```

### Transfer: 18O(p,t)16O — two-neutron pickup

```
reset
REACTION: 18O(P,T)16O(0+ 0.000) ELAB=25
PARAMETERSET dpsb r0target
lstep=1 lmin=0 lmax=30 asymptopia=50

PROJECTILE
nodes=0 l=0 r0=1.25 a=0.65 rc0=1.3
;
TARGET
JBIGA=0+
nodes=0 l=0 jp=0
r0=1.25 a=.65
rc0=1.3
;
INCOMING
v=50.0 r0=1.17 a=0.75 vi=2.0 ri0=1.32 ai=0.53
vsi=8.0 rsi0=1.32 asi=0.53 vso=6.0 rso0=1.01 aso=0.75 rc0=1.25
;
OUTGOING
v=155.0 r0=1.20 a=0.72 vi=20.0 ri0=1.40 ai=0.88 rc0=1.30
;
anglemin=0 anglemax=180 anglestep=1
;
end
```

### Transfer: 16O(d,3He)15N — proton pickup

```
reset
REACTION: 16O(D,3HE)15N(1/2- 0.000) ELAB=40
PARAMETERSET dpsb r0target
lstep=1 lmin=0 lmax=30 asymptopia=50

PROJECTILE
nodes=0 l=0 r0=1.25 a=0.65 rc0=1.3
;
TARGET
JBIGA=0+
nodes=0 l=1 jp=1/2
r0=1.25 a=.65
vso=6 rso0=1.10 aso=.65
rc0=1.3
;
INCOMING
v=89.0 r0=1.15 a=0.75 vi=2.3 ri0=1.34 ai=0.60
vsi=10.2 rsi0=1.39 asi=0.69 vso=3.6 rso0=0.97 aso=1.01 rc0=1.30
;
OUTGOING
v=155.0 r0=1.20 a=0.72 vi=20.0 ri0=1.40 ai=0.88 rc0=1.40
;
anglemin=0 anglemax=180 anglestep=1
;
end
```

### Transfer: 12C(d,n)13N — proton stripping

```
reset
REACTION: 12C(D,N)13N(1/2- 0.000) ELAB=20
PARAMETERSET dpsb r0target
lstep=1 lmin=0 lmax=30 asymptopia=50

PROJECTILE
nodes=0 l=0 r0=1.25 a=0.65 rc0=1.3
;
TARGET
JBIGA=0+
nodes=0 l=1 jp=1/2
r0=1.25 a=.65
vso=6 rso0=1.10 aso=.65
rc0=1.3
;
INCOMING
v=89.0 r0=1.15 a=0.75 vi=2.3 ri0=1.34 ai=0.60
vsi=10.2 rsi0=1.39 asi=0.69 vso=3.6 rso0=0.97 aso=1.01 rc0=1.30
;
OUTGOING
v=50.0 r0=1.17 a=0.75 vi=2.0 ri0=1.32 ai=0.53
vsi=8.0 rsi0=1.32 asi=0.53 vso=6.0 rso0=1.01 aso=0.75
;
anglemin=0 anglemax=180 anglestep=1
;
end
```

### Transfer: 16O(12C,13C)15O — heavy-ion one-neutron

```
reset
REACTION: 16O(12C,13C)15O(1/2- 0.000) ELAB=120
PARAMETERSET dpsb r0target
lstep=1 lmin=0 lmax=50 asymptopia=30

PROJECTILE
nodes=0 l=0 r0=1.25 a=0.65 rc0=1.3
;
TARGET
JBIGA=0+
nodes=0 l=1 jp=1/2
r0=1.25 a=.65
vso=6 rso0=1.10 aso=.65
rc0=1.3
;
INCOMING
v=300.0 r0=1.24 a=0.68 vi=50.0 ri0=1.43 ai=0.87 rc0=1.30
;
OUTGOING
v=300.0 r0=1.24 a=0.68 vi=50.0 ri0=1.43 ai=0.87 rc0=1.30
;
anglemin=0 anglemax=60 anglestep=0.5
;
end
```

### Transfer: 12C(8Li,a)16N — heavy-ion four-nucleon

```
PARAMETERSET dpsb
REACTION: 12C(8LI(2+), A)16N(3- 0.5) JX=2- MXCX=25.93 ELAB=28
PROJECTILE
NODES=0 L=1
R0=1.20 A=.65 RC0=1.28
;
TARGET
NODES=0 L=2
R0=1.28 A=.76
;
INCOMING
V=152 r0=1.24 a=.685
w=23 ri0=1.432 ai=.87
RC0=1.3
;
OUTGOING
v=55 r0=1.24 a=.65
w=7 ri0=1.24 ai=.65
RC0=1.3
;
ANGLEMAX=180 anglestep=1
;
end
```

### Elastic scattering: 48Ca + 16O

```
CHANNEL: 16O + 48CA
ELAB=56
V=36.9 A=.4175 R0=1.35 RC0=1.2
VI=78 RI0=1.27 AI=.28
ELASTIC SCATTERING
;
RETURN
```

### Inelastic: 208Pb(p,p')208Pb* 2+

```
HEADER: 208Pb(p,p') 2+ at 4 MeV
JBIGA=0+
REACTION: 208PB(P, P)208PB(2+ 4.000) ELAB=20
BELX = 0.10
PARAMETERSET INELOCA3
INCOMING
V=53.0 R0=1.17 A=0.75 VI=2.0 RI0=1.32 AI=0.53
VSI=8.0 RSI0=1.32 ASI=0.53 VSO=6.0 RSO0=1.01 ASO=0.75 RC0=1.25
;
OUTGOING
;
;
RETURN
```

### Inelastic: 40Ca(a,a')40Ca* 2+

```
HEADER: 40Ca(a,a') 2+ at 4 MeV
JBIGA=0+
REACTION: 40CA(A, A)40CA(2+ 4.000) ELAB=80
BELX = 0.10
PARAMETERSET INELOCA3
INCOMING
V=120.0 R0=1.20 A=0.72 VI=25.0 RI0=1.40 AI=0.52 RC0=1.40
;
OUTGOING
;
;
RETURN
```

### Inelastic: 44Ca(16O,16O')44Ca* 2+

```
HEADER: 44Ca inelastic excitation
REACTION: 44CA(16O, 16O)44CA(2+ 1.156) ELAB=60
BELX = .0473
PARAMETERSET INELOCA1
INCOMING
V=36.9 VI=78
A=.4175 R0=1.35 RI0=1.27 AI=.28
RC0=1.2
;
OUTGOING
;
;
RETURN
```

### Elastic fit with data

```
CHANNEL: 16O + 208PB
HEADER: One-energy fit
PARAMETERSET EL1 ASYMPTOPIA=20
FIT ( R0=RI0=RC0  AI  V )
R0=1.3 A=.5 AI=.5 V=40 VI=15
DATA ( ELAB=104 ANGLE SIGMATORU PERCENTER
  26.87  1.120  15.
  32.21  0.994   5.
  37.54  1.005   5.
  42.84  1.006   5.
)
NSCATALOG
;
RETURN
```

### Multiple reactions in one file

```
$ First reaction
reset
REACTION: 16O(D,P)17O(5/2+ 0.000) ELAB=20
PARAMETERSET dpsb r0target
...
;

$ Second reaction
reset
REACTION: 16O(D,P)17O(1/2+ 0.871) ELAB=20
PARAMETERSET dpsb r0target
...
;

end
```

---

## Notes

- For inelastic with `OUTGOING ;` (empty), incoming potentials are copied
  to the outgoing channel automatically.
- The `NODES` keyword is the number of radial nodes (n_r), not the
  principal quantum number. For 1p: nodes=0, for 2s: nodes=1, for 1d:
  nodes=0, for 2p: nodes=1, for 1f: nodes=0.
- For (d,p) and (d,t), use `wavefunction av18` in the PROJECTILE block
  for the deuteron wave function. For other composite projectiles
  (3He, triton), specify `nodes`, `l`, `r0`, `a` explicitly.
- Spin-0 projectiles (alpha) do not need spin-orbit potentials.
- The `W` keyword is an alias for `VI` (imaginary volume depth).
