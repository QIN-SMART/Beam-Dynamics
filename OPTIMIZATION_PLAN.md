
# UED Beam Dynamics Simulator Optimization Plan

## Project Background

This project aims to build a physics-based simulator for ultrafast electron diffraction (UED) beam dynamics.

Current implementation:

- RMS beam envelope model
- 3 transverse/longitudinal beam sizes:
    - σx
    - σy
    - σz
- slopes:
    - νx=dσx/dz
    - νy=dσy/dz
    - νz=dσz/dz
- relativistic acceleration:
    - γ evolution
- external fields:
    - solenoid lens
    - RF cavity
    - electrostatic acceleration
- space charge:
    - Gaussian model (Kelisani 2023)
    - Ellipsoid model (Luiten 2004)

Current state:

The code is a good RMS envelope framework, but it is not yet a complete 6D UED beam dynamics simulator.

The main missing physics:

1. Correct space charge implementation
2. Solenoid x-y coupling
3. Longitudinal phase space:
    - z-δ correlation
    - RF chirp
    - bunch compression
4. Full 6D covariance matrix tracking


---

# Phase 1: Fix Current Physics Model

Goal:

Make the existing RMS envelope model physically reliable.

Do NOT rewrite the whole architecture.

Only modify and validate existing modules.


---

# 1. Correct Space Charge Force Sign

## Problem

Current implementation:

```

Fs_x = -common * Mx/(σxσz)
Fs_y = -common * My/(σyσz)
Fs_z = -common * Mz/(σxσy)

```

This produces a focusing effect.

However, electron beam space charge is repulsive.

The envelope equation:

\[
\sigma_u''=F_u^e+F_u^{sc}+F_u^\epsilon
\]

requires:

\[
F^{sc}>0
\]


because space charge increases beam size.


## Required modification

Check and correct:

```

_ellipsoid_space_charge_forces()

```

The sign convention should be consistent with:

\[
F_x^{sc}\propto+\frac{K}{\sigma_x}
\]


Verify also:

```

_space_charge_forces()

```

Gaussian Kelisani model.


## Validation

Create a test:

No external fields:

- emittance = 0
- initial ν=0

Expected:

\[
\sigma_x(z),\sigma_y(z),\sigma_z(z)
\]

must monotonically increase.


---

# 2. Improve Space Charge Numerical Stability


Current:

alpha interpolation:

```

kx_range=(1e-2,1e2)
ky_range=(1e-2,1e2)

```


Need:

## Add clipping:

```

kx=np.clip(kx,1e-2,1e2)
ky=np.clip(ky,1e-2,1e2)

```


Avoid extrapolation instability.


## Add comparison test:

Compare:

Gaussian model

vs

Uniform ellipsoid model


For:

- round beam
- long bunch
- pancake beam


Check physical consistency.


---

# 3. Add Solenoid x-y Coupling


## Current limitation

Current model:

```

σx
σy

```

evolve independently.


But real solenoid:

\[
x-y
\]

coupled.

The solenoid introduces Larmor rotation:


\[
k_s=\frac{eB_z}{2p}
\]


and transport matrix:


\[
R_{solenoid}
\]


contains:

\[
x\leftrightarrow y
\]


## Required modification


Add a transverse covariance model:


Instead of only:

\[
\sigma_x,\sigma_y
\]


introduce:

\[
\langle xy\rangle
\]


or transverse 4D matrix:


\[
X=(x,x',y,y')
\]


with:


\[
\Sigma_4=
\langle XX^T\rangle
\]


At minimum:

add:

- xy coupling
- Larmor rotation angle


---

# Phase 2: Add UED Longitudinal Physics

Goal:

Enable realistic UED bunch compression simulation.


Current model cannot simulate:

- RF longitudinal lens
- velocity bunching
- temporal focus


because it lacks:


\[
\delta=\frac{\Delta p}{p}
\]


---

# 4. Add Longitudinal Energy Spread Variable


Current:


```

σz
νz

```


is insufficient.


Need:


\[
(z,\delta)
\]


Introduce:


\[
\sigma_\delta
\]


and correlation:


\[
C_{z\delta}
=
\langle z\delta\rangle
\]


State becomes:


\[
(\sigma_z,\sigma_\delta,C_{z\delta})
\]


---

# 5. Implement RF Longitudinal Lens (LL)


## Physical model


RF cavity gives:


\[
\Delta p_z=q\int E_z(t)dt
\]


Linearized:


\[
\delta
=
\delta_0+h z
\]


where:


\[
h=
\frac{1}{p_0}
\frac{dp}{dz}
\]


or:


\[
\Delta\delta=-hz
\]


Equivalent to thin lens:


\[
\begin{pmatrix}
z\\
\delta
\end{pmatrix}
\rightarrow
\begin{pmatrix}
1&0\\
-h&1
\end{pmatrix}
\begin{pmatrix}
z\\
\delta
\end{pmatrix}
\]


---

# Required implementation


Create:


```

longitudinal_lens.py

```


Containing:


```

RFLongitudinalLens()

```


Parameters:


- RF frequency
- voltage
- phase
- electron energy


Output:


chirp:

\[
h
\]


---

# 6. Implement Drift Compression


After LL:


particles evolve:


\[
z_f=z_i+R_{56}\delta
\]


where:


\[
R_{56}
=
\frac{L}{\gamma^2}
\]


approximately.


Compression condition:


\[
1+hR_{56}=0
\]


Need implement:


```

drift_longitudinal_transport()

```


---

# Phase 3: Upgrade to True 6D RMS Model


Goal:

Research-grade UED simulator.


---

# 7. Replace Envelope Model With 6D Covariance Tracking


Current:

Only:


\[
\sigma_x,\sigma_y,\sigma_z
\]


Upgrade:


phase space:


\[
X=
(x,x',y,y',z,\delta)
\]


Covariance:


\[
\Sigma=
\langle XX^T\rangle
\]


36 elements.


Evolution:


\[
\frac{d\Sigma}{dz}
=
A\Sigma+\Sigma A^T
\]


where:


\[
A
\]

is beam transport matrix.


---

# 8. Implement Beam Transport Matrices


Need modules:


```

transport_matrix.py

```


Elements:


## Drift

\[
R_{drift}
\]


## Solenoid

Include:

- focusing
- rotation
- coupling


## RF cavity


Include:


- transverse RF focusing
- longitudinal chirp


## Quadrupole


Include:


- x focusing
- y defocusing


---

# 9. Space Charge in 6D Model


Use self-consistent approach:


Option 1:

RMS KV/envelope approximation


Option 2:

Particle-in-cell / macroparticle tracking


For UED:

Start with RMS.


---

# Phase 4: Validation Against Literature


Before claiming physical accuracy:


Reproduce:


## Case 1

100 keV UED beam


Parameters:


- energy
- bunch charge
- initial size
- emittance


Compare:


- σx(z)
- σz(z)


---

## Case 2

RF bunch compression


Compare:


Before compression:

\[
\sigma_z
\]


After LL:


\[
\sigma_z
\]


---

## Case 3

Temporal resolution


Calculate:


\[
\Delta t
=
\sqrt{
\Delta t_L^2+
\Delta t_e^2
}
\]


where:


\[
\Delta t_e
=
\frac{\sigma_z}{v}
\]


---

# Development Rules


## Important

Do not blindly rewrite the project.

Follow incremental development:


Priority order:


1. Fix space charge sign
2. Validate solenoid focusing
3. Add longitudinal chirp variable δ
4. Implement LL thin lens
5. Add drift compression
6. Upgrade covariance matrix


Every modification must include:


- physical equation
- code location
- unit checking
- validation example


---

# Required Documentation

For every new module:

Provide:

1. Physical principle
2. Mathematical derivation
3. Numerical implementation
4. Test case


The final goal:

A physically interpretable UED beam dynamics simulator,
not only a numerical fitting program.



