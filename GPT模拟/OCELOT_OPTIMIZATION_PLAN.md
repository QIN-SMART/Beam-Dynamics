
# OCELOT UED Particle Tracking Model Optimization Plan

## Project Purpose

This project implements an independent particle tracking model for UED beam dynamics using OCELOT.

The purpose is NOT to replace the Analytical Gaussian (AG) envelope model.

The purpose is:

1. Build a General Particle Tracker (GPT-like) model.
2. Use it as an independent benchmark.
3. Compare particle tracking results with analytical Gaussian envelope results.

The final goal:

Validate the physical correctness of the analytical model.

---

# Model Architecture

Two independent simulation approaches:

## Model A: Analytical Gaussian Envelope

State:

\[
\sigma_x,\sigma_y,\sigma_z
\]

Physics:

- envelope equation
- space charge
- emittance
- external focusing


Purpose:

Fast optimization.


---

## Model B: OCELOT Particle Tracking

Particle coordinates:

\[
X_i=
(x,y,z,p_x,p_y,p_z)
\]


Each electron/macroparticle is tracked individually.

Physics:

\[
\frac{d\vec p}{dt}
=
q(\vec E+\vec v\times \vec B)
\]


Purpose:

High-fidelity validation.

---

# Current Code Status

Current file:

```

ocelot_beamline.py

```


Implemented:

- OCELOT lattice
- Drift
- Solenoid
- Cavity placeholder
- Space charge interface
- Macroparticle tracking


Current limitations:

1. Initial beam distribution incomplete.
2. RF longitudinal compression not implemented.
3. Solenoid strength not physically calibrated.
4. Missing beam diagnostics.
5. Missing benchmark tests.

---

# Development Rule

IMPORTANT:

Do NOT blindly increase model complexity.

Do NOT add covariance matrix model.

Do NOT merge with Analytical Gaussian code.

This project should remain an independent particle tracking benchmark.

Every modification must include:

1. Physical equation.
2. Unit consistency.
3. Validation test.

---

# Phase 1: Improve Initial Beam Distribution

## Goal

Generate a realistic UED electron bunch.

Current:

Only:

\[
\sigma_x,\sigma_y,\sigma_t
\]


Need to include:


## 1. Transverse emittance


Electron beam has angular spread:

\[
x'
\]


Relationship:

\[
\epsilon_x
=
\sigma_x\sigma_{x'}
\]


Include:

- normalized emittance
- transverse momentum distribution


Example:

Input:

```

epsilon_n

````

Calculate:

\[
\epsilon
=
\frac{\epsilon_n}{\beta\gamma}
\]


Then:


\[
\sigma_{x'}
=
\frac{\epsilon}{\sigma_x}
\]


Generate:


\[
x'
\sim N(0,\sigma_{x'})
\]


---

## 2. Energy spread


Add:

\[
\delta=
\frac{\Delta E}{E}
\]


Initial distribution:


\[
\delta
\sim N(0,\sigma_\delta)
\]


This is necessary for:

- chromatic effects
- RF compression


---

# Phase 2: Add Beam Diagnostics

The tracker must output more than beam size.

Add:

---

## 1. RMS beam size


Already exists:


\[
\sigma_x,\sigma_y,\sigma_z
\]


Keep.

---

## 2. Transverse emittance


Calculate:


\[
\epsilon_x
=
\sqrt{
<x^2><x'^2>
-
<xx'>^2
}
\]


Same for y.


---

## 3. Energy spread


Output:


\[
\sigma_\delta
\]


---

## 4. Longitudinal phase space


Generate plots:


\[
z-\delta
\]


before and after RF.


This is essential for UED bunch compression.

---

# Phase 3: Correct Solenoid Implementation

## Problem

Current:

```python
Solenoid(k=k_tl)
````

uses estimated k.

Need physical calibration.

---

## Required input

Use experimental magnetic field:

[
B_z
]

Calculate:

[
k_s=
\frac{eB_z}{2p}
]

where:

[
p=\gamma mv
]

Then convert to OCELOT parameter.

---

## Validation

Test:

Only one solenoid.

Expected:

* beam waist
* focusing
* symmetric x/y evolution

Compare:

OCELOT

vs

Analytical Gaussian model.

---

# Phase 4: Implement Real RF Longitudinal Lens

## Current problem

Current:

```python
Cavity(v=0)
```

means no energy modulation.

Therefore:

No bunch compression.

---

## Required physics

RF field:

[
E_z=E_0\cos(\omega t+\phi)
]

Energy gain:

[
\Delta E
========

eV\cos(\phi+kz)
]

Linear approximation:

[
\delta(z)
=========

hz
]

where:

[
h=
\frac{1}{E_0}
\frac{dE}{dz}
]

---

## Implementation

Enable:

```
Cavity(v=V_RF)
```

Parameters:

* frequency
* voltage
* phase

Study:

[
\sigma_z
]

before and after compression.

---

# Phase 5: Space Charge Validation

## Current implementation

Uses:

```python
SpaceCharge()
```

Need verification.

---

## Benchmark

Run:

### Case A

No external fields.

Only space charge.

Expected:

[
\sigma_x(z)
]

monotonically increases.

---

Compare:

OCELOT

vs

Analytical Gaussian SC model

---

# Phase 6: Benchmark Tests

Before full UED simulation:

## Test 1: Drift

No SC.

No lens.

Analytical:

[
\sigma_x(z)
===========

\sqrt{
\sigma_0^2+
\sigma_{x'}^2z^2
}
]

OCELOT should reproduce.

---

## Test 2: Solenoid

Compare:

beam waist position.

---

## Test 3: Space charge

Compare:

beam expansion rate.

---

## Test 4: RF compression

Compare:

minimum bunch length.

---

# Final UED Beamline

After validation:

```
Photocathode

↓

initial bunch distribution

↓

Solenoid TL

↓

RF longitudinal lens

↓

Drift compression

↓

Solenoid focusing

↓

Sample

↓

Detector
```

Outputs:

* beam size
* emittance
* bunch length
* time resolution

---

# Scientific Goal

The final project should establish:

Analytical Gaussian Model

```
    |

    | benchmark

    ↓
```

OCELOT Particle Tracking Model

Agreement between two models validates the physical assumptions.



