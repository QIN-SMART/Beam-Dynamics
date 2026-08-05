
# Phase 1.9 — Physical Solenoid Calibration

## Objective

Upgrade the current solenoid benchmark from an artificial focusing strength parameter k_sol to a physically calibrated solenoid model.

Do NOT add:
- space charge
- RF cavity
- longitudinal compression
- covariance model

Keep the current benchmark simple:
Drift + Solenoid + Drift.

The goal is to establish the physical relation:

B-field → electron momentum → solenoid focusing strength → beam waist


---

# Current status

Already validated:

1. Drift propagation
2. Emittance conservation
3. Solenoid focusing behavior

Current limitation:

The solenoid strength is manually set:

k_sol = 1.5

This is not directly connected to experimental parameters.


---

# Required modification 1: Use magnetic field as input

Replace:

```python
k_sol
````

with:

```python
B_sol
```

Input:

* magnetic field B [Tesla]
* solenoid length L [m]
* electron kinetic energy [keV]

---

# Required modification 2: Calculate solenoid focusing strength

Implement:

Electron relativistic parameters:

[
\gamma=1+\frac{E_k}{m_ec^2}
]

[
\beta=\sqrt{1-\frac{1}{\gamma^2}}
]

Momentum:

[
p=\gamma m_e\beta c
]

Solenoid Larmor focusing strength:

[
k_s=\frac{eB_z}{2p}
]

Check units carefully.

Document all conversions.

---

# Required modification 3: Compare three models

Generate comparison:

## Model A

OCELOT particle tracking

## Model B

4×4 transfer matrix envelope model

## Model C

Thin lens approximation

Thin lens:

[
\frac1f=k_s^2L
]

Compare:

* focal position
* beam waist
* minimum beam size

---

# Required output

Generate:

## Figure 1

Beam size evolution:

[
\sigma_x(z),\sigma_y(z)
]

Three curves:

* OCELOT
* transfer matrix
* thin lens

---

## Figure 2

Focal length comparison

Print:

```
B field =
Electron energy =

k_s =

Thin lens focal length =

OCELOT waist position =
Matrix waist position =
Thin lens prediction =
```

---

# Required validation

Test at least three magnetic fields:

Example:

B = 0.02 T

B = 0.05 T

B = 0.10 T

Expected behavior:

Stronger B:

[
k_s \propto B
]

Therefore:

[
f\propto \frac1{B^2}
]

Beam waist should move closer.

---

# Coding requirement

Do not refactor the whole project.

Only modify:

* solenoid parameter module
* benchmark script
* diagnostics

Keep compatibility with previous Phase 1.5 and Phase 1.8 results.

Create:

benchmark_solenoid_physical.py

Do not overwrite the old benchmark.




