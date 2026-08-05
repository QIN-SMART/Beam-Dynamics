
# Phase 2A — Space Charge Drift Benchmark

## Background

The OCELOT UED particle tracking framework has completed:

- Phase 1: initial beam generation
- Phase 1.5: drift benchmark
- Phase 1.8: solenoid benchmark
- Phase 1.9: physical solenoid calibration

These benchmarks verified:

- particle initialization
- emittance calculation
- linear transport
- solenoid focusing physics

Now we start Space Charge validation.

IMPORTANT:
Do NOT add:
- solenoid
- RF cavity
- longitudinal compression
- chromatic effects

This phase only studies:

Electron bunch self-repulsion during free drift.

The purpose is to validate the Space Charge module independently.

---

# Physics objective

Compare:

## Case A: No space charge

Linear drift:

\[
x(z)=x_0+x'_0z
\]

Expected:

\[
\sigma_x(z)
=
\sqrt{
\sigma_{x0}^2+
\sigma_{x'}^2z^2
}
\]

Emittance:

\[
\epsilon_x=constant
\]


---

## Case B: With space charge

Electron bunch generates Coulomb repulsion.

Expected:

\[
\sigma_x(z)
\]

grows faster.


The transverse equation:

\[
\sigma_x''=
F_{sc,x}
+
F_{\epsilon}
\]


where:

\[
F_{sc,x}
\]

is the space charge defocusing term.


---

# Implementation requirements


## 1. Keep current beam parameters

Do not modify:

- energy
- charge
- initial spot size
- emittance
- particle number


Use the same parameters as Phase 1.5 and Phase 1.9.

The benchmark must be directly comparable.


---

## 2. Create a new benchmark file

Do NOT overwrite previous files.

Create:

```

benchmark_space_charge_drift.py

```


Keep:

```

benchmark_drift.py
benchmark_solenoid.py
benchmark_solenoid_physical.py

```

unchanged.


---

# 3. Implement two simulations


## Simulation A

Space charge OFF


Equivalent to Phase 1.5 drift.


## Simulation B

Space charge ON


Use OCELOT native SpaceCharge element/module.


Document:

- which OCELOT SC algorithm is used
- mesh parameters
- step size
- assumptions


---

# 4. Required diagnostics


Generate:


## Figure 1

Transverse beam size:

\[
\sigma_x(z),\sigma_y(z)
\]


Two curves:

- No SC
- With SC


Expected:

SC curve should expand faster.


---

## Figure 2

Longitudinal bunch length:

\[
\sigma_z(z)
\]


Compare:

No SC

vs

SC


Because space charge is a 3D effect.


---

## Figure 3

Emittance evolution:


\[
\epsilon_x(z)
\]


\[
\epsilon_y(z)
\]


Expected:


No SC:

\[
\epsilon=constant
\]


SC:

emittance growth may appear because of nonlinear Coulomb forces.


---

## Figure 4

Phase space:

x-x'


at:

z=0

and

final position


Compare:

No SC

vs

SC


---

# 5. Add charge scan

Run at least:


\[
Q=10 fC
\]


\[
Q=50 fC
\]


\[
Q=100 fC
\]


Study:


\[
\sigma_x(z,Q)
\]


Expected:

Higher charge:

stronger space charge expansion.


---

# 6. Add validation output


Print:


```

=============================
Space Charge Benchmark

Energy:
Charge:
Particles:

SC:
ON/OFF

Initial sigma_x:
Final sigma_x:

Initial emittance:
Final emittance:

# Emittance growth:

```


---

# 7. Important physics checks


Before completing, verify:


## Check 1

When:

```

charge = 0

```

SC result should approach no-SC result.


---

## Check 2

When:

```

SC=False

```

Result should reproduce Phase 1.5 drift benchmark.


---

## Check 3

Increasing charge should monotonically increase beam expansion.


---

# Coding philosophy

Do not optimize performance yet.

Prioritize:

1. physical correctness
2. reproducibility
3. clear comparison with analytical Gaussian model


The final goal is:

OCELOT particle tracking
        vs
Analytical Gaussian envelope model

under identical beam parameters.


