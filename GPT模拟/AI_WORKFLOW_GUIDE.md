# AI-Assisted Scientific Computing Workflow Guide

## Purpose

This document defines the required workflow for AI-assisted development of scientific simulation software.

The goal is to ensure that AI-generated code remains:

- physically interpretable
- scientifically traceable
- numerically verifiable
- maintainable

AI should not be treated as a black-box code generator.

AI acts as a computational physics assistant responsible for:

1. implementation
2. documentation
3. testing support

The human researcher remains responsible for:

1. physical assumptions
2. model selection
3. validation criteria
4. interpretation of results


---

# Core Principle

## Every physics implementation must have a corresponding Physics Note.

No new physics module should be created without documentation.

For every:

- force model
- beamline element
- numerical solver
- coordinate transformation
- diagnostic quantity

create a corresponding:

```

xxx_physics_note.md

```

Example:

```

solenoid.py

solenoid_physics_note.md

```

---

# Development Workflow

All AI coding tasks must follow the sequence:

```

Physical definition

```
    ↓
```

Mathematical model

```
    ↓
```

Numerical implementation

```
    ↓
```

Unit validation

```
    ↓
```

Integration test

```
    ↓
```

Documentation update

```


Do not directly jump from:

"implement feature"

to

"modify code".

---

# Before Modifying Code

AI must first analyze:

## 1. Existing implementation

Answer:

- Does this function already exist?
- Is there previous validated code?
- Can existing modules be reused?

Never duplicate existing physics.


---

## 2. Physical definition

Before coding, provide:

### Purpose

What physical phenomenon is modeled?


Example:

"Model transverse focusing caused by solenoid magnetic field."


---

### Coordinate system

Clearly define:

Example:

Transverse:

\[
(x,x')
\]


Longitudinal:

\[
(z,\delta)
\]


or particle coordinates:

\[
(x,p_x,y,p_y,z,\delta)
\]


---

### Variables and units

Every variable must have:

|Variable|Meaning|Unit|
|-|-|-|
|B|magnetic field|Tesla|
|E|electric field|V/m|
|sigma_x|beam RMS size|m|
|delta|relative momentum deviation|-|


No undefined variables.

---

# Physics Note Template

Every module must generate:

```

xxx_physics_note.md

```

with the following structure.


---

# [Module Name] Physics Note


## 1. Physical Purpose

Describe:

- what physical component is modeled
- what effect it produces


Example:

"The solenoid lens focuses electron beams through the azimuthal magnetic field and Larmor rotation."


---

## 2. Physical Model


Write governing equations.

Example:

Solenoid focusing:

\[
k_s=\frac{eB_z}{2p}
\]


Envelope equation:

\[
\sigma''+k_s^2\sigma=0
\]


Explain every term.


---

## 3. Assumptions


Explicitly list:

Example:

- paraxial approximation
- hard-edge magnetic field
- Gaussian beam distribution
- neglect fringe field


---

## 4. Coordinate Definition


Must specify:


Input coordinates:

\[
(x,x')
\]


Output coordinates:

\[
(x_f,x_f')
\]


Transformation:

\[
X_f=MX_i
\]


---

## 5. Numerical Implementation


Explain:

- algorithm
- integration method
- step size
- interpolation
- approximations


Example:

"Fourth-order Runge-Kutta integration is used for envelope propagation."


---

## 6. Input Parameters


Example:

```

Bz:
length:
energy:
charge:

```


---

## 7. Output Quantities


Example:

```

sigma_x(z)

sigma_y(z)

emittance

phase space

```


---

## 8. Validation


Every module requires validation.


Include:

### Analytical limit

Example:

Drift:

\[
\sigma_x(z)
=
\sqrt{
\sigma_{x0}^2+
\sigma_{x'}^2z^2
}
\]


### Numerical comparison

Compare:

simulation

vs

analytical solution


### Expected behavior

Example:

Increasing solenoid strength should reduce beam size.


---

## 9. Known Limitations


Must explicitly state:


Example:

Current model:

- ignores fringe field
- ignores nonlinear aberration
- assumes Gaussian distribution


---

# Modification Report Requirement

Every AI code modification must produce:

```

CHANGELOG_xxx.md

```


Format:


# Modification

## Date


## Reason

Why was this modification needed?


---

## Previous Problem

Example:

"RF compression did not produce longitudinal chirp."


---

## Physical Cause

Example:

"Tau coordinate unit mismatch."


---

## Modification

Describe code changes.


---

## Physics Impact

Explain:

Before:

\[
\sigma_z=...
\]


After:

\[
\sigma_z=...
\]


---

## Validation Result


Include:

- figures
- numerical comparison
- tests passed


---

# Debugging Rules


When simulation result is wrong:

Do NOT immediately modify parameters.


Follow:


## Step 1

Check units.


Examples:

- mm vs m
- eV vs Joule
- time vs length


---

## Step 2

Check coordinate definition.


Examples:

Is delta:

\[
\frac{\Delta p}{p_0}
\]

or:

\[
\frac{\Delta E}{E_0}
\]


---

## Step 3

Check limiting cases.


Examples:

Space charge OFF:

Should:

\[
\epsilon=constant
\]


No field:

Should:

\[
beam=drift
\]


---

## Step 4

Compare with analytical solution.


---

# Integration Rules


A validated module should not be rewritten during integration.


Example:


Already validated:

```

solenoid.py
rf.py
space_charge.py

```


Beamline assembly should only:

- import
- configure
- connect


Never:

- copy formulas
- duplicate force calculation


---

# Model Comparison Rules


When multiple simulation models exist:


Example:

AG model:

RMS envelope


GPT/OCELOT:

Particle tracking


They should share:


```

beam_parameters

beamline_definition

physical constants

```


Comparison outputs must use identical definitions:


```

sigma_x

sigma_y

sigma_z

emittance

energy spread

time resolution

chirp

```


---

# AI Behavior Rules


AI must:


1. Explain physics before coding.

2. Identify assumptions.

3. Reuse existing verified modules.

4. Generate documentation.

5. Provide validation.


AI must NOT:


1. Rewrite working modules unnecessarily.

2. Introduce unexplained formulas.

3. Hide unit conversions.

4. Optimize parameters to force agreement.

5. Claim physical correctness only because code runs.


---

# Final Research Principle


A simulation result is trustworthy only when:

\[
\boxed{
Physical\ Model
+
Numerical\ Implementation
+
Validation
+
Documentation
}
\]

are all consistent.


A beautiful plot without physical explanation is not a scientific result.

A simple model with clear assumptions and validation is more valuable than a complicated black-box simulation.

=
