"""UED Simulation Framework — validation layer.

Manages the two INDEPENDENT physics backends:
  - AG_model        : AG/  (Kelisani 6D envelope ODE)
  - GPT_OCELOT_model: GPT模拟/ (OCELOT macroparticle tracking)

This layer only *drives* the existing backends through their public APIs and
reads shared parameters (shared/beamline_config.yaml).  It re-implements NO
physics kernel (Drift/Solenoid/RF are provided by the backends).
"""
