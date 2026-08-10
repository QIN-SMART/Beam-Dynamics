"""
OCELOT ParticleArray coordinate access (v0.13) — semantic names for the
native 6-row layout, eliminating magic indices in business code.

Rows (OCELOT native, unchanged):
  0 x      [m]      1 px (=x')     2 y [m]      3 py (=y')
  4 tau    [m] (c·t) 5 p_oc = ΔE/(c·p0)

Conversion (project standard, R56 audit classification B):
  p_oc = β0 · δ_p          (into OCELOT)
  δ_p  = p_oc / β0         (out of OCELOT)
"""

I_X, I_PX, I_Y, I_PY, I_TAU, I_P = 0, 1, 2, 3, 4, 5


# ── generic accessors ─────────────────────────────────────────────────────
def get_x(p):        return p.rparticles[I_X]
def get_px(p):       return p.rparticles[I_PX]
def get_y(p):        return p.rparticles[I_Y]
def get_py(p):       return p.rparticles[I_PY]
def get_tau(p):      return p.rparticles[I_TAU]
def get_p_oc(p):     return p.rparticles[I_P]


# ── writers (used at the adapter boundary only) ────────────────────────────
def set_x(p, v):        p.rparticles[I_X, :] = v
def set_px(p, v):       p.rparticles[I_PX, :] = v
def set_y(p, v):        p.rparticles[I_Y, :] = v
def set_py(p, v):       p.rparticles[I_PY, :] = v
def set_tau(p, v):      p.rparticles[I_TAU, :] = v
def set_p_oc(p, v):     p.rparticles[I_P, :] = v
def add_p_oc(p, v):     p.rparticles[I_P, :] += v
def add_px(p, v):       p.rparticles[I_PX, :] += v
def add_py(p, v):       p.rparticles[I_PY, :] += v
