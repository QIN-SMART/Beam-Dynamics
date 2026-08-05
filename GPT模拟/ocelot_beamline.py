#!/usr/bin/env python3
"""
OCELOT 束线模拟 — 与 beam_dynamics.py 对应束线

运行：
  python3 ocelot_beamline.py
  python3 ocelot_beamline.py --no-sc     # 无空间电荷
  python3 ocelot_beamline.py --ne 500000  # 修改电子数

首次运行需安装 OCELOT：
  pip install git+https://github.com/ocelot-collab/ocelot.git
"""

import sys, numpy as np

print("加载 OCELOT（约 5—10 秒）...")
import ocelot
from ocelot.cpbd.elements import Drift, Solenoid, Cavity, Marker
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.beam import generate_parray
from ocelot.cpbd.navi import Navigator
from ocelot.cpbd.track import tracking_step

# ═══════════════════════════════════════════════════════════════════════
# 参数
# ═══════════════════════════════════════════════════════════════════════

use_sc = '--no-sc' not in sys.argv
Ne = 100_000
spot_rms   = 85e-6
sig_z0     = 300e-6
E_keV      = 100.0
epsilon_n  = 0.08e-6             # normalized emittance [m·rad]
sigma_delta = 1.0e-4              # relative energy spread ΔE/E

for i, a in enumerate(sys.argv):
    if a == '--ne' and i + 1 < len(sys.argv):
        Ne = int(float(sys.argv[i + 1]))
    elif a == '--epsn' and i + 1 < len(sys.argv):
        epsilon_n = float(sys.argv[i + 1])
    elif a == '--sigd' and i + 1 < len(sys.argv):
        sigma_delta = float(sys.argv[i + 1])

# 元件 z 坐标 (m)
z_tl1  = 0.100
z_ap1  = 0.270
z_ll1  = 0.400
z_tl2  = 0.495
z_ap2  = 0.615
z_samp = 0.777
z_det  = 1.015

# 透镜参数（需手动标定使 σ_x 匹配 PhaseSpace 模型）
k_tl1 = 1.5          # 对应 I=0.90 A
k_tl2 = 3.5          # 对应 I=2.70 A

# 构建束线（跳过 Aperture + Cavity V=0 以避免 RF 横向效应）
lat = MagneticLattice([
    Drift(l=z_tl1, eid='D1'),
    Solenoid(l=0.06, k=k_tl1, eid='TL1'),
    Drift(l=z_ap1 - z_tl1, eid='D2'),
    Drift(l=z_ll1 - z_ap1, eid='D3'),
    Cavity(l=0.022, freq=2.856e9, v=0, phi=0, eid='LL1'),
    Drift(l=z_tl2 - z_ll1, eid='D4'),
    Solenoid(l=0.06, k=k_tl2, eid='TL2'),
    Drift(l=z_ap2 - z_tl2, eid='D5'),
    Drift(l=z_samp - z_ap2, eid='D6'),
    Drift(l=z_det - z_samp, eid='D7'),
])
lat.update_transfer_maps()

# ═══════════════════════════════════════════════════════════════════════
# 束团
# ═══════════════════════════════════════════════════════════════════════

# relativistic parameters (needed for tau=c·t initialisation)
E_rest = 511.0              # keV
gamma = 1.0 + E_keV / E_rest
beta = np.sqrt(1.0 - 1.0 / gamma**2)
beta_gamma = beta * gamma

p = generate_parray(
    sigma_x=spot_rms, sigma_y=spot_rms,
    sigma_tau=sig_z0 / beta,          # OCELOT tau = c·t [m]; σ_tau = σ_z/β
    energy=(E_keV + 511.0) * 1e-6,    # TOTAL energy in GeV (E_kin+mc²)
    charge=Ne,
)

# ＝＝ Phase 1: improved initial beam distribution ＝＝

# 1. transverse emittance → angular spread x', y'
epsilon_geom = epsilon_n / beta_gamma         # geometric emittance [m·rad]
sigma_xp = epsilon_geom / spot_rms             # σ_x' [rad]
sigma_yp = epsilon_geom / spot_rms             # σ_y' [rad] (round beam)

np.random.seed(42)
n_particles = p.rparticles.shape[1]
p.rparticles[1, :] = np.random.normal(0.0, sigma_xp, n_particles)   # px (= x')
p.rparticles[3, :] = np.random.normal(0.0, sigma_yp, n_particles)   # py (= y')

# 2. energy spread δ = ΔE/E
p.rparticles[5, :] = np.random.normal(0.0, sigma_delta, n_particles)

# verify initial distribution
x  = p.x()
xp = p.px()
y  = p.y()
yp = p.py()
dd = p.p()
eps_x_initial = np.sqrt(np.mean(x**2) * np.mean(xp**2) - np.mean(x * xp)**2) * 1e6
eps_y_initial = np.sqrt(np.mean(y**2) * np.mean(yp**2) - np.mean(y * yp)**2) * 1e6
print(f"\n  初始分布:")
print(f"    ε_n  = {epsilon_n*1e6:.4f} mm·mrad  (input)")
print(f"    ε_geo = {epsilon_geom*1e6:.4f} mm·mrad  (ε_n/βγ, βγ={beta_gamma:.4f})")
print(f"    ε_x = {eps_x_initial:.4f} mm·mrad  (from particles)")
print(f"    ε_y = {eps_y_initial:.4f} mm·mrad  (from particles)")
print(f"    σ_δ  = {np.std(dd)*1e3:.2f} e-3  (from particles)")
print(f"    σ_x' = {sigma_xp*1e3:.3f} mrad,  σ_y' = {sigma_yp*1e3:.3f} mrad")

# ═══════════════════════════════════════════════════════════════════════
# 追踪
# ═══════════════════════════════════════════════════════════════════════

navi = Navigator(lat)
dz = 0.005
total = sum(e.l if hasattr(e, 'l') else 0 for e in lat.sequence)

# SC 模块
sc = None
if use_sc:
    try:
        from ocelot.cpbd.sc import SpaceCharge
        sc = SpaceCharge(step=1)
    except Exception as e:
        print(f"SC 加载失败: {e}")

# 记录面 z 坐标
probes = {
    'AccelExit': z_tl1,          # 加速段出口 ≈ TL1 位置（此处以来束）
    'TL1':      z_tl1,
    'AP1':      z_ap1,
    'LL1':      z_ll1,
    'TL2':      z_tl2,
    'AP2':      z_ap2,
    'Sample':   z_samp,
    'Detector': z_det,
}
results = {}

for step_i in range(int(total / dz)):
    z_before = navi.z0
    tracking_step(lat, p, dz, navi)
    if sc is not None:
        sc.apply(lat, p, dz)
    z_after = navi.z0

    for name, z_p in list(probes.items()):
        if name not in results and z_before <= z_p < z_after:
            xa  = p.x()
            ya  = p.y()
            ta  = p.tau()
            xpa = p.px()
            ypa = p.py()
            da  = p.p()
            sx = np.std(xa) * 1e6
            sy = np.std(ya) * 1e6
            st = np.std(ta) * 1e12
            ex = np.sqrt(np.mean(xa**2)*np.mean(xpa**2) - np.mean(xa*xpa)**2) * 1e6
            ey = np.sqrt(np.mean(ya**2)*np.mean(ypa**2) - np.mean(ya*ypa)**2) * 1e6
            sd = np.std(da) * 1e3
            results[name] = (z_p * 1e3, sx, sy, st, ex, ey, sd)

sx_f = np.std(p.x()) * 1e6
sy_f = np.std(p.y()) * 1e6
st_f = np.std(p.tau()) * 1e12
xa_f = p.x(); xpa_f = p.px()
ya_f = p.y(); ypa_f = p.py()
da_f = p.p()
ex_f = np.sqrt(np.mean(xa_f**2)*np.mean(xpa_f**2) - np.mean(xa_f*xpa_f)**2) * 1e6
ey_f = np.sqrt(np.mean(ya_f**2)*np.mean(ypa_f**2) - np.mean(ya_f*ypa_f)**2) * 1e6
sd_f = np.std(da_f) * 1e3

# ═══════════════════════════════════════════════════════════════════════
# 输出
# ═══════════════════════════════════════════════════════════════════════

print(f"\n{'='*80}")
print(f"  OCELOT 束线模拟  Ne={Ne}  SC={'ON' if use_sc else 'OFF'}")
print(f"{'='*80}")
print(f"  {'位置':>10s}  {'z(mm)':>7s}  {'σ_x(μm)':>9s}  {'σ_y(μm)':>9s}  {'σ_t(fs)':>9s}  {'ε_x(mm·mrad)':>13s}  {'ε_y(mm·mrad)':>13s}  {'σ_δ(e-3)':>9s}")
print(f"  {'-'*78}")

for name in ['AccelExit', 'TL1', 'AP1', 'LL1', 'TL2', 'AP2', 'Sample', 'Detector']:
    if name in results:
        z_mm, sx, sy, st, ex, ey, sd = results[name]
        print(f"  {name:>10s}  {z_mm:7.0f}  {sx:9.1f}  {sy:9.1f}  {st:9.0f}  {ex:13.4f}  {ey:13.4f}  {sd:9.1f}")

print(f"\n  最终: σ_x={sx_f:.1f} μm,  σ_y={sy_f:.1f} μm,  σ_t={st_f:.1f} fs")
print(f"        ε_x={ex_f:.4f} mm·mrad,  ε_y={ey_f:.4f} mm·mrad,  σ_δ={sd_f:.1f} e-3")
print(f"  传输率: {p.n} 粒子")

print(f"\n  对比 Python PhaseSpace (Ne={Ne}):")
print(f"    Sample σ_x ≈ {'94' if Ne==100000 else '~100'} μm")
print(f"    差异来源：")
print(f"      • Solenoid.k 值需用磁场标定（当前为估计值）")
print(f"      • RF 腔设为 V=0（避免横向 defocusing，匹配 PhaseSpace 薄透镜近似）")
print(f"      • 光阑、空间电荷未启用（启用请去掉 --no-sc）")
print(f"      • OCELOT 的 RF 横向效应是 PhaseSpace 未建模的额外物理")
print(f"\n  命令行参数:")
print(f"    --ne N     电子数 (默认 {100_000})")
print(f"    --epsn X   归一化发射度 [m·rad] (默认 {0.08e-6})")
print(f"    --sigd X   相对能散 ΔE/E (默认 {1.0e-4})")
print(f"    --no-sc    关闭空间电荷")
