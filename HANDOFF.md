# UED Simulation — AI Session Handoff

Date: 2026-08-10 · 由上一 session 生成，供下一个 session 快速接管
工作目录: /Users/qin/Desktop/shuyan/Beam_dynamics_simu
解释器: /opt/anaconda3/bin/python3（系统 python3 无 yaml/ocelot，勿用）

## 1. 当前可信版本

branch: **main**
commit: **6ac4ae6**（同步脚本 A 级快照；v0.14 代码提交为 2b3d6e3）
latest tag: **v0.13-preSC-maintainability**（v0.14 未打 tag——SC 修复后需打 v0.14 tag 或并入 v0.14.1）

当前阶段:
**v0.14 SC audit 完成**（P0 已修，P1 未修）
下一阶段:
**v0.14.1 SC integration hardening**（scheduler 等价性 + AG 电荷语义 + SC 状态契约）

禁止跳到:
full SC validation / GUI

## 2. 项目目标

AG: **6D RMS envelope**（Kelisani 2023，Gaussian α 系数 / Luiten 2004 ellipsoid form factors）
OCELOT: **macroparticle tracking**（SecondTM + PIC SpaceCharge）
最终目标: 两个独立 backend + unified result（BeamResult）+ GUI（未来）

## 3. 当前架构

```
shared/beamline_config.yaml (唯一参数源, config SHA 见 §10)
   ↓ shared/params.py (parse/derived) + shared/beam_physics.py (BeamReference γ/β/p0)
   ↓ lattice.elements (唯一几何源)
   ├─ validation/backend.py::run_ag      (envelope ODE + thin-lens RF)
   ├─ validation/backend.py::run_ocelot  (SecondTM + RF kick + SC)
   └─ GPT模拟/ued_beamline_v2.py         (主路由, 同源 lattice)
   ↓ BeamResult (统一语义: σx σy σz σt σ_δ_p εnx εny + meta/provenance)
   ↓ validation/test_*.py → reports + baseline
```

## 4. 已解决的重大物理问题

| 问题 | 结论 |
|---|---|
| SOLENOID | 圆束降阶 Larmor 耦合必须 OFF（精确 4×4 ⇒ σ_xy≡0）；AG 适配层修正，σ_x 0.40% |
| RF | 薄透镜纵向 kick（H=eV·k·cosφ/(β²E₀)）；**RF 横向 kick 默认 OFF**（开关未独立验证）|
| R56 | 统一 δ_p=Δp/p0；OCELOT 原生 p_oc=ΔE/(c·p0)；适配 p_oc=β0·δ_p（只存在于 backend 边界）|
| Lattice | lattice.elements 全项目唯一几何源（主路由/验证路由/AG 同源，等价性测试证明）|
| RNG | config random.seed=42 配置化；generate_parray 用 seed、px/py/δ_p 用 seed+1；同 seed **位级可复现** |

## 5. 当前 SC 状态

v0.14 审计发现 **P0**:
```
OCELOT tracking_step() 只应用 transfer maps，PhysProc (SpaceCharge)
只在 track() 的 counter 机制中触发 → 生产代码 SC 从未执行
```
Fix（已实施，validation/backend.py + ued_beamline_v2.py + sc_audit_diagnostics.py）:
```
循环内复刻 counter：sc.counter-=1; 若≤0 → sc.z0=navi.z0; sc.apply(p, step*dz); counter=step
SC OFF 路径与 v0.13 位级一致（hash e041d6ae9fb7a0d2）
```
Evidence（500 fC 纯漂移 0.5 m）:
- σx: 725.3 → 2436.1 µm（**+236%**）· σz +428% · εnx +49%
- charge monotonicity: 0-1000 fC 单调（σx 725→3439 µm）PASS
- convergence: N(1e4-1e5)±0.3%、mesh(33-127)±0.6%、SC step≤5 PASS
- OCELOT SC 算法: NGP 沉积 + FFT Green 卷积 + 三线性插值 + 含 Ez + 束团系 Poisson（sc.py 全行号见 reports/SC_ocelot_source_audit.md）

## 6. 当前未解决问题

1. **manual PhysProc scheduler 尚未证明与 OCELOT native `track()` 严格等价**（P1）——下一步必须做 A/B 对照
2. **AG SC 电荷语义错误**（P1）：AG 用 Ne=beam.n_particles(50000)→Ne·e≈8 fC，而 config charge_fC=100 fC——弱 12.5 倍。AG 应使用 **Ne_phys = Q/e**
3. **SC runtime state contract 未实现**（requested/configured/available/attached/effective + HARD FAIL）——当前 ImportError 仍 silent fallback
4. **RF transverse physics 未独立验证**（开关默认 OFF，K_trans=-2.68 m⁻¹ 公式来自 AG，未对照文献/Panofsky-Wenzel）
5. OCELOT SC `prepare()` 固定 `np.random.seed(10)`（random_mesh=False 时影响小，需记录）
6. backend.py adapter+physics 混叠、主路由单体（结构性债务，SC/GUI 阶段再拆）

## 7. 禁止修改的东西

AG 核心公式（beam_dynamics_6d.py / external_forces.py / beamline_sim.py）· OCELOT 安装源码 · R56 公式与 p_oc 转换 · RF 纵向 kick 方程 · solenoid 已验证设置（圆束耦合 OFF）· baseline tags（v0.10/v0.11/v0.13）· 测试 tolerance · lattice 几何与 config 物理值

## 8. 关键文件

```
shared/beamline_config.yaml   唯一参数源（含 space_charge/random/physics_switches）
shared/params.py              parse/derived/lattice helpers
shared/beam_physics.py        BeamReference（γ/β/p0/velocity 单一派生）
shared/ocelot_coords.py       rparticles 语义访问（I_X..I_P, set/add_*）
shared/constants.py           物理常数单一来源
validation/backend.py         run_ag / run_ocelot / _ocelot_rf_kick / _provenance
validation/beam_result.py     BeamResult 容器
validation/config_check.py    Level-1 只读配置一致性检查
validation/sc_audit_diagnostics.py  SC 诊断（smoke/charge/convergence）
validation/test_{drift,solenoid,rf,full_beamline,gpt_route_equivalence,r56_convention,config_consistency}.py
GPT模拟/ued_beamline_v2.py    主路由（build_lattice_from_shared + counter SC）
AG/run_shared.py              AG 统一输出适配
scripts/sync_to_phone.py      iCloud 同步（--gpt-review 自动归档 A 级到 versions/<ts>/）
```

## 9. 当前数据语义（必须遵守）

```
δ_p = Δp/p0                    项目标准动量偏差
p_oc = ΔE/(c·p0)               OCELOT 原生第六坐标（只在 adapter 边界出现）
p_oc = β0 · δ_p                双向转换
σ_z = β0 · std(tau)            OCELOT 输出换算
σ_δ_p = std(p_oc)/β0           报告值恒为 δ_p
beam.n_particles               宏粒子数（OCELOT 用）
beam.charge_fC                 物理束团电荷（OCELOT generate_parray charge=Q）
AG: 必须用 Q/e 作为物理电子数   ← P1 待修
```

## 10. 当前验证基线

```
config_schema PASS · drift PASS · solenoid PASS · rf PASS
full_beamline PASS · gpt_route PASS · R56 PASS（全部 7 项）
SC OFF 必须保持位级不变（AG 位级 + OCELOT 同 seed 位级）
config SHA（v0.13 后新增 random 键）: 以当前 load_config()+config_sha() 为准
AG 关键参考值: σ_x=1984.191, σ_z=477.001（v0.10/0.12/0.13/0.14 逐位不变）
```

## 11. 下一任务 v0.14.1

A. **scheduler equivalence**：manual counter 复刻 vs OCELOT native `track()`（同一 seed/束/lattice，比较位级或浮点容差）——若不等价需修
B. **AG charge semantics**：AG 的 Ne 改为 Q/e（注意：这会改变 SC ON 数值——是修正，需先报告再实施，并更新 SC 相关测试）
C. **SC runtime state contract**：sc_requested/configured/available/attached/effective + HARD FAIL（禁 silent fallback）
完成前禁止正式 AG/OCELOT SC comparison。

## 12. STOP CONDITIONS

出现以下任意情况立即停止并报告（不连续打补丁）：
- 任何 no-SC regression 变化（位级）
- R56 表征变化
- AG no-SC 数组变化
- config lattice / 物理参数变化
- 单位/坐标语义不明确
- manual scheduler 与 OCELOT native 不一致
- 为使两模型吻合而准备调参

## 13. 审核包与同步

- iCloud: UED_Sync/gpt_review/（全量最新）+ versions/<时间戳>/（A 级精选，14 文件）
- 每次任务收尾自动跑 `python3 scripts/sync_to_phone.py --gpt-review`（AGENTS.md 约定）
- 代码自动 push 到 github.com/QIN-SMART/Beam-Dynamics（HTTP/1.1 + 清代理规避网络超时）
- 完整历史: validation/CHECKPOINTS.md（每次决策）+ validation/reports/v0.1[2-4]_*.md
