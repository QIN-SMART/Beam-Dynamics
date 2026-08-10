# SC Integration Audit (v0.14)

Date: 2026-08-10 · 先审计，最小修复需经本报告批准后实施

## 1. OCELOT SC 源码证据
见 `SC_ocelot_source_audit.md`（sc.py 全部关键行号）。

## 2. AG SC 公式
见 `SC_AG_model_audit.md`（Gaussian Kelisani + Ellipsoid Luiten，公式与代码位置）。

## 3. 两模型假设差异
| 维度 | AG (Gaussian/Ellipsoid) | OCELOT (PIC) |
|---|---|---|
| 表示 | 6D RMS 包络 / 均匀椭球 | 宏粒子 + 场求解 |
| SC 力 | 解析 α 系数 / form factors，只依赖 σ 矩 | 逐粒子：束团静止系 Poisson（Green 卷积 FFT）+ 三线性插值 |
| 分布假设 | 高斯/椭球 | 任意（无假设）|
| emittance | 输入常量，**天然不增长** | 投影发射度可因非线性 SC 增长 |
| 非线性 | 无 | 有（近点非线性、halo）|
| 数值噪声 | 无 | 有（NGP 沉积、网格）|

## 4. shared config → SC 数据流
```
config space_charge.{enabled, mesh, step}
   ├─ OCELOT backend:  SpaceCharge(step=1) ← ① mesh 未传(默认63³，恰好=config)
   ├─ GPT 主路由:      SpaceCharge(step=config.step) ← ② mesh 未传
   └─ AG backend:      Ne → fb (Ne·e 总电荷)  ← ③ 与 config charge_fC 脱钩
```
**③ AG SC 强度 bug（P1）**：AG 用 `Ne·e` = 50000·1.6e-19 = **8 fC 等效总电荷**，而 config 是 **100 fC**——差 12.5 倍（AG SC 偏弱）。修复会改变 SC ON 数值 → 按任务十二先报告暂停，本阶段不实施。

## 5. SC state 控制（设计）
现状：无状态机。设计（本阶段仅记录）：
```
sc_requested = (sc_enabled 入参 or config.enabled) and config.enabled
sc_configured = SpaceCharge(step=config.step, nmesh_xyz=config.mesh)
sc_available = import 成功（当前 OK）
sc_attached = navi.add_physics_proc 成功
sc_effective = apply 至少被调用一次（需 get_next_step 机制）
规则：sc_requested 且 (¬available or ¬attached) → HARD FAIL（禁 silent fallback）
```

## 6. mesh 真实状态
从运行对象读取（诊断脚本已验证）：`nmesh_xyz=[63,63,63], step=1, random_seed=10, random_mesh=False` —— 对象属性正常，但 **apply 从未执行**（见 P0）。

## 7. 当前发现的 bug
| 级别 | 描述 | 状态 |
|---|---|---|
| **P0** | `tracking_step()` 不触发 physics procs：SpaceCharge.apply 从未被调用（track.py:476-479 vs 406-425）。生产代码 backend.py:392 / ued_beamline_v2.py:202 均用 tracking_step → **SC 完全无效** | 本报告后实施最小修复 |
| P1 | mesh 未从 config 传入（值恰好=默认，63³；若改 config 将失效）| 随 P0 一起修 |
| P1 | AG SC 电荷语义：Ne·e vs config charge_fC（8 vs 100 fC，弱 12.5 倍）| 报告暂停 |
| P1 | SC state 无双源检查：ued_beamline_v2 step>=4 与 config.enabled 双来源 | 报告 |
| P2 | silent fallback（ImportError→pass）| 报告 |

## 8. 建议修复（最小，数值不变或仅 SC 路径）
1. **backend.py / ued_beamline_v2.py**：tracking 循环改用
   `for t_maps, dz_step, proc_list, phys_steps in navi.get_next_step():`
   （与 OCELOT track() 相同机制），SC OFF 时 proc_list 为空 → no-SC 行为不变。
2. SpaceCharge 构造传 `nmesh_xyz=config.mesh, step=config.step`（当前值相同，改后数值不变）。
3. AG 侧与 SC 状态机：留待 SC 正式验证阶段。

## 9-11. 诊断结果（修复后：SC 生效）
smoke（500 fC 纯漂移 0.5m，同 seed）：SC OFF vs ON
- σx: 725.3 → 2436.1 µm（**+235.9%**）· σz: 302.0 → 1594.9 µm（+428%）· εnx: 0.080 → 0.120（+49%）· SC effect **PRESENT**
charge scan（0/10/50/100/500/1000 fC）：单调递增
- σx: 725→787→1004→1230→2436→3439 µm；σz: 302→357→543→726→1595→2269 µm；εnx: 0.080→0.149
convergence（500 fC）：
- N: 1e4/5e4/1e5 → σx 2444/2436/2438 µm（±0.3%，**收敛**）
- mesh: 33/63/127 → 2425/2436/2439 µm（±0.6%，**收敛**）
- SC step: 1/5/10 → 2436/2472/2620 µm（**step≥10 偏离 7.5%**，建议 step≤5）
实际 SC attrs：nmesh_xyz=[63,63,63], step=1, random_seed=10, random_mesh=False ✓

## 12. no-SC 回归
修复后 6/6 测试 PASS + R56 表征不变；full_beamline OCELOT 1996.205/474.022 与 v0.13 **逐位一致**（hash e041d6ae9fb7a0d2）；AG 位级不变。step4 运行时 σ_δ 2.94→4.65e-3（SC 生效的直接证据）。

## 13. SC 正式验证方案（修复后）
1. smoke（纯 drift SC off/on，同 seed）→ 验证 Δσ ≠ 0
2. charge monotonicity（0..1000 fC）
3. convergence（N/mesh/step）
4. AG vs OCELOT 四级比较（方向→趋势→量级→定量解释）

## 14. 结论
**P0 已定位**：不是 SC 公式或参数问题，而是**调用机制缺失**（tracking_step 不触发 PhysProc）。修复是纯数据流改动，不触碰任何 SC 物理。
