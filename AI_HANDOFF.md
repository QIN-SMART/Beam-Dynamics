# AI_HANDOFF — UED 项目当前状态唯一入口

> 本文档是**当前状态入口**，不是完整历史记录。新 session 从本节开始建立认知；
> 需要细节时再读 validation/CHECKPOINTS.md、validation/reports/、validation/baselines/。
> 生成时间：2026-08-10 · 本文档在 handoff 提交后需同步更新 HEAD 号。

## 0. 版本与位置
- 分支：`main`；HEAD：`ad43c69`（v0.14.1 task 1 提交；本文档提交后 HEAD 再变，以 git log 为准）
- 最近可信标签：**`v0.13-preSC-maintainability`**（上一冻结基线，可回退点）
- 其他标签：`v0.10-noSC-longitudinal-validated`（首次纵向冻结）、`v0.11-noSC-single-source-lattice`
- config SHA：以 `git show HEAD:shared/beamline_config.yaml` 实际计算为准；历史上 v0.13 加 `random.seed` 键导致 SHA 变化过一次，**物理参数值从未改变**
- 运行环境：macOS 13.7.1 (x86_64) / Python 3.12.2（conda-forge）/ OCELOT 26.06.1 / NumPy 1.26.4 / SciPy 1.13.1 / Matplotlib 3.9.2 / PyYAML 6.0.1
- **重要**：所有命令用 `/opt/anaconda3/bin/python3`；系统 `/usr/bin/python3`（3.6）缺 numpy/ocelot/yaml，不可用
- 远程：GitHub `QIN-SMART/Beam-Dynamics`（public，token 已存钥匙串；push 若遇 HTTP 400/超时用 `git -c http.version=HTTP/1.1 -c http.proxy= push`）

## 1. 当前阶段 / 下一阶段
- 当前阶段：**v0.14.1 SC 正式 beamline 验证准备**——task 1（scheduler equivalence）已完成：manual counter 经 characterization 证明非一般等价（T5–T7），production 已迁移到 OCELOT native `get_next_step()` 调度（SC ON），SC OFF 位级不变；SC 覆盖修正为 cathode→sample（0→0.777m，stop anchor 来自 lattice.elements 的 sample marker）；**SC 数值仍未经正式 beamline 验证，未声称 fully validated**
- 当前阶段：**v0.14.1 三个代码硬化任务全部完成**（scheduler equivalence / AG charge semantics / SC runtime state contract，见 §10/§11/§12/§14）；**待冻结 `v0.14.1-SC-integration-hardened`**；v0.15（native SC 数值验证、收敛固化、AG vs OCELOT SC 物理对比）未开始
- 长期目标：SC 验证完成后才有资格讨论 AG/OCELOT 的 SC 物理对比、GUI、优化
- 禁止事项：修改 AG/OCELOT 核心、RF 方程、R56、lattice 几何、测试阈值、随机策略；未获批准不得实施 v0.14.1 task 2/3 或任何 SC 参数调整

## 2. 项目架构
```
shared/           参数解析 params.py；常量 constants.py（单一来源）；
                  γ/β/p0 派生 beam_physics.py；OCELOT 坐标封装 ocelot_coords.py；
                  输出 schema output_schema.py；双路由对比 compare.py
AG/               核心：beam_dynamics_6d.py（6D 包络 ODE 1289 行）、
                  external_forces.py（螺线管/RF 力）、beamline_sim.py（出图）；
                  适配：run_shared.py（统一输出 JSON）
GPT模拟/          主路由 ued_beamline_v2.py（lattice 单一来源，可 import 测试）
validation/       backend.py（双后端适配，最高风险）、beam_result.py、
                  reference.py（解析参考）、config_check.py（Level-1 检查）、
                  test_*.py、reports/、baselines/
scripts/          sync_to_phone.py（iCloud 同步 + bundle + gpt_review 版本归档）
```
- **唯一几何真源**：`shared/beamline_config.yaml → lattice.elements`；主路由与验证路由均从它构建（v0.11 达成，有 `test_gpt_route_equivalence.py` 保证等价）
- 统一输出：`BeamResult`（σ_x/y/z_um、σ_t、σ_δ_p、ε_nx/ny_mm_mrad + meta provenance）
- provenance meta 含：git_commit、config_sha、lattice_hash、timestamp、python、random_seed、坐标约定

## 3. 两个 backend 职责
- **AG（解析包络，确定性）**：`validation/backend.py::run_ag` → `make_beam_100keV`（Ne、能量、σ 矩、发射度）→ 按 lattice 建 `ExtFieldRegion` 力区（螺线管 Bz、RF 横向力）→ `propagate()`（solve_ivp RK45 积分 6D 包络 ODE）→ 输出 σ 矩。RF 纵向用薄透镜 `apply_rf_thin_lens(H)` 分段施加。SC 通过 `space_charge_forces(beam, model)`（gaussian/ellipsoid）进入 ODE。**AG 是确定性算法，无随机数，SC OFF 时结果位级固定**
- **OCELOT（宏粒子）**：`run_ocelot` → `generate_parray`（配置 seed）→ Drift/Solenoid（SecondTM）→ 解析 RF kick → 可选 SpaceCharge（counter 机制触发）→ 粒子统计 → BeamResult
- 两者比较原则：不要求逐点重合；按"方向→趋势→量级→定量解释"四级推进（AG 是 RMS 解析，OCELOT 是 PIC，模型能力不同）

## 4. 坐标约定（最重要，改代码前必读）
- **项目标准变量：δ_p = Δp/p0**（动量偏差，无量纲）
- **OCELOT 原生第六坐标：p_oc = ΔE/(c·p0)**（能量偏差归一化，**不是** Δp/p0；源码证据 `ocelot/cpbd/beam/generator.py:51`、`beam/particle.py:20`）
- 两者换算（仅限 adapter 边界）：**p_oc = β0·δ_p**（进入 OCELOT 时），**δ_p = p_oc/β0**（离开时）
- 允许出现 p_oc 的位置只有两处：`validation/backend.py`（run_ocelot 束生成、_ocelot_rf_kick）和 `GPT模拟/ued_beamline_v2.py`（apply_rf_kick）——除此之外任何代码出现 p_oc 都是违规（已登记例外：只读诊断脚本 `validation/sc_audit_diagnostics.py:62` 用 set_p_oc 镜像束生成，属 manual-scheduler characterization 保留代码）
- tau = c·t [m]（OCELOT 原生，rparticles[4]）；共动空间坐标 z = −β0·c·Δt（头为正）；输出 σ_z = β0·σ_tau；σ_t = σ_z/(β0c)
- **BeamResult.sigma_delta_e3 恒表示 δ_p**（OCELOT 输出已 ÷β0），禁止直接暴露 raw p_oc
- 业务代码禁止直接写 `rparticles[5]` 等魔法索引，必须用 `shared/ocelot_coords.py` 的语义访问（I_X/I_PX/I_Y/I_PY/I_TAU/I_P、set_*/add_*）
- 禁止在 AG 核心出现 p_oc（已用 grep 验证 0 处）
- 磁矩语义：γ/β/p0 统一从 `shared/beam_physics.py::BeamReference.from_energy_keV` 派生（v0.13），禁止各处重复手写 1+E/511

## 5. R56 最终结论
- 分类：**B — 输入 δ 变量约定不匹配**（不是 OCELOT bug，也不是 AG bug；两模型各自内部自洽）
- 关键证据（冻结表征 `test_r56_convention.py`，禁止修改）：OCELOT raw 斜率 = 传输矩阵 R56_tm = **−1.163624**，与 `−L/(β²γ²)` 偏差 1.25e-5，而与精确 c·t 斜率 `−L/(βγ²)=−0.6379` 差 82.4%（=1/β−1）；形式转换（δ_p=p_oc/β0，Δz=−β0·Δτ）闭合残差 0.2%
- 修复（v0.13 前完成）：adapter 输入 ×β0、RF kick ×β0、输出 ÷β0；**OCELOT 原生 R56 未改**
- 结论：σ_z 压缩曲线在样品面一致（<1%）；深压缩区腰附近有几 µm 的高阶数值残差（τ 空间 vs z 空间 RMS 演化差异），已解释，非 bug

## 6. Solenoid 最终结论
- AG 降阶 Larmor 耦合项（dν_x⊃−2k_s·ν_y、dν_y⊃+2k_s·ν_x、dσ_xy⊃2k_s(σ_x²−σ_y²)）对**圆束**产生虚假 σ_xy 并欠聚焦（z=160mm 处 σ_x=147µm vs 精确 97µm，x-y 对称破坏 σ_y=2469）
- 精确证据：硬边 4×4 Brown-Chao 矩阵（= OCELOT SolenoidTM）对圆束给出 σ_xy≡0（解析可证），OCELOT 结果与之一致
- 修复：适配层关闭耦合（`run_ag(..., solenoid_coupling=False)`；`AG/run_shared.py` 默认 False）；**AG 核心未改**，修正理由写入 CHECKPOINTS
- 验证：coupling OFF 后 AG vs OCELOT σ_x 偏差 0.40%（样品面），x-y 对称恢复

## 7. RF 状态
- **纵向（已解决）**：两后端统一薄透镜模型；chirp H = eV·k·cosφ/(β²E₀) = **−9.78 m⁻¹**；kick δ_p = (V/(β²E_tot))·sin(φ+k·z_phys)，z_phys=−β0·tau；σ_δ_p 一致性 0.3–1%（阈值 2%）；kick 语义测试残差 6e-13
- **横向（open item）**：K_trans = −e·k·V/(2γβ·m_e·c²) = **−2.68 m⁻¹**（Panofsky-Wenzel），默认 **OFF**（`physics_switches.rf_transverse_kick: false`）；公式尚未独立验证，仅两后端自洽。开启时两路一致（AG σ_x 1119→542、OCELOT 1120→535）
- kick 时机按每个 rf 实例的 z_start（多实例支持），每次只踢一次；drift/solenoid 节禁止误踢（有 routing 测试防回归 {0,0,1,1}）
- 参数来源：`lattice.elements[rf*].parameters`（frequency_GHz/voltage_kV/phase_rad）

## 8. Lattice single source 状态
- 已达成：`build_lattice_from_shared(cfg, active_types)` 是主路由唯一几何构建入口；验证路由 run_ocelot 从同一 `_lattice_elements` 构建
- step 语义（**测试语义，不是物理路由**）：1{drift} / 2{+solenoid} / 3{+rf} / 4{+rf+SC}；未激活元件保留等长 Drift，总长恒等于样品位（777mm）
- 多螺线管、多 RF 支持（每实例独立 B/V/f/φ）；几何等价由 `test_gpt_route_equivalence.py` 验证（名称/类型/长度/顺序/RF 位置精确一致，样品面 dev<2%）
- 配置一致性由 `validation/config_check.py`（Level-1）只读检查：energy↔γ↔β↔p0、RF 频率合理、位置单调、sample==总长、schema 合法

## 9. RNG Policy
- 配置：`config.random.seed = 42`（v0.13 起）
- 生成序列：`seed` → generate_parray（x/y/tau 确定）；`seed+1` → px/py/δ_p（独立序列，**避免 x 与 px 用同序列导致虚假相关的陷阱**——这是 v0.13 之前发现的坑，禁止回到单 seed 覆盖）
- **同 seed 两次运行位级一致**（OCELOT full 验证 hash e041d6ae9fb7a0d2）——历史"首束随机"问题已解决
- 陷阱：OCELOT `SpaceCharge.prepare()` 无条件 `np.random.seed(10)`（sc.py:106）——若在束生成前 prepare 会污染全局 RNG；当前流程 prepare 在束生成后，影响有限，但记录在案
- 测试受控集用 `default_rng(seed+2)`（test_rf），与主序列隔离

## 10. SC v0.14 审计最终结论
- **P0 根因（已修复）**：OCELOT 的 `tracking_step()` 只应用 transfer maps；PhysicsProcess（含 SpaceCharge）只在 `track()` 主循环的 `get_next_step()` 的 counter 机制中触发（`ocelot/cpbd/track.py:476-479`）。**历史所有代码（backend、主路由、benchmark）都用 tracking_step，因此 SC 从未真正运行**——这就是 v0.12 之前"SC ON/OFF 几乎重合"疑点的根因
- **v0.14 修复（manual counter，已退役）**：曾在 `run_ocelot`、`ued_beamline_v2.run_beamline`、`sc_audit_diagnostics` 的循环内复刻 counter 逻辑（`sc.counter -= 1; if counter<=0: ...`）；SpaceCharge 构造从 config 传 `step` 与 `nmesh_xyz`
- **v0.14.1 task 1（已落地）**：manual counter 经 characterization（`validation/test_sc_scheduler_equivalence.py`）证明非一般等价——T5 尾段丢失、T6 区间外触发、T7 生产 lattice 下覆盖 [0,0.422) vs [0,0.777)；**production 已迁移到 OCELOT native `get_next_step()` 调度**（== track() 核心，backend.py::run_ocelot 与 ued_beamline_v2.py::run_beamline 的 SC ON 分支）；SC OFF 保持原 tracking_step 循环，位级不变（数组 hash 7790fd9c2a2b）
- **stop anchor**：SC ON 时 runtime sequence 保留 lattice.elements 中的 cathode/sample 零长 marker，PhysProc attach 为 cathode→sample（覆盖 [0, sample.z_start=0.777)，无硬编码）；SC OFF sequence 不变（带/不带 marker 已验证位级一致）
- **metadata（SC ON，只读）**：sc_scheduler="ocelot_native"、sc_apply_count、sc_coverage_start_m、sc_coverage_stop_m、sc_events（z,zstep 列表）
- **⚠️ SC 数值未验证**：v0.14 smoke/charge/收敛数值是 **manual-scheduler 行为**（保留于 sc_audit_diagnostics.py 作为 manual characterization），native-scheduler 的 SC 数值验证属 v0.15。**禁止把 SC 描述为 fully validated**
- 修复后诊断（`sc_audit_diagnostics.py`）：smoke 500 fC 纯漂移 0.5m：σx 725→2436µm（+236%）、σz 302→1595µm（+428%）、εnx 0.080→0.120（+49%）；charge 0→1000 fC 单调（σx 725→3439、σz 302→2269、εnx→0.149）；N(1e4/5e4/1e5)±0.3%、mesh(33/63/127)±0.6% 收敛；SC step≤5 建议（step=10 偏 7.5%）
- OCELOT SC 算法（源码证据 sc.py）：NGP 沉积（:191-193）→ 自由空间 Green 卷积 + FFT 解 Poisson（:85-168，ASTRA 同款）→ 三线性插值（:202-204）→ 含纵向 Ez（:200,248）→ 束团静止系解场、实验室系 kick（横向乘 1−β0²，:246-248）；电荷来自 q_array（:241）；默认 low_order_kick=True、random_mesh=False、random_seed=10

## 11. AG charge semantics bug（已修复，v0.14.1 task 2）
- 根因：AG SC 强度 `fb = η·Ne·e/(8π√π·ε₀)`（beam_dynamics_6d.py:379）中 **Ne = config beam.n_particles = 50000** → 等效总电荷 Ne·e = **8 fC**，比物理束团电荷（`beam.charge_fC = 100`）弱 **12.5 倍**
- 语义区分：**beam.n_particles 对 OCELOT 是宏粒子数值计数（数值精度参数）；beam.charge_fC 是物理束团电荷**
- **修复（2026-08-12，仅 adapter 层）**：`validation/backend.py::run_ag` 与 `AG/run_shared.py` 注入 `Ne_phys = abs(Q_C)/E_SI`（≈6.24e5 @100 fC）；AG 核心与 SC 公式未改；SC OFF 路径 Ne 强制 0 → 位级不变
- 验证：`validation/test_ag_charge_semantics.py`——charge 语义（ag_ne_phys==Q/e，σx(Q) 单调 2202→12106 µm）与 n_particles 不变性（Q 固定时 n=1e4/5e4/1e5 SC ON 位级一致；SC OFF sample σx=1984.191/σz=477.001 与 v0.13 基线精确一致）；run_all 6/6 + r56 不变
- 注意：SC ON 数值较旧 bug 时代明显增大（Q=100 fC 时 σx≈3657 µm vs 旧 ~1 mm 量级）——**预期变化，不是回归**

## 12. SC runtime state contract（已实施，v0.14.1 task 3）
- **统一状态机**（`shared/sc_state.py::SCState`）：`sc_requested → sc_available → sc_configured → sc_attached → sc_apply_count → sc_effective`；`sc_effective = (sc_apply_count > 0)` 是"SC 真正运行"的唯一判据
- **HARD FAIL**（requested=True 时）：import/构造/attach 失败、tracking 后 apply_count==0、coverage != cathode→sample 全部直接 raise；`except ImportError: pass` 已删除（backend.py 与 ued_beamline_v2.py）
- **双状态源消除**：`step` 只定义 route capability（1-3 无 SC），config.enabled 定义 requested；GPT 主路由 `sc_requested_from(cfg,step)=config.enabled AND step>=4`，main 单次运行（r_uni 第二跑已删除）→ 终端显示与 saved sc_enabled 不可能不一致
- **AG meta 对齐**：physical_charge_C / physical_electron_number（正式 contract）+ sc_requested/sc_effective；ag_ne_phys 保留 alias（task 2 兼容）
- 测试：`validation/test_sc_runtime_state.py` A-F 全 PASS；回归 6/6 + r56 + canonical hash 7790fd9c2a2b 不变

## 13. no-SC regression baseline（必须保持）
- 测试命令：`/opt/anaconda3/bin/python3 validation/run_all.py`（6 项）+ `validation/test_r56_convention.py`
- **canonical no-SC regression hash（v0.14.1 起正式定义）**：
  `SHA1( contiguous bytes of [z_mm, sigma_x_um, sigma_y_um, sigma_z_um, eps_nx_mm_mrad, eps_ny_mm_mrad, sigma_delta_e3] )[:12]`
  （run_ocelot SC OFF，config seed 42，N=5e4，dz=0.001，full beamline）→ **`7790fd9c2a2b`**
  旧 `e041d6ae9fb7a0d2` 不再追踪：其输入对象/算法未记录，无法严格复现
- 6 项全 PASS：config_schema（Level-1 一致性）、drift（σ_x 0.2%、σ_δ_p 语义 0.022%）、solenoid（0.40%）、rf（σ_δ_p <1%、kick 语义 6e-13）、full_beamline（七项量化：σ_x 0.60%、σ_y 0.36%、σ_δ_p 0.46%、ε_nx 0.53%、ε_ny 0.26%、σ_z 0.63%、σ_t 0.63%、腰 Δz=0.4mm）、gpt_route（样品面 <1%）
- full_beamline 样品面参考值：AG σ_x=1984.191µm/σ_z=477.001µm（位级固定）；OCELOT σ_x=1996.205/σ_z=474.022（同 seed 位级，hash e041d6ae9fb7a0d2）
- r56_convention：raw slope=−1.163624、闭合残差 2.11e-3、naive 0.826 —— 冻结不变
- AG 数组位级不变；OCELOT 同 seed 位级可复现

## 14. 下一阶段 v0.14.1 明确任务（SC 正式 beamline 验证）
1. **✅ 已完成（2026-08-12）**：manual counter vs OCELOT native scheduler 等价性证明。结论：非一般等价（T5–T7 差异），production 已迁移到 native `get_next_step()`（SC ON），SC OFF 位级不变；SC 覆盖修正为 cathode→sample。证据：`validation/test_sc_scheduler_equivalence.py`（T1–T7 + acceptance A–F 全 PASS）、run_all 6/6、r56 不变
2. **✅ 已完成（2026-08-12）**：AG charge semantics——`run_ag`/`run_shared` 注入 `Ne_phys = Q/e`（≈6.24e5）；SC OFF 位级不变（1984.191/477.001）；n_particles 不变性测试全 PASS（详见 CHECKPOINTS v0.14.1-ag-charge-semantics）
3. **✅ 已完成（2026-08-12）**：SC 状态机（六态 + HARD FAIL，shared/sc_state.py）与双状态源统一（step=capability、config=requested、effective=apply_count>0）；测试 test_sc_runtime_state.py A-F 全 PASS
4. 固化 SC 收敛参数：N=5e4、mesh=63³、SC step≤5（文档化到 config 注释与报告）——**属 v0.15**
5. AG vs OCELOT SC 四级比较：方向（SC ON/OFF 一致）→ 趋势（charge 单调一致）→ 量级 → 定量解释；**禁止逐点强求一致**；ε 增长差异（AG 天然不增长 vs OCELOT PIC 投影增长）是模型能力差异，**不可调参消除**——**属 v0.15**
7. **物理证据链审计（v0.15 新增，来源：GPT 审核意见 2026-08-16）**：为每个 AG SC 公式建立「公式 → 原始论文 → 代码实现行号 → 独立验证」证据链。已知外部锚点：Kelisani et al. 2023 PRA（6D envelope + Gaussian SC，f_b ∝ Ne）、Luiten et al. 2004 PRL（均匀椭球 SC，线性内部场）、van Oudheusden et al. 2010 PRL（95 keV/200 fC RF 压缩，UED 同量级基准）。目标：排除"测试全 PASS 但公式源于 AI 生成"的风险（测试只证明代码符合已规定行为，不证明规定本身正确）。验证方式：复现论文关键图 + 与 OCELOT PIC 独立对比
6. 每步回归：no-SC 6 项 + r56 保持 PASS；SC ON 结果附完整 provenance（含 seed/mesh/step/N）
- **禁止**：声称 SC fully validated；修改 AG/OCELOT 核心；调整测试阈值；进入 GUI；优化参数让两模型吻合

## 15. 关键文件（职责速查）
- `validation/backend.py`：双后端适配 + 力组装 + RF kick + p_oc 转换 + provenance（**最高风险，改前必读 §4**）
- `shared/params.py`：解析 + derived() + lattice helpers；`shared/beam_physics.py`：BeamReference；`shared/constants.py`：物理常数；`shared/ocelot_coords.py`：rparticles 语义访问
- `GPT模拟/ued_beamline_v2.py`：主路由（可 import；main() 守卫）；`AG/run_shared.py`：AG 统一输出
- `validation/beam_result.py`：统一容器；`validation/reference.py`：解析参考；`validation/config_check.py`：只读一致性
- 测试：`test_config_consistency / test_drift / test_solenoid / test_rf / test_full_beamline / test_gpt_route_equivalence`；**`test_r56_convention.py` 冻结禁改**；`test_sc_scheduler_equivalence.py`（SC 调度 characterization + production acceptance A–F，mesh 33³ 加速）；`sc_audit_diagnostics.py`（manual-scheduler SC characterization 只读脚本）

## 16. 关键报告
- `validation/reports/SC_integration_audit.md`：P0 根因 + 修复 + 诊断数值 + 状态机设计（§8 附 v0.14.1 manual→native 实施状态更新）
- `validation/reports/SC_ocelot_source_audit.md`：OCELOT SC 源码逐条行号证据
- `validation/reports/SC_AG_model_audit.md`：AG SC 公式与模型假设
- `validation/reports/R56_convention_resolution.md`：B 类结论全推导
- `validation/reports/v0.12_dataflow_audit.md` / `v0.12_data_contract.md` / `v0.12_maintainability_risk.md` / `v0.13_architecture_audit.md`：架构/数据流/契约/风险
- `validation/baselines/v0.10-noSC-longitudinal-validated/`：冻结基线（manifest/environment/physics-status）

## 17. 禁止修改区域（红线）
- AG 核心三文件（beam_dynamics_6d.py / external_forces.py / beamline_sim.py）
- OCELOT 安装源码（/opt/anaconda3/.../ocelot/）
- `test_r56_convention.py`（冻结表征）
- 测试阈值；baseline 目录；`shared/beamline_config.yaml` 物理参数值
- 随机策略（除非作为已验证的接口改动）；RF 纵/横向方程；R56；solenoid 方程

## 18. STOP CONDITIONS（触发即停，报告第一个失败点，不连续打补丁）
- 任何原有测试从 PASS 变 FAIL；AG 确定性数组非预期变化；config SHA 意外改变
- 元件数量/顺序不一致；RF kick 数量错误；step1/2 出现 RF kick
- 出现 NaN、负方差、发射度≈0；为使两模型吻合而准备调参；无法确认某变量定义或单位
- SC 相关：SC ON 无任何效应；manual counter 与 native scheduler 对照不等价且无法解释

## 19. 交接后第一件事
1. `git log --oneline -5` + `git status` 确认工作树干净
2. 读本文件 §0-§14；如需细节读 validation/CHECKPOINTS.md 最新条目（注：CHANGELOG/CHECKPOINTS 目录在 validation/ 下）
3. 跑 `validation/run_all.py` + `test_r56_convention.py` 确认基线
4. 按 §14 任务书开始 v0.14.1 的第一步（scheduler 等价性证明），每步带回归

## 20. 阶段历史沿革（速览，细节在 validation/CHECKPOINTS.md）
- v0.10（基线冻结）：drift/solenoid/RF 分节验证通过，冻结标签 + baselines/ 目录；当时发现并修了 OCELOT 能量单位 bug（energy 应是总能量 GeV）、ε_nz 一致性、螺线管耦合、RF 薄透镜统一
- v0.11（lattice 单一来源）：主路由硬编码几何（0.100/0.240/0.355/0.777）清除，`build_lattice_from_shared` + 多实例 + step 语义；新增 `test_gpt_route_equivalence`
- v0.12（架构治理）：常量单一来源 constants.py、provenance meta、Level-1 config_check、7 份 v0.12 审计报告
- v0.13（结构债务收敛）：BeamReference 统一 γ/β/p0、ocelot_coords 语义封装、random.seed 配置化（同 seed 位级可复现）；标签 `v0.13-preSC-maintainability`
- v0.14（SC 接入审计）：定位 P0（tracking_step 不触发 PhysProc）、manual counter 修复、SC 诊断（smoke/charge/收敛）、三份 SC 报告；**SC 未正式验证**

## 21. 常见陷阱与调试要点
- **单位**：OCELOT energy 参数是总能量 GeV（`(E_keV+511)*1e-6`）；tau 是米（c·t）；σ_z 输出 = β0·std(tau)；σ_δ_p = std(p_oc)/β0。任何"看起来差 1.8 倍"先检查 β0 是否漏乘/漏除
- **seed**：束生成必须 seed→generate_parray、seed+1→px/py/δ_p；禁止在 generate_parray 前后重复同 seed（x-px 相关陷阱）；SC 的 prepare() 会 seed(10)（影响全局 RNG）
- **rparticles**：一律用 ocelot_coords 语义名；禁止裸索引
- **SC 无效症状**：SC ON/OFF 逐位相同 → 检查 counter 触发逻辑是否在循环内；若用 get_next_step() 注意无 proc 时按元素大步走（破坏细采样）
- **AG SC 偏弱**：Ne 语义是宏粒子数，AG SC 强度 ∝ Ne·e，与 charge_fC 无关（12.5× 偏弱，未修）
- **调试顺序**（AGENTS.md 规范）：先查单位 → 再查坐标定义 → 再查极限情形（SC OFF 应 ε=const、无场应 drift）→ 再对照解析解
- 物理常数一律从 shared/constants.py 取；γ/β/p0 一律从 BeamReference 取

## 22. AG SC 公式速查（详情见 SC_AG_model_audit.md）
- Gaussian（Kelisani 2023）：主项 Fs_x = fb·αx/(β²γ³·σx·σz)，fb = η·Ne·e/(8π√π·ε₀)，η=e/(mc²)；α 系数由 kx=σx/(γσz)、ky=σy/(γσz) 查插值表；另有 O(1/γ) 括号修正项（beam_dynamics_6d.py:384-431）
- Ellipsoid（Luiten 2004）：Fs_x = Ne·r_e·Mx/(γ³β²·σx·σz)，Mx/My/Mz 为均匀椭球 form factors（球对称解析 / 三轴 Carlson 积分），Mx+My+Mz=1（:506-545）
- 两模型 ε_nx 均为输入常量 → **天然不产生发射度增长**（这是模型能力限制，不是 bug；与 OCELOT PIC 的投影 ε 增长对比时禁止据此判 AG 错）
- SC 与发射度项在包络方程中线性相加：dν_u = Fe_u + Fs_u + F_eps_u + 耦合（envelope_ode:997-999）
- SC OFF 语义：`Ne=0` → fb=0 → 无 SC 力（no-SC 基线的标准做法）

## 23. 关键数值速查（100 keV 参考束）
- γ0=1.1957、β0=0.5482、βγ=0.6556、p0=1.790e-22 kg·m/s、v0=1.644e8 m/s
- ε_nx=ε_ny=0.08 mm·mrad、ε_nz=0.02 mm·mrad（=βγ·σ_z·σ_δ 自洽值）、σ_x0=σ_y0=85µm、σ_z0=300µm、σ_δ_p=1e-4、Q=100 fC、N=5e4
- 螺线管 k_s=22.38 m⁻¹（k_s²=500.65 m⁻²）、RF H=−9.78 m⁻¹、K_trans=−2.68 m⁻¹、k_rf=59.86 m⁻¹
- 样品面（no-SC）：AG σ_x=1984µm/σ_z=477µm；OCELOT σ_x=1996µm/σ_z=474µm
- SC 诊断（500 fC 漂移）：OCELOT σx 725→2436µm；AG 侧尚未做同等 SC 对比（Ne 语义未修，不可直接比）

## 24. 测试命令与协作约定
- 一键回归：`/opt/anaconda3/bin/python3 validation/run_all.py`（约 30-40 分钟，输出 6 项 PASS/FAIL + 样品面数值）
- 冻结表征：`... validation/test_r56_convention.py`（结果必须逐位不变）
- SC 调度：`... validation/test_sc_scheduler_equivalence.py`（native vs manual 事件对照 + production 验收 A–F，约 3 分钟）
- SC 诊断（manual characterization）：`... validation/sc_audit_diagnostics.py`（smoke + charge scan + convergence，只读）
- iCloud 同步：`... scripts/sync_to_phone.py --gpt-review`（任务收尾自动跑；A 级归档到 gpt_review/versions/<时间戳>/，iCloud 侧一律 .txt）
- 任务完成输出规范：按根目录 AGENTS.md 的审核输出格式（一句话结论、修改文件表、关键代码、数据流、验收数值表、风险、审核材料、可复制 GPT 摘要）；每轮任务结束自动跑 gpt-review 同步
- 停止条件：任何测试失败 / AG 位级变化 / config SHA 意外变化 / 元件顺序或长度不一致 / RF kick 数量异常 / NaN / 为拟合而调参——立即停止并报告第一个失败点
- 与用户协作：物理假设、模型选择、验收标准由用户决定；AI 负责实现、文档、测试；禁止未经批准修改红线区

## 25. 交接自检清单（下个 session 开场必做）
- [ ] `git log --oneline -3` 确认 HEAD 与本文 §0 一致（或已更新）
- [ ] `git status --short` 工作树干净（handoff 提交后应无未提交改动）
- [ ] `git rev-parse v0.13-preSC-maintainability` 可回退点存在
- [ ] 读 §4 坐标约定（δ_p/p_oc/tau/z 语义）——所有改动的前提
- [ ] 读 §10 SC 结论——**禁止声称 SC fully validated**
- [ ] 读 §14 v0.14.1 任务书——按序执行，第一步是 scheduler 等价性证明
- [ ] 跑 run_all + r56 确认基线，与 §13 数值对照
- [ ] 任何新代码遵守红线（§17）与 AGENTS.md 规范

## 26. 与旧 session 的边界
- 上一 session 的最后结论：SC 审计完成、P0 修复落地（manual counter）、SC 诊断通过（smoke +236% 等）、no-SC 六项回归 + R56 全部 PASS、AG 位级不变
- 上一 session 的已知债务（v0.14.1 处理）：manual counter 等价性未证、AG Ne 语义未修、SC 状态机未建、SC 收敛参数未固化
- **v0.14.1 task 1/2/3（本 session）全部完成**：scheduler 迁移（§10）、AG charge semantics（§11）、SC runtime state contract（§12）
- 本 session 未做的事：冻结标签 v0.14.1-SC-integration-hardened（待用户确认）、SC 正式 beamline 对比（v0.15）、GUI、优化、任何 SC 参数调整
