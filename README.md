# UED Beam Dynamics Simulator

超快电子衍射（UED）束流动力学模拟器 —— 基于物理的、可互验的双路线束流模拟框架。

- **AG 路线**：Kelisani 6D 束流包络方程（RMS 解析模型，完全确定性）
- **OCELOT 路线**：宏粒子追踪（SecondTM 传输矩阵 + 可选 PIC 空间电荷）

两条独立路线共享**唯一参数源** `shared/beamline_config.yaml`（lattice 单一来源），
并用一套统一的验证框架进行交叉校验（AG vs OCELOT）。

---

## 项目目标

以"物理可解释、数值可验证、过程可追溯"为原则，构建一个 UED 束流动力学模拟器：

- 横向：螺线管聚焦、漂移、束团尺寸演化
- 纵向：RF 纵向透镜（chirp）与束团压缩、时间分辨
- 空间电荷：高斯模型（Kelisani 2023）与均匀椭球模型（Luiten 2004）
- 双路线互验：AG（解析包络）与 OCELOT（宏粒子）在相同参数下交叉校验

> 信任一个模拟结果的前提是「物理模型 + 数值实现 + 验证 + 文档」四者一致。
> 一个物理假设清晰、经过验证的简单模型，胜过复杂的黑箱模拟。
> —— 详见 [AGENTS.md](AGENTS.md)（AI 辅助科研编码规范）

---

## 架构总览

```
beamline_config.yaml (唯一参数真源)
        │
        ├──────────────────────────────┐
        ▼                              ▼
   shared/params.py              (两条路线共用)
        │
   ┌────┴─────────┬───────────────────┐
   ▼              ▼                   ▼
 AG 路线        OCELOT/GPT 路线      validation 框架
 (解析包络)     (宏粒子追踪)          (交叉验证/回归)
```

### 两个独立 backend

| | AG（解析包络） | OCELOT（宏粒子） |
|---|---|---|
| 方法 | 6D 包络 ODE（`solve_ivp` RK45） | 传输矩阵（SecondTM）+ 粒子统计 |
| 模型 | Kelisani 2023 | OCELOT `generate_parray` + `Navigator` |
| 空间电荷 | Gaussian / Ellipsoid 解析模型 | PIC（NGP 沉积 + FFT 泊松求解） |
| 随机性 | 确定性算法，位级固定 | 配置化 seed，同 seed 位级可复现 |
| 发射度 | 输入常量（不产生发射度增长） | PIC 投影发射度（可能增长） |

比较原则：不要求逐点重合，按「方向 → 趋势 → 量级 → 定量解释」四级推进。

---

## 目录结构

```
├── shared/                       # 共享层（单一来源）
│   ├── beamline_config.yaml      # ★ 唯一参数真源（beam + lattice + SC）
│   ├── params.py                 # 参数解析与派生量
│   ├── constants.py              # 物理常数（单一来源）
│   ├── beam_physics.py           # γ/β/p0 派生（BeamReference）
│   ├── ocelot_coords.py          # rparticles 语义访问（禁止裸索引）
│   ├── output_schema.py          # 统一输出 schema
│   └── compare.py                # 双路线对比工具
│
├── AG/                           # AG 路线（解析包络核心，勿改）
│   ├── beam_dynamics_6d.py       # 6D 包络 ODE（1289 行，核心）
│   ├── external_forces.py        # 螺线管 / RF 力模型
│   ├── beamline_sim.py           # 出图
│   └── run_shared.py             # AG 统一输出适配层（薄适配器）
│
├── GPT模拟/                      # OCELOT 路线（主路由）
│   ├── ued_beamline_v2.py        # 主路由（lattice 单一来源，可 import）
│   └── {drift,solenoid,rf,space_charge}/   # 分节 benchmark
│
├── validation/                   # 验证框架
│   ├── backend.py                # 双后端适配（最高风险）
│   ├── beam_result.py            # 统一结果容器
│   ├── run_all.py                # 一键回归（6 项）
│   ├── test_*.py                 # 分节/整线/路由等价测试
│   ├── CHECKPOINTS.md            # 验证检查点日志
│   ├── reports/                  # 验证报告
│   └── baselines/                # 冻结基线
│
└── scripts/
    └── sync_to_phone.py          # iCloud 同步 + gpt-review 打包
```

---

## 环境要求

- macOS / Linux
- Python 3.10+
- 依赖：`numpy`、`scipy`、`matplotlib`、`PyYAML`、`ocelot`（≥ 26.x）

```bash
pip install numpy scipy matplotlib pyyaml ocelot
```

> 注意：OCELOT 的 `energy` 参数是**总能量**（GeV），即 `(E_keV + 511) * 1e-6`，
> 这是本项目曾经踩过的单位坑（见 `validation/CHECKPOINTS.md`）。

---

## 快速开始

### 1. 运行 OCELOT 主路由（GPT 路线）

```bash
# 支持 step 语义：1=drift, 2=+solenoid, 3=+RF, 4=+空间电荷
python3 "GPT模拟/ued_beamline_v2.py" --step 3
```

### 2. 运行 AG 路线（解析包络）

```bash
python3 AG/run_shared.py
```

### 3. 一键回归验证（约 30–40 分钟）

```bash
python3 validation/run_all.py            # 6 项测试 + 样品面数值
python3 validation/test_r56_convention.py # 冻结表征（结果必须逐位不变）
```

两条路线运行后都会写入统一结果文件（`shared/results/`），
可用 `shared/compare.py` 做 AG vs OCELOT 对比。

---

## 配置文件

所有物理参数的唯一来源是 `shared/beamline_config.yaml`。核心结构：

```yaml
beam:
  energy_keV: 100.0        # 束流动能 [keV]
  charge_fC: 100.0         # 束团电荷 [fC]
  n_particles: 50000       # 宏粒子数 (OCELOT) / Ne (AG)

initial_distribution:
  spot_rms_um: 85.0        # 横向 RMS 尺寸 [μm]
  bunch_length_um: 300.0   # 纵向 RMS 长度 [μm]
  epsilon_n_mm_mrad: 0.08  # 归一化横向发射度 [mm·mrad]
  sigma_delta: 1.0e-4      # 动量偏差 δ_p = Δp/p₀

space_charge:
  enabled: false           # 空间电荷默认关闭
  mesh: [63, 63, 63]       # PIC 网格

lattice:
  elements:                # ★ 唯一几何定义（z_start/length 单位：m）
    - {name: cathode, type: cathode, z_start: 0.000, length: 0.000}
    - {name: drift1,  type: drift,   z_start: 0.000, length: 0.100}
    - {name: solenoid1, type: solenoid, z_start: 0.100, length: 0.060,
       parameters: {B_field_T: 0.05}}
    - {name: rf1, type: rf_cavity, z_start: 0.400, length: 0.022,
       parameters: {frequency_GHz: 2.856, voltage_kV: 30.0, phase_rad: 3.1416}}
    - {name: sample, type: sample, z_start: 0.777, length: 0.000}
```

元件类型：`cathode`（阴极）、`drift`（漂移）、`solenoid`（螺线管）、
`rf_cavity`（RF 腔）、`sample`（样品面/诊断面）。支持多螺线管、多 RF 实例。

> 物理开关 `physics_switches`：`rf_longitudinal_kick`（默认开）、
> `rf_transverse_kick`（默认关，Panofsky-Wenzel 横向踢，公式尚未独立验证）。

---

## 坐标约定（最重要）

项目的一切改动都以坐标约定为前提，详见 [AI_HANDOFF.md](AI_HANDOFF.md) §4。

| 符号 | 含义 | 单位 |
|---|---|---|
| δ_p = Δp/p₀ | 项目标准动量偏差（无量纲） | - |
| p_oc = ΔE/(c·p₀) | OCELOT 原生第六坐标（**不是** Δp/p₀） | - |
| τ = c·t | OCELOT 共动时间坐标 | m |
| z = −β₀·c·Δt | 空间坐标（头为正） | m |

转换规则（仅限 adapter 边界）：

```
进入 OCELOT : p_oc = β₀ · δ_p
离开 OCELOT : δ_p  = p_oc / β₀
```

`BeamResult.sigma_delta_e3` 恒表示 δ_p。禁止在业务代码直接写
`rparticles[5]` 等魔法索引，必须用 `shared/ocelot_coords.py` 的语义访问。

---

## 物理模型

### 螺线管聚焦

聚焦强度：

```
k_s = e·B_z / (2·p₀)
```

本项目参考束 k_s = 22.38 m⁻¹（k_s² = 500.65 m⁻²）。
对于**圆束、无关联**的束团，AG 核心的降阶 Larmor 耦合项会产生虚假 σ_xy 并欠聚焦，
已在适配层关闭（`solenoid_coupling=False`），AG 核心未改，与精确硬边 4×4 传输一致。

### RF 纵向透镜（薄透镜模型）

两后端统一：

```
chirp  H = e·V·k·cosφ / (β²·E₀)      # = −9.78 m⁻¹
kick   δ_p += (V/(β²E_total))·sin(φ + k·z_phys)
```

φ = π 时产生最大线性 chirp（速度聚束）。每次每个 RF 实例在自身 z_start 踢一次。

### 空间电荷（SC）

- **AG**：Gaussian（Kelisani 2023，α 系数插值）或 Ellipsoid（Luiten 2004，三轴 Carlson 积分）
- **OCELOT**：PIC —— NGP 沉积 → 自由空间 Green 卷积 + FFT 解泊松 → 三线性插值 → 含纵向 E_z

> ⚠️ 当前状态（v0.14.1）：SC 调度已迁移到 OCELOT native `get_next_step()`，
> 但 **SC 数值尚未正式通过 beamline 验证**，不声称 fully validated。
> 已知未修问题：AG 的 SC 强度 ∝ Ne·e（宏粒子数），与物理束团电荷
> `charge_fC` 语义不一致（偏弱 12.5 倍，见 v0.14.1 task 2）。

---

## 验证与回归

一键回归 `validation/run_all.py` 共 6 项，全部 PASS 为基线：

| 测试 | 内容 | 验收 |
|---|---|---|
| config_schema | Level-1 配置一致性 | 能量↔γ↔β↔p₀ 自洽 |
| drift | 漂移传输 | σ_x 0.2%（解析参考） |
| solenoid | 螺线管聚焦（耦合关闭） | 0.40% |
| rf | RF 薄透镜 | σ_δ <1%、kick 语义 6e-13 |
| full_beamline | 整线 7 项量化 + 腰位置 | σ_x 0.60%、σ_z 0.63%、腰 Δz=0.4mm |
| gpt_route | 主路由 vs 验证路由 lattice 等价 | 样品面 <1% |

当前冻结数值（100 keV 参考束，no-SC）：

```
AG    : σ_x = 1984.191 μm   σ_z = 477.001 μm   （位级固定）
OCELOT: σ_x = 1996.205 μm   σ_z = 474.022 μm   （同 seed 位级可复现）
```

RNG 策略（v0.13）：`seed → x/y/tau`；`seed+1 → px/py/δ_p`（独立序列，
避免 x–px 虚假相关）。

**停止条件**：任何测试从 PASS 变 FAIL / AG 位级变化 / config SHA 意外变化 /
元件数量或顺序不一致 / RF kick 数量错误 / 出现 NaN——立即停止并报告第一个失败点。

---

## 版本状态

| 版本 | 内容 | 标签 |
|---|---|---|
| v0.10 | drift/solenoid/RF 分节验证通过，冻结基线 | `v0.10-noSC-longitudinal-validated` |
| v0.11 | lattice 单一来源（清除硬编码几何） | `v0.11-noSC-single-source-lattice` |
| v0.12 | 架构治理（常量单一来源、provenance） | - |
| v0.13 | 结构债务收敛（BeamReference、seed 配置化） | `v0.13-preSC-maintainability` |
| v0.14 | SC 接入审计：定位 P0（tracking_step 不触发 PhysProc）并修复 | - |
| v0.14.1 | SC 调度迁移到 OCELOT native `get_next_step()` | - |

当前阶段：v0.14.1 SC 正式 beamline 验证准备。长期目标：SC 验证完成后
才有资格讨论 AG/OCELOT 的 SC 物理对比、GUI、优化。

---

## 文档导航

| 文档 | 内容 |
|---|---|
| [AGENTS.md](AGENTS.md) | AI 辅助科研编码规范（物理笔记、验证、红线） |
| [AI_HANDOFF.md](AI_HANDOFF.md) | 项目当前状态唯一入口（版本、坐标、结论、任务书） |
| [validation/CHECKPOINTS.md](validation/CHECKPOINTS.md) | 验证检查点详细日志 |
| [validation/reports/](validation/reports/) | 各专题验证报告（SC/R56/架构/数据流） |
| [REFERENCES.md](REFERENCES.md) | 公开文献清单（PDF 不入库） |
| [OPTIMIZATION_PLAN.md](OPTIMIZATION_PLAN.md) | 物理模型优化路线图 |
| CHANGELOG_*.md | 各次修改报告（原因→修改→物理影响→验证） |

---

## 红线（禁止修改）

- AG 核心三文件（`beam_dynamics_6d.py` / `external_forces.py` / `beamline_sim.py`）
- OCELOT 安装源码
- `validation/test_r56_convention.py`（冻结表征）
- 测试阈值、baseline 目录、`shared/beamline_config.yaml` 物理参数值
- 随机策略、RF 方程、R56、螺线管方程

---

## 许可证

[MIT](LICENSE) · Copyright (c) 2026 QIN-SMART
