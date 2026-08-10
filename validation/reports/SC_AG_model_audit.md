# SC — AG Model Audit

Date: 2026-08-10

## 模型假设 vs 代码实现（明确区分）

### 模型假设（Kelisani 2023 / Luiten 2004）
- **Gaussian**：6D 相关高斯分布；SC 力由 9 个 α 系数（形状因子）解析表达。
- **Ellipsoid**：均匀椭球（Luiten 2004）；内部电场线性，form factors Mx+My+Mz=1。
- 两者都是 **RMS/包络级** 模型：SC 力只依赖 σ 矩，不追踪粒子间非线性相互作用。

### 代码实现

**1. SC 数学公式（beam_dynamics_6d.py）**
- Gaussian（`:384-431`）：
  - 主项 `Fs_x = fb·αx / (β²γ³·σx·σz)`（`:392`），fb = η·Ne·e/(8π√π·ε₀)（`:379`），η=e/(mc²)
  - 高阶括号项 A_x/B_x（`:395-404`，含 ν²、θ²、p0² 修正，O(1/γ) 项 `-fb/γ·(A+B)`）
  - Fs_y 对称（`:407-417`），Fs_z = fb·αz/(β²γ³·σx·σy)（`:420`）+ 括号项
- Ellipsoid（`:506-545`）：
  - `common = Ne·r_e/(γ³β²)`（`:537`）
  - `Fs_x = +common·Mx/(σx·σz)`，`Fs_z = +common·Mz/(σx·σy)`（`:541-543`）
  - form factors：球对称解析（`:462-479`）/ 三轴 Carlson R_D（`:480-501`）

**2. generalized perveance / charge 依赖**：两者都线性 ∝ Ne（高斯 fb∝Ne；椭球 common∝Ne）。100 fC/50k 粒子 → Ne·e 即总电荷。

**3. αx αy αz 计算**：`compute_all_alpha(kx,ky)`（`:77`）数值积分（s 代换 t/(1-t)，`_integrate_alpha`），kx=σx/(γσz)、ky=σy/(γσz)；表驱动插值 `alpha_at`（缓存 `get_alpha_interpolators`）。

**4. Gaussian/ellipsoid 假设**：两模型各自假设，通过 `space_charge_forces(beam, model)` 分派（`:548-565`）。

**5/6. 横向/纵向项**：见上 Fs_x/y/z；纵向 Fs_z 存在且非零。

**7. γ/β 依赖**：两模型都 ∝ 1/(β²γ³)（主项）。

**8. 与 emittance 项共同进入 ODE**：`envelope_ode`（`:965`）`Fs_x,Fs_y,Fs_z = space_charge_forces(beam, model)`；`dν_u = Fe_u + Fs_u + F_eps_u + 耦合`（`:997-999`）——SC 与发射度压力线性相加。

**9. SC 是否产生 emittance 增长**：**RMS 模型天然不产生**——ε_nx 是输入常量（theta_x 恒用 ε_nx/p0），SC 只改变 σ 演化；**这是模型能力限制**（见十一节原则）。

**10. 模型无法描述的效应**：非线性 SC、相空间畸变、PIC 噪声、投影发射度增长、非椭球/非高斯分布演化、halo。

## 代码质量备注
- Gaussian 括号项（A/B）带 `(1-p0²)` 因子——p0=γβ 归一化动量，注意与 δ_p 区分。
- 分母 clip `max(...,1e-30)` 防 NaN（`:385-388`）。
- ellipsoid 与 gaussian 的 SC 符号：ellipsoid 显式 `+`（发散，`:541-543`）；gaussian 主项 `+fb·αx/denom` 且 `-fb/γ(A+B)`——α>0，主项为正（发散）。
