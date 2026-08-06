# 文献清单 / Reference Literature

本清单列出项目依赖的公开文献。PDF 原件**不入库**（保持仓库轻量），
仅此清单随仓库同步。文献按主题分组；完整元数据以 PDF 为准。

## 包络方程 / Envelope equations

1. **Kelisani, Barzegar, Craievich, Doebert** —
   "Six-Dimensional Beam-Envelope Equations: An Ultrafast Computational
   Approach for Interactive Modeling of Accelerator Structures",
   Phys. Rev. Applied **19**, 054011 (2023).
   → AG 6D 包络 ODE 的核心理论（`AG/beam_dynamics_6d.py`）。
   File: `AG/文献/Kelisani_2023_Six-Dimensional_Beam-Envelope_Equations.pdf`

## 均匀椭球束团 / Uniform ellipsoidal bunches

2. **Luiten, van der Geer, Moors, de Loos** —
   "How to Realize Uniform Three-Dimensional Ellipsoidal Electron Bunches",
   Phys. Rev. Lett. **93**, 094802 (2004).
   → AG 空间电荷椭球模型（`sc_model='ellipsoid'`）。
   File: `AG/文献/Luiten_2004_Uniform_Ellipsoidal_Bunches_PRL93.pdf`

## 螺线管与 RF 耦合发射度 / Coupled transverse dynamics

3. **Dowell** —
   "Exact cancellation of emittance growth due to coupled transverse
   dynamics in solenoids and RF linacs",
   Phys. Rev. Accel. Beams **21**, 010101 (2018).
   → 螺线管/射频耦合与发射度补偿（AG 横向耦合修正的背景）。
   File: `AG/文献/Dowell_2018_Exact_Cancellation_Emittance_Coupled_Transverse.pdf`

## 空间电荷 / Space charge

4. **Stupakov** —
   "Space charge effects in an accelerated beam",
   Phys. Rev. ST Accel. Beams **11** (2008) [文章编号以 PDF 为准].
   → 加速束团空间电荷效应（AG/OCELOT SC 模块背景）。
   File: `AG/文献/Stupakov_2008_Space_Charge_Accelerated_Beam_PRSTAB11.pdf`

5. **Kim** —
   "RF and Space-Charge Effects in Laser-Driven RF Electron Guns",
   Nucl. Instrum. Methods Phys. Res. A **275**, 201–218 (1989).
   → RF 腔纵向/横向效应与空间电荷的经典处理（含中文翻译工程）。
   File: `AG/文献/RF-SC/kim1989_lecture.pdf`
   (中文翻译: `AG/文献/RF-SC/translation_project/`)

## UED 束流控制 / UED beam control

6. **Williams et al.** —
   "Active control of bright electron beams with RF optics for
   femtosecond electron diffraction",
   (2017).
   → UED RF 光学束流控制（实验应用背景）。
   File: `AG/文献/Williams et al. - 2017 - Active control of bright electron beams with RF op.pdf`

## 内部笔记（不属公开文献，随项目文档保留）

- `AG/文献/physics_principles.tex/.pdf` — 物理原理笔记
- `AG/文献/code_update_log.tex/.pdf` — 代码更新日志
- `AG/文献/main.tex` — 汇总文档
- `AG/文献/RF-SC/RF-SC.pdf` — RF 与空间电荷综述（来源与上述文献相关）
