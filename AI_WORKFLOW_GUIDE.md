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

# 任务结束时的审核输出规范

任务完成后，不要只输出“已完成”或文件列表。

必须在最终回复中输出一个完整的「给用户和 GPT 审核的摘要」，即使代码和报告已经写入文件，也必须把最关键的信息直接写在终端回复中。

## 一、最终回复固定结构

### 1. 一句话结论

用一句话说明：

- 本轮到底做了什么；
- 是否修改物理；
- 是否通过验收；
- 是否仍有开放问题。

示例：

> 本轮仅重构 GPT 主路由的 lattice 数据流，未修改任何物理公式；主路由与 validation 路由几何完全一致，四项回归测试全部通过。

---

### 2. 修改范围

列出：

- 修改的文件；
- 每个文件修改的作用；
- 明确说明哪些核心文件没有修改。

格式：

| 文件 | 修改内容 | 风险等级 |
|---|---|---|
| `xxx.py` | …… | 低/中/高 |

同时单独列出：

```text
确认未修改：
- AG 核心……
- OCELOT 核心……
- RF 方程……

3. 最关键的代码变化
不要粘贴完整文件。
只粘贴最关键的 2～5 段代码，每段控制在 10～40 行，并标注：
* 文件路径；
* 函数名；
* 修改前逻辑；
* 修改后逻辑；
* 为什么这样改。
优先展示容易产生数据流错误、单位错误或重复调用的代码。
必须包括：
1. 参数从哪里读取；
2. 物理量在哪里转换；
3. 元件在哪里组装；
4. 开关在哪里生效；
5. 输出量在哪里计算。
若代码较长，使用精简伪代码或 diff，不要省略关键乘除因子、单位和条件判断。

4. 数据流与调用链
用纯文本画出本轮实际调用链，例如：
beamline_config.yaml
        ↓
shared.params
        ↓
validation.backend.run_ocelot
        ↓
lattice builder
        ↓
OCELOT tracking
        ↓
BeamResult
        ↓
test_full_beamline
必须特别指出：
* 是否存在重复参数源；
* 是否存在旧路径仍被调用；
* 是否存在未使用的历史脚本；
* 是否存在某个测试绕过了正式入口。

5. 验收结果
用表格给出关键数值，不能只说 PASS。
至少包括：
* 修改前；
* 修改后；
* 参考值或另一后端；
* 相对偏差；
* 阈值；
* PASS/FAIL。
示例：
指标	修改前	修改后	参考	偏差	判定
若有随机宏粒子噪声，必须说明：
* 随机种子策略；
* 预期噪声范围；
* 本轮差异是否在噪声范围内。

6. 图像摘要
如果任务产生多张图，额外生成一张：
validation/reports/review_summary.png
要求：
* 将最关键的 3～6 幅结果放进一张合并图；
* 每幅图有清晰标题；
* 标出元件位置、关键差异和样品面；
* 分辨率适合手机查看；
* 不要只给本地文件路径。
最终回复中说明这张图包含什么，以及最值得人工检查的位置。
如果当前终端支持 Markdown 图片显示，输出：
![审核汇总图](validation/reports/review_summary.png)
如果终端不显示图片，则明确告诉用户：
最值得上传给 GPT 审核的图片是 validation/reports/review_summary.png。

7. 风险与开放问题
必须分成三类：
已解决
列出已经有测试和证据支持的问题。
尚未解决
列出还不能下结论的问题。
本轮可能引入的新风险
例如：
* 默认参数隐藏；
* 多实例遗漏；
* 元件重复调用；
* 坐标重复转换；
* 随机数相关性；
* 测试和正式入口走不同路径；
* 输出名称与内部变量含义不一致。
不要为了显得任务成功而省略风险。

8. 最值得发给 GPT Plus 审核的材料
按优先级列出最多 5 项：
A 级：强烈建议发送
只包括最容易出现严重错误、仅靠文字难以验证的文件或图片。
B 级：出现异常时发送
用于进一步排查。
C 级：通常不需要发送
仅作存档。
示例：
A级：
1. validation/backend.py
   原因：包含双后端接口、RF/R56 转换和元件组装，风险最高。
2. validation/reports/review_summary.png
   原因：一张图可同时判断曲线、腰位置和异常跳变。

B级：
3. shared/beamline_config.yaml
   原因：怀疑参数源或元件顺序时再发。

C级：
4. CHECKPOINTS.md
   原因：用于追溯历史，不是本轮代码正确性的直接证据。

9. 给 GPT 复制的最小审核文本
在最终回复最后，生成一段可以直接复制给 GPT 的文字，包含：
* 项目背景一句话；
* 本轮目标；
* 修改文件；
* 关键公式；
* 修改前后结果；
* 尚未解决问题；
* 希望 GPT 重点审核的三件事。
控制在 800～1500 字以内。
标题固定为：
===== 可直接复制给 GPT 审核 =====

二、停止条件
出现以下任意情况时，不要继续自动补丁：
* 原有通过测试变为失败；
* AG 确定性数组发生非预期变化；
* config SHA 意外改变；
* 元件数量或顺序不一致；
* RF kick 数量错误；
* 出现 NaN、负方差或发射度接近零；
* 为使两模型吻合而准备修改参数；
* 无法确认某变量的定义或单位。
应立即停止，并在最终回复中给出第一个失败点、相关文件和建议上传给 GPT 审核的材料。
---
