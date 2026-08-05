
# UED Beamline Integration and Validation Task

## 项目背景

这是一个超快电子衍射(UED) beamline模拟项目。

目前已经完成：

1. Drift元件模型
2. Solenoid磁透镜模型
3. RF cavity纵向压缩模型
4. Space charge模块初步实现

这些模块已经经过独立验证。

当前任务不是重新开发物理模型。

当前任务：

> 将已有验证模块组装成为完整beamline，并保证不同模型之间参数一致、物理定义一致。

---

# 第一原则（非常重要）

禁止重新实现已经存在的模块。

禁止：

- 新写beam generation
- 新写solenoid公式
- 新写RF kick公式
- 新写emittance定义
- 新写space charge模型


必须：

先搜索已有代码：

```

find .
grep -r "Solenoid"
grep -r "RF"
grep -r "generate"
grep -r "emittance"

```

找到已经验证的实现。

然后：

import并调用。

不要复制代码。

---

# 当前已有验证结果

## 已验证模块

### Drift

验证目标：

自由漂移：

\[
x=x_0+x'_0z
\]

beam size:

\[
\sigma_x(z)
\]

符合解析结果。


---

### Solenoid

已验证：

输入：

B field

length

particle energy


输出：

beam focusing

beam waist

符合理论趋势。


注意：

solenoid内部不能直接使用普通projected emittance。

原因：

存在canonical angular momentum:

\[
P_x=p_x+eA_x
\]


因此：

solenoid内部emittance只作为参考。

---

### RF

目前需要重新检查。

不要认为RF模块正确。

必须单独benchmark。

验证：

输入：

(z)

输出：

(delta)


应该出现：

linear chirp:

\[
\delta=h z
\]


---

# 当前发现的问题

当前beamline结果：

## 问题1

SC OFF 和 SC ON几乎完全重合。

说明：

space charge可能没有真正作用。


需要检查：

1.
SpaceCharge是否加入Navigator

2.
physics process范围是否覆盖整个beamline

3.
step size是否合理


不要直接调参数。


---

## 问题2

RF后纵向phase space没有明显chirp。


检查：

RF kick定义。


注意OCELOT：

delta定义：

\[
\delta=\frac{p-p_0}{p_0}
\]


不是：

\[
\frac{\Delta E}{E}
\]


需要确认：

energy kick

转换

是否正确。


---

## 问题3

不要手动覆盖OCELOT longitudinal tracking。


检查是否存在类似：

```

tau += ...

```

或者：

```

z = ...

```

人为推进。


OCELOT已经负责：

longitudinal transport。


禁止重复。


---

# Beamline正确架构

应该：

```

Beam parameters
|
v
Generate ParticleArray
|
v
Navigator
|
|
+ Drift
|
+ Solenoid
|
+ Drift
|
+ RF cavity
|
+ Drift
|
v
Diagnostics

```


---

# 参数管理要求

所有参数必须来自：

beamline_config.yaml


包括：

beam:

- energy
- charge
- bunch length
- emittance


solenoid:

- B
- length
- position


RF:

- voltage
- frequency
- phase


禁止：

代码内部出现：

```

100e3
0.2
10e6

```


这种magic number。


---

# Debug流程（必须遵守）

不要一次修改多个地方。

按照：

## Step 1

只运行：

Cathode → Drift

比较：

AG模型

OCELOT


要求：

sigma_x误差 <5%


---

## Step 2

加入Solenoid


比较：

beam waist位置


---

## Step 3

加入RF

只看：

(z,delta)


不看beam size。


---

## Step 4

加入Space charge


比较：

SC OFF

SC ON


检查：

emittance增长。


---

# 每次修改必须输出

1.
修改原因

2.
修改文件

3.
修改前后物理意义

4.
验证结果


禁止：

直接大规模重构。


---

# 代码质量要求

如果发现已有函数：

例如：

```

apply_rf_kick()
emit()
generate_beam()

```

必须先判断：

是否已经存在。


如果存在：

修改已有函数。


不要复制。


---

# 最终目标

实现两个模型：

AG Gaussian model

和

OCELOT particle tracking model


共享：

beam_parameters.yaml


beamline.yaml


输出：

sigma_x(z)

sigma_z(z)

epsilon(z)

phase space


用于论文结果交叉验证。


```



