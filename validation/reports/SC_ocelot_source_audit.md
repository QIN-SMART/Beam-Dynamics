# SC — OCELOT Source Audit (installed 26.06.1)

Date: 2026-08-10 · 基于当前安装源码，非历史注释

## 1. SpaceCharge 类定义位置
`/opt/anaconda3/lib/python3.12/site-packages/ocelot/cpbd/sc.py:76` `class SpaceCharge(PhysProc)`

## 2. 构造参数 (sc.py:92-102)
```python
def __init__(self, step=1, **kwargs):
    self.step = step                     # in Navigator.unit_step
    self.nmesh_xyz = kwargs.get("nmesh_xyz", [63, 63, 63])
    self.low_order_kick = kwargs.get("low_order_kick", True)
    self.random_mesh = kwargs.get("random_mesh", False)
    self.random_seed = 10                # 固定! prepare() 中 np.random.seed(10)
```

## 3. step 含义 (sc.py:81, docstring + apply 调用)
"step of the Space Charge kick applying **in Navigator.unit_step**" —— 每 `step × unit_step` 米施加一次 kick（Navigator 每步调一次 physics process，sc.apply 内部按计数触发）。

## 4/5. mesh 参数名与默认值
`nmesh_xyz`（kwargs 键），默认 `[63, 63, 63]`（sc.py:95）。**没有独立的 "mesh" 键**——若传 `SpaceCharge(mesh=...)` 无效，只有 `nmesh_xyz`。

## 6. 粒子电荷来源
`p_array.q_array`（sc.py:182,193,241）。q_array 由 `generate_parray(charge=Q_total)` 设置，每粒子电荷 = Q_total/N。

## 7. 总电荷进入 Poisson 的方式
`el_field()` 内 `q = np.bincount(inds, Q, ...)`（sc.py:193）——直接把每粒子电荷沉积到网格；`potential()` 用 FFT 卷积 Green 函数后除以 `4πε₀·hx·hy·hz`（sc.py:167）。电荷单位：库仑（generate_parray charge 为 C）。

## 8. Deposition 方法
**NGP（最近网格点）**：`Xi = floor(X)+1`，`np.bincount(inds, Q)`（sc.py:191-193）。不是 CIC。

## 9. Poisson solver 类型
**积分表示 + FFT 卷积**（自由空间 Green 函数，ASTRA 同款，sc.py:85-90,135-168）。不是差分/谱求解器；无边界条件（自由空间）。

## 10. Field interpolation
`scipy.ndimage.map_coordinates(..., order=1)` —— **三线性插值**（sc.py:202-204）。

## 11. 纵向电场
**包含**：Ez 由势差分得到（sc.py:200），kick 更新 `xp[5] += cdT·Exyz[:,2]`（sc.py:248）。

## 12. Relativistic 变换
束团静止系解 Poisson：`X[:,2] *= gamma` 洛伦兹拉伸（sc.py:172）；动量→MAD 坐标（`xxstg_2_xp_mad`，sc.py:221）；束团平均速度系变换（sc.py:223-239）；回实验室系 kick：横向乘 `(1-β0²)` 因子（sc.py:246-247），纵向不乘（sc.py:248）。

## 13. Kick 更新方式
在 MAD 动量坐标上更新 `xp[3:6]`（= rparticles 的 px/py/p 行），再 `xp_2_xxstg_mad` 转回（sc.py:251）。即**直接改 rparticles[1,3,5]**——p_oc 与 px/py 被 SC 更新。

## 14. Navigator 调用频率
`apply(p_array, zstep)` 由 Navigator 的 physics-process 机制在每个 tracking step 调用；`step` 参数控制每隔几次执行（sc.py:81 语义）。

## ⚠️ 发现
- **random_seed=10 固定**（sc.py:102），`prepare()` 无条件 `np.random.seed(10)`（sc.py:106-107）——若 Navigator 在束生成前 prepare，会污染全局 RNG 状态（当前流程 prepare 在束生成后，影响有限，但需记录）。
- `low_order_kick=True` 默认——存在低阶 kick 分支（未深入）。
