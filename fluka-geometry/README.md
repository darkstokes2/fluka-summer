# FLUKA Geometry Parser

解析 FLUKA `.inp` 几何定义文件（GEOBEGIN/GEOEND），评估 CSG 区域表达式，计算任意切面上的材料分布。

## 功能

- 解析 body 卡片：SPH、ZCC、PLA、XYP、XZP、YZP
- 解析 region CSG 表达式（`+body` / `-body` / `|` / `&` / `()`）
- 解析 ASSIGNMA 材料到区域的映射
- 计算切面（X/Y/Z = const）上的材料分布网格
- 提取解析边界曲线，可叠加到 USRBIN 热图上

## 依赖

```
numpy
```

## 用法

### 解析 .inp 文件

```python
from fluka_geometry import parse_inp

geom = parse_inp("simulation.inp")
# geom.bodies: {name: Body}
# geom.regions: list[Region]
# geom.material_to_regions: {material: [region_names]}
```

### 计算切面材料分布

```python
from fluka_geometry import cut_plane
import numpy as np

free_edges_1 = np.linspace(-5.0, 5.0, 101)  # x 边界
free_edges_2 = np.linspace(-5.0, 5.0, 101)  # y 边界

mgrid = cut_plane(geom, axis="z", value=0.0,
                  free_edges_1=free_edges_1,
                  free_edges_2=free_edges_2)

# 查询材料掩码
mask = mgrid.material_mask("COPPER")

# 导出 CSV
mgrid.to_csv("material_distribution.csv")
```

### 提取边界曲线

```python
from fluka_geometry import collect_analytic_boundaries

curves = collect_analytic_boundaries(geom, axis="z", value=0.0)

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
for c in curves:
    ax.plot(c["free1"], c["free2"], "k-", linewidth=0.5)
ax.set_aspect("equal")
```

## 支持的 Body 类型

| 类型 | 卡片格式 | 说明 |
|------|---------|------|
| SPH | `SPH name cx cy cz r` | 球 |
| ZCC | `ZCC name cx cy r` | Z 轴无限圆柱 |
| XYP | `XYP name z0` | 平面 z = z0，法向 +z |
| XZP | `XZP name y0` | 平面 y = y0，法向 +y |
| YZP | `YZP name x0` | 平面 x = x0，法向 +x |
| PLA | `PLA name Vx Vy Vz X1 Y1 Z1` | 任意平面 |

**符号约定：** `+body` = 内部（法向量反方向半空间）；`-body` = 外部（法向量所指方向半空间）。

## 扩展

在 `BODY_PARSERS` 字典中注册新的 body 类型即可扩展。预注册但未实现的类型槽位：XCC、YCC、RCC、BOX、REC、TRC、CCC、ARB、WED、RAW。

## 独立使用

本模块不依赖任何其他模块，可独立 import。也被 `usrbin-reader` 可选导入，用于在热图上叠加材料边界。
