# USRBIN Post-Processor

读取 FLUKA USRBIN `usbrea` `.lis` 输出文件，生成纵向剖面图和二维截面热图。

## 功能

- 支持一个 `.lis` 文件中包含多个 USRBIN Detector
- 按 Detector 编号和名称关键词选择目标
- 读取三维评分矩阵和百分比误差矩阵
- 将三维矩阵沿 X-Y 横截面约化为 Z 向一维分布（平均 / 面积积分 / 求和）
- 生成 X-Y、X-Z、Y-Z 二维截面热图，支持多文件合并
- 支持电流缩放、人工缩放、阶梯图/折线、误差棒
- 统计不足数据屏蔽、线性/对数坐标
- 指定 Z 范围积分或求和
- 输出 PNG、PDF、SVG、曲线 CSV 和积分结果 CSV
- 可选叠加材料区域背景（Z 轴条带或任意区域函数）

## 依赖

```
matplotlib
numpy
```

可选：`fluka-geometry` 模块（用于在热图上叠加几何边界）。

## 用法

编辑脚本顶部的配置区，然后运行：

```bash
python usrbin_postprocessor.py
```

### 关键配置

```python
# 输入文件列表
CASES = [
    ("output.lis", r"$\mathrm{Case\ A}$", 1.0，1.0，),
    # (文件名, 图例标签, 束流电流 mA，人工缩放倍率)
]

# Detector 选择
DETECTOR_NUMBER = 1
DETECTOR_NAME_CONTAINS = "Edeposit"

# 横截面约化模式："average" / "area_integral" / "sum"
TRANSVERSE_REDUCTION_MODE = "area_integral"

# 电流缩放
APPLY_CURRENT_SCALING = False
REFERENCE_CURRENT_MA = 1.0

# 误差棒
SHOW_ERROR_BARS = False

# 二维热图
GENERATE_XY_HEATMAPS = False
HEATMAP_Z_POSITIONS_CM = [0.0]
```

### 输出文件

| 文件 | 内容 |
|------|------|
| `{stem}.png/pdf/svg` | 纵向剖面图 |
| `{stem}.csv` | Z-bin 数据（原始/缩放/误差） |
| `{stem}_totals.csv` | 积分/求和结果 |
| `{stem}_xy_heatmaps/` | 二维热图目录 |

## 与 fluka-geometry 配合

本脚本会尝试 import `fluka_geometry` 模块。如果可用，可在热图上叠加几何边界曲线。将 `fluka_geometry.py` 放在同目录或加入 PYTHONPATH 即可。

## 材料区域背景

支持两种模式：

- `MATERIAL_REGION_MODE = "z_bands"`：简单的 Z 轴条带
- `MATERIAL_REGION_MODE = "region_function"`：任意空间谓词函数（径向、角向、嵌套壳等）
