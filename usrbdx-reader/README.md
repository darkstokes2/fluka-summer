# USRBDX Spectrum Reader

解析 FLUKA USRBDX `tab.lis` 输出文件，生成微分能谱图、误差棒、CSV 导出和能量积分。

## 功能

- 解析 `tab.lis` 四列表格（Emin, Emax, 微分值, 相对误差%）
- 支持多文件叠加绘图，对比不同工况
- 连续折线或阶梯图（推荐阶梯图用于分箱数据）
- bin 中心标记、误差棒、统计不足数据屏蔽
- 线性/对数坐标、自动坐标范围
- 电流缩放、人工 Y 轴缩放
- 能量积分或求和（支持部分 bin 重叠计算）
- 输出 PNG、PDF、SVG、曲线 CSV 和积分结果 CSV

## 依赖

```
matplotlib
numpy
```

## 用法

编辑脚本顶部的配置区，然后运行：

```bash
python usrbdx_spectrum_reader.py
```

### 关键配置

```python
# 输入文件列表
CASES = [
    ("spectrum_tab.lis", r"$\mathrm{Case\ A}$", 1.0, 1.0),
    # (文件名, 图例标签, 束流电流 mA, 人工 Y 缩放倍率)
]

# Detector 选择
DETECTOR_NUMBER = 1
DETECTOR_NAME_CONTAINS = None

# 电流缩放
APPLY_CURRENT_SCALING = False
REFERENCE_CURRENT_MA = 1.0

# 曲线形式
USE_CONTINUOUS_LINE = True  # True=折线, False=阶梯图

# 误差棒
SHOW_ERROR_BARS = True

# 坐标轴
USE_LOG_X = True
USE_LOG_Y = True

# 能量积分
CALCULATE_TOTAL = True
TOTAL_MODE = "integral"  # "integral" 或 "sum"
```

### 输出文件

| 文件 | 内容 |
|------|------|
| `{stem}.png/pdf/svg` | 能谱图（600 DPI） |
| `{stem}.csv` | 能量 bin、原始/缩放能谱、误差 |
| `{stem}_totals.csv` | 每个工况的积分/求和结果 |

## 能量 bin 中心计算

- 对数 bin（所有 Emin/Emax > 0）：几何平均 `sqrt(Emin * Emax)`
- 线性 bin：算术平均 `(Emin + Emax) / 2`

## 对数 Y 轴误差棒

对数 Y 轴下，下误差棒被截断以防止穿过零点（仅影响显示，不影响 CSV 导出的误差值）。

## 注意事项

- 只读取每个 Detector 的第一张四列表格（积分立体角能谱），自动跳过双重微分分布
- 多文件绘图时要求所有文件使用相同的能量 bin 网格
- 阶梯图模式要求相邻 bin 连续（`Emax[i] == Emin[i+1]`）
- 脚本不额外除以或乘以面积，直接使用 `tab.lis` 第三列值
