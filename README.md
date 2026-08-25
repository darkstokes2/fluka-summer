# fluka-summer

一套 FLUKA 蒙特卡洛模拟后处理工具集，涵盖几何解析、USRBIN/USRBDX 批处理合并、数据读取与可视化绘图。

## 功能概览

| 模块 | 目录 | 功能 |
|------|------|------|
| FLUKA 几何解析器 | [`fluka-geometry/`](fluka-geometry/) | 解析 `.inp` 几何定义（body/region CSG 表达式），计算切面材料分布，提取解析边界曲线用于热图叠加 |
| USRBIN 读图程序 | [`usrbin-reader/`](usrbin-reader/) | 读取 `usbrea` `.lis` 输出，生成 Z 向一维分布剖面图、X-Y/X-Z/Y-Z 二维截面热图，支持误差棒、电流缩放、积分 |
| USRBIN 批处理 | [`usrbin-merge/`](usrbin-merge/) | 用 `usbsuw` 合并 `_fort.NN` 文件，`usbrea` 转 `.lis`，支持自动发现文件夹和命令行参数 |
| USRBDX 批处理 | [`usrbdx-merge/`](usrbdx-merge/) | 用 `usxsuw` 合并 `_fort.NN` 文件，自动生成 `.bnn` + `_sum.lis` + `_tab.lis` |
| USRBDX 能谱读取 | [`usrbdx-reader/`](usrbdx-reader/) | 解析 `tab.lis` 四列表格，绘制微分能谱图（阶梯图/折线）、误差棒、积分/求和、CSV 导出 |

## 后处理流水线

```
_fort.NN 文件 ──→ [usrbin-merge / usrbdx-merge] ──→ .lis / tab.lis 文件
                                                          │
              [fluka-geometry] ←── .inp 几何文件 ────────┤
                    │                                     │
                    └── 材料边界叠加 ──→ [usrbin-reader / usrbdx-reader] ──→ 图表 + CSV
```

- `fluka-geometry` 和 `usrbin-reader` 配套使用：geometry 提供的材料边界可直接叠加到 USRBIN 热图上
- `fluka-geometry` 也可独立使用，解析输入卡的几何和材料映射
- 批处理脚本（merge）是流水线第一步，生成 `.lis` 供 reader 程序读取

## 快速开始

### 环境要求

- Python 3.8+
- FLUKA 工具链（`usbsuw`、`usbrea`、`usxsuw`）需要在 PATH 中（仅批处理脚本需要）

### 安装依赖

```bash
pip install numpy matplotlib pexpect
```

- `numpy`：所有模块
- `matplotlib`：reader 模块（usrbin-reader、usrbdx-reader）
- `pexpect`：merge 模块（usrbin-merge、usrbdx-merge）

### 使用示例

```bash
# 1. 批处理合并 USRBIN fort 文件
cd usrbin-merge
python fluka_batch_processor.py --base /path/to/simulation --fort-unit 42

# 2. 读取并绘图 USRBIN 结果
cd usrbin-reader
# 编辑 usrbin_postprocessor.py 顶部的 CASES 配置，指定 .lis 文件
python usrbin_postprocessor.py

# 3. 批处理合并 USRBDX fort 文件
cd usrbdx-merge
python usrbdx_batch_processor.py --base /path/to/simulation --fort-units 50

# 4. 读取并绘图 USRBDX 能谱
cd usrbdx-reader
# 编辑 usrbdx_spectrum_reader.py 顶部的 CASES 配置，指定 tab.lis 文件
python usrbdx_spectrum_reader.py

# 5. 解析几何（独立使用或配合 usrbin-reader）
cd fluka-geometry
python -c "from fluka_geometry import parse_inp; g = parse_inp('simulation.inp'); print(len(g.bodies), 'bodies')"
```

各模块的详细用法见对应目录下的 README。

## 目录结构

```
fluka-summer/
├── README.md              ← 你在这里
├── LICENSE                ← GPL-3.0
├── .gitignore
├── fluka-geometry/        ← 几何解析器
│   ├── README.md
│   └── fluka_geometry.py
├── usrbin-reader/         ← USRBIN 读图程序
│   ├── README.md
│   └── usrbin_postprocessor.py
├── usrbin-merge/          ← USRBIN 批处理
│   ├── README.md
│   └── fluka_batch_processor.py
├── usrbdx-merge/          ← USRBDX 批处理
│   ├── README.md
│   └── usrbdx_batch_processor.py
└── usrbdx-reader/         ← USRBDX 能谱读取
    ├── README.md
    └── usrbdx_spectrum_reader.py
```

## 许可证

本项目采用 [GPL-3.0](LICENSE) 许可证。

这意味着你可以自由使用、修改和分发本项目，但**衍生作品必须同样以 GPL-3.0 开源**，且必须保留原始版权声明。

## 引用与致谢

如果本项目对你的研究有帮助，欢迎引用：

```
fluka-summer — FLUKA post-processing toolkit
https://github.com/darkstokes2/fluka-summer
```

欢迎在 GitHub 提 Issue 和 PR。

## 联系方式

- zhangpb22@mails.tsinghua.edu.cn
- darkstokes@163.com

## 致谢

This toolkit was developed with assistance from AI coding tools (Trae). All physics logic, feature design, and engineering decisions were made by the author. The code is released under GPL-3.0.
