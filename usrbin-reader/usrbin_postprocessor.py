#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# fluka-summer — FLUKA post-processing toolkit
# Copyright (C) 2026  Pengbin Zhang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
FLUKA USRBIN usbrea 格式化输出读取、处理和绘图程序。

主要功能：
1. 支持一个 usbrea .lis 文件中包含多个 USRBIN Detector；
2. 按 Detector 编号和名称关键词选择目标 Detector；
3. 同时读取三维评分矩阵和对应百分比误差矩阵；
4. 将三维矩阵沿 X-Y 横截面约化为 Z 向一维分布；
5. 保留横截面平均（通量/评分密度）和横截面积分
   （径迹长度密度 dL/dz、能量沉积 dE/dz 等）功能；
6. 支持电流缩放、人工缩放、阶梯图、折线、标记、误差棒；
7. 支持统计不足数据屏蔽、线性/对数坐标和自动坐标范围；
8. 支持指定 Z 范围积分或求和；
9. 输出 PNG、PDF、SVG、曲线 CSV 和积分结果 CSV；
10. 保留原程序中的材料区域背景和薄层引线标注功能；
11. 支持按实际 z 坐标选择 Z-bin，并输出 X-Y 二维截面热图及截面 CSV；
12. 支持多维度截面热图（X-Y、X-Z、Y-Z），可按维度组合选择；
13. 支持多文件热力图合并绘制，以子图形式排列于同一张图中。

说明：
- usbrea 输出的三维数组按 Fortran 顺序 A(ix, iy, iz) 展开，
  即 ix 变化最快、iy 次之、iz 最慢。本程序使用 order="F" 还原。
- 将不同 X-Y bin 合并时，百分比误差按各 bin 独立的假设进行
  方差传播：sigma_profile = sqrt(sum(sigma_bin^2)) × reduction_factor。
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize
from matplotlib.transforms import blended_transform_factory

# FLUKA .inp 几何解析器（同目录模块）
try:
    from fluka_geometry import (
        parse_inp as _fluka_parse_inp,
        collect_analytic_boundaries as _fluka_collect_boundaries,
        _collect_body_refs as _fluka_collect_body_refs,
    )
    _FLUKA_GEOMETRY_AVAILABLE = True
except ImportError:
    _FLUKA_GEOMETRY_AVAILABLE = False


# ============================================================
# 1. 输入文件
# ============================================================
BASE_DIR = Path(__file__).resolve().parent

# CASES 每项格式：
#
# (
#     文件名，
#     图例名称，
#     束流电流/mA，
#     用户人工 Y 轴缩放倍率，
# )
#
# 示例：
#
# (
#     "your_energy_deposition.lis",
#     r"$\mathrm{Case\ A}$",
#     1.0,
#     1.0，
# ),
#
# 请填入你自己的 *.lis 文件名。
CASES = [
    # (
    #     "your_energy_deposition.lis",
    #     r"$\mathrm{Case\ A}$",
    #     1.0,
    #     1.0,
    # ),
]


# ============================================================
# 2. 用户配置区
# ============================================================

# ------------------------------------------------------------
# A. USRBIN Detector 选择
# ------------------------------------------------------------
# 一个 usbrea .lis 文件可能包含多个 USRBIN Detector。
#
# 按 Detector 编号选择；不按编号筛选时设为 None。
DETECTOR_NUMBER = None

# 同时按 Detector 名称中的关键词筛选。
# 例如：DETECTOR_NAME_CONTAINS = "Edeposit"
# 不需要名称筛选时设为 None。
DETECTOR_NAME_CONTAINS = None

# 如果百分比误差矩阵缺失：
# True  -> 立即报错；
# False -> 允许继续，误差数组以 NaN 填充，不绘制误差棒。
REQUIRE_ERROR_MATRIX = False


# ------------------------------------------------------------
# B. X-Y 横截面约化方式
# ------------------------------------------------------------
# 可选：
#
# "average"：
#     profile(z) = mean_xy(score)
#
#     保留原 USRBIN score 的量纲。
#     常用于横截面平均通量、平均能量密度、平均 DPA 等。
#
# "area_integral"：
#     profile(z) = sum_xy(score) * dx * dy
#
#     对横截面积分。
#     若原 score 为 track-length fluence（cm^-2 primary^-1），
#     则结果可解释为 dL/dz；再沿 z 积分得到总径迹长度。
#     若原 score 为能量沉积密度（GeV cm^-3 primary^-1），
#     则结果可解释为 dE_dep/dz。
#
# "sum"：
#     profile(z) = sum_xy(score)
#
#     仅求和，不乘 dx*dy。只有确认物理量本身已是每 bin 总量时使用。
TRANSVERSE_REDUCTION_MODE = "area_integral"

# 兼容原代码中的开关：
# None  -> 使用上面的 TRANSVERSE_REDUCTION_MODE；
# True  -> 强制使用 "average"；
# False -> 强制使用 "area_integral"。
#
# 因此原来的 average_x_and_y=True/False 功能仍然保留。
average_x_and_y = None

# True：要求所有工况的 X、Y、Z 三维网格均与第一个文件相同。
# False：只要求 Z 网格一致；不同 X-Y 网格仍可生成各自的约化曲线。
REQUIRE_IDENTICAL_XY_GRID = True


# ------------------------------------------------------------
# C. 束流电流缩放
# ------------------------------------------------------------
# False：保留 FLUKA 原始 per-primary 结果。
# True ：按 current_mA / REFERENCE_CURRENT_MA 缩放。
APPLY_CURRENT_SCALING = False
REFERENCE_CURRENT_MA = 1.0


# ------------------------------------------------------------
# D. 每根曲线的人工缩放倍率
# ------------------------------------------------------------
# 在 CASES 每项的最后一个数中设置：
# 1.0 不缩放；2.0 放大 2 倍；0.2 缩小至 20%。
SHOW_USER_SCALE_IN_LEGEND = True


# ------------------------------------------------------------
# E. 曲线形式
# ------------------------------------------------------------
# True ：连接 Z-bin 中心，绘制普通折线。
# False：按 Z-bin 边界绘制阶梯图，更符合分箱结果。
USE_CONTINUOUS_LINE = False


# ------------------------------------------------------------
# F. bin 中心标记
# ------------------------------------------------------------
SHOW_BIN_MARKERS = False
BIN_MARKER = "o"
BIN_MARKER_SIZE = 3.0
MARK_EVERY = 1


# ------------------------------------------------------------
# G. 阶梯图辅助连接线
# ------------------------------------------------------------
ADD_GUIDE_LINE = False
GUIDE_LINE_WIDTH = 0.75
GUIDE_LINE_ALPHA = 0.55


# ------------------------------------------------------------
# H. 统计误差棒
# ------------------------------------------------------------
# USRBIN 的每个三维 bin 都有相对误差。
# 合并 X-Y bin 后，本程序按独立误差进行方差传播。
SHOW_ERROR_BARS = False
ERROR_BAR_LINEWIDTH = 0.7
ERROR_BAR_CAPSIZE = 2.0
ERROR_BAR_ALPHA = 0.60
ERROR_BAR_EVERY = 1


# ------------------------------------------------------------
# I. 统计不足数据处理
# ------------------------------------------------------------
INSUFFICIENT_ERROR_PERCENT = 99.0

# True：隐藏约化后一维误差达到或超过 99% 的非零数据点。
MASK_INSUFFICIENT_STATISTICS = False

# True：误差棒不显示误差达到或超过 99% 的点。
SKIP_INSUFFICIENT_ERROR_BARS = True


# ------------------------------------------------------------
# J. X 轴（Z 坐标）
# ------------------------------------------------------------
X_AXIS_LABEL = r"$z$ (cm)"
USE_LOG_X = False

# None 表示使用完整 Z 范围。
PLOT_Z_MIN = -5.0
PLOT_Z_MAX = 10.0


# ------------------------------------------------------------
# K. Y 轴标题和范围
# ------------------------------------------------------------
# 当前示例：横截面积分后的能量沉积分布 dE/dz。
Y_AXIS_LABEL = (
    r"Energy deposition per unit $z$, "
    r"$\mathrm{d}E_{\mathrm{dep}}/\mathrm{d}z$ "
    r"(GeV cm$^{-1}$ primary$^{-1}$)"
)

# 其他示例：
#
# 横截面平均通量：
# Y_AXIS_LABEL = (
#     r"Cross-section-averaged fluence "
#     r"(cm$^{-2}$ primary$^{-1}$)"
# )
#
# 径迹长度密度：
# Y_AXIS_LABEL = (
#     r"Track length per unit $z$, "
#     r"$\mathrm{d}L/\mathrm{d}z$ "
#     r"(cm cm$^{-1}$ primary$^{-1}$)"
# )
#
# 横截面平均 DPA：
# Y_AXIS_LABEL = r"Cross-section-averaged DPA (primary$^{-1}$)"

USE_LOG_Y = False
Y_HEADROOM_FACTOR = 1.18
LOG_Y_BOTTOM_FACTOR = 0.7
USER_Y_MIN = None
USER_Y_MAX = None


# ------------------------------------------------------------
# L. Z 方向积分/求和
# ------------------------------------------------------------
CALCULATE_TOTAL = True

# "integral"：sum(profile_i * 与积分范围重叠的 dz_i)
# "sum"     ：直接对与范围相交的 Z-bin 求和。
TOTAL_MODE = "integral"

# True：对缩放后、图中显示的数据计算。
# False：对 FLUKA 原始约化 profile 计算。
TOTAL_USE_SCALED_PROFILE = True

TOTAL_NAME = "Total deposited energy"
TOTAL_UNIT = "GeV/primary"

# None 表示完整 Z 范围。
INTEGRAL_Z_MIN = None
INTEGRAL_Z_MAX = None


# ------------------------------------------------------------
# M. CSV 输出
# ------------------------------------------------------------
EXPORT_RAW_AND_SCALED = True
EXPORT_RELATIVE_ERROR = True
EXPORT_ABSOLUTE_ERROR = True


# ------------------------------------------------------------
# N. 图标题
# ------------------------------------------------------------
PLOT_TITLE = None


# ------------------------------------------------------------
# O. 输出文件名称
# ------------------------------------------------------------
OUTPUT_STEM = BASE_DIR / "usrbin_longitudinal_profile"


# ------------------------------------------------------------
# P. 材料区域背景（可选功能，默认关闭）
# ------------------------------------------------------------
# 开启后会在纵向曲线底部绘制材料区域彩色背景。
SHOW_MATERIAL_REGIONS = False
USE_THIN_REGION_CALLOUTS = True

# 每项格式：(左边界 cm, 右边界 cm, 名称, 颜色)
# 空列表表示不标注任何材料区域；请根据你的 FLUKA 几何填入。
MATERIAL_REGIONS: list[tuple[float, float, str, str]] = []

# 材料边界位置列表（cm），用于在曲线上绘制竖直虚线。
# 空元组表示不绘制边界线；请根据你的几何填入。
MATERIAL_BOUNDARIES: tuple[float, ...] = ()


# ------------------------------------------------------------
# Q. 二维截面热图
# ------------------------------------------------------------
# 总开关，False 时程序行为与原版本一致，不生成热图。
GENERATE_XY_HEATMAPS = True

# 需要生成热图的 CASE 编号，按 CASES 中从 1 开始的序号填写。
# 例如 [1, 6] 表示第 1 和第 6 个工况；None 表示全部工况。
HEATMAP_CASE_INDICES = None

# 需要查看的实际 z 坐标，单位 cm。可同时填写多个位置。
HEATMAP_Z_POSITIONS_CM = [1.0]  # 示例值，请根据你的几何填入实际 z 坐标

# z 位置选择方式：
# "containing_bin"：选择包含该 z 坐标的 Z-bin；若正好位于内部边界，选择右侧 bin。
# "nearest_center"：选择中心坐标距离给定 z 最近的 Z-bin。
HEATMAP_Z_SELECTION_MODE = "containing_bin"

# True：热图使用与纵向曲线相同的电流缩放和人工倍率。
# False：热图使用 FLUKA 原始三维 score。
HEATMAP_USE_SCALED_SCORE = True

# True：屏蔽百分比误差达到或超过 INSUFFICIENT_ERROR_PERCENT 的像素。
HEATMAP_MASK_INSUFFICIENT_STATISTICS = False

# 色标设置。
HEATMAP_USE_LOG_COLOR = True
HEATMAP_COLORMAP = "viridis"
HEATMAP_COLORBAR_LABEL = "Energy Score"

# None 表示根据每张热图的数据自动确定。
# 若希望多个工况使用完全相同的颜色范围，可手动设置这两个值。
HEATMAP_VMIN = None
HEATMAP_VMAX = None

# True：显示 X-Y bin 边界；网格很密时建议设为 False。
HEATMAP_SHOW_BIN_EDGES = False

# True：X、Y 坐标采用相同比例，保证几何形状不失真。
# 仅对 X-Y 横截面（俯视图）有物理意义；Y-Z/X-Z 纵向剖面
# 不同物理方向长度差异巨大，等比例会把细网格方向压扁，
# 建议关闭以让 Z 方向自然拉伸、细网格得以呈现。
HEATMAP_EQUAL_ASPECT = True
HEATMAP_FIGSIZE = (8.2, 6.0)
HEATMAP_DPI = 2400

# 输出格式只能从 png、pdf、svg 中选择。
HEATMAP_OUTPUT_FORMATS = ("png", "pdf", "svg")

# True：同时输出长表格式 CSV，包含每个 X-Y bin 的坐标、值和误差。
HEATMAP_EXPORT_CSV = True

# None：自动生成标题；字符串：所有热图统一使用该标题。
HEATMAP_TITLE = None

# 热图输出到独立目录，避免与纵向曲线文件混在一起。
HEATMAP_OUTPUT_DIR = Path(
    str(OUTPUT_STEM) + "_xy_heatmaps"
)


# ------------------------------------------------------------
# Q2. 热图维度选择
# ------------------------------------------------------------
# 需要绘制的截面维度组合，可多选：
#
# "xy"：固定 Z 坐标，绘制 X-Y 截面（原功能）；
# "xz"：固定 Y 坐标，绘制 X-Z 截面；
# "yz"：固定 X 坐标，绘制 Y-Z 截面。
#
# 例如 ("xy", "xz") 同时绘制 X-Y 和 X-Z 截面。
HEATMAP_DIMENSIONS = ("xy",)

# "xz" 维度所需的固定 Y 坐标列表（cm）。
# 用法与 HEATMAP_Z_POSITIONS_CM 相同。
HEATMAP_Y_POSITIONS_CM = [0.0]

# "yz" 维度所需的固定 X 坐标列表（cm）。
HEATMAP_X_POSITIONS_CM = [0.0]


# ------------------------------------------------------------
# Q3. 多文件热力图合并
# ------------------------------------------------------------
# True：将多个 CASE 的同一维度、同一固定坐标的热图
#       放在同一张图中，以子图（subplot）形式排列。
# False：每个 CASE 独立生成热图（原功能）。
HEATMAP_MERGE_FILES = True

# 合并模式下是否对不同文件使用不同视觉标识。
# True：不同文件使用不同色图和透明度，便于区分。
# False：所有文件使用统一色图（一般用法）。
HEATMAP_MERGE_DIFFERENTIATE = False

# 合并模式下不同文件使用的色图列表，按 CASE 顺序循环。
HEATMAP_MERGE_COLORMAPS = (
    "viridis",
    "plasma",
    "inferno",
    "magma",
    "cividis",
)

# 合并模式下不同文件使用的透明度列表，按 CASE 顺序循环。
HEATMAP_MERGE_ALPHAS = (
    1.0,
    0.6,
    0.6,
    0.6,
)

# 合并图排列列数。
HEATMAP_MERGE_NCOLS = 3


# ------------------------------------------------------------
# Q4. FLUKA 几何边界叠加（2D 热图专用）
# ------------------------------------------------------------
# 从 .inp 文件解析几何体，在 2D 热图上叠加材料区域边界线。
# 仅画边界线，不填色，不影响热图自身 colormap。
# 1D 纵向 profile 不受影响（保留原有 Z 条带）。
# 需要 fluka_geometry.py 模块（同目录）。

# 总开关。False 时完全不解析几何，热图与原版一致。
FLUKA_GEOMETRY_ENABLED = True

# .inp 文件路径（相对于本脚本所在目录，或绝对路径）。
# 所有 CASE 共享同一几何。
FLUKA_GEOMETRY_INP = ""  # 填入你的 .inp 文件名（相对于脚本目录或绝对路径）

# 截面位置模式：
#   None  = 自动跟随热图切片（与 HEATMAP_*_POSITIONS_CM 一致）
#   dict  = 独立指定，如 {"x": 0.0, "y": 0.0, "z": 0.05}
#           仅指定的轴生效，未指定的轴仍自动跟随
FLUKA_GEOMETRY_CUT_PLANE = None

# 边界线样式（单一颜色，不做材料色映射）。
FLUKA_GEOMETRY_LINE_COLOR = "#202020"
FLUKA_GEOMETRY_LINE_WIDTH = 0.3
FLUKA_GEOMETRY_LINE_STYLE = "-"

# 跳过这些材料的区域边界（避免真空/黑体边界干扰热图）。
FLUKA_GEOMETRY_SKIP_MATERIALS = {"VACUUM", "BLCKHOLE"}

# True：在热图目录下额外导出材料分布 CSV（长格式）。
FLUKA_GEOMETRY_EXPORT_CSV = True


# ------------------------------------------------------------
# Q5. 局部放大 inset（2D 热图专用）
# ------------------------------------------------------------
# 在大图角落画一个放大子图，并在大图上用虚线框标注放大区域。
# 仅画线，不改变热图数据或 colormap。

# 总开关。False 时无 inset。
HEATMAP_INSET_ENABLED = False

# 放大区域（物理坐标，cm）。
# 横轴范围（对应热图 h 轴：yz->Z, xz->X, xy->X）
HEATMAP_INSET_H_RANGE = (-1.0, 1.0)
# 纵轴范围（对应热图 v 轴：yz->Y, xz->Z, xy->Y）
HEATMAP_INSET_V_RANGE = (-1.0, 1.0)

# inset 在大图中的位置（相对坐标 0-1）。
# 格式 (left, bottom, width, height)。
# 相对于每个子图（不是整个 figure）。
HEATMAP_INSET_LOC = (0.13, 0.62, 0.15, 0.15)

# merged 模式下是否每个子图都画 inset。
# True: 每个子图都画（推荐）
# False: 只在第一个子图画
HEATMAP_INSET_PER_SUBPLOT = True

# 大图上虚线框样式
HEATMAP_INSET_BBOX_COLOR = "#E00000"
HEATMAP_INSET_BBOX_LINESTYLE = "--"
HEATMAP_INSET_BBOX_LINEWIDTH = 1.0

# 从虚线框到 inset 的连接线（指示放大关系）
# 单线模式：只从虚线框右上角到 inset 左下角
HEATMAP_INSET_CONNECT_LINES = True
HEATMAP_INSET_CONNECT_COLOR = "#888888"
HEATMAP_INSET_CONNECT_LINESTYLE = ":"
HEATMAP_INSET_CONNECT_LINEWIDTH = 0.6

# inset 边框样式
HEATMAP_INSET_EDGE_COLOR = "#666666"
HEATMAP_INSET_EDGE_LINEWIDTH = 0.8

# inset 标题（None 则自动生成 "Zoom: h=[..], v=[..]"）
HEATMAP_INSET_TITLE = None

# inset 是否使用与大图相同的 colormap/norm（True=一致，False=自动适配放大区数据）
HEATMAP_INSET_SHARE_NORM = True

# inset 是否也叠加 FLUKA 几何边界线
HEATMAP_INSET_SHOW_GEOMETRY = True


# ------------------------------------------------------------
# Q6. 3D 能量沉积导出（CSV + npz）
# ------------------------------------------------------------
# 总开关。False 时不导出 3D 数据。
# 输出 score_3d 全网格 (nx, ny, nz) 到专用文件夹，便于老师/合作者二次分析。
EXPORT_3D_SCORE_ENABLED = True

# 专用输出目录（脚本所在目录的子目录，带空格的显示名）
EXPORT_3D_SCORE_DIR = BASE_DIR / "figure of energy deposition in 3D"

# 是否输出 CSV（含坐标/raw/scaled/error 列）
EXPORT_3D_SCORE_CSV = True

# 是否输出 npz（推荐，二进制紧凑；float32 ~65MB，float64 ~130MB）
EXPORT_3D_SCORE_NPZ = True

# npz 精度：True=float32(省空间)，False=float64(全精度)
EXPORT_3D_SCORE_NPZ_FLOAT32 = True

# CSV 是否同时输出 raw + scaled 两列（与 2D CSV 一致；False 则只输出 raw）
EXPORT_3D_SCORE_CSV_SCALED = True

# CSV 是否含百分比误差列（无误差矩阵时填 NaN）
EXPORT_3D_SCORE_CSV_ERROR = True


# ============================================================
# 3. 正则表达式与基础工具
# ============================================================
FLOAT_RE = re.compile(
    r"[-+]?(?:\d+\.\d*|\.\d+|\d+)"
    r"(?:[EeDd][-+]?\d+)?"
)

# 示例：
# Cartesian binning n.   1  "Edeposit  " , generalized particle n.  208
BINNING_HEADER_RE = re.compile(
    r"^\s*(Cartesian|Cylindrical)\s+binning\s+n\.\s*"
    r"(\d+)\s*(?:\"([^\"]*)\")?\s*,?\s*(.*?)\s*$",
    re.IGNORECASE | re.MULTILINE,
)

AXIS_RE = re.compile(
    r"([XYZ])\s+coordinate:\s*"
    r"from\s*([-+0-9.EeDd]+)\s+to\s*([-+0-9.EeDd]+)\s+cm,\s*"
    r"(\d+)\s+bins\s*"
    r"\(\s*([-+0-9.EeDd]+)\s+cm\s+wide\s*\)",
    re.IGNORECASE,
)

PARTICLE_NUMBER_RE = re.compile(
    r"(?:generalized\s+particle|particle)\s+n\.\s*(\d+)",
    re.IGNORECASE,
)

DATA_MARKER_RE = re.compile(
    r"Data\s+follow\s+in\s+a\s+matrix",
    re.IGNORECASE,
)

ERROR_MARKER_RE = re.compile(
    r"Percentage\s+errors\s+follow\s+in\s+a\s+matrix",
    re.IGNORECASE,
)


def to_float(value: str) -> float:
    """将 Fortran 的 E/D 科学计数格式转换为 Python 浮点数。"""
    return float(value.replace("D", "E").replace("d", "e"))


def resolve_case_path(filename: str | Path) -> Path:
    """将 CASES 中的字符串或 Path 转换为实际路径。"""
    path = Path(filename)

    if not path.is_absolute():
        path = BASE_DIR / path

    return path


def get_reduction_mode() -> str:
    """解析新版模式和旧版 average_x_and_y 兼容开关。"""
    if average_x_and_y is None:
        return TRANSVERSE_REDUCTION_MODE

    return "average" if average_x_and_y else "area_integral"


def extract_first_n_floats(
    text: str,
    expected: int,
    context: str,
) -> np.ndarray:
    """
    从给定文本中顺序读取前 expected 个浮点数。

    只读取指定数量，避免把矩阵后的其他数字或下一段说明误读为数据。
    """
    tokens = FLOAT_RE.findall(text)

    if len(tokens) < expected:
        raise ValueError(
            f"{context}的数据数量不足："
            f"实际找到 {len(tokens)} 个，预期 {expected} 个。"
        )

    return np.asarray(
        [to_float(token) for token in tokens[:expected]],
        dtype=float,
    )


# ============================================================
# 4. 配置检查
# ============================================================
def validate_configuration() -> None:
    """检查用户配置是否有效。"""
    reduction_mode = get_reduction_mode()

    if reduction_mode not in {
        "average",
        "area_integral",
        "sum",
    }:
        raise ValueError(
            "TRANSVERSE_REDUCTION_MODE 必须为 "
            "'average'、'area_integral' 或 'sum'。"
        )

    if REFERENCE_CURRENT_MA <= 0:
        raise ValueError(
            "REFERENCE_CURRENT_MA 必须大于 0。"
        )

    if MARK_EVERY < 1:
        raise ValueError(
            "MARK_EVERY 必须大于或等于 1。"
        )

    if ERROR_BAR_EVERY < 1:
        raise ValueError(
            "ERROR_BAR_EVERY 必须大于或等于 1。"
        )

    if Y_HEADROOM_FACTOR <= 1.0:
        raise ValueError(
            "Y_HEADROOM_FACTOR 建议设置为大于 1。"
        )

    if LOG_Y_BOTTOM_FACTOR <= 0:
        raise ValueError(
            "LOG_Y_BOTTOM_FACTOR 必须大于 0。"
        )

    if TOTAL_MODE not in {
        "integral",
        "sum",
    }:
        raise ValueError(
            "TOTAL_MODE 必须为 'integral' 或 'sum'。"
        )

    if (
        PLOT_Z_MIN is not None
        and PLOT_Z_MAX is not None
        and PLOT_Z_MAX <= PLOT_Z_MIN
    ):
        raise ValueError(
            "PLOT_Z_MAX 必须大于 PLOT_Z_MIN。"
        )

    if (
        INTEGRAL_Z_MIN is not None
        and INTEGRAL_Z_MAX is not None
        and INTEGRAL_Z_MAX <= INTEGRAL_Z_MIN
    ):
        raise ValueError(
            "INTEGRAL_Z_MAX 必须大于 INTEGRAL_Z_MIN。"
        )

    if (
        USE_LOG_X
        and PLOT_Z_MIN is not None
        and PLOT_Z_MIN <= 0
    ):
        raise ValueError(
            "对数 X 轴要求 PLOT_Z_MIN 大于 0。"
        )

    if not CASES:
        raise ValueError(
            "CASES 不能为空。"
        )


    if GENERATE_XY_HEATMAPS:
        # --- 维度选择验证 ---
        valid_dimensions = {"xy", "xz", "yz"}

        if not HEATMAP_DIMENSIONS:
            raise ValueError(
                "启用热图时，HEATMAP_DIMENSIONS 不能为空。"
            )

        for dim in HEATMAP_DIMENSIONS:
            if dim not in valid_dimensions:
                raise ValueError(
                    f"HEATMAP_DIMENSIONS 中的 '{dim}' 无效，"
                    "必须为 'xy'、'xz' 或 'yz'。"
                )

        # 按维度验证固定轴坐标列表。
        if "xy" in HEATMAP_DIMENSIONS:
            if not HEATMAP_Z_POSITIONS_CM:
                raise ValueError(
                    "选择 'xy' 维度时，"
                    "HEATMAP_Z_POSITIONS_CM 不能为空。"
                )

            for z_position in HEATMAP_Z_POSITIONS_CM:
                if not np.isfinite(float(z_position)):
                    raise ValueError(
                        "HEATMAP_Z_POSITIONS_CM 中的 z 坐标"
                        "必须为有限数值。"
                    )

        if "xz" in HEATMAP_DIMENSIONS:
            if not HEATMAP_Y_POSITIONS_CM:
                raise ValueError(
                    "选择 'xz' 维度时，"
                    "HEATMAP_Y_POSITIONS_CM 不能为空。"
                )

            for y_position in HEATMAP_Y_POSITIONS_CM:
                if not np.isfinite(float(y_position)):
                    raise ValueError(
                        "HEATMAP_Y_POSITIONS_CM 中的 y 坐标"
                        "必须为有限数值。"
                    )

        if "yz" in HEATMAP_DIMENSIONS:
            if not HEATMAP_X_POSITIONS_CM:
                raise ValueError(
                    "选择 'yz' 维度时，"
                    "HEATMAP_X_POSITIONS_CM 不能为空。"
                )

            for x_position in HEATMAP_X_POSITIONS_CM:
                if not np.isfinite(float(x_position)):
                    raise ValueError(
                        "HEATMAP_X_POSITIONS_CM 中的 x 坐标"
                        "必须为有限数值。"
                    )

        if HEATMAP_Z_SELECTION_MODE not in {
            "containing_bin",
            "nearest_center",
        }:
            raise ValueError(
                "HEATMAP_Z_SELECTION_MODE 必须为 "
                "'containing_bin' 或 'nearest_center'。"
            )

        if HEATMAP_CASE_INDICES is not None:
            if not HEATMAP_CASE_INDICES:
                raise ValueError(
                    "HEATMAP_CASE_INDICES 不能是空列表；"
                    "如需全部工况请设为 None。"
                )

            for case_index in HEATMAP_CASE_INDICES:
                if (
                    not isinstance(case_index, int)
                    or isinstance(case_index, bool)
                    or case_index < 1
                    or case_index > len(CASES)
                ):
                    raise ValueError(
                        "HEATMAP_CASE_INDICES 中的编号必须是 "
                        f"1 到 {len(CASES)} 之间的整数。"
                    )

        if (
            len(HEATMAP_FIGSIZE) != 2
            or HEATMAP_FIGSIZE[0] <= 0
            or HEATMAP_FIGSIZE[1] <= 0
        ):
            raise ValueError(
                "HEATMAP_FIGSIZE 必须包含两个大于 0 的数值。"
            )

        if HEATMAP_DPI <= 0:
            raise ValueError(
                "HEATMAP_DPI 必须大于 0。"
            )

        allowed_formats = {"png", "pdf", "svg"}
        selected_formats = {
            str(item).lower()
            for item in HEATMAP_OUTPUT_FORMATS
        }

        if not selected_formats:
            raise ValueError(
                "HEATMAP_OUTPUT_FORMATS 不能为空。"
            )

        if not selected_formats.issubset(allowed_formats):
            raise ValueError(
                "HEATMAP_OUTPUT_FORMATS 只能包含 "
                "png、pdf 和 svg。"
            )

        if HEATMAP_VMIN is not None:
            if not np.isfinite(float(HEATMAP_VMIN)):
                raise ValueError(
                    "HEATMAP_VMIN 必须为有限数值或 None。"
                )

            if HEATMAP_USE_LOG_COLOR and HEATMAP_VMIN <= 0:
                raise ValueError(
                    "对数热图要求 HEATMAP_VMIN 大于 0。"
                )

        if HEATMAP_VMAX is not None:
            if not np.isfinite(float(HEATMAP_VMAX)):
                raise ValueError(
                    "HEATMAP_VMAX 必须为有限数值或 None。"
                )

            if HEATMAP_USE_LOG_COLOR and HEATMAP_VMAX <= 0:
                raise ValueError(
                    "对数热图要求 HEATMAP_VMAX 大于 0。"
                )

        if (
            HEATMAP_VMIN is not None
            and HEATMAP_VMAX is not None
            and HEATMAP_VMAX <= HEATMAP_VMIN
        ):
            raise ValueError(
                "HEATMAP_VMAX 必须大于 HEATMAP_VMIN。"
            )

        # 提前验证 Matplotlib 是否认识该色图名称。
        plt.get_cmap(HEATMAP_COLORMAP)

        # --- 合并参数验证 ---
        if HEATMAP_MERGE_NCOLS < 1:
            raise ValueError(
                "HEATMAP_MERGE_NCOLS 必须大于或等于 1。"
            )

        if HEATMAP_MERGE_DIFFERENTIATE:
            for cmap_name in (
                HEATMAP_MERGE_COLORMAPS
            ):
                plt.get_cmap(cmap_name)

            for alpha in (
                HEATMAP_MERGE_ALPHAS
            ):
                if not (
                    0.0
                    < float(alpha)
                    <= 1.0
                ):
                    raise ValueError(
                        "HEATMAP_MERGE_ALPHAS 中的每个值"
                        "必须在 (0, 1] 范围内。"
                    )


def validate_case_filenames() -> None:
    """检查 CASES 中的文件名、电流和人工缩放倍率。"""
    empty_cases: list[int] = []

    for index, case in enumerate(
        CASES,
        start=1,
    ):
        if len(case) != 4:
            raise ValueError(
                f"CASES 第 {index} 项必须包含四个元素："
                "文件名、图例、电流和人工缩放倍率。"
            )

        (
            filename,
            _,
            current_ma,
            user_y_scale,
        ) = case

        if not str(filename).strip():
            empty_cases.append(index)

        if float(current_ma) <= 0:
            raise ValueError(
                f"CASES 第 {index} 项的电流必须大于 0。"
            )

        if not np.isfinite(
            float(user_y_scale)
        ):
            raise ValueError(
                f"CASES 第 {index} 项的人工缩放倍率必须有限。"
            )

    if empty_cases:
        case_text = ", ".join(
            str(index)
            for index in empty_cases
        )

        raise ValueError(
            "以下 CASES 尚未填写 USRBIN 文件名："
            f"{case_text}"
        )


# ============================================================
# 5. USRBIN usbrea .lis 文件解析
# ============================================================
def parse_usrbin_detectors(
    path: Path,
) -> list[dict]:
    """
    读取一个 usbrea .lis 文件中的全部 USRBIN Detector。

    每个 Detector 返回：
        number
        name
        binning_type
        header_tail
        particle_number
        track_length_binning
        axes
        score
        error_percent
        has_error_matrix
    """
    if not path.is_file():
        raise FileNotFoundError(
            f"找不到输入文件：\n{path}"
        )

    text = path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    header_matches = list(
        BINNING_HEADER_RE.finditer(text)
    )

    # 某些简化输出可能缺少完整的 binning 标题。
    # 此时将整个文件作为一个 Detector 块尝试解析。
    if not header_matches:
        block_specs = [
            {
                "start": 0,
                "end": len(text),
                "binning_type": "Cartesian",
                "number": 1,
                "name": path.stem,
                "header_tail": "",
            }
        ]

    else:
        block_specs = []

        for index, match in enumerate(
            header_matches
        ):
            block_specs.append({
                "start": match.start(),
                "end": (
                    header_matches[
                        index + 1
                    ].start()
                    if index + 1
                    < len(header_matches)
                    else len(text)
                ),
                "binning_type": (
                    match.group(1)
                    .strip()
                    .title()
                ),
                "number": int(
                    match.group(2)
                ),
                "name": (
                    match.group(3) or ""
                ).strip(),
                "header_tail": (
                    match.group(4) or ""
                ).strip(),
            })

    detectors: list[dict] = []

    for spec in block_specs:
        block = text[
            spec["start"]:spec["end"]
        ]

        if (
            spec["binning_type"].lower()
            != "cartesian"
        ):
            raise ValueError(
                f"{path.name} 中 Detector "
                f"{spec['number']} 为 "
                f"{spec['binning_type']} binning。"
                "当前程序的 Z 向剖面功能仅支持 "
                "Cartesian X-Y-Z 网格。"
            )

        axes: dict[str, dict] = {}

        for match in AXIS_RE.finditer(
            block
        ):
            axis = match.group(1).upper()

            axes[axis] = {
                "min": to_float(
                    match.group(2)
                ),
                "max": to_float(
                    match.group(3)
                ),
                "n": int(
                    match.group(4)
                ),
                "width": to_float(
                    match.group(5)
                ),
            }

        if set(axes) != {
            "X",
            "Y",
            "Z",
        }:
            raise ValueError(
                f"无法从 {path.name} 的 Detector "
                f"{spec['number']} 中完整识别 "
                "X、Y、Z 网格信息。"
            )

        for axis_name, axis in axes.items():
            if axis["n"] < 1:
                raise ValueError(
                    f"{path.name} Detector "
                    f"{spec['number']} 的 "
                    f"{axis_name} bin 数必须大于 0。"
                )

            if axis["max"] <= axis["min"]:
                raise ValueError(
                    f"{path.name} Detector "
                    f"{spec['number']} 的 "
                    f"{axis_name} 最大值必须大于最小值。"
                )

            if axis["width"] <= 0:
                raise ValueError(
                    f"{path.name} Detector "
                    f"{spec['number']} 的 "
                    f"{axis_name} bin 宽度必须大于 0。"
                )

        nx = int(
            axes["X"]["n"]
        )

        ny = int(
            axes["Y"]["n"]
        )

        nz = int(
            axes["Z"]["n"]
        )

        expected = nx * ny * nz

        data_marker = DATA_MARKER_RE.search(
            block
        )

        if data_marker is None:
            raise ValueError(
                f"无法在 {path.name} 的 Detector "
                f"{spec['number']} 中找到 "
                "'Data follow in a matrix'。"
            )

        data_line_end = block.find(
            "\n",
            data_marker.end(),
        )

        if data_line_end < 0:
            raise ValueError(
                f"{path.name} Detector "
                f"{spec['number']} 的数据标志行不完整。"
            )

        error_marker = ERROR_MARKER_RE.search(
            block,
            data_line_end + 1,
        )

        if error_marker is None:
            data_section = block[
                data_line_end + 1:
            ]
        else:
            data_section = block[
                data_line_end + 1:
                error_marker.start()
            ]

        values = extract_first_n_floats(
            data_section,
            expected,
            f"{path.name} Detector "
            f"{spec['number']} 评分矩阵",
        )

        score = values.reshape(
            (nx, ny, nz),
            order="F",
        )

        if error_marker is None:
            if REQUIRE_ERROR_MATRIX:
                raise ValueError(
                    f"{path.name} Detector "
                    f"{spec['number']} 缺少百分比误差矩阵。"
                )

            error_percent = np.full(
                score.shape,
                np.nan,
                dtype=float,
            )

            has_error_matrix = False

        else:
            error_line_end = block.find(
                "\n",
                error_marker.end(),
            )

            if error_line_end < 0:
                raise ValueError(
                    f"{path.name} Detector "
                    f"{spec['number']} 的误差标志行不完整。"
                )

            error_section = block[
                error_line_end + 1:
            ]

            error_values = extract_first_n_floats(
                error_section,
                expected,
                f"{path.name} Detector "
                f"{spec['number']} 百分比误差矩阵",
            )

            error_percent = error_values.reshape(
                (nx, ny, nz),
                order="F",
            )

            has_error_matrix = True

        particle_match = PARTICLE_NUMBER_RE.search(
            spec["header_tail"]
            + "\n"
            + block[:data_marker.start()]
        )

        particle_number = (
            int(
                particle_match.group(1)
            )
            if particle_match is not None
            else None
        )

        track_length_binning = (
            "this is a track-length binning"
            in block.lower()
        )

        detectors.append({
            "number": spec["number"],
            "name": spec["name"],
            "binning_type": (
                spec["binning_type"]
            ),
            "header_tail": (
                spec["header_tail"]
            ),
            "particle_number": (
                particle_number
            ),
            "track_length_binning": (
                track_length_binning
            ),
            "axes": axes,
            "score": score,
            "error_percent": (
                error_percent
            ),
            "has_error_matrix": (
                has_error_matrix
            ),
        })

    if not detectors:
        raise ValueError(
            "无法在以下文件中找到 USRBIN Detector：\n"
            f"{path}"
        )

    return detectors


def select_detector(
    detectors: list[dict],
    path: Path,
) -> dict:
    """根据编号和名称关键词选择唯一 Detector。"""
    matched: list[dict] = []

    for detector in detectors:
        number_matches = (
            DETECTOR_NUMBER is None
            or detector["number"]
            == DETECTOR_NUMBER
        )

        if DETECTOR_NAME_CONTAINS is None:
            name_matches = True
        else:
            name_matches = (
                DETECTOR_NAME_CONTAINS.lower()
                in detector["name"].lower()
            )

        if (
            number_matches
            and name_matches
        ):
            matched.append(
                detector
            )

    available_text = "\n".join(
        "  Detector "
        f"{item['number']}: "
        f"{item['name'] or '(unnamed)'}"
        for item in detectors
    )

    if len(matched) == 1:
        return matched[0]

    if not matched:
        raise ValueError(
            "在以下文件中没有找到符合条件的 "
            "USRBIN Detector：\n"
            f"{path}\n\n"
            "当前文件中的 Detector 为：\n"
            f"{available_text}"
        )

    raise ValueError(
        "Detector 筛选条件匹配到多个结果：\n"
        f"{path}\n\n"
        "请设置 DETECTOR_NUMBER 和/或 "
        "DETECTOR_NAME_CONTAINS，"
        "以唯一确定 Detector。\n\n"
        f"匹配结果：\n{available_text}"
    )


def read_usrbin_lis(
    path: Path,
) -> dict:
    """读取并返回用户指定的单个 USRBIN Detector。"""
    detectors = parse_usrbin_detectors(
        path
    )

    detector = select_detector(
        detectors,
        path,
    )

    # 返回副本，避免后续处理意外修改解析结果。
    return {
        **detector,
        "axes": {
            key: value.copy()
            for key, value
            in detector["axes"].items()
        },
        "score": (
            detector["score"].copy()
        ),
        "error_percent": (
            detector["error_percent"].copy()
        ),
    }


# ============================================================
# 6. 三维矩阵约化为 Z 向一维 profile
# ============================================================
def get_axis_edges(
    axis: dict,
) -> np.ndarray:
    """根据轴上下限和 bin 数生成精确边界。"""
    return np.linspace(
        float(axis["min"]),
        float(axis["max"]),
        int(axis["n"]) + 1,
        dtype=float,
    )


def build_longitudinal_profile(
    axes: dict,
    score: np.ndarray,
    error_percent_3d: np.ndarray,
) -> dict:
    """
    将 A(ix,iy,iz) 沿 X-Y 约化为 Z 向 profile，并传播误差。

    返回：
        z_min, z_max, z_center, dz
        raw_profile
        raw_absolute_error
        profile_error_percent
        reduction_mode
        transverse_factor
    """
    nx = int(
        axes["X"]["n"]
    )

    ny = int(
        axes["Y"]["n"]
    )

    nz = int(
        axes["Z"]["n"]
    )

    if score.shape != (
        nx,
        ny,
        nz,
    ):
        raise ValueError(
            f"score 形状 {score.shape} 与网格 "
            f"({nx}, {ny}, {nz}) 不一致。"
        )

    if (
        error_percent_3d.shape
        != score.shape
    ):
        raise ValueError(
            "误差矩阵形状与评分矩阵不一致。"
        )

    x_edges = get_axis_edges(
        axes["X"]
    )

    y_edges = get_axis_edges(
        axes["Y"]
    )

    z_edges = get_axis_edges(
        axes["Z"]
    )

    dx = float(
        np.mean(
            np.diff(x_edges)
        )
    )

    dy = float(
        np.mean(
            np.diff(y_edges)
        )
    )

    dz_values = np.diff(
        z_edges
    )

    z_min = z_edges[:-1]
    z_max = z_edges[1:]

    z_center = (
        z_min + z_max
    ) / 2.0

    reduction_mode = get_reduction_mode()

    if reduction_mode == "average":
        transverse_factor = (
            1.0 / (nx * ny)
        )

    elif reduction_mode == "area_integral":
        transverse_factor = dx * dy

    else:
        transverse_factor = 1.0

    if np.any(
        ~np.isfinite(score)
    ):
        raise ValueError(
            "USRBIN 评分矩阵包含 NaN 或无穷值。"
        )

    raw_profile = (
        np.sum(
            score,
            axis=(0, 1),
        )
        * transverse_factor
    )

    # 将三维相对误差转换为三维绝对标准差。
    valid_error = (
        np.isfinite(error_percent_3d)
        & (error_percent_3d >= 0)
    )

    absolute_error_3d = np.where(
        valid_error,
        np.abs(score)
        * error_percent_3d
        / 100.0,
        0.0,
    )

    # 假设各 X-Y bin 独立，方差相加。
    raw_absolute_error = (
        np.sqrt(
            np.sum(
                absolute_error_3d ** 2,
                axis=(0, 1),
            )
        )
        * abs(transverse_factor)
    )

    # 只要某个有效 score bin 没有误差信息，
    # 该 z 层误差就标为未知。
    complete_error_for_z = np.all(
        valid_error,
        axis=(0, 1),
    )

    raw_absolute_error = (
        raw_absolute_error.astype(float)
    )

    raw_absolute_error[
        ~complete_error_for_z
    ] = np.nan

    profile_error_percent = np.full(
        raw_profile.shape,
        np.nan,
        dtype=float,
    )

    nonzero_profile = (
        np.isfinite(raw_profile)
        & (raw_profile != 0)
    )

    valid_profile_error = (
        nonzero_profile
        & np.isfinite(
            raw_absolute_error
        )
    )

    profile_error_percent[
        valid_profile_error
    ] = (
        raw_absolute_error[
            valid_profile_error
        ]
        / np.abs(
            raw_profile[
                valid_profile_error
            ]
        )
        * 100.0
    )

    zero_with_known_error = (
        (raw_profile == 0)
        & np.isfinite(
            raw_absolute_error
        )
        & (raw_absolute_error == 0)
    )

    profile_error_percent[
        zero_with_known_error
    ] = 0.0

    return {
        "z_min": z_min,
        "z_max": z_max,
        "z_center": z_center,
        "dz": dz_values,
        "raw_profile": raw_profile,
        "raw_absolute_error": (
            raw_absolute_error
        ),
        "profile_error_percent": (
            profile_error_percent
        ),
        "reduction_mode": (
            reduction_mode
        ),
        "transverse_factor": (
            transverse_factor
        ),
        "dx": dx,
        "dy": dy,
    }


# ============================================================
# 7. 绘图数据处理
# ============================================================
def get_plot_range_mask(
    z_min: np.ndarray,
    z_max: np.ndarray,
) -> np.ndarray:
    """返回与用户绘图 Z 范围相交的 bin。"""
    mask = np.ones(
        z_min.shape,
        dtype=bool,
    )

    if PLOT_Z_MIN is not None:
        mask &= (
            z_max > PLOT_Z_MIN
        )

    if PLOT_Z_MAX is not None:
        mask &= (
            z_min < PLOT_Z_MAX
        )

    if USE_LOG_X:
        mask &= (
            z_max > 0
        )

    return mask


def prepare_profile_for_plot(
    z_min: np.ndarray,
    z_max: np.ndarray,
    profile: np.ndarray,
    error_percent: np.ndarray,
) -> np.ndarray:
    """生成绘图用 profile 副本并应用范围、对数轴和统计屏蔽。"""
    plotted = (
        profile.copy()
        .astype(float)
    )

    plotted[
        ~np.isfinite(plotted)
    ] = np.nan

    range_mask = get_plot_range_mask(
        z_min,
        z_max,
    )

    plotted[
        ~range_mask
    ] = np.nan

    if USE_LOG_Y:
        plotted[
            plotted <= 0
        ] = np.nan

    if MASK_INSUFFICIENT_STATISTICS:
        insufficient_mask = (
            np.isfinite(error_percent)
            & (
                error_percent
                >= INSUFFICIENT_ERROR_PERCENT
            )
            & np.isfinite(profile)
            & (profile != 0)
        )

        plotted[
            insufficient_mask
        ] = np.nan

    return plotted


def get_error_bar_mask(
    profile_for_plot: np.ndarray,
    error_percent: np.ndarray,
) -> np.ndarray:
    """返回实际绘制误差棒的数据点掩码。"""
    mask = (
        np.isfinite(profile_for_plot)
        & np.isfinite(error_percent)
        & (error_percent >= 0)
    )

    if SKIP_INSUFFICIENT_ERROR_BARS:
        mask &= (
            error_percent
            < INSUFFICIENT_ERROR_PERCENT
        )

    interval_mask = np.zeros(
        profile_for_plot.shape,
        dtype=bool,
    )

    interval_mask[
        np.arange(
            profile_for_plot.size
        )[::ERROR_BAR_EVERY]
    ] = True

    return mask & interval_mask


def get_visible_y_limits(
    z_min: np.ndarray,
    z_max: np.ndarray,
    profile: np.ndarray,
    error_percent: np.ndarray,
) -> tuple[float, float]:
    """获取一条曲线实际显示的最小值和最大值。"""
    plotted = prepare_profile_for_plot(
        z_min,
        z_max,
        profile,
        error_percent,
    )

    finite_mask = np.isfinite(
        plotted
    )

    if not np.any(finite_mask):
        return np.inf, -np.inf

    finite_values = plotted[
        finite_mask
    ]

    if USE_LOG_Y:
        positive_values = finite_values[
            finite_values > 0
        ]

        minimum_value = (
            float(
                np.min(
                    positive_values
                )
            )
            if positive_values.size > 0
            else np.inf
        )

    else:
        minimum_value = float(
            np.min(
                finite_values
            )
        )

    maximum_value = float(
        np.max(
            finite_values
        )
    )

    if SHOW_ERROR_BARS:
        error_mask = get_error_bar_mask(
            plotted,
            error_percent,
        )

        if np.any(error_mask):
            y_values = plotted[
                error_mask
            ]

            absolute_error = (
                np.abs(y_values)
                * error_percent[
                    error_mask
                ]
                / 100.0
            )

            upper_values = (
                y_values
                + absolute_error
            )

            lower_values = (
                y_values
                - absolute_error
            )

            finite_upper = upper_values[
                np.isfinite(
                    upper_values
                )
            ]

            if finite_upper.size > 0:
                maximum_value = max(
                    maximum_value,
                    float(
                        np.max(
                            finite_upper
                        )
                    ),
                )

            if USE_LOG_Y:
                valid_lower = lower_values[
                    np.isfinite(
                        lower_values
                    )
                    & (lower_values > 0)
                ]

            else:
                valid_lower = lower_values[
                    np.isfinite(
                        lower_values
                    )
                ]

            if valid_lower.size > 0:
                minimum_value = min(
                    minimum_value,
                    float(
                        np.min(
                            valid_lower
                        )
                    ),
                )

    return (
        minimum_value,
        maximum_value,
    )


def plot_one_profile(
    ax,
    z_min: np.ndarray,
    z_max: np.ndarray,
    z_center: np.ndarray,
    profile: np.ndarray,
    error_percent: np.ndarray,
    label: str,
    color,
) -> None:
    """绘制一条 USRBIN Z 向 profile。"""
    profile_for_plot = prepare_profile_for_plot(
        z_min,
        z_max,
        profile,
        error_percent,
    )

    if not np.any(
        np.isfinite(
            profile_for_plot
        )
    ):
        return

    if USE_CONTINUOUS_LINE:
        ax.plot(
            z_center,
            profile_for_plot,
            color=color,
            linewidth=2.0,
            label=label,
            zorder=3,
        )

    else:
        z_edges = np.concatenate([
            z_min[:1],
            z_max,
        ])

        ax.stairs(
            values=profile_for_plot,
            edges=z_edges,
            color=color,
            linewidth=2.0,
            label=label,
            zorder=3,
        )

    if SHOW_BIN_MARKERS:
        marker_mask = np.zeros(
            z_center.shape,
            dtype=bool,
        )

        marker_mask[
            np.arange(
                z_center.size
            )[::MARK_EVERY]
        ] = True

        marker_mask &= np.isfinite(
            profile_for_plot
        )

        ax.plot(
            z_center[
                marker_mask
            ],
            profile_for_plot[
                marker_mask
            ],
            linestyle="none",
            marker=BIN_MARKER,
            markersize=(
                BIN_MARKER_SIZE
            ),
            color=color,
            label="_nolegend_",
            zorder=4,
        )

    if (
        not USE_CONTINUOUS_LINE
        and ADD_GUIDE_LINE
    ):
        ax.plot(
            z_center,
            profile_for_plot,
            color=color,
            linewidth=(
                GUIDE_LINE_WIDTH
            ),
            alpha=(
                GUIDE_LINE_ALPHA
            ),
            solid_capstyle="round",
            label="_nolegend_",
            zorder=4,
        )

    if SHOW_ERROR_BARS:
        error_mask = get_error_bar_mask(
            profile_for_plot,
            error_percent,
        )

        if np.any(error_mask):
            x_values = z_center[
                error_mask
            ]

            y_values = profile_for_plot[
                error_mask
            ]

            relative_errors = (
                error_percent[
                    error_mask
                ]
            )

            absolute_error = (
                np.abs(y_values)
                * relative_errors
                / 100.0
            )

            if USE_LOG_Y:
                lower_error = np.minimum(
                    absolute_error,
                    y_values
                    * (1.0 - 1.0e-12),
                )

                y_error = np.vstack([
                    lower_error,
                    absolute_error,
                ])

            else:
                y_error = (
                    absolute_error
                )

            ax.errorbar(
                x_values,
                y_values,
                yerr=y_error,
                fmt="none",
                ecolor=color,
                elinewidth=(
                    ERROR_BAR_LINEWIDTH
                ),
                capsize=(
                    ERROR_BAR_CAPSIZE
                ),
                alpha=(
                    ERROR_BAR_ALPHA
                ),
                label="_nolegend_",
                zorder=2,
            )


# ============================================================
# 8. Z 向积分或求和
# ============================================================
def calculate_profile_total(
    z_min: np.ndarray,
    z_max: np.ndarray,
    profile: np.ndarray,
) -> tuple[float, float, float]:
    """在指定 Z 范围内对 profile 积分或求和。"""
    lower_limit = (
        float(
            np.min(z_min)
        )
        if INTEGRAL_Z_MIN is None
        else float(
            INTEGRAL_Z_MIN
        )
    )

    upper_limit = (
        float(
            np.max(z_max)
        )
        if INTEGRAL_Z_MAX is None
        else float(
            INTEGRAL_Z_MAX
        )
    )

    if upper_limit <= lower_limit:
        raise ValueError(
            "积分 Z 上限必须大于下限。"
        )

    overlap_left = np.maximum(
        z_min,
        lower_limit,
    )

    overlap_right = np.minimum(
        z_max,
        upper_limit,
    )

    overlap_width = np.maximum(
        overlap_right
        - overlap_left,
        0.0,
    )

    valid_mask = (
        np.isfinite(profile)
        & (overlap_width > 0)
    )

    if not np.any(valid_mask):
        raise ValueError(
            "指定 Z 范围内没有有效 USRBIN 数据。"
        )

    if TOTAL_MODE == "integral":
        total_value = float(
            np.sum(
                profile[
                    valid_mask
                ]
                * overlap_width[
                    valid_mask
                ]
            )
        )

    else:
        total_value = float(
            np.sum(
                profile[
                    valid_mask
                ]
            )
        )

    return (
        total_value,
        lower_limit,
        upper_limit,
    )


# ============================================================
# 9. 网格检查与 CSV 名称
# ============================================================
def axes_are_equal(
    axes_a: dict,
    axes_b: dict,
    include_xy: bool,
) -> bool:
    """检查两个 Detector 的空间网格是否一致。"""
    axis_names = (
        ("X", "Y", "Z")
        if include_xy
        else ("Z",)
    )

    for axis_name in axis_names:
        a = axes_a[
            axis_name
        ]

        b = axes_b[
            axis_name
        ]

        if int(a["n"]) != int(b["n"]):
            return False

        if not np.isclose(
            a["min"],
            b["min"],
            rtol=1.0e-8,
            atol=1.0e-12,
        ):
            return False

        if not np.isclose(
            a["max"],
            b["max"],
            rtol=1.0e-8,
            atol=1.0e-12,
        ):
            return False

    return True


def make_clean_csv_name(
    path: Path,
    case_index: int,
) -> str:
    """生成适合 CSV 列名的工况名称。"""
    clean_stem = (
        path.stem
        .replace(",", "_")
        .replace(" ", "_")
        .replace("(", "")
        .replace(")", "")
        .replace("-", "_")
    )

    clean_stem = re.sub(
        r"_+",
        "_",
        clean_stem,
    ).strip("_")

    return (
        f"case_{case_index:02d}_"
        f"{clean_stem}"
    )


# ============================================================
# 10. 二维截面热图
# ============================================================

# 各维度组合的配置信息。
# free_axes：热图水平轴和垂直轴对应的 USRBIN 轴名；
# fixed_axis：被固定的轴名（沿该轴选取截面）；
# xlabel / ylabel：热图坐标轴标签。
HEATMAP_DIMENSION_INFO = {
    "xy": {
        "free_axes": ("X", "Y"),
        "fixed_axis": "Z",
        "xlabel": r"$x$ (cm)",
        "ylabel": r"$y$ (cm)",
    },
    "xz": {
        "free_axes": ("X", "Z"),
        "fixed_axis": "Y",
        "xlabel": r"$x$ (cm)",
        "ylabel": r"$z$ (cm)",
    },
    "yz": {
        # 物理习惯：Z（束流方向）水平，Y（横向）垂直。
        "free_axes": ("Z", "Y"),
        "fixed_axis": "X",
        "xlabel": r"$z$ (cm)",
        "ylabel": r"$y$ (cm)",
    },
}


# ============================================================
# 10.2 FLUKA 几何边界叠加辅助函数
# ============================================================
# 模块级缓存：.inp 只解析一次
_FLUKA_GEOMETRY_CACHE = None


def _liang_barsky(x0, y0, x1, y1, xmin, ymin, xmax, ymax):
    """Liang-Barsky 直线裁剪算法。

    将直线 (x0,y0)->(x1,y1) 裁剪到矩形 [xmin,xmax]×[ymin,ymax]。
    返回 (nx0, ny0, nx1, ny1) 或 None（直线在矩形外）。
    """
    dx = x1 - x0
    dy = y1 - y0
    p = [-dx, dx, -dy, dy]
    q = [x0 - xmin, xmax - x0, y0 - ymin, ymax - y0]
    u1, u2 = 0.0, 1.0
    for i in range(4):
        if abs(p[i]) < 1e-12:
            if q[i] < 0:
                return None
        else:
            t = q[i] / p[i]
            if p[i] < 0:
                if t > u2:
                    return None
                if t > u1:
                    u1 = t
            else:
                if t < u1:
                    return None
                if t < u2:
                    u2 = t
    if u1 > u2:
        return None
    return (x0 + u1 * dx, y0 + u1 * dy,
            x0 + u2 * dx, y0 + u2 * dy)


def _body_outward_normal(body, x, y, z, eps=1e-6):
    """计算 body 在 (x,y,z) 处的外法向量（signed_distance 增加方向）。"""
    sd_xp = body.signed_distance(x + eps, y, z)
    sd_xm = body.signed_distance(x - eps, y, z)
    sd_yp = body.signed_distance(x, y + eps, z)
    sd_ym = body.signed_distance(x, y - eps, z)
    sd_zp = body.signed_distance(x, y, z + eps)
    sd_zm = body.signed_distance(x, y, z - eps)
    dx = sd_xp - sd_xm
    dy = sd_yp - sd_ym
    dz = sd_zp - sd_zm
    norm = np.sqrt(dx**2 + dy**2 + dz**2)
    norm = np.where(norm < 1e-12, 1.0, norm)
    return dx / norm, dy / norm, dz / norm


def _get_fluka_geometry():
    """解析 .inp 文件，结果缓存。返回 Geometry 或 None。"""
    global _FLUKA_GEOMETRY_CACHE
    if not FLUKA_GEOMETRY_ENABLED or not _FLUKA_GEOMETRY_AVAILABLE:
        return None
    if _FLUKA_GEOMETRY_CACHE is None:
        inp_path = Path(FLUKA_GEOMETRY_INP)
        if not inp_path.is_absolute():
            inp_path = BASE_DIR / inp_path
        if not inp_path.exists():
            print(f"  [FLUKA几何] 警告: .inp 文件不存在: {inp_path}")
            return None
        try:
            _FLUKA_GEOMETRY_CACHE = _fluka_parse_inp(inp_path)
            print(f"  [FLUKA几何] 已解析: {len(_FLUKA_GEOMETRY_CACHE.bodies)} bodies, "
                  f"{len(_FLUKA_GEOMETRY_CACHE.regions)} regions")
        except Exception as e:
            print(f"  [FLUKA几何] 解析失败: {e}")
            return None
    return _FLUKA_GEOMETRY_CACHE


def _get_geometry_cut_value(dimension, requested_pos_cm):
    """根据维度和热图切片位置，返回 (cut_axis, cut_value)。"""
    dim_info = HEATMAP_DIMENSION_INFO[dimension]
    fixed_axis = dim_info["fixed_axis"]
    cut_axis = fixed_axis.lower()

    if FLUKA_GEOMETRY_CUT_PLANE is not None:
        if cut_axis in FLUKA_GEOMETRY_CUT_PLANE:
            return cut_axis, float(FLUKA_GEOMETRY_CUT_PLANE[cut_axis])

    return cut_axis, float(requested_pos_cm)


def _draw_geometry_on_heatmap(
    ax,
    dimension: str,
    requested_pos_cm: float,
    h_edges,
    v_edges,
    output_stem=None,
    n_samples=2000,
    eps=1e-4,
) -> None:
    """在热图 ax 上叠加 FLUKA 几何边界线（仅画线，不涂色）。

    边界被裁剪到 region 范围内：对每条 body 边界曲线采样后，
    沿 body 法向量偏移 ±eps 检查两侧是否在 region 内，
    任一侧在 region 内则画该点。
    """
    geom = _get_fluka_geometry()
    if geom is None:
        return

    cut_axis, cut_value = _get_geometry_cut_value(dimension, requested_pos_cm)

    dim_info = HEATMAP_DIMENSION_INFO[dimension]
    h_axis, v_axis = dim_info["free_axes"]
    axis_to_name = {"X": "x", "Y": "y", "Z": "z"}
    h_name = axis_to_name[h_axis]
    v_name = axis_to_name[v_axis]

    # cut_plane 内部 free1/free2 的物理轴名
    cut_free1_name, cut_free2_name = {
        "x": ("y", "z"),
        "y": ("x", "z"),
        "z": ("x", "y"),
    }[cut_axis]

    # 映射物理名 -> edges
    edges_map = {h_name: h_edges, v_name: v_edges}
    free_edges_1 = edges_map[cut_free1_name]
    free_edges_2 = edges_map[cut_free2_name]

    h_lo, h_hi = float(h_edges[0]), float(h_edges[-1])
    v_lo, v_hi = float(v_edges[0]), float(v_edges[-1])

    # 切面固定坐标
    if cut_axis == "x":
        x_fixed = cut_value
    elif cut_axis == "y":
        y_fixed = cut_value
    else:
        z_fixed = cut_value

    # 遍历每个 region，对 region 引用的每个 body 画裁剪后的边界
    for region in geom.regions:
        if region.material in FLUKA_GEOMETRY_SKIP_MATERIALS:
            continue

        body_names = _fluka_collect_body_refs(region.expr_tree)
        for bn in body_names:
            if bn not in geom.bodies:
                continue
            body = geom.bodies[bn]
            cs = body.analytic_boundary_on_plane(cut_axis, cut_value)
            for c in cs:
                f1 = np.asarray(c["free1"], dtype=float)
                f2 = np.asarray(c["free2"], dtype=float)

                # 映射到画图坐标
                if cut_free1_name == h_name:
                    h_raw, v_raw = f1, f2
                else:
                    h_raw, v_raw = f2, f1

                # 两点直线：先 Liang-Barsky 裁剪到画图范围，再密集采样
                if len(h_raw) == 2:
                    clipped = _liang_barsky(
                        float(h_raw[0]), float(v_raw[0]),
                        float(h_raw[1]), float(v_raw[1]),
                        h_lo, v_lo, h_hi, v_hi
                    )
                    if clipped is None:
                        continue
                    ch0, cv0, ch1, cv1 = clipped
                    t = np.linspace(0.0, 1.0, n_samples)
                    h_data = ch0 + t * (ch1 - ch0)
                    v_data = cv0 + t * (cv1 - cv0)
                else:
                    # 多点曲线：直接 mask
                    h_data = h_raw
                    v_data = v_raw
                    range_mask = (h_data >= h_lo) & (h_data <= h_hi) & \
                                 (v_data >= v_lo) & (v_data <= v_hi)
                    h_data = h_data[range_mask]
                    v_data = v_data[range_mask]
                    if len(h_data) < 2:
                        continue

                # 构造 (x, y, z) 坐标用于 region.inside 检查
                pts_count = len(h_data)
                if cut_free1_name == h_name:
                    f1_s = h_data
                    f2_s = v_data
                else:
                    f1_s = v_data
                    f2_s = h_data

                if cut_axis == "x":
                    X = np.full(pts_count, x_fixed)
                    if cut_free1_name == "y":
                        Y, Z = f1_s, f2_s
                    else:
                        Y, Z = f2_s, f1_s
                elif cut_axis == "y":
                    Y = np.full(pts_count, y_fixed)
                    if cut_free1_name == "x":
                        X, Z = f1_s, f2_s
                    else:
                        X, Z = f2_s, f1_s
                else:  # z
                    Z = np.full(pts_count, z_fixed)
                    if cut_free1_name == "x":
                        X, Y = f1_s, f2_s
                    else:
                        X, Y = f2_s, f1_s

                # 沿 body 法向量偏移 ±eps，检查两侧
                nx, ny, nz = _body_outward_normal(body, X, Y, Z)
                X_in = X - eps * nx
                Y_in = Y - eps * ny
                Z_in = Z - eps * nz
                X_out = X + eps * nx
                Y_out = Y + eps * ny
                Z_out = Z + eps * nz

                mask_in = region.inside(geom.bodies, X_in, Y_in, Z_in)
                mask_out = region.inside(geom.bodies, X_out, Y_out, Z_out)
                region_mask = mask_in | mask_out

                if region_mask.sum() >= 2:
                    ax.plot(
                        h_data[region_mask],
                        v_data[region_mask],
                        linestyle=FLUKA_GEOMETRY_LINE_STYLE,
                        color=FLUKA_GEOMETRY_LINE_COLOR,
                        linewidth=FLUKA_GEOMETRY_LINE_WIDTH,
                        alpha=0.9,
                        zorder=5,
                    )

    # 可选：导出材料分布 CSV
    if FLUKA_GEOMETRY_EXPORT_CSV and output_stem is not None:
        try:
            from fluka_geometry import cut_plane as _fluka_cut_plane
            mg = _fluka_cut_plane(geom, cut_axis, cut_value,
                                   free_edges_1, free_edges_2)
            csv_name = f"{output_stem.name}_materials_{cut_axis}{cut_value:.4f}.csv"
            csv_path = output_stem.with_name(csv_name)
            mg.to_csv(csv_path)
        except Exception as e:
            print(f"  [FLUKA几何] CSV 导出失败: {e}")


def _draw_inset_on_heatmap(
    parent_ax,
    figure,
    mesh,
    norm,
    masked_score_2d,
    h_edges,
    v_edges,
    h_axis_name: str,
    v_axis_name: str,
    dimension: str,
    requested_pos_cm: float,
    xlabel: str,
    ylabel: str,
) -> None:
    """在大图角落画 inset 放大子图 + 大图上画虚线框。

    parent_ax: 大图 axes
    figure: 大图 figure
    mesh: 大图 pcolormesh 返回的 QuadMesh（用于共享 cmap）
    norm: 大图 norm
    masked_score_2d: 大图 masked_score（转置前，即原 (n_h, n_v)）
    h_edges/v_edges: bin 边界
    h_axis_name/v_axis_name: 横纵轴物理名（'x'/'y'/'z'）
    dimension: 'xy'/'xz'/'yz'
    requested_pos_cm: 切面位置
    """
    if not HEATMAP_INSET_ENABLED:
        return

    h_lo, h_hi = float(HEATMAP_INSET_H_RANGE[0]), float(HEATMAP_INSET_H_RANGE[1])
    v_lo, v_hi = float(HEATMAP_INSET_V_RANGE[0]), float(HEATMAP_INSET_V_RANGE[1])

    # 限制到画图范围内
    h_lo = max(h_lo, float(h_edges[0]))
    h_hi = min(h_hi, float(h_edges[-1]))
    v_lo = max(v_lo, float(v_edges[0]))
    v_hi = min(v_hi, float(v_edges[-1]))
    if h_hi <= h_lo or v_hi <= v_lo:
        return

    # ---- 1. 大图上画虚线框 ----
    parent_ax.plot(
        [h_lo, h_hi, h_hi, h_lo, h_lo],
        [v_lo, v_lo, v_hi, v_hi, v_lo],
        linestyle=HEATMAP_INSET_BBOX_LINESTYLE,
        color=HEATMAP_INSET_BBOX_COLOR,
        linewidth=HEATMAP_INSET_BBOX_LINEWIDTH,
        zorder=6,
    )

    # ---- 2. 创建 inset axes（相对 parent_ax 子图坐标）----
    inset_ax = parent_ax.inset_axes(HEATMAP_INSET_LOC)

    # 裁剪数据到 inset 范围
    h_centers = 0.5 * (h_edges[:-1] + h_edges[1:])
    v_centers = 0.5 * (v_edges[:-1] + v_edges[1:])
    h_idx = np.where((h_centers >= h_lo) & (h_centers <= h_hi))[0]
    v_idx = np.where((v_centers >= v_lo) & (v_centers <= v_hi))[0]
    if len(h_idx) < 2 or len(v_idx) < 2:
        inset_ax.set_visible(False)
        return

    inset_h_edges = h_edges[h_idx[0]: h_idx[-1] + 2]
    inset_v_edges = v_edges[v_idx[0]: v_idx[-1] + 2]
    inset_data = masked_score_2d[np.ix_(h_idx, v_idx)]

    # inset norm
    if HEATMAP_INSET_SHARE_NORM and norm is not None:
        inset_norm = norm
    else:
        inset_norm = build_heatmap_norm(inset_data)

    # ---- 3. inset 上画 pcolormesh ----
    inset_mesh = inset_ax.pcolormesh(
        inset_h_edges,
        inset_v_edges,
        inset_data.T,
        shading="flat",
        cmap=HEATMAP_COLORMAP,
        norm=inset_norm,
        edgecolors="none",
        linewidth=0.0,
        rasterized=True,
    )

    inset_ax.set_xlim(float(inset_h_edges[0]), float(inset_h_edges[-1]))
    inset_ax.set_ylim(float(inset_v_edges[0]), float(inset_v_edges[-1]))
    inset_ax.set_aspect("equal", adjustable="box")
    inset_ax.set_xlabel("")
    inset_ax.set_ylabel("")
    inset_ax.tick_params(direction="in", top=True, right=True, labelsize=7,
                         length=3, pad=2)

    # inset 边框
    for spine in inset_ax.spines.values():
        spine.set_edgecolor(HEATMAP_INSET_EDGE_COLOR)
        spine.set_linewidth(HEATMAP_INSET_EDGE_LINEWIDTH)

    # 无标题

    # ---- 4. inset 上叠加几何边界 ----
    if HEATMAP_INSET_SHOW_GEOMETRY:
        _draw_geometry_on_heatmap(
            inset_ax,
            dimension,
            requested_pos_cm,
            inset_h_edges,
            inset_v_edges,
            None,  # 不导出 CSV
        )

    # ---- 5. 连接线（单线：虚线框右上角 -> inset 底部中心）----
    if HEATMAP_INSET_CONNECT_LINES:
        inset_left, inset_bottom, inset_w, inset_h = HEATMAP_INSET_LOC
        # inset 底部中心（axes 坐标，相对 parent_ax）
        inset_corner_axes = (inset_left + inset_w / 2, inset_bottom)
        # 虚线框右上角（data 坐标 -> axes 坐标）
        bbox_corner_data = (h_hi, v_hi)
        # 用 blended transform: x 用 data, y 用 data
        from matplotlib.transforms import blended_transform_factory
        transform = blended_transform_factory(parent_ax.transData, parent_ax.transData)
        # 在 parent_ax 上画线：从 data 点到 axes 点
        # 用 ConnectionPatch 更简单
        from matplotlib.patches import ConnectionPatch
        con = ConnectionPatch(
            xyA=bbox_corner_data, coordsA=parent_ax.transData,
            xyB=inset_corner_axes, coordsB=parent_ax.transAxes,
            linestyle=HEATMAP_INSET_CONNECT_LINESTYLE,
            color=HEATMAP_INSET_CONNECT_COLOR,
            linewidth=HEATMAP_INSET_CONNECT_LINEWIDTH,
            zorder=10,
            clip_on=False,
        )
        parent_ax.add_artist(con)


def should_generate_heatmap_for_case(
    case_index: int,
) -> bool:
    """判断当前 CASE 是否需要输出截面热图。"""
    if not GENERATE_XY_HEATMAPS:
        return False

    if HEATMAP_CASE_INDICES is None:
        return True

    return case_index in HEATMAP_CASE_INDICES


def get_heatmap_positions(
    dimension: str,
) -> list[float]:
    """返回指定维度对应的固定轴坐标列表。"""
    fixed_axis = (
        HEATMAP_DIMENSION_INFO[dimension][
            "fixed_axis"
        ]
    )

    if fixed_axis == "Z":
        return [
            float(z)
            for z in HEATMAP_Z_POSITIONS_CM
        ]

    if fixed_axis == "Y":
        return [
            float(y)
            for y in HEATMAP_Y_POSITIONS_CM
        ]

    return [
        float(x)
        for x in HEATMAP_X_POSITIONS_CM
    ]


def select_heatmap_bin(
    axes: dict,
    axis_name: str,
    requested_position_cm: float,
) -> dict:
    """根据实际坐标选择指定轴的 bin。"""
    axis_edges = get_axis_edges(
        axes[axis_name]
    )

    bin_centers = (
        axis_edges[:-1]
        + axis_edges[1:]
    ) / 2.0

    requested_position_cm = float(
        requested_position_cm
    )

    if HEATMAP_Z_SELECTION_MODE == "nearest_center":
        bin_index = int(
            np.argmin(
                np.abs(
                    bin_centers
                    - requested_position_cm
                )
            )
        )

    else:
        tolerance = max(
            1.0,
            abs(float(axis_edges[0])),
            abs(float(axis_edges[-1])),
        ) * 1.0e-12

        if (
            requested_position_cm
            < axis_edges[0] - tolerance
            or requested_position_cm
            > axis_edges[-1] + tolerance
        ):
            raise ValueError(
                f"请求的热图位置 "
                f"{axis_name}="
                f"{requested_position_cm:g} cm "
                f"超出 USRBIN {axis_name} 网格范围 "
                f"[{axis_edges[0]:g}, "
                f"{axis_edges[-1]:g}] cm。"
            )

        # 处理浮点误差，并使最右端点归入最后一个 bin。
        clipped_pos = min(
            max(
                requested_position_cm,
                float(axis_edges[0]),
            ),
            float(axis_edges[-1]),
        )

        if np.isclose(
            clipped_pos,
            axis_edges[-1],
            rtol=0.0,
            atol=tolerance,
        ):
            bin_index = len(bin_centers) - 1
        else:
            bin_index = int(
                np.searchsorted(
                    axis_edges,
                    clipped_pos,
                    side="right",
                ) - 1
            )

    return {
        "index": bin_index,
        "requested_cm": (
            requested_position_cm
        ),
        "min_cm": float(
            axis_edges[bin_index]
        ),
        "max_cm": float(
            axis_edges[bin_index + 1]
        ),
        "center_cm": float(
            bin_centers[bin_index]
        ),
    }


def prepare_2d_heatmap_data(
    score_3d: np.ndarray,
    error_percent_3d: np.ndarray,
    dimension: str,
    bin_index: int,
    effective_scale: float,
) -> dict:
    """提取一个固定轴 bin 的二维截面数据，并应用缩放和统计屏蔽。"""
    if dimension == "xy":
        raw_score = score_3d[
            :, :, bin_index
        ]
        error_percent = (
            error_percent_3d[
                :, :, bin_index
            ]
        )
    elif dimension == "xz":
        raw_score = score_3d[
            :, bin_index, :
        ]
        error_percent = (
            error_percent_3d[
                :, bin_index, :
            ]
        )
    else:  # "yz"
        # free_axes=("Z","Y")，切片 score_3d[bin,:,:]
        # 原始排列为 (Y, Z)，转置为 (Z, Y) 与 free_axes 一致。
        raw_score = score_3d[
            bin_index, :, :
        ].T
        error_percent = (
            error_percent_3d[
                bin_index, :, :
            ].T
        )

    raw_score = (
        raw_score
        .copy()
        .astype(float)
    )

    error_percent = (
        error_percent
        .copy()
        .astype(float)
    )

    heatmap_scale = (
        float(effective_scale)
        if HEATMAP_USE_SCALED_SCORE
        else 1.0
    )

    plotted_score = (
        raw_score
        * heatmap_scale
    )

    invalid_mask = (
        ~np.isfinite(plotted_score)
    )

    if HEATMAP_MASK_INSUFFICIENT_STATISTICS:
        invalid_mask |= (
            np.isfinite(error_percent)
            & (
                error_percent
                >= INSUFFICIENT_ERROR_PERCENT
            )
        )

    if HEATMAP_USE_LOG_COLOR:
        invalid_mask |= (
            plotted_score <= 0
        )

    masked_score = np.ma.array(
        plotted_score,
        mask=invalid_mask,
        copy=False,
    )

    if masked_score.count() == 0:
        raise ValueError(
            "所选截面在应用统计筛选和色标条件后"
            "没有可绘制数据。"
        )

    return {
        "raw_score": raw_score,
        "plotted_score": plotted_score,
        "masked_score": masked_score,
        "error_percent": error_percent,
        "heatmap_scale": heatmap_scale,
    }


def build_heatmap_norm(
    masked_score: np.ma.MaskedArray,
):
    """根据热图数据和用户配置构造线性或对数颜色归一化。"""
    finite_values = np.asarray(
        masked_score.compressed(),
        dtype=float,
    )

    if finite_values.size == 0:
        raise ValueError(
            "热图没有可用于确定颜色范围的有效数据。"
        )

    if HEATMAP_USE_LOG_COLOR:
        positive_values = finite_values[
            finite_values > 0
        ]

        if positive_values.size == 0:
            raise ValueError(
                "对数热图要求至少存在一个大于 0 的数据。"
            )

        vmin = (
            float(np.min(positive_values))
            if HEATMAP_VMIN is None
            else float(HEATMAP_VMIN)
        )

        vmax = (
            float(np.max(positive_values))
            if HEATMAP_VMAX is None
            else float(HEATMAP_VMAX)
        )

        if vmin <= 0 or vmax <= 0:
            raise ValueError(
                "对数热图的颜色上下限必须大于 0。"
            )

        if vmax <= vmin:
            if np.isclose(vmax, vmin):
                vmax = vmin * (1.0 + 1.0e-6)
            else:
                raise ValueError(
                    "热图颜色上限必须大于下限。"
                )

        return LogNorm(
            vmin=vmin,
            vmax=vmax,
        )

    vmin = (
        float(np.min(finite_values))
        if HEATMAP_VMIN is None
        else float(HEATMAP_VMIN)
    )

    vmax = (
        float(np.max(finite_values))
        if HEATMAP_VMAX is None
        else float(HEATMAP_VMAX)
    )

    if vmax <= vmin:
        if np.isclose(vmax, vmin):
            padding = max(
                abs(vmin) * 1.0e-6,
                1.0e-15,
            )
            vmin -= padding
            vmax += padding
        else:
            raise ValueError(
                "热图颜色上限必须大于下限。"
            )

    return Normalize(
        vmin=vmin,
        vmax=vmax,
    )


def format_coordinate_for_filename(
    value: float,
) -> str:
    """将坐标转换为不含小数点和正负号歧义的文件名片段。"""
    text_value = f"{float(value):+.8g}"

    return (
        text_value
        .replace("+", "p")
        .replace("-", "m")
        .replace(".", "p")
    )


def export_2d_heatmap_csv(
    output_path: Path,
    axes: dict,
    dimension: str,
    bin_info: dict,
    heatmap_data: dict,
) -> None:
    """输出每个 bin 的坐标、评分值和百分比误差。"""
    dim_info = (
        HEATMAP_DIMENSION_INFO[dimension]
    )

    h_axis, v_axis = dim_info["free_axes"]
    fixed_axis = dim_info["fixed_axis"]

    h_edges = get_axis_edges(
        axes[h_axis]
    )

    v_edges = get_axis_edges(
        axes[v_axis]
    )

    h_min_grid, v_min_grid = np.meshgrid(
        h_edges[:-1],
        v_edges[:-1],
        indexing="ij",
    )

    h_max_grid, v_max_grid = np.meshgrid(
        h_edges[1:],
        v_edges[1:],
        indexing="ij",
    )

    h_center_grid = (
        h_min_grid + h_max_grid
    ) / 2.0

    v_center_grid = (
        v_min_grid + v_max_grid
    ) / 2.0

    h_lower = h_axis.lower()
    v_lower = v_axis.lower()
    f_lower = fixed_axis.lower()

    table = np.column_stack([
        h_min_grid.ravel(order="C"),
        h_max_grid.ravel(order="C"),
        h_center_grid.ravel(order="C"),
        v_min_grid.ravel(order="C"),
        v_max_grid.ravel(order="C"),
        v_center_grid.ravel(order="C"),
        np.full(
            h_min_grid.size,
            bin_info["requested_cm"],
        ),
        np.full(
            h_min_grid.size,
            bin_info["index"],
        ),
        np.full(
            h_min_grid.size,
            bin_info["min_cm"],
        ),
        np.full(
            h_min_grid.size,
            bin_info["max_cm"],
        ),
        np.full(
            h_min_grid.size,
            bin_info["center_cm"],
        ),
        heatmap_data[
            "raw_score"
        ].ravel(order="C"),
        heatmap_data[
            "plotted_score"
        ].ravel(order="C"),
        heatmap_data[
            "error_percent"
        ].ravel(order="C"),
    ])

    np.savetxt(
        output_path,
        table,
        delimiter=",",
        header=(
            f"{h_lower}_min_cm,"
            f"{h_lower}_max_cm,"
            f"{h_lower}_center_cm,"
            f"{v_lower}_min_cm,"
            f"{v_lower}_max_cm,"
            f"{v_lower}_center_cm,"
            f"requested_{f_lower}_cm,"
            f"{f_lower}_bin_index_zero_based,"
            f"{f_lower}_min_cm,"
            f"{f_lower}_max_cm,"
            f"{f_lower}_center_cm,"
            "raw_score,plotted_score,"
            "error_percent"
        ),
        comments="",
        fmt="%.10e",
    )


def export_3d_score(
    output_dir: Path,
    case_name: str,
    axes: dict,
    score: np.ndarray,
    error_percent: np.ndarray,
    has_error_matrix: bool,
    effective_scale: float,
) -> None:
    """导出 3D score 全网格到 CSV + npz（专用文件夹）。

    axes: detector["axes"]，含 X/Y/Z 的 min/max/n/width
    score: shape (nx, ny, nz)，FLUKA 原始 per-primary 值
    error_percent: 同 shape，百分比误差（无误差矩阵时填 NaN）
    effective_scale: current_scale * user_y_scale，用于 scaled 列
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    x_edges = get_axis_edges(axes["X"])
    y_edges = get_axis_edges(axes["Y"])
    z_edges = get_axis_edges(axes["Z"])
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])

    nx, ny, nz = score.shape
    if (nx, ny, nz) != (len(x_centers), len(y_centers), len(z_centers)):
        raise ValueError(
            f"3D 导出维度不匹配: score={score.shape}, "
            f"centers=({len(x_centers)},{len(y_centers)},{len(z_centers)})"
        )

    scaled_score = score * float(effective_scale)

    # ---- npz：紧凑二进制 ----
    if EXPORT_3D_SCORE_NPZ:
        dtype = np.float32 if EXPORT_3D_SCORE_NPZ_FLOAT32 else np.float64
        npz_path = output_dir / f"{case_name}_3d_score.npz"
        # 轴边界/中心保留 float64（精度），score/error 按 dtype 压缩
        import json
        axes_meta = json.dumps(
            {k: {kk: float(vv) if kk != "n" else int(vv)
                 for kk, vv in v.items()}
             for k, v in axes.items()},
            ensure_ascii=False,
        )
        np.savez_compressed(
            npz_path,
            score_raw=score.astype(dtype),
            score_scaled=scaled_score.astype(dtype),
            error_percent=error_percent.astype(dtype),
            x_edges=x_edges,
            y_edges=y_edges,
            z_edges=z_edges,
            x_centers=x_centers,
            y_centers=y_centers,
            z_centers=z_centers,
            axes_meta=axes_meta,
            effective_scale=float(effective_scale),
            has_error_matrix=bool(has_error_matrix),
        )
        print(f"  [3D 导出] npz -> {npz_path.name} ({npz_path.stat().st_size / 1e6:.1f} MB)")

    # ---- CSV：长表，逐 bin 输出 ----
    if EXPORT_3D_SCORE_CSV:
        csv_path = output_dir / f"{case_name}_3d_score.csv"
        # 用 meshgrid 生成 3D 中心坐标（indexing="ij" 与 score 维度顺序一致）
        X, Y, Z = np.meshgrid(x_centers, y_centers, z_centers, indexing="ij")

        # 拼接列（C order ravel，与 meshgrid 默认一致）
        cols = [X.ravel(order="C"),
                Y.ravel(order="C"),
                Z.ravel(order="C"),
                score.ravel(order="C")]
        header_cols = ["x_center", "y_center", "z_center", "raw_score"]

        if EXPORT_3D_SCORE_CSV_SCALED:
            cols.append(scaled_score.ravel(order="C"))
            header_cols.append("scaled_score")

        if EXPORT_3D_SCORE_CSV_ERROR:
            cols.append(error_percent.ravel(order="C"))
            header_cols.append("error_percent")

        table = np.column_stack(cols)
        np.savetxt(
            csv_path,
            table,
            delimiter=",",
            header=",".join(header_cols),
            comments="",
            fmt="%.10e",
        )
        print(f"  [3D 导出] csv -> {csv_path.name} ({csv_path.stat().st_size / 1e6:.1f} MB)")


def generate_heatmaps_for_case(
    case_index: int,
    path: Path,
    label: str,
    detector: dict,
    effective_scale: float,
) -> list[Path]:
    """为一个工况生成用户指定维度的全部热图。"""
    if not should_generate_heatmap_for_case(
        case_index
    ):
        return []

    HEATMAP_OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    axes = detector["axes"]

    clean_case_name = make_clean_csv_name(
        path,
        case_index,
    )

    output_paths: list[Path] = []

    for dimension in HEATMAP_DIMENSIONS:
        dim_info = (
            HEATMAP_DIMENSION_INFO[
                dimension
            ]
        )

        h_axis, v_axis = (
            dim_info["free_axes"]
        )

        axis_to_name = {"X": "x", "Y": "y", "Z": "z"}

        fixed_axis = (
            dim_info["fixed_axis"]
        )

        h_edges = get_axis_edges(
            axes[h_axis]
        )

        v_edges = get_axis_edges(
            axes[v_axis]
        )

        positions = get_heatmap_positions(
            dimension
        )

        for requested_pos_cm in positions:
            bin_info = select_heatmap_bin(
                axes,
                fixed_axis,
                float(
                    requested_pos_cm
                ),
            )

            heatmap_data = (
                prepare_2d_heatmap_data(
                    detector["score"],
                    detector[
                        "error_percent"
                    ],
                    dimension,
                    bin_info["index"],
                    effective_scale,
                )
            )

            norm = build_heatmap_norm(
                heatmap_data[
                    "masked_score"
                ]
            )

            figure, heatmap_ax = (
                plt.subplots(
                    figsize=(
                        HEATMAP_FIGSIZE
                    )
                )
            )

            if HEATMAP_SHOW_BIN_EDGES:
                edgecolors = "black"
                linewidth = 0.15
            else:
                edgecolors = "none"
                linewidth = 0.0

            mesh = heatmap_ax.pcolormesh(
                h_edges,
                v_edges,
                heatmap_data[
                    "masked_score"
                ].T,
                shading="flat",
                cmap=HEATMAP_COLORMAP,
                norm=norm,
                edgecolors=edgecolors,
                linewidth=linewidth,
                rasterized=True,
            )

            colorbar = figure.colorbar(
                mesh,
                ax=heatmap_ax,
                pad=0.025,
            )

            colorbar.set_label(
                HEATMAP_COLORBAR_LABEL,
                fontsize=13,
                labelpad=10,
            )

            heatmap_ax.set_xlabel(
                dim_info["xlabel"],
                fontsize=14,
                labelpad=8,
            )

            heatmap_ax.set_ylabel(
                dim_info["ylabel"],
                fontsize=14,
                labelpad=8,
            )

            if HEATMAP_EQUAL_ASPECT:
                heatmap_ax.set_aspect(
                    "equal",
                    adjustable="box",
                )

            heatmap_ax.set_xlim(
                float(h_edges[0]),
                float(h_edges[-1]),
            )

            heatmap_ax.set_ylim(
                float(v_edges[0]),
                float(v_edges[-1]),
            )

            if HEATMAP_TITLE is None:
                detector_name = (
                    detector["name"]
                    or "unnamed"
                )

                title = (
                    f"{label}\n"
                    f"Detector "
                    f"{detector['number']} "
                    f"({detector_name}), "
                    f"requested "
                    f"{fixed_axis}="
                    f"{bin_info['requested_cm']:.6g}"
                    " cm; "
                    f"bin {bin_info['index']} "
                    f"[{bin_info['min_cm']:.6g}, "
                    f"{bin_info['max_cm']:.6g}] cm"
                )
            else:
                title = HEATMAP_TITLE

            heatmap_ax.set_title(
                title,
                fontsize=13,
                pad=12,
            )

            heatmap_ax.tick_params(
                direction="in",
                top=True,
                right=True,
            )

            # 叠加 inset 放大图
            _draw_inset_on_heatmap(
                heatmap_ax,
                figure,
                mesh,
                norm,
                heatmap_data["masked_score"],
                h_edges,
                v_edges,
                axis_to_name[h_axis],
                axis_to_name[v_axis],
                dimension,
                float(requested_pos_cm),
                dim_info["xlabel"],
                dim_info["ylabel"],
            )

            figure.tight_layout()

            pos_name = (
                format_coordinate_for_filename(
                    bin_info[
                        "requested_cm"
                    ]
                )
            )

            output_stem = (
                HEATMAP_OUTPUT_DIR
                / (
                    f"{clean_case_name}_"
                    f"{dimension}_"
                    f"{fixed_axis.lower()}_"
                    f"{pos_name}_bin_"
                    f"{bin_info['index']:04d}"
                )
            )

            # 叠加 FLUKA 几何边界线
            _draw_geometry_on_heatmap(
                heatmap_ax,
                dimension,
                float(requested_pos_cm),
                h_edges,
                v_edges,
                output_stem,
            )

            for output_format in (
                HEATMAP_OUTPUT_FORMATS
            ):
                output_format = str(
                    output_format
                ).lower()

                image_path = (
                    output_stem
                    .with_suffix(
                        f".{output_format}"
                    )
                )

                save_kwargs = {
                    "bbox_inches": (
                        "tight"
                    ),
                }

                if (
                    output_format
                    == "png"
                ):
                    save_kwargs["dpi"] = (
                        HEATMAP_DPI
                    )

                figure.savefig(
                    image_path,
                    **save_kwargs,
                )

                output_paths.append(
                    image_path
                )

            if HEATMAP_EXPORT_CSV:
                csv_path = (
                    output_stem
                    .with_suffix(".csv")
                )

                export_2d_heatmap_csv(
                    csv_path,
                    axes,
                    dimension,
                    bin_info,
                    heatmap_data,
                )

                output_paths.append(
                    csv_path
                )

            print(
                f"  {dimension} heatmap: "
                f"requested {fixed_axis}="
                f"{bin_info['requested_cm']:.8e}"
                " cm, "
                f"selected bin="
                f"{bin_info['index']}, "
                f"range=["
                f"{bin_info['min_cm']:.8e}, "
                f"{bin_info['max_cm']:.8e}"
                "] cm, "
                f"center="
                f"{bin_info['center_cm']:.8e}"
                " cm"
            )

    return output_paths


def generate_merged_heatmaps(
    case_infos: list[dict],
) -> list[Path]:
    """将多个 CASE 的热图合并到同一张图中。"""
    if not case_infos:
        return []

    HEATMAP_OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_paths: list[Path] = []

    for dimension in HEATMAP_DIMENSIONS:
        dim_info = (
            HEATMAP_DIMENSION_INFO[
                dimension
            ]
        )

        h_axis, v_axis = (
            dim_info["free_axes"]
        )

        fixed_axis = (
            dim_info["fixed_axis"]
        )

        positions = get_heatmap_positions(
            dimension
        )

        for requested_pos_cm in positions:
            # 收集所有 CASE 在该维度和位置
            # 的截面数据。
            slice_data: list[dict] = []

            for info in case_infos:
                detector = info[
                    "detector"
                ]

                axes = detector["axes"]

                try:
                    bin_info = (
                        select_heatmap_bin(
                            axes,
                            fixed_axis,
                            float(
                                requested_pos_cm
                            ),
                        )
                    )
                except ValueError:
                    continue

                heatmap_data = (
                    prepare_2d_heatmap_data(
                        detector[
                            "score"
                        ],
                        detector[
                            "error_percent"
                        ],
                        dimension,
                        bin_info["index"],
                        info[
                            "effective_scale"
                        ],
                    )
                )

                slice_data.append({
                    "case_index": info[
                        "case_index"
                    ],
                    "label": info[
                        "label"
                    ],
                    "path": info[
                        "path"
                    ],
                    "axes": axes,
                    "bin_info": bin_info,
                    "heatmap_data": (
                        heatmap_data
                    ),
                    "detector_name": (
                        detector["name"]
                        or "unnamed"
                    ),
                    "detector_number": (
                        detector[
                            "number"
                        ]
                    ),
                })

            if not slice_data:
                continue

            # 用全部 CASE 的有效数据确定
            # 统一颜色范围。
            all_values = np.concatenate([
                np.asarray(
                    d[
                        "heatmap_data"
                    ][
                        "masked_score"
                    ].compressed(),
                    dtype=float,
                )
                for d in slice_data
            ])

            all_masked = np.ma.array(
                all_values,
                mask=(
                    ~np.isfinite(
                        all_values
                    )
                ),
            )

            norm = build_heatmap_norm(
                all_masked
            )

            n_cases = len(slice_data)
            ncols = min(
                HEATMAP_MERGE_NCOLS,
                n_cases,
            )

            nrows = (
                n_cases + ncols - 1
            ) // ncols

            figure, axes_arr = (
                plt.subplots(
                    nrows,
                    ncols,
                    figsize=(
                        HEATMAP_FIGSIZE[0]
                        * ncols,
                        HEATMAP_FIGSIZE[1]
                        * nrows,
                    ),
                    squeeze=False,
                )
            )

            for idx, data in enumerate(
                slice_data
            ):
                row = idx // ncols
                col = idx % ncols
                ax = axes_arr[row][col]

                h_edges = (
                    get_axis_edges(
                        data["axes"][
                            h_axis
                        ]
                    )
                )

                v_edges = (
                    get_axis_edges(
                        data["axes"][
                            v_axis
                        ]
                    )
                )

                if HEATMAP_SHOW_BIN_EDGES:
                    edgecolors = (
                        "black"
                    )
                    linewidth = 0.15
                else:
                    edgecolors = "none"
                    linewidth = 0.0

                cmap = (
                    HEATMAP_COLORMAP
                )
                alpha = 1.0

                if (
                    HEATMAP_MERGE_DIFFERENTIATE
                ):
                    cmap = (
                        HEATMAP_MERGE_COLORMAPS[
                            idx
                            % len(
                                HEATMAP_MERGE_COLORMAPS
                            )
                        ]
                    )

                    alpha = (
                        HEATMAP_MERGE_ALPHAS[
                            idx
                            % len(
                                HEATMAP_MERGE_ALPHAS
                            )
                        ]
                    )

                mesh = ax.pcolormesh(
                    h_edges,
                    v_edges,
                    data[
                        "heatmap_data"
                    ][
                        "masked_score"
                    ].T,
                    shading="flat",
                    cmap=cmap,
                    norm=norm,
                    alpha=alpha,
                    edgecolors=edgecolors,
                    linewidth=linewidth,
                    rasterized=True,
                )

                # 记录第一个 mesh，循环结束后再创建色标，
                # 避免与 tight_layout 的子图重排冲突。
                if idx == 0:
                    first_mesh = mesh

                ax.set_xlabel(
                    dim_info[
                        "xlabel"
                    ],
                    fontsize=14,
                    labelpad=8,
                )

                ax.set_ylabel(
                    dim_info[
                        "ylabel"
                    ],
                    fontsize=14,
                    labelpad=8,
                )

                if HEATMAP_EQUAL_ASPECT:
                    ax.set_aspect(
                        "equal",
                        adjustable="box",
                    )

                ax.set_xlim(
                    float(
                        h_edges[0]
                    ),
                    float(
                        h_edges[-1]
                    ),
                )

                ax.set_ylim(
                    float(
                        v_edges[0]
                    ),
                    float(
                        v_edges[-1]
                    ),
                )

                bin_info = data[
                    "bin_info"
                ]

                if HEATMAP_TITLE is None:
                    title = (
                        f"{data['label']}\n"
                        f"requested "
                        f"{fixed_axis}="
                        f"{bin_info['requested_cm']:.6g}"
                        " cm; "
                        f"bin "
                        f"{bin_info['index']} "
                        f"["
                        f"{bin_info['min_cm']:.6g}"
                        ", "
                        f"{bin_info['max_cm']:.6g}"
                        "] cm"
                    )
                else:
                    title = (
                        HEATMAP_TITLE
                    )

                ax.set_title(
                    title,
                    fontsize=13,
                    pad=12,
                )

                ax.tick_params(
                    direction="in",
                    top=True,
                    right=True,
                )

                # 叠加 FLUKA 几何边界线
                _draw_geometry_on_heatmap(
                    ax,
                    dimension,
                    float(requested_pos_cm),
                    h_edges,
                    v_edges,
                    None,  # merged 模式 CSV 在下方统一导出
                )

            # 隐藏多余的子图。
            for idx in range(
                n_cases,
                nrows * ncols,
            ):
                row = idx // ncols
                col = idx % ncols
                axes_arr[row][col].set_visible(
                    False
                )

            # 叠加 inset 放大图
            if slice_data:
                axis_to_name_local = {"X": "x", "Y": "y", "Z": "z"}
                if HEATMAP_INSET_PER_SUBPLOT:
                    # 每个子图都画 inset
                    for idx, data in enumerate(slice_data):
                        row = idx // ncols
                        col = idx % ncols
                        ax_sub = axes_arr[row][col]
                        h_edges_sub = get_axis_edges(data["axes"][h_axis])
                        v_edges_sub = get_axis_edges(data["axes"][v_axis])
                        _draw_inset_on_heatmap(
                            ax_sub,
                            figure,
                            None,  # mesh 不单独传
                            norm,
                            data["heatmap_data"]["masked_score"],
                            h_edges_sub,
                            v_edges_sub,
                            axis_to_name_local[h_axis],
                            axis_to_name_local[v_axis],
                            dimension,
                            float(requested_pos_cm),
                            dim_info["xlabel"],
                            dim_info["ylabel"],
                        )
                else:
                    # 只在第一个子图画
                    first_data = slice_data[0]
                    first_h_edges = get_axis_edges(first_data["axes"][h_axis])
                    first_v_edges = get_axis_edges(first_data["axes"][v_axis])
                    _draw_inset_on_heatmap(
                        axes_arr[0][0],
                        figure,
                        first_mesh,
                        norm,
                        first_data["heatmap_data"]["masked_score"],
                        first_h_edges,
                        first_v_edges,
                        axis_to_name_local[h_axis],
                        axis_to_name_local[v_axis],
                        dimension,
                        float(requested_pos_cm),
                        dim_info["xlabel"],
                        dim_info["ylabel"],
                    )

            figure.tight_layout()

            # 子图布局确定后再创建色标，
            # 关联所有可见子图，使其从右侧统一预留空间。
            visible_axes = [
                ax
                for ax in axes_arr.ravel()
                if ax.get_visible()
            ]

            colorbar = figure.colorbar(
                first_mesh,
                ax=visible_axes,
                pad=0.02,
                shrink=0.95,
            )

            colorbar.set_label(
                HEATMAP_COLORBAR_LABEL,
                fontsize=13,
                labelpad=10,
            )

            pos_name = (
                format_coordinate_for_filename(
                    float(
                        requested_pos_cm
                    )
                )
            )

            output_stem = (
                HEATMAP_OUTPUT_DIR
                / (
                    f"merged_"
                    f"{dimension}_"
                    f"{fixed_axis.lower()}_"
                    f"{pos_name}"
                )
            )

            # merged 模式：统一导出一次材料分布 CSV
            if FLUKA_GEOMETRY_ENABLED and FLUKA_GEOMETRY_EXPORT_CSV:
                try:
                    from fluka_geometry import cut_plane as _fluka_cut_plane
                    geom = _get_fluka_geometry()
                    if geom is not None:
                        cut_axis, cut_value = _get_geometry_cut_value(
                            dimension, float(requested_pos_cm)
                        )
                        axis_to_name = {"X": "x", "Y": "y", "Z": "z"}
                        h_name = axis_to_name[h_axis]
                        v_name = axis_to_name[v_axis]
                        cut_free1_name, cut_free2_name = {
                            "x": ("y", "z"),
                            "y": ("x", "z"),
                            "z": ("x", "y"),
                        }[cut_axis]
                        h_edges_0 = get_axis_edges(slice_data[0]["axes"][h_axis])
                        v_edges_0 = get_axis_edges(slice_data[0]["axes"][v_axis])
                        edges_map = {h_name: h_edges_0, v_name: v_edges_0}
                        mg = _fluka_cut_plane(
                            geom, cut_axis, cut_value,
                            edges_map[cut_free1_name],
                            edges_map[cut_free2_name],
                        )
                        csv_name = f"{output_stem.name}_materials_{cut_axis}{cut_value:.4f}.csv"
                        mg.to_csv(output_stem.with_name(csv_name))
                except Exception as e:
                    print(f"  [FLUKA几何] merged CSV 导出失败: {e}")

            for output_format in (
                HEATMAP_OUTPUT_FORMATS
            ):
                output_format = str(
                    output_format
                ).lower()

                image_path = (
                    output_stem
                    .with_suffix(
                        f".{output_format}"
                    )
                )

                save_kwargs = {
                    "bbox_inches": (
                        "tight"
                    ),
                }

                if (
                    output_format
                    == "png"
                ):
                    save_kwargs["dpi"] = (
                        HEATMAP_DPI
                    )

                figure.savefig(
                    image_path,
                    **save_kwargs,
                )

                output_paths.append(
                    image_path
                )

            print(
                f"  merged {dimension}"
                f" heatmap: {n_cases}"
                " cases, requested "
                f"{fixed_axis}="
                f"{float(requested_pos_cm):.8e}"
                " cm"
            )

    return output_paths


# ============================================================
# 11. 材料区域背景
# ============================================================
def add_material_regions(
    ax,
) -> None:
    """添加原程序中的 Z 向材料分区和标注。"""
    if not SHOW_MATERIAL_REGIONS:
        return

    for (
        left,
        right,
        _,
        color,
    ) in MATERIAL_REGIONS:
        if (
            left is None
            or right is None
            or right <= left
        ):
            continue

        ax.axvspan(
            left,
            right,
            color=color,
            alpha=0.78,
            zorder=0,
        )

    for boundary in MATERIAL_BOUNDARIES:
        ax.axvline(
            boundary,
            color="0.45",
            linestyle="--",
            linewidth=1.0,
            alpha=0.85,
            zorder=1,
        )

    transform = blended_transform_factory(
        ax.transData,
        ax.transAxes,
    )

    plot_left = (
        float(PLOT_Z_MIN)
        if PLOT_Z_MIN is not None
        else float(
            ax.get_xlim()[0]
        )
    )

    plot_right = (
        float(PLOT_Z_MAX)
        if PLOT_Z_MAX is not None
        else float(
            ax.get_xlim()[1]
        )
    )

    # 根据 MATERIAL_REGIONS 动态生成标签。
    # 宽区域用居中横排文字，窄区域用引线标注或竖排文字。
    # 阈值：宽度低于此值的区域视为"窄区域"。
    thin_threshold = 0.5  # cm

    for (
        left,
        right,
        name,
        _color,
    ) in MATERIAL_REGIONS:
        if (
            left is None
            or right is None
            or right <= left
        ):
            continue

        width = right - left
        center = (left + right) / 2.0

        if width >= thin_threshold:
            # 宽区域：居中横排文字
            fontsize = 12 if width >= 1.0 else 11
            ax.text(
                center,
                0.965,
                name,
                transform=transform,
                ha="center",
                va="top",
                fontsize=fontsize,
                zorder=5,
            )
        elif USE_THIN_REGION_CALLOUTS:
            # 窄区域 + 引线模式：用 annotate 指向区域中心
            ax.annotate(
                name,
                xy=(center, 0.98),
                xycoords=(
                    "data",
                    "axes fraction",
                ),
                xytext=(center, 1.08),
                textcoords=(
                    "data",
                    "axes fraction",
                ),
                ha="center",
                va="bottom",
                fontsize=10,
                annotation_clip=False,
                arrowprops={
                    "arrowstyle": "-",
                    "color": "0.30",
                    "linewidth": 0.9,
                },
            )
        else:
            # 窄区域 + 竖排模式：在区域中心放竖排文字
            ax.text(
                center,
                0.965,
                name,
                transform=transform,
                ha="center",
                va="top",
                rotation=90,
                fontsize=7,
                zorder=5,
            )


# ============================================================
# 12. 主程序
# ============================================================
def main() -> None:
    validate_configuration()
    validate_case_filenames()

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": [
            "Times New Roman",
            "STIXGeneral",
            "DejaVu Serif",
            "SimSun",
            "Noto Serif CJK SC",
        ],
        "mathtext.fontset": "stix",
        "font.size": 12,
        "axes.linewidth": 1.2,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
    })

    OUTPUT_STEM.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fig, ax = plt.subplots(
        figsize=(12.0, 7.0)
    )

    add_material_regions(
        ax
    )

    color_cycle = (
        plt.rcParams[
            "axes.prop_cycle"
        ]
        .by_key()
        .get(
            "color",
            ["C0"],
        )
    )

    csv_columns: list[
        np.ndarray
    ] = []

    csv_names: list[
        str
    ] = []

    total_results: list[
        dict
    ] = []

    heatmap_output_paths: list[
        Path
    ] = []

    heatmap_case_infos: list[
        dict
    ] = []

    reference_axes = None
    reference_z_min = None
    reference_z_max = None
    reference_z_center = None

    global_ymin = np.inf
    global_ymax = -np.inf

    reduction_mode = get_reduction_mode()

    print()
    print(
        "USRBIN longitudinal-profile processing"
    )
    print("-" * 78)

    print(
        "Transverse reduction mode "
        f"= {reduction_mode}"
    )

    if APPLY_CURRENT_SCALING:
        print(
            "Current scaling is enabled: "
            "current_mA / "
            f"{REFERENCE_CURRENT_MA:g} mA"
        )

    if SHOW_ERROR_BARS:
        print(
            "Statistical error bars are enabled; "
            "3D bin errors are propagated "
            "to the Z profile."
        )

    print("-" * 78)

    for case_index, (
        filename,
        label,
        current_ma,
        user_y_scale,
    ) in enumerate(
        CASES,
        start=1,
    ):
        path = resolve_case_path(
            filename
        )

        detector = read_usrbin_lis(
            path
        )

        axes = detector[
            "axes"
        ]

        profile_data = build_longitudinal_profile(
            axes,
            detector["score"],
            detector[
                "error_percent"
            ],
        )

        z_min = profile_data[
            "z_min"
        ]

        z_max = profile_data[
            "z_max"
        ]

        z_center = profile_data[
            "z_center"
        ]

        raw_profile = profile_data[
            "raw_profile"
        ]

        raw_absolute_error = profile_data[
            "raw_absolute_error"
        ]

        error_percent = profile_data[
            "profile_error_percent"
        ]

        if APPLY_CURRENT_SCALING:
            current_scale = (
                float(current_ma)
                / REFERENCE_CURRENT_MA
            )

        else:
            current_scale = 1.0

        effective_scale = (
            current_scale
            * float(user_y_scale)
        )

        plotted_profile = (
            raw_profile
            * effective_scale
        )

        plotted_absolute_error = (
            raw_absolute_error
            * abs(effective_scale)
        )

        # ---- 3D 能量沉积导出（独立开关，不受 heatmap case 筛选影响）----
        if EXPORT_3D_SCORE_ENABLED:
            export_3d_score(
                output_dir=EXPORT_3D_SCORE_DIR,
                case_name=make_clean_csv_name(path, case_index),
                axes=axes,
                score=detector["score"],
                error_percent=detector["error_percent"],
                has_error_matrix=detector.get("has_error_matrix", False),
                effective_scale=effective_scale,
            )

        if should_generate_heatmap_for_case(
            case_index
        ):
            heatmap_case_infos.append({
                "case_index": (
                    case_index
                ),
                "label": label,
                "path": path,
                "detector": detector,
                "effective_scale": (
                    effective_scale
                ),
            })

        # ----------------------------------------------------
        # 检查空间网格
        # ----------------------------------------------------
        if reference_axes is None:
            reference_axes = axes

            reference_z_min = (
                z_min.copy()
            )

            reference_z_max = (
                z_max.copy()
            )

            reference_z_center = (
                z_center.copy()
            )

            csv_columns.extend([
                reference_z_min,
                reference_z_max,
                reference_z_center,
            ])

            csv_names.extend([
                "z_min_cm",
                "z_max_cm",
                "z_center_cm",
            ])

        else:
            if not axes_are_equal(
                axes,
                reference_axes,
                include_xy=(
                    REQUIRE_IDENTICAL_XY_GRID
                ),
            ):
                checked_axes = (
                    "X、Y、Z"
                    if REQUIRE_IDENTICAL_XY_GRID
                    else "Z"
                )

                raise ValueError(
                    f"{path.name} 的 "
                    f"{checked_axes} 网格"
                    "与第一个文件不同。"
                )

        # ----------------------------------------------------
        # 图例名称
        # ----------------------------------------------------
        display_label = label

        if (
            SHOW_USER_SCALE_IN_LEGEND
            and not np.isclose(
                float(user_y_scale),
                1.0,
            )
        ):
            display_label += (
                rf" ($\times "
                rf"{float(user_y_scale):g}$)"
            )

        # ----------------------------------------------------
        # 自动 Y 轴范围
        # ----------------------------------------------------
        (
            local_ymin,
            local_ymax,
        ) = get_visible_y_limits(
            z_min,
            z_max,
            plotted_profile,
            error_percent,
        )

        if np.isfinite(
            local_ymin
        ):
            global_ymin = min(
                global_ymin,
                local_ymin,
            )

        if np.isfinite(
            local_ymax
        ):
            global_ymax = max(
                global_ymax,
                local_ymax,
            )

        # ----------------------------------------------------
        # 绘图
        # ----------------------------------------------------
        curve_color = color_cycle[
            (case_index - 1)
            % len(color_cycle)
        ]

        plot_one_profile(
            ax=ax,
            z_min=z_min,
            z_max=z_max,
            z_center=z_center,
            profile=plotted_profile,
            error_percent=error_percent,
            label=display_label,
            color=curve_color,
        )

        # ----------------------------------------------------
        # CSV 输出
        # ----------------------------------------------------
        clean_name = make_clean_csv_name(
            path,
            case_index,
        )

        if EXPORT_RAW_AND_SCALED:
            csv_columns.append(
                raw_profile
            )

            csv_names.append(
                clean_name
                + "_raw"
            )

            csv_columns.append(
                plotted_profile
            )

            csv_names.append(
                clean_name
                + "_scaled"
            )

        else:
            csv_columns.append(
                plotted_profile
            )

            csv_names.append(
                clean_name
                + "_plotted"
            )

        if EXPORT_RELATIVE_ERROR:
            csv_columns.append(
                error_percent
            )

            csv_names.append(
                clean_name
                + "_relative_error_percent"
            )

        if EXPORT_ABSOLUTE_ERROR:
            csv_columns.append(
                plotted_absolute_error
            )

            csv_names.append(
                clean_name
                + "_plotted_absolute_error"
            )

        # ----------------------------------------------------
        # 积分或求和
        # ----------------------------------------------------
        total_value = None
        total_z_min = None
        total_z_max = None

        if CALCULATE_TOTAL:
            if TOTAL_USE_SCALED_PROFILE:
                profile_for_total = (
                    plotted_profile
                )

                total_source = "scaled"

            else:
                profile_for_total = (
                    raw_profile
                )

                total_source = "raw"

            (
                total_value,
                total_z_min,
                total_z_max,
            ) = calculate_profile_total(
                z_min,
                z_max,
                profile_for_total,
            )

            total_results.append({
                "file": path.name,
                "label": label,
                "detector_number": (
                    detector["number"]
                ),
                "detector_name": (
                    detector["name"]
                ),
                "particle_number": (
                    detector[
                        "particle_number"
                    ]
                ),
                "track_length_binning": (
                    detector[
                        "track_length_binning"
                    ]
                ),
                "reduction_mode": (
                    reduction_mode
                ),
                "current_mA": (
                    current_ma
                ),
                "current_scale": (
                    current_scale
                ),
                "user_y_scale": (
                    user_y_scale
                ),
                "effective_scale": (
                    effective_scale
                ),
                "total_mode": (
                    TOTAL_MODE
                ),
                "total_source": (
                    total_source
                ),
                "z_min_cm": (
                    total_z_min
                ),
                "z_max_cm": (
                    total_z_max
                ),
                "total": (
                    total_value
                ),
                "unit": (
                    TOTAL_UNIT
                ),
            })

        # ----------------------------------------------------
        # 终端输出
        # ----------------------------------------------------
        finite_raw = raw_profile[
            np.isfinite(
                raw_profile
            )
        ]

        peak_raw_value = (
            float(
                np.max(
                    finite_raw
                )
            )
            if finite_raw.size > 0
            else np.nan
        )

        insufficient_bin_count = int(
            np.count_nonzero(
                np.isfinite(
                    error_percent
                )
                & (
                    error_percent
                    >= INSUFFICIENT_ERROR_PERCENT
                )
                & np.isfinite(
                    raw_profile
                )
                & (raw_profile != 0)
            )
        )

        print(
            path.name
        )

        print(
            "  detector = "
            f"{detector['number']}"
        )

        print(
            "  detector name = "
            f"{detector['name'] or '(unnamed)'}"
        )

        print(
            "  particle number = "
            f"{detector['particle_number']}"
        )

        print(
            "  track-length binning = "
            f"{detector['track_length_binning']}"
        )

        print(
            "  grid = "
            f"{axes['X']['n']} × "
            f"{axes['Y']['n']} × "
            f"{axes['Z']['n']}"
        )

        print(
            "  number of Z bins = "
            f"{z_center.size}"
        )

        print(
            "  X-Y reduction mode = "
            f"{reduction_mode}"
        )

        print(
            "  current scaling factor = "
            f"{current_scale:g}"
        )

        print(
            "  user Y scaling factor = "
            f"{float(user_y_scale):g}"
        )

        print(
            "  effective scaling factor = "
            f"{effective_scale:g}"
        )

        print(
            "  peak raw profile value = "
            f"{peak_raw_value:.8e}"
        )

        print(
            "  bins with insufficient statistics = "
            f"{insufficient_bin_count}"
        )

        if not detector[
            "has_error_matrix"
        ]:
            print(
                "  warning: percentage-error "
                "matrix was not found"
            )

        if CALCULATE_TOTAL:
            print(
                f"  {TOTAL_NAME} = "
                f"{total_value:.8e} "
                f"{TOTAL_UNIT}"
            )

            print(
                "  total Z range = "
                f"[{total_z_min:.8e}, "
                f"{total_z_max:.8e}] cm"
            )

        print("-" * 78)

    if (
        reference_z_min is None
        or reference_z_max is None
    ):
        raise ValueError(
            "没有读取到任何 USRBIN 数据。"
        )

    # ========================================================
    # 生成热图
    # ========================================================
    if heatmap_case_infos:
        if HEATMAP_MERGE_FILES:
            heatmap_output_paths = (
                generate_merged_heatmaps(
                    heatmap_case_infos
                )
            )
        else:
            for info in (
                heatmap_case_infos
            ):
                heatmap_output_paths.extend(
                    generate_heatmaps_for_case(
                        case_index=(
                            info[
                                "case_index"
                            ]
                        ),
                        path=info[
                            "path"
                        ],
                        label=info[
                            "label"
                        ],
                        detector=info[
                            "detector"
                        ],
                        effective_scale=(
                            info[
                                "effective_scale"
                            ]
                        ),
                    )
                )

    # ========================================================
    # 坐标轴设置
    # ========================================================
    if USE_LOG_X:
        if np.any(
            reference_z_min <= 0
        ):
            raise ValueError(
                "对数 X 轴要求所有显示的 "
                "Z 坐标大于 0。"
            )

        ax.set_xscale(
            "log"
        )

    else:
        ax.set_xscale(
            "linear"
        )

    x_min = (
        float(
            np.min(
                reference_z_min
            )
        )
        if PLOT_Z_MIN is None
        else float(
            PLOT_Z_MIN
        )
    )

    x_max = (
        float(
            np.max(
                reference_z_max
            )
        )
        if PLOT_Z_MAX is None
        else float(
            PLOT_Z_MAX
        )
    )

    if x_max <= x_min:
        raise ValueError(
            "PLOT_Z_MAX 必须大于 PLOT_Z_MIN。"
        )

    if (
        USE_LOG_X
        and x_min <= 0
    ):
        raise ValueError(
            "对数 X 轴的 PLOT_Z_MIN "
            "必须大于 0。"
        )

    ax.set_xlim(
        x_min,
        x_max,
    )

    if not np.isfinite(
        global_ymax
    ):
        raise ValueError(
            "实际绘图范围内没有有效 USRBIN 数据。"
        )

    if USE_LOG_Y:
        if (
            global_ymax <= 0
            or not np.isfinite(
                global_ymin
            )
        ):
            raise ValueError(
                "对数 Y 轴模式下没有找到大于 0 的数据。"
            )

        ax.set_yscale(
            "log",
            nonpositive="clip",
        )

        if USER_Y_MIN is None:
            y_min = (
                global_ymin
                * LOG_Y_BOTTOM_FACTOR
            )

        else:
            if USER_Y_MIN <= 0:
                raise ValueError(
                    "对数 Y 轴的 USER_Y_MIN "
                    "必须大于 0。"
                )

            y_min = float(
                USER_Y_MIN
            )

        y_max = (
            global_ymax
            * Y_HEADROOM_FACTOR
            if USER_Y_MAX is None
            else float(
                USER_Y_MAX
            )
        )

    else:
        ax.set_yscale(
            "linear"
        )

        if USER_Y_MIN is None:
            # 与原程序一致：全为非负时默认从 0 开始；
            # 若存在负值，则留出少量下边距。
            if global_ymin >= 0:
                y_min = 0.0

            else:
                y_min = (
                    global_ymin
                    * Y_HEADROOM_FACTOR
                )

        else:
            y_min = float(
                USER_Y_MIN
            )

        if USER_Y_MAX is None:
            if global_ymax > 0:
                y_max = (
                    global_ymax
                    * Y_HEADROOM_FACTOR
                )

            elif global_ymax < 0:
                y_max = (
                    global_ymax
                    / Y_HEADROOM_FACTOR
                )

            else:
                y_max = 1.0

        else:
            y_max = float(
                USER_Y_MAX
            )

    if y_max <= y_min:
        raise ValueError(
            "最终 Y 轴上限必须大于下限。"
        )

    ax.set_ylim(
        y_min,
        y_max,
    )

    ax.set_xlabel(
        X_AXIS_LABEL,
        fontsize=14,
        labelpad=10,
    )

    ax.set_ylabel(
        Y_AXIS_LABEL,
        fontsize=14,
        labelpad=8,
    )

    if PLOT_TITLE:
        ax.set_title(
            PLOT_TITLE,
            fontsize=15,
            pad=14,
        )

    # ========================================================
    # 网格与图例
    # ========================================================
    ax.grid(
        True,
        which="major",
        linewidth=0.6,
        alpha=0.25,
        zorder=0,
    )

    if USE_LOG_X or USE_LOG_Y:
        ax.grid(
            True,
            which="minor",
            linewidth=0.35,
            alpha=0.12,
            zorder=0,
        )

    (
        handles,
        legend_labels,
    ) = ax.get_legend_handles_labels()

    fig.legend(
        handles,
        legend_labels,
        loc="lower center",
        bbox_to_anchor=(
            0.5,
            0.018,
        ),
        ncol=3,
        frameon=False,
        fontsize=11,
        handlelength=3.0,
        columnspacing=1.8,
    )

    fig.subplots_adjust(
        left=0.115,
        right=0.985,
        top=(
            0.90
            if SHOW_MATERIAL_REGIONS
            else 0.94
        ),
        bottom=0.255,
    )

    # ========================================================
    # 保存图片
    # ========================================================
    png_path = OUTPUT_STEM.with_suffix(
        ".png"
    )

    pdf_path = OUTPUT_STEM.with_suffix(
        ".pdf"
    )

    svg_path = OUTPUT_STEM.with_suffix(
        ".svg"
    )

    fig.savefig(
        png_path,
        dpi=600,
        bbox_inches="tight",
    )

    fig.savefig(
        pdf_path,
        bbox_inches="tight",
    )

    fig.savefig(
        svg_path,
        bbox_inches="tight",
    )

    # ========================================================
    # 保存曲线 CSV
    # ========================================================
    curve_csv_path = (
        OUTPUT_STEM.with_suffix(
            ".csv"
        )
    )

    table = np.column_stack(
        csv_columns
    )

    np.savetxt(
        curve_csv_path,
        table,
        delimiter=",",
        header=",".join(
            csv_names
        ),
        comments="",
        fmt="%.10e",
    )

    # ========================================================
    # 保存积分结果 CSV
    # ========================================================
    total_csv_path = Path(
        str(OUTPUT_STEM)
        + "_totals.csv"
    )

    if CALCULATE_TOTAL:
        with total_csv_path.open(
            "w",
            newline="",
            encoding="utf-8-sig",
        ) as file:
            writer = csv.DictWriter(
                file,
                fieldnames=[
                    "file",
                    "label",
                    "detector_number",
                    "detector_name",
                    "particle_number",
                    "track_length_binning",
                    "reduction_mode",
                    "current_mA",
                    "current_scale",
                    "user_y_scale",
                    "effective_scale",
                    "total_mode",
                    "total_source",
                    "z_min_cm",
                    "z_max_cm",
                    "total",
                    "unit",
                ],
            )

            writer.writeheader()

            writer.writerows(
                total_results
            )

    # ========================================================
    # 最终终端输出
    # ========================================================
    print()
    print(
        "Output files:"
    )
    print(
        png_path
    )
    print(
        pdf_path
    )
    print(
        svg_path
    )
    print(
        curve_csv_path
    )

    if CALCULATE_TOTAL:
        print(
            total_csv_path
        )

    if heatmap_output_paths:
        print()
        print(
            "Heatmap output files:"
        )

        for heatmap_path in heatmap_output_paths:
            print(
                heatmap_path
            )

    plt.show()


if __name__ == "__main__":
    main()
