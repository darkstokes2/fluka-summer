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

from __future__ import annotations

import csv
import re
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np


# ============================================================
# 1. 输入文件
# ============================================================
BASE_DIR = Path(__file__).resolve().parent

# CASES每项格式：
#
# (
#     文件名,
#     图例名称,
#     束流电流/mA,
#     用户人工Y轴缩放倍率,
# )
#
# 示例：
#
# (
#     "your_spectrum_tab.lis",
#     r"$\mathrm{Case\ A}$",
#     1.0,
#     1.0,
# ),
#
# 请填入你自己的 *_tab.lis 文件名。
CASES = [
    # (
    #     "your_spectrum_tab.lis",
    #     r"$\mathrm{Case\ A}$",
    #     1.0,
    #     1.0,
    # ),
]


# ============================================================
# 2. 用户配置区
# ============================================================

# ------------------------------------------------------------
# A. USRBDX Detector选择
# ------------------------------------------------------------
# 一个tab.lis文件可能包含多个Detector。
#
# 按Detector编号选择：
DETECTOR_NUMBER = 1

# 可选：同时按照Detector名称中的关键词筛选。
#
# 例如：
#
# DETECTOR_NAME_CONTAINS = "piFluenUD"
#
# 不需要名称筛选时设置为None。
DETECTOR_NAME_CONTAINS = None


# ------------------------------------------------------------
# B. 束流电流缩放
# ------------------------------------------------------------
# False：
#     保留FLUKA原始per-primary结果。
#
# True：
#     按相对于REFERENCE_CURRENT_MA的电流比例缩放。
APPLY_CURRENT_SCALING = False

REFERENCE_CURRENT_MA = 1.0


# ------------------------------------------------------------
# C. 每根曲线的人工缩放倍率
# ------------------------------------------------------------
# 在CASES中最后一个数设置：
#
# 1.0：不缩放
# 2.0：放大2倍
# 0.2：缩小至20%
SHOW_USER_SCALE_IN_LEGEND = True


# ------------------------------------------------------------
# D. 曲线形式
# ------------------------------------------------------------
# True：
#     连接能量bin中心，绘制普通折线。
#
# False：
#     按能量bin边界绘制阶梯图。
#     对USRBDX分箱结果，更推荐使用阶梯图。
USE_CONTINUOUS_LINE = True


# ------------------------------------------------------------
# E. bin中心标记
# ------------------------------------------------------------
SHOW_BIN_MARKERS = True

BIN_MARKER = "o"
BIN_MARKER_SIZE = 3.0

# 每隔几个bin显示一个marker。
# 1表示每个bin都显示。
MARK_EVERY = 1


# ------------------------------------------------------------
# F. 阶梯图辅助连接线
# ------------------------------------------------------------
# 只有USE_CONTINUOUS_LINE=False时生效。
#
# 辅助线只连接bin中心，不进行平滑或插值，
# 且不会进入图例。
ADD_GUIDE_LINE = False

GUIDE_LINE_WIDTH = 0.75
GUIDE_LINE_ALPHA = 0.55


# ------------------------------------------------------------
# G. Error bar开关
# ------------------------------------------------------------
# True：
#     根据USRBDX第四列的相对统计误差绘制error bar。
#
# False：
#     不绘制error bar。
SHOW_ERROR_BARS = True

ERROR_BAR_LINEWIDTH = 0.7
ERROR_BAR_CAPSIZE = 2.0
ERROR_BAR_ALPHA = 0.60

# 每隔几个能量bin绘制一次error bar。
#
# 1：每个bin均显示
# 2：每隔一个bin显示
# 5：每隔四个bin显示
ERROR_BAR_EVERY = 1


# ------------------------------------------------------------
# H. 统计不足数据处理
# ------------------------------------------------------------
# FLUKA约定：
#
# 当统计量不足以计算标准差时，
# 相对统计误差通常输出为99%。
#
# 对于零通量，统计误差通常也为零。

INSUFFICIENT_ERROR_PERCENT = 99.0

# True：
#     整个曲线中隐藏误差达到或超过99%的非零数据点。
#
# False：
#     曲线仍保留这些数据点。
MASK_INSUFFICIENT_STATISTICS = False

# True：
#     error bar不显示误差达到或超过99%的点。
#
# False：
#     99%误差点也显示error bar。
SKIP_INSUFFICIENT_ERROR_BARS = True


# ------------------------------------------------------------
# I. X轴设置
# ------------------------------------------------------------
# USRBDX tab.lis中的能量通常以GeV给出。
X_AXIS_LABEL = r"particle energy, $E$ (GeV)"

USE_LOG_X = True

# 用户强制设置绘图能量范围。
# None表示使用完整能量范围。
PLOT_ENERGY_MIN = None
PLOT_ENERGY_MAX = None


# ------------------------------------------------------------
# J. Y轴标题
# ------------------------------------------------------------
# 微分通量并对立体角积分时，可以使用：
Y_AXIS_LABEL = (
    r"Differential neutron fluence"
    r"(cm$^{-2}$ GeV$^{-1}$ primary$^{-1}$)"
)

# 若读取的是未按面积归一化的particle current，
# 可以自行修改为：
#
# Y_AXIS_LABEL = (
#     r"Differential particle current "
#     r"(GeV$^{-1}$ primary$^{-1}$)"
# )
#
# 本程序不会额外除以面积或乘以面积，
# 而是直接使用tab.lis第三列的结果。


# ------------------------------------------------------------
# K. Y轴设置
# ------------------------------------------------------------
USE_LOG_Y = True

# 自动Y轴上限留白。
Y_HEADROOM_FACTOR = 1.5

# 对数Y轴下，自动下限相对于最小正值的倍率。
LOG_Y_BOTTOM_FACTOR = 0.7

# 用户强制设置Y轴范围。
# None表示自动。
USER_Y_MIN = None
USER_Y_MAX = None


# ------------------------------------------------------------
# L. 积分/求和
# ------------------------------------------------------------
CALCULATE_TOTAL = True

# "integral"：
#
#     sum(value_i * Delta_E_i)
#
#     适用于第三列为dPhi/dE、dN/dE等微分量。
#
# "sum"：
#
#     sum(value_i)
#
#     仅适用于第三列本身已经是每个能量bin的积分量。
TOTAL_MODE = "integral"

# True：
#     对最终绘图的缩放后结果计算。
#
# False：
#     对FLUKA原始结果计算。
TOTAL_USE_SCALED_PROFILE = True

TOTAL_NAME = "Energy-integrated USRBDX result"

# 对微分面通量积分后的示例单位：
TOTAL_UNIT = "cm^-2 primary^-1"

# 如果第三列是未按面积归一化的particle current，
# 可改为：
#
# TOTAL_UNIT = "primary^-1"

# 积分能量范围，单位通常为GeV。
#
# None表示完整能量范围。
INTEGRAL_ENERGY_MIN = None
INTEGRAL_ENERGY_MAX = None


# ------------------------------------------------------------
# M. CSV输出
# ------------------------------------------------------------
# True：
#     同时输出原始结果与缩放后结果。
#
# False：
#     只输出最终绘图结果。
EXPORT_RAW_AND_SCALED = True

# 是否将USRBDX第四列的相对统计误差输出到CSV。
EXPORT_RELATIVE_ERROR = True

# 是否将绘图使用的绝对误差输出到CSV。
#
# absolute_error =
# abs(plotted_value) * relative_error_percent / 100
EXPORT_ABSOLUTE_ERROR = True


# ------------------------------------------------------------
# N. 图标题
# ------------------------------------------------------------
# 不需要图标题时设为None。
PLOT_TITLE = None

# 示例：
#
# PLOT_TITLE = "USRBDX particle spectrum"


# ------------------------------------------------------------
# O. 输出文件名称
# ------------------------------------------------------------
OUTPUT_STEM = BASE_DIR / "usrbdx_energy_spectrum_neutron"


# ============================================================
# 3. USRBDX tab.lis文件解析
# ============================================================

FLOAT_RE = re.compile(
    r"[-+]?(?:\d+\.\d*|\.\d+|\d+)"
    r"(?:[EeDd][-+]?\d+)?"
)

DETECTOR_RE = re.compile(
    r"^\s*#\s*Detector\s+n:\s*"
    r"(\d+)\s+(.*?)\s*$",
    re.IGNORECASE,
)

ENERGY_INTERVAL_RE = re.compile(
    r"^\s*#\s*N\.\s*of\s*energy\s+intervals\s+"
    r"(\d+)\s*$",
    re.IGNORECASE,
)


def to_float(value: str) -> float:
    """
    将Fortran的E/D科学计数格式转换为Python浮点数。
    """
    return float(
        value.replace("D", "E").replace("d", "e")
    )


def validate_configuration():
    """
    检查用户配置是否有效。
    """
    if REFERENCE_CURRENT_MA <= 0:
        raise ValueError(
            "REFERENCE_CURRENT_MA必须大于0。"
        )

    if MARK_EVERY < 1:
        raise ValueError(
            "MARK_EVERY必须大于或等于1。"
        )

    if ERROR_BAR_EVERY < 1:
        raise ValueError(
            "ERROR_BAR_EVERY必须大于或等于1。"
        )

    if Y_HEADROOM_FACTOR <= 1.0:
        raise ValueError(
            "Y_HEADROOM_FACTOR建议设置为大于1。"
        )

    if LOG_Y_BOTTOM_FACTOR <= 0:
        raise ValueError(
            "LOG_Y_BOTTOM_FACTOR必须大于0。"
        )

    if TOTAL_MODE not in {
        "integral",
        "sum",
    }:
        raise ValueError(
            "TOTAL_MODE必须为"
            "'integral'或'sum'。"
        )

    if (
        PLOT_ENERGY_MIN is not None
        and PLOT_ENERGY_MAX is not None
        and PLOT_ENERGY_MAX
        <= PLOT_ENERGY_MIN
    ):
        raise ValueError(
            "PLOT_ENERGY_MAX必须大于"
            "PLOT_ENERGY_MIN。"
        )

    if (
        INTEGRAL_ENERGY_MIN is not None
        and INTEGRAL_ENERGY_MAX is not None
        and INTEGRAL_ENERGY_MAX
        <= INTEGRAL_ENERGY_MIN
    ):
        raise ValueError(
            "INTEGRAL_ENERGY_MAX必须大于"
            "INTEGRAL_ENERGY_MIN。"
        )


def validate_case_filenames():
    """
    检查CASES中的输入文件名是否已经填写。
    """
    empty_cases = []

    for index, (
        filename,
        _,
        _,
        _,
    ) in enumerate(
        CASES,
        start=1,
    ):
        if not str(filename).strip():
            empty_cases.append(index)

    if empty_cases:
        case_text = ", ".join(
            str(index)
            for index in empty_cases
        )

        raise ValueError(
            "以下CASES尚未填写USRBDX文件名："
            f"{case_text}\n\n"
            "请在CASES每个元组的第一项中填写"
            "*_tab.lis文件名。"
        )


def parse_usrbdx_detectors(path: Path):
    """
    读取USRBDX tab.lis中的所有四列表格。

    每个Detector的积分立体角能谱通常为：

        # Detector n: 1 detectorName
        # N. of energy intervals 50

        Emin  Emax  differential_value  error_percent

    本函数只读取每个Detector的第一张四列表格。

    在读取到指定数量的能量区间之后立即停止，
    因此不会把后面的double differential distributions
    当作积分立体角能谱读取。

    返回：
        detectors

    每个Detector保存为：

        {
            "number": Detector编号,
            "description": Detector说明,
            "n_intervals": 能量bin数,
            "energy_min": ndarray,
            "energy_max": ndarray,
            "value": ndarray,
            "error_percent": ndarray,
        }
    """
    if not path.is_file():
        raise FileNotFoundError(
            f"找不到输入文件：\n{path}"
        )

    text = path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    lines = text.splitlines()

    detectors = []
    line_index = 0

    while line_index < len(lines):
        detector_match = DETECTOR_RE.match(
            lines[line_index]
        )

        if detector_match is None:
            line_index += 1
            continue

        detector_number = int(
            detector_match.group(1)
        )

        detector_description = (
            detector_match.group(2).strip()
        )

        interval_count = None
        interval_line_index = None

        search_index = line_index + 1

        while search_index < len(lines):
            current_line = lines[search_index]

            if DETECTOR_RE.match(current_line):
                break

            if (
                "double differential distributions"
                in current_line.lower()
            ):
                break

            interval_match = (
                ENERGY_INTERVAL_RE.match(
                    current_line
                )
            )

            if interval_match is not None:
                interval_count = int(
                    interval_match.group(1)
                )

                interval_line_index = (
                    search_index
                )

                break

            search_index += 1

        if (
            interval_count is None
            or interval_line_index is None
        ):
            raise ValueError(
                f"{path.name}中Detector "
                f"{detector_number}后未找到"
                "能量区间数量。"
            )

        rows = []
        data_index = interval_line_index + 1

        while (
            data_index < len(lines)
            and len(rows) < interval_count
        ):
            current_line = lines[
                data_index
            ].strip()

            if not current_line:
                data_index += 1
                continue

            if current_line.startswith("#"):
                if (
                    "double differential distributions"
                    in current_line.lower()
                ):
                    break

                data_index += 1
                continue

            tokens = FLOAT_RE.findall(
                current_line
            )

            if len(tokens) >= 4:
                rows.append([
                    to_float(tokens[0]),
                    to_float(tokens[1]),
                    to_float(tokens[2]),
                    to_float(tokens[3]),
                ])

            data_index += 1

        if len(rows) != interval_count:
            raise ValueError(
                f"{path.name}中Detector "
                f"{detector_number}的数据数量不正确："
                f"实际读取{len(rows)}行，"
                f"预期{interval_count}行。"
            )

        table = np.asarray(
            rows,
            dtype=float,
        )

        energy_min = table[:, 0]
        energy_max = table[:, 1]
        values = table[:, 2]
        error_percent = table[:, 3]

        if np.any(
            energy_max <= energy_min
        ):
            raise ValueError(
                f"{path.name}中Detector "
                f"{detector_number}存在"
                "Emax <= Emin的能量区间。"
            )

        if np.any(
            np.diff(energy_min) < 0
        ):
            raise ValueError(
                f"{path.name}中Detector "
                f"{detector_number}的能量bin"
                "没有按照从低到高排列。"
            )

        detectors.append({
            "number": detector_number,
            "description": (
                detector_description
            ),
            "n_intervals": (
                interval_count
            ),
            "energy_min": energy_min,
            "energy_max": energy_max,
            "value": values,
            "error_percent": (
                error_percent
            ),
        })

        line_index = data_index

    if not detectors:
        raise ValueError(
            "无法在以下文件中找到USRBDX "
            "Detector四列表格：\n"
            f"{path}"
        )

    return detectors


def select_detector(
    detectors,
    path: Path,
):
    """
    根据DETECTOR_NUMBER和DETECTOR_NAME_CONTAINS
    选择需要处理的Detector。
    """
    matched = []

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
                in detector[
                    "description"
                ].lower()
            )

        if number_matches and name_matches:
            matched.append(detector)

    available_text = "\n".join(
        "  Detector "
        f"{item['number']}: "
        f"{item['description']}"
        for item in detectors
    )

    if len(matched) == 1:
        return matched[0]

    if not matched:
        raise ValueError(
            "在以下文件中没有找到符合条件的"
            f"Detector：\n{path}\n\n"
            "当前文件中的Detector为：\n"
            f"{available_text}"
        )

    raise ValueError(
        "Detector筛选条件匹配到多个结果：\n"
        f"{path}\n\n"
        "请同时设置DETECTOR_NUMBER和"
        "DETECTOR_NAME_CONTAINS，"
        "以唯一确定Detector。\n\n"
        f"匹配结果：\n{available_text}"
    )


def read_usrbdx_tab_lis(path: Path):
    """
    读取用户指定的USRBDX Detector。

    返回：
        detector_number
        detector_description
        energy_min
        energy_max
        energy_center
        raw_spectrum
        error_percent
    """
    detectors = parse_usrbdx_detectors(
        path
    )

    detector = select_detector(
        detectors,
        path,
    )

    energy_min = detector[
        "energy_min"
    ].copy()

    energy_max = detector[
        "energy_max"
    ].copy()

    raw_spectrum = detector[
        "value"
    ].copy()

    error_percent = detector[
        "error_percent"
    ].copy()

    # 对数能量bin使用几何中心。
    # 若存在非正能量，则使用算术中心。
    if np.all(
        (energy_min > 0)
        & (energy_max > 0)
    ):
        energy_center = np.sqrt(
            energy_min * energy_max
        )
    else:
        energy_center = (
            energy_min + energy_max
        ) / 2.0

    return (
        detector["number"],
        detector["description"],
        energy_min,
        energy_max,
        energy_center,
        raw_spectrum,
        error_percent,
    )


# ============================================================
# 4. 能量网格与数据显示
# ============================================================

def check_contiguous_energy_bins(
    energy_min,
    energy_max,
    path: Path,
):
    """
    检查相邻能量bin是否连续。

    阶梯图要求：

        Emax[i] == Emin[i + 1]
    """
    if energy_min.size <= 1:
        return

    if not np.allclose(
        energy_max[:-1],
        energy_min[1:],
        rtol=1.0e-7,
        atol=1.0e-14,
    ):
        raise ValueError(
            f"{path.name}的能量bin不连续，"
            "无法使用标准阶梯图。\n"
            "请检查tab.lis文件，或者设置：\n"
            "USE_CONTINUOUS_LINE = True"
        )


def get_energy_edges(
    energy_min,
    energy_max,
):
    """
    将连续的Emin和Emax转换为阶梯图边界数组。
    """
    return np.concatenate([
        energy_min[:1],
        energy_max,
    ])


def get_plot_range_mask(
    energy_min,
    energy_max,
):
    """
    返回与用户绘图能量范围存在重叠的bin。
    """
    mask = np.ones(
        energy_min.shape,
        dtype=bool,
    )

    if PLOT_ENERGY_MIN is not None:
        mask &= (
            energy_max
            > PLOT_ENERGY_MIN
        )

    if PLOT_ENERGY_MAX is not None:
        mask &= (
            energy_min
            < PLOT_ENERGY_MAX
        )

    if USE_LOG_X:
        mask &= (
            energy_max > 0
        )

    return mask


def prepare_profile_for_plot(
    energy_min,
    energy_max,
    profile,
    error_percent,
):
    """
    生成绘图用数据副本。

    处理内容：
        非有限值设为NaN；
        绘图能量范围之外设为NaN；
        对数Y轴下非正值设为NaN；
        可选隐藏统计误差达到99%的非零bin。
    """
    plotted = (
        profile.copy().astype(float)
    )

    plotted[
        ~np.isfinite(plotted)
    ] = np.nan

    range_mask = get_plot_range_mask(
        energy_min,
        energy_max,
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
    profile_for_plot,
    error_percent,
):
    """
    返回实际绘制error bar的数据点掩码。
    """
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

    mask &= interval_mask

    return mask


def get_visible_y_limits(
    energy_min,
    energy_max,
    profile,
    error_percent,
):
    """
    获取一条曲线实际显示的Y轴最小正值和最大值。

    当SHOW_ERROR_BARS=True时，
    自动范围同时考虑error bar的上下边界。
    """
    plotted = prepare_profile_for_plot(
        energy_min,
        energy_max,
        profile,
        error_percent,
    )

    finite_mask = np.isfinite(
        plotted
    )

    if not np.any(finite_mask):
        return np.inf, 0.0

    finite_values = plotted[
        finite_mask
    ]

    if USE_LOG_Y:
        positive_values = finite_values[
            finite_values > 0
        ]

        if positive_values.size == 0:
            minimum_positive = np.inf
        else:
            minimum_positive = float(
                np.min(positive_values)
            )
    else:
        minimum_positive = float(
            np.min(finite_values)
        )

    maximum_value = float(
        np.max(finite_values)
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

            finite_upper = upper_values[
                np.isfinite(upper_values)
            ]

            if finite_upper.size > 0:
                maximum_value = max(
                    maximum_value,
                    float(
                        np.max(finite_upper)
                    ),
                )

            lower_values = (
                y_values
                - absolute_error
            )

            if USE_LOG_Y:
                valid_lower = lower_values[
                    np.isfinite(lower_values)
                    & (lower_values > 0)
                ]
            else:
                valid_lower = lower_values[
                    np.isfinite(lower_values)
                ]

            if valid_lower.size > 0:
                minimum_positive = min(
                    minimum_positive,
                    float(
                        np.min(valid_lower)
                    ),
                )

    return (
        minimum_positive,
        maximum_value,
    )


# ============================================================
# 5. 绘制单条USRBDX能谱
# ============================================================

def plot_one_spectrum(
    ax,
    energy_min,
    energy_max,
    energy_center,
    profile,
    error_percent,
    label,
    color,
):
    """
    绘制一条USRBDX能谱。

    支持：
        连续bin中心折线；
        能量bin阶梯图；
        bin中心marker；
        阶梯图辅助连接线；
        相对统计误差error bar；
        统计不足bin隐藏；
        线性或对数坐标。
    """
    profile_for_plot = (
        prepare_profile_for_plot(
            energy_min,
            energy_max,
            profile,
            error_percent,
        )
    )

    if not np.any(
        np.isfinite(profile_for_plot)
    ):
        return

    # --------------------------------------------------------
    # 主曲线
    # --------------------------------------------------------
    if USE_CONTINUOUS_LINE:
        ax.plot(
            energy_center,
            profile_for_plot,
            color=color,
            linewidth=2.0,
            label=label,
            zorder=3,
        )

    else:
        energy_edges = get_energy_edges(
            energy_min,
            energy_max,
        )

        ax.stairs(
            values=profile_for_plot,
            edges=energy_edges,
            color=color,
            linewidth=2.0,
            label=label,
            zorder=3,
        )

    # --------------------------------------------------------
    # bin中心marker
    # --------------------------------------------------------
    if SHOW_BIN_MARKERS:
        marker_mask = np.zeros(
            energy_center.shape,
            dtype=bool,
        )

        marker_mask[
            np.arange(
                energy_center.size
            )[::MARK_EVERY]
        ] = True

        marker_mask &= np.isfinite(
            profile_for_plot
        )

        ax.plot(
            energy_center[marker_mask],
            profile_for_plot[marker_mask],
            linestyle="none",
            marker=BIN_MARKER,
            markersize=BIN_MARKER_SIZE,
            color=color,
            label="_nolegend_",
            zorder=4,
        )

    # --------------------------------------------------------
    # 阶梯模式下的bin中心辅助连接线
    # --------------------------------------------------------
    if (
        not USE_CONTINUOUS_LINE
        and ADD_GUIDE_LINE
    ):
        ax.plot(
            energy_center,
            profile_for_plot,
            color=color,
            linewidth=GUIDE_LINE_WIDTH,
            alpha=GUIDE_LINE_ALPHA,
            solid_capstyle="round",
            label="_nolegend_",
            zorder=4,
        )

    # --------------------------------------------------------
    # Error bar
    # --------------------------------------------------------
    if SHOW_ERROR_BARS:
        error_mask = get_error_bar_mask(
            profile_for_plot,
            error_percent,
        )

        if np.any(error_mask):
            x_values = energy_center[
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

            # 对数坐标下，error bar下端不能小于或等于0。
            #
            # 这里只限制绘图时的下误差长度，
            # 不修改原始误差值或CSV输出。
            if USE_LOG_Y:
                lower_error = np.minimum(
                    absolute_error,
                    y_values
                    * (1.0 - 1.0e-12),
                )

                upper_error = (
                    absolute_error
                )

                y_error = np.vstack([
                    lower_error,
                    upper_error,
                ])

            else:
                y_error = absolute_error

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
                alpha=ERROR_BAR_ALPHA,
                label="_nolegend_",
                zorder=2,
            )


# ============================================================
# 6. 能量积分或求和
# ============================================================

def calculate_spectrum_total(
    energy_min,
    energy_max,
    profile,
):
    """
    对指定能量范围内的能谱积分或求和。

    integral模式：
        假设每个能量bin内的微分值为常数，
        根据积分范围与bin的实际重叠宽度计算：

            sum(value_i * overlap_width_i)

        当积分边界位于某个bin内部时，
        只计算该bin与积分范围重叠的部分。

    sum模式：
        对与指定能量范围有重叠的bin直接求和。
    """
    lower_limit = (
        float(np.min(energy_min))
        if INTEGRAL_ENERGY_MIN is None
        else float(INTEGRAL_ENERGY_MIN)
    )

    upper_limit = (
        float(np.max(energy_max))
        if INTEGRAL_ENERGY_MAX is None
        else float(INTEGRAL_ENERGY_MAX)
    )

    if upper_limit <= lower_limit:
        raise ValueError(
            "积分能量上限必须大于下限。"
        )

    overlap_left = np.maximum(
        energy_min,
        lower_limit,
    )

    overlap_right = np.minimum(
        energy_max,
        upper_limit,
    )

    overlap_width = np.maximum(
        overlap_right - overlap_left,
        0.0,
    )

    valid_mask = (
        np.isfinite(profile)
        & (overlap_width > 0)
    )

    if not np.any(valid_mask):
        raise ValueError(
            "指定能量范围内没有有效USRBDX数据。"
        )

    if TOTAL_MODE == "integral":
        total_value = float(
            np.sum(
                profile[valid_mask]
                * overlap_width[valid_mask]
            )
        )

    else:
        total_value = float(
            np.sum(
                profile[valid_mask]
            )
        )

    return (
        total_value,
        lower_limit,
        upper_limit,
    )


# ============================================================
# 7. CSV列名处理
# ============================================================

def make_clean_csv_name(
    path: Path,
    case_index: int,
):
    """
    生成适合CSV列名的工况名称。
    """
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
# 8. 主程序
# ============================================================

def main():
    validate_configuration()
    validate_case_filenames()
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": [
            "Times New Roman",
            "STIXGeneral",
            "SimSun",
            "DejaVu Serif",
            "Noto Serif CJK SC",              
        ],
        "axes.unicode_minus": False,
        "mathtext.fontset": "cm",
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

    color_cycle = (
        plt.rcParams[
            "axes.prop_cycle"
        ].by_key().get(
            "color",
            ["C0"],
        )
    )

    csv_columns = []
    csv_names = []
    total_results = []

    reference_energy_min = None
    reference_energy_max = None
    reference_energy_center = None

    selected_detector_number = None
    selected_detector_description = None

    global_ymax = 0.0
    global_ymin = np.inf

    print()
    print(
        "USRBDX differential-spectrum processing"
    )
    print("-" * 76)

    if APPLY_CURRENT_SCALING:
        print(
            "Current scaling is enabled:"
        )
        print(
            "  current scale = "
            "current_mA / "
            f"{REFERENCE_CURRENT_MA:g} mA"
        )
        print("-" * 76)

    if SHOW_ERROR_BARS:
        print(
            "Statistical error bars are enabled."
        )
        print(
            "  Error values are read from "
            "the fourth USRBDX column."
        )
        print("-" * 76)

    if TOTAL_MODE == "sum":
        print(
            "Warning: TOTAL_MODE='sum'直接对"
            "第三列数据求和。"
        )
        print(
            "只有当第三列本身已经是每个能量bin的"
            "积分量时，该操作才具有正确物理意义。"
        )
        print("-" * 76)

    for case_index, (
        filename,
        label,
        current_ma,
        user_y_scale,
    ) in enumerate(
        CASES,
        start=1,
    ):
        path = BASE_DIR / filename

        (
            detector_number,
            detector_description,
            energy_min,
            energy_max,
            energy_center,
            raw_spectrum,
            error_percent,
        ) = read_usrbdx_tab_lis(
            path
        )

        if not USE_CONTINUOUS_LINE:
            check_contiguous_energy_bins(
                energy_min,
                energy_max,
                path,
            )

        # ----------------------------------------------------
        # 电流缩放
        # ----------------------------------------------------
        if APPLY_CURRENT_SCALING:
            current_scale = (
                current_ma
                / REFERENCE_CURRENT_MA
            )
        else:
            current_scale = 1.0

        effective_scale = (
            current_scale
            * user_y_scale
        )

        plotted_spectrum = (
            raw_spectrum
            * effective_scale
        )

        absolute_error_plotted = (
            np.abs(plotted_spectrum)
            * error_percent
            / 100.0
        )

        # ----------------------------------------------------
        # 检查能量网格
        # ----------------------------------------------------
        if reference_energy_min is None:
            reference_energy_min = (
                energy_min.copy()
            )

            reference_energy_max = (
                energy_max.copy()
            )

            reference_energy_center = (
                energy_center.copy()
            )

            selected_detector_number = (
                detector_number
            )

            selected_detector_description = (
                detector_description
            )

            csv_columns.extend([
                reference_energy_min,
                reference_energy_max,
                reference_energy_center,
            ])

            csv_names.extend([
                "energy_min_GeV",
                "energy_max_GeV",
                "energy_center_GeV",
            ])

        else:
            same_shape = (
                energy_min.shape
                == reference_energy_min.shape
            )

            same_min = (
                same_shape
                and np.allclose(
                    energy_min,
                    reference_energy_min,
                    rtol=1.0e-7,
                    atol=1.0e-14,
                )
            )

            same_max = (
                same_shape
                and np.allclose(
                    energy_max,
                    reference_energy_max,
                    rtol=1.0e-7,
                    atol=1.0e-14,
                )
            )

            if not (
                same_min and same_max
            ):
                raise ValueError(
                    f"{path.name}的能量网格"
                    "与第一个文件不同。\n"
                    "当前统一CSV输出要求所有工况"
                    "使用相同的能量bin。"
                )

        # ----------------------------------------------------
        # 图例名称
        # ----------------------------------------------------
        display_label = label

        if (
            SHOW_USER_SCALE_IN_LEGEND
            and not np.isclose(
                user_y_scale,
                1.0,
            )
        ):
            display_label += (
                rf" ($\times"
                rf"{user_y_scale:g}$)"
            )

        # ----------------------------------------------------
        # 更新自动Y轴范围
        # ----------------------------------------------------
        (
            local_ymin,
            local_ymax,
        ) = get_visible_y_limits(
            energy_min,
            energy_max,
            plotted_spectrum,
            error_percent,
        )

        if np.isfinite(local_ymin):
            global_ymin = min(
                global_ymin,
                local_ymin,
            )

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

        plot_one_spectrum(
            ax=ax,
            energy_min=energy_min,
            energy_max=energy_max,
            energy_center=energy_center,
            profile=plotted_spectrum,
            error_percent=error_percent,
            label=display_label,
            color=curve_color,
        )

        # ----------------------------------------------------
        # CSV输出
        # ----------------------------------------------------
        clean_name = make_clean_csv_name(
            path,
            case_index,
        )

        if EXPORT_RAW_AND_SCALED:
            csv_columns.append(
                raw_spectrum
            )

            csv_names.append(
                clean_name + "_raw"
            )

            csv_columns.append(
                plotted_spectrum
            )

            csv_names.append(
                clean_name + "_scaled"
            )

        else:
            csv_columns.append(
                plotted_spectrum
            )

            csv_names.append(
                clean_name + "_plotted"
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
                absolute_error_plotted
            )

            csv_names.append(
                clean_name
                + "_plotted_absolute_error"
            )

        # ----------------------------------------------------
        # 积分或求和
        # ----------------------------------------------------
        total_value = None
        total_energy_min = None
        total_energy_max = None

        if CALCULATE_TOTAL:
            if TOTAL_USE_SCALED_PROFILE:
                profile_for_total = (
                    plotted_spectrum
                )

                total_source = "scaled"
            else:
                profile_for_total = (
                    raw_spectrum
                )

                total_source = "raw"

            (
                total_value,
                total_energy_min,
                total_energy_max,
            ) = calculate_spectrum_total(
                energy_min,
                energy_max,
                profile_for_total,
            )

            total_results.append({
                "file": path.name,
                "label": label,
                "detector_number": (
                    detector_number
                ),
                "detector_description": (
                    detector_description
                ),
                "current_mA": current_ma,
                "current_scale": (
                    current_scale
                ),
                "user_y_scale": (
                    user_y_scale
                ),
                "effective_scale": (
                    effective_scale
                ),
                "total_mode": TOTAL_MODE,
                "total_source": total_source,
                "energy_min_GeV": (
                    total_energy_min
                ),
                "energy_max_GeV": (
                    total_energy_max
                ),
                "total": total_value,
                "unit": TOTAL_UNIT,
            })

        # ----------------------------------------------------
        # 终端输出
        # ----------------------------------------------------
        finite_raw = raw_spectrum[
            np.isfinite(raw_spectrum)
        ]

        if finite_raw.size > 0:
            peak_raw_value = float(
                np.max(finite_raw)
            )
        else:
            peak_raw_value = np.nan

        insufficient_bin_count = int(
            np.count_nonzero(
                np.isfinite(error_percent)
                & (
                    error_percent
                    >= INSUFFICIENT_ERROR_PERCENT
                )
                & np.isfinite(raw_spectrum)
                & (raw_spectrum != 0)
            )
        )

        print(path.name)

        print(
            "  detector "
            f"= {detector_number}"
        )

        print(
            "  detector description "
            f"= {detector_description}"
        )

        print(
            "  number of energy bins "
            f"= {energy_min.size}"
        )

        print(
            "  current scaling factor "
            f"= {current_scale:g}"
        )

        print(
            "  user Y scaling factor "
            f"= {user_y_scale:g}"
        )

        print(
            "  effective scaling factor "
            f"= {effective_scale:g}"
        )

        print(
            "  peak raw differential value "
            f"= {peak_raw_value:.8e}"
        )

        print(
            "  bins with insufficient "
            "statistics "
            f"= {insufficient_bin_count}"
        )

        if CALCULATE_TOTAL:
            print(
                f"  {TOTAL_NAME} "
                f"= {total_value:.8e} "
                f"{TOTAL_UNIT}"
            )

            print(
                "  total energy range "
                f"= [{total_energy_min:.8e}, "
                f"{total_energy_max:.8e}] GeV"
            )

        print("-" * 76)

    # ========================================================
    # 坐标轴设置
    # ========================================================
    if reference_energy_min is None:
        raise ValueError(
            "没有读取到任何USRBDX数据。"
        )

    # --------------------------------------------------------
    # X轴
    # --------------------------------------------------------
    if USE_LOG_X:
        if np.any(
            reference_energy_min <= 0
        ):
            raise ValueError(
                "对数X轴要求所有显示能量大于0。"
            )

        ax.set_xscale("log")

    else:
        ax.set_xscale("linear")

    if PLOT_ENERGY_MIN is None:
        x_min = float(
            np.min(reference_energy_min)
        )
    else:
        x_min = float(
            PLOT_ENERGY_MIN
        )

    if PLOT_ENERGY_MAX is None:
        x_max = float(
            np.max(reference_energy_max)
        )
    else:
        x_max = float(
            PLOT_ENERGY_MAX
        )

    if x_max <= x_min:
        raise ValueError(
            "PLOT_ENERGY_MAX必须大于"
            "PLOT_ENERGY_MIN。"
        )

    if USE_LOG_X and x_min <= 0:
        raise ValueError(
            "对数X轴的PLOT_ENERGY_MIN"
            "必须大于0。"
        )

    ax.set_xlim(
        x_min,
        x_max,
    )

    # --------------------------------------------------------
    # Y轴
    # --------------------------------------------------------
    if global_ymax <= 0.0:
        raise ValueError(
            "实际绘图范围内没有有效的正USRBDX数据。"
        )

    if USE_LOG_Y:
        if not np.isfinite(
            global_ymin
        ):
            raise ValueError(
                "对数Y轴模式下没有找到"
                "大于0的USRBDX数据。"
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
                    "对数Y轴的USER_Y_MIN"
                    "必须大于0。"
                )

            y_min = float(
                USER_Y_MIN
            )

        if USER_Y_MAX is None:
            y_max = (
                global_ymax
                * Y_HEADROOM_FACTOR
            )
        else:
            y_max = float(
                USER_Y_MAX
            )

        ax.set_ylim(
            y_min,
            y_max,
        )

    else:
        ax.set_yscale("linear")

        y_min = (
            0.0
            if USER_Y_MIN is None
            else float(USER_Y_MIN)
        )

        y_max = (
            global_ymax
            * Y_HEADROOM_FACTOR
            if USER_Y_MAX is None
            else float(USER_Y_MAX)
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
    # 网格
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

    # ========================================================
    # 图例
    # ========================================================
    handles, legend_labels = (
        ax.get_legend_handles_labels()
    )

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
        top=0.94,
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
    # 保存曲线CSV
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
    # 保存积分结果CSV
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
                    "detector_description",
                    "current_mA",
                    "current_scale",
                    "user_y_scale",
                    "effective_scale",
                    "total_mode",
                    "total_source",
                    "energy_min_GeV",
                    "energy_max_GeV",
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
    print("Selected USRBDX detector:")

    print(
        "  number      = "
        f"{selected_detector_number}"
    )

    print(
        "  description = "
        f"{selected_detector_description}"
    )

    print()
    print("Output files:")
    print(png_path)
    print(pdf_path)
    print(svg_path)
    print(curve_csv_path)

    if CALCULATE_TOTAL:
        print(total_csv_path)

    plt.show()


if __name__ == "__main__":
    main()
