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
FLUKA .inp 几何解析器（最小集）。

支持范围：
  - Body 类型: SPH, ZCC, PLA, XYP（预留 XCC/YCC/XZP/YZP/RPP/BOX/RCC 接口）
  - Region 表达式: +body / -body / |（或）/ & 或空格（与）/ ()
  - 材料映射: ASSIGNMA
  - 切面材料分布: cut_plane() 返回 MaterialGrid

符号约定（已与用户确认）:
  +body  = 内部 / 法向量反方向半空间
  -body  = 外部 / 法向量所指方向半空间
  平面法向量指向 - 区域；XYP/YZP/XZP 法向分别为 +z/+x/+y。

模块可独立 import，不依赖 usrbin 读图程序。
未来可作为 fluka-ai-agent 后处理 tool 使用。
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

import numpy as np


# ============================================================
# 1. Body 基类与最小集实现
# ============================================================

class Body:
    """几何体基类。子类需实现 inside() 和 analytic_boundary_on_plane()。"""

    name: str

    def inside(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        """返回与 x/y/z 同形状的 bool 数组，True = 在 body 内部（+body 半空间）。"""
        raise NotImplementedError

    def signed_distance(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        """有符号距离：<=0 表示 +body 内部，>0 表示 -body 外部。默认抛错。"""
        raise NotImplementedError

    def analytic_boundary_on_plane(
        self, axis: str, value: float
    ) -> list[dict]:
        """在切面 axis=value 上返回解析边界曲线列表。

        每条曲线为 dict: {"free1": np.ndarray, "free2": np.ndarray}
        其中 free1/free2 对应切面剩余两个轴的坐标。
        无解析式时返回 []（由调用方回退到 contour）。
        """
        return []


class Sph(Body):
    """SPH name cx cy cz r —— 球"""

    def __init__(self, name: str, cx: float, cy: float, cz: float, r: float):
        self.name = name
        self.cx, self.cy, self.cz, self.r = float(cx), float(cy), float(cz), float(r)

    def signed_distance(self, x, y, z):
        return np.sqrt((x - self.cx) ** 2 + (y - self.cy) ** 2 + (z - self.cz) ** 2) - self.r

    def inside(self, x, y, z):
        return self.signed_distance(x, y, z) <= 0

    def analytic_boundary_on_plane(self, axis, value):
        # 球被任一切面截，交线为圆
        if axis == "x":
            d2 = self.r ** 2 - (value - self.cx) ** 2
        elif axis == "y":
            d2 = self.r ** 2 - (value - self.cy) ** 2
        elif axis == "z":
            d2 = self.r ** 2 - (value - self.cz) ** 2
        else:
            return []
        if d2 <= 0:
            return []
        r_cut = np.sqrt(d2)
        t = np.linspace(0, 2 * np.pi, 361)
        # 切面剩余两轴：x切->(y,z), y切->(x,z), z切->(x,y)
        if axis == "x":
            f1 = self.cy + r_cut * np.sin(t)   # y
            f2 = self.cz + r_cut * np.cos(t)    # z
        elif axis == "y":
            f1 = self.cx + r_cut * np.cos(t)    # x
            f2 = self.cz + r_cut * np.sin(t)    # z
        else:  # z
            f1 = self.cx + r_cut * np.cos(t)    # x
            f2 = self.cy + r_cut * np.sin(t)    # y
        return [{"free1": f1, "free2": f2}]


class Zcc(Body):
    """ZCC name cx cy r —— 轴平行 z 的无限圆柱"""

    def __init__(self, name: str, cx: float, cy: float, r: float):
        self.name = name
        self.cx, self.cy, self.r = float(cx), float(cy), float(r)

    def signed_distance(self, x, y, z):
        return np.sqrt((x - self.cx) ** 2 + (y - self.cy) ** 2) - self.r

    def inside(self, x, y, z):
        return self.signed_distance(x, y, z) <= 0

    def analytic_boundary_on_plane(self, axis, value):
        # z 切面：圆
        if axis == "z":
            t = np.linspace(0, 2 * np.pi, 361)
            return [{"free1": self.cx + self.r * np.cos(t),
                     "free2": self.cy + self.r * np.sin(t)}]
        # x 切面：两条水平直线 y = cy ± r（z 自由）
        elif axis == "x":
            d2 = self.r ** 2 - (value - self.cx) ** 2
            if d2 <= 0:
                return []
            yc = self.cy + np.sqrt(d2)
            yc2 = self.cy - np.sqrt(d2)
            return [
                {"free1": np.array([yc, yc]), "free2": np.array([-1e6, 1e6])},
                {"free1": np.array([yc2, yc2]), "free2": np.array([-1e6, 1e6])},
            ]
        # y 切面：两条垂直直线 x = cx ± r（z 自由）
        elif axis == "y":
            d2 = self.r ** 2 - (value - self.cy) ** 2
            if d2 <= 0:
                return []
            xc = self.cx + np.sqrt(d2)
            xc2 = self.cx - np.sqrt(d2)
            return [
                {"free1": np.array([xc, xc]), "free2": np.array([-1e6, 1e6])},
                {"free1": np.array([xc2, xc2]), "free2": np.array([-1e6, 1e6])},
            ]
        return []


class Plane(Body):
    """平面类基类。法向量指向 - 区域（外部）。"""

    def analytic_boundary_on_plane(self, axis, value):
        # 平面在切面上通常为一条直线；具体子类决定
        return []


class Xyp(Plane):
    """XYP name z0 —— 平面 z = z0，法向 +z"""

    def __init__(self, name: str, z0: float):
        self.name = name
        self.z0 = float(z0)

    def signed_distance(self, x, y, z):
        # 法向 +z 指向 - 区域，所以 + 区域是 z <= z0，signed_distance = z - z0
        return z - self.z0

    def inside(self, x, y, z):
        return z <= self.z0

    def analytic_boundary_on_plane(self, axis, value):
        # XYP 是 z=z0 平面（法向 +z）。
        # 在非 z 轴切面上，交线为 z=z0 的水平线（自由轴方向延伸）。
        # 在 z 轴切面上，仅当 value==z0 时重合（整面），否则不相交。
        if axis == "z":
            return []
        # axis='x' -> free1=y, free2=z；axis='y' -> free1=x, free2=z
        # 交线：free2 固定为 z0，free1 自由
        big = 1.0e6
        return [{"free1": np.array([-big, big]), "free2": np.array([self.z0, self.z0])}]


class Xzp(Plane):
    """XZP name y0 —— 平面 y = y0，法向 +y"""

    def __init__(self, name: str, y0: float):
        self.name = name
        self.y0 = float(y0)

    def signed_distance(self, x, y, z):
        return y - self.y0

    def inside(self, x, y, z):
        return y <= self.y0

    def analytic_boundary_on_plane(self, axis, value):
        # XZP 是 y=y0 平面（法向 +y）。
        # 在非 y 轴切面上，交线为 y=y0 的水平线。
        # 在 y 轴切面上，仅当 value==y0 时重合，否则不相交。
        if axis == "y":
            return []
        # axis='x' -> free1=y, free2=z；axis='z' -> free1=x, free2=y
        # 交线：free1 固定为 y0（axis='x'）或 free2 固定为 y0（axis='z'）
        big = 1.0e6
        if axis == "x":
            return [{"free1": np.array([self.y0, self.y0]),
                     "free2": np.array([-big, big])}]
        else:  # axis == "z"
            return [{"free1": np.array([-big, big]),
                     "free2": np.array([self.y0, self.y0])}]


class Yzp(Plane):
    """YZP name x0 —— 平面 x = x0，法向 +x"""

    def __init__(self, name: str, x0: float):
        self.name = name
        self.x0 = float(x0)

    def signed_distance(self, x, y, z):
        return x - self.x0

    def inside(self, x, y, z):
        return x <= self.x0

    def analytic_boundary_on_plane(self, axis, value):
        # YZP 是 x=x0 平面（法向 +x）。
        # 在非 x 轴切面上，交线为 x=x0 的垂直线。
        # 在 x 轴切面上，仅当 value==x0 时重合，否则不相交。
        if axis == "x":
            return []
        # axis='y' -> free1=x, free2=z；axis='z' -> free1=x, free2=y
        # 交线：free1 固定为 x0
        big = 1.0e6
        return [{"free1": np.array([self.x0, self.x0]),
                 "free2": np.array([-big, big])}]


class Pla(Plane):
    """PLA name Vx Vy Vz X1 Y1 Z1 —— 任意平面"""

    def __init__(self, name: str, vx, vy, vz, x1, y1, z1):
        self.name = name
        self.v = np.array([float(vx), float(vy), float(vz)], dtype=float)
        self.p0 = np.array([float(x1), float(y1), float(z1)], dtype=float)

    def signed_distance(self, x, y, z):
        # V·(p - p0)，<=0 为 + 区域
        return (self.v[0] * (x - self.p0[0])
                + self.v[1] * (y - self.p0[1])
                + self.v[2] * (z - self.p0[2]))

    def inside(self, x, y, z):
        return self.signed_distance(x, y, z) <= 0

    def analytic_boundary_on_plane(self, axis, value):
        # 平面 V·(p-p0)=0 与切面 axis=value 的交线为一条直线
        # 解析表达：在切面上，一个自由轴随另一个线性变化
        if axis == "x":
            # x 固定，求 y(z) 或 z(y)
            # Vx(x-x0) + Vy(y-y0) + Vz(z-z0) = 0
            if abs(self.v[1]) < 1e-12 and abs(self.v[2]) < 1e-12:
                return []
            # 取 z 作自由轴 f2，y 作 f1：y = y0 - (Vx(x-x0) + Vz(z-z0)) / Vy
            if abs(self.v[1]) >= 1e-12:
                z_lin = np.array([-1e6, 1e6])
                y_lin = (self.p0[1]
                         - (self.v[0] * (value - self.p0[0])
                            + self.v[2] * (z_lin - self.p0[2])) / self.v[1])
                return [{"free1": y_lin, "free2": z_lin}]
            else:
                # Vy=0，平面是 x=const 类型，在 x 切面上要么全有要么全无
                return []
        elif axis == "y":
            if abs(self.v[0]) < 1e-12 and abs(self.v[2]) < 1e-12:
                return []
            if abs(self.v[0]) >= 1e-12:
                z_lin = np.array([-1e6, 1e6])
                x_lin = (self.p0[0]
                         - (self.v[1] * (value - self.p0[1])
                            + self.v[2] * (z_lin - self.p0[2])) / self.v[0])
                return [{"free1": x_lin, "free2": z_lin}]
            else:
                return []
        elif axis == "z":
            if abs(self.v[0]) < 1e-12 and abs(self.v[1]) < 1e-12:
                return []
            if abs(self.v[0]) >= 1e-12:
                y_lin = np.array([-1e6, 1e6])
                x_lin = (self.p0[0]
                         - (self.v[1] * (y_lin - self.p0[1])
                            + self.v[2] * (value - self.p0[2])) / self.v[0])
                return [{"free1": x_lin, "free2": y_lin}]
            else:
                return []
        return []


# ============================================================
# 2. Body 注册器（扩展点：新增 body 类型在此登记）
# ============================================================

# 每个 parser 接收 tokens（list[str]）和 name，返回 Body 实例
BODY_PARSERS: dict[str, Callable[[str, list[str]], Body]] = {
    "SPH": lambda name, t: Sph(name, t[0], t[1], t[2], t[3]),
    "ZCC": lambda name, t: Zcc(name, t[0], t[1], t[2]),
    "XYP": lambda name, t: Xyp(name, t[0]),
    "XZP": lambda name, t: Xzp(name, t[0]),
    "YZP": lambda name, t: Yzp(name, t[0]),
    "PLA": lambda name, t: Pla(name, t[0], t[1], t[2], t[3], t[4], t[5]),
    # 扩展示例（未实现）:
    # "XCC": lambda name, t: Xcc(name, t[0], t[1], t[2]),
    # "RPP": lambda name, t: Rpp(name, *t[:6]),
}


# ============================================================
# 3. CSG 表达式解析（递归下降）
# ============================================================

# AST 节点
@dataclass
class BodyRef:
    name: str
    positive: bool  # True = +body, False = -body

    def eval(self, bodies: dict, x, y, z) -> np.ndarray:
        mask = bodies[self.name].inside(x, y, z)
        return mask if self.positive else ~mask


@dataclass
class AndNode:
    children: list

    def eval(self, bodies, x, y, z):
        result = self.children[0].eval(bodies, x, y, z)
        for c in self.children[1:]:
            result = result & c.eval(bodies, x, y, z)
        return result


@dataclass
class OrNode:
    children: list

    def eval(self, bodies, x, y, z):
        result = self.children[0].eval(bodies, x, y, z)
        for c in self.children[1:]:
            result = result | c.eval(bodies, x, y, z)
        return result


class CSGParser:
    """递归下降解析 FLUKA region 表达式。

    文法:
      or_expr   := and_expr ('|' and_expr)*
      and_expr  := term (& term | term)*    # & 或 空格分隔
      term      := '+' body | '-' body | '(' or_expr ')'
                 | '+' '(' or_expr ')'       # +(...) = 括号本身
                 | '-' '(' or_expr ')'       # -(...) = 取反
    """

    def __init__(self, tokens: list[str]):
        self.tokens = tokens
        self.pos = 0

    def peek(self) -> Optional[str]:
        return self.tokens[self.pos] if self.pos < len(self.tokens) else None

    def next(self) -> str:
        t = self.tokens[self.pos]
        self.pos += 1
        return t

    def parse(self):
        node = self.parse_or()
        return node

    def parse_or(self):
        children = [self.parse_and()]
        while self.peek() == "|":
            self.next()
            children.append(self.parse_and())
        return OrNode(children) if len(children) > 1 else children[0]

    def parse_and(self):
        children = [self.parse_term()]
        while True:
            t = self.peek()
            if t is None:
                break
            if t == "|":
                break
            if t == ")":
                break
            if t == "&":
                self.next()
                children.append(self.parse_term())
            elif t in ("+", "-"):
                # 直接跟一个 body 或括号，无显式 &，FLUKA 默认 and
                children.append(self.parse_term())
            else:
                # 裸 body 名（无前导符号，默认 +）
                children.append(self.parse_term())
        return AndNode(children) if len(children) > 1 else children[0]

    def parse_term(self):
        t = self.peek()
        if t is None:
            raise ValueError("表达式意外结束")
        if t == "+":
            self.next()
            nxt = self.peek()
            if nxt == "(":
                self.next()
                node = self.parse_or()
                assert self.next() == ")"
                return node  # +(...) = 括号本身
            else:
                name = self.next()
                return BodyRef(name, True)
        elif t == "-":
            self.next()
            nxt = self.peek()
            if nxt == "(":
                self.next()
                node = self.parse_or()
                assert self.next() == ")"
                # 取反
                return NegateNode(node)
            else:
                name = self.next()
                return BodyRef(name, False)
        elif t == "(":
            self.next()
            node = self.parse_or()
            assert self.next() == ")"
            return node
        else:
            # 裸 body 名，默认 +
            name = self.next()
            return BodyRef(name, True)


@dataclass
class NegateNode:
    child: object

    def eval(self, bodies, x, y, z):
        return ~self.child.eval(bodies, x, y, z)


def tokenize_expr(expr: str) -> list[str]:
    """将 region 表达式切分为 tokens。

    处理: + - | & ( ) 以及 body 名。
    注意 body 名可能含字母数字下划线。
    """
    # 在符号两侧加空格再分割
    s = expr
    for ch in ("+", "-", "|", "&", "(", ")"):
        s = s.replace(ch, f" {ch} ")
    raw = s.split()
    # 合并：但要注意 +- 紧跟 body 名时已被分开，OK
    tokens = []
    i = 0
    while i < len(raw):
        tokens.append(raw[i])
        i += 1
    return tokens


# ============================================================
# 4. Region 与 Geometry
# ============================================================

@dataclass
class Region:
    name: str
    expr_tree: object           # AST 根节点
    material: Optional[str] = None
    material_alias: Optional[str] = None  # ASSIGNMA 第二字段（材料别名）

    def inside(self, bodies: dict, x, y, z) -> np.ndarray:
        return self.expr_tree.eval(bodies, x, y, z)


@dataclass
class Geometry:
    bodies: dict = field(default_factory=dict)   # name -> Body
    regions: list = field(default_factory=list)   # list[Region]
    material_to_regions: dict = field(default_factory=dict)  # material -> [region_names]


# ============================================================
# 5. .inp 文件解析
# ============================================================

# Body 行：类型 + name + 参数
# 例: "SPH blkbody    0.0 0.0 0.0 10000.0"
BODY_TYPE_RE = re.compile(r"^\s*(SPH|ZCC|XCC|YCC|RCC|BOX|RPP|REC|TRC|PLA|XYP|XZP|YZP|CCC|ARB|WED|RAW)\b", re.IGNORECASE)


def _to_float(s: str) -> float:
    """转换 FLUKA 数值（支持 D/d 科学记数法）。"""
    s = s.strip()
    if not s:
        raise ValueError("空数值")
    s = s.replace("D", "E").replace("d", "e")
    return float(s)


def parse_inp(path) -> Geometry:
    """解析 .inp 文件，返回 Geometry。

    流程:
      1. 定位 GEOBEGIN ... GEOEND
      2. 解析 body 卡片（SPH/ZCC/...）
      3. 解析 region 卡片（CSG 表达式）
      4. 定位 ASSIGNMA，关联 region -> material
    """
    path = Path(path)
    text = path.read_text(encoding="utf-8", errors="replace")

    # 标准化行
    lines = text.splitlines()

    # ---- 1. 截取 GEOBEGIN..GEOEND ----
    geobegin_idx = None
    geoend_idx = None
    for i, line in enumerate(lines):
        up = line.upper()
        if "GEOBEGIN" in up and geobegin_idx is None:
            geobegin_idx = i
        if "GEOEND" in up and geobegin_idx is not None and geoend_idx is None:
            geoend_idx = i
            break
    if geobegin_idx is None or geoend_idx is None:
        raise ValueError("未找到 GEOBEGIN/GEOEND")

    # COMBNAME 模式：GEOBEGIN 行末含 COMBNAME，body/region 用名称
    combname = "COMBNAME" in lines[geobegin_idx].upper()

    geo_lines = lines[geobegin_idx + 1: geoend_idx]

    # ---- 2. 分离 body 段与 region 段 ----
    # body 段以第一个 "END" 结束（在 COMBNAME 模式下）
    # region 段以第二个 "END" 结束
    end_positions = [i for i, l in enumerate(geo_lines) if l.strip().upper() == "END"]
    if len(end_positions) < 2:
        raise ValueError(f"未找到两个 END 分隔符，找到 {len(end_positions)} 个")

    body_lines = geo_lines[:end_positions[0]]
    region_lines = geo_lines[end_positions[0] + 1: end_positions[1]]

    geom = Geometry()

    # ---- 3. 解析 body ----
    for line in body_lines:
        s = line.strip()
        if not s or s.startswith("*"):
            continue
        # 跳过 continuation 标记行（如 "    0    0"）
        # body 行: TYPE name params...
        m = BODY_TYPE_RE.match(line)
        if not m:
            continue
        body_type = m.group(1).upper()
        rest = line[m.end():].split()
        if not rest:
            continue
        name = rest[0]
        param_tokens = rest[1:]
        # 转 float
        try:
            params = [_to_float(t) for t in param_tokens]
        except ValueError:
            continue

        if body_type not in BODY_PARSERS:
            raise NotImplementedError(f"Body 类型 {body_type} 暂未实现（在 BODY_PARSERS 注册即可扩展）")
        body = BODY_PARSERS[body_type](name, params)
        geom.bodies[name] = body

    # ---- 4. 解析 region ----
    for line in region_lines:
        s = line.strip()
        if not s or s.startswith("*") or s.upper() == "END":
            continue
        # region 行: NAME [view] expr
        # FLUKA region 卡片可能含可选的 view 编号（整数），需跳过
        parts = s.split(None, 2)
        if len(parts) < 2:
            continue
        region_name = parts[0]
        # 若第二个 token 是纯整数，视为 view 编号并跳过
        if len(parts) >= 3 and parts[1].lstrip("-").isdigit():
            expr_str = parts[2]
        else:
            # 第二个 token 不是整数 -> 是表达式的一部分
            expr_str = parts[1]
            # 但 parts[1] 此时只含一个 token，需取完整剩余
            full = s.split(None, 1)[1]
            expr_str = full
        tokens = tokenize_expr(expr_str)
        parser = CSGParser(tokens)
        ast = parser.parse()
        geom.regions.append(Region(name=region_name, expr_tree=ast))

    # ---- 5. 解析 ASSIGNMA ----
    assignma_re = re.compile(r"^\s*ASSIGNMA\b", re.IGNORECASE)
    region_by_name = {r.name: r for r in geom.regions}
    for line in lines:
        if not assignma_re.match(line):
            continue
        # ASSIGNMA  material  region  [material_alias]
        fields = line.split()
        if len(fields) < 3:
            continue
        material = fields[1]
        region_name = fields[2]
        alias = fields[3] if len(fields) >= 4 else None
        if region_name in region_by_name:
            region_by_name[region_name].material = material
            region_by_name[region_name].material_alias = alias
            geom.material_to_regions.setdefault(material, []).append(region_name)

    return geom


# ============================================================
# 6. MaterialGrid 与切面
# ============================================================

@dataclass
class MaterialGrid:
    """切面上的材料分布网格。

    region_id[i,j] = 区域索引（1..N），0 = 无归属
    free_edges_1/2 与热图 h_edges/v_edges 对齐。
    """
    region_id: np.ndarray            # (n1, n2) int
    free_edges_1: np.ndarray         # (n1+1,) bin 边界
    free_edges_2: np.ndarray         # (n2+1,) bin 边界
    regions: list                    # list[Region]，索引 0 对应 region_id=1
    cut_axis: str                    # 'x'/'y'/'z'
    cut_value: float

    def material_mask(self, material: str) -> np.ndarray:
        """返回材料对应的 bool mask。"""
        mask = np.zeros_like(self.region_id, dtype=bool)
        for i, r in enumerate(self.regions, start=1):
            if r.material == material:
                mask |= (self.region_id == i)
        return mask

    def material_name_grid(self) -> np.ndarray:
        """返回 (n1, n2) 字符串数组（调试用）。"""
        names = np.empty(self.region_id.shape, dtype=object)
        names[:] = "VOID"
        for i, r in enumerate(self.regions, start=1):
            names[self.region_id == i] = r.material or "UNASSIGNED"
        return names

    def to_csv(self, path) -> None:
        """导出长格式 CSV: free1,free2,region_id,region_name,material_name"""
        path = Path(path)
        # 像素中心
        c1 = 0.5 * (self.free_edges_1[:-1] + self.free_edges_1[1:])
        c2 = 0.5 * (self.free_edges_2[:-1] + self.free_edges_2[1:])
        C1, C2 = np.meshgrid(c1, c2, indexing="ij")
        lines = ["free1,free2,region_id,region_name,material_name"]
        for i in range(self.region_id.shape[0]):
            for j in range(self.region_id.shape[1]):
                rid = int(self.region_id[i, j])
                if rid == 0:
                    rn, mn = "NONE", "VOID"
                else:
                    r = self.regions[rid - 1]
                    rn, mn = r.name, r.material or "UNASSIGNED"
                lines.append(f"{C1[i,j]:.6f},{C2[i,j]:.6f},{rid},{rn},{mn}")
        path.write_text("\n".join(lines), encoding="utf-8")


def cut_plane(
    geom: Geometry,
    axis: str,
    value: float,
    free_edges_1: np.ndarray,
    free_edges_2: np.ndarray,
) -> MaterialGrid:
    """在 axis=value 切面上计算材料分布。

    切面轴约定:
      axis='x' -> free1=y, free2=z
      axis='y' -> free1=x, free2=z
      axis='z' -> free1=x, free2=y
    free_edges_1/2 为 bin 边界（与热图对齐）。
    """
    # 像素中心
    c1 = 0.5 * (free_edges_1[:-1] + free_edges_1[1:])
    c2 = 0.5 * (free_edges_2[:-1] + free_edges_2[1:])
    C1, C2 = np.meshgrid(c1, c2, indexing="ij")

    if axis == "x":
        X = np.full_like(C1, value)
        Y, Z = C1, C2
    elif axis == "y":
        X = C1
        Y = np.full_like(C1, value)
        Z = C2
    elif axis == "z":
        X, Y = C1, C2
        Z = np.full_like(C1, value)
    else:
        raise ValueError(f"未知切面轴: {axis}")

    # 按列表顺序，先匹配先得
    region_id = np.zeros(C1.shape, dtype=np.int32)
    for i, region in enumerate(geom.regions, start=1):
        mask = region.inside(geom.bodies, X, Y, Z)
        region_id[mask & (region_id == 0)] = i

    return MaterialGrid(
        region_id=region_id,
        free_edges_1=free_edges_1,
        free_edges_2=free_edges_2,
        regions=geom.regions,
        cut_axis=axis,
        cut_value=value,
    )


# ============================================================
# 7. 边界绘制辅助
# ============================================================

def collect_analytic_boundaries(
    geom: Geometry,
    axis: str,
    value: float,
    skip_materials: Optional[set] = None,
) -> list[dict]:
    """收集所有 body 在切面上的解析边界曲线。

    用于在热图上叠加边界线。skip_materials 中的材料对应的 region 跳过。
    返回 [{"free1": ..., "free2": ..., "body": name, "region": name, "material": ...}, ...]
    """
    skip_materials = skip_materials or set()
    curves = []
    for region in geom.regions:
        if region.material in skip_materials:
            continue
        # 提取 region 引用的所有 body 名
        body_names = _collect_body_refs(region.expr_tree)
        for bn in body_names:
            if bn not in geom.bodies:
                continue
            body = geom.bodies[bn]
            cs = body.analytic_boundary_on_plane(axis, value)
            for c in cs:
                curves.append({
                    "free1": c["free1"],
                    "free2": c["free2"],
                    "body": bn,
                    "region": region.name,
                    "material": region.material,
                })
    return curves


def _collect_body_refs(node) -> set:
    """递归收集 AST 中引用的所有 body 名。"""
    refs = set()
    if isinstance(node, BodyRef):
        refs.add(node.name)
    elif isinstance(node, (AndNode, OrNode)):
        for c in node.children:
            refs |= _collect_body_refs(c)
    elif isinstance(node, NegateNode):
        refs |= _collect_body_refs(node.child)
    return refs


# ============================================================
# 模块自检（直接运行本文件时执行）
# ============================================================

if __name__ == "__main__":
    print("fluka_geometry.py 模块加载成功")
    print(f"已注册 body 类型: {list(BODY_PARSERS.keys())}")
