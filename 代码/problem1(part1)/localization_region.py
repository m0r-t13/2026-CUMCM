"""由示向度构造定位区域（半平面求交），作为旋转卡壳的输入。

纯 Python，无第三方依赖。

几何设定
--------
检测点 S_i 处测得示向度 theta_i，误差界 alpha（本题 alpha = 1 度）。则该次观测
给出的可行域是**角楔**

    W_i = { S_i + t*(cos phi, sin phi) : t >= 0, |phi - theta_i| <= alpha }

记 u = (cos theta_i, sin theta_i)、t = tan(alpha)、d = P - S_i，则

    P in W_i  <=>  d.u >= 0  且  |d x u| <= t * (d.u)

拆开即 **3 个半平面** a*x + b*y <= c：
  两条来自角楔的两条边界射线所在直线，第三条来自过顶点、垂直于中轴的
  "前方"半平面。**第三条不可省**：只用两条边界会同时保留反向的对顶楔。

定位区域 R = 交集 W_i，是凸多边形。其直径即问题 1 前半问所求。
"""

from __future__ import annotations

import math
from itertools import combinations

import rotating_calipers as rc

__all__ = [
    "wedge_halfplanes",
    "normals_positively_span",
    "intersect_halfplanes",
    "region_from_bearings",
    "disk_halfplanes",
    "disk_polygon",
    "clip_halfplane",
]


def wedge_halfplanes(sx: float, sy: float, theta_deg: float, alpha_deg: float = 1.0):
    """把一个角楔写成 3 个半平面 (a, b, c)，含义 a*x + b*y <= c，(a,b) 为外法向。"""
    th = math.radians(theta_deg)
    al = math.radians(alpha_deg)
    ux, uy = math.cos(th), math.sin(th)
    t = math.tan(al)

    cuS = ux * sy - uy * sx          # cross(u, S)
    duS = ux * sx + uy * sy          # dot(u, S)

    # (A)  cross(u,d) - t*(u.d) <= 0
    # (B) -cross(u,d) - t*(u.d) <= 0
    # (C) -(u.d) <= 0                  （只保留前方，排除对顶楔）
    return [
        (-uy - t * ux, ux - t * uy, cuS - t * duS),
        (uy - t * ux, -ux - t * uy, -cuS - t * duS),
        (-ux, -uy, -duS),
    ]


def normals_positively_span(hps) -> bool:
    """外法向是否正张成 R^2（等价于半平面交集有界）。

    判据：外法向按幅角排序后，相邻夹角的最大值 < pi。
    若存在 >= pi 的空隙，则所有法向落在某个闭半平面内，交集沿该方向无界。
    """
    angs = sorted(math.atan2(b, a) for a, b, _ in hps)
    n = len(angs)
    if n < 2:
        return False
    gaps = [angs[i + 1] - angs[i] for i in range(n - 1)]
    gaps.append(angs[0] + 2 * math.pi - angs[-1])
    return max(gaps) < math.pi - 1e-12


def _line_intersect(l1, l2):
    a1, b1, c1 = l1
    a2, b2, c2 = l2
    det = a1 * b2 - a2 * b1
    if abs(det) < 1e-14:
        return None
    return ((c1 * b2 - c2 * b1) / det, (a1 * c2 - a2 * c1) / det)


def intersect_halfplanes(hps, scale: float = 1.0):
    """半平面交集 {P : a.P <= c}。

    返回 ('ok', 顶点列表) / ('empty', []) / ('unbounded', [])。

    做法（顶点枚举）：有界多胞体的每个顶点都是两条边界直线的交点，故枚举全部
    线对交点、保留满足所有半平面者，再求凸包即可。该方法无需人为包围盒。
    """
    if not hps:
        return "unbounded", []
    if not normals_positively_span(hps):
        return "unbounded", []

    tol = 1e-9 * max(1.0, scale)
    cand = []
    for l1, l2 in combinations(hps, 2):
        p = _line_intersect(l1, l2)
        if p is None:
            continue
        if all(a * p[0] + b * p[1] - c <= tol for a, b, c in hps):
            cand.append(p)
    if not cand:
        return "empty", []
    poly = rc.simplify_polygon(rc.convex_hull(cand))
    return ("ok", poly) if poly else ("empty", [])


def disk_halfplanes(cx: float, cy: float, radius: float, n_sides: int = 720):
    """用内接正 n 边形近似圆盘，返回其半平面（供可选约束使用）。

    径向亏空 = radius*(1 - cos(pi/n))；r=1800、n=720 时约 0.017 m。
    """
    out = []
    for i in range(n_sides):
        ang = 2 * math.pi * i / n_sides
        a, b = math.cos(ang), math.sin(ang)
        out.append((a, b, a * cx + b * cy + radius * math.cos(math.pi / n_sides)))
    return out


def disk_polygon(cx: float, cy: float, radius: float, n_sides: int = 720):
    """圆盘的内接正 n 边形顶点（逆时针）。"""
    return [
        (cx + radius * math.cos(2 * math.pi * i / n_sides),
         cy + radius * math.sin(2 * math.pi * i / n_sides))
        for i in range(n_sides)
    ]


def clip_halfplane(poly, hp, tol: float = 1e-9):
    """用单个半平面 a.x + b.y <= c 裁剪多边形（Sutherland-Hodgman）。"""
    a, b, c = hp
    n = len(poly)
    if n == 0:
        return []
    out = []
    for i in range(n):
        p = poly[i]
        q = poly[(i + 1) % n]
        dp = a * p[0] + b * p[1] - c
        dq = a * q[0] + b * q[1] - c
        pin, qin = dp <= tol, dq <= tol
        if pin:
            out.append(p)
        if pin != qin:
            t = dp / (dp - dq)
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return out


def region_from_bearings(
    detections,
    alpha_deg: float = 1.0,
    region_disk_radius=None,
    reception_max=None,
    n_sides: int = 720,
):
    """由检测点与示向度构造定位区域。

    detections : [(x, y, bearing_deg), ...]

    可选约束（问题 1 基础模型不需要，问题 3/4 中有用）：
      region_disk_radius : 叠加"干扰源位于以原点为心、该半径的圆域内"
      reception_max      : 叠加"干扰源到各检测点距离 <= reception_max"

    返回 (status, polygon, notes)，polygon 为逆时针严格凸顶点列表。

    实现分两条路：
      * 无约束（问题 1 的基础模型）——顶点枚举法，O(k^3)，精确且无需包围盒；
      * 有圆盘约束——从圆盘多边形出发逐条裁剪，避免 720*k 条半平面的枚举开销。
    """
    notes = []
    scale = 1.0
    for x, y, _ in detections:
        scale = max(scale, abs(x), abs(y))

    wedge_hps = []
    for x, y, th in detections:
        wedge_hps.extend(wedge_halfplanes(x, y, th, alpha_deg))

    if not region_disk_radius and not reception_max:
        status, poly = intersect_halfplanes(wedge_hps, scale)
        if status == "unbounded":
            notes.append(
                "角楔交集无界：各示向度方向锥存在公共方向，没有有限的定位区域直径。"
                "需增加一个不同方位的检测点，或叠加圆域/接收距离约束。"
            )
        return status, poly, notes

    # ---- 有约束：从圆盘多边形出发裁剪 ----
    seed = None
    if region_disk_radius:
        seed = disk_polygon(0.0, 0.0, region_disk_radius, n_sides)
        notes.append(
            f"叠加圆域约束 r <= {region_disk_radius:g} m"
            f"（内接 {n_sides} 边形，径向亏空 "
            f"{region_disk_radius * (1 - math.cos(math.pi / n_sides)):.4f} m）"
        )
    if reception_max:
        for x, y, _ in detections:
            disk = disk_polygon(x, y, reception_max, n_sides)
            seed = disk if seed is None else _clip_by_convex(seed, disk, scale)
            if not seed:
                break
        notes.append(f"叠加接收距离约束 <= {reception_max:g} m")

    if not seed:
        return "empty", [], notes

    tol = 1e-9 * max(1.0, scale)
    poly = seed
    for hp in wedge_hps:
        poly = clip_halfplane(poly, hp, tol)
        if not poly:
            break
    if not poly:
        return "empty", [], notes
    poly = rc.simplify_polygon(rc.convex_hull(poly))
    return ("ok", poly, notes) if poly else ("empty", [], notes)


def _clip_by_convex(subject, clipper, scale: float = 1.0):
    """用逆时针凸多边形 clipper 裁剪 subject。"""
    tol = 1e-9 * max(1.0, scale)
    clip = rc.rotate_to_ccw(list(clipper))
    m = len(clip)
    if m < 3:
        return []
    poly = list(subject)
    for i in range(m):
        p, q = clip[i], clip[(i + 1) % m]
        ex, ey = q[0] - p[0], q[1] - p[1]
        a, b = ey, -ex                      # 逆时针多边形的外法向
        c = a * p[0] + b * p[1]
        poly = clip_halfplane(poly, (a, b, c), tol)
        if not poly:
            break
    return poly
