#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2026 CUMCM B题 - 问题3 机器狗自动搜索/定位/清除程序

策略：
1) 自适应7站保证搜索：中心固定，外圈站允许在“仍可由剩余预算完成1000m全覆盖”的前提下向高价值定位点偏移；
   若自适应偏移不能保持覆盖证书，则自动回退到原六边形骨架；
2) 搜索阶段只扫描尚未发现的频道；
3) 对已发现源维护保守可行多边形（场地外接多边形 ∩ 每次±1°角楔 ∩ 1500m外接圆盘）；
4) 在搜索骨架点顺便进行有价值的二次/多次测向；
5) 搜索完成后，用“几何指标粗筛 + 角域交会预测精筛”选择下一测点；
6) 当保守可行域的最小包围圆半径 <= READY_RADIUS 时，到圆心 /clear；
7) 全程串行 HTTP；网络故障重试同一动作时复用同一个 request_id；记录 JSONL 日志。

仅针对问题3（全部为全向干扰源）。
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import random
import socket
import sys
import time
import uuid
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

Point = Tuple[float, float]

# ----------------------------
# 可调参数
# ----------------------------
ARENA_R = 1800.0
MIN_RECEIVE_R = 1000.0
MAX_RECEIVE_R = 1500.0
BEARING_ERR_DEG = 1.0
NEAR_R = 5.0
CLEAR_R = 20.0
READY_RADIUS = 19.0          # 留1m数值安全余量
COVER_RING_R = 1140.0        # 7点覆盖：中心 + 外圈6点；1140m仍有约7m严格覆盖裕量，路线更短
CIRCLE_SIDES = 180           # 外接多边形近似圆；始终是圆的超集，保证保守性
COARSE_STEP = 125.0          # 主动定位粗网格
FINE_STEP = 30.0             # 主动定位局部精网格
GEOM_KEEP_RATIO = 0.35       # 几何粗筛保留比例
NEAR_OPT_RATIO = 1.15        # 精确预测近优集合；在其中选最顺路点
MAX_LOCALIZE_MEASURES = 7
MIN_ACTIVE_BASELINE = 80.0    # 主动新增测点尽量与既有测点拉开，避免同地重复测向
STATION_PROBE_IMPROVEMENT = 0.78  # 搜索站顺便测已知源时，预测直径至少缩到当前的78%
MAX_REAL_MARGIN_S = 20.0

# 自适应保证搜索：最多仍使用7个“扫描未知频道”的站点。
# 由于1000/1800略低于6圆覆盖阈值，若仍有未知频道，严格全覆盖至少需要7站；
# 因此这里优化的是站位/顺序与定位复用，而不是无条件把7站压成5~6站。
MAX_SEARCH_STATIONS = 7
ADAPTIVE_SOURCE_LIMIT = 4
ADAPTIVE_SHIFT_FACTORS = (0.05, 0.10, 0.15, 0.20)
ADAPTIVE_COVER_TOL = 1e-6

# 自适应站评分（量纲统一为“等效米”）：低分优先。
# 覆盖计划剩余站数优先级最高，其次兼顾当前移动、覆盖半径下降和主频道定位收益。
ADAPTIVE_REMAINING_STATION_PENALTY = 1800.0
ADAPTIVE_COVER_RADIUS_WEIGHT = 0.55
ADAPTIVE_LOCALIZE_BONUS = 650.0
ADAPTIVE_JOINT_BONUS = 90.0
ADAPTIVE_CLEAR_BONUS = 1200.0      # 已就绪源中心若可兼作搜索站，可省掉未来一次专程访问
ADAPTIVE_READY_CLEAR_LIMIT = 6

# ---- 虚拟时间模型（附件2默认参数） ----
# 目标函数应直接按“秒”优化，而不是只追求几何精度。移动 500m = 100s，远大于一次测向5s。
MOVE_SPEED_MPS = 5.0
MEASURE_ACTION_S = 5.0
SWITCH_CHANNEL_S = 1.0
CLEAR_ACTION_S = 5.0

# 主动定位的时间型目标：允许牺牲少量单次几何精度，换取显著更短的移动。
TIME_RESIDUAL_WEIGHT = 0.85
TIME_TO_CENTER_WEIGHT = 0.38
TIME_QUALITY_CAP_EARLY = 1.50   # 前两次观测：允许比几何最优直径差50%，换取少移动
TIME_QUALITY_CAP_LATE = 1.18    # 3次及以后：必须接近几何最优，防止反复低价值测向

# 到达任意定位点/清除点后，顺手给其他源测向。一次顺手测量最多约6s，
# 只要能避免后续几十/几百米的专程移动就非常划算。
JOINT_PROBE_RATIO_FIRST = 0.88
JOINT_PROBE_RATIO_LATER = 0.74
MAX_JOINT_PROBES_PER_STOP = 5
JOINT_GAIN_VALUE_S = 95.0
TASK_LOOKAHEAD_WEIGHT = 0.18

# 搜索阶段若一个已就绪源的清除点几乎位于当前->下一搜索站的路上，就顺路清掉。
SEARCH_CLEAR_MAX_DETOUR_M = 220.0
# 主动定位后，若清除点离当前位置不远，立即清除；否则交给全局动作调度，避免折返。
IMMEDIATE_CLEAR_DISTANCE_M = 260.0
AUX_SCAN_MIN_SPACING_M = 1200.0  # 回归更优：更积极让定位/清除点替代后备覆盖站

# ---- 第二轮路线/局部清除优化 ----
REGION_ROUTE_2OPT_PASSES = 8
EARLY_CENTER_PROGRESS_WEIGHT = 1.05   # 前两次测向更偏向“朝源区域前进”，减少V形折返
LOCAL_CLEAR_MAX_POINTS = 6            # 小可行域至多用6个保证清除点替代下一次远距离补测
LOCAL_CLEAR_CELL_RADIUS = 18.5        # 每个网格单元半对角线<=18.5m，给20m清除半径留余量
LOCAL_CLEAR_ANGLES_DEG = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0)
LOCAL_CLEAR_FUTURE_MARGIN_S = 45.0   # 扫清成本略高于一次测向也可接受，因为它直接完成该源

EPS = 1e-9


# ----------------------------
# 基础几何
# ----------------------------
def add(a: Point, b: Point) -> Point:
    return (a[0] + b[0], a[1] + b[1])


def sub(a: Point, b: Point) -> Point:
    return (a[0] - b[0], a[1] - b[1])


def mul(a: Point, k: float) -> Point:
    return (a[0] * k, a[1] * k)


def dot(a: Point, b: Point) -> float:
    return a[0] * b[0] + a[1] * b[1]


def cross(a: Point, b: Point) -> float:
    return a[0] * b[1] - a[1] * b[0]


def dist(a: Point, b: Point) -> float:
    return math.hypot(a[0] - b[0], a[1] - b[1])


def centroid(poly: Sequence[Point]) -> Point:
    if not poly:
        return (0.0, 0.0)
    # 对凸多边形，用顶点均值做稳定的搜索中心即可
    return (sum(p[0] for p in poly) / len(poly), sum(p[1] for p in poly) / len(poly))


def bearing_deg(a: Point, b: Point) -> float:
    return math.degrees(math.atan2(b[1] - a[1], b[0] - a[0])) % 360.0


def angle_diff_deg(a: float, b: float) -> float:
    d = (a - b + 180.0) % 360.0 - 180.0
    return abs(d)


def regular_circumscribed_polygon(center: Point, radius: float, n: int) -> List[Point]:
    """返回包含真实圆盘的正n边形（内切圆半径=radius）。"""
    rc = radius / math.cos(math.pi / n)
    return [
        (center[0] + rc * math.cos(2 * math.pi * k / n),
         center[1] + rc * math.sin(2 * math.pi * k / n))
        for k in range(n)
    ]


def clip_halfplane(poly: Sequence[Point], f) -> List[Point]:
    """Sutherland-Hodgman；保留 f(p)>=0 的半平面。"""
    if not poly:
        return []
    out: List[Point] = []
    prev = poly[-1]
    fp = f(prev)
    prev_in = fp >= -1e-10
    for cur in poly:
        fc = f(cur)
        cur_in = fc >= -1e-10
        if cur_in != prev_in:
            den = fp - fc
            if abs(den) > 1e-15:
                t = fp / den
                inter = (prev[0] + t * (cur[0] - prev[0]),
                         prev[1] + t * (cur[1] - prev[1]))
                out.append(inter)
        if cur_in:
            out.append(cur)
        prev, fp, prev_in = cur, fc, cur_in
    return dedup_polygon(out)


def dedup_polygon(poly: Sequence[Point]) -> List[Point]:
    if not poly:
        return []
    out: List[Point] = []
    for p in poly:
        if not out or dist(p, out[-1]) > 1e-7:
            out.append(p)
    if len(out) >= 2 and dist(out[0], out[-1]) <= 1e-7:
        out.pop()
    return out


def clip_wedge(poly: Sequence[Point], station: Point, measured_deg: float, alpha_deg: float = BEARING_ERR_DEG) -> List[Point]:
    """保留从station看过去方位落在 measured±alpha 内的点。"""
    lo = math.radians(measured_deg - alpha_deg)
    hi = math.radians(measured_deg + alpha_deg)
    ulo = (math.cos(lo), math.sin(lo))
    uhi = (math.cos(hi), math.sin(hi))

    p1 = clip_halfplane(poly, lambda p: cross(ulo, sub(p, station)))
    p2 = clip_halfplane(p1, lambda p: -cross(uhi, sub(p, station)))
    return p2


def clip_disk_outer(poly: Sequence[Point], center: Point, radius: float, n: int = CIRCLE_SIDES) -> List[Point]:
    """用圆的外接正n边形约束：真实圆盘完全包含在结果允许区域中（保守）。"""
    out = list(poly)
    for k in range(n):
        a = 2 * math.pi * k / n
        normal = (math.cos(a), math.sin(a))
        out = clip_halfplane(out, lambda p, nrm=normal: radius - dot(nrm, sub(p, center)))
        if not out:
            break
    return out


def polygon_diameter(poly: Sequence[Point]) -> float:
    if len(poly) <= 1:
        return 0.0
    best = 0.0
    for i in range(len(poly)):
        for j in range(i + 1, len(poly)):
            best = max(best, dist(poly[i], poly[j]))
    return best


def sample_poly(poly: Sequence[Point], max_vertices: int = 16) -> List[Point]:
    if not poly:
        return []
    pts: List[Point] = []
    n = len(poly)
    if n <= max_vertices:
        pts.extend(poly)
    else:
        for i in range(max_vertices):
            pts.append(poly[(i * n) // max_vertices])
    pts.append(centroid(poly))
    # 加一些边中点，减小只看顶点的偏差
    m = min(8, n)
    for i in range(m):
        a = poly[(i * n) // m]
        b = poly[((i * n) // m + 1) % n]
        pts.append(((a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5))
    return pts


def robust_receive(point: Point, poly: Sequence[Point], radius: float = MIN_RECEIVE_R) -> bool:
    """凸集到point的最远距离在顶点取到；若所有顶点<=1000，则真实源必能收到。"""
    return bool(poly) and max(dist(point, v) for v in poly) <= radius + 1e-7


# ----------------------------
# 最小包围圆（Nayuki/Welzl风格的随机增量实现）
# ----------------------------
@dataclass
class Circle:
    c: Point
    r: float


def is_in_circle(c: Circle, p: Point) -> bool:
    return dist(c.c, p) <= c.r + 1e-7


def circle_two(a: Point, b: Point) -> Circle:
    c = ((a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0)
    return Circle(c, dist(a, b) / 2.0)


def circle_three(a: Point, b: Point, c: Point) -> Optional[Circle]:
    ax, ay = a
    bx, by = b
    cx, cy = c
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if abs(d) < 1e-12:
        return None
    ux = ((ax * ax + ay * ay) * (by - cy) +
          (bx * bx + by * by) * (cy - ay) +
          (cx * cx + cy * cy) * (ay - by)) / d
    uy = ((ax * ax + ay * ay) * (cx - bx) +
          (bx * bx + by * by) * (ax - cx) +
          (cx * cx + cy * cy) * (bx - ax)) / d
    cc = (ux, uy)
    return Circle(cc, dist(cc, a))


def minimum_enclosing_circle(points: Sequence[Point]) -> Circle:
    pts = list(points)
    if not pts:
        return Circle((0.0, 0.0), 0.0)
    rnd = random.Random(20260911)
    rnd.shuffle(pts)
    c: Optional[Circle] = None
    for i, p in enumerate(pts):
        if c is None or not is_in_circle(c, p):
            c = Circle(p, 0.0)
            for j in range(i):
                q = pts[j]
                if not is_in_circle(c, q):
                    c = circle_two(p, q)
                    for k in range(j):
                        r = pts[k]
                        if not is_in_circle(c, r):
                            cc = circle_three(p, q, r)
                            if cc is not None:
                                c = cc
                            else:
                                # 共线：取三对中能覆盖三点且半径最小者
                                opts = [circle_two(p, q), circle_two(p, r), circle_two(q, r)]
                                opts = [o for o in opts if all(is_in_circle(o, z) for z in (p, q, r))]
                                c = min(opts, key=lambda o: o.r)
    assert c is not None
    return c


# ----------------------------
# 源状态
# ----------------------------
@dataclass
class Observation:
    pos: Point
    svd_deg: float


@dataclass
class SourceState:
    channel: int
    poly: List[Point]
    observations: List[Observation] = field(default_factory=list)
    cleared: bool = False
    clear_attempts: int = 0

    def add_direction(self, pos: Point, svd_deg: float) -> None:
        p = clip_wedge(self.poly, pos, svd_deg)
        # direction 返回意味着真实距离>5m。对 ±1° 角楔内任一点，
        # 沿测得中心方向的投影至少应超过 5*cos(1°)。这是线性、保守且凸的约束。
        a = math.radians(svd_deg)
        u = (math.cos(a), math.sin(a))
        min_forward = NEAR_R * math.cos(math.radians(BEARING_ERR_DEG))
        p = clip_halfplane(p, lambda x: dot(u, sub(x, pos)) - min_forward)
        p = clip_disk_outer(p, pos, MAX_RECEIVE_R)
        if not p:
            raise RuntimeError(f"频道{self.channel}加入观测后可行域为空；请检查示向角/实现。")
        self.poly = p
        self.observations.append(Observation(pos, svd_deg))

    @property
    def diameter(self) -> float:
        return polygon_diameter(self.poly)

    @property
    def mec(self) -> Circle:
        return minimum_enclosing_circle(self.poly)

    @property
    def ready(self) -> bool:
        return bool(self.poly) and self.mec.r <= READY_RADIUS


# ----------------------------
# 主动选点：几何粗筛 + 角域预测精筛
# ----------------------------
def candidate_bbox(poly: Sequence[Point], radius: float = MIN_RECEIVE_R) -> Optional[Tuple[float, float, float, float]]:
    if not poly:
        return None
    xmin = max(v[0] - radius for v in poly)
    xmax = min(v[0] + radius for v in poly)
    ymin = max(v[1] - radius for v in poly)
    ymax = min(v[1] + radius for v in poly)
    if xmin > xmax or ymin > ymax:
        return None
    return xmin, xmax, ymin, ymax


def grid_points(bounds: Tuple[float, float, float, float], step: float) -> List[Point]:
    xmin, xmax, ymin, ymax = bounds
    pts: List[Point] = []
    # 保证即使区间很窄也含中心
    nx = max(1, int(math.floor((xmax - xmin) / step)) + 1)
    ny = max(1, int(math.floor((ymax - ymin) / step)) + 1)
    for ix in range(nx + 1):
        x = min(xmax, xmin + ix * step)
        for iy in range(ny + 1):
            y = min(ymax, ymin + iy * step)
            pts.append((x, y))
    pts.append(((xmin + xmax) / 2.0, (ymin + ymax) / 2.0))
    # 去重
    seen = set()
    out = []
    for p in pts:
        key = (round(p[0], 6), round(p[1], 6))
        if key not in seen:
            seen.add(key)
            out.append(p)
    return out


def geometric_score(q: Point, state: SourceState) -> float:
    """问题2中的U型指标用于粗筛。以第一次观测点为参考站。"""
    if not state.observations:
        return float("inf")
    s1 = state.observations[0].pos
    worst = 0.0
    for g in sample_poly(state.poly, 12):
        d1 = max(dist(s1, g), 1e-6)
        d2 = max(dist(q, g), 1e-6)
        v1 = sub(s1, g)
        v2 = sub(q, g)
        den = d1 * d2
        sin_g = abs(cross(v1, v2)) / max(den, 1e-12)
        u = d1 * d2 / max(sin_g, 1e-4)
        worst = max(worst, u)
    return worst


def predicted_worst_diameter(q: Point, state: SourceState) -> float:
    """
    预测若在q测一次后最坏定位直径。
    对可能真值G抽样，并考虑新的测向误差中心相对真实方位偏 -1/0/+1 度。
    这是选点启发式；实际定位区域仍由真实返回值精确更新，因此不会破坏最终保守性。
    """
    worst = 0.0
    for g in sample_poly(state.poly, 14):
        if dist(q, g) <= NEAR_R:
            d_here = 0.0
        else:
            true_b = bearing_deg(q, g)
            d_here = 0.0
            for e in (-BEARING_ERR_DEG, 0.0, BEARING_ERR_DEG):
                p = clip_wedge(state.poly, q, true_b + e)
                # q 已通过 robust_receive：当前整个保守可行域距 q 均 <=1000m，
                # 因而 1500m 接收上界不会进一步收紧预测区域，无需重复裁圆。
                if not p:
                    continue
                d_here = max(d_here, polygon_diameter(p))
        worst = max(worst, d_here)
    return worst


def _time_score_measure_candidate(q: Point, predicted_d: float, state: SourceState, current_pos: Point) -> float:
    """按官方虚拟时间模型给主动测点打分（单位：秒，越小越好）。

    核心思想：移动耗时 = 距离/5。几何上“最漂亮”的测点如果多走800m，
    就会多花160s；相比之下，多做一次测向通常只要5~6s。
    因此这里把移动、测后到当前MEC中心的距离、以及残余不确定性统一换算成秒。
    """
    move_s = dist(current_pos, q) / MOVE_SPEED_MPS
    center_leg_s = dist(q, state.mec.c) / MOVE_SPEED_MPS
    residual_s = predicted_d / MOVE_SPEED_MPS
    center_weight = (EARLY_CENTER_PROGRESS_WEIGHT
                     if len(state.observations) <= 2 else TIME_TO_CENTER_WEIGHT)
    return (move_s + MEASURE_ACTION_S
            + center_weight * center_leg_s
            + TIME_RESIDUAL_WEIGHT * residual_s)


def _filter_measure_candidates(points: Sequence[Point], state: SourceState) -> List[Point]:
    """保留保证能收到信号、且与旧测点有足够基线的候选点。"""
    raw = [p for p in points if robust_receive(p, state.poly)]
    good = [p for p in raw
            if all(dist(p, obs.pos) >= MIN_ACTIVE_BASELINE for obs in state.observations)]
    if good:
        return good
    return [p for p in raw if all(dist(p, obs.pos) >= 1.0 for obs in state.observations)]


def _append_unique_point(dst: List[Point], p: Point) -> None:
    key = (round(p[0], 6), round(p[1], 6))
    if not any((round(q[0], 6), round(q[1], 6)) == key for q in dst):
        dst.append(p)


def _early_progress_candidates(state: SourceState, current_pos: Point) -> List[Point]:
    """前两次观测时显式加入“朝MEC中心推进”的候选点。

    这些点仍必须经过 robust_receive 与基线检查；因此只是路线偏好，不改变保证性。
    """
    if len(state.observations) > 2 or not state.poly:
        return []
    c = state.mec.c
    pts: List[Point] = []
    for t in (0.35, 0.55, 0.75, 1.0):
        q = blend_point(current_pos, c, t)
        if inside_arena(q) and robust_receive(q, state.poly):
            if all(dist(q, obs.pos) >= MIN_ACTIVE_BASELINE for obs in state.observations):
                _append_unique_point(pts, q)
    return pts


def choose_next_measure_point(state: SourceState, current_pos: Point) -> Point:
    """选择“预计虚拟时间最短”的下一测点，而不是纯几何最优点。

    仍然只在 robust_receive=True 的点中选择，因此不会牺牲全向源的保证接收性质。
    """
    bounds = candidate_bbox(state.poly)
    if bounds is None:
        return state.mec.c

    coarse = _filter_measure_candidates(grid_points(bounds, COARSE_STEP), state)
    if not coarse:
        c = state.mec
        if c.r <= MIN_RECEIVE_R + 1e-6:
            return c.c
        raise RuntimeError(f"频道{state.channel}找不到保证接收的下一检测点")

    # 几何粗筛仍保留，但不再只允许best_d的1.15倍以内；最终直接按预计秒数比较。
    scored = sorted((geometric_score(p, state), p) for p in coarse)
    keep_n = max(10, int(len(scored) * max(GEOM_KEEP_RATIO, 0.45)))
    keep = [p for _, p in scored[:keep_n]]

    # 强制加入“离当前位置最近”“离MEC中心最近”的候选，防止几何粗筛把省路点全部删掉。
    keep.extend(sorted(coarse, key=lambda p: dist(current_pos, p))[:6])
    keep.extend(sorted(coarse, key=lambda p: dist(state.mec.c, p))[:6])
    keep.extend(_early_progress_candidates(state, current_pos))
    uniq: Dict[Tuple[int, int], Point] = {}
    for p in keep:
        uniq[(round(p[0] * 1000), round(p[1] * 1000))] = p
    keep = list(uniq.values())

    exact = []
    for p in keep:
        pd = predicted_worst_diameter(p, state)
        exact.append((_time_score_measure_candidate(p, pd, state, current_pos), pd, p))
    best_pd = min(x[1] for x in exact)
    quality_cap = TIME_QUALITY_CAP_EARLY if len(state.observations) <= 2 else TIME_QUALITY_CAP_LATE
    admissible = [x for x in exact if x[1] <= best_pd * quality_cap + 1e-8]
    admissible.sort(key=lambda x: (x[0], x[1]))
    coarse_best = admissible[0][2]

    # 围绕时间最优粗点局部精化。
    bx = (coarse_best[0] - 2 * COARSE_STEP, coarse_best[0] + 2 * COARSE_STEP,
          coarse_best[1] - 2 * COARSE_STEP, coarse_best[1] + 2 * COARSE_STEP)
    xmin, xmax, ymin, ymax = bounds
    fine_bounds = (max(xmin, bx[0]), min(xmax, bx[1]), max(ymin, bx[2]), min(ymax, bx[3]))
    fine = _filter_measure_candidates(grid_points(fine_bounds, FINE_STEP), state)
    if not fine:
        return coarse_best

    fs = sorted((geometric_score(p, state), p) for p in fine)
    fkeep = [p for _, p in fs[:max(14, int(len(fs) * 0.55))]]
    fkeep.extend(sorted(fine, key=lambda p: dist(current_pos, p))[:8])
    fkeep.extend(sorted(fine, key=lambda p: dist(state.mec.c, p))[:8])
    fkeep.extend(_early_progress_candidates(state, current_pos))
    uniq.clear()
    for p in fkeep:
        uniq[(round(p[0] * 1000), round(p[1] * 1000))] = p

    fitems = []
    for p in uniq.values():
        pd = predicted_worst_diameter(p, state)
        fitems.append((_time_score_measure_candidate(p, pd, state, current_pos), pd, p))
    if not fitems:
        return coarse_best
    fbest_pd = min(x[1] for x in fitems)
    quality_cap = TIME_QUALITY_CAP_EARLY if len(state.observations) <= 2 else TIME_QUALITY_CAP_LATE
    fadmissible = [x for x in fitems if x[1] <= fbest_pd * quality_cap + 1e-8]
    fadmissible.sort(key=lambda x: (x[0], x[1]))
    return fadmissible[0][2]


# ----------------------------
# 7点覆盖
# ----------------------------
def coverage_points(phi_deg: float = 0.0) -> List[Point]:
    pts: List[Point] = [(0.0, 0.0)]
    for k in range(6):
        a = math.radians(phi_deg + 60.0 * k)
        pts.append((COVER_RING_R * math.cos(a), COVER_RING_R * math.sin(a)))
    return pts


def verify_coverage(points: Sequence[Point], samples_r: int = 90, samples_a: int = 720) -> float:
    worst = 0.0
    for ir in range(samples_r + 1):
        r = ARENA_R * ir / samples_r
        for ia in range(samples_a):
            a = 2 * math.pi * ia / samples_a
            p = (r * math.cos(a), r * math.sin(a))
            nearest = min(dist(p, q) for q in points)
            worst = max(worst, nearest)
    return worst


def guaranteed_cover_radius(points: Sequence[Point]) -> Tuple[float, Point]:
    """计算圆形场地内“到最近搜索站的最大距离”。

    对有限站点的最近距离函数，其最大值只可能出现在：
      1) 场地内部的Voronoi顶点（三站等距点）；
      2) Voronoi边与场地圆周的交点（两站等距点）；
      3) 某站Voronoi圆弧上的距离极大点（该站反方向的圆周点）。
    枚举这些候选即可得到圆域上的覆盖半径；额外加入圆心/四个轴向边界点处理退化情形。

    返回 (最大最近站距离, 对应最坏点)。该值<=MIN_RECEIVE_R 时，
    才允许把剩余未知频道严格判空。
    """
    pts = list(points)
    if not pts:
        return float("inf"), (0.0, 0.0)

    candidates: List[Point] = [
        (0.0, 0.0),
        (ARENA_R, 0.0), (-ARENA_R, 0.0),
        (0.0, ARENA_R), (0.0, -ARENA_R),
    ]

    # 圆周上相对每个站最远的方向；若该点属于该站Voronoi圆弧，
    # 它就是该圆弧内部可能出现的极大点。
    for sx, sy in pts:
        rr = math.hypot(sx, sy)
        if rr > 1e-12:
            candidates.append((-ARENA_R * sx / rr, -ARENA_R * sy / rr))

    # 两站等距线（垂直平分线）与场地圆的交点。
    for a, b in itertools.combinations(pts, 2):
        nx = b[0] - a[0]
        ny = b[1] - a[1]
        nn = nx * nx + ny * ny
        if nn <= 1e-18:
            continue
        rhs = (b[0] * b[0] + b[1] * b[1] - a[0] * a[0] - a[1] * a[1]) * 0.5
        px = nx * rhs / nn
        py = ny * rhs / nn
        p2 = px * px + py * py
        if p2 <= ARENA_R * ARENA_R + 1e-7:
            h = math.sqrt(max(0.0, ARENA_R * ARENA_R - p2))
            invn = 1.0 / math.sqrt(nn)
            tx, ty = -ny * invn, nx * invn
            candidates.append((px + h * tx, py + h * ty))
            candidates.append((px - h * tx, py - h * ty))

    # 三站对应的Voronoi顶点候选（外接圆圆心）。
    for a, b, c in itertools.combinations(pts, 3):
        cc = circle_three(a, b, c)
        if cc is None:
            continue
        x, y = cc.c
        if x * x + y * y <= ARENA_R * ARENA_R + 1e-7:
            candidates.append(cc.c)

    worst = -1.0
    worst_p = (0.0, 0.0)
    for p in candidates:
        if p[0] * p[0] + p[1] * p[1] > ARENA_R * ARENA_R + 1e-6:
            continue
        nearest = min(dist(p, q) for q in pts)
        if nearest > worst:
            worst = nearest
            worst_p = p
    return worst, worst_p


def minimum_completion_subset(visited: Sequence[Point], remaining: Sequence[Tuple[int, Point]],
                              max_additional: int) -> Optional[List[Tuple[int, Point]]]:
    """从剩余保证骨架点中找最少的补点集合，使 visited ∪ subset 严格覆盖整个场地。

    只枚举最多6个外圈骨架点，因此2^6规模很小。返回None表示在给定剩余站数预算内
    无法由当前骨架完成覆盖；自适应候选遇到这种情况必须拒绝。
    """
    rem = list(remaining)
    cap = min(max_additional, len(rem))
    for k in range(cap + 1):
        for idxs in itertools.combinations(range(len(rem)), k):
            subset = [rem[i] for i in idxs]
            radius, _ = guaranteed_cover_radius(list(visited) + [p for _, p in subset])
            if radius <= MIN_RECEIVE_R + ADAPTIVE_COVER_TOL:
                return subset
    return None


def inside_arena(p: Point, margin: float = 1e-7) -> bool:
    return p[0] * p[0] + p[1] * p[1] <= (ARENA_R + margin) ** 2


def blend_point(a: Point, b: Point, t: float) -> Point:
    return (a[0] * (1.0 - t) + b[0] * t,
            a[1] * (1.0 - t) + b[1] * t)


# ----------------------------
# HTTP客户端
# ----------------------------
class SimulatorClient:
    def __init__(self, base_url: str, robot_id: str, log_path: Path, timeout_s: float = 5.0, max_retries: int = 5):
        self.base_url = base_url.rstrip("/")
        self.robot_id = robot_id
        self.timeout_s = timeout_s
        self.max_retries = max_retries
        self.log_path = log_path
        self.counter = 0
        self.current_pos: Point = (0.0, 0.0)
        self.current_channel = 1
        self.last_virtual_time = 0.0
        self.remaining_real_duration_s: Optional[float] = None
        self.real_start_monotonic: Optional[float] = None
        # 仅用于演练统计，不参与控制逻辑。
        self.total_move_distance_m = 0.0
        self.measure_count = 0
        self.channel_switch_count = 0
        self.clear_count = 0
        self.clear_success_count = 0
        self.clear_fail_count = 0

    def _new_id(self, prefix: str) -> str:
        self.counter += 1
        return f"{prefix}-{self.counter:05d}-{uuid.uuid4().hex[:8]}"

    def _base(self, rid: str) -> dict:
        return {"arena_id": "default", "robot_id": self.robot_id, "request_id": rid}

    def _log(self, obj: dict) -> None:
        obj = dict(obj)
        obj["local_time"] = datetime.now().isoformat(timespec="milliseconds")
        with self.log_path.open("a", encoding="utf-8") as f:
            f.write(json.dumps(obj, ensure_ascii=False) + "\n")

    def real_time_left(self) -> float:
        if self.remaining_real_duration_s is None or self.real_start_monotonic is None:
            return float("inf")
        return self.remaining_real_duration_s - (time.monotonic() - self.real_start_monotonic)

    def _post(self, path: str, payload: dict) -> dict:
        data = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
        last_err = None
        for attempt in range(1, self.max_retries + 1):
            req = Request(
                self.base_url + path,
                data=data,
                headers={"Content-Type": "application/json"},
                method="POST",
            )
            try:
                with urlopen(req, timeout=self.timeout_s) as resp:
                    status = resp.getcode()
                    raw = resp.read().decode("utf-8")
                    body = json.loads(raw)
                self._log({"event": "http", "path": path, "attempt": attempt,
                           "request": payload, "http_status": status, "response": body})
                if status != 200:
                    raise RuntimeError(f"HTTP {status}: {body}")
                if body.get("accepted") is not True:
                    raise RuntimeError(f"动作未被模拟器接受: {path} {body}")
                if "virtual_time_s" in body:
                    self.last_virtual_time = float(body["virtual_time_s"])
                return body
            except HTTPError as e:
                try:
                    raw = e.read().decode("utf-8")
                    body = json.loads(raw)
                except Exception:
                    body = {"raw": str(e)}
                self._log({"event": "http_error", "path": path, "attempt": attempt,
                           "request": payload, "http_status": e.code, "response": body})
                # 400/404/405/409/413/415 基本都是程序/参数错误，不盲目重试
                if e.code in (400, 404, 405, 409, 413, 415):
                    raise RuntimeError(f"HTTP {e.code}: {body}") from e
                last_err = e
            except (URLError, socket.timeout, TimeoutError, ConnectionError) as e:
                self._log({"event": "network_error", "path": path, "attempt": attempt,
                           "request": payload, "error": repr(e)})
                last_err = e
            # 关键：重试使用完全相同 payload/request_id
            if attempt < self.max_retries:
                time.sleep(min(0.5 * (2 ** (attempt - 1)), 3.0))
        raise RuntimeError(f"请求失败，已重试{self.max_retries}次: {path}: {last_err}")

    def enter(self) -> dict:
        rid = self._new_id("enter")
        body = self._post("/enter", self._base(rid))
        self.remaining_real_duration_s = float(body.get("remaining_real_duration_s", 1200.0))
        self.real_start_monotonic = time.monotonic()
        self.current_pos = (0.0, 0.0)
        self.current_channel = 1
        return body

    def measure(self, pos: Point, channel: int) -> dict:
        rid = self._new_id("measure")
        payload = self._base(rid)
        payload["position"] = {"x": float(pos[0]), "y": float(pos[1])}
        payload["channel"] = int(channel)
        old_pos = self.current_pos
        old_channel = self.current_channel
        body = self._post("/measure", payload)
        self.total_move_distance_m += dist(old_pos, pos)
        self.measure_count += 1
        if old_channel != channel:
            self.channel_switch_count += 1
        self.current_pos = pos
        self.current_channel = channel
        return body

    def clear(self, pos: Point, channel: int) -> dict:
        rid = self._new_id("clear")
        payload = self._base(rid)
        payload["position"] = {"x": float(pos[0]), "y": float(pos[1])}
        payload["channel"] = int(channel)
        old_pos = self.current_pos
        body = self._post("/clear", payload)
        self.total_move_distance_m += dist(old_pos, pos)
        self.clear_count += 1
        if body.get("clear_result") == "success":
            self.clear_success_count += 1
        else:
            self.clear_fail_count += 1
        self.current_pos = pos
        # /clear 不改变 current_channel
        return body

    def exit(self) -> Optional[dict]:
        rid = self._new_id("exit")
        try:
            return self._post("/exit", self._base(rid))
        except Exception as e:
            self._log({"event": "exit_error", "error": repr(e)})
            return None


# ----------------------------
# 控制器
# ----------------------------
class Problem3Controller:
    def __init__(self, client: SimulatorClient):
        self.client = client
        self.sources: Dict[int, SourceState] = {}
        self.unknown = set(range(1, 21))
        self.empty = set()
        self.cleared = set()
        self.arena_poly = regular_circumscribed_polygon((0.0, 0.0), ARENA_R, CIRCLE_SIDES)
        self.cover_pts = coverage_points(0.0)
        self.station_scanned: Dict[int, set] = {}
        self.search_scan_points: List[Point] = []
        # 对问题3全向源，发现前的no_signal不是废信息：
        # 若N点无信号、S点有信号，则同一个未知接收半径R满足 d(S,G)<=R<d(N,G)，
        # 因而真实源必位于“更靠近S而不是N”的垂直平分线半平面内。
        self.no_signal_points: Dict[int, List[Point]] = {}

    def _ensure_real_time(self) -> None:
        left = self.client.real_time_left()
        if left <= MAX_REAL_MARGIN_S:
            raise TimeoutError(f"现实运行时间只剩 {left:.1f}s，停止新动作")

    def _apply_signal_vs_nosignal_constraints(self, st: SourceState, signal_pos: Point) -> None:
        """利用历史no_signal免费收紧问题3全向源的可行域。

        同一干扰源有效接收半径R固定。若signal_pos收到信号而旧点n无信号，则
        |G-signal_pos| <= R < |G-n|，所以 |G-signal_pos| < |G-n|。
        两边平方后是严格线性半平面；用<=边界的闭半平面保守包含真实源。
        """
        for n in self.no_signal_points.get(st.channel, []):
            if dist(n, signal_pos) < 1e-7:
                continue
            sx, sy = signal_pos
            nx, ny = n
            const = nx * nx + ny * ny - sx * sx - sy * sy
            dx = nx - sx
            dy = ny - sy
            st.poly = clip_halfplane(
                st.poly,
                lambda x, c=const, ax=dx, ay=dy: c - 2.0 * (ax * x[0] + ay * x[1])
            )
            if not st.poly:
                raise RuntimeError(f"频道{st.channel}应用signal/no_signal约束后可行域为空")

    def _apply_nosignal_vs_signals(self, st: SourceState, no_signal_pos: Point) -> None:
        """已知源在新点无信号时，用所有历史有信号点构造免费半平面约束。"""
        nx, ny = no_signal_pos
        for obs in st.observations:
            sx, sy = obs.pos
            if dist(obs.pos, no_signal_pos) < 1e-7:
                continue
            # d(G,S) < d(G,N)  <=>  |G-N|^2 - |G-S|^2 > 0
            const = nx * nx + ny * ny - sx * sx - sy * sy
            dx = nx - sx
            dy = ny - sy
            st.poly = clip_halfplane(
                st.poly,
                lambda x, c=const, ax=dx, ay=dy: c - 2.0 * (ax * x[0] + ay * x[1])
            )
            if not st.poly:
                raise RuntimeError(f"频道{st.channel}应用no_signal/signal约束后可行域为空")

    def _new_source(self, channel: int, pos: Point, svd_deg: float) -> SourceState:
        st = SourceState(channel=channel, poly=list(self.arena_poly))
        st.add_direction(pos, svd_deg)
        self._apply_signal_vs_nosignal_constraints(st, pos)
        self.sources[channel] = st
        self.unknown.discard(channel)
        print(f"  [发现] 频道{channel}: svd={svd_deg:.2f}°, 初始直径≈{st.diameter:.1f}m")
        return st

    def _handle_measure(self, channel: int, pos: Point, resp: dict) -> None:
        result = resp.get("measure_result")
        if result == "direction":
            svd = float(resp["svd_deg"])
            if channel not in self.sources:
                st = self._new_source(channel, pos, svd)
            else:
                st = self.sources[channel]
                if st.cleared:
                    return
                st.add_direction(pos, svd)
                self._apply_signal_vs_nosignal_constraints(st, pos)
                print(f"  [测向] 频道{channel}: svd={svd:.2f}°, 直径≈{st.diameter:.1f}m, MEC半径≈{st.mec.r:.1f}m")
        elif result == "near":
            print(f"  [近距] 频道{channel} 距离<=5m，直接清除")
            self.unknown.discard(channel)
            if channel not in self.sources:
                self.sources[channel] = SourceState(channel=channel, poly=[], observations=[])
            cr = self.client.clear(pos, channel)
            if cr.get("clear_result") == "success":
                self.sources[channel].cleared = True
                self.cleared.add(channel)
                print(f"  [清除成功] 频道{channel}")
            else:
                print(f"  [警告] near后清除失败: 频道{channel}, {cr}")
        elif result == "no_signal":
            # 问题3全向源：在频道尚未发现时保存无信号点。未来一旦在别处收到信号，
            # 可用同一未知接收半径构造“signal点比no_signal点更近”的免费半平面约束。
            pts = self.no_signal_points.setdefault(channel, [])
            if not any(dist(pos, q) < 1e-6 for q in pts):
                pts.append(pos)
            st = self.sources.get(channel)
            if st is not None and (not st.cleared) and st.poly and st.observations:
                self._apply_nosignal_vs_signals(st, pos)
                print(f"  [无信号约束] 频道{channel}: 直径≈{st.diameter:.1f}m, MEC半径≈{st.mec.r:.1f}m")
        else:
            raise RuntimeError(f"未知 measure_result={result}")

    def scan_channels_at(self, pos: Point, channels: Sequence[int], station_idx: Optional[int] = None) -> None:
        chs = [c for c in channels if c not in self.cleared]
        if not chs:
            return
        # 当前频道优先，减少一次切换
        if self.client.current_channel in chs:
            chs.remove(self.client.current_channel)
            chs.insert(0, self.client.current_channel)
        for c in chs:
            self._ensure_real_time()
            resp = self.client.measure(pos, c)
            if station_idx is not None:
                self.station_scanned.setdefault(station_idx, set()).add(c)
            self._handle_measure(c, pos, resp)

    def _joint_probe_candidates(self, pos: Point, primary: Optional[int] = None,
                                search_mode: bool = False) -> List[Tuple[float, int]]:
        """找出在当前位置“顺手测一次”最划算的其他频道。

        当前位置已经到达，因此额外成本只有5s检测 + 最多1s换台；候选仍要求
        robust_receive，保证问题3全向源一定有信号。
        """
        candidates: List[Tuple[float, int]] = []
        for c, st in self.sources.items():
            if c == primary or st.cleared or not st.poly or st.ready:
                continue
            if any(dist(obs.pos, pos) < 1.0 for obs in st.observations):
                continue
            if not robust_receive(pos, st.poly):
                continue
            cur_d = max(st.diameter, 1e-6)
            pred = predicted_worst_diameter(pos, st)
            ratio = pred / cur_d
            if search_mode:
                limit = 1.0 if len(st.observations) == 1 else STATION_PROBE_IMPROVEMENT
            else:
                limit = JOINT_PROBE_RATIO_FIRST if len(st.observations) == 1 else JOINT_PROBE_RATIO_LATER
            if ratio <= limit + 1e-12:
                candidates.append((ratio, c))
        candidates.sort()
        return candidates

    def _joint_gain_at(self, pos: Point, primary: Optional[int] = None) -> float:
        """估计一个位置可同时服务多个源的价值；用于全局动作排序。"""
        gain = 0.0
        for ratio, _ in self._joint_probe_candidates(pos, primary=primary, search_mode=False):
            gain += max(0.0, 1.0 - ratio)
        return min(gain, 2.5)

    def _probe_known_at_point(self, pos: Point, primary: Optional[int] = None,
                              search_mode: bool = False) -> int:
        """已到达pos后批量测其他源；不产生额外移动。返回实际顺手测量次数。"""
        candidates = self._joint_probe_candidates(pos, primary=primary, search_mode=search_mode)
        if not search_mode:
            candidates = candidates[:MAX_JOINT_PROBES_PER_STOP]
        n = 0
        for ratio, c in candidates:
            if c in self.cleared:
                continue
            st = self.sources.get(c)
            if st is None or st.cleared or st.ready or not st.poly:
                continue
            # 前面某次顺手测量可能改变了几何区域，再做一次便宜的安全检查。
            if any(dist(obs.pos, pos) < 1.0 for obs in st.observations) or not robust_receive(pos, st.poly):
                continue
            print(f"  [同点复用] 频道{c}: 预计直径比例≈{ratio:.2f}")
            self._ensure_real_time()
            resp = self.client.measure(pos, c)
            self._handle_measure(c, pos, resp)
            n += 1
        return n

    def _predicted_nosignal_diameter(self, pos: Point, st: SourceState) -> float:
        """若在pos得到no_signal，利用已有signal点可保证得到的直径上界。"""
        p = list(st.poly)
        nx, ny = pos
        for obs in st.observations:
            sx, sy = obs.pos
            if dist(obs.pos, pos) < 1e-7:
                continue
            const = nx * nx + ny * ny - sx * sx - sy * sy
            dx = nx - sx
            dy = ny - sy
            p = clip_halfplane(
                p,
                lambda x, c=const, ax=dx, ay=dy: c - 2.0 * (ax * x[0] + ay * x[1])
            )
            if not p:
                return 0.0
        return polygon_diameter(p)

    def maybe_probe_known_at_station(self, pos: Point) -> None:
        """搜索/覆盖站移动成本已支付，只补做“6秒内信息价值足够高”的已知源检测。

        对保证接收点评估direction后的最坏直径；对非保证接收点评估no_signal
        能带来的垂直平分线收缩。避免旧版那种大量“无信号但区域完全没变”的空测。
        """
        candidates = []
        for c, st in self.sources.items():
            if st.cleared or not st.poly or st.ready:
                continue
            if any(dist(obs.pos, pos) < 1.0 for obs in st.observations):
                continue
            cur_d = max(st.diameter, 1e-6)
            if robust_receive(pos, st.poly):
                pred = predicted_worst_diameter(pos, st)
                ratio = pred / cur_d
                # 第2条有效方位通常很值钱；之后要求至少有约7%的保证收缩。
                limit = 1.0 if len(st.observations) == 1 else 0.93
                if ratio <= limit + 1e-12:
                    candidates.append((ratio, c, "direction"))
            else:
                ns_d = self._predicted_nosignal_diameter(pos, st)
                ratio = ns_d / cur_d
                # 非保证接收点若返回no_signal也必须能明显缩域，否则这6秒大概率不值。
                if ratio <= 0.88 + 1e-12:
                    candidates.append((ratio + 0.04, c, "range-order"))

        candidates.sort()
        for _, c, mode in candidates:
            st = self.sources.get(c)
            if st is None or st.cleared or not st.poly or st.ready:
                continue
            if any(dist(obs.pos, pos) < 1.0 for obs in st.observations):
                continue
            print(f"  [覆盖站补测] 频道{c} ({mode})")
            self._ensure_real_time()
            resp = self.client.measure(pos, c)
            self._handle_measure(c, pos, resp)

    def try_clear_source(self, channel: int) -> bool:
        st = self.sources[channel]
        if st.cleared:
            return True
        if not st.poly:
            return False
        mec = st.mec
        if mec.r > READY_RADIUS:
            return False
        self._ensure_real_time()
        st.clear_attempts += 1
        print(f"  [清除尝试] 频道{channel} -> ({mec.c[0]:.2f},{mec.c[1]:.2f}), 保守MEC半径={mec.r:.2f}m")
        resp = self.client.clear(mec.c, channel)
        if resp.get("clear_result") == "success":
            st.cleared = True
            self.cleared.add(channel)
            print(f"  [清除成功] 频道{channel}")
            return True
        print(f"  [清除失败] 频道{channel}，继续增加测向；响应={resp.get('clear_result')}")
        return False

    def _clear_ready_along_route(self, destination: Point) -> int:
        """去下一搜索站途中顺路清除已就绪源。

        只接受额外路程 <= SEARCH_CLEAR_MAX_DETOUR_M 的清除点，因此不会为了清源
        大幅破坏7站搜索路线；而被顺路清掉的源以后无需专程再访问。
        """
        done = 0
        while True:
            cur = self.client.current_pos
            direct = dist(cur, destination)
            choices = []
            for c, st in self.sources.items():
                if st.cleared or not st.poly or not st.ready:
                    continue
                p = st.mec.c
                extra = dist(cur, p) + dist(p, destination) - direct
                if extra <= SEARCH_CLEAR_MAX_DETOUR_M + 1e-7:
                    choices.append((extra, dist(cur, p), c, p))
            if not choices:
                break
            choices.sort()
            _, _, c, p = choices[0]
            print(f"  [搜索途中顺路清除] 频道{c}，相对直达额外路程≈{choices[0][0]:.1f}m")
            if not self.try_clear_source(c):
                break
            done += 1
            # 清除点也当作免费定位站复用。
            self._probe_known_at_point(p, primary=c, search_mode=False)
        return done

    def _adaptive_source_candidates(self, remaining_skeleton: Sequence[Tuple[int, Point]]) -> List[Tuple[Point, Optional[int], str]]:
        """生成“搜索+定位复用”的候选站。

        除原保证骨架点外，对当前最难定位的若干源求其主动定位理想点，并把附近骨架点
        朝理想点偏移5/10/15/20%。偏移点只有在后续仍能用剩余站数完成严格覆盖时才会被采用。
        """
        candidates: List[Tuple[Point, Optional[int], str]] = []
        for idx, p in remaining_skeleton:
            candidates.append((p, None, f"骨架{idx}"))

        # 若某源已经定位到可清除，把它的MEC中心也作为“清除+搜索”二合一候选站。
        # 只要该点通过后续覆盖完成性证书，它就可以替代一个专门的搜索骨架站。
        ready_states = [st for st in self.sources.values()
                        if (not st.cleared and st.poly and st.ready)]
        ready_states.sort(key=lambda st: dist(self.client.current_pos, st.mec.c))
        for st in ready_states[:ADAPTIVE_READY_CLEAR_LIMIT]:
            q = st.mec.c
            if inside_arena(q):
                candidates.append((q, st.channel, f"频道{st.channel}清除+搜索"))

        active = [st for st in self.sources.values()
                  if (not st.cleared and st.poly and not st.ready and st.observations)]
        active.sort(key=lambda st: st.diameter, reverse=True)
        for st in active[:ADAPTIVE_SOURCE_LIMIT]:
            try:
                ideal = choose_next_measure_point(st, self.client.current_pos)
            except Exception:
                continue
            if inside_arena(ideal) and robust_receive(ideal, st.poly):
                candidates.append((ideal, st.channel, f"频道{st.channel}理想定位点"))

            # 只向距离理想点最近的两个剩余骨架站做偏移，控制候选数量。
            near_skeleton = sorted(remaining_skeleton, key=lambda z: dist(z[1], ideal))[:2]
            for idx, base in near_skeleton:
                for t in ADAPTIVE_SHIFT_FACTORS:
                    q = blend_point(base, ideal, t)
                    if not inside_arena(q):
                        continue
                    if not robust_receive(q, st.poly):
                        continue
                    if any(dist(q, obs.pos) < 1.0 for obs in st.observations):
                        continue
                    candidates.append((q, st.channel, f"骨架{idx}→频道{st.channel}偏移{int(t*100)}%"))

        # 坐标去重；若同一点既是骨架又是某源候选，优先保留带主频道的版本。
        unique: Dict[Tuple[int, int], Tuple[Point, Optional[int], str]] = {}
        for q, ch, label in candidates:
            key = (round(q[0] * 1000), round(q[1] * 1000))
            old = unique.get(key)
            if old is None or (old[1] is None and ch is not None):
                unique[key] = (q, ch, label)
        return list(unique.values())

    def _count_joint_localization_opportunities(self, q: Point, primary: Optional[int]) -> int:
        n = 0
        for c, st in self.sources.items():
            if c == primary or st.cleared or not st.poly or st.ready:
                continue
            if any(dist(q, obs.pos) < MIN_ACTIVE_BASELINE for obs in st.observations):
                continue
            if robust_receive(q, st.poly):
                n += 1
        return n

    def _choose_adaptive_search_station(self, remaining_skeleton: Sequence[Tuple[int, Point]]) -> Tuple[Point, Optional[int], str, List[Tuple[int, Point]], float]:
        """在不破坏“最多7站完成严格覆盖”的前提下选下一搜索站。"""
        cur = self.client.current_pos
        visited = list(self.search_scan_points)
        current_radius, _ = guaranteed_cover_radius(visited)
        budget_after = MAX_SEARCH_STATIONS - (len(visited) + 1)
        if budget_after < 0:
            raise RuntimeError("自适应搜索站数超过MAX_SEARCH_STATIONS")

        scored = []
        for q, primary, label in self._adaptive_source_candidates(remaining_skeleton):
            if any(dist(q, old) < 1e-6 for old in visited):
                continue

            # 若q本身就是某个未访问骨架点，则从后续骨架候选里移除，避免重复计数。
            rem_after = [(idx, p) for idx, p in remaining_skeleton if dist(q, p) > 1e-6]
            completion = minimum_completion_subset(visited + [q], rem_after, budget_after)
            if completion is None:
                continue  # 关键安全闸：这个自适应点会让7站预算无法保证全覆盖

            cover_radius, _ = guaranteed_cover_radius(visited + [q])
            move_cost = dist(cur, q)
            loc_gain = 0.0
            clear_gain = 0.0
            if primary is not None and primary in self.sources:
                st = self.sources[primary]
                if st.ready and dist(q, st.mec.c) <= 1e-4:
                    clear_gain = 1.0
                else:
                    cur_d = max(st.diameter, 1e-6)
                    try:
                        pred = predicted_worst_diameter(q, st)
                        loc_gain = max(0.0, min(1.0, 1.0 - pred / cur_d))
                    except Exception:
                        loc_gain = 0.0
            joint = self._count_joint_localization_opportunities(q, primary)

            # completion越少越优先；之后在同等可完成性下兼顾移动、覆盖收缩、定位与顺路清除收益。
            score = (
                len(completion) * ADAPTIVE_REMAINING_STATION_PENALTY
                + move_cost
                + cover_radius * ADAPTIVE_COVER_RADIUS_WEIGHT
                - loc_gain * ADAPTIVE_LOCALIZE_BONUS
                - clear_gain * ADAPTIVE_CLEAR_BONUS
                - joint * ADAPTIVE_JOINT_BONUS
            )
            scored.append((score, move_cost, cover_radius, -loc_gain, q, primary, label, completion, joint))

        if not scored:
            # 理论上只要原7点骨架未被破坏就不应发生；仍给出明确失败而不是偷偷失去覆盖保证。
            raise RuntimeError("找不到能在7站预算内保持全覆盖保证的下一搜索点")

        scored.sort(key=lambda x: (x[0], x[1]))
        _, _, cover_radius, neg_gain, q, primary, label, completion, joint = scored[0]
        print(f"  [自适应决策] {label}: 当前覆盖半径≈{current_radius:.1f}m -> 访问后≈{cover_radius:.1f}m, "
              f"后续最少骨架站={len(completion)}, 定位收益≈{-neg_gain:.1%}, 可顺带源={joint}")
        return q, primary, label, completion, cover_radius

    def _execute_search_station(self, q: Point, station_idx: int, primary: Optional[int] = None) -> None:
        """在q扫描全部尚未知频道，并在不增加移动的情况下利用该站定位已知源。"""
        print(f"\n[自适应搜索站 {station_idx}] ({q[0]:.1f}, {q[1]:.1f})，未知频道数={len(self.unknown)}")
        self.scan_channels_at(q, sorted(self.unknown), station_idx=station_idx)
        self.search_scan_points.append(q)

        # 若这个搜索站本身就是某个ready源的MEC中心，扫描完未知频道后原地清除。
        if primary is not None and primary in self.sources and primary not in self.cleared:
            pst = self.sources[primary]
            if pst.poly and pst.ready and dist(q, pst.mec.c) <= 1e-3:
                self.try_clear_source(primary)

        # 若该站是为某个未ready主频道偏移出来的，兑现这次定位收益。
        if primary is not None and primary in self.sources and primary not in self.cleared:
            st = self.sources[primary]
            if (st.poly and not st.ready and
                    all(dist(q, obs.pos) >= 1.0 for obs in st.observations) and
                    robust_receive(q, st.poly)):
                print(f"  [自适应主测] 频道{primary}")
                resp = self.client.measure(q, primary)
                self._handle_measure(primary, q, resp)

        self.maybe_probe_known_at_station(q)

    def choose_coverage_orientation(self) -> float:
        """选择后备六边形整体旋转角。

        第二轮优化优先最小化“当前已知源区域锚点 + 6个覆盖点”的开放路线长度；
        几何交会质量只作轻量次级项。旋转不改变1000m覆盖保证。
        """
        active = [st for st in self.sources.values() if (not st.cleared and st.poly and st.observations)]
        best = (float("inf"), 0.0)
        cur = self.client.current_pos
        for phi in range(0, 60, 5):
            outer = coverage_points(float(phi))[1:]
            nodes: List[Tuple[str, int, Point]] = []
            for i, q in enumerate(outer, start=1):
                nodes.append(("cover", i, q))
            for st in active:
                nodes.append(("source", st.channel, st.mec.c))

            route = self._optimize_open_route(cur, nodes) if hasattr(self, "_optimize_open_route") else nodes
            if route:
                route_m = dist(cur, route[0][2]) + sum(
                    dist(route[i - 1][2], route[i][2]) for i in range(1, len(route))
                )
            else:
                route_m = 0.0

            # 次级几何项只用于打破路线接近时的平局，避免完全忽略顺路形成好交会的机会。
            geom = 0.0
            for st in active:
                vals = [geometric_score(q, st) for q in outer
                        if robust_receive(q, st.poly) and
                        all(dist(q, obs.pos) >= MIN_ACTIVE_BASELINE for obs in st.observations)]
                if vals:
                    geom += min(vals) / 1e6
                else:
                    geom += 1e6
            score = route_m + 0.02 * geom
            if score < best[0]:
                best = (score, float(phi))
        return best[1]

    def global_search(self) -> None:
        print("\n=== 阶段A：自适应7站保证搜索 ===")

        # 第0站固定在原点。之后仍以一个完整六边形骨架作为“安全后备”，
        # 但每一步允许把骨架站向高价值定位点偏移；偏移必须通过覆盖完成性证书。
        p0 = (0.0, 0.0)
        self._execute_search_station(p0, station_idx=0, primary=None)

        if len(self.sources) >= 16:
            print("已发现16个频道，达到题目源数量上界，提前结束全域搜索。")
            self.empty |= self.unknown
            self.unknown.clear()
            return

        phi = self.choose_coverage_orientation()
        self.cover_pts = coverage_points(phi)
        remaining_skeleton: List[Tuple[int, Point]] = [(i, self.cover_pts[i]) for i in range(1, 7)]
        print(f"保证骨架采用旋转角 phi={phi:.1f}°；后续站允许自适应偏移，但始终保留7站覆盖证书。")

        station_no = 1
        while self.unknown:
            if len(self.sources) >= 16:
                print("已发现16个频道，达到题目源数量上界，提前结束全域搜索。")
                self.empty |= self.unknown
                self.unknown.clear()
                break

            radius, worst_p = guaranteed_cover_radius(self.search_scan_points)
            if radius <= MIN_RECEIVE_R + ADAPTIVE_COVER_TOL:
                print(f"\n覆盖证书完成：最大最近搜索站距离≈{radius:.3f}m <= {MIN_RECEIVE_R:.1f}m。")
                print(f"剩余未知频道可严格判空: {sorted(self.unknown)}")
                self.empty |= self.unknown
                self.unknown.clear()
                break

            if len(self.search_scan_points) >= MAX_SEARCH_STATIONS:
                raise RuntimeError(
                    f"已使用{MAX_SEARCH_STATIONS}个搜索站但覆盖证书仍为{radius:.3f}m，最坏点={worst_p}"
                )

            q, primary, label, completion, _ = self._choose_adaptive_search_station(remaining_skeleton)
            # 去下一保证搜索站之前，若有已定位源几乎在线路上，先顺路清除。
            self._clear_ready_along_route(q)
            self._execute_search_station(q, station_idx=station_no, primary=primary)
            station_no += 1

            # 若本次访问与某个原骨架点重合，则把它从后备集合移除。
            remaining_skeleton = [(idx, p) for idx, p in remaining_skeleton if dist(q, p) > 1e-6]

        # 最终再做一次证书输出，便于演练日志核查。
        if self.search_scan_points:
            radius, worst_p = guaranteed_cover_radius(self.search_scan_points)
            print(f"搜索阶段结束：共使用{len(self.search_scan_points)}个扫描站，覆盖半径证书≈{radius:.3f}m，最坏点≈({worst_p[0]:.1f},{worst_p[1]:.1f})")

    def _optimize_open_route(self, start: Point, nodes: Sequence[Tuple[str, int, Point]]) -> List[Tuple[str, int, Point]]:
        """开放式最近邻 + 2-opt。路线无需回到起点，适合机器狗剩余任务排序。"""
        rem = list(nodes)
        route: List[Tuple[str, int, Point]] = []
        cur = start
        while rem:
            k = min(range(len(rem)), key=lambda i: dist(cur, rem[i][2]))
            node = rem.pop(k)
            route.append(node)
            cur = node[2]
        if len(route) < 3:
            return route

        def path_len(rt: Sequence[Tuple[str, int, Point]]) -> float:
            if not rt:
                return 0.0
            total = dist(start, rt[0][2])
            for i in range(1, len(rt)):
                total += dist(rt[i - 1][2], rt[i][2])
            return total

        best_len = path_len(route)
        for _ in range(REGION_ROUTE_2OPT_PASSES):
            improved = False
            n = len(route)
            for i in range(n - 1):
                prev = start if i == 0 else route[i - 1][2]
                for j in range(i + 1, n):
                    nxt = route[j + 1][2] if j + 1 < n else None
                    old = dist(prev, route[i][2])
                    new = dist(prev, route[j][2])
                    if nxt is not None:
                        old += dist(route[j][2], nxt)
                        new += dist(route[i][2], nxt)
                    if new + 1e-7 < old:
                        cand = route[:i] + list(reversed(route[i:j + 1])) + route[j + 1:]
                        cand_len = path_len(cand)
                        if cand_len + 1e-7 < best_len:
                            route = cand
                            best_len = cand_len
                            improved = True
            if not improved:
                break
        return route

    def _region_anchor_nodes(self, completion: Sequence[Tuple[int, Point]]) -> List[Tuple[str, int, Point]]:
        """把每个源压缩为稳定MEC锚点，与必须补的覆盖点统一规划。"""
        nodes: List[Tuple[str, int, Point]] = []
        for c, st in self.sources.items():
            if st.cleared or not st.poly:
                continue
            nodes.append(("source", c, st.mec.c))
        for idx, p in completion:
            nodes.append(("cover", idx, p))
        return nodes

    def _plan_region_route(self, completion: Sequence[Tuple[int, Point]]) -> List[Tuple[str, int, Point]]:
        return self._optimize_open_route(self.client.current_pos, self._region_anchor_nodes(completion))

    def _route_aware_completion(self, remaining_skeleton: Sequence[Tuple[int, Point]]) -> Optional[List[Tuple[int, Point]]]:
        """在“最少补点数”相同的方案中，选与剩余源区域联合路线最短的一组覆盖点。"""
        rem = list(remaining_skeleton)
        for k in range(len(rem) + 1):
            feasible: List[Tuple[float, List[Tuple[int, Point]]]] = []
            for idxs in itertools.combinations(range(len(rem)), k):
                subset = [rem[i] for i in idxs]
                radius, _ = guaranteed_cover_radius(self.search_scan_points + [p for _, p in subset])
                if radius > MIN_RECEIVE_R + ADAPTIVE_COVER_TOL:
                    continue
                route = self._plan_region_route(subset)
                if route:
                    length = dist(self.client.current_pos, route[0][2]) + sum(
                        dist(route[i - 1][2], route[i][2]) for i in range(1, len(route))
                    )
                else:
                    length = 0.0
                feasible.append((length, subset))
            if feasible:
                feasible.sort(key=lambda x: x[0])
                return feasible[0][1]
        return None

    def _dynamic_route_aware_completion(self) -> Tuple[Optional[List[Tuple[int, Point]]], float]:
        """动态旋转后备六边形。

        已有定位/清除复用扫描点发生变化后，重新枚举phi=0..55°：
        先最小化仍需补的覆盖点数量，再在同数量方案中最小化与所有源区域锚点的开放路线。
        最终安全性仍由 guaranteed_cover_radius 的精确证书决定。
        """
        best_k: Optional[int] = None
        best_len = float("inf")
        best_subset: Optional[List[Tuple[int, Point]]] = None
        best_phi = 0.0
        visited = list(self.search_scan_points)

        for phi in range(0, 60, 5):
            outer = [(i, p) for i, p in enumerate(coverage_points(float(phi))[1:], start=1)
                     if not any(dist(p, v) < 1e-6 for v in visited)]
            max_k = len(outer) if best_k is None else min(len(outer), best_k)
            found_this_phi = False
            for k in range(max_k + 1):
                if best_k is not None and k > best_k:
                    break
                for idxs in itertools.combinations(range(len(outer)), k):
                    subset = [outer[i] for i in idxs]
                    radius, _ = guaranteed_cover_radius(visited + [p for _, p in subset])
                    if radius > MIN_RECEIVE_R + ADAPTIVE_COVER_TOL:
                        continue
                    found_this_phi = True
                    route = self._plan_region_route(subset)
                    if route:
                        length = dist(self.client.current_pos, route[0][2]) + sum(
                            dist(route[i - 1][2], route[i][2]) for i in range(1, len(route))
                        )
                    else:
                        length = 0.0
                    if best_k is None or k < best_k or (k == best_k and length < best_len):
                        best_k = k
                        best_len = length
                        best_subset = subset
                        best_phi = float(phi)
                # 对同一phi，一旦k可完成覆盖，就无需枚举更大的k。
                if found_this_phi:
                    break
            if best_k == 0:
                break
        return best_subset, best_phi

    def _local_clear_grid(self, st: SourceState) -> Optional[List[Point]]:
        """若当前保守可行域很小，构造至多6个保证覆盖它的清除点。

        对多个旋转角建立外接矩形网格；每个单元半对角线<=18.5m，因此单元内任一点
        到单元中心都小于20m。覆盖整个外接矩形自然覆盖真实可行域，保证最终有一格能清除。
        """
        if not st.poly:
            return None
        origin = centroid(st.poly)
        side = LOCAL_CLEAR_CELL_RADIUS * math.sqrt(2.0)
        best: Optional[List[Point]] = None
        ox, oy = origin
        for deg in LOCAL_CLEAR_ANGLES_DEG:
            a = math.radians(deg)
            ca, sa = math.cos(a), math.sin(a)

            def to_local(p: Point) -> Point:
                x, y = p[0] - ox, p[1] - oy
                return (ca * x + sa * y, -sa * x + ca * y)

            def to_world(p: Point) -> Point:
                x, y = p
                return (ox + ca * x - sa * y, oy + sa * x + ca * y)

            rp = [to_local(v) for v in st.poly]
            xmin = min(x for x, _ in rp)
            xmax = max(x for x, _ in rp)
            ymin = min(y for _, y in rp)
            ymax = max(y for _, y in rp)
            w = max(0.0, xmax - xmin)
            h = max(0.0, ymax - ymin)
            nx = max(1, int(math.ceil(w / side)))
            ny = max(1, int(math.ceil(h / side)))
            count = nx * ny
            if count > LOCAL_CLEAR_MAX_POINTS:
                continue
            cw = w / nx if w > 1e-12 else 0.0
            ch = h / ny if h > 1e-12 else 0.0
            if 0.5 * math.hypot(cw, ch) > LOCAL_CLEAR_CELL_RADIUS + 1e-7:
                continue
            pts: List[Point] = []
            for ix in range(nx):
                x = xmin + (ix + 0.5) * cw if w > 1e-12 else xmin
                for iy in range(ny):
                    y = ymin + (iy + 0.5) * ch if h > 1e-12 else ymin
                    q = to_world((x, y))
                    if not inside_arena(q):
                        pts = []
                        break
                    pts.append(q)
                if not pts and count > 0:
                    break
            if len(pts) != count:
                continue
            ordered_nodes = self._optimize_open_route(
                self.client.current_pos, [("clearcell", i, q) for i, q in enumerate(pts)]
            )
            ordered = [q for _, _, q in ordered_nodes]
            if best is None or len(ordered) < len(best):
                best = ordered
        return best

    def _local_clear_sweep_cost_s(self, points: Sequence[Point]) -> float:
        if not points:
            return float("inf")
        cur = self.client.current_pos
        move = 0.0
        for q in points:
            move += dist(cur, q)
            cur = q
        # 最坏情况：前n-1次失败各3s，最后一次成功5s。
        return move / MOVE_SPEED_MPS + max(0, len(points) - 1) * 3.0 + CLEAR_ACTION_S

    def _task_action_for_channel(self, c: int) -> Optional[Tuple[str, int, Point, float]]:
        st = self.sources.get(c)
        if st is None or st.cleared or not st.poly:
            return None
        cur = self.client.current_pos
        if st.ready:
            p = st.mec.c
            return ("clear", c, p, dist(cur, p) / MOVE_SPEED_MPS + CLEAR_ACTION_S)

        q = choose_next_measure_point(st, cur)
        pred = predicted_worst_diameter(q, st)
        measure_s = _time_score_measure_candidate(q, pred, st, cur)
        if self.client.current_channel != c:
            measure_s += SWITCH_CHANNEL_S

        sweep = self._local_clear_grid(st)
        if sweep:
            sweep_s = self._local_clear_sweep_cost_s(sweep)
            if sweep_s <= measure_s + LOCAL_CLEAR_FUTURE_MARGIN_S:
                return ("sweep", c, sweep[0], sweep_s)
        return ("measure", c, q, measure_s)

    def _execute_local_clear_sweep(self, c: int) -> bool:
        st = self.sources[c]
        points = self._local_clear_grid(st)
        if not points:
            return False
        print(f"  [保证式局部扫清] 频道{c}: {len(points)}个清除格")
        for i, q in enumerate(points, start=1):
            self._ensure_real_time()
            st.clear_attempts += 1
            print(f"    清除格 {i}/{len(points)} -> ({q[0]:.2f},{q[1]:.2f})")
            resp = self.client.clear(q, c)
            if resp.get("clear_result") == "success":
                st.cleared = True
                self.cleared.add(c)
                print(f"  [局部扫清成功] 频道{c}")
                return True
        # 理论上可行域包含真实源且网格覆盖整个外接矩形时不应发生。
        raise RuntimeError(f"频道{c}保证式局部扫清全部失败；可行域/清除半径假设可能不一致")

    def next_task_action(self) -> Optional[Tuple[str, int, Point, float]]:
        """全局选择下一动作：ready源可清除，未ready源可继续测向。

        与旧版“只选最近频道”不同，这里直接按预计虚拟秒数排序，并加入：
        1) 测后靠近MEC中心的收益；2) 当前位置可顺手服务其他源的收益；
        3) 很轻的下一步路程前瞻，减少最近邻策略造成的折返。
        """
        cur = self.client.current_pos
        raw = []
        for c, st in self.sources.items():
            if st.cleared:
                continue
            if not st.poly:
                continue
            if st.ready:
                p = st.mec.c
                base_s = dist(cur, p) / MOVE_SPEED_MPS + CLEAR_ACTION_S
                kind = "clear"
            else:
                p = choose_next_measure_point(st, cur)
                pred = predicted_worst_diameter(p, st)
                base_s = _time_score_measure_candidate(p, pred, st, cur)
                if self.client.current_channel != c:
                    base_s += SWITCH_CHANNEL_S
                kind = "measure"
            raw.append({"kind": kind, "channel": c, "point": p, "base": base_s})

        if not raw:
            return None

        # 先用便宜的base分筛前几名，避免对所有n^2组合反复做角域预测。
        raw.sort(key=lambda x: x["base"])
        shortlist = raw[:min(7, len(raw))]
        all_points = [(x["channel"], x["point"]) for x in raw]

        best = None
        for item in shortlist:
            p = item["point"]
            c = item["channel"]
            joint_gain = self._joint_gain_at(p, primary=c)

            # 一步前瞻：两个base接近时，优先选择完成后离下一任务也近的点。
            future = [dist(p, q) / MOVE_SPEED_MPS for oc, q in all_points if oc != c]
            lookahead_s = min(future) if future else 0.0
            score = (item["base"]
                     + TASK_LOOKAHEAD_WEIGHT * lookahead_s
                     - JOINT_GAIN_VALUE_S * joint_gain)
            cand = (score, item["base"], item["kind"], c, p, joint_gain)
            if best is None or cand[:2] < best[:2]:
                best = cand

        assert best is not None
        score, base_s, kind, c, p, joint_gain = best
        print(f"  [时间调度] {kind} 频道{c}: base≈{base_s:.1f}s, "
              f"同点复用收益指标≈{joint_gain:.2f}, 综合分≈{score:.1f}s")
        return kind, c, p, score

    def localize_and_clear_all(self) -> None:
        print("\n=== 阶段B：时间最优的联合定位与清除 ===")
        guard = 0
        while True:
            guard += 1
            if guard > 260:
                raise RuntimeError("任务循环超过260次，防止死循环")

            action = self.next_task_action()
            if action is None:
                break
            kind, c, q, _ = action
            st = self.sources[c]

            if kind == "clear":
                if self.try_clear_source(c):
                    # 清除点本身也已经付了移动成本，把它当共享定位站使用。
                    self._probe_known_at_point(self.client.current_pos, primary=c, search_mode=False)
                continue

            if len(st.observations) >= MAX_LOCALIZE_MEASURES:
                print(f"  [告警] 频道{c}已测向{len(st.observations)}次仍未达到{READY_RADIUS}m，继续尝试")

            print(f"\n[时间型主动测向] 频道{c} -> ({q[0]:.2f},{q[1]:.2f}), "
                  f"当前直径={st.diameter:.1f}m, MEC={st.mec.r:.1f}m")
            self._ensure_real_time()
            resp = self.client.measure(q, c)
            self._handle_measure(c, q, resp)

            # 关键优化：到达一次测点后，最多顺手服务若干其他源，不再“一次移动只服务一个源”。
            self._probe_known_at_point(q, primary=c, search_mode=False)

            # 不再无条件立刻跑向清除点；只有清除点较近时才立即清，远时交回全局调度，
            # 这样能避免“测完A->清A->再折返去B”的锯齿路径。
            st2 = self.sources.get(c)
            if st2 is not None and (not st2.cleared) and st2.poly and st2.ready:
                mec = st2.mec
                if dist(self.client.current_pos, mec.c) <= IMMEDIATE_CLEAR_DISTANCE_M:
                    if self.try_clear_source(c):
                        self._probe_known_at_point(self.client.current_pos, primary=c, search_mode=False)

    def _scan_unknown_aux(self, pos: Point, force: bool = False, label: str = "复用点") -> bool:
        """在已经到达的定位/清除点顺便扫描未知频道，并把该点加入覆盖证书。

        不为扫描额外移动；若离既有扫描点太近则跳过，避免用大量5~6s检测换来几乎重复的覆盖。
        """
        if not self.unknown:
            return False
        if len(self.sources) >= 16:
            self.empty |= self.unknown
            self.unknown.clear()
            return False
        if not force and self.search_scan_points:
            nearest = min(dist(pos, q) for q in self.search_scan_points)
            if nearest < AUX_SCAN_MIN_SPACING_M:
                return False
        print(f"  [搜索复用] {label} ({pos[0]:.1f},{pos[1]:.1f}) 扫描未知频道{len(self.unknown)}个")
        self.scan_channels_at(pos, sorted(self.unknown), station_idx=None)
        if not any(dist(pos, q) < 1e-6 for q in self.search_scan_points):
            self.search_scan_points.append(pos)
        return True

    def _execute_interleaved_action(self, action: Tuple[str, int, Point, float]) -> None:
        """执行一次源任务，并把到达位置复用为未知频道搜索点。"""
        kind, c, q, _ = action
        st = self.sources[c]

        if kind == "clear":
            if self.try_clear_source(c):
                p = self.client.current_pos
                self._scan_unknown_aux(p, force=False, label=f"频道{c}清除点")
                self._probe_known_at_point(p, primary=c, search_mode=False)
            return

        if kind == "sweep":
            if self._execute_local_clear_sweep(c):
                p = self.client.current_pos
                self._scan_unknown_aux(p, force=False, label=f"频道{c}局部扫清点")
                self._probe_known_at_point(p, primary=c, search_mode=False)
            return

        if len(st.observations) >= MAX_LOCALIZE_MEASURES:
            print(f"  [告警] 频道{c}已测向{len(st.observations)}次仍未达到{READY_RADIUS}m，继续尝试")
        print(f"\n[交织主动测向] 频道{c} -> ({q[0]:.2f},{q[1]:.2f}), "
              f"当前直径={st.diameter:.1f}m, MEC={st.mec.r:.1f}m")
        self._ensure_real_time()
        resp = self.client.measure(q, c)
        self._handle_measure(c, q, resp)

        # 先利用这个新位置扩大未知频道覆盖，再做其他已知源的同点复用。
        self._scan_unknown_aux(q, force=False, label=f"频道{c}定位点")
        self._probe_known_at_point(q, primary=c, search_mode=False)

        st2 = self.sources.get(c)
        if st2 is not None and (not st2.cleared) and st2.poly and st2.ready:
            mec = st2.mec
            if dist(self.client.current_pos, mec.c) <= IMMEDIATE_CLEAR_DISTANCE_M:
                if self.try_clear_source(c):
                    p = self.client.current_pos
                    self._scan_unknown_aux(p, force=False, label=f"频道{c}即时清除点")
                    self._probe_known_at_point(p, primary=c, search_mode=False)

    def interleaved_search_localize_clear(self) -> None:
        """问题3主流程：让搜索站与定位/清除路径互相复用。

        正确性仍由两道闸保证：
        1) 已发现源只在保守MEC<=READY_RADIUS时做保证清除；
        2) 未发现频道只有在所有扫描点形成1000m全域覆盖证书，或已发现16源时，才判空。
        """
        print("\n=== 交织阶段：搜索 + 定位 + 清除 ===")
        p0 = (0.0, 0.0)
        self._scan_unknown_aux(p0, force=True, label="原点")

        phi = self.choose_coverage_orientation()
        self.cover_pts = coverage_points(phi)
        print(f"初始后备骨架参考 phi={phi:.1f}°；后续每轮会动态旋转，只补真正必要的覆盖洞。")

        guard = 0
        while True:
            guard += 1
            if guard > 320:
                raise RuntimeError("交织任务循环超过320次，防止死循环")

            # 已发现达到上界，无需再为不存在的频道做覆盖搜索。
            if len(self.sources) >= 16 and self.unknown:
                self.empty |= self.unknown
                self.unknown.clear()

            # 若仍有未知频道，先看当前复用扫描点是否已经严格覆盖整个圆域。
            if self.unknown and self.search_scan_points:
                radius, worst_p = guaranteed_cover_radius(self.search_scan_points)
                if radius <= MIN_RECEIVE_R + ADAPTIVE_COVER_TOL:
                    print(f"  [覆盖完成] 半径证书≈{radius:.3f}m，剩余频道严格判空。")
                    self.empty |= self.unknown
                    self.unknown.clear()

            pending = [st for st in self.sources.values() if not st.cleared]
            if not self.unknown and not pending:
                break

            # 若未知频道尚未被严格覆盖，求“最少补点数 + 与剩余源联合路线最短”的后备点集合。
            completion: List[Tuple[int, Point]] = []
            completion_phi = 0.0
            if self.unknown:
                comp, completion_phi = self._dynamic_route_aware_completion()
                if comp is None:
                    raise RuntimeError("动态后备覆盖无法完成1000m覆盖证书")
                completion = comp
                if not completion:
                    radius, _ = guaranteed_cover_radius(self.search_scan_points)
                    if radius <= MIN_RECEIVE_R + ADAPTIVE_COVER_TOL:
                        self.empty |= self.unknown
                        self.unknown.clear()
                        continue

            # 区域级规划：源MEC中心 + 必须补的覆盖点统一做开放式TSP/2-opt。
            # 只执行整条短路线的第一个区域，动作后根据新观测重新规划。
            route = self._plan_region_route(completion)
            if route:
                node_kind, node_id, anchor = route[0]
                if node_kind == "cover":
                    idx, q = node_id, anchor
                    print(f"\n[区域TSP覆盖点 phi={completion_phi:.0f}° / {idx}] ({q[0]:.1f},{q[1]:.1f})")
                    self._scan_unknown_aux(q, force=True, label=f"动态后备phi={completion_phi:.0f}°点{idx}")
                    # 搜索移动已经支付：direction/no_signal都可能继续严格缩域。
                    self.maybe_probe_known_at_station(q)
                    continue

                if node_kind == "source":
                    action = self._task_action_for_channel(node_id)
                    if action is not None:
                        print(f"  [区域TSP] 下一源区域=频道{node_id}, 锚点=({anchor[0]:.1f},{anchor[1]:.1f})")
                        self._execute_interleaved_action(action)
                        continue

            # 兜底：若区域路线暂时为空，但仍有源任务，使用原时间调度。
            if pending:
                action = self.next_task_action()
                if action is not None:
                    self._execute_interleaved_action(action)
                    continue
            break

        if self.search_scan_points:
            radius, worst_p = guaranteed_cover_radius(self.search_scan_points)
            print(f"交织流程结束：扫描点={len(self.search_scan_points)}，覆盖证书≈{radius:.3f}m，"
                  f"清除={len(self.cleared)}/{len(self.sources)}")

    def summary(self) -> dict:
        found_count = len(self.sources)
        cleared_count = len(self.cleared)
        real_elapsed = None
        if self.client.real_start_monotonic is not None:
            real_elapsed = time.monotonic() - self.client.real_start_monotonic
        return {
            "found_channels": sorted(self.sources.keys()),
            "found_count": found_count,
            "cleared_channels": sorted(self.cleared),
            "cleared_count": cleared_count,
            "clear_rate_vs_found": (cleared_count / found_count) if found_count else 0.0,
            "average_virtual_localize_clear_time_s": (self.client.last_virtual_time / cleared_count) if cleared_count else None,
            "empty_channels": sorted(self.empty),
            "virtual_time_s": self.client.last_virtual_time,
            "local_real_elapsed_s": real_elapsed,
            "route_metrics": {
                "move_distance_m": self.client.total_move_distance_m,
                "estimated_move_time_s": self.client.total_move_distance_m / MOVE_SPEED_MPS,
                "measure_count": self.client.measure_count,
                "channel_switch_count": self.client.channel_switch_count,
                "clear_count": self.client.clear_count,
                "clear_success_count": self.client.clear_success_count,
                "clear_fail_count": self.client.clear_fail_count,
            },
            "source_details": {
                str(c): {
                    "cleared": st.cleared,
                    "observations": len(st.observations),
                    "diameter_m": st.diameter if st.poly else None,
                    "mec_radius_m": st.mec.r if st.poly else None,
                    "clear_attempts": st.clear_attempts,
                }
                for c, st in self.sources.items()
            },
        }

    def run(self) -> dict:
        self.interleaved_search_localize_clear()
        return self.summary()


# ----------------------------
# 自检
# ----------------------------
def self_test() -> None:
    print("运行本地自检……")
    pts = coverage_points(0.0)
    worst = verify_coverage(pts, samples_r=60, samples_a=720)
    exact_worst, exact_p = guaranteed_cover_radius(pts)
    print(f"7点覆盖离散验证：目标圆内最大最近站距离 ≈ {worst:.3f} m")
    print(f"7点覆盖证书：最大最近站距离 ≈ {exact_worst:.3f} m, 最坏点≈({exact_p[0]:.1f},{exact_p[1]:.1f})")
    if exact_worst > MIN_RECEIVE_R + ADAPTIVE_COVER_TOL:
        raise AssertionError("7点覆盖失败")

    # 自适应安全闸测试：中心之后，完整6个外圈骨架应能在剩余6站预算内完成覆盖；
    # 但只允许再加5站时，不应错误声称已经有保证完成方案。
    rem = list(enumerate(pts[1:], start=1))
    comp6 = minimum_completion_subset([pts[0]], rem, 6)
    comp5 = minimum_completion_subset([pts[0]], rem, 5)
    if comp6 is None or len(comp6) != 6:
        raise AssertionError("自适应覆盖完成性检查失败：原7点骨架应可完成")
    if comp5 is not None:
        raise AssertionError("自适应覆盖完成性检查失败：错误地认为中心+5外圈可保证覆盖")
    print("自适应搜索安全闸：通过（7站可完成，中心+5外圈不可误判为全覆盖）")

    # 角楔交会测试：两条近90°方向应得到小区域
    poly = regular_circumscribed_polygon((0.0, 0.0), ARENA_R, CIRCLE_SIDES)
    poly = clip_wedge(poly, (-500.0, 0.0), 0.0)
    poly = clip_disk_outer(poly, (-500.0, 0.0), MAX_RECEIVE_R)
    poly = clip_wedge(poly, (0.0, -500.0), 90.0)
    poly = clip_disk_outer(poly, (0.0, -500.0), MAX_RECEIVE_R)
    d = polygon_diameter(poly)
    mec = minimum_enclosing_circle(poly)
    print(f"几何测试：交会区域顶点={len(poly)}, 直径≈{d:.3f}m, MEC半径≈{mec.r:.3f}m")
    if not poly or d <= 0 or not math.isfinite(d):
        raise AssertionError("角域交会测试失败")

    # 新增：保证式局部扫清网格测试。40x20m矩形应能由少量清除格完整覆盖。
    dummy_client = SimulatorClient("http://127.0.0.1:1", "self-test", Path(os.devnull))
    dummy_ctl = Problem3Controller(dummy_client)
    test_poly = [(-20.0, -10.0), (20.0, -10.0), (20.0, 10.0), (-20.0, 10.0)]
    test_st = SourceState(channel=99, poly=test_poly)
    clear_grid = dummy_ctl._local_clear_grid(test_st)
    if not clear_grid or len(clear_grid) > LOCAL_CLEAR_MAX_POINTS:
        raise AssertionError("局部保证扫清网格生成失败")
    if any(min(dist(v, q) for q in clear_grid) > CLEAR_R + 1e-7 for v in test_poly):
        raise AssertionError("局部保证扫清网格未覆盖测试多边形顶点")
    print(f"局部保证扫清测试：{len(clear_grid)}个清除格，覆盖通过")

    # request JSON字段不在这里联网测试
    print("自检通过。")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="2026高教社杯B题 问题3 机器狗程序")
    p.add_argument("--robot-id", help="当前登录参赛队号；也可设置环境变量 ROBOT_ID")
    p.add_argument("--base-url", default=os.environ.get("SIM_BASE_URL", "http://127.0.0.1:2026"),
                   help="模拟器接口地址，默认 http://127.0.0.1:2026")
    p.add_argument("--log-dir", default="logs", help="本地日志目录，默认 ./logs")
    p.add_argument("--self-test", action="store_true", help="只做本地几何自检，不连接模拟器")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if args.self_test:
        self_test()
        return 0

    robot_id = args.robot_id or os.environ.get("ROBOT_ID")
    if not robot_id:
        print("错误：请使用 --robot-id <参赛队号> 或设置环境变量 ROBOT_ID", file=sys.stderr)
        return 2

    log_dir = Path(args.log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = log_dir / f"problem3_{stamp}.jsonl"
    summary_path = log_dir / f"problem3_{stamp}_summary.json"

    print(f"模拟器地址: {args.base_url}")
    print(f"robot_id: {robot_id}")
    print(f"本地日志: {log_path}")
    print("注意：请先在模拟器中登录并启动【问题3演练测试/正式测试】，等界面显示接口已就绪后再运行本程序。")

    client = SimulatorClient(args.base_url, robot_id, log_path)
    controller = Problem3Controller(client)
    entered = False
    summary = None
    try:
        er = client.enter()
        entered = True
        print(f"/enter 成功，本局剩余现实运行时间: {er.get('remaining_real_duration_s')} s")
        summary = controller.run()
        print("\n=== 任务完成 ===")
        print(json.dumps(summary, ensure_ascii=False, indent=2))
    except KeyboardInterrupt:
        print("\n用户中断程序。")
    except TimeoutError as e:
        print(f"\n[时间保护退出] {e}")
    except Exception as e:
        print(f"\n[程序异常] {type(e).__name__}: {e}", file=sys.stderr)
        client._log({"event": "fatal", "error_type": type(e).__name__, "error": repr(e)})
    finally:
        if entered and client.real_time_left() > 2.0:
            xr = client.exit()
            if xr:
                print(f"/exit: {xr.get('exit_reason')}")
        if summary is None:
            try:
                summary = controller.summary()
            except Exception:
                summary = {"virtual_time_s": client.last_virtual_time}
        summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"总结已写入: {summary_path}")

    # 清除数若少于发现数，返回非0方便演练识别
    if summary and summary.get("cleared_count", 0) == summary.get("found_count", -1) and summary.get("found_count", 0) >= 10:
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
