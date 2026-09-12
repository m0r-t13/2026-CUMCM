#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2026 CUMCM B题 - 问题4 机器狗自动搜索/定位/清除程序

问题4特点：目标区域同时存在全向与定向干扰源；定向源有效覆盖角为180°，方向未知。

策略：
1) 保留“中心 + 12个内环点 + 12个外环点”的25点严格保证网，但将READY/可局部扫清源动态插入剩余搜索路线，搜索与清除交织。
   内环半径950m；外环为外切目标圆的正12边形顶点，半径=1800/cos(15°)。
   所有剖分三角形边长均<1000m，因此任一源都位于一个三个顶点均距源<1000m的三角形内；
   又因源位于三个顶点的凸包中，任意经过源的180°覆盖半平面至少包含一个顶点，故定向源必被发现。
2) 搜索阶段只扫描尚未发现的频道；25个保证点一个不少（除非已发现16源），因此剩余频道判空证明不变。
3) 首次发现后保存示向信息。定位阶段不再假设“1000m内必有信号”，而采用方向感知的
   “左右对称探测 + 失败缩步”策略：从最近一次成功测向点沿示向方向前进，同时向左右偏置，
   若两侧都无信号则按1/2缩小步长，直到重新进入定向覆盖半平面。
4) 每次成功测向仍使用问题1的±1°角楔交会更新保守可行域；MEC半径<=19m时直接/clear。
   若探测进入5m近距区，则直接清除。
5) 全程串行HTTP；网络重试复用同一request_id；记录JSONL日志。

本程序专用于问题4，同时兼容全向源。
"""

from __future__ import annotations

import argparse
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
READY_RADIUS = 19.0          # 清除半径20m，留1m安全余量
CIRCLE_SIDES = 180           # 用外接正多边形保守近似1500m圆盘
MAX_REAL_MARGIN_S = 20.0

# 问题4保证搜索网：25点三角剖分
P4_INNER_R = 950.0
P4_OUTER_R = ARENA_R / math.cos(math.radians(15.0))
P4_RING_N = 12

# 定向源定位：左右对称探测。满尺度位移约284m；双失败后逐次减半。
P4_PROBE_FORWARD = 180.0
P4_PROBE_LATERAL = 30.0
P4_MIN_FORWARD = 40.0
P4_MIN_LATERAL = 12.0

# 搜索站同点复用：不增加移动，只在预计“仍在定向覆盖半平面内且几何收益足够大”时补测。
P4_OPPORTUNISTIC_MAX_PER_STATION = 4
P4_OPPORTUNISTIC_MIN_FORWARD = 35.0
P4_OPPORTUNISTIC_MAX_DISP = 900.0
P4_OPPORTUNISTIC_CONE_BIAS = 95.0
P4_OPPORTUNISTIC_CONE_SLOPE = 0.55
P4_OPPORTUNISTIC_MIN_REL_GAIN = 0.10
P4_OPPORTUNISTIC_INLINE_CLEAR_DIST = 70.0

# 小可行域直接做“定向包围盒网格”保证清除：网格单元中心到任一点最大距离 <20m。
# 与固定7点圆扫相比，细长角楔通常只需2~4个clear点。
P4_CLEAR_GRID_STEP = 26.0          # 单元对角半径 <= 26/sqrt(2) = 18.385m < 20m
P4_CLEAR_GRID_MAX_POINTS = 6       # 点数再多时继续测向通常更划算
P4_CLEAR_GRID_MARGIN = 0.5         # 数值安全余量，仅用于自检说明

MOVE_SPEED_MPS = 5.0
P4_MIN_PROBE_SCALE = 1.0 / 64.0
P4_MAX_SUCCESSFUL_PROBES = 12
P4_MAX_TOTAL_PROBES_PER_SOURCE = 60

# ---- 搜索/清除联合路线（-88s版） ----
# 已READY或可用少量格点保证扫清的源，直接作为真实节点插入剩余保证搜索点路线。
# 未READY源的统计交会位置只作为“软引导”，绝不改变保守可行域或搜索保证。
P4_ROUTE_MULTISTART = 10
P4_FUTURE_ROUTE_WEIGHT = 0.12
P4_FUTURE_ROUTE_CAP_M = 700.0
P4_FUTURE_MIN_OBS = 2
P4_TRIAL_CLEAR_MAX_MEC = 70.0

# mixed-route 仍保留 future 软引导，但不允许为了软目标付出过大的确定性绕路。
# 120m=24s移动时间；只有超过该上限时才退回同一候选集合中的纯距离最短路线。
P4_MIXED_ROUTE_MAX_EXTRA_M = 120.0

# ---- 少干扰源专项模式（仅完整保证搜索后确认11~15源时开启） ----
# 16源场景在阶段A达到上界后仍按原逻辑提前退出；阶段B完全不启用这些专项行为。
P4_SPARSE_MIN_SOURCES = 11
P4_SPARSE_MAX_SOURCES = 15
P4_SPARSE_DP_MAX_N = 15
P4_SPARSE_CROSS_MAX_PER_STOP = 2
P4_SPARSE_CROSS_MIN_SCORE = 0.30
P4_SPARSE_CROSS_MIN_REL_GAIN = 0.16
P4_SPARSE_CROSS_RETRY_DIST_M = 15.0

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


def point_segment_distance(p: Point, a: Point, b: Point) -> float:
    ab = sub(b, a)
    den = dot(ab, ab)
    if den <= 1e-18:
        return dist(p, a)
    t = max(0.0, min(1.0, dot(sub(p, a), ab) / den))
    q = add(a, mul(ab, t))
    return dist(p, q)


def polyline_distance(p: Point, route: Sequence[Point], start: Point) -> float:
    if not route:
        return dist(p, start)
    best = dist(p, start)
    prev = start
    for q in route:
        best = min(best, point_segment_distance(p, prev, q))
        prev = q
    return best


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
# 问题4：25点定向保证搜索网与路径
# ----------------------------
def p4_guarantee_search_points() -> List[Point]:
    """
    返回25个保证搜索点：中心1点、半径950m内环12点、外切目标圆的外环12点。

    三角剖分：
      - 中心与相邻两个内环点构成12个三角形；
      - 每个30°扇区的内外环四边形用一条对角线拆为2个三角形。
    这些三角形覆盖外切正12边形，而目标圆完全包含于该12边形。
    所有三角形边长均小于1000m。
    """
    pts: List[Point] = [(0.0, 0.0)]
    for k in range(P4_RING_N):
        a = math.radians(30.0 * k)
        pts.append((P4_INNER_R * math.cos(a), P4_INNER_R * math.sin(a)))
    for k in range(P4_RING_N):
        a = math.radians(15.0 + 30.0 * k)
        pts.append((P4_OUTER_R * math.cos(a), P4_OUTER_R * math.sin(a)))
    return pts


def open_route_length(route: Sequence[Point], start: Point = (0.0, 0.0)) -> float:
    cur = start
    total = 0.0
    for p in route:
        total += dist(cur, p)
        cur = p
    return total


def nearest_neighbor_route(points: Sequence[Point], start: Point = (0.0, 0.0)) -> List[Point]:
    rem = list(points)
    route: List[Point] = []
    cur = start
    while rem:
        idx = min(range(len(rem)), key=lambda i: dist(cur, rem[i]))
        cur = rem.pop(idx)
        route.append(cur)
    return route


def two_opt_open_route(route: Sequence[Point], start: Point = (0.0, 0.0)) -> List[Point]:
    """固定起点、终点自由的2-opt。"""
    r = list(route)
    improved = True
    while improved:
        improved = False
        n = len(r)
        for i in range(n - 1):
            a = start if i == 0 else r[i - 1]
            b = r[i]
            for j in range(i + 1, n):
                c = r[j]
                dnext = None if j == n - 1 else r[j + 1]
                old = dist(a, b) + (0.0 if dnext is None else dist(c, dnext))
                new = dist(a, c) + (0.0 if dnext is None else dist(b, dnext))
                if new + 1e-9 < old:
                    r[i:j + 1] = reversed(r[i:j + 1])
                    improved = True
                    break
            if improved:
                break
    return r


def p4_search_route() -> List[Point]:
    """从原点开始；原点先扫描，其余24点用NN+2-opt排一个较短开放路径。"""
    pts = p4_guarantee_search_points()
    center = pts[0]
    rest = pts[1:]
    route = two_opt_open_route(nearest_neighbor_route(rest, center), center)
    return [center] + route


def p4_mesh_max_edge() -> float:
    inner_adj = 2.0 * P4_INNER_R * math.sin(math.radians(15.0))
    outer_adj = 2.0 * P4_OUTER_R * math.sin(math.radians(15.0))
    inner_outer = math.sqrt(
        P4_INNER_R ** 2 + P4_OUTER_R ** 2
        - 2.0 * P4_INNER_R * P4_OUTER_R * math.cos(math.radians(15.0))
    )
    return max(P4_INNER_R, inner_adj, outer_adj, inner_outer)


def p4_discrete_directional_gap(samples_r: int = 45, samples_a: int = 360) -> float:
    """
    离散验证：对目标圆内采样点，取1000m内搜索点的方位角，计算最大相邻角缺口。
    若最大缺口<=180°，则任意闭180°覆盖扇区都至少包含一个搜索点。
    """
    pts = p4_guarantee_search_points()
    worst = 0.0
    for ir in range(samples_r + 1):
        rr = ARENA_R * ir / samples_r
        for ia in range(samples_a):
            a = 2.0 * math.pi * ia / samples_a
            g = (rr * math.cos(a), rr * math.sin(a))
            dirs: List[float] = []
            coincident = False
            for q in pts:
                d = dist(g, q)
                if d <= 1e-8:
                    coincident = True
                    break
                if d <= MIN_RECEIVE_R + 1e-8:
                    dirs.append(math.atan2(q[1] - g[1], q[0] - g[0]) % (2.0 * math.pi))
            if coincident:
                gap = 0.0
            else:
                if not dirs:
                    return 360.0
                dirs.sort()
                gaps = [
                    (dirs[(i + 1) % len(dirs)] - dirs[i]) % (2.0 * math.pi)
                    for i in range(len(dirs))
                ]
                gap = math.degrees(max(gaps))
            worst = max(worst, gap)
    return worst

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
        # 运行统计：用于直接核验优化是否真的省下虚拟时间。
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
        if int(channel) != int(old_channel):
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
# 问题4控制器
# ----------------------------
class Problem4Controller:
    def __init__(self, client: SimulatorClient):
        self.client = client
        self.sources: Dict[int, SourceState] = {}
        self.unknown = set(range(1, 21))
        self.empty = set()
        self.cleared = set()
        self.arena_poly = regular_circumscribed_polygon((0.0, 0.0), ARENA_R, CIRCLE_SIDES)
        self.search_route = p4_search_route()
        self.station_scanned: Dict[int, set] = {i: set() for i in range(len(self.search_route))}
        self.search_route_length_m = open_route_length(self.search_route[1:], self.search_route[0])
        self.probe_no_signal_count = 0
        self.opportunistic_measure_count = 0
        self.opportunistic_hit_count = 0
        self.opportunistic_no_signal_count = 0
        self.local_clear_sweep_count = 0
        self.local_clear_sweep_attempts = 0
        self.interleaved_clear_count = 0
        self.trial_clear_count = 0
        self.trial_clear_success_count = 0

        # 少源专项状态：只能在25点保证搜索完成后由_activate_sparse_mode开启。
        self.sparse_mode = False
        self.sparse_cross_measure_count = 0
        self.sparse_cross_hit_count = 0
        self.sparse_cross_no_signal_count = 0
        self.sparse_route_dp_count = 0
        self.mixed_route_guard_count = 0
        self.sparse_cross_attempt_positions: Dict[int, List[Point]] = {}

    def _activate_sparse_mode(self) -> None:
        """完整保证搜索结束后，根据最终发现源数决定是否启用少源专项阶段B。"""
        n = len(self.sources)
        self.sparse_mode = P4_SPARSE_MIN_SOURCES <= n <= P4_SPARSE_MAX_SOURCES
        if self.sparse_mode:
            print(
                f"\n[少源专项模式] 最终发现{n}个干扰源：启用阶段B跨频道同点顺测 + 精确开放路径DP"
            )
        else:
            print(f"\n[常规模式] 最终发现{n}个干扰源：保持原阶段B策略")

    def _ensure_real_time(self) -> None:
        left = self.client.real_time_left()
        if left <= MAX_REAL_MARGIN_S:
            raise TimeoutError(f"现实运行时间只剩 {left:.1f}s，停止新动作")

    def _new_source(self, channel: int, pos: Point, svd_deg: float) -> SourceState:
        st = SourceState(channel=channel, poly=list(self.arena_poly))
        st.add_direction(pos, svd_deg)
        self.sources[channel] = st
        self.unknown.discard(channel)
        print(f"  [发现] 频道{channel}: svd={svd_deg:.2f}°, 初始直径≈{st.diameter:.1f}m")
        return st

    def _handle_measure(self, channel: int, pos: Point, resp: dict) -> str:
        result = resp.get("measure_result")
        if result == "direction":
            svd = float(resp["svd_deg"])
            if channel not in self.sources:
                self._new_source(channel, pos, svd)
            else:
                st = self.sources[channel]
                if not st.cleared:
                    # 避免同一坐标重复叠加完全相同观测；通常不会发生，只做保险。
                    if not st.observations or dist(st.observations[-1].pos, pos) > 1e-8:
                        st.add_direction(pos, svd)
                    print(f"  [测向] 频道{channel}: svd={svd:.2f}°, 直径≈{st.diameter:.1f}m, MEC半径≈{st.mec.r:.1f}m")
            return result

        if result == "near":
            print(f"  [近距] 频道{channel} 距离<=5m，直接清除")
            self.unknown.discard(channel)
            if channel not in self.sources:
                self.sources[channel] = SourceState(channel=channel, poly=[], observations=[])
            self._ensure_real_time()
            cr = self.client.clear(pos, channel)
            if cr.get("clear_result") == "success":
                self.sources[channel].cleared = True
                self.cleared.add(channel)
                print(f"  [清除成功] 频道{channel}")
            else:
                print(f"  [警告] near后清除失败: 频道{channel}, {cr}")
            return result

        if result == "no_signal":
            return result
        raise RuntimeError(f"未知 measure_result={result}")

    def scan_unknown_at(self, pos: Point, station_idx: int) -> None:
        chs = sorted(self.unknown)
        if self.client.current_channel in chs:
            chs.remove(self.client.current_channel)
            chs.insert(0, self.client.current_channel)
        for c in chs:
            self._ensure_real_time()
            resp = self.client.measure(pos, c)
            self.station_scanned[station_idx].add(c)
            self._handle_measure(c, pos, resp)

    def _predict_poly_at_station(self, st: SourceState, pos: Point, measured_deg: float) -> List[Point]:
        """只用于搜索站“同点顺测”的收益估计，不修改真实状态。"""
        p = clip_wedge(st.poly, pos, measured_deg)
        if not p:
            return []
        a = math.radians(measured_deg)
        u = (math.cos(a), math.sin(a))
        min_forward = NEAR_R * math.cos(math.radians(BEARING_ERR_DEG))
        p = clip_halfplane(p, lambda x: dot(u, sub(x, pos)) - min_forward)
        if not p:
            return []
        return clip_disk_outer(p, pos, MAX_RECEIVE_R)

    def _opportunistic_score(self, st: SourceState, pos: Point) -> Optional[Tuple[float, float, float]]:
        """评估搜索站对已发现源的“零额外移动”补测价值。

        返回 (score, predicted_mec_r, relative_gain)。
        这里不把“1000m内”误当成定向源必有信号；只把它作为距离层面的安全项，
        再用“从最近成功点沿示向前进、侧偏较小”的锥形判据过滤定向半平面风险。
        """
        if st.cleared or not st.poly or st.ready or not st.observations:
            return None
        if any(dist(pos, obs.pos) < 1.0 for obs in st.observations):
            return None

        obs = st.observations[-1]
        a = math.radians(obs.svd_deg)
        u = (math.cos(a), math.sin(a))
        v = (-math.sin(a), math.cos(a))
        delta = sub(pos, obs.pos)
        forward = dot(delta, u)
        lateral = abs(dot(delta, v))
        displacement = math.hypot(forward, lateral)

        mec = st.mec
        est_forward = dot(sub(mec.c, obs.pos), u)
        distance_safe = robust_receive(pos, st.poly, MIN_RECEIVE_R)
        near_current_region = dist(pos, mec.c) <= MIN_RECEIVE_R + 50.0
        in_forward_cone = (
            forward >= P4_OPPORTUNISTIC_MIN_FORWARD
            and displacement <= P4_OPPORTUNISTIC_MAX_DISP
            and lateral <= P4_OPPORTUNISTIC_CONE_BIAS + P4_OPPORTUNISTIC_CONE_SLOPE * forward
            and forward <= max(P4_OPPORTUNISTIC_MIN_FORWARD + 20.0, est_forward + 180.0)
        )
        if not (distance_safe or in_forward_cone or near_current_region):
            return None

        # 用当前MEC中心作为“最可能真源”预测示向。这里只决定值不值得尝试；
        # 真正写入可行域的仍然是模拟器返回的svd。
        pred_deg = bearing_deg(pos, mec.c)
        pp = self._predict_poly_at_station(st, pos, pred_deg)
        if not pp:
            return None
        pred_r = minimum_enclosing_circle(pp).r
        cur_r = max(mec.r, 1e-6)
        rel_gain = max(0.0, min(1.0, 1.0 - pred_r / cur_r))
        if rel_gain < P4_OPPORTUNISTIC_MIN_REL_GAIN and pred_r > READY_RADIUS:
            return None

        # 预计直接READY时强奖励；1000m距离保证也略加权。
        score = rel_gain + (1.15 if pred_r <= READY_RADIUS else 0.0) + (0.12 if distance_safe else 0.0)
        return score, pred_r, rel_gain

    def opportunistic_known_at_search_station(self, pos: Point) -> None:
        """搜索移动已经发生：在同一坐标最多补测少量高价值已发现频道。

        这一步的关键是“只增加measure/switch成本，不增加移动成本”。如果一次顺测能替代
        后续一个180m级定向探测，通常净省约30s以上；多个搜索站累计后会大幅削减阶段B。
        """
        cand: List[Tuple[float, int, float, float]] = []
        for c, st in self.sources.items():
            val = self._opportunistic_score(st, pos)
            if val is None:
                continue
            score, pred_r, rel_gain = val
            cand.append((score, c, pred_r, rel_gain))
        cand.sort(reverse=True)
        cand = cand[:P4_OPPORTUNISTIC_MAX_PER_STATION]
        if not cand:
            return

        # 当前频道优先可以少一次切台；其余仍按收益降序。
        ordered = cand
        for i, item in enumerate(cand):
            if item[1] == self.client.current_channel:
                ordered = [item] + cand[:i] + cand[i + 1:]
                break

        for score, c, pred_r, rel_gain in ordered:
            st = self.sources.get(c)
            if st is None or st.cleared or st.ready:
                continue
            # 前面的同点观测不会改变本频道之外的状态，但保险起见重新做一次判定。
            val = self._opportunistic_score(st, pos)
            if val is None:
                continue
            print(
                f"  [同点顺测] 频道{c}: 预计MEC {st.mec.r:.1f}->{pred_r:.1f}m, "
                f"rel_gain≈{rel_gain:.2f}, score≈{score:.2f}"
            )
            self._ensure_real_time()
            self.opportunistic_measure_count += 1
            resp = self.client.measure(pos, c)
            result = resp.get("measure_result")
            if result == "no_signal":
                self.opportunistic_no_signal_count += 1
                print(f"    [顺测无信号] 频道{c}，不改变可行域")
                continue
            self.opportunistic_hit_count += 1
            self._handle_measure(c, pos, resp)
            if c in self.sources and self.sources[c].cleared:
                continue

    def sparse_opportunistic_known_at_current_position(self, exclude_channel: Optional[int]) -> None:
        """少源阶段B：复用已经到达的probe/clear位置，原地顺测其它高价值频道。

        只在完整25点搜索后确认11~15源时启用，因此16源及其它场景不会进入此分支。
        仍复用_opportunistic_score的定向半平面/几何收益过滤，并增加更严格阈值与
        位置去重，避免用大量measure换取很小的定位收益。
        """
        if not self.sparse_mode:
            return
        pos = self.client.current_pos
        cand: List[Tuple[float, int, float, float]] = []
        for c, st in self.sources.items():
            if c == exclude_channel or st.cleared or not st.poly or st.ready:
                continue
            attempted = self.sparse_cross_attempt_positions.get(c, [])
            if any(dist(pos, q) < P4_SPARSE_CROSS_RETRY_DIST_M for q in attempted):
                continue
            val = self._opportunistic_score(st, pos)
            if val is None:
                continue
            score, pred_r, rel_gain = val
            # 预计一次就READY时无条件保留；否则要求比搜索站顺测更高的收益。
            if pred_r > READY_RADIUS and (
                score < P4_SPARSE_CROSS_MIN_SCORE or rel_gain < P4_SPARSE_CROSS_MIN_REL_GAIN
            ):
                continue
            cand.append((score, c, pred_r, rel_gain))

        if not cand:
            return
        cand.sort(reverse=True)
        cand = cand[:P4_SPARSE_CROSS_MAX_PER_STOP]
        # 当前频道优先，尽量减少切台。
        for i, item in enumerate(cand):
            if item[1] == self.client.current_channel:
                cand = [item] + cand[:i] + cand[i + 1:]
                break

        for score, c, pred_r, rel_gain in cand:
            st = self.sources.get(c)
            if st is None or st.cleared or st.ready:
                continue
            val = self._opportunistic_score(st, pos)
            if val is None:
                continue
            score2, pred_r2, rel_gain2 = val
            if pred_r2 > READY_RADIUS and (
                score2 < P4_SPARSE_CROSS_MIN_SCORE or rel_gain2 < P4_SPARSE_CROSS_MIN_REL_GAIN
            ):
                continue
            self._ensure_real_time()
            self.sparse_cross_attempt_positions.setdefault(c, []).append(pos)
            self.sparse_cross_measure_count += 1
            print(
                f"    [少源跨频道顺测] 频道{c}: 预计MEC {st.mec.r:.1f}->{pred_r2:.1f}m, "
                f"rel_gain≈{rel_gain2:.2f}"
            )
            resp = self.client.measure(pos, c)
            result = resp.get("measure_result")
            if result == "no_signal":
                self.sparse_cross_no_signal_count += 1
                continue
            self.sparse_cross_hit_count += 1
            self._handle_measure(c, pos, resp)

    def _oriented_clear_grid(self, st: SourceState, order_from: Optional[Point] = None) -> List[Point]:
        """用可行域的主轴方向构造保证覆盖的矩形网格。

        先取凸多边形直径端点作为u轴，再取垂直v轴；多边形必包含在这个
        旋转包围矩形中。每个网格单元边长<=P4_CLEAR_GRID_STEP，clear点放在
        单元中心，因此矩形内任一点到某个clear点的距离<=step/sqrt(2)<20m。
        这给出严格的局部清除覆盖，同时对细长角楔比固定7点圆扫省很多动作。
        """
        poly = st.poly
        if not poly:
            return []
        if len(poly) == 1:
            return [poly[0]]

        # 直径端点（顶点数通常很小，O(n^2)足够且稳定）。
        ia, ib, best = 0, 1, -1.0
        for i in range(len(poly)):
            for j in range(i + 1, len(poly)):
                d = dist(poly[i], poly[j])
                if d > best:
                    ia, ib, best = i, j, d
        if best <= 1e-9:
            return [centroid(poly)]

        a, b = poly[ia], poly[ib]
        u = ((b[0] - a[0]) / best, (b[1] - a[1]) / best)
        v = (-u[1], u[0])
        pu = [dot(p, u) for p in poly]
        pv = [dot(p, v) for p in poly]
        u0, u1 = min(pu), max(pu)
        v0, v1 = min(pv), max(pv)
        lu, lv = u1 - u0, v1 - v0

        nu = max(1, int(math.ceil(lu / P4_CLEAR_GRID_STEP - 1e-12)))
        nv = max(1, int(math.ceil(lv / P4_CLEAR_GRID_STEP - 1e-12)))
        if nu * nv > P4_CLEAR_GRID_MAX_POINTS:
            return []

        du = lu / nu if nu else 0.0
        dv = lv / nv if nv else 0.0
        pts: List[Point] = []
        for i in range(nu):
            cu = (u0 + u1) * 0.5 if lu <= 1e-12 else u0 + (i + 0.5) * du
            for j in range(nv):
                cv = (v0 + v1) * 0.5 if lv <= 1e-12 else v0 + (j + 0.5) * dv
                pts.append((u[0] * cu + v[0] * cv, u[1] * cu + v[1] * cv))

        startp = self.client.current_pos if order_from is None else order_from
        return two_opt_open_route(nearest_neighbor_route(pts, startp), startp)

    def try_local_clear_sweep(self, channel: int) -> bool:
        """当当前可行域可被<=6个20m清除圆保证覆盖时，直接局部扫清。

        这里按“实际可行多边形”而不是MEC圆盘铺点；角楔交会后的区域通常很细长，
        因而常见2~4点即可保证覆盖。若所需点数>阈值则返回False，继续定向测量。
        """
        st = self.sources[channel]
        if st.cleared or not st.poly:
            return st.cleared
        if st.ready:
            return self.try_clear_source(channel)

        pts = self._oriented_clear_grid(st)
        if not pts:
            return False

        self.local_clear_sweep_count += 1
        print(
            f"  [自适应保证清除] 频道{channel}: 直径≈{st.diameter:.2f}m, "
            f"网格点={len(pts)}（单元边长<={P4_CLEAR_GRID_STEP:.1f}m）"
        )
        for i, q in enumerate(pts, 1):
            self._ensure_real_time()
            st.clear_attempts += 1
            self.local_clear_sweep_attempts += 1
            resp = self.client.clear(q, channel)
            if resp.get("clear_result") == "success":
                st.cleared = True
                self.cleared.add(channel)
                print(f"    [局部清除成功] 频道{channel}, 第{i}/{len(pts)}点")
                self.sparse_opportunistic_known_at_current_position(exclude_channel=channel)
                return True
            self.sparse_opportunistic_known_at_current_position(exclude_channel=channel)
        # 理论上不会到这里；保留测向回退，避免浮点/模拟器边界差异使程序中断。
        print(f"    [局部清除未命中] 频道{channel}，回退到定向追踪")
        return False

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
            self.sparse_opportunistic_known_at_current_position(exclude_channel=channel)
            return True
        self.sparse_opportunistic_known_at_current_position(exclude_channel=channel)
        print(f"  [清除失败] 频道{channel}，继续缩小定位区域；响应={resp.get('clear_result')}")
        return False

    def _least_squares_source_point(self, st: SourceState) -> Optional[Point]:
        """由>=2条测向中心线做最小二乘交会；仅用于路线软引导/廉价试探，不参与保守定位。"""
        if len(st.observations) < P4_FUTURE_MIN_OBS:
            return None
        # 对每条线，法向n=(-uy,ux)，最小化 sum (n·(x-p))^2
        a11 = a12 = a22 = b1 = b2 = 0.0
        for obs in st.observations:
            ang = math.radians(obs.svd_deg)
            ux, uy = math.cos(ang), math.sin(ang)
            nx, ny = -uy, ux
            rhs = nx * obs.pos[0] + ny * obs.pos[1]
            a11 += nx * nx
            a12 += nx * ny
            a22 += ny * ny
            b1 += nx * rhs
            b2 += ny * rhs
        det = a11 * a22 - a12 * a12
        if abs(det) < 1e-8:
            return None
        x = (b1 * a22 - b2 * a12) / det
        y = (a11 * b2 - a12 * b1) / det
        q = (x, y)
        rr = math.hypot(x, y)
        if rr > ARENA_R:
            q = (x * ARENA_R / rr, y * ARENA_R / rr)
        return q

    def _source_action_anchor(self, channel: int) -> Optional[Tuple[str, Point]]:
        """返回当前可立即执行的源任务锚点；None表示继续依赖搜索站免费补测。"""
        st = self.sources[channel]
        if st.cleared or not st.poly:
            return None
        if st.ready:
            return ("clear", st.mec.c)
        grid = self._oriented_clear_grid(st, self.client.current_pos)
        if grid:
            return ("sweep", grid[0])
        return None

    def _future_prediction_points(self) -> List[Point]:
        pts: List[Point] = []
        for c, st in self.sources.items():
            if st.cleared or not st.poly or st.ready:
                continue
            if self._oriented_clear_grid(st, self.client.current_pos):
                continue
            q = self._least_squares_source_point(st)
            if q is not None:
                pts.append(q)
        return pts

    def _mixed_route_score(self, nodes: Sequence[Tuple[str, int, Point]], start: Point) -> float:
        pts = [x[2] for x in nodes]
        score = open_route_length(pts, start)
        future = self._future_prediction_points()
        if future and pts:
            penalty = 0.0
            for q in future:
                penalty += min(P4_FUTURE_ROUTE_CAP_M, polyline_distance(q, pts, start))
            score += P4_FUTURE_ROUTE_WEIGHT * penalty
        return score

    def _two_opt_nodes(self, nodes: Sequence[Tuple[str, int, Point]], start: Point) -> List[Tuple[str, int, Point]]:
        r = list(nodes)
        improved = True
        while improved:
            improved = False
            n = len(r)
            for i in range(n - 1):
                a = start if i == 0 else r[i - 1][2]
                b = r[i][2]
                for j in range(i + 1, n):
                    c = r[j][2]
                    dnext = None if j == n - 1 else r[j + 1][2]
                    old = dist(a, b) + (0.0 if dnext is None else dist(c, dnext))
                    new = dist(a, c) + (0.0 if dnext is None else dist(b, dnext))
                    if new + 1e-9 < old:
                        r[i:j + 1] = reversed(r[i:j + 1])
                        improved = True
                        break
                if improved:
                    break
        return r

    def _select_mixed_candidate(
        self,
        candidates: Sequence[List[Tuple[str, int, Point]]],
        start: Point,
    ) -> List[Tuple[str, int, Point]]:
        """future软引导 + 硬绕路保护。

        先按原_mixed_route_score选择；若它相对同一候选集合的纯距离最短路线
        多走超过P4_MIXED_ROUTE_MAX_EXTRA_M，则退回纯距离路线。
        """
        if not candidates:
            return []
        shortest = min(candidates, key=lambda r: open_route_length([x[2] for x in r], start))
        guided = min(candidates, key=lambda r: self._mixed_route_score(r, start))
        shortest_len = open_route_length([x[2] for x in shortest], start)
        guided_len = open_route_length([x[2] for x in guided], start)
        if guided_len > shortest_len + P4_MIXED_ROUTE_MAX_EXTRA_M + 1e-9:
            self.mixed_route_guard_count += 1
            return shortest
        return guided

    def _plan_mixed_route(self, search_nodes: Sequence[Tuple[int, Point]]) -> List[Tuple[str, int, Point]]:
        """把剩余保证搜索点和当前可清源统一做开放式多起点2-opt。"""
        cur = self.client.current_pos
        nodes: List[Tuple[str, int, Point]] = [("search", idx, p) for idx, p in search_nodes]
        for c in sorted(self.sources):
            act = self._source_action_anchor(c)
            if act is not None:
                kind, q = act
                nodes.append((kind, c, q))
        if not nodes:
            return []

        # NN基线
        rem = list(nodes)
        nn: List[Tuple[str, int, Point]] = []
        cp = cur
        while rem:
            k = min(range(len(rem)), key=lambda i: dist(cp, rem[i][2]))
            item = rem.pop(k)
            nn.append(item)
            cp = item[2]
        candidates = [self._two_opt_nodes(nn, cur)]

        # 强制多个不同首节点，多起点能显著减少2-opt局部最优。
        firsts = sorted(nodes, key=lambda x: dist(cur, x[2]))[:min(P4_ROUTE_MULTISTART, len(nodes))]
        for first in firsts:
            rem = [x for x in nodes if x is not first]
            route = [first]
            cp = first[2]
            while rem:
                k = min(range(len(rem)), key=lambda i: dist(cp, rem[i][2]))
                item = rem.pop(k)
                route.append(item)
                cp = item[2]
            candidates.append(self._two_opt_nodes(route, cur))
        return self._select_mixed_candidate(candidates, cur)

    def _maybe_trial_clear(self, channel: int) -> bool:
        """统计交会点的廉价试探清除；失败不改变可行域，成功可省后续移动。"""
        st = self.sources[channel]
        if st.cleared or not st.poly or len(st.observations) < 2 or st.mec.r > P4_TRIAL_CLEAR_MAX_MEC:
            return False
        q = self._least_squares_source_point(st)
        if q is None:
            return False
        # 只在估计点位于保守MEC附近时尝试，防止病态交会导致大绕路。
        if dist(q, st.mec.c) > max(80.0, 1.8 * st.mec.r):
            return False
        self._ensure_real_time()
        self.trial_clear_count += 1
        st.clear_attempts += 1
        print(f"  [交会点试探清除] 频道{channel} -> ({q[0]:.1f},{q[1]:.1f}), MEC≈{st.mec.r:.1f}m")
        resp = self.client.clear(q, channel)
        if resp.get("clear_result") == "success":
            st.cleared = True
            self.cleared.add(channel)
            self.trial_clear_success_count += 1
            print(f"  [试探清除成功] 频道{channel}")
            self.sparse_opportunistic_known_at_current_position(exclude_channel=channel)
            return True
        self.sparse_opportunistic_known_at_current_position(exclude_channel=channel)
        return False

    def guaranteed_directional_search(self) -> None:
        print("\n=== 阶段A：25点保证搜索 + 顺路清除联合路线 ===")
        print(
            f"搜索网最大剖分边长≈{p4_mesh_max_edge():.1f}m < {MIN_RECEIVE_R:.0f}m，"
            f"原固定搜索路线长度≈{self.search_route_length_m:.1f}m"
        )

        # 原点固定先扫；剩余24点一个不少，严格保证完全保留。
        center = self.search_route[0]
        self.scan_unknown_at(center, 0)
        self.opportunistic_known_at_search_station(center)
        remaining: List[Tuple[int, Point]] = [(i, p) for i, p in enumerate(self.search_route[1:], start=1)]

        guard = 0
        while remaining:
            guard += 1
            if guard > 200:
                raise RuntimeError("问题4联合搜索循环超过200次")
            if len(self.sources) >= 16:
                print("已发现16个频道，达到题目源数量上界；剩余固定搜索点无需再扫未知频道。")
                self.empty |= self.unknown
                self.unknown.clear()
                break

            route = self._plan_mixed_route(remaining)
            if not route:
                break
            kind, ident, q = route[0]

            if kind == "search":
                idx = ident
                print(f"\n[联合保证搜索点] idx={idx}, ({q[0]:.1f},{q[1]:.1f}), 尚未知={len(self.unknown)}")
                self.scan_unknown_at(q, idx)
                self.opportunistic_known_at_search_station(q)
                remaining = [(j, p) for j, p in remaining if j != idx]
                continue

            c = ident
            st = self.sources.get(c)
            if st is None or st.cleared:
                continue
            self.interleaved_clear_count += 1
            print(f"\n[搜索途中顺路处理] 频道{c}, 类型={kind}, ({q[0]:.1f},{q[1]:.1f})")
            if kind == "clear":
                self.try_clear_source(c)
            else:
                self.try_local_clear_sweep(c)

        # 若因为达到16源提前退出，unknown已判空；否则必须完整走完25点才能判空。
        if self.unknown:
            if remaining:
                # 只有未达到16源时才会到这里；理论上remaining应为空。
                raise RuntimeError("问题4保证搜索尚未完成却尝试判空")
            print(f"\n25点保证搜索完成，以下频道始终未发现，判定为空频道: {sorted(self.unknown)}")
            self.empty |= self.unknown
            self.unknown.clear()

        # 只有此时才能严格知道最终源数；此前不能为了“疑似少源”改变16源搜索行为。
        self._activate_sparse_mode()

    def _probe_offsets(self, st: SourceState, scale: float) -> Tuple[float, float]:
        obs = st.observations[-1]
        a = math.radians(obs.svd_deg)
        u = (math.cos(a), math.sin(a))
        est_forward = dot(sub(st.mec.c, obs.pos), u) if st.poly else P4_PROBE_FORWARD
        # 远距离保持180m前进；接近源时按估计前向距离自动缩短，减少越过真源导致的双侧无信号。
        f0 = min(P4_PROBE_FORWARD, max(P4_MIN_FORWARD, 0.55 * max(est_forward, 0.0)))
        # 侧偏只负责制造交会基线；窄侧偏比原±180m显著省路，同时更不易跨出定向半平面。
        b0 = min(P4_PROBE_LATERAL, max(P4_MIN_LATERAL, 0.18 * f0))
        return f0 * scale, b0 * scale

    def _paired_probe_points(self, st: SourceState, scale: float) -> List[Point]:
        obs = st.observations[-1]
        a = math.radians(obs.svd_deg)
        u = (math.cos(a), math.sin(a))
        v = (-math.sin(a), math.cos(a))
        f, b = self._probe_offsets(st, scale)
        base = add(obs.pos, mul(u, f))
        return [add(base, mul(v, b)), add(base, mul(v, -b))]

    def _choose_next_localize_channel(self) -> Optional[int]:
        """阶段B按所有源“下一实际动作点”做多起点2-opt，而不是单纯最近邻。"""
        cur = self.client.current_pos
        nodes: List[Tuple[str, int, Point]] = []
        for c, st in self.sources.items():
            if st.cleared:
                continue
            if st.poly and st.ready:
                q = st.mec.c
            elif st.poly:
                grid = self._oriented_clear_grid(st, cur)
                if grid:
                    q = grid[0]
                elif st.observations:
                    q = min(self._paired_probe_points(st, 1.0), key=lambda x: dist(cur, x))
                else:
                    q = cur
            elif st.observations:
                q = st.observations[-1].pos
            else:
                q = cur
            nodes.append(("source", c, q))
        if not nodes:
            return None
        # 重用混合路线的多起点2-opt；这里没有search节点。
        route = self._plan_nodes_only(nodes, cur)
        return route[0][1]

    def _exact_open_route_dp(
        self,
        nodes: Sequence[Tuple[str, int, Point]],
        start: Point,
    ) -> List[Tuple[str, int, Point]]:
        """Held-Karp精确求固定起点、终点自由的开放Hamilton路径。

        少源阶段最多15个待处理源，状态数n*2^n可控；只用于选择当前全局更好的
        下一源，源状态变化后会重新规划。
        """
        n = len(nodes)
        if n <= 1:
            return list(nodes)
        full = (1 << n) - 1
        states = (1 << n) * n
        inf = float("inf")
        dp = [inf] * states
        parent = [-1] * states
        start_dist = [dist(start, nodes[j][2]) for j in range(n)]
        dmat = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                d = dist(nodes[i][2], nodes[j][2])
                dmat[i][j] = d
                dmat[j][i] = d
        for j in range(n):
            dp[(1 << j) * n + j] = start_dist[j]

        for mask in range(1, full + 1):
            base = mask * n
            rem_mask_all = full ^ mask
            if rem_mask_all == 0:
                continue
            bits = mask
            while bits:
                bit_j = bits & -bits
                j = bit_j.bit_length() - 1
                bits ^= bit_j
                cur = dp[base + j]
                if cur == inf:
                    continue
                rem = rem_mask_all
                while rem:
                    bit_k = rem & -rem
                    k = bit_k.bit_length() - 1
                    rem ^= bit_k
                    nmask = mask | bit_k
                    idx = nmask * n + k
                    nd = cur + dmat[j][k]
                    if nd + 1e-12 < dp[idx]:
                        dp[idx] = nd
                        parent[idx] = j

        full_base = full * n
        end = min(range(n), key=lambda j: dp[full_base + j])
        order_rev: List[int] = []
        mask = full
        while end >= 0:
            order_rev.append(end)
            prev = parent[mask * n + end]
            mask ^= 1 << end
            end = prev
        order_rev.reverse()
        return [nodes[i] for i in order_rev]

    def _plan_nodes_only(self, nodes: Sequence[Tuple[str, int, Point]], start: Point) -> List[Tuple[str, int, Point]]:
        if not nodes:
            return []
        if self.sparse_mode and len(nodes) <= P4_SPARSE_DP_MAX_N:
            self.sparse_route_dp_count += 1
            return self._exact_open_route_dp(nodes, start)
        candidates: List[List[Tuple[str, int, Point]]] = []
        firsts = sorted(nodes, key=lambda x: dist(start, x[2]))[:min(P4_ROUTE_MULTISTART, len(nodes))]
        for first in firsts:
            rem = [x for x in nodes if x is not first]
            r = [first]
            cp = first[2]
            while rem:
                k = min(range(len(rem)), key=lambda i: dist(cp, rem[i][2]))
                item = rem.pop(k)
                r.append(item)
                cp = item[2]
            candidates.append(self._two_opt_nodes(r, start))
        return min(candidates, key=lambda r: open_route_length([x[2] for x in r], start))

    def localize_one_source(self, channel: int) -> None:
        st = self.sources[channel]
        if st.cleared:
            return
        if st.poly and st.ready and self.try_clear_source(channel):
            return
        if st.poly and self.try_local_clear_sweep(channel):
            return
        if st.poly and self._maybe_trial_clear(channel):
            return
        if not st.observations:
            raise RuntimeError(f"频道{channel}没有可用于定向追踪的成功测向记录")

        successful = max(0, len(st.observations) - 1)
        total_probes = 0
        print(f"\n[定向追踪] 频道{channel}: 从最近成功测向点开始左右对称探测")

        while not st.cleared:
            if total_probes >= P4_MAX_TOTAL_PROBES_PER_SOURCE:
                raise RuntimeError(f"频道{channel}定向追踪探测次数超过{P4_MAX_TOTAL_PROBES_PER_SOURCE}")
            if successful >= P4_MAX_SUCCESSFUL_PROBES:
                print(f"  [告警] 频道{channel}已有{len(st.observations)}次成功测向仍未READY，继续追踪")

            # 如果此前已经达到安全清除半径，优先清除；失败则继续增加观测。
            if st.poly and st.ready:
                if self.try_clear_source(channel):
                    return
            if st.poly and self.try_local_clear_sweep(channel):
                return
            if st.poly and self._maybe_trial_clear(channel):
                return

            scale = 1.0
            got_signal = False
            while scale >= P4_MIN_PROBE_SCALE - 1e-12:
                qlist = self._paired_probe_points(st, scale)
                qlist.sort(key=lambda q: dist(self.client.current_pos, q))
                f, b = self._probe_offsets(st, scale)
                print(f"  [成对探测] scale={scale:.5f}, 前进≈{f:.1f}m, 侧偏≈±{b:.1f}m")

                for q in qlist:
                    self._ensure_real_time()
                    total_probes += 1
                    resp = self.client.measure(q, channel)
                    result = resp.get("measure_result")
                    if result == "no_signal":
                        self.probe_no_signal_count += 1
                        print(f"    [无信号] ({q[0]:.2f},{q[1]:.2f})，尝试对称点/缩步")
                        self.sparse_opportunistic_known_at_current_position(exclude_channel=channel)
                        continue

                    self._handle_measure(channel, q, resp)
                    self.sparse_opportunistic_known_at_current_position(exclude_channel=channel)
                    got_signal = True
                    if result == "near" or self.sources[channel].cleared:
                        return
                    successful += 1
                    st = self.sources[channel]
                    if st.poly and st.ready and self.try_clear_source(channel):
                        return
                    if st.poly and self.try_local_clear_sweep(channel):
                        return
                    break

                if got_signal:
                    break
                scale *= 0.5
                print(f"    [双侧均无信号] 缩小探测尺度至 {scale:.5f}")

            if not got_signal:
                # 理论上最后成功点仍在覆盖半平面内，足够小的对称扰动应能重新获得信号。
                # 若模拟器数值边界导致失败，回到最后成功点重测一次以恢复状态，再继续。
                obs = st.observations[-1]
                print("  [恢复探测] 极小尺度仍无信号，回到最近成功点重新测量")
                self._ensure_real_time()
                resp = self.client.measure(obs.pos, channel)
                result = resp.get("measure_result")
                if result == "no_signal":
                    raise RuntimeError(f"频道{channel}在此前成功测向点重新测量却无信号，状态不一致")
                self._handle_measure(channel, obs.pos, resp)
                self.sparse_opportunistic_known_at_current_position(exclude_channel=channel)
                if result == "near" or self.sources[channel].cleared:
                    return

    def localize_and_clear_all(self) -> None:
        print("\n=== 阶段B：方向感知定位与清除 ===")
        guard = 0
        while True:
            guard += 1
            if guard > 100:
                raise RuntimeError("问题4任务循环超过100次，防止死循环")
            c = self._choose_next_localize_channel()
            if c is None:
                break
            self.localize_one_source(c)

    def summary(self) -> dict:
        found_count = len(self.sources)
        cleared_count = len(self.cleared)
        real_elapsed = None
        if self.client.real_start_monotonic is not None:
            real_elapsed = time.monotonic() - self.client.real_start_monotonic
        return {
            "problem": 4,
            "found_channels": sorted(self.sources.keys()),
            "found_count": found_count,
            "cleared_channels": sorted(self.cleared),
            "cleared_count": cleared_count,
            "clear_rate_vs_found": (cleared_count / found_count) if found_count else 0.0,
            "average_virtual_localize_clear_time_s": (self.client.last_virtual_time / cleared_count) if cleared_count else None,
            "empty_channels": sorted(self.empty),
            "search_strategy": {
                "mode": "25_point_triangulation_directional_guarantee",
                "search_points": len(self.search_route),
                "mesh_max_edge_m": p4_mesh_max_edge(),
                "fixed_route_length_m": self.search_route_length_m,
            },
            "localize_strategy": {
                "mode": "symmetric_pair_probe_with_halving_and_wedge_MEC",
                "probe_forward_m": P4_PROBE_FORWARD,
                "probe_lateral_m": P4_PROBE_LATERAL,
                "no_signal_probe_count": self.probe_no_signal_count,
                "interleaved_clear_count": self.interleaved_clear_count,
                "trial_clear_count": self.trial_clear_count,
                "trial_clear_success_count": self.trial_clear_success_count,
                "sparse_mode": self.sparse_mode,
                "sparse_cross_measure_count": self.sparse_cross_measure_count,
                "sparse_cross_hit_count": self.sparse_cross_hit_count,
                "sparse_cross_no_signal_count": self.sparse_cross_no_signal_count,
                "sparse_route_dp_count": self.sparse_route_dp_count,
                "mixed_route_guard_count": self.mixed_route_guard_count,
            },
            "route_metrics": {
                "move_distance_m": self.client.total_move_distance_m,
                "estimated_move_time_s": self.client.total_move_distance_m / MOVE_SPEED_MPS,
                "measure_count": self.client.measure_count,
                "channel_switch_count": self.client.channel_switch_count,
                "clear_count": self.client.clear_count,
                "clear_success_count": self.client.clear_success_count,
                "clear_fail_count": self.client.clear_fail_count,
            },
            "virtual_time_s": self.client.last_virtual_time,
            "local_real_elapsed_s": real_elapsed,
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
        self.guaranteed_directional_search()
        self.localize_and_clear_all()
        return self.summary()


# ----------------------------
# 自检
# ----------------------------
def self_test() -> None:
    print("运行问题4本地自检……")

    pts = p4_guarantee_search_points()
    if len(pts) != 25:
        raise AssertionError(f"问题4保证搜索点数量异常: {len(pts)}")
    max_edge = p4_mesh_max_edge()
    print(f"25点三角网：点数={len(pts)}, 最大剖分边长≈{max_edge:.3f}m")
    if max_edge >= MIN_RECEIVE_R - 1e-9:
        raise AssertionError("问题4三角网边长未严格小于1000m")

    gap = p4_discrete_directional_gap(samples_r=40, samples_a=360)
    print(f"定向覆盖离散验证：1000m邻域方向最大角缺口≈{gap:.6f}°")
    if gap > 180.000001:
        raise AssertionError("问题4定向180°保证覆盖离散验证失败")

    route = p4_search_route()
    rlen = open_route_length(route[1:], route[0])
    print(f"固定保证搜索路线长度≈{rlen:.1f}m")

    # 保留问题1角域交会测试。
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

    # 对称探测几何参数检查：满尺度候选相对成功点位移不应过大。
    disp = math.hypot(P4_PROBE_FORWARD, P4_PROBE_LATERAL)
    print(f"定向追踪测试：满尺度单次试探位移≈{disp:.1f}m，最小尺度={P4_MIN_PROBE_SCALE:.5f}")
    if disp >= MIN_RECEIVE_R:
        raise AssertionError("成对探测位移过大")

    print("自检通过。")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="2026高教社杯B题 问题4 全向/定向干扰源机器狗程序")
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
    log_path = log_dir / f"problem4_{stamp}.jsonl"
    summary_path = log_dir / f"problem4_{stamp}_summary.json"

    print(f"模拟器地址: {args.base_url}")
    print(f"robot_id: {robot_id}")
    print(f"本地日志: {log_path}")
    print("注意：请先在模拟器中登录并启动【问题4演练测试/正式测试】，等界面显示接口已就绪后再运行本程序。")

    client = SimulatorClient(args.base_url, robot_id, log_path)
    controller = Problem4Controller(client)
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
