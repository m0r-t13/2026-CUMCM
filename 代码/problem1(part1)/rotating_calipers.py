"""旋转卡壳（rotating calipers）求凸多边形直径

核心结论:
  定理 A  凸集的直径必在一对顶点（极点）上取到。
  定理 B  该顶点对必为对踵点对（存在平行支撑线的顶点对）。
  引理 C  沿多边形"旋转"卡尺时，对踵点是单调前进的，故全部对踵点对只有 O(n) 个。

"""

from __future__ import annotations

import math

__all__ = [
    "cross",
    "convex_hull",
    "simplify_polygon",
    "polygon_area",
    "polygon_perimeter",
    "min_edge_length",
    "rotating_calipers_diameter",
    "antipodal_pairs",
    "brute_force_diameter",
    "diameter_pairs",
    "rotate_to_ccw",
]

EPS = 1e-12


# 基础工具


def cross(ox: float, oy: float, ax: float, ay: float, bx: float, by: float) -> float:
    """向量 (oa) × (ob) 的 z 分量。>0 表示 o->a->b 为左转。"""
    return (ax - ox) * (by - oy) - (ay - oy) * (bx - ox)


def dist(p, q) -> float:
    return math.hypot(q[0] - p[0], q[1] - p[1])


def _scale_of(poly) -> float:
    m = 1.0
    for x, y in poly:
        m = max(m, abs(x), abs(y))
    return m


def rotate_to_ccw(poly):
    """把多边形规范为逆时针（有向面积为正）。"""
    s = 0.0
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return list(poly) if s >= 0 else list(reversed(poly))


# 凸包与顶点净化


def convex_hull(points):
    """Andrew 单调链求凸包。

    返回**逆时针**、**严格凸**（去掉共线与重复顶点）的顶点列表。
    旋转卡壳要求输入严格凸，故这一步不能省。
    """
    pts = sorted({(float(x), float(y)) for x, y in points})
    if len(pts) <= 2:
        return pts

    def half(seq):
        out = []
        for p in seq:
            while len(out) >= 2 and cross(
                out[-2][0], out[-2][1], out[-1][0], out[-1][1], p[0], p[1]
            ) <= 0:
                out.pop()
            out.append(p)
        return out

    lower = half(pts)
    upper = half(reversed(pts))
    return lower[:-1] + upper[:-1]


def simplify_polygon(poly, rtol: float = 1e-9):
    """按坐标量级的相对容差合并近重合顶点、去掉近共线顶点。

    必要性：半平面求交得到的候选顶点在退化几何下会出现相差 1e-13 m 量级的
    重合点（例如角楔的两条边界与过顶点的法向直线三线共点）。这类极短边会让
    以顶点为基准的几何量（如张角）失去意义，必须净化。
    """
    pts = [(float(x), float(y)) for x, y in poly]
    n = len(pts)
    if n < 3:
        return pts
    eps = rtol * _scale_of(pts)

    # 1) 合并近重合顶点（环形）
    keep = []
    for p in pts:
        if keep and dist(p, keep[-1]) <= eps:
            continue
        keep.append(p)
    if len(keep) >= 2 and dist(keep[0], keep[-1]) <= eps:
        keep.pop()
    pts = keep
    if len(pts) < 3:
        return pts

    # 2) 迭代去掉近共线顶点
    for _ in range(len(pts)):
        m = len(pts)
        out = []
        for i in range(m):
            a, b, c = pts[(i - 1) % m], pts[i], pts[(i + 1) % m]
            lac = dist(a, c)
            if dist(a, b) <= eps or dist(b, c) <= eps or lac <= eps:
                continue
            h = abs(cross(a[0], a[1], c[0], c[1], b[0], b[1])) / lac
            if h <= eps:
                continue  # b 到直线 ac 的距离可忽略，视为共线
            out.append(b)
        if len(out) == m:
            break
        pts = out
        if len(pts) < 3:
            break
    return pts


def polygon_area(poly) -> float:
    n = len(poly)
    if n < 3:
        return 0.0
    s = 0.0
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return abs(s) / 2.0


def polygon_perimeter(poly) -> float:
    n = len(poly)
    if n < 2:
        return 0.0
    return sum(dist(poly[i], poly[(i + 1) % n]) for i in range(n))


def min_edge_length(poly) -> float:
    n = len(poly)
    if n < 2:
        return 0.0
    return min(dist(poly[i], poly[(i + 1) % n]) for i in range(n))


# 旋转卡壳


def _area2(poly, i: int, j: int, k: int) -> float:
    """三角形 (p_i, p_j, p_k) 面积的两倍。逆时针凸多边形上恒为非负。"""
    ax, ay = poly[i]
    bx, by = poly[j]
    cx, cy = poly[k]
    return abs((bx - ax) * (cy - ay) - (by - ay) * (cx - ax))


def rotating_calipers_diameter(poly):
    """旋转卡壳 O(n) 求凸多边形直径。

    参数
    ----
    poly : 逆时针、严格凸的顶点序列（用 convex_hull 处理过的结果）

    返回
    ----
    (D, i, j) —— 直径长度与达到直径的两个顶点下标

    算法
    ----
    沿每条边 (p_i, p_{i+1}) 维护它的"最远顶点" k，使 p_k 到该边所在直线的距离
    最大（等价于三角形面积最大）。这正是该边对应的对踵顶点。面积随 k 单调
    变化，故 k 只前进不后退，全程推进总次数 <= 2n。
    """
    poly = rotate_to_ccw(list(poly))
    n = len(poly)
    if n == 0:
        return 0.0, -1, -1
    if n == 1:
        return 0.0, 0, 0
    if n == 2:
        return dist(poly[0], poly[1]), 0, 1

    best = -1.0
    bi = bj = 0
    k = 1
    advances = 0
    limit = 2 * n + 4  # 安全上限：k 单调前进，总步数不超过 2n

    for i in range(n):
        ni = (i + 1) % n
        while advances <= limit:
            nk = (k + 1) % n
            if _area2(poly, i, ni, nk) > _area2(poly, i, ni, k):
                k = nk
                advances += 1
            else:
                break
        # 该边的两个候选对踵点对
        for p, q in ((i, k), (ni, k)):
            if p == q:
                continue
            d = dist(poly[p], poly[q])
            if d > best:
                best, bi, bj = d, p, q

    return best, bi, bj


def antipodal_pairs(poly):
    """用旋转卡壳枚举**全部**对踵点对，返回排序去重后的 (i, j) 列表。

    该函数不用于求直径（直径用 rotating_calipers_diameter 即可），而是用于
    验证"对踵点对只有 O(n) 个"这一引理：返回值长度 <= 3n。

    平行边情形：若 p_k 与 p_{k+1} 到边 (p_i,p_{i+1}) 的距离相同（面积相等），
    则 (i,k),(i+1,k),(i,k+1),(i+1,k+1) 都是对踵点对。
    """
    poly = rotate_to_ccw(list(poly))
    n = len(poly)
    if n < 2:
        return []
    if n == 2:
        return [(0, 1)]

    tol = 1e-12 * max(1.0, _scale_of(poly)) ** 2
    pairs = set()
    k = 1
    advances = 0
    limit = 2 * n + 4

    for i in range(n):
        ni = (i + 1) % n
        while advances <= limit:
            nk = (k + 1) % n
            if _area2(poly, i, ni, nk) > _area2(poly, i, ni, k) + tol:
                k = nk
                advances += 1
            else:
                break
        a = _area2(poly, i, ni, k)
        for j in (i, ni):
            pairs.add((min(j, k), max(j, k)))
        # 平行边（等面积）情形
        nk = (k + 1) % n
        if abs(_area2(poly, i, ni, nk) - a) <= tol:
            for j in (i, ni):
                pairs.add((min(j, nk), max(j, nk)))
    return sorted(pairs)


# 交叉校验用：暴力枚举


def brute_force_diameter(poly):
    """O(n^2) 枚举全部顶点对。由定理 A 这是精确的，用作旋转卡壳的对照。"""
    n = len(poly)
    if n == 0:
        return 0.0, -1, -1
    if n == 1:
        return 0.0, 0, 0
    best = -1.0
    bi = bj = 0
    for i in range(n):
        for j in range(i + 1, n):
            d = dist(poly[i], poly[j])
            if d > best:
                best, bi, bj = d, i, j
    return best, bi, bj


def diameter_pairs(poly, tol: float = 1e-9):
    """全部达到最大距离的顶点对（处理直径不唯一的情形）。"""
    dmax, _, _ = brute_force_diameter(poly)
    if dmax <= 0:
        return []
    n = len(poly)
    atol = tol * max(1.0, dmax)
    return [
        (i, j)
        for i in range(n)
        for j in range(i + 1, n)
        if abs(dist(poly[i], poly[j]) - dmax) <= atol
    ]


def _self_test() -> int:
    print("自检：旋转卡壳 vs 暴力枚举")
    pts = [(0, 0), (4, 0), (5, 3), (2, 5), (-1, 3)]
    poly = convex_hull(pts)
    d1, i1, j1 = rotating_calipers_diameter(poly)
    d2, i2, j2 = brute_force_diameter(poly)
    ap = antipodal_pairs(poly)
    print(f"  凸包 {len(poly)} 顶点: {poly}")
    print(f"  卡壳 D = {d1:.12f}  端点 {i1},{j1}")
    print(f"  暴力 D = {d2:.12f}  端点 {i2},{j2}")
    print(f"  对踵点对 {len(ap)} 个（<= 3n = {3 * len(poly)}）: {ap}")
    ok = abs(d1 - d2) < 1e-9 and len(ap) <= 3 * len(poly)
    print("  自检", "通过" if ok else "失败")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(_self_test())
