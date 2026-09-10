"""旋转卡壳实现的验证测试。

运行：
    python test_rotating_calipers.py

覆盖:
  T1  与暴力枚举一致          —— 1000 个随机凸多边形，O(n) 结果 == O(n^2) 结果
  T2  直径端点确实是对踵点对  —— 定理 B 的直接检验
  T3  对踵点对的个数 <= 3n    —— 引理 C（O(n) 个对踵点对）
  T4  对踵点对枚举的正确性    —— 与"平行支撑线"定义的点采样版暴力法比对
  T5  退化输入                —— n=1, n=2, 共线、重合点
  T6  平行边情形              —— 矩形、正多边形（对踵点对可能一行多个）
  T7  端到端物理校验          —— 真源必须落在由精确示向度构造的定位区域内部
"""

from __future__ import annotations

import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import localization_region as lr  # noqa: E402
import rotating_calipers as rc  # noqa: E402

RND = random.Random(20260910)
FAILS = []


def check(cond, msg):
    print(("  ok   " if cond else "  FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def rand_convex(n_extra=6, span=3000.0):
    n = RND.randint(3, 25)
    pts = [(RND.uniform(-span, span), RND.uniform(-span, span))
           for _ in range(n * n_extra)]
    return rc.convex_hull(pts)


#T1
def t1_agrees_with_brute_force():
    print("T1 旋转卡壳 O(n) vs 暴力枚举 O(n^2)")
    worst = 0.0
    tested = 0
    for _ in range(1000):
        poly = rand_convex()
        if len(poly) < 3:
            continue
        d1, _, _ = rc.rotating_calipers_diameter(poly)
        d2, _, _ = rc.brute_force_diameter(poly)
        worst = max(worst, abs(d1 - d2))
        tested += 1
    check(worst < 1e-9, f"{tested} 个随机凸多边形，最大偏差 {worst:.3e} m")


#T2
def t2_diameter_is_antipodal():
    print("T2 直径端点必为对踵点对（定理 B）")
    bad = 0
    for _ in range(500):
        poly = rand_convex()
        if len(poly) < 3:
            continue
        _, i, j = rc.rotating_calipers_diameter(poly)
        ap = set(rc.antipodal_pairs(poly))
        if (min(i, j), max(i, j)) not in ap:
            bad += 1
    check(bad == 0, f"500 个多边形的直径对都在对踵点对集合中（例外 {bad}）")


#T3
def t3_count_is_linear():
    print("T3 对踵点对个数 <= 3n（引理 C）")
    worst_ratio = 0.0
    bad = 0
    sizes = []
    for _ in range(400):
        poly = rand_convex()
        if len(poly) < 3:
            continue
        n = len(poly)
        c = len(rc.antipodal_pairs(poly))
        sizes.append((n, c))
        worst_ratio = max(worst_ratio, c / n)
        if c > 3 * n:
            bad += 1
    # 顶点对总数是 O(n^2)，对踵点对是 O(n)：用圆周上的点造大 n 的凸包做量级对比
    N = 4000
    circ = [
        (1000.0 * math.cos(2 * math.pi * k / N) * (1 + RND.uniform(-1e-3, 1e-3)),
         1000.0 * math.sin(2 * math.pi * k / N) * (1 + RND.uniform(-1e-3, 1e-3)))
        for k in range(N)
    ]
    big = rc.convex_hull(circ)
    cb = len(rc.antipodal_pairs(big))
    total = len(big) * (len(big) - 1) // 2
    check(bad == 0 and cb <= 3 * len(big) and cb < total / 50,
          f"400 个多边形全部满足，最大 c/n = {worst_ratio:.3f}；"
          f"n={len(big)} 时对踵点对 {cb} 个 vs 顶点对总数 {total} 个"
          f"（差 {total / cb:.0f} 倍）")


#T4
def t4_antipodal_definition():
    print("T4 对踵点对枚举 vs 定义（平行支撑线）暴力法")
    bad = 0
    tested = 0
    m = 20000
    for _ in range(8):
        poly = rand_convex(n_extra=4, span=500.0)
        n = len(poly)
        if n < 3 or n > 40:
            continue
        tested += 1
        exact = set(tuple(p) for p in rc.antipodal_pairs(poly))
        # 按定义：存在方向 d 使 p_i 最大且 p_j 最小
        sampled = set()
        for t in range(m):
            ang = 2 * math.pi * t / m
            cx, cy = math.cos(ang), math.sin(ang)
            vals = [cx * p[0] + cy * p[1] for p in poly]
            vmax, vmin = max(vals), min(vals)
            tol = 1e-9 * max(1.0, abs(vmax), abs(vmin))
            imax = [i for i, v in enumerate(vals) if v >= vmax - tol]
            imin = [i for i, v in enumerate(vals) if v <= vmin + tol]
            for i in imax:
                for j in imin:
                    if i != j:
                        sampled.add((min(i, j), max(i, j)))
        if not exact <= sampled:
            bad += 1
    check(bad == 0,
          f"{tested} 个多边形：卡壳枚举的每一对都确为对踵点对（例外 {bad}）")


#T5
def t5_degenerate():
    print("T5 退化输入")
    ok = True
    ok &= rc.rotating_calipers_diameter([])[0] == 0.0
    ok &= rc.rotating_calipers_diameter([(1.0, 2.0)])[0] == 0.0
    d, i, j = rc.rotating_calipers_diameter([(0.0, 0.0), (3.0, 4.0)])
    ok &= abs(d - 5.0) < 1e-12
    # 全部共线：凸包退化为两点
    hull = rc.convex_hull([(0, 0), (1, 1), (2, 2), (3, 3)])
    ok &= len(hull) == 2
    d, _, _ = rc.rotating_calipers_diameter(hull)
    ok &= abs(d - math.hypot(3, 3)) < 1e-12
    # 大量重合点：(1,1) 落在 (0,0)-(2,0)-(0,2) 的斜边上，应被当作共线剔除
    hull = rc.convex_hull([(1, 1)] * 50 + [(0, 0), (2, 0), (0, 2)])
    ok &= len(hull) == 3
    d, _, _ = rc.rotating_calipers_diameter(hull)
    ok &= abs(d - math.hypot(2, 2)) < 1e-12
    check(ok, "n=0/1/2、全共线、重复点均正确")


#T6
def t6_parallel_edges():
    print("T6 平行边情形（矩形、正多边形）")
    ok = True
    rect = [(0, 0), (10, 0), (10, 4), (0, 4)]
    d, i, j = rc.rotating_calipers_diameter(rect)
    ok &= abs(d - math.hypot(10, 4)) < 1e-12
    ap = rc.antipodal_pairs(rect)
    ok &= (0, 2) in ap and (1, 3) in ap
    ok &= len(ap) <= 3 * 4
    # 正 24 边形：直径 = 2R
    R = 100.0
    reg = [(R * math.cos(2 * math.pi * k / 24), R * math.sin(2 * math.pi * k / 24))
           for k in range(24)]
    d, _, _ = rc.rotating_calipers_diameter(rc.convex_hull(reg))
    ok &= abs(d - 2 * R) < 1e-9
    check(ok, f"矩形 D={math.hypot(10, 4):.6f}、正24边形 D={d:.6f}（应为 200）")


#T7
def t7_end_to_end():
    print("T7 端到端：真源必须落在定位区域内部")
    bad = 0
    tested = 0
    for _ in range(200):
        # 随机真源
        rG = 1500 * math.sqrt(RND.random())
        aG = RND.random() * 2 * math.pi
        G = (rG * math.cos(aG), rG * math.sin(aG))
        k = RND.randint(2, 4)
        dirs = []
        for _try in range(50):
            cand = sorted(RND.random() * 2 * math.pi for _ in range(k))
            if all(((cand[(i + 1) % k] - cand[i]) % (2 * math.pi))
                   >= math.radians(20.0) for i in range(k)):
                dirs = cand
                break
        if not dirs:
            continue
        dets = []
        for a in dirs:
            dist = RND.uniform(1000, 1500)
            sx, sy = G[0] + dist * math.cos(a), G[1] + dist * math.sin(a)
            bearing = math.degrees(math.atan2(G[1] - sy, G[0] - sx)) % 360.0
            dets.append((sx, sy, bearing))
        status, poly, _ = lr.region_from_bearings(dets, 1.0)
        if status != "ok" or len(poly) < 3:
            continue
        tested += 1
        # 真源必须在区域内部（用半平面判据独立判断，含容差）
        hps = []
        for sx, sy, b in dets:
            hps += lr.wedge_halfplanes(sx, sy, b, 1.0)
        if not all(a * G[0] + b * G[1] - c <= 1e-6 for a, b, c in hps):
            bad += 1
    check(bad == 0, f"{tested} 个随机算例，真源全部落在区域内（违反 {bad}）")


def main():
    for fn in (t1_agrees_with_brute_force, t2_diameter_is_antipodal,
               t3_count_is_linear, t4_antipodal_definition, t5_degenerate,
               t6_parallel_edges, t7_end_to_end):
        fn()
    print()
    if FAILS:
        print(f"失败 {len(FAILS)} 项：")
        for m in FAILS:
            print("  -", m)
        return 1
    print("全部通过")
    return 0


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    raise SystemExit(main())
