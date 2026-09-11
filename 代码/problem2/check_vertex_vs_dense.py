"""核查：把源不确定区只取 4 个顶点，会不会漏掉真正的最坏情况？

我的版本对"定位误差"取 max 时只用了源可行域的 4 个顶点，依据是
"凸函数最大值在极点取到"。但那个依据只对**探测约束**（|S2-G| 关于 G 凸）严格成立；
对**误差函数 D(G)** 不一定成立（D 不是 G 的凸函数）。

所以必须实测：把源区稠密采样，看最大值是否比 4 顶点更大。
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from screen_min_error import (  # noqa: E402
    region_diameter, bearing_from, source_vertices, t_arena_along, unit,
    ARENA, R_MIN, T_RECV, T_NEAR, ALPHA,
)


def ds(S1, th, r_cap, n_eps=21, n_r=60):
    """把源不确定区稠密采样（含角楔厚度）。"""
    S1 = np.asarray(S1, float)
    pts = []
    for e in np.linspace(-ALPHA, ALPHA, n_eps):
        phi = th + e
        u = unit(phi)
        t_ar = t_arena_along(S1, phi)
        r_hi = max(0.0, min(r_cap, t_ar if t_ar else 0.0))
        if r_hi <= T_NEAR:
            continue
        for r in np.linspace(T_NEAR, r_hi, n_r):
            pts.append(S1 + r * u)
    return np.array(pts)


def worst_D(S1, th1, S2, Gs):
    worst, arg = 0.0, None
    for G in Gs:
        d = region_diameter(S1, bearing_from(S1, G), S2, bearing_from(S2, G))
        if d is None:
            return None, G
        if d > worst:
            worst, arg = d, G
    return worst, arg


def edge_source(S1, th, r_cap, n_edge=80):
    """只沿两条楔边采样（含近/远端角点）——比全稠密便宜得多。"""
    S1 = np.asarray(S1, float)
    pts = []
    for s in (+1.0, -1.0):
        phi = th + s * ALPHA
        u = unit(phi)
        t_ar = t_arena_along(S1, phi)
        r_hi = max(0.0, min(r_cap, t_ar if t_ar else 0.0))
        if r_hi <= T_NEAR:
            continue
        for r in np.linspace(T_NEAR, r_hi, n_edge):
            pts.append(S1 + r * u)
    return np.array(pts)


def main():
    CASES = [
        ("(a) 沿半径朝外", (0.0, 0.0), 0.0),
        ("(b) 斜穿圆域", (1200.0, 700.0), 120.0),
        ("(c) 近切向斜穿", (1500.0, 900.0), 120.0),
    ]
    # 每个场景取几个有代表性的候选 S2（含最优附近与两侧）
    PROBES = [(0.6, 0.75), (0.9, 0.4), (0.4, 0.9), (0.8, 0.2), (0.3, 0.35)]

    print("=" * 108)
    print("核查：三种源区采样方式给出的最坏 D —— 4顶点 / 只沿两条楔边 / 全稠密")
    print("=" * 108)
    bad4 = bad_edge = 0
    gap4 = gap_edge = 0.0
    for tag, S1, th in CASES:
        S1 = np.asarray(S1, float)
        u1, n1 = unit(th), unit(th + 90.0)
        V4 = source_vertices(S1, th, T_RECV)
        VE = edge_source(S1, th, T_RECV, 80)
        Dd = ds(S1, th, T_RECV)
        Gd = np.vstack([Dd, V4])
        print(f"\n{tag}   4顶点={len(V4)}  楔边采样={len(VE)}  全稠密={len(Gd)}")
        print(f"   {'探针 (a,b)':<16}{'D(4顶点)':>10}{'D(楔边)':>10}{'D(全稠密)':>11}"
              f"{'4顶点漏?':>10}{'楔边漏?':>9}")
        for fa, fb in PROBES:
            a = fa * (T_NEAR + R_MIN)
            b = fb * R_MIN
            S2 = S1 + a * u1 + b * n1
            d4, _ = worst_D(S1, th, S2, V4)
            de, _ = worst_D(S1, th, S2, VE)
            dd, _ = worst_D(S1, th, S2, Gd)
            if None in (d4, de, dd):
                continue
            g4 = (dd - d4) / max(d4, 1e-9)
            ge = (dd - de) / max(dd, 1e-9)
            gap4 = max(gap4, g4)
            gap_edge = max(gap_edge, ge)
            if dd > d4 * 1.001:
                bad4 += 1
            if dd > de * 1.001:
                bad_edge += 1
            print(f"   ({a:6.0f},{b:6.0f}){d4:>10.2f}{de:>10.2f}{dd:>11.2f}"
                  f"{('**是**' if dd > d4*1.001 else '否'):>10}"
                  f"{('**是**' if dd > de*1.001 else '否'):>9}")

    print("\n" + "=" * 108)
    print(f"4 顶点采样：漏掉最坏情况的探针 {bad4} 个，最大相对差 {gap4:.2%}")
    print(f"楔边采样  ：漏掉最坏情况的探针 {bad_edge} 个，最大相对差 {gap_edge:.2%}")
    print("=" * 108)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    main()
