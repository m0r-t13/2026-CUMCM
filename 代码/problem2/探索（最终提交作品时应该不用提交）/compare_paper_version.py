"""把《问题2_求解_论文版.md》的算法原样实现，与我的版本在同一把尺子上对比。

论文版的算法（按其 §16 伪代码 / §17 Python）
-------------------------------------------
源不确定区 𝒢：ε ∈ linspace(-1°, 1°, 21)，φ = θ1 + ε，
              r ∈ [0, r_max(φ)]（r_max 由目标圆域二次方程解出），步长 20 m。
              —— 注意：**没有**近场 5 m 阈值，**没有**对源施加 1500 m 接收半径上限。
误差指标 U(S1,S2,G) = d1·d2 / |sin∠S1GS2| = (d1 d2)² / |v1×v2|   （量纲 m²）
约束：|S2 - G| ≤ R 对所有 G ∈ 𝒢（R 固定取 1000 或 1500）
目标：min_{S2} max_{G∈𝒢} U
搜索：先 50 m 全域粗搜，再在最优点周围 ±100 m、5 m 精搜。

我的版本
--------
源不确定区：𝒢 ∩ {5 < r ≤ 1500}（近场阈值 + 接收半径帽）
误差指标：定位区域直径 D（米），用问题一的精确方法算
约束：|S2 - G| ≤ max(r1, 1000)
目标：min_{S2} max D

公平对比的做法
--------------
两者各自给出 S2* 后，用**同一套尺子**评价：
  (i)  正确的源区上最坏 D（我的精确方法）
  (ii) 正确的源区上最坏完成时间无关量 U
  (iii) 各自的可行性是否在正确模型下成立
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from screen_min_error import (  # noqa: E402
    region_diameter, bearing_from, source_vertices, edge_source,
    t_arena_along, unit, ARENA, R_MIN, T_RECV, T_NEAR, ALPHA,
)

R_THEIRS = 1000.0     # 论文版的保守取值


# ==================================================== 论文版的源不确定区
def source_paper(S1, th, d_eps=0.1, d_r=20.0):
    """论文版 𝒢：ε∈[-1°,1°] 21 档，r∈[0, r_max] 步长 d_r。无近场阈值、无 1500 帽。"""
    S1 = np.asarray(S1, float)
    pts = []
    for e in np.arange(-ALPHA, ALPHA + 1e-9, d_eps):
        phi = th + e
        u = unit(phi)
        b = float(S1 @ u)
        C = ARENA * ARENA - float(S1 @ S1)
        disc = b * b + C
        if disc <= 0:
            continue
        r_max = -b + math.sqrt(disc)
        if r_max <= 0:
            continue
        for r in np.arange(0.0, r_max, d_r):
            pts.append(S1 + r * u)
    return np.array(pts)


def U_metric(V, S1, S2):
    """U = (d1 d2)² / |v1 × v2|，向量化。V 为 (N,2)。"""
    v1 = S1 - V
    v2 = S2 - V
    d1 = np.hypot(v1[:, 0], v1[:, 1])
    d2 = np.hypot(v2[:, 0], v2[:, 1])
    cross = np.abs(v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0])
    return (d1 * d2) ** 2 / np.maximum(cross, 1e-12), d1, d2, cross


def solve_paper(S1, th, G, R=R_THEIRS, coarse=50.0, fine=5.0, halfwin=100.0):
    """论文版：50 m 粗搜 + 局部 5 m 精搜。返回 (S2*, worst_U)。

    退化点处理：论文版 𝒢 含 r = 0（即 G = S1 处），此时 d1 = 0、U 发散。
    按其伪代码里 `if d1 < 1e-6 or d2 < 1e-6: continue` 的做法**跳过**这些点。
    """
    S1 = np.asarray(S1, float)
    d1_all = np.linalg.norm(S1 - G, axis=1)
    keep = d1_all > 1e-6                    # 跳过退化点
    G = G[keep]

    def scan(xs, ys):
        best_S2, best_U = None, np.inf
        for x in xs:
            for y in ys:
                S2 = np.array([float(x), float(y)])
                if np.linalg.norm(S2) > ARENA:
                    continue
                _, d1, d2, cross = U_metric(G, S1, S2)
                if np.any(d2 > R + 1e-9):
                    continue
                if np.any(cross < 1e-8):
                    continue
                U = float(np.max((d1 * d2) ** 2 / np.maximum(cross, 1e-12)))
                if U < best_U:
                    best_U, best_S2 = U, S2
        return best_S2, best_U

    xs = np.arange(-ARENA, ARENA + 1e-9, coarse)
    s0, u0 = scan(xs, xs)
    if s0 is None:
        return None, None
    xs2 = np.arange(s0[0] - halfwin, s0[0] + halfwin + 1e-9, fine)
    ys2 = np.arange(s0[1] - halfwin, s0[1] + halfwin + 1e-9, fine)
    s1, u1 = scan(xs2, ys2)
    return (s1, u1) if s1 is not None else (s0, u0)


# ==================================================== 我的版本
def solve_mine(S1, th, na=130, nb=130, n_edge=100):
    """我的版本：约束用 4 顶点（严格），误差沿楔边稠密采样。返回最优 S2 与最坏 D。"""
    S1 = np.asarray(S1, float)
    u1, n1 = unit(th), unit(th + 90.0)
    Vc = source_vertices(S1, th, R_MIN)
    Vo = edge_source(S1, th, T_RECV, n_edge)

    a = np.linspace(-ARENA, T_NEAR + R_MIN, na)
    b = np.linspace(-R_MIN, R_MIN, nb)
    A, B = np.meshgrid(a, b)
    S2x = S1[0] + A * u1[0] + B * n1[0]
    S2y = S1[1] + A * u1[1] + B * n1[1]
    ok = np.ones_like(A, dtype=bool)
    for V in Vc:
        ok &= ((S2x - V[0]) ** 2 + (S2y - V[1]) ** 2) <= R_MIN ** 2
    ok &= (S2x ** 2 + S2y ** 2) <= ARENA ** 2

    best, bestS2 = np.inf, None
    for i, j in np.argwhere(ok):
        S2 = np.array([S2x[i, j], S2y[i, j]])
        worst = 0.0
        bad = False
        for V in Vo:
            d = region_diameter(S1, bearing_from(S1, V), S2, bearing_from(S2, V))
            if d is None:
                bad = True
                break
            worst = max(worst, d)
        if bad:
            continue
        if worst < best:
            best, bestS2 = worst, S2
    return bestS2, best


def evaluate(S2, S1, G_correct, G_paper, Vo_edge):
    """用同一把尺子评价某个 S2：正确的源区上的最坏 D 与最坏 U。"""
    d_worst, bad = 0.0, False
    for V in Vo_edge:
        d = region_diameter(S1, bearing_from(S1, V), S2, bearing_from(S2, V))
        if d is None:
            bad = True
            break
        d_worst = max(d_worst, d)
    U_c, d1c, d2c, _ = U_metric(G_correct, S1, S2)
    U_p, _, d2p, _ = U_metric(G_paper, S1, S2)
    return dict(
        D=float("inf") if bad else d_worst, D_unbounded=bad,
        U_correct=float(np.max(U_c)), U_paper=float(np.max(U_p)),
        max_r2_correct=float(np.max(d2c)), max_r2_paper=float(np.max(d2p)),
    )


def main():
    CASES = [
        ("(a) 视线沿半径朝外", (0.0, 0.0), 0.0),
        ("(b) 斜穿圆域", (1200.0, 700.0), 120.0),
        ("(c) 近切向斜穿", (1500.0, 900.0), 120.0),
    ]
    print("=" * 112)
    print("论文版算法 vs 本仓库版本 —— 在同一把尺子上对比")
    print("=" * 112)

    rows = []
    for tag, S1, th in CASES:
        S1 = np.asarray(S1, float)
        Gp = source_paper(S1, th)
        Vc = source_vertices(S1, th, R_MIN)
        Vo_full = edge_source(S1, th, T_RECV, 100)
        Gc = np.vstack([Vo_full, Vc])       # 正确源区（含近场阈值与 1500 帽）

        print(f"\n{tag}   S1={S1}, θ1={th}°")
        print(f"   论文版 𝒢：{len(Gp)} 点，r 最大 "
              f"{np.max(np.linalg.norm(Gp - S1, axis=1)):.0f} m")
        print(f"   正确源区：r ∈ ({T_NEAR:.0f}, {T_RECV:.0f}] m")

        Sp, Up = solve_paper(S1, th, Gp)
        Sm, Dm = solve_mine(S1, th)
        print(f"   论文版最优 S2* = ({Sp[0]:8.1f},{Sp[1]:8.1f})   最坏 U = {Up:.3e}")
        print(f"   本版  最优 S2* = ({Sm[0]:8.1f},{Sm[1]:8.1f})   最坏 D = {Dm:.2f} m")

        ep = evaluate(Sp, S1, Gc, Gp, Vo_full)
        em = evaluate(Sm, S1, Gc, Gp, Vo_full)

        def fmt(e):
            d = "无界" if e["D_unbounded"] else f"{e['D']:.2f} m"
            return d, f"{e['U_correct']:.3e}", f"{e['U_paper']:.3e}", \
                   f"{e['max_r2_correct']:.0f}", f"{e['max_r2_paper']:.0f}"

        dp, upc, upp, r2c_p, r2p_p = fmt(ep)
        dm, umc, ump, r2c_m, r2p_m = fmt(em)
        print(f"   {'':22}{'正确源区最坏 D':>16}{'正确源区最坏 U':>16}"
              f"{'论文源区最坏 U':>16}{'到源最远(正确)':>15}")
        print(f"   {'论文版 S2*':<22}{dp:>16}{upc:>16}{upp:>16}{r2c_p:>15}")
        print(f"   {'本版 S2*':<22}{dm:>16}{umc:>16}{ump:>16}{r2c_m:>15}")

        feasible_p = ep["max_r2_correct"] <= R_MIN + 1e-6
        feasible_m = em["max_r2_correct"] <= R_MIN + 1e-6
        print(f"   在正确源区上满足 |S2-G|<=1000 ？  论文版 {'是' if feasible_p else '**否**'}"
              f"   本版 {'是' if feasible_m else '否'}")
        rows.append((tag, Sp, Sm, ep, em, Dm))

    # -------- 汇总 --------
    print("\n" + "=" * 112)
    print("汇总：用「正确源区上的最坏定位误差 D」作为统一尺子")
    print("=" * 112)
    print(f"{'场景':<22}{'论文版 D':>14}{'本版 D':>12}{'差距':>10}"
          f"{'论文版是否可行':>16}")
    for tag, Sp, Sm, ep, em, Dm in rows:
        dp = float("inf") if ep["D_unbounded"] else ep["D"]
        gap = "—" if not np.isfinite(dp) else f"{(dp - Dm) / Dm:.1%}"
        ok = ep["max_r2_correct"] <= R_MIN + 1e-6
        print(f"{tag:<22}{('无界' if not np.isfinite(dp) else f'{dp:.2f} m'):>14}"
              f"{f'{Dm:.2f} m':>12}{gap:>10}{('是' if ok else '**否**'):>16}")


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    main()
