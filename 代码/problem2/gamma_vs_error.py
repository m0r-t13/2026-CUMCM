"""最坏交会角 γ_min 越大就越好吗？—— 用真实定位误差来检验。

疑问：
前面一直把"最大化最坏交会角 γ_min"当作第二问的目标。但交会角只是
代理指标：真正的定位误差还取决于两个距离。

本脚本直接对比两个目标：
    目标 A：最大化 γ_min（最坏交会角）
    目标 B：最小化最坏定位误差（用问题一的精确方法算定位区域直径 D）

做法
----
对场景 A(S1=(0,0), θ1=0°)在 (a,b) 网格上：
  1. 约束：对源可行域每个顶点，|S2-V| ≤ R_min=1000；且 |S2| ≤ 1800
  2. 目标 A:γ_min(S2) = 90° − max_V |γ(V;S2) − 90°|
  3. 目标 B:D_worst(S2) = max_V [ 源在 V 时两站交会定位区域的直径 ]
     定位区域用问题一的半平面求交 + 顶点对枚举精确计算
  4. 分别取最优，比较两个最优点是否一致
"""

from __future__ import annotations

import math
import os
import sys
from itertools import combinations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
P1 = os.path.join(os.path.dirname(HERE), "problem1(part1)")
if os.path.isdir(P1):
    sys.path.insert(0, P1)

ARENA = 1800.0
R_MIN = 1000.0
T_RECV = 1500.0
T_NEAR = 5.0
ALPHA = 1.0


# ---------------------------------------------------------------- 几何工具
def unit(deg):
    r = math.radians(deg)
    return np.array([math.cos(r), math.sin(r)])


def wedge_hps(sx, sy, th_deg, alpha=ALPHA):
    th, al = math.radians(th_deg), math.radians(alpha)
    ux, uy = math.cos(th), math.sin(th)
    t = math.tan(al)
    cuS = ux * sy - uy * sx
    duS = ux * sx + uy * sy
    return [(-uy - t * ux, ux - t * uy, cuS - t * duS),
            (uy - t * ux, -ux - t * uy, -cuS - t * duS),
            (-ux, -uy, -duS)]


def _hull(pts):
    pts = sorted({(float(x), float(y)) for x, y in pts})
    if len(pts) <= 2:
        return np.array(pts, float).reshape(-1, 2)

    def half(seq):
        out = []
        for p in seq:
            while len(out) >= 2:
                (ox, oy), (ax, ay) = out[-2], out[-1]
                if (ax - ox) * (p[1] - oy) - (ay - oy) * (p[0] - ox) <= 0:
                    out.pop()
                else:
                    break
            out.append(p)
        return out

    lo, up = half(pts), half(pts[::-1])
    return np.array(lo[:-1] + up[:-1], float).reshape(-1, 2)


def region_diameter(S1, th1, S2, th2, scale=3000.0):
    """两站交会定位区域的直径（问题一方法）。无界或空返回 None。"""
    hps = wedge_hps(S1[0], S1[1], th1) + wedge_hps(S2[0], S2[1], th2)
    angs = sorted(math.atan2(b, a) for a, b, _ in hps)
    gaps = [angs[i + 1] - angs[i] for i in range(len(angs) - 1)]
    gaps.append(angs[0] + 2 * math.pi - angs[-1])
    if max(gaps) >= math.pi - 1e-12:
        return None                       # 无界
    tol = 1e-8 * scale
    cand = []
    for l1, l2 in combinations(hps, 2):
        a1, b1, c1 = l1
        a2, b2, c2 = l2
        det = a1 * b2 - a2 * b1
        if abs(det) < 1e-14:
            continue
        x = (c1 * b2 - c2 * b1) / det
        y = (a1 * c2 - a2 * c1) / det
        if all(a * x + b * y - c <= tol for a, b, c in hps):
            cand.append((x, y))
    if not cand:
        return 0.0
    poly = _hull(cand)
    if len(poly) < 2:
        return 0.0
    best = 0.0
    for i in range(len(poly)):
        for j in range(i + 1, len(poly)):
            best = max(best, float(np.linalg.norm(poly[i] - poly[j])))
    return best


def t_arena_along(S1, phi_deg, arena=ARENA):
    u = unit(phi_deg)
    b0 = float(S1 @ u)
    disc = b0 * b0 + arena * arena - float(S1 @ S1)
    return None if disc <= 0 else -b0 + math.sqrt(disc)


def source_vertices(S1, th, r_cap=T_RECV, n_edge=0):
    """源可行域 W1 ∩ 圆域 ∩ {r ≤ r_cap} 的采样点。

    n_edge = 0 时只返回 4 个顶点（对探测约束充分，但对误差函数不充分）；
    n_edge > 0 时沿两条楔边各采 n_edge 个点（推荐，用于误差评估）。
    """
    out = []
    for s in (+1.0, -1.0):
        phi = th + s * ALPHA
        u = unit(phi)
        t_ar = t_arena_along(S1, phi)
        r_far = max(0.0, min(r_cap, t_ar if t_ar else 0.0))
        if r_far <= 1e-9:
            continue
        if n_edge <= 0:
            out.append(S1 + T_NEAR * u)
            out.append(S1 + r_far * u)
        else:
            for r in np.linspace(T_NEAR, r_far, n_edge):
                out.append(S1 + r * u)
    return np.array(out)


def bearing_from(S, G):
    return math.degrees(math.atan2(G[1] - S[1], G[0] - S[0])) % 360.0


def cross_angle(V, S1, S2):
    v1, v2 = S1 - V, S2 - V
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 < 1e-9 or n2 < 1e-9:
        return 90.0
    c = float(np.clip(v1 @ v2 / (n1 * n2), -1.0, 1.0))
    return math.degrees(math.acos(c))


# ---------------------------------------------------------------- 主实验
def main():
    S1 = np.array([0.0, 0.0])
    TH1 = 0.0
    u1, n1 = unit(TH1), unit(TH1 + 90.0)

    Vc = source_vertices(S1, TH1, R_MIN)              # 约束顶点（帽 1000）
    Vo = source_vertices(S1, TH1, T_RECV, n_edge=60)  # 误差评估：沿楔边采样
    print(f"源可行域（约束用，R_cap=1000）：{len(Vc)} 个顶点")
    print(f"源可行域（误差评估用，R_cap=1500）：沿两条楔边采样 {len(Vo)} 点，"
          f"r ∈ [{np.linalg.norm(Vo[0]-S1):.0f}, "
          f"{np.linalg.norm(Vo[-1]-S1):.0f}] m")

    # 网格
    aa = np.linspace(-200.0, T_NEAR + R_MIN, 90)
    bb = np.linspace(20.0, R_MIN - 1e-6, 90)
    print(f"\n扫描网格：a {len(aa)} × b {len(bb)} = {len(aa)*len(bb)} 个候选点")

    rows = []
    for a in aa:
        for b in bb:
            S2 = S1 + a * u1 + b * n1
            if np.linalg.norm(S2) > ARENA:
                continue
            if any(np.linalg.norm(S2 - V) > R_MIN + 1e-9 for V in Vc):
                continue
            # 目标 A：最坏交会角
            gmin = 90.0 - max(abs(cross_angle(V, S1, S2) - 90.0) for V in Vo)
            # 目标 B：最坏定位区域直径
            dmax = 0.0
            ok = True
            for V in Vo:
                th1v = bearing_from(S1, V)      # 注意：用真实示向度，不是 0°
                th2 = bearing_from(S2, V)
                d = region_diameter(S1, th1v, S2, th2)
                if d is None:
                    ok = False
                    break
                dmax = max(dmax, d)
            if not ok:
                continue
            rows.append((a, b, gmin, dmax))

    arr = np.array(rows)
    print(f"可行候选点 {len(arr)} 个")

    iA = int(np.argmax(arr[:, 2]))
    iB = int(np.argmin(arr[:, 3]))
    aA, bA, gA, dA = arr[iA]
    aB, bB, gB, dB = arr[iB]

    print("\n" + "=" * 88)
    print("两个目标各自的最优解")
    print("=" * 88)
    print(f"{'目标':<26}{'a':>8}{'b':>9}{'γ_min':>10}{'最坏定位直径 D':>16}")
    print("-" * 88)
    print(f"{'A: 最大化 γ_min':<26}{aA:>8.1f}{bA:>9.1f}{gA:>9.2f}°{dA:>15.2f} m")
    print(f"{'B: 最小化最坏 D':<26}{aB:>8.1f}{bB:>9.1f}{gB:>9.2f}°{dB:>15.2f} m")
    print()
    print(f"  两点相距 = {math.hypot(aA-aB, bA-bB):.1f} m")
    print(f"  用 A 的做法，定位直径比最优差 {(dA-dB)/dB:.2%}")
    print(f"  用 B 的做法，交会角比最优差 {gA-gB:.2f}°")

    # γ_min 与 D 的相关性
    g = arr[:, 2]
    d = arr[:, 3]
    corr = float(np.corrcoef(g, d)[0, 1])
    print(f"\n  γ_min 与最坏定位直径 D 的相关系数 = {corr:.3f}")
    print(f"  （若接近 −1 说明两者完全一致；越接近 0 越说明 γ_min 不可靠）")

    # 按 γ_min 分箱看 D 的分布
    print("\n  按 γ_min 分箱，看对应的最坏定位直径 D：")
    print(f"  {'γ_min 区间':<16}{'样本数':>8}{'D 中位数':>12}{'D 最小':>10}{'D 最大':>10}")
    edges = [0, 20, 30, 35, 38, 40, 41, 42]
    for i in range(len(edges) - 1):
        m = (g >= edges[i]) & (g < edges[i + 1])
        if m.sum() == 0:
            continue
        print(f"  {f'{edges[i]}–{edges[i+1]}°':<16}{int(m.sum()):>8}"
              f"{np.median(d[m]):>11.2f} m{d[m].min():>9.2f} m{d[m].max():>9.2f} m")

    out = os.path.join(HERE, "交会角与定位误差关系.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("a,b,gamma_min,D_worst_m\n")
        for a, b, gg, dd in arr:
            fh.write(f"{a:.4f},{b:.4f},{gg:.6f},{dd:.6f}\n")
    print("\n-> 逐点数据:", out)

    # ---- 诊断：D 最大的那个点，看是哪个源顶点造成的 ----
    iMax = int(np.argmax(arr[:, 3]))
    aM, bM, gM, dM = arr[iMax]
    S2M = S1 + aM * u1 + bM * n1
    print("\n" + "=" * 88)
    print(f"诊断：最坏定位直径最大的候选点  a={aM:.1f}, b={bM:.1f}, "
          f"γ_min={gM:.2f}°, D={dM:.1f} m")
    print("=" * 88)
    best = None
    for i, V in enumerate(Vo):
        th1v = bearing_from(S1, V)
        th2 = bearing_from(S2M, V)
        d = region_diameter(S1, th1v, S2M, th2)
        if d is None:
            continue
        if best is None or d > best[0]:
            best = (d, V, cross_angle(V, S1, S2M))
    if best:
        d, V, gam = best
        print(f"  决定最坏误差的源位置：({V[0]:.2f}, {V[1]:.2f})，"
              f"r1={np.linalg.norm(V - S1):.1f} m，"
              f"该处交会角={gam:.3f}°，定位直径={d:.1f} m")
        print(f"  说明：交会角 $\\gamma<2\\alpha=2°$ 时两方向锥重叠、区域无界；"
              f"此处 {gam:.2f}° 仅略大于 2°，故区域被拉得极长。")


    # ---- 图：γ_min 与真实定位误差的关系 ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for fam in ("Microsoft YaHei", "SimHei", "SimSun", "DejaVu Sans"):
        try:
            font_manager.findfont(fam, fallback_to_default=False)
            plt.rcParams["font.sans-serif"] = [fam]
            break
        except Exception:
            continue
    plt.rcParams["axes.unicode_minus"] = False

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16.5, 6.6))
    gg, dd = arr[:, 2], arr[:, 3]
    ax1.semilogy(gg, dd, ".", color="tab:blue", ms=4, alpha=0.35,
                 label=f"每个可行候选点（{len(arr)} 个）")
    centers, meds = [], []
    for i in range(len(edges) - 1):
        m = (gg >= edges[i]) & (gg < edges[i + 1])
        if m.sum() == 0:
            continue
        centers.append(0.5 * (edges[i] + edges[i + 1]))
        meds.append(float(np.median(dd[m])))
    ax1.semilogy(centers, meds, "o-", color="tab:red", lw=2.4, ms=9,
                 label="各区间中位数")
    ax1.semilogy(gA, dA, "*", color="tab:green", ms=24, mec="k", mew=1.2,
                 label=f"A 最大化 γ：{gA:.1f}°, {dA:.0f} m")
    ax1.semilogy(gB, dB, "P", color="tab:orange", ms=15, mec="k", mew=1.2,
                 label=f"B 最小化 D：{gB:.1f}°, {dB:.0f} m")
    ax1.axvline(2 * ALPHA, color="k", ls="--", lw=1.6)
    ax1.annotate(f"$\\gamma=2\\alpha$={2*ALPHA:.0f}° 以下\n两个方向锥重叠\n区域无界",
                 (2 * ALPHA, max(dd) * 0.25), textcoords="offset points",
                 xytext=(10, 0), fontsize=10)
    ax1.set_xlabel("最坏交会角 $\\gamma_{min}$ / 度")
    ax1.set_ylabel("最坏定位区域直径 $D$ / m（对数轴）")
    ax1.grid(alpha=0.25, which="both")
    ax1.legend(fontsize=9.5, loc="upper right")
    ax1.set_title(f"(a) $\\gamma_{{min}}$ 与真实定位误差："
                  f"中位数单调改善，但同区间内散布极大（相关系数 {corr:.3f}）",
                  fontsize=12)

    ax2.plot(centers, meds, "o-", color="tab:red", lw=2.4, ms=10)
    for cx, my in zip(centers, meds):
        ax2.annotate(f"{my:.0f}", (cx, my), textcoords="offset points",
                     xytext=(0, 9), ha="center", fontsize=10)
    ax2.axhline(dB, color="tab:orange", ls="--", lw=1.6)
    ax2.annotate(f"真正的最优 D = {dB:.0f} m", (centers[0], dB),
                 textcoords="offset points", xytext=(6, 6),
                 color="tab:orange", fontsize=11)
    ax2.set_yscale("log")
    ax2.set_xlabel("最坏交会角 $\\gamma_{min}$ / 度")
    ax2.set_ylabel("该区间内最坏定位直径的中位数 / m（对数轴）")
    ax2.grid(alpha=0.25, which="both")
    ax2.set_title("(b) 作为粗筛是有效的：中位数从 475 m 降到 131 m，\n"
                  "但降幅在 $\\gamma_{min}>35°$ 后明显变平", fontsize=12)

    fig.suptitle("「最坏交会角越大越好」只对了一半：它是有效粗筛，"
                 "但不是精确目标——两个目标的最优点相距 "
                 f"{math.hypot(aA-aB, bA-bB):.0f} m", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    p = os.path.join(HERE, "交会角与定位误差关系.png")
    fig.savefig(p, dpi=145)
    plt.close(fig)
    print("->", p)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    main()
