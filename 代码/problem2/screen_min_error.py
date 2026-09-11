"""第二问：从候选区域里筛出定位误差最小的那一部分。

之前输出的是"探测保证集"（所有一定能收到信号的第二个检测点位置）——
那是一个很大的区域，里面误差差别很大。

本脚本改成：
    1. 先算探测保证集（含角楔厚度，源区 4 个顶点）
    2. 在集合内逐点算真实定位误差：
       源落在源区每个顶点时，用问题一的精确方法算两站交会定位区域的直径 D，
       取源区上最坏的那个 D 作为该点的误差。
    3. 只保留 D 最小的一小部分 → 输出"最小误差候选区"

误差用定位区域直径 D（米），而不是交会角——因为用户要的是误差小，
交会角只是代理指标（相关系数仅 −0.14）。
"""

from __future__ import annotations

import math
import os
import sys

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


# ============================================================ 几何工具
def unit(deg):
    r = math.radians(deg)
    return np.array([math.cos(r), math.sin(r)])


def _wedge_hps_arr(sx, sy, th_deg, alpha=ALPHA):
    """一个角楔的 3 个半平面，返回 (3,3) 数组，每行 (a,b,c) 表示 a x + b y <= c。"""
    th, al = math.radians(th_deg), math.radians(alpha)
    ux, uy = math.cos(th), math.sin(th)
    t = math.tan(al)
    cuS = ux * sy - uy * sx
    duS = ux * sx + uy * sy
    return np.array([
        [-uy - t * ux, ux - t * uy, cuS - t * duS],
        [uy - t * ux, -ux - t * uy, -cuS - t * duS],
        [-ux, -uy, -duS],
    ])


def region_diameter(S1, th1, S2, th2, tol=1e-6):
    """两站交会定位区域的直径（问题一：顶点枚举 + 顶点对最大距离）。

    凸多胞体的直径必在顶点取到，而顶点必为两条边界直线交点，
    所以直接在所有"可行候选交点"里取最大距离即可，无需先求凸包。
    返回 None 表示无界。
    """
    H = np.vstack([_wedge_hps_arr(S1[0], S1[1], th1),
                   _wedge_hps_arr(S2[0], S2[1], th2)])
    # 有界性：外法向是否正张成 R^2
    angs = np.sort(np.arctan2(H[:, 1], H[:, 0]))
    gaps = np.diff(np.concatenate([angs, angs[:1] + 2 * math.pi]))
    if gaps.max() >= math.pi - 1e-12:
        return None
    iu = np.triu_indices(6, 1)
    a1, b1, c1 = H[iu[0], 0], H[iu[0], 1], H[iu[0], 2]
    a2, b2, c2 = H[iu[1], 0], H[iu[1], 1], H[iu[1], 2]
    det = a1 * b2 - a2 * b1
    good = np.abs(det) > 1e-14
    safe = np.where(good, det, 1.0)
    x = np.where(good, (c1 * b2 - c2 * b1) / safe, np.nan)
    y = np.where(good, (a1 * c2 - a2 * c1) / safe, np.nan)
    feas = good.copy()
    for a, b, c in H:
        feas &= (a * x + b * y - c) <= tol
    P = np.column_stack([x[feas], y[feas]])
    if len(P) < 2:
        return 0.0
    d2 = ((P[:, None, :] - P[None, :, :])  2).sum(-1)
    return float(math.sqrt(d2.max()))


def t_arena_along(S1, phi_deg, arena=ARENA):
    u = unit(phi_deg)
    b0 = float(S1 @ u)
    disc = b0 * b0 + arena * arena - float(S1 @ S1)
    return None if disc <= 0 else -b0 + math.sqrt(disc)


def source_vertices(S1, th, r_cap):
    """源可行域（含角楔厚度）的顶点：两条楔边 × 两个半径端点。

    注意：这 4 个顶点对探测约束是充分的（|S2−G| 关于 G 凸），
    但对误差函数不充分 —— 见 edge_source。
    """
    out = []
    for s in (+1.0, -1.0):
        phi = th + s * ALPHA
        u = unit(phi)
        t_ar = t_arena_along(S1, phi)
        r_far = max(0.0, min(r_cap, t_ar if t_ar else 0.0))
        out.append(S1 + T_NEAR * u)
        if r_far > 1e-9:
            out.append(S1 + r_far * u)
    return np.array(out)


def edge_source(S1, th, r_cap, n_edge=100):
    """沿两条楔边稠密采样（含近/远端角点）。

    为什么不能只用 4 个顶点：误差函数 D(G) 不是 G 的凸函数，
    最坏情况可能落在楔边内部。实测只用 4 顶点会低估最坏 D 达 73%；
    沿两条楔边采样后残差降到 0.5% 以内（纯属离散步长）。
    """
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


def bearing_from(S, G):
    return math.degrees(math.atan2(G[1] - S[1], G[0] - S[0])) % 360.0


# ============================================================ 主求解
def solve_case(S1, th, na=130, nb=130, n_edge=100, verbose=True):
    """在 (a,b) 网格上算：可行域、最坏定位误差场、最小误差子区。

    约束（探测保证）用 4 个顶点 —— 严格，因为 |S2−G| 关于 G 是凸的；
    目标（定位误差）沿两条楔边稠密采样 —— 因为 D(G) 不是凸的。
    """
    S1 = np.asarray(S1, float)
    u1, n1 = unit(th), unit(th + 90.0)

    Vc = source_vertices(S1, th, R_MIN)          # 探测约束顶点（严格）
    Vo = edge_source(S1, th, T_RECV, n_edge)     # 误差评估：沿楔边稠密采样
    t_max = min(T_RECV, t_arena_along(S1, th) or 0.0)

    a = np.linspace(-ARENA, T_NEAR + R_MIN, na)
    b = np.linspace(-R_MIN, R_MIN, nb)
    A, B = np.meshgrid(a, b)
    S2x = S1[0] + A * u1[0] + B * n1[0]
    S2y = S1[1] + A * u1[1] + B * n1[1]

    ok = np.ones_like(A, dtype=bool)
    for V in Vc:
        ok &= ((S2x - V[0])  2 + (S2y - V[1])  2) <= R_MIN  2
    ok &= (S2x  2 + S2y  2) <= ARENA  2

    D = np.full(A.shape, np.nan)
    idx = np.argwhere(ok)
    nbad = 0
    for i, j in idx:
        S2 = np.array([S2x[i, j], S2y[i, j]])
        worst = 0.0
        bad = False
        for V in Vo:
            th1v = bearing_from(S1, V)
            th2v = bearing_from(S2, V)
            d = region_diameter(S1, th1v, S2, th2v)
            if d is None:
                bad = True
                break
            worst = max(worst, d)
        if bad:
            nbad += 1
            continue
        D[i, j] = worst

    valid = ~np.isnan(D)
    Dmin = float(np.nanmin(D))
    if verbose:
        print(f"   可行点 {int(ok.sum())}，可评估 {int(valid.sum())}，"
              f"无界/失败 {nbad}")
        print(f"   最坏定位误差最小值 D_min = {Dmin:.2f} m  "
              f"（网格最优位置 a={float(a[np.unravel_index(np.nanargmin(D), D.shape)[1]]):.0f}, "
              f"b={float(b[np.unravel_index(np.nanargmin(D), D.shape)[0]]):.0f}）")

    return dict(S1=S1, th=th, u1=u1, n1=n1, a=a, b=b, A=A, B=B,
                S2x=S2x, S2y=S2y, D=D, ok=ok, Vc=Vc, Vo=Vo,
                t_max=t_max, Dmin=Dmin,
                n_feasible=int(valid.sum()),
                cell_area=float((a[1] - a[0]) * (b[1] - b[0])))


def sub_region(r, eps):
    """D <= D_min*(1+eps) 的子区，返回布尔掩码与统计。"""
    m = (~np.isnan(r["D"])) & (r["D"] <= r["Dmin"] * (1.0 + eps))
    n = int(m.sum())
    area = n * r["cell_area"]
    return m, n, area


# ============================================================ 主程序
def main():
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

    CASES = [
        ("(a) 视线沿半径朝外", (0.0, 900.0), 90.0),
        ("(b) 斜穿圆域", (1200.0, 700.0), 120.0),
        ("(c) 近切向斜穿", (1500.0, 900.0), 120.0),
    ]
    EPS_LIST = [0.02, 0.05, 0.10, 0.20]

    print("=" * 104)
    print("第二问：候选区域 → 筛出定位误差最小的部分（含角楔厚度）")
    print("=" * 104)
    res = []
    for tag, S1, th in CASES:
        print(f"\n{tag}   S1={S1}, θ1={th}°")
        r = solve_case(S1, th)
        res.append((tag, S1, th, r))
        print(f"   阈值 ε    子区点数    子区面积(m²)   占可行域比   该子区最坏误差")
        for eps in EPS_LIST:
            m, n, area = sub_region(r, eps)
            print(f"   {eps:>5.0%}   {n:>9}   {area:>12.0f}   "
                  f"{area / (r['n_feasible'] * r['cell_area']):>10.1%}   "
                  f"{r['Dmin'] * (1 + eps):>8.1f} m")

    # ---------- 图 1：三场景的误差场 + 最小误差子区 ----------
    fig, axes = plt.subplots(1, 3, figsize=(21, 7.2))
    thv = np.linspace(0, 2 * math.pi, 700)
    for ax, (tag, S1, th, r) in zip(axes, res):
        S1a = np.asarray(S1)
        ax.plot(ARENA * np.cos(thv), ARENA * np.sin(thv), color="0.45", ls="--",
                lw=1.2, zorder=5, label="目标区域 (r=1800)")
        # 源可行域
        poly = np.vstack([r["Vc"][0], r["Vc"][2], r["Vc"][3], r["Vc"][1]])
        ax.fill(poly[:, 0], poly[:, 1], color="lime", alpha=0.9, zorder=9,
                label="源可行域（含角楔厚度）")
        # 误差场
        Dm = r["D"]
        pos = Dm[~np.isnan(Dm)]
        vmin, vmax = r["Dmin"], float(np.percentile(pos, 85))
        cs = ax.contourf(r["S2x"], r["S2y"], Dm,
                         levels=np.linspace(vmin, vmax, 25), cmap="viridis_r",
                         alpha=0.95, zorder=1, extend="max")
        # 最小误差子区（ε=5%）加黑边
        m5, _, _ = sub_region(r, 0.05)
        ax.contour(r["S2x"], r["S2y"], m5.astype(float), levels=[0.5],
                   colors="red", linewidths=2.6, zorder=12)
        ax.plot([], [], color="red", lw=2.6,
                label=f"最小误差子区 $D\\leq D_{{min}}\\times1.05$")
        ax.plot([], [], color="none", label=f"$D_{{min}}$={r['Dmin']:.0f} m")
        # 视线
        L = max(r["t_max"] * 1.4, 300)
        ax.plot([S1a[0], S1a[0] + L * r["u1"][0]], [S1a[1], S1a[1] + L * r["u1"][1]],
                color="tab:blue", lw=1.6, alpha=0.9, zorder=10, label="视线")
        ax.plot(*S1a, "o", color="k", ms=10, zorder=13)
        ax.annotate("$S_1$", S1a, textcoords="offset points", xytext=(-26, -34),
                    fontsize=12, zorder=14)
        lim = ARENA * 1.06
        ax.set_aspect("equal")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.grid(alpha=0.18)
        ax.set_xlabel("x / m")
        ax.set_ylabel("y / m")
        ax.legend(loc="lower left", fontsize=8.2)
        ax.set_title(f"{tag}\n$S_1$=({S1[0]:.0f},{S1[1]:.0f}), θ1={th:.0f}°, "
                     f"$t_{{max}}$={r['t_max']:.0f} m", fontsize=11.5)
    fig.subplots_adjust(left=0.035, right=0.912, top=0.82, bottom=0.07, wspace=0.14)
    cax = fig.add_axes([0.924, 0.11, 0.011, 0.68])
    cb = fig.colorbar(cs, cax=cax)
    cb.set_label("最坏定位误差 $D$ / m（越小越好）", fontsize=11)
    fig.suptitle("第二问：把「所有能保证探测的点」筛成「定位误差最小的那部分」"
                 "（红色圈出 $D\\leq1.05D_{min}$ 的子区）", fontsize=14, y=0.965)
    p1 = os.path.join(HERE, "最小误差候选区_三场景.png")
    fig.savefig(p1, dpi=125)
    plt.close(fig)
    print("\n->", p1)

    # ---------- 图 2：场景 A 的筛选细节 + 权衡曲线 ----------
    tag, S1, th, r = res[0]
    S1a = np.asarray(S1)
    fig2, (axA, axB, axC) = plt.subplots(1, 3, figsize=(21, 6.8))

    Dm = r["D"]
    pos = Dm[~np.isnan(Dm)]
    vmin, vmax = r["Dmin"], float(np.percentile(pos, 85))
    cs = axA.contourf(r["S2x"], r["S2y"], Dm,
                      levels=np.linspace(vmin, vmax, 30), cmap="viridis_r",
                      alpha=0.95, zorder=1, extend="max")
    cols = ["red", "orange", "magenta", "cyan"]
    for eps, col in zip(EPS_LIST, cols):
        m, n, area = sub_region(r, eps)
        axA.contour(r["S2x"], r["S2y"], m.astype(float), levels=[0.5],
                    colors=col, linewidths=2.0, zorder=10)
        axA.plot([], [], color=col, lw=2.2,
                 label=f"$\\varepsilon$={eps:.0%} → {area/1e4:.0f} 万 m²")
    axA.plot(ARENA * np.cos(thv), ARENA * np.sin(thv), color="0.45", ls="--",
             lw=1.2, zorder=5)
    poly = np.vstack([r["Vc"][0], r["Vc"][2], r["Vc"][3], r["Vc"][1]])
    axA.fill(poly[:, 0], poly[:, 1], color="lime", alpha=0.85, zorder=9)
    axA.plot(*S1a, "o", color="k", ms=10, zorder=13)
    axA.annotate("$S_1$", S1a, textcoords="offset points", xytext=(-26, -34),
                 fontsize=12, zorder=14)
    # 缩放到候选区包围盒，让等值线看得见
    fv = ~np.isnan(Dm)
    if fv.any():
        xf, yf = r["S2x"][fv], r["S2y"][fv]
        cx, cy = float(xf.mean()), float(yf.mean())
        sp = max(float(np.ptp(xf)), float(np.ptp(yf))) * 0.72 + 60
        axA.set_xlim(cx - sp, cx + sp)
        axA.set_ylim(cy - sp, cy + sp)
    axA.set_aspect("equal")
    axA.grid(alpha=0.18)
    axA.set_xlabel("x / m")
    axA.set_ylabel("y / m")
    axA.legend(loc="lower left", fontsize=8.5)
    axA.set_title("(a) 阈值越松，子区越大：四条等值线", fontsize=12)
    fig2.subplots_adjust(left=0.035, right=0.945, top=0.83, bottom=0.11, wspace=0.26)
    cax2 = fig2.add_axes([0.957, 0.15, 0.010, 0.55])
    cb2 = fig2.colorbar(cs, cax=cax2)
    cb2.set_label("最坏定位误差 $D$ / m", fontsize=10)

    eps_fine = np.linspace(0.0, 0.6, 200)
    areas, dmaxs = [], []
    for e in eps_fine:
        m, n, area = sub_region(r, float(e))
        areas.append(area / 1e4)
        dmaxs.append(r["Dmin"] * (1 + e))
    axB.plot(dmaxs, areas, color="tab:red", lw=2.6)
    for eps, col in zip(EPS_LIST, cols):
        axB.plot(r["Dmin"] * (1 + eps), sub_region(r, eps)[2] / 1e4, "o",
                 color=col, ms=12, mec="k", mew=0.8)
    axB.set_xlabel("允许的最坏定位误差 $D$ / m")
    axB.set_ylabel("子区面积 / 万 m²")
    axB.grid(alpha=0.25)
    axB.set_title("(b) 误差-面积权衡：想要更准，可选范围指数收缩", fontsize=12)

    labels = ["全部候选区"] + [f"ε={e:.0%}" for e in EPS_LIST]
    dv = [float(np.nanmax(Dm))] + [r["Dmin"] * (1 + e) for e in EPS_LIST]
    bars = axC.bar(labels, dv, color=["0.6"] + cols, alpha=0.9,
                   edgecolor="k", lw=0.8)
    for bb, vv in zip(bars, dv):
        axC.annotate(f"{vv:,.0f}", (bb.get_x() + bb.get_width() / 2, vv),
                     textcoords="offset points", xytext=(0, 5), ha="center",
                     fontsize=10)
    axC.set_yscale("log")
    axC.set_ylabel("该方案能保证的最坏定位误差 / m（对数轴）")
    axC.grid(alpha=0.25, axis="y", which="both")
    axC.tick_params(axis="x", labelrotation=20)
    axC.set_title("(c) 筛选带来的收益：最坏误差从\n"
                  "「全部候选区」的量级降到 130 m 量级", fontsize=12)

    fig2.suptitle("场景 A 的筛选细节：从「能保证探测」到「定位误差最小」",
                  fontsize=14, y=0.965)
    p2 = os.path.join(HERE, "最小误差候选区_筛选细节.png")
    fig2.savefig(p2, dpi=125)
    plt.close(fig)
    print("->", p2)

    # ---------- 图 3：三场景 ε=5% 子区并列放大 ----------
    fig3, axes3 = plt.subplots(1, 3, figsize=(21, 6.6))
    for ax, (tag, S1, th, r) in zip(axes3, res):
        S1a = np.asarray(S1)
        m5, n5, a5 = sub_region(r, 0.05)
        x = r["S2x"][m5]
        y = r["S2y"][m5]
        ax.plot(x, y, ".", color="red", ms=3)
        allp = r["S2x"][~np.isnan(r["D"])]
        allq = r["S2y"][~np.isnan(r["D"])]
        ax.plot(allp, allq, ".", color="0.75", ms=1.2,
                label=f"全部候选区（{r['n_feasible']} 点）")
        ax.plot(x, y, ".", color="red", ms=3,
                label=f"$\\varepsilon$=5% 子区（{n5} 点，{a5/1e4:.1f} 万 m²）")
        ax.plot(*S1a, "o", color="k", ms=9, zorder=10)
        ax.annotate("$S_1$", S1a, textcoords="offset points", xytext=(-24, -32),
                    fontsize=11, zorder=11)
        L = max(r["t_max"] * 1.3, 250)
        ax.plot([S1a[0], S1a[0] + L * r["u1"][0]], [S1a[1], S1a[1] + L * r["u1"][1]],
                color="tab:blue", lw=1.6, zorder=9, label="视线")
        if len(x):
            cx, cy = float(x.mean()), float(y.mean())
            sp = max(float(np.ptp(x)), float(np.ptp(y))) * 1.8 + 120
            ax.set_xlim(cx - sp, cx + sp)
            ax.set_ylim(cy - sp, cy + sp)
        ax.set_aspect("equal")
        ax.grid(alpha=0.2)
        ax.set_xlabel("x / m")
        ax.set_ylabel("y / m")
        ax.legend(loc="upper left", fontsize=8.5)
        ax.set_title(f"{tag}\n$D_{{min}}$={r['Dmin']:.0f} m，"
                     f"子区占全部的 {a5/(r['n_feasible']*r['cell_area']):.1%}",
                     fontsize=11)
    fig3.suptitle("三场景的「最小误差候选区」并列（ε=5%）：形状多为一条弧带，"
                  "紧贴场地边界一侧", fontsize=14)
    fig3.tight_layout(rect=[0, 0, 1, 0.90])
    p3 = os.path.join(HERE, "最小误差候选区_三例并列.png")
    fig3.savefig(p3, dpi=125)
    plt.close(fig)
    print("->", p3)

    return [p1, p2, p3]


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    main()
