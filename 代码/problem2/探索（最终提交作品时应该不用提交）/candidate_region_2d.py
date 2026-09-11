"""第二问：把角楔当成**有厚度的二维区域**时的候选区域求解。

与 candidate_region.py（把楔形压成线段）的区别
-----------------------------------------------
线段版：源可行集 = 一条线段，端点 2 个 → 2 个约束圆。
二维版：源可行集 = 角楔 ∩ 圆域 ∩ 半径帽，是一个**有厚度的细长多边形**。

核心思路（凸性）
----------------
"S2 必须对所有可能的源位置都保证探测" 即

        max_{G in T} |S2 - G| <= R_min

被 max 的 |S2-G| 关于 G 是**凸函数**；凸函数在紧凸集上的最大值必在**极点**
取到。于是无穷多条约束退化成"每个顶点一条"：

        ∩_{G in T} D(G,R) = ∩_{V in 顶点(T)} D(V,R)

角楔∩圆域∩半径帽 的顶点只有 4 个：
    近端 2 个 (r = 5, 夹角 ±alpha)   —— 相距 2*5*sin(alpha) ≈ 0.17 m, 几乎重合
    远端 2 个 (r = min(帽, t_arena)) —— 横向各偏 ±r*sin(alpha)  ← "厚度"的全部影响

目标函数同理：最坏交会角也在**源区域的顶点**取到（用的是完整源区域，
半径帽取 1500 而不是 1000）。
"""

from __future__ import annotations

import math
import os

import numpy as np
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

HERE = os.path.dirname(os.path.abspath(__file__))

ARENA = 1800.0
R_MIN = 1000.0
T_RECV = 1500.0
T_NEAR = 5.0
ALPHA = 1.0


def unit(deg):
    r = math.radians(deg)
    return np.array([math.cos(r), math.sin(r)])


def t_arena_along(S1, u, arena=ARENA):
    """从 S1 沿方向 u 走多远穿出以原点为心、半径 arena 的圆。"""
    b0 = float(S1 @ u)
    disc = b0 * b0 + arena * arena - float(S1 @ S1)
    return -b0 + math.sqrt(max(0.0, disc))


def region_vertices(S1, theta1, r_cap, t_near=T_NEAR, alpha=ALPHA, arena=ARENA):
    """角楔 ∩ 圆域 ∩ {r <= r_cap} 的顶点（世界坐标）。

    顶点只有 4 个：两条楔边界各自与 {r = t_near} 和 {r = min(r_cap,t_arena)} 相交。
    返回 (近端两顶点, 远端两顶点)。
    """
    near, far = [], []
    for s in (+1.0, -1.0):
        u = unit(theta1 + s * alpha)
        t_ar = t_arena_along(S1, u, arena)
        r_far = max(0.0, min(r_cap, t_ar))
        near.append(S1 + t_near * u)
        if r_far > 1e-9:
            far.append(S1 + r_far * u)
    return np.array(near), np.array(far)


def crossing_angle(V, S1, S2):
    """在源位置 V 处，两条视线 S1-V 与 S2-V 的夹角（度）。"""
    v1 = S1 - V
    v2 = S2 - V
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 < 1e-9 or n2 < 1e-9:
        return 90.0
    c = float(np.clip(v1 @ v2 / (n1 * n2), -1.0, 1.0))
    return math.degrees(math.acos(c))


def solve_2d(S1, theta1, arena=ARENA, r_min=R_MIN, t_recv=T_RECV,
             t_near=T_NEAR, alpha=ALPHA, na=1400, nb=1400):
    """在 (a,b) 网格上搜索。

    目标不是"最大化交会角的最小值"（那会把 174° 这种近乎共线的坏几何当成好），
    而是**最小化"偏离 90° 的最大程度"**：
        worst_dev(S2) = max_{V in 顶点(源区域)} |gamma(V; S2) - 90°|
    """
    S1 = np.asarray(S1, float)
    u1 = unit(theta1)
    n1 = np.array([-u1[1], u1[0]])

    near_c, far_c = region_vertices(S1, theta1, r_min, t_near, alpha, arena)
    near_o, far_o = region_vertices(S1, theta1, t_recv, t_near, alpha, arena)
    Vc = np.vstack([near_c, far_c])
    Vo = np.vstack([near_o, far_o]) if len(far_o) else near_o
    t_max = min(t_recv, t_arena_along(S1, u1, arena))

    a = np.linspace(-arena, t_near + r_min, na)
    b = np.linspace(-r_min, r_min, nb)
    A, B = np.meshgrid(a, b)
    S2x = S1[0] + A * u1[0] + B * n1[0]
    S2y = S1[1] + A * u1[1] + B * n1[1]

    # 约束一：对每个约束顶点都要在 R_min 圆内
    ok = np.ones_like(A, dtype=bool)
    for V in Vc:
        ok &= ((S2x - V[0]) ** 2 + (S2y - V[1]) ** 2) <= r_min ** 2
    # 约束二：S2 自身在目标圆域内
    ok &= (S2x ** 2 + S2y ** 2) <= arena ** 2

    # 目标：完整源区域各顶点处，"偏离 90° 的程度"的最大值
    dev = np.zeros_like(A)
    for V in Vo:
        v1 = S1 - V
        n1v = np.linalg.norm(v1)
        v2x = S2x - V[0]
        v2y = S2y - V[1]
        n2v = np.hypot(v2x, v2y)
        c = np.clip((v1[0] * v2x + v1[1] * v2y) / (n1v * n2v), -1.0, 1.0)
        g = np.degrees(np.arccos(c))
        dev = np.maximum(dev, np.abs(g - 90.0))

    dev = np.where(ok, dev, np.nan)
    k = np.nanargmin(dev)
    i, j = np.unravel_index(k, dev.shape)
    a_star, b_star = float(a[j]), abs(float(b[i]))
    dev_star = float(dev[i, j])

    # 用于显示的 gamma_min 场
    gmin = np.full_like(A, 1e9)
    for V in Vo:
        v1 = S1 - V
        n1v = np.linalg.norm(v1)
        v2x = S2x - V[0]
        v2y = S2y - V[1]
        n2v = np.hypot(v2x, v2y)
        c = np.clip((v1[0] * v2x + v1[1] * v2y) / (n1v * n2v), -1.0, 1.0)
        gmin = np.minimum(gmin, np.degrees(np.arccos(c)))
    gmin = np.where(ok, gmin, np.nan)

    return dict(
        S1=S1, u1=u1, n1=n1, theta1=theta1,
        near_c=near_c, far_c=far_c, near_o=near_o, far_o=far_o,
        Vc=Vc, Vo=Vo, t_max=t_max,
        a=a, b=b, A=A, B=B, S2x=S2x, S2y=S2y,
        gam=gmin, dev=dev,
        a_star=a_star, b_star=b_star, dev_star=dev_star,
        gamma_star=90.0 - dev_star,
        S2p=S1 + a_star * u1 + b_star * n1,
        S2m=S1 + a_star * u1 - b_star * n1,
    )


def main():
    SCEN = [
        ("A", (0.0, 0.0), 0.0),
        ("B", (0.0, 1100.0), 90.0),
        ("C", (0.0, 1500.0), 90.0),
    ]
    print("=" * 100)
    print("第二问（角楔有厚度版）：约束顶点 → 有限个圆")
    print("=" * 100)
    res = {}
    for tag, S1, th in SCEN:
        r = solve_2d(S1, th)
        res[tag] = (S1, th, r)
        print(f"\n场景 {tag}:  S1={S1}, θ1={th}°")
        print(f"  约束顶点 (R_min={R_MIN:.0f} m)：")
        for nm, V in (("近+", r["near_c"][0]), ("近-", r["near_c"][1]),
                      ("远+", r["far_c"][0]), ("远-", r["far_c"][1])):
            lam = np.linalg.norm(V - np.asarray(S1))
            print(f"    {nm}: ({V[0]:9.3f}, {V[1]:9.3f})   r={lam:8.3f}")
        d_far = float(np.linalg.norm(r["far_c"][0] - r["far_c"][1]))
        print(f"  远端两顶点间距 = {d_far:.3f} m  ← 这就是角楔的'厚度'")
        print(f"  目标顶点数 = {len(r['Vo'])}（半径帽 {T_RECV:.0f} m）")
        print(f"  最优: a*={r['a_star']:.1f}, b*={r['b_star']:.1f}  "
              f"γ_min={r['gamma_star']:.2f}°（最大偏离 90°={r['dev_star']:.2f}°）")
        print(f"        S2+ = ({r['S2p'][0]:8.1f}, {r['S2p'][1]:8.1f})  "
              f"|S2+|={np.linalg.norm(r['S2p']):7.1f}")

    # ---------- 图 ----------
    fig, axes = plt.subplots(1, 3, figsize=(21, 7.6))
    levels = np.linspace(0, 90, 46)
    for ax, (tag, S1, th) in zip(axes, SCEN):
        r = res[tag][2]
        thv = np.linspace(0, 2 * math.pi, 600)
        ax.plot(ARENA * np.cos(thv), ARENA * np.sin(thv), color="0.45", ls="--",
                lw=1.3, zorder=6, label="目标区域 (r=1800 m)")

        u1, n1 = r["u1"], r["n1"]
        S1a = np.asarray(S1)

        # 楔（放大显示厚度）
        AX = 9.0
        L = max(r["t_max"] * 1.35, 200.0)
        for s in (-1, 1):
            uu = unit(th + s * AX)
            ax.plot([S1a[0], S1a[0] + L * uu[0]], [S1a[1], S1a[1] + L * uu[1]],
                    color="tab:blue", lw=1.0, alpha=0.7, zorder=8)
        ax.plot([S1a[0], S1a[0] + L * u1[0]], [S1a[1], S1a[1] + L * u1[1]],
                color="tab:blue", lw=1.4, alpha=0.9, zorder=8)

        # 真实可观测的源区域（用真实 alpha，填充）
        rr = np.linspace(T_NEAR, r["t_max"], 60)
        for s in (-1, 1):
            uu = unit(th + s * ALPHA)
            P = S1a[None, :] + rr[:, None] * uu[None, :]
            ax.plot(P[:, 0], P[:, 1], color="lime", lw=0.8, alpha=0.8, zorder=9)
        # 用细长四边形示意（真实厚度已按真实 alpha 画）
        P1 = S1a + T_NEAR * unit(th + ALPHA)
        P2 = S1a + T_NEAR * unit(th - ALPHA)
        poly = np.array([P1, r["far_c"][0], r["far_c"][1], P2])
        ax.fill(poly[:, 0], poly[:, 1], color="lime", alpha=0.85, zorder=9,
                label="源可行域（角楔∩圆域∩r≤1000）")

        # 约束顶点
        for V in r["Vc"]:
            ax.add_patch(plt.Circle(V, R_MIN, color="0.1", fill=False, ls=":",
                                    lw=1.2, zorder=7))
        ax.plot(r["Vc"][:, 0], r["Vc"][:, 1], "x", color="k", ms=9, mew=2, zorder=12)
        ax.plot([], [], color="0.1", ls=":", lw=1.2,
                label=f"顶点约束圆 (R={R_MIN:.0f} m)")

        # 候选区域（对每个格子算 γ_min）
        G = r["gam"]
        ax.contourf(r["S2x"], r["S2y"], np.where(np.isnan(G), np.nan, G),
                    levels=levels, cmap="turbo", vmin=0, vmax=90, alpha=0.95,
                    zorder=1)
        try:
            cs = ax.contour(r["S2x"], r["S2y"], G, levels=[20, 30, 40, 50, 60, 70],
                            colors="k", linewidths=0.7, alpha=0.5, zorder=6)
            ax.clabel(cs, fmt="%d°", fontsize=8)
        except Exception:
            pass

        for S2, lab in ((r["S2p"], "$S_2^+$"), (r["S2m"], "$S_2^-$")):
            ax.plot(*S2, "*", color="w", ms=27, zorder=13)
            ax.plot(*S2, "*", color="red", ms=20, zorder=14)
            ax.annotate(lab, S2, textcoords="offset points", xytext=(12, 8),
                        color="red", fontsize=13, zorder=15)
        ax.plot(*S1a, "o", color="k", ms=9, zorder=13)
        ax.annotate("$S_1$", S1a, textcoords="offset points", xytext=(-28, -34),
                    fontsize=13, zorder=15)

        lim = ARENA * 1.06
        ax.set_aspect("equal")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.grid(alpha=0.18)
        ax.set_xlabel("x / m")
        ax.set_ylabel("y / m")
        d_far = float(np.linalg.norm(r["far_c"][0] - r["far_c"][1]))
        ax.set_title(f"({tag}) $S_1$=({S1[0]:.0f},{S1[1]:.0f}), $\\theta_1$={th:.0f}°\n"
                     f"源区远端厚度 {d_far:.1f} m → $\\gamma_{{min}}$="
                     f"{r['gamma_star']:.1f}°", fontsize=11.5)
    axes[0].legend(loc="lower left", fontsize=8.5, framealpha=0.93)
    sm = plt.cm.ScalarMappable(cmap="turbo", norm=plt.Normalize(0, 90))
    sm.set_array([])
    fig.subplots_adjust(left=0.035, right=0.905, top=0.80, bottom=0.07, wspace=0.13)
    cax = fig.add_axes([0.917, 0.12, 0.011, 0.66])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("最坏交会角 $\\gamma_{min}$ / 度", fontsize=11)
    fig.suptitle("角楔有厚度时的候选区域：把无穷多个圆退化成「每个顶点一个圆」\n"
                 "★ = 最优点；远端两顶点横向错开，使候选区域略微收窄",
                 fontsize=14, y=0.965)
    out = os.path.join(HERE, "候选区域_有厚度版.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("\n->", out)


if __name__ == "__main__":
    import sys
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    main()
