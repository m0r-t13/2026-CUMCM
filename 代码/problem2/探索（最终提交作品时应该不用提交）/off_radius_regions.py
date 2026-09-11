"""三个「可行段是否在半径上」例子的第二问候选区域。

复用 candidate_region_2d 的顶点法：
    源可行域 = 角楔 ∩ 目标圆域 ∩ {r ≤ 半径帽}，顶点只有 4 个
    "对每个顶点各画一个半径 R_min 的圆" = 保证探测的可行集
    目标 = 最小化"偏离 90° 的最大程度"

与之前不同的地方：射向不再沿半径，左右**不对称**，所以要
分别求 b>0 侧和 b<0 侧各自的最优点。
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from candidate_region_2d import (  # noqa: E402
    unit, t_arena_along, region_vertices, ARENA, R_MIN, T_RECV, T_NEAR, ALPHA,
)

for fam in ("Microsoft YaHei", "SimHei", "SimSun", "DejaVu Sans"):
    try:
        font_manager.findfont(fam, fallback_to_default=False)
        plt.rcParams["font.sans-serif"] = [fam]
        break
    except Exception:
        continue
plt.rcParams["axes.unicode_minus"] = False

HERE = os.path.dirname(os.path.abspath(__file__))


def solve_case(S1, th, na=1300, nb=1300):
    S1 = np.asarray(S1, float)
    u1 = unit(th)
    n1 = np.array([-u1[1], u1[0]])

    near_c, far_c = region_vertices(S1, th, R_MIN)
    near_o, far_o = region_vertices(S1, th, T_RECV)
    Vc = np.vstack([near_c, far_c])
    Vo = np.vstack([near_o, far_o]) if len(far_o) else near_o
    t_max = min(T_RECV, t_arena_along(S1, u1))

    a = np.linspace(-ARENA, T_NEAR + R_MIN, na)
    b = np.linspace(-R_MIN, R_MIN, nb)
    A, B = np.meshgrid(a, b)
    S2x = S1[0] + A * u1[0] + B * n1[0]
    S2y = S1[1] + A * u1[1] + B * n1[1]

    ok = np.ones_like(A, dtype=bool)
    for V in Vc:
        ok &= ((S2x - V[0]) ** 2 + (S2y - V[1]) ** 2) <= R_MIN ** 2
    ok &= (S2x ** 2 + S2y ** 2) <= ARENA ** 2

    dev = np.zeros_like(A)
    gmin = np.full_like(A, 1e9)
    for V in Vo:
        v1 = S1 - V
        n1v = np.linalg.norm(v1)
        v2x = S2x - V[0]
        v2y = S2y - V[1]
        n2v = np.hypot(v2x, v2y)
        c = np.clip((v1[0] * v2x + v1[1] * v2y) / (n1v * n2v), -1.0, 1.0)
        g = np.degrees(np.arccos(c))
        dev = np.maximum(dev, np.abs(g - 90.0))
        gmin = np.minimum(gmin, g)

    dev = np.where(ok, dev, np.nan)
    gmin = np.where(ok, gmin, np.nan)

    out = {}
    for sign, key in ((+1, "plus"), (-1, "minus")):
        mask = (B > 0) if sign > 0 else (B < 0)
        d = np.where(mask, dev, np.nan)
        if np.all(np.isnan(d)):
            out[key] = None
            continue
        k = np.nanargmin(d)
        i, j = np.unravel_index(k, d.shape)
        a_s, b_s = float(a[j]), float(b[i])
        out[key] = dict(a=a_s, b=b_s, gamma=90.0 - float(dev[i, j]),
                        dev=float(dev[i, j]),
                        S2=S1 + a_s * u1 + b_s * n1)
    out.update(S1=S1, u1=u1, n1=n1, th=th, Vc=Vc, Vo=Vo, t_max=t_max,
               A=A, B=B, S2x=S2x, S2y=S2y, gmin=gmin, dev=dev,
               near_c=near_c, far_c=far_c)
    return out


def main():
    CASES = [
        ("(a) 视线沿半径朝外", (0.0, 900.0), 90.0),
        ("(b) 斜穿圆域", (1200.0, 700.0), 120.0),
        ("(c) 近切向斜穿", (1500.0, 900.0), 120.0),
    ]

    print("=" * 104)
    print("三例的第二问候选区域与最优第二检测点")
    print("=" * 104)
    print(f"{'情形':<20}{'d':>8}{'t_max':>7}{'角点间距':>9}"
          f"{'b+':>8}{'γ_min+':>9}{'b-':>8}{'γ_min-':>9}{'两侧差':>9}{'朝心侧':>8}")
    print("-" * 104)

    res = []
    for tag, S1, th in CASES:
        r = solve_case(S1, th)
        res.append((tag, S1, th, r))
        S1a = np.asarray(S1)
        d = abs(float(S1a[0] * r["u1"][1] - S1a[1] * r["u1"][0]))
        far_sep = float(np.linalg.norm(r["far_c"][0] - r["far_c"][1]))
        gp = r["plus"]["gamma"] if r["plus"] else float("nan")
        gm = r["minus"]["gamma"] if r["minus"] else float("nan")
        bp = r["plus"]["b"] if r["plus"] else float("nan")
        bm = r["minus"]["b"] if r["minus"] else float("nan")
        # "朝心侧"：+n1 方向与"指向圆心"的方向是否同向
        to_center = -S1a / max(np.linalg.norm(S1a), 1e-9)
        side = "+n1 侧" if float(to_center @ r["n1"]) > 0 else "-n1 侧"
        print(f"{tag:<20}{d:>8.1f}{r['t_max']:>7.1f}{far_sep:>9.1f}"
              f"{bp:>8.1f}{gp:>8.2f}°{bm:>8.1f}{gm:>8.2f}°"
              f"{abs(gp-gm):>8.2f}°{side:>8}")

    # ---------------- 画图：2 行 × 3 列 ----------------
    fig, axes = plt.subplots(2, 3, figsize=(22, 14.5))
    thv = np.linspace(0, 2 * math.pi, 800)
    levels = np.linspace(0, 90, 46)

    for col, (tag, S1, th, r) in enumerate(res):
        S1a = np.asarray(S1)
        u1, n1 = r["u1"], r["n1"]
        t_bar = 0.5 * (T_NEAR + r["t_max"])
        x0 = S1a + t_bar * u1

        for row in (0, 1):
            ax = axes[row, col]
            ax.plot(ARENA * np.cos(thv), ARENA * np.sin(thv), color="0.45",
                    ls="--", lw=1.2, zorder=5, label="目标区域 (r=1800)")

            # 源可行域（真实厚度）
            poly = np.vstack([r["near_c"][0], r["far_c"][0],
                              r["far_c"][1], r["near_c"][1]])
            ax.fill(poly[:, 0], poly[:, 1], color="lime", alpha=0.85, zorder=9,
                    label="源可行域（角楔∩圆域∩r≤1000）")

            # 视线 + 楔（放大显示）
            L = max(r["t_max"] * 1.5, 300.0)
            AX = 8.0
            for s in (-1, 1):
                uu = unit(th + s * AX)
                ax.plot([S1a[0], S1a[0] + L * uu[0]],
                        [S1a[1], S1a[1] + L * uu[1]],
                        color="tab:blue", lw=0.9, alpha=0.45, zorder=6)
            ax.plot([S1a[0], S1a[0] + L * u1[0]], [S1a[1], S1a[1] + L * u1[1]],
                    color="tab:blue", lw=1.5, alpha=0.9, zorder=6, label="视线")

            # 约束顶点与它们的圆
            for V in r["Vc"]:
                ax.add_patch(plt.Circle(V, R_MIN, color="0.15", fill=False,
                                        ls=":", lw=1.1, zorder=7))
            ax.plot(r["Vc"][:, 0], r["Vc"][:, 1], "x", color="k", ms=9, mew=2,
                    zorder=12)
            ax.plot([], [], color="0.15", ls=":", lw=1.1,
                    label=f"顶点约束圆 (R={R_MIN:.0f})")

            # 候选区域
            G = r["gmin"]
            ax.contourf(r["S2x"], r["S2y"], G, levels=levels, cmap="turbo",
                        vmin=0, vmax=90, alpha=0.95, zorder=1)
            try:
                cs = ax.contour(r["S2x"], r["S2y"], G,
                                levels=[30, 40, 50, 60, 70],
                                colors="k", linewidths=0.7, alpha=0.45, zorder=6)
                ax.clabel(cs, fmt="%d°", fontsize=7)
            except Exception:
                pass

            # 过 x0 的"垂直于视线"参考线（用来看候选区有没有歪）
            perp = 900.0
            ax.plot([x0[0] - perp * n1[0], x0[0] + perp * n1[0]],
                    [x0[1] - perp * n1[1], x0[1] + perp * n1[1]],
                    color="w", ls="-.", lw=1.6, zorder=13,
                    label="过 $x_0$ 垂直视线（对称时的候选线）")
            ax.plot(*x0, "s", color="k", ms=8, zorder=13)

            # 两个最优解
            for key, lab, mk in (("plus", "$S_2^+$", "*"),
                                 ("minus", "$S_2^-$", "P")):
                rr = r[key]
                if not rr:
                    continue
                ax.plot(*rr["S2"], mk, color="w",
                        ms=26 if mk == "*" else 16, zorder=14)
                ax.plot(*rr["S2"], mk, color="red",
                        ms=19 if mk == "*" else 12, zorder=15)
                ax.annotate(f"{lab}\n$\\gamma$={rr['gamma']:.1f}°", rr["S2"],
                            textcoords="offset points", xytext=(13, 8),
                            color="red", fontsize=10, zorder=16,
                            bbox=dict(fc="white", alpha=0.8, ec="none", pad=1.2))

            ax.plot(*S1a, "o", color="k", ms=10, zorder=14)
            ax.annotate("$S_1$", S1a, textcoords="offset points", xytext=(-26, -34),
                        fontsize=12, zorder=15)

            # 从 S1 指向圆心：提示哪一侧"有空间"
            tc = -S1a / max(np.linalg.norm(S1a), 1e-9)
            ax.annotate("", xy=(S1a[0] + 620 * tc[0], S1a[1] + 620 * tc[1]),
                        xytext=(S1a[0] + 150 * tc[0], S1a[1] + 150 * tc[1]),
                        arrowprops=dict(arrowstyle="-|>", color="k", lw=2.0,
                                        alpha=0.65), zorder=16)
            ax.annotate("朝圆心\n(有空间)", (S1a[0] + 640 * tc[0], S1a[1] + 640 * tc[1]),
                        color="k", fontsize=10, ha="center", va="center", zorder=17,
                        bbox=dict(fc="white", alpha=0.8, ec="none", pad=1.2))

            if row == 0:
                lim = ARENA * 1.06
                ax.set_xlim(-lim, lim)
                ax.set_ylim(-lim, lim)
                ax.set_title(f"{tag}\n$S_1$=({S1[0]:.0f},{S1[1]:.0f}), θ1={th:.0f}°, "
                             f"$t_{{max}}$={r['t_max']:.0f} m", fontsize=11.5)
            else:
                pts = [rr["S2"] for k in ("plus", "minus")
                       if (rr := r[k]) is not None]
                if pts:
                    mid = np.mean(pts, axis=0)
                    span = max(float(np.linalg.norm(p - mid)) for p in pts)
                    span = span * 1.5 + 80.0
                else:
                    mid = np.array([float(np.nanmean(r["S2x"])),
                                    float(np.nanmean(r["S2y"]))])
                    span = 480.0
                ax.set_xlim(mid[0] - span, mid[0] + span)
                ax.set_ylim(mid[1] - span, mid[1] + span)
                gl = r["plus"]["gamma"] if r["plus"] else float("nan")
                gr = r["minus"]["gamma"] if r["minus"] else float("nan")
                ax.set_title(f"{tag} —— 放大：两侧最优点 $\\gamma$ = "
                             f"{gl:.1f}° vs {gr:.1f}°", fontsize=11.5)
            ax.set_aspect("equal")
            ax.grid(alpha=0.18)
            ax.set_xlabel("x / m")
            ax.set_ylabel("y / m")
            if row == 0 and col == 0:
                ax.legend(loc="lower left", fontsize=8.0)

    sm = plt.cm.ScalarMappable(cmap="turbo", norm=plt.Normalize(0, 90))
    sm.set_array([])
    fig.subplots_adjust(left=0.035, right=0.925, top=0.92, bottom=0.055,
                        wspace=0.13, hspace=0.22)
    cax = fig.add_axes([0.937, 0.10, 0.009, 0.72])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("最坏交会角 $\\gamma_{min}$ / 度", fontsize=11)
    fig.suptitle("可行段「在不在半径上」→ 候选区域「正不正」：\n"
                 "(a) 沿半径，候选区关于视线对称；(b)(c) 斜穿，候选区相对视线歪斜、"
                 "两侧最优点不再镜像", fontsize=14, y=0.975)
    out = os.path.join(HERE, "三例候选区域对比.png")
    fig.savefig(out, dpi=125)
    plt.close(fig)
    print("\n->", out)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    main()
