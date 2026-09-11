"""正文用图：问题二模型总览（三面板示意图）。

(a) 一次观测 → 源可行域是一个有厚度的角楔（不是线段）
(b) 探测保证 → 无穷多个圆退化成「每个顶点一个圆」
(c) 两条选址策略 → ①垂直于视线 ②朝场地中心一侧
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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from screen_min_error import (  # noqa: E402
    unit, t_arena_along, source_vertices, edge_source,
    ARENA, R_MIN, T_RECV, T_NEAR, ALPHA,
)

for fam in ("Microsoft YaHei", "SimHei", "SimSun", "DejaVu Sans"):
    try:
        font_manager.findfont(fam, fallback_to_default=False)
        plt.rcParams["font.sans-serif"] = [fam]
        break
    except Exception:
        continue
plt.rcParams["axes.unicode_minus"] = False

THV = np.linspace(0, 2 * math.pi, 800)


def arena(ax, lab="目标区域 (R=1800 m)"):
    ax.plot(ARENA * np.cos(THV), ARENA * np.sin(THV), color="0.45", ls="--",
            lw=1.4, zorder=4, label=lab)


def main():
    S1 = np.array([0.0, 0.0])
    TH = 0.0
    u1, n1 = unit(TH), unit(TH + 90.0)

    fig, axes = plt.subplots(1, 3, figsize=(22, 7.4))

    # ================= (a) 源可行域：有厚度的角楔 =================
    ax = axes[0]
    arena(ax)
    L = 1500.0
    # 夸张显示厚度的对照
    AXBIG = 7.0
    for s in (-1, 1):
        ub = unit(TH + s * AXBIG)
        ax.plot([S1[0], S1[0] + L * ub[0]], [S1[1], S1[1] + L * ub[1]],
                color="0.75", lw=1.0, ls="--", zorder=3)
    ax.fill([0, L * math.cos(math.radians(-AXBIG)), L * math.cos(math.radians(AXBIG))],
            [0, L * math.sin(math.radians(-AXBIG)), L * math.sin(math.radians(AXBIG))],
            color="0.85", alpha=0.45, zorder=2)
    ax.annotate(f"（虚线为把 $\\pm1°$ 夸大到 $\\pm{AXBIG:.0f}°$ 的示意）",
                (700, -1130), fontsize=10, color="0.4", ha="center")

    # 真实可行域（加黑边，突出它是有厚度的细长条）
    Vc4 = source_vertices(S1, TH, T_RECV)
    poly = np.array([Vc4[0], Vc4[2], Vc4[3], Vc4[1]])
    ax.fill(poly[:, 0], poly[:, 1], color="lime", alpha=0.95, zorder=9,
            label="源可行域（真实 $\\pm1°$ 厚度）")
    ax.plot([poly[0, 0], poly[1, 0]], [poly[0, 1], poly[1, 1]],
            color="k", lw=1.8, zorder=10)
    ax.plot([poly[2, 0], poly[3, 0]], [poly[2, 1], poly[3, 1]],
            color="k", lw=1.8, zorder=10)
    # 夸大版可行域，便于看出"厚度"
    xb = [0.0, T_RECV * math.cos(math.radians(AXBIG)),
          T_RECV * math.cos(math.radians(AXBIG)), 0.0]
    yb = [0.0, T_RECV * math.sin(math.radians(AXBIG)),
          T_RECV * math.sin(math.radians(-AXBIG)), 0.0]
    ax.fill(xb, yb, color="lime", alpha=0.22, zorder=8,
            label=f"同上，厚度放大到 $\\pm{AXBIG:.0f}°$（仅示意）")

    ax.plot(*S1, "o", color="k", ms=11, zorder=12)
    ax.annotate("$S_1$", S1, textcoords="offset points", xytext=(-26, -36),
                fontsize=14, zorder=13)
    ax.annotate("", xy=(1100, 0), xytext=(120, 0),
                arrowprops=dict(arrowstyle="-|>", color="tab:blue", lw=2.2),
                zorder=11)
    ax.annotate("示向度 $\\theta_1$", (1000, 70), color="tab:blue", fontsize=12)
    for V in (Vc4[2], Vc4[3]):
        ax.plot(*V, "x", color="k", ms=11, mew=2.6, zorder=13)
    ax.annotate("远端两角点横向错开 $\\rho\\sin\\alpha$\n（典型 10~60 m）"
                "——这就是「厚度」",
                (760, 620), fontsize=11, ha="center", zorder=14,
                bbox=dict(fc="yellow", alpha=0.9, ec="k", lw=0.9))
    ax.annotate("", xy=(1499, 40), xytext=(1010, 560),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=1.6), zorder=13)
    ax.set_xlim(-450, 1900)
    ax.set_ylim(-1250, 1250)
    ax.set_title("(a) 一次观测 → 源可行域是「有厚度的角楔」", fontsize=13)
    ax.legend(loc="lower left", fontsize=9)

    # ================= (b) 顶点法 =================
    S1b = np.array([0.0, 900.0])
    THb = 90.0
    ax = axes[1]
    arena(ax)
    u1b, n1b = unit(THb), unit(THb + 90.0)
    Vcb = source_vertices(S1b, THb, R_MIN)
    polyb = np.array([Vcb[0], Vcb[2], Vcb[3], Vcb[1]])
    ax.fill(polyb[:, 0], polyb[:, 1], color="lime", alpha=0.9, zorder=8,
            label="源可行域（$r\\leq R_{min}$ 部分）")
    # 顶点圆
    ang = np.linspace(0, 2 * math.pi, 400)
    for k, V in enumerate(Vcb):
        ax.plot(V[0] + R_MIN * np.cos(ang), V[1] + R_MIN * np.sin(ang),
                color="tab:red" if k in (2, 3) else "tab:purple",
                ls=":", lw=1.6, alpha=0.9, zorder=5)
    ax.plot([], [], color="tab:purple", ls=":", lw=1.8,
            label="近端顶点圆（2 个，相距仅 0.17 m）")
    ax.plot([], [], color="tab:red", ls=":", lw=1.8,
            label="远端顶点圆（2 个，绝不能合并）")
    ax.plot(Vcb[:, 0], Vcb[:, 1], "x", color="k", ms=11, mew=2.6, zorder=12)
    ax.annotate("对「每个顶点」各画一个 $R_{min}$ 圆\n"
                "→ 无穷多条约束退化成有限个",
                (Vcb[2][0], Vcb[2][1]), textcoords="offset points",
                xytext=(-70, 70), fontsize=11, zorder=14,
                bbox=dict(fc="yellow", alpha=0.85, ec="k", lw=0.9))
    ax.plot(*S1b, "o", color="k", ms=11, zorder=12)
    ax.annotate("$S_1$", S1b, textcoords="offset points", xytext=(-28, -36),
                fontsize=14, zorder=13)
    ax.set_xlim(-1750, 1750)
    ax.set_ylim(-1750, 1750)
    ax.set_aspect("equal")
    ax.set_title("(b) 探测保证 → 「每个顶点一个圆」", fontsize=13)
    ax.legend(loc="lower left", fontsize=9)

    # ================= (c) 两条策略 =================
    S1c = np.array([1200.0, 700.0])
    THc = 120.0
    ax = axes[2]
    arena(ax)
    u1c, n1c = unit(THc), unit(THc + 90.0)
    Vcc = source_vertices(S1c, THc, R_MIN)
    polyc = np.array([Vcc[0], Vcc[2], Vcc[3], Vcc[1]])
    ax.fill(polyc[:, 0], polyc[:, 1], color="lime", alpha=0.9, zorder=8,
            label="源可行域")
    t_bar = 0.5 * (T_NEAR + min(T_RECV, t_arena_along(S1c, THc)))
    x0 = S1c + t_bar * u1c
    LP = 820.0
    ax.plot([S1c[0] - LP * u1c[0], S1c[0] + LP * u1c[0]],
            [S1c[1] - LP * u1c[1], S1c[1] + LP * u1c[1]],
            color="tab:blue", lw=2.0, alpha=0.9, zorder=10, label="视线（方向已定）")
    ax.plot([x0[0] - LP * n1c[0], x0[0] + LP * n1c[0]],
            [x0[1] - LP * n1c[1], x0[1] + LP * n1c[1]],
            color="tab:red", lw=2.4, zorder=11,
            label="策略①：垂直于视线（$90°$ 交会）")
    ax.plot(*x0, "s", color="tab:red", ms=9, zorder=12)
    # 朝心箭头
    tc = -S1c / np.linalg.norm(S1c)
    ax.annotate("", xy=(S1c[0] + 560 * tc[0], S1c[1] + 560 * tc[1]),
                xytext=(S1c[0] + 130 * tc[0], S1c[1] + 130 * tc[1]),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=2.6), zorder=13)
    ax.annotate("策略②：选「朝场地中心」的一侧\n"
                "（背离圆心的一侧被边界挤死）",
                (S1c[0] + 600 * tc[0], S1c[1] + 600 * tc[1]),
                textcoords="offset points", xytext=(-20, -66), fontsize=11,
                ha="center", zorder=14,
                bbox=dict(fc="yellow", alpha=0.88, ec="k", lw=0.9))
    # 两侧标注（数值取自本场景实算：S1=(1200,700), θ1=120°）
    for s, txt, col in ((+1, "＋侧（朝心）\n$\\gamma_{min}$≈54.6°", "tab:green"),
                        (-1, "−侧（背心）\n$\\gamma_{min}$≈28.4°", "tab:orange")):
        P = x0 + s * 430 * n1c
        ax.plot(*P, "o", color=col, ms=13, mec="k", mew=1.0, zorder=12)
        ax.annotate(txt, P, textcoords="offset points", xytext=(0, 22),
                    fontsize=10.5, color=col, ha="center", zorder=13)
    ax.plot(*S1c, "o", color="k", ms=11, zorder=12)
    ax.annotate("$S_1$", S1c, textcoords="offset points", xytext=(-30, -36),
                fontsize=14, zorder=13)
    ax.plot(0, 0, "P", color="k", ms=11, zorder=12)
    ax.annotate("圆心", (0, 0), textcoords="offset points", xytext=(-18, -36),
                fontsize=12, zorder=13)
    ax.set_xlim(-1750, 1750)
    ax.set_ylim(-1750, 1750)
    ax.set_aspect("equal")
    ax.set_title("(c) 两条选址策略", fontsize=13)
    ax.legend(loc="lower left", fontsize=9)

    for ax in axes:
        ax.grid(alpha=0.2)
        ax.set_xlabel("x / m")
        ax.set_ylabel("y / m")

    fig.suptitle("问题二模型总览：源可行域（含厚度）→ 探测保证（顶点法）→ 选址策略",
                 fontsize=15)
    fig.subplots_adjust(left=0.04, right=0.99, top=0.87, bottom=0.08, wspace=0.20)
    out = os.path.join(HERE, "图0_模型总览.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("->", out)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    main()
