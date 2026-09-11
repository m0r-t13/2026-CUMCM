"""源可行段"在不在半径上"的对比。

定义
----
源可行段 = 从 S1 沿视线方向 u1 出发、落在目标圆域内、且在接收半径内的一段。

"在半径上" ⟺ 过 S1 沿 u1 的**整条直线穿过圆心**  ⟺ 直线到圆心的垂距 d = 0
                                                      （d = |S1 × u1|）

之前所有例子都取 S1 在 y 轴上、视线朝正北，此时 d = 0，左右完全对称。
本脚本给三个例子：

  (a) 沿半径朝外      S1=(0,900),  θ1=90°   d=0      左右对称
  (b) 斜穿圆域        S1=(1200,700), θ1=120° d>0      楔两条边出圆域的距离不同
  (c) 近切向斜穿      S1=(1500,900), θ1=120° d≈|S1|   不对称最大

关键后果
--------
d ≠ 0 时：
  * S1 不在弦的中点（弦中点到 S1 的距离 = |S1·u1|）
  * 楔的两条边界穿出圆域的距离**不相等** → 远端两个角点不等距
  * 于是两个"远端约束圆"的错位不再对称 → 候选区域相对视线**倾斜**
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
R_CAP = 1000.0
T_NEAR = 5.0
ALPHA = 1.0


def unit(deg):
    r = math.radians(deg)
    return np.array([math.cos(r), math.sin(r)])


def t_arena_along(S1, phi_deg, arena=ARENA):
    u = unit(phi_deg)
    b0 = float(S1 @ u)
    disc = b0 * b0 + arena * arena - float(S1 @ S1)
    if disc <= 0:
        return None
    return -b0 + math.sqrt(disc)


def analyse(S1, th, cap=R_CAP, alpha=ALPHA, t_recv=1500.0):
    S1 = np.asarray(S1, float)
    u1 = unit(th)
    d = abs(float(S1[0] * u1[1] - S1[1] * u1[0]))      # 直线到圆心的垂距
    b0 = float(S1 @ u1)                                # 弦中点到 S1 的有向距离
    t_ar = t_arena_along(S1, th)
    t_lo = t_arena_along(S1, th - alpha)
    t_hi = t_arena_along(S1, th + alpha)
    t_max = min(t_recv, t_ar)
    return dict(S1=S1, th=th, u1=u1, d=d, b0=b0, t_arena=t_ar,
                t_hi=t_hi, t_lo=t_lo, t_max=t_max,
                chord_half=math.sqrt(max(0.0, ARENA**2 - d**2)),
                asym=abs(t_hi - t_lo))


def main():
    cases = [
        ("(a) 视线沿半径朝外", (0.0, 900.0), 90.0),
        ("(b) 斜穿圆域", (1200.0, 700.0), 120.0),
        ("(c) 近切向斜穿", (1500.0, 900.0), 120.0),
    ]
    R = []
    print("=" * 100)
    print("源可行段「在不在半径上」的对比")
    print("=" * 100)
    print(f"{'情形':<22}{'|S1|':>8}{'垂距 d':>9}{'弦半长':>9}"
          f"{'t_arena':>9}{'左边':>8}{'右边':>8}{'不对称':>9}{'t_max':>8}")
    print("-" * 100)
    for tag, S1, th in cases:
        r = analyse(S1, th)
        R.append((tag, r))
        print(f"{tag:<22}{np.linalg.norm(S1):>8.1f}{r['d']:>9.2f}"
              f"{r['chord_half']:>9.1f}{r['t_arena']:>9.1f}"
              f"{r['t_lo']:>8.1f}{r['t_hi']:>8.1f}"
              f"{r['asym']:>9.2f}{r['t_max']:>8.1f}")
        print(f"{'':<22}弦中点到 S1 的距离 = |S1·u1| = {abs(r['b0']):.2f} m"
              f"   （S1 偏离弦中点这么多）")

    print("\n判别：")
    for tag, r in R:
        verdict = "在半径上（d=0，左右对称）" if r["d"] < 1e-9 else \
                  f"**不在半径上**（d={r['d']:.1f} m）"
        print(f"  {tag:<22}{verdict}")

    # ---------------- 画图 ----------------
    fig, axes = plt.subplots(1, 4, figsize=(25, 6.6))
    thv = np.linspace(0, 2 * math.pi, 800)

    for ax, (tag, r) in zip(axes[:3], R):
        S1, u1, th = r["S1"], r["u1"], r["th"]
        ax.plot(ARENA * np.cos(thv), ARENA * np.sin(thv), color="0.45",
                ls="--", lw=1.4, zorder=3, label="目标区域 (r=1800 m)")

        # 弦所在直线 + 弦
        u_perp = np.array([-u1[1], u1[0]])
        L = 2600.0
        ax.plot([S1[0] - L * u1[0], S1[0] + L * u1[0]],
                [S1[1] - L * u1[1], S1[1] + L * u1[1]],
                color="0.6", lw=0.8, ls="-", zorder=2)
        A = S1 - r["b0"] * u1 - r["chord_half"] * u1
        B = S1 - r["b0"] * u1 + r["chord_half"] * u1
        ax.plot([A[0], B[0]], [A[1], B[1]], color="0.25", lw=4, alpha=0.55,
                zorder=4, label="圆域内的弦")
        mid = S1 - r["b0"] * u1
        ax.plot(*mid, "s", color="0.25", ms=9, zorder=9, label="弦中点（圆心垂足）")
        # 垂足到圆心的虚线
        ax.plot([0, mid[0]], [0, mid[1]], color="0.25", ls=":", lw=1.2, zorder=4)
        ax.annotate(f"$d$={r['d']:.1f} m", (mid[0] / 2, mid[1] / 2),
                    textcoords="offset points", xytext=(6, 8),
                    color="0.25", fontsize=10)

        # 楔（放大显示）
        AX = 9.0
        for s in (-1, 1):
            uu = unit(th + s * AX)
            ax.plot([S1[0], S1[0] + 950 * uu[0]], [S1[1], S1[1] + 950 * uu[1]],
                    color="tab:blue", lw=0.9, alpha=0.5, zorder=5)
        # 两条真实楔边界（到出圆域为止）
        for s, t_exit, col, off in ((-1, r["t_lo"], "tab:purple", (-8, -22)),
                                    (+1, r["t_hi"], "tab:orange", (10, 10))):
            uu = unit(th + s * ALPHA)
            P = S1 + t_exit * uu
            ax.plot([S1[0], P[0]], [S1[1], P[1]], color=col, lw=1.8, zorder=7)
            ax.plot(*P, "o", color=col, ms=9, mec="k", mew=0.7, zorder=10)
            ax.annotate(f"{t_exit:.0f} m", P, textcoords="offset points",
                        xytext=off, color=col, fontsize=10, zorder=11,
                        ha="right" if off[0] < 0 else "left",
                        bbox=dict(fc="white", alpha=0.75, ec="none", pad=1.0))
        ax.plot([], [], color="tab:purple", lw=1.8, label="下侧楔边（真实 $\\alpha$）")
        ax.plot([], [], color="tab:orange", lw=1.8, label="上侧楔边（真实 $\\alpha$）")

        # 可行段
        P_lo = S1 + T_NEAR * u1
        P_hi = S1 + r["t_max"] * u1
        ax.plot([P_lo[0], P_hi[0]], [P_lo[1], P_hi[1]], color="lime", lw=6,
                solid_capstyle="butt", zorder=11, label=f"源可行段 [5, {r['t_max']:.0f}] m")
        ax.plot(*P_lo, "o", color="lime", ms=8, mec="k", mew=0.6, zorder=12)
        ax.plot(*P_hi, "o", color="lime", ms=8, mec="k", mew=0.6, zorder=12)

        ax.plot(*S1, "o", color="k", ms=11, zorder=12)
        ax.annotate("$S_1$", S1, textcoords="offset points", xytext=(-26, -34),
                    fontsize=13, zorder=13)
        ax.plot(0, 0, "P", color="k", ms=10, zorder=12)
        ax.annotate("圆心", (0, 0), textcoords="offset points", xytext=(-16, -32),
                    fontsize=11, zorder=13)

        lim = ARENA * 1.08
        ax.set_aspect("equal")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.grid(alpha=0.2)
        ax.set_xlabel("x / m")
        ax.set_ylabel("y / m")
        ax.legend(loc="lower left", fontsize=8.0)
        ax.set_title(f"{tag}\n$S_1$=({S1[0]:.0f},{S1[1]:.0f}), θ1={th:.0f}°, "
                     f"$d$={r['d']:.1f} m\n"
                     f"两楔边出圆域距离差 {r['asym']:.1f} m", fontsize=11)

    # (d) 不平衡度量对比
    ax = axes[3]
    labels = [t.split(" ")[0] for t, _ in R]
    asym = [r["asym"] for _, r in R]
    dd = [r["d"] for _, r in R]
    x = np.arange(len(labels))
    b1 = ax.bar(x - 0.2, asym, width=0.4, color="tab:blue",
                label="两楔边出圆域距离差 (m)")
    ax.set_ylabel("两楔边出圆域距离差 / m", color="tab:blue")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.grid(alpha=0.25, axis="y")
    ax2 = ax.twinx()
    ax2.plot(x, dd, "o-", color="red", lw=2, ms=10, label="直线到圆心垂距 $d$ (m)")
    ax2.set_ylabel("直线到圆心垂距 $d$ / m", color="red")
    for xi, (a_, d_) in enumerate(zip(asym, dd)):
        ax.annotate(f"{a_:.1f}", (xi - 0.2, a_), textcoords="offset points",
                    xytext=(0, 5), ha="center", fontsize=11, color="tab:blue")
        ax2.annotate(f"{d_:.0f}", (xi, d_), textcoords="offset points",
                     xytext=(0, -20), ha="center", fontsize=11, color="red")
    ax.set_title("(d) 越偏离半径，两侧越不对称", fontsize=11)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=9)

    fig.suptitle("源可行段「在不在半径上」：$d=0$ 时左右对称；$d>0$（斜穿）时"
                 "楔两条边出圆域的距离不同，远端两个角点不再等距", fontsize=14)
    fig.subplots_adjust(left=0.04, right=0.99, top=0.84, bottom=0.09, wspace=0.22)
    out = os.path.join(HERE, "可行段是否在半径上.png")
    fig.savefig(out, dpi=135)
    plt.close(fig)
    print("\n->", out)


if __name__ == "__main__":
    import sys
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    main()
