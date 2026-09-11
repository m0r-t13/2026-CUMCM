"""远端截线到底长什么样：单弧 / 双弧+角点 / 多边形近似的折线。

背景
----
第二问里源可行域 = 角楔 ∩ 目标圆域 ∩ {r ≤ R_cap}（外加减掉 r ≤ 5 的小圆）。
能"截住远端"的只有两个圆形约束：
    目标圆域  |P| <= 1800        （圆心在原点）
    接收半径  |P - S1| <= R_cap  （圆心在 S1）
角楔的两条边界线是**沿视线方向**的，只约束两侧，约束不到远端。
所以远端只可能是圆弧。

但两圆谁"先到"是随方向变的：t_arena(phi) 与 R_cap 比较，在角楔张角范围内
可能一半被圆域截、一半被接收半径截 —— 于是远端由**两条不同半径的圆弧**拼成，
中间出现一个**角点**。这就是"不是单一圆弧"的唯一情形。

本脚本画三种情况：
  (a) 单弧   ：整段都被接收半径帽截        S1 在圆心附近
  (b) 双弧+角点：一半被帽截、一半被圆域截   S1 靠近圆域边界且方向合适
  (c) 折线   ：把圆域换成多边形近似（数值实现的常见做法）
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


def t_arena_along(S1, phi_rad, arena=ARENA):
    """从 S1 沿方向 phi 走多远穿出以原点为心、半径 arena 的圆。"""
    u = np.array([math.cos(phi_rad), math.sin(phi_rad)])
    b0 = float(S1 @ u)
    disc = b0 * b0 + arena * arena - float(S1 @ S1)
    if disc <= 0:
        return None
    return -b0 + math.sqrt(disc)


def t_poly_along(S1, phi_rad, poly):
    """沿方向 phi 与凸多边形 poly 的边界交点：返回 (距离, 命中的边序号)。"""
    u = np.array([math.cos(phi_rad), math.sin(phi_rad)])
    best, best_i = None, None
    m = len(poly)
    for i in range(m):
        a, b = poly[i], poly[(i + 1) % m]
        e = b - a
        det = u[0] * (-e[1]) - u[1] * (-e[0])
        if abs(det) < 1e-14:
            continue
        rhs = a - S1
        t = (rhs[0] * (-e[1]) - rhs[1] * (-e[0])) / det
        s = (u[0] * rhs[1] - u[1] * rhs[0]) / det
        if t > 0 and -1e-9 <= s <= 1 + 1e-9:
            if best is None or t < best:
                best, best_i = t, i
    return (best, best_i) if best is not None else (None, None)


def boundary(S1, theta_deg, cap=R_CAP, alpha=ALPHA, arena=ARENA,
             poly=None, n=600):
    """返回远端边界点、每点由谁限制、以及近端边界点。

    poly 为 None 时 who 取 "cap"/"arena"；
    给了 poly 时 who 取 "edge{k}"（第 k 条边），便于检出折线的角点。
    """
    lo, hi = theta_deg - alpha, theta_deg + alpha
    phis = np.linspace(math.radians(lo), math.radians(hi), n)
    far, who, near = [], [], []
    for ph in phis:
        u = np.array([math.cos(ph), math.sin(ph)])
        if poly is not None:
            t_ar, ei = t_poly_along(S1, ph, poly)
        else:
            t_ar, ei = t_arena_along(S1, ph, arena), None
        if t_ar is None or t_ar <= T_NEAR:
            continue
        r = min(cap, t_ar)
        far.append(S1 + r * u)
        near.append(S1 + T_NEAR * u)
        if cap < t_ar - 1e-9:
            who.append("cap")
        elif poly is not None:
            who.append(f"edge{ei}")
        else:
            who.append("arena")
    return np.array(far), who, np.array(near)


def find_kinks(who):
    """找出限制来源切换的位置（角点）。"""
    idx = [i for i in range(1, len(who)) if who[i] != who[i - 1]]
    return idx


def main():
    S1a = np.array([0.0, 0.0])
    S1b = np.array([0.0, 1500.0])

    cases = []

    # (a) 单弧：整段都被接收半径帽截
    far, who, near = boundary(S1a, 0.0)
    cases.append(dict(tag="(a) 单弧 —— 被接收半径帽截", S1=S1a, th=0.0,
                      far=far, who=who, near=near, kind="cap",
                      note="帽 1000 m < 圆域 1800 m → 远端是**以 $S_1$ 为心**的弧"))

    # (b) 单弧：整段都被目标圆域截
    far, who, near = boundary(S1b, 90.0)
    cases.append(dict(tag="(b) 单弧 —— 被目标圆域截", S1=S1b, th=90.0,
                      far=far, who=who, near=near, kind="arena",
                      note="视线 300 m 就出圆域，比帽近 → 远端是**以原点为心**的弧"))

    # (c) 双弧 + 角点
    far, who, near = boundary(S1b, 0.0)
    cases.append(dict(tag="(c) 两条弧 + 角点", S1=S1b, th=0.0,
                      far=far, who=who, near=near, kind="mixed",
                      note="一半被**帽**截、一半被**圆域**截 → 两种弧，中间一个角点"))

    # (d) 圆域用多边形近似 → 折线
    #     注意：多边形顶点是按「以原点为心」的角度排的，而角楔是按「以 S1 为心」
    #     的 2° 张开的。所以要先算出远端那段边界相对原点的方位，再把一个多边形
    #     顶点正好转到那个方位上，折线的角点才会落进楔内。
    N = 24
    far_c, _, _ = boundary(S1b, 0.0)
    ang_mean = math.atan2(float(far_c[:, 1].mean()), float(far_c[:, 0].mean()))
    step = 2 * math.pi / N
    rot = ang_mean - round(ang_mean / step) * step
    ang = np.linspace(0, 2 * math.pi, N, endpoint=False) + rot
    poly = np.column_stack([ARENA * np.cos(ang), ARENA * np.sin(ang)])
    far, who, near = boundary(S1b, 0.0, poly=poly)
    cases.append(dict(tag=f"(d) 折线 —— 圆域换成 {N} 边形", S1=S1b, th=0.0,
                      far=far, who=who, near=near, kind="poly", poly=poly,
                      note=f"用 {N} 边形近似圆域 → 远端变成**折线**（数值实现的产物）"))
    print(f"  [d] 多边形顶点方位对准 {math.degrees(ang_mean):.2f}°，"
          f"旋转量 {math.degrees(rot):.2f}°")

    print("=" * 92)
    print("远端截线形态统计")
    print("=" * 92)
    for c in cases:
        who = c["who"]
        n_cap = who.count("cap")
        n_ar = who.count("arena")
        print(f"\n{c['tag']}   S1=({c['S1'][0]:.0f},{c['S1'][1]:.0f}), θ1={c['th']}°")
        kinds = sorted(set(who))
        if c["kind"] == "poly":
            print(f"  远端采样点 {len(who)} 个，落在 {len(kinds)} 条多边形边上")
        else:
            print(f"  远端采样点 {len(who)} 个：被帽截 {who.count('cap')} 个，"
                  f"被圆域截 {who.count('arena')} 个")
        print(f"  {c['note']}")
        kk = find_kinks(who)
        if kk:
            print(f"  角点（限制来源切换）位置：")
            for i in kk:
                P = c["far"][i]
                r = np.linalg.norm(P - c["S1"])
                print(f"     ({P[0]:9.3f}, {P[1]:9.3f})   r={r:8.3f}  "
                      f"{who[i-1]} → {who[i]}")
        else:
            print(f"  无角点 → 单条圆弧")

    # 弧的"平坦程度"：弦高
    print("\n" + "=" * 92)
    print("圆弧其实非常平（这是它看起来像直线的原因）")
    print("=" * 92)
    f0 = cases[0]["far"]
    chord = float(np.linalg.norm(f0[0] - f0[-1]))
    sag_cap = R_CAP - math.sqrt(R_CAP**2 - (chord / 2) ** 2)
    print(f"  角楔在 r={R_CAP:.0f} 处的弦长 = {chord:.2f} m")
    print(f"  以 S1 为心的帽弧，弦高 = {sag_cap:.4f} m")
    print(f"  以原点为心的圆域弧（同弦长）弦高 = "
          f"{ARENA - math.sqrt(ARENA**2 - (chord/2)**2):.4f} m")
    print("  → 弦高只有几厘米到十几厘米，所以图上看着像直线，但它是弧。")

    # ---------------- 画图 ----------------
    fig, axes = plt.subplots(1, 4, figsize=(26, 6.8))
    for ax, c in zip(axes, cases):
        S1 = c["S1"]
        far, near = c["far"], c["near"]
        ax.fill(np.concatenate([near[:, 0], far[::-1, 0]]),
                np.concatenate([near[:, 1], far[::-1, 1]]),
                color="lime", alpha=0.55, zorder=3)

        # 远端边界按"谁限制"着色
        who = c["who"]
        for i in range(len(far) - 1):
            col = "red" if who[i] == "cap" else "blue"
            ax.plot(far[i:i + 2, 0], far[i:i + 2, 1], color=col, lw=3.2, zorder=8,
                    solid_capstyle="round")
        ax.plot([], [], color="red", lw=3,
                label=f"被**帽**截（圆心 $S_1$, R={R_CAP:.0f}）")
        ax.plot([], [], color="blue", lw=3, label="被**圆域**截（圆心原点, R=1800）")

        # 两个圆（局部画）
        thv = np.linspace(0, 2 * math.pi, 500)
        ax.plot(S1[0] + R_CAP * np.cos(thv), S1[1] + R_CAP * np.sin(thv),
                color="red", ls=":", lw=1.1, alpha=0.75, zorder=2)
        if c["kind"] == "poly":
            pp = np.vstack([c["poly"], c["poly"][:1]])
            ax.plot(pp[:, 0], pp[:, 1], color="blue", ls=":", lw=1.1, alpha=0.75,
                    zorder=2)
        elif np.linalg.norm(S1) < ARENA:
            ax.plot(ARENA * np.cos(thv), ARENA * np.sin(thv),
                    color="blue", ls=":", lw=1.1, alpha=0.75, zorder=2)

        # 角点
        kk = find_kinks(who)
        for i in kk:
            P = far[i]
            ax.plot(*P, "o", color="k", ms=13, mfc="yellow", mew=2, zorder=12)
            ax.annotate("角点", P, textcoords="offset points", xytext=(16, 10),
                        fontsize=11, zorder=13,
                        bbox=dict(fc="yellow", alpha=0.85, ec="k", lw=0.8))

        ax.plot(*S1, "o", color="k", ms=10, zorder=11)
        ax.annotate("$S_1$", S1, textcoords="offset points", xytext=(-18, -32),
                    fontsize=13, zorder=13)

        # 弦（两端点连线）：用来对照"弧偏离直线有多少"
        ax.plot([far[0, 0], far[-1, 0]], [far[0, 1], far[-1, 1]],
                color="0.35", ls="--", lw=1.2, zorder=6,
                label="两端点连线（弦）")
        d = far[-1] - far[0]
        nrm = np.array([-d[1], d[0]]) / max(np.linalg.norm(d), 1e-12)
        sag = float(np.max(np.abs((far - far[0]) @ nrm)))   # 弧到弦直线的最大垂距
        chord = float(np.linalg.norm(d))

        # 统一放大到远端附近（约 70 m 见方）
        cx, cy = far[:, 0].mean(), far[:, 1].mean()
        SP = 36.0
        ax.set_xlim(cx - SP, cx + SP)
        ax.set_ylim(cy - SP, cy + SP)
        ax.set_aspect("equal")
        ax.grid(alpha=0.22)
        ax.set_xlabel("x / m")
        ax.set_ylabel("y / m")
        ax.legend(loc="upper left", fontsize=8.0)
        ax.set_title(f"{c['tag']}\n$S_1$=({S1[0]:.0f},{S1[1]:.0f}), θ1={c['th']:.0f}°\n"
                     f"{c['note']}\n弦长 {chord:.1f} m，弧偏离弦 {sag:.3f} m",
                     fontsize=10.5)

    fig.suptitle("远端截线的四种形态（各图统一放大到远端附近约 70 m 见方）："
                 "前两种是单条圆弧，看着像直线；第三种是「两条弧 + 角点」；"
                 "折线只来自把圆域近似成多边形", fontsize=13.5)
    fig.tight_layout(rect=[0, 0, 1, 0.86])
    out = os.path.join(HERE, "远端截线的四种形态.png")
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
