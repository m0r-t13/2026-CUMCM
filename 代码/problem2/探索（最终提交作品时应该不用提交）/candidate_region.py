"""第二问：第二个检测点的选择算法（修正版）+ 候选区域可视化。

修正要点（相对最初版本）
------------------------
1. 探测半径下界是 max(r1, R_min) 而非恒为 R_min；
   对 t > R_min 的约束被 t = R_min 那条蕴含，
   故"保证集"是  D(G(t_min), R_min) ∩ D(G(t~), R_min)，  t~ = min(t_max, R_min)。
   —— 两个**等半径**圆盘，圆心都在视线上（不是 G(t_max)！）。
2. 目标函数最坏点：gamma(t) 关于 t 单调递减，
   偏离 90 度的程度 = arctan( M(a) / |b| )，  M(a) = max(|a-t_min|, |a-t_max|)。
3. 中垂线点 (a=(t_min+t_max)/2) 只在 t_max <= R_min 时严格最优；
   t_max > R_min 时最优 a 会偏离中垂线，需要一维搜索。
4. 追加约束 |S2| <= 1800（检测点须在目标区域内）。

输出：三场景对比图 + 场景 A 的约束/目标曲线图。
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

# ---------------- 常数 ----------------
ARENA = 1800.0
R_MIN = 1000.0
T_NEAR = 5.0
T_RECV = 1500.0
ALPHA = 1.0


# ============================================================
# 算法
# ============================================================
def solve(S1, theta1_deg, arena=ARENA, r_min=R_MIN,
          t_recv=T_RECV, t_near=T_NEAR, gamma_req=None,
          n_scan=4001):
    """返回结果字典。

    S1          : (x, y)
    theta1_deg  : 示向度（度）
    gamma_req   : 要求的交会角下限；None 表示取最优 gamma_min
    """
    S1 = np.asarray(S1, dtype=float)
    th = math.radians(theta1_deg)
    u1 = np.array([math.cos(th), math.sin(th)])
    n1 = np.array([-u1[1], u1[0]])

    # ① 视线穿出目标圆域的距离
    b0 = float(S1 @ u1)
    disc = b0 * b0 + arena * arena - float(S1 @ S1)
    t_arena = -b0 + math.sqrt(max(0.0, disc))

    # ② 可行段
    t_min = t_near
    t_max = min(t_recv, t_arena)
    feasible = t_max > t_min

    # ③ 关键：透镜远端用 t~ = min(t_max, r_min)
    t_til = min(t_max, r_min)

    c1 = float(S1 @ u1)      # S1 在 u1 方向的分量
    c2 = float(S1 @ n1)      # S1 在 n1 方向的分量

    def b_lens(a):
        """探测保证（透镜）给出的 |b| 上界。"""
        m = max((a - t_min) ** 2, (a - t_til) ** 2)
        v = r_min * r_min - m
        return math.sqrt(v) if v > 0 else -1.0

    def b_arena(a):
        """|S2| <= arena 给出的 |b| 上界（取正侧）。"""
        K = float(S1 @ S1) + a * a + 2 * a * c1 - arena * arena
        d = c2 * c2 - K
        if d < 0:
            return -1.0
        return -c2 + math.sqrt(d)

    def b_max(a):
        bl, ba = b_lens(a), b_arena(a)
        if bl <= 0 or ba <= 0:
            return 0.0
        return min(bl, ba)

    def M(a):
        return max(abs(a - t_min), abs(a - t_max))

    def gamma_min_of(a):
        bm = b_max(a)
        if bm <= 1e-9:
            return -1.0
        return 90.0 - math.degrees(math.atan2(M(a), bm))

    # ④ 一维搜索最优 a
    a_lo = max(t_til - r_min, -arena) + 1e-6
    a_hi = t_min + r_min - 1e-6
    aa = np.linspace(a_lo, a_hi, n_scan)
    gg = np.array([gamma_min_of(float(a)) for a in aa])
    best_i = int(np.argmax(gg))
    a_star = float(aa[best_i])
    gamma_opt = float(gg[best_i])
    b_star = b_max(a_star)

    # 解析解对照（仅 t_max <= r_min 时严格成立）
    if t_max <= r_min:
        a_ana = 0.5 * (t_min + t_max)
        g_ana = gamma_min_of(a_ana)
    else:
        a_ana = 0.5 * (t_min + t_max)     # 中垂线点
        g_ana = gamma_min_of(a_ana)

    # ⑤ 候选区域用的阈值
    greq = gamma_opt if gamma_req is None else gamma_req

    res = dict(
        S1=S1, u1=u1, n1=n1, theta1=theta1_deg,
        t_arena=t_arena, t_min=t_min, t_max=t_max, t_til=t_til,
        feasible=feasible,
        a_star=a_star, b_star=b_star, gamma_opt=gamma_opt,
        a_perp=a_ana, gamma_perp=g_ana,
        gamma_req=greq,
        S2p=S1 + a_star * u1 + b_star * n1,
        S2m=S1 + a_star * u1 - b_star * n1,
        b_max=b_max, M=M, gamma_min_of=gamma_min_of,
    )
    return res


def region_grid(r, n=900):
    """在 (a, b) 网格上求：可行集 / 优质集。"""
    a_lo = max(r["t_til"] - R_MIN, -ARENA)
    a_hi = r["t_min"] + R_MIN
    b_hi = R_MIN
    a = np.linspace(a_lo, a_hi, n)
    b = np.linspace(-b_hi, b_hi, n)
    A, B = np.meshgrid(a, b)

    # 探测保证（透镜）
    lens = ((A - r["t_min"]) ** 2 + B**2 <= R_MIN**2) & \
           ((A - r["t_til"]) ** 2 + B**2 <= R_MIN**2)
    # 目标圆域
    X = r["S1"][0] + A * r["u1"][0] + B * r["n1"][0]
    Y = r["S1"][1] + A * r["u1"][1] + B * r["n1"][1]
    inside = X**2 + Y**2 <= ARENA**2
    adm = lens & inside

    with np.errstate(divide="ignore", invalid="ignore"):
        Mv = np.maximum(np.abs(A - r["t_min"]), np.abs(A - r["t_max"]))
        gam = 90.0 - np.degrees(np.arctan2(Mv, np.abs(B)))
    good = adm & (gam >= r["gamma_req"])
    return A, B, X, Y, adm, good, gam


# ============================================================
# 绘图
# ============================================================
SCEN = [
    dict(tag="A", S1=(0.0, 0.0), th=0.0, note="视线几乎横穿圆域"),
    dict(tag="B", S1=(0.0, 1100.0), th=90.0, note="视线朝外，圆域在 700 m 截断"),
    dict(tag="C", S1=(0.0, 1500.0), th=90.0, note="贴近边界，仅 300 m 可用"),
]


def draw_scene(ax, r, title, levels):
    # 目标圆域
    th = np.linspace(0, 2 * math.pi, 600)
    ax.plot(ARENA * np.cos(th), ARENA * np.sin(th), color="0.45", ls="--", lw=1.3,
            label="目标区域 (r=1800 m)", zorder=5)

    # 视线与角楔边界（放大显示；不填充，避免冲淡色彩场）
    L = max(r["t_max"] * 1.35, 200.0)
    AX = 9.0
    for sgn in (-1, 1):
        al = math.radians(r["theta1"] + sgn * AX)
        ax.plot([r["S1"][0], r["S1"][0] + L * math.cos(al)],
                [r["S1"][1], r["S1"][1] + L * math.sin(al)],
                color="tab:blue", lw=1.0, alpha=0.75, zorder=9,
                label=(f"$\\pm{ALPHA:.0f}^\\circ$ 角楔（放大到 $\\pm{AX:.0f}^\\circ$ 显示）"
                       if sgn > 0 else None))
    ax.plot([r["S1"][0], r["S1"][0] + L * r["u1"][0]],
            [r["S1"][1], r["S1"][1] + L * r["u1"][1]],
            color="tab:blue", lw=1.4, alpha=0.9, zorder=9)

    # 可行源段
    P_lo = r["S1"] + r["t_min"] * r["u1"]
    P_hi = r["S1"] + r["t_max"] * r["u1"]
    ax.plot([P_lo[0], P_hi[0]], [P_lo[1], P_hi[1]], color="lime", lw=6,
            alpha=0.95, solid_capstyle="butt", zorder=10,
            label=f"源可行段 [{r['t_min']:.0f}, {r['t_max']:.0f}] m")
    for P in (P_lo, P_hi):
        ax.plot(*P, "o", color="lime", ms=8, mec="k", mew=0.6, zorder=11)

    # 两个等半径约束圆
    for P in (P_lo, r["S1"] + r["t_til"] * r["u1"]):
        ax.add_patch(plt.Circle(P, R_MIN, color="0.15", fill=False, ls=":",
                                lw=1.3, zorder=8))
        ax.plot(*P, "x", color="0.15", ms=8, mew=1.8, zorder=11)
    ax.plot([], [], color="0.15", ls=":", lw=1.3,
            label=f"探测保证圆 (R={R_MIN:.0f} m)")

    # —— 候选区域：以 gamma_min 着色的色彩场 ——
    _, _, X, Y, adm, good, gam = region_grid(r, n=560)
    G = np.where(adm, gam, np.nan)
    ax.contourf(X, Y, G, levels=levels, cmap="turbo", vmin=0, vmax=90,
                alpha=0.96, zorder=1)
    # 候选区边界（加粗）
    ax.contour(X, Y, adm.astype(float), levels=[0.5], colors="k",
               linewidths=1.8, zorder=7)
    # 等值线
    try:
        cs = ax.contour(X, Y, G, levels=[20, 30, 40, 50, 60, 70, 80],
                        colors="k", linewidths=0.7, alpha=0.5, zorder=6)
        ax.clabel(cs, fmt="%d°", fontsize=8)
    except Exception:
        pass

    # 最优两点
    for S2, lab in ((r["S2p"], "$S_2^+$"), (r["S2m"], "$S_2^-$")):
        ax.plot(*S2, "*", color="w", ms=27, zorder=13)
        ax.plot(*S2, "*", color="red", ms=20, zorder=14)
        ax.annotate(lab, S2, textcoords="offset points", xytext=(12, 8),
                    color="red", fontsize=13, zorder=15)

    ax.plot(*r["S1"], "o", color="k", ms=9, zorder=13)
    ax.annotate("$S_1$", r["S1"], textcoords="offset points", xytext=(-28, -34),
                fontsize=13, zorder=15)

    lim = ARENA * 1.06
    ax.set_aspect("equal")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.grid(alpha=0.18)
    ax.set_xlabel("x / m")
    ax.set_ylabel("y / m")
    ax.set_title(title, fontsize=11.5)


def main():
    outs = []
    print("=" * 96)
    print("第二问：第二检测点选择（修正版）")
    print("=" * 96)
    hdr = (f"{'场景':<4}{'S1':<16}{'θ1':>6}{'t_arena':>9}{'t_max':>8}{'t~':>7}"
           f"{'a*(最优)':>10}{'a(中垂线)':>11}{'b*':>8}{'γ_opt':>8}{'γ_中垂线':>9}")
    print(hdr)
    print("-" * 96)

    results = []
    for sc in SCEN:
        r = solve(sc["S1"], sc["th"])
        results.append((sc, r))
        print(f"{sc['tag']:<4}{str(tuple(int(v) for v in sc['S1'])):<16}"
              f"{sc['th']:>6.0f}{r['t_arena']:>9.1f}{r['t_max']:>8.1f}"
              f"{r['t_til']:>7.0f}{r['a_star']:>10.1f}{r['a_perp']:>11.1f}"
              f"{r['b_star']:>8.1f}{r['gamma_opt']:>7.2f}°{r['gamma_perp']:>8.2f}°")
        print(f"     S2+ = ({r['S2p'][0]:8.1f}, {r['S2p'][1]:8.1f})   "
              f"|S2+| = {np.linalg.norm(r['S2p']):7.1f} m   "
              f"({'在圆域内' if np.linalg.norm(r['S2p']) <= ARENA else '出圆域'})")

    # ---------------- 图 1：三场景对比 ----------------
    levels = np.linspace(0, 90, 46)
    fig, axes = plt.subplots(1, 3, figsize=(21, 7.6))
    for ax, (sc, r) in zip(axes, results):
        draw_scene(
            ax, r,
            f"({sc['tag']}) $S_1$=({sc['S1'][0]:.0f},{sc['S1'][1]:.0f}), "
            f"$\\theta_1$={sc['th']:.0f}°\n"
            f"$t_{{arena}}$={r['t_arena']:.0f} → $t_{{max}}$={r['t_max']:.0f}, "
            f"$\\tilde t$={r['t_til']:.0f} → $\\gamma_{{\\min}}$={r['gamma_opt']:.1f}°"
            f"   （{sc['note']}）", levels)
    axes[0].legend(loc="lower left", fontsize=8.5, framealpha=0.93)
    sm = plt.cm.ScalarMappable(cmap="turbo", norm=plt.Normalize(0, 90))
    sm.set_array([])
    fig.subplots_adjust(left=0.035, right=0.905, top=0.80, bottom=0.07, wspace=0.13)
    cax = fig.add_axes([0.917, 0.12, 0.011, 0.66])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("该位置能做到的最坏交会角 $\\gamma_{min}$ / 度", fontsize=11)
    fig.suptitle("第二问候选区域：彩色区 = 可保证探测的位置（黑色粗线为边界），颜色 = 能做到的最坏交会角\n"
                 "★ = 最优点；视线在目标圆域内还能走多远，决定了定位能有多准",
                 fontsize=14, y=0.965)
    p1 = os.path.join(HERE, "候选区域_三场景对比.png")
    fig.savefig(p1, dpi=140)
    plt.close(fig)
    outs.append(p1)
    print("\n->", p1)

    # ---------------- 图 2：场景 A 的约束与目标曲线 ----------------
    scA, rA = results[0]
    fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6.4))

    aa = np.linspace(max(rA["t_til"] - R_MIN, -ARENA) + 1e-6,
                     rA["t_min"] + R_MIN - 1e-6, 3000)
    bl = np.array([rA["b_max"](float(a)) for a in aa])
    blens = np.array([rA["b_max"](float(a)) for a in aa])
    # 单独算透镜与圆域两条界
    def b_lens_only(a):
        m = max((a - rA["t_min"]) ** 2, (a - rA["t_til"]) ** 2)
        v = R_MIN * R_MIN - m
        return math.sqrt(v) if v > 0 else np.nan
    def b_arena_only(a):
        c1 = float(rA["S1"] @ rA["u1"]); c2 = float(rA["S1"] @ rA["n1"])
        K = float(rA["S1"] @ rA["S1"]) + a * a + 2 * a * c1 - ARENA * ARENA
        d = c2 * c2 - K
        return (-c2 + math.sqrt(d)) if d >= 0 else np.nan
    bL = np.array([b_lens_only(float(a)) for a in aa])
    bA = np.array([b_arena_only(float(a)) for a in aa])

    ax1.plot(aa, bL, color="tab:blue", lw=2.2, label="透镜约束 (探测保证)")
    ax1.plot(aa, bA, color="tab:purple", lw=2.2, label="目标圆域约束 $|S_2|\\leq 1800$")
    ax1.plot(aa, blens, color="k", lw=1.2, ls="--", label="$b_{max}=\\min$（实际）")
    ax1.axvline(rA["a_star"], color="red", lw=2, ls="-.",
                label=f"最优 $a^*$={rA['a_star']:.0f}")
    ax1.axvline(rA["a_perp"], color="tab:green", lw=1.6, ls=":",
                label=f"中垂线点 {rA['a_perp']:.0f}")
    ax1.set_xlabel("$a$（沿视线方向的偏移）/ m")
    ax1.set_ylabel("$|b|$ 上界 / m")
    ax1.set_ylim(0, 1900)
    ax1.grid(alpha=0.25)
    ax1.legend(fontsize=9, loc="upper right")
    ax1.annotate("圆域约束始终在上方\n（场景 A 中不起作用）",
                 (aa[len(aa) // 6], bA[len(aa) // 6]),
                 textcoords="offset points", xytext=(10, 10),
                 color="tab:purple", fontsize=9)
    ax1.set_title("(a) 候选点距视线能偏多远——两条约束谁在起作用", fontsize=12)

    gg = np.array([rA["gamma_min_of"](float(a)) for a in aa])
    ax2.plot(aa, gg, color="tab:red", lw=2.4)
    ax2.plot(rA["a_star"], rA["gamma_opt"], "*", color="red", ms=22,
             label=f"最优 $\\gamma_{{min}}$={rA['gamma_opt']:.2f}°")
    ax2.plot(rA["a_perp"], rA["gamma_perp"], "o", color="tab:green", ms=9,
             label=f"中垂线点 $\\gamma$={rA['gamma_perp']:.2f}°")
    ax2.axhline(90, color="0.5", ls=":", lw=1)
    ax2.annotate("$90^\\circ$（源的先验位置处）", (aa[0], 90),
                 textcoords="offset points", xytext=(8, 6), color="0.4", fontsize=10)
    ax2.set_xlabel("$a$（沿视线方向的偏移）/ m")
    ax2.set_ylabel("最坏交会角 $\\gamma_{min}$ / 度")
    ax2.grid(alpha=0.25)
    ax2.legend(fontsize=10)
    ax2.set_title("(b) 最坏交会角随 $a$ 变化——最优点略偏中垂线", fontsize=12)

    fig2.suptitle(f"场景 A（$S_1$=(0,0), $\\theta_1$=0°）：候选点是怎么被两条约束夹出来的",
                  fontsize=14)
    fig2.tight_layout(rect=[0, 0, 1, 0.94])
    p2 = os.path.join(HERE, "候选区域_约束与目标曲线.png")
    fig2.savefig(p2, dpi=140)
    plt.close(fig2)
    outs.append(p2)
    print("->", p2)

    return outs


if __name__ == "__main__":
    import sys
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    main()
