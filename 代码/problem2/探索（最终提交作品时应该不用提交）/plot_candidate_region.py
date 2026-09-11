"""第二问候选区域的几何示意图。

算例：S1 = (0,0)，示向度 theta1 = 0°（正东），目标区域半径 1800 m，
      接收半径下界 R_min = 1000 m，近场阈值 5 m。

输出：候选区域示意图.png（双面板）
  左：整体几何 —— 目标区域、视线、源距离可行段、先验位置 x0、
      候选区域（两条垂直线段）、两个"最坏探测圆"
  右：交会角为什么会有最坏值 —— 源在段中点时 gamma = 90°，
      在两端点时 gamma 偏离到 gamma_min 与其补角
"""

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

# ---------------- 参数 ----------------
ARENA = 1800.0
R_MIN = 1000.0
ALPHA = 1.0
T_NEAR = 5.0
T_RECV = 1500.0

S1 = np.array([0.0, 0.0])
TH1 = 0.0
u1 = np.array([math.cos(math.radians(TH1)), math.sin(math.radians(TH1))])
n1 = np.array([-u1[1], u1[0]])

# 视线穿出圆域的距离（S1 在原点时即 ARENA）
b = float(S1 @ u1)
t_arena = -b + math.sqrt(max(0.0, b * b + ARENA**2 - float(S1 @ S1)))
t_max = min(T_RECV, t_arena)
t_min = T_NEAR

Delta = (t_max - t_min) / 2.0
t_bar = (t_max + t_min) / 2.0
x0 = S1 + t_bar * u1

feasible = Delta <= R_MIN
d_max = math.sqrt(max(0.0, R_MIN**2 - Delta**2))
gamma_min = math.degrees(math.acos(min(1.0, Delta / R_MIN)))
GAMMA_REQ = 40.0
d_min = Delta * math.tan(math.radians(GAMMA_REQ))

print(f"t_arena={t_arena:.4f}  t_max={t_max:.4f}  t_min={t_min:.4f}")
print(f"Delta={Delta:.4f}  t_bar={t_bar:.4f}  x0={x0}")
print(f"d_max={d_max:.4f}  gamma_min={gamma_min:.4f} deg  d_min={d_min:.4f}")
for s in (+1, -1):
    P = x0 + s * d_max * n1
    print(f"S2({'+-'[s<0]}) = ({P[0]:.4f}, {P[1]:.4f})  |P|={np.linalg.norm(P):.4f}")

# ---------------- 画图 ----------------
fig, (axL, axR) = plt.subplots(1, 2, figsize=(16.5, 8.4))

# ============ 左：整体几何 ============
th = np.linspace(0, 2 * math.pi, 600)
axL.plot(ARENA * np.cos(th), ARENA * np.sin(th), color="0.55", ls="--", lw=1.4,
         label=f"目标区域（半径 {ARENA:.0f} m）", zorder=1)

# 视线
L = t_max * 1.28
axL.plot([S1[0], S1[0] + L * u1[0]], [S1[1], S1[1] + L * u1[1]],
         color="tab:blue", lw=1.4, zorder=2,
         label="第一条视线（示向度 $\\theta_1$）")

# ±alpha 角楔（放大显示）
AX = 9.0
for sgn in (-1, 1):
    a = math.radians(TH1 + sgn * AX)
    axL.plot([S1[0], S1[0] + L * math.cos(a)], [S1[1], S1[1] + L * math.sin(a)],
             color="tab:blue", lw=0.9, alpha=0.45, zorder=2)
axL.fill_between([S1[0], S1[0] + L * u1[0]],
                 [S1[1], S1[1] + L * math.tan(math.radians(-AX))],
                 [S1[1], S1[1] + L * math.tan(math.radians(+AX))],
                 color="tab:blue", alpha=0.07, zorder=1)
axL.text(1450, -360,
         f"$\\pm{ALPHA:.0f}^\\circ$ 角楔（为看清已放大）",
         color="tab:blue", fontsize=9.5)

# 源距离可行段
P_lo = S1 + t_min * u1
P_hi = S1 + t_max * u1
axL.plot([P_lo[0], P_hi[0]], [P_lo[1], P_hi[1]],
         color="tab:green", lw=6, alpha=0.75, solid_capstyle="butt", zorder=4,
         label=f"源距离可行段 $t\\in[{t_min:.0f},\\,{t_max:.0f}]$ m")
for P in (P_lo, P_hi):
    axL.plot(*P, "o", color="tab:green", ms=7, zorder=6)
axL.annotate(f"$t_{{max}}={t_max:.0f}$", P_hi, textcoords="offset points",
             xytext=(-10, 30), color="tab:green", fontsize=10, ha="right")
axL.annotate(f"$t_{{min}}={t_min:.0f}$", P_lo, textcoords="offset points",
             xytext=(6, 34), color="tab:green", fontsize=10, ha="left")

# 两个"最坏探测圆"：以可行段两端为心、R_min 为半径
for P in (P_lo, P_hi):
    c = plt.Circle(P, R_MIN, color="0.45", fill=False, ls=":", lw=1.2, zorder=2)
    axL.add_patch(c)
axL.plot([], [], color="0.45", ls=":", lw=1.2,
         label=f"最坏探测圆（半径 $R_{{min}}={R_MIN:.0f}$ m）")

# 过 x0 的垂线
axL.plot([x0[0], x0[0]], [-860, 860], color="tab:orange", ls="-.", lw=1.3, zorder=3,
         label="过 $x_0$ 且垂直视线（短轴方向）")

# 候选区域
for s in (-1, 1):
    A = x0 + s * d_min * n1
    B = x0 + s * d_max * n1
    axL.plot([A[0], B[0]], [A[1], B[1]], color="tab:red", lw=7,
             solid_capstyle="round", zorder=7,
             label="候选区域" if s > 0 else None)
    axL.plot(*B, "o", color="tab:red", ms=9, zorder=8)

axL.plot(*x0, "s", color="k", ms=8, zorder=8, label="$x_0$（源的先验估计位置 $\\bar t$）")
axL.plot(*S1, "o", color="k", ms=9, zorder=8)
axL.annotate("$S_1$", S1, textcoords="offset points", xytext=(-30, -40), fontsize=12)
axL.annotate("$x_0$", x0, textcoords="offset points", xytext=(-26, -38), fontsize=12)
Pn = x0 + d_max * n1
axL.annotate(f"$S_2^+$ ({Pn[0]:.0f}, {Pn[1]:.0f})\n$d^*={d_max:.0f}$ m",
             Pn, textcoords="offset points", xytext=(16, 8),
             color="tab:red", fontsize=10)
Pn2 = x0 - d_max * n1
axL.annotate("$S_2^-$", Pn2, textcoords="offset points", xytext=(16, -16),
             color="tab:red", fontsize=11)

# 探测保证的交集说明
axL.annotate(f"两圆交点即候选点\n$\\gamma_{{min}}=\\arccos(\\Delta/R_{{min}})={gamma_min:.1f}^\\circ$",
             xy=(x0[0], d_max), xytext=(-1850, 1120), fontsize=11, color="tab:red",
             arrowprops=dict(arrowstyle="->", color="tab:red", lw=1.2))

axL.set_aspect("equal")
axL.set_xlim(-1950, 2050)
axL.set_ylim(-1950, 1950)
axL.grid(alpha=0.22)
axL.set_xlabel("x / m（东）")
axL.set_ylabel("y / m（北）")
axL.set_title("(a) 整体几何：一次观测后的可行域与第二检测点候选区域", fontsize=12)
axL.legend(loc="lower left", fontsize=8.5, framealpha=0.92)

# ============ 右：交会角 ============
G_mid = np.array([t_bar, 0.0])
G_lo = np.array([t_min, 0.0])
G_hi = np.array([t_max, 0.0])
S2p = x0 + d_max * n1

axR.plot([t_min - 120, t_max + 120], [0, 0], color="tab:green", lw=5,
         alpha=0.6, solid_capstyle="butt", zorder=3)
axR.plot([t_min, t_max], [0, 0], color="tab:green", lw=6, zorder=4)
for P, lab, off in ((G_lo, "$G(t_{min})$", (34, -34)),
                    (G_mid, "$G(\\bar t)$", (0, -40)),
                    (G_hi, "$G(t_{max})$", (30, -34))):
    axR.plot(*P, "o", color="tab:green", ms=9, zorder=6)
    axR.annotate(lab, P, textcoords="offset points", xytext=off,
                 color="tab:green", fontsize=11, ha="center")
axR.plot(*S1, "o", color="k", ms=9, zorder=6)
axR.annotate("$S_1$", S1, textcoords="offset points", xytext=(-30, 14), fontsize=12)
axR.plot(*S2p, "o", color="tab:red", ms=11, zorder=7)
axR.annotate("$S_2^+$", S2p, textcoords="offset points", xytext=(10, 6),
             color="tab:red", fontsize=12)

for P in (G_lo, G_mid, G_hi):
    axR.plot([P[0], S2p[0]], [P[1], S2p[1]], color="tab:red", lw=1.4, ls="--",
             alpha=0.8, zorder=5)
    axR.plot([P[0], S1[0]], [P[1], S1[1]], color="tab:blue", lw=1.4, ls="--",
             alpha=0.8, zorder=5)

# 角弧
def draw_angle(ax, V, A, B, r, color, label, lab_off=(0, 0)):
    a1 = math.atan2(A[1] - V[1], A[0] - V[0])
    a2 = math.atan2(B[1] - V[1], B[0] - V[0])
    d = (a2 - a1 + math.pi) % (2 * math.pi) - math.pi
    ts = np.linspace(0, d, 80)
    ax.plot(V[0] + r * np.cos(a1 + ts), V[1] + r * np.sin(a1 + ts),
            color=color, lw=2.0, zorder=8)
    am = a1 + d / 2
    ax.annotate(label, (V[0] + r * 1.55 * math.cos(am), V[1] + r * 1.55 * math.sin(am)),
                color=color, fontsize=11, ha="center", va="center",
                textcoords="offset points", xytext=lab_off, zorder=9)

draw_angle(axR, G_mid, S1, S2p, 250, "tab:purple", "$\\gamma=90^\\circ$")
draw_angle(axR, G_hi, S1, S2p, 250, "tab:orange",
           f"$\\gamma_{{min}}={gamma_min:.1f}^\\circ$", lab_off=(0, -14))
draw_angle(axR, G_lo, S1, S2p, 250, "tab:brown",
           f"$180^\\circ-\\gamma_{{min}}={180-gamma_min:.1f}^\\circ$", lab_off=(0, 12))

axR.set_aspect("equal")
axR.set_xlim(-560, 2020)
axR.set_ylim(-340, 1200)
axR.grid(alpha=0.22)
axR.set_xlabel("x / m（东）")
axR.set_ylabel("y / m（北）")
axR.set_title("(b) 交会角随源位置变化：中点处 90°，两端点最坏", fontsize=12)

fig.suptitle(
    f"第二问候选区域几何 —— $S_1$=({S1[0]:.0f},{S1[1]:.0f}), "
    f"$\\theta_1$={TH1:.0f}°, $\\Delta$={Delta:.1f} m, "
    f"$d^*$={d_max:.1f} m, $\\gamma_{{min}}$={gamma_min:.1f}°",
    fontsize=13.5)
fig.tight_layout(rect=[0, 0, 1, 0.955])

out = os.path.join(HERE, "候选区域示意图.png")
fig.savefig(out, dpi=155)
plt.close(fig)
print("->", out)
