"""问题 1 前半问端到端求解：示向度 → 定位区域 → 旋转卡壳求直径。

用法
----
    python solve_problem1_diameter.py --list
    python solve_problem1_diameter.py --demo all
    python solve_problem1_diameter.py --demo case1
    python solve_problem1_diameter.py --json mycase.json
    python solve_problem1_diameter.py --benchmark          # 卡壳 vs 暴力 的规模对比

自定义算例 JSON：
    {
      "name": "my case",
      "alpha_deg": 1.0,
      "detections": [{"x": 1350, "y": -250, "bearing_deg": 180.0}, ...]
    }
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys

import localization_region as lr
import rotating_calipers as rc

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
OUTDIR = os.path.join(HERE, "results")
ARENA_RADIUS = 1800.0


# 求解


def solve(detections, alpha_deg: float = 1.0, name: str = "case", **kw):
    """返回一个结果字典。detections = [(x, y, bearing_deg), ...]"""
    status, poly, notes = lr.region_from_bearings(detections, alpha_deg, **kw)
    res = dict(
        name=name,
        status=status,
        alpha_deg=alpha_deg,
        detections=[list(d) for d in detections],
        notes=notes,
        region=[list(p) for p in poly],
        n_vertices=len(poly),
        area=rc.polygon_area(poly) if poly else 0.0,
        perimeter=rc.polygon_perimeter(poly) if poly else 0.0,
        min_edge=rc.min_edge_length(poly) if poly else 0.0,
        diameter=0.0,
        diameter_pair=None,
        diameter_pair_brute=None,
        calipers_vs_brute=0.0,
        antipodal_pairs=[],
        n_antipodal=0,
        tied_pairs=[],
    )
    if status != "ok" or len(poly) < 2:
        return res

    d_cal, i_cal, j_cal = rc.rotating_calipers_diameter(poly)
    d_bf, i_bf, j_bf = rc.brute_force_diameter(poly)
    ap = rc.antipodal_pairs(poly)

    res.update(
        diameter=d_bf,
        diameter_pair=[int(i_cal), int(j_cal)],
        diameter_pair_brute=[int(i_bf), int(j_bf)],
        calipers_vs_brute=abs(d_cal - d_bf),
        antipodal_pairs=[[int(a), int(b)] for a, b in ap],
        n_antipodal=len(ap),
        tied_pairs=[[int(a), int(b)] for a, b in rc.diameter_pairs(poly)],
    )
    return res


def report(res) -> None:
    print("=" * 76)
    print(f"算例：{res['name']}")
    print("=" * 76)
    print(f"检测点数 k = {len(res['detections'])}，角楔半角 α = {res['alpha_deg']:g}°")
    for i, (x, y, th) in enumerate(res["detections"], 1):
        print(f"    S{i}: ({x:12.4f}, {y:12.4f})   示向度 = {th:.4f}°")
    for n in res["notes"]:
        print(f"  [约束] {n}")

    if res["status"] != "ok":
        print(f"\n区域状态：{res['status']} —— 无法给出有限直径")
        return

    print(f"\n定位区域：{res['n_vertices']} 顶点，面积 {res['area']:.4f} m²，"
          f"周长 {res['perimeter']:.4f} m，最短边 {res['min_edge']:.6f} m")
    for i, (x, y) in enumerate(res["region"]):
        print(f"    V{i}: ({x:14.6f}, {y:14.6f})")

    print(f"\n【旋转卡壳】")
    i, j = res["diameter_pair"]
    print(f"    直径 D = {res['diameter']:.9f} m")
    print(f"    直径端点 V{i} ---- V{j}")
    print(f"    枚举到的对踵点对 {res['n_antipodal']} 个"
          f"（上界 3n = {3 * res['n_vertices']}），"
          f"顶点对总数 C(n,2) = {res['n_vertices'] * (res['n_vertices'] - 1) // 2}")
    print(f"    对踵点对：{res['antipodal_pairs']}")
    print(f"\n【交叉校验】暴力枚举 O(n²) 给出 "
          f"D = {res['diameter']:.9f} m，端点 {res['diameter_pair_brute']}")
    print(f"    两种方法偏差 = {res['calipers_vs_brute']:.3e} m")
    if len(res["tied_pairs"]) > 1:
        print(f"    达到直径的顶点对不唯一：{res['tied_pairs']}")


# 内置算例


def dfs(G, pts):
    return [
        math.degrees(math.atan2(G[1] - py, G[0] - px)) % 360.0 for px, py in pts
    ]


def demo_cases():
    cases = {}

    def make(key, name, G, pts, **kw):
        bs = dfs(G, pts)
        cases[key] = dict(
            name=name,
            alpha_deg=kw.pop("alpha_deg", 1.0),
            detections=[[px, py, b] for (px, py), b in zip(pts, bs)],
            true_source=list(G),
            **kw,
        )

    G = (150.0, -250.0)

    make("case1", "case1 双站正交交会（交会角 90°）", G,
         [(G[0] + 1200.0, G[1]), (G[0], G[1] + 1200.0)])

    make("case2", "case2 三站交会", G,
         [(G[0] + 1200.0, G[1]),
          (G[0], G[1] + 1200.0),
          (G[0] + 1300.0 * math.cos(math.radians(210.0)),
           G[1] + 1300.0 * math.sin(math.radians(210.0)))])

    make("case3", "case3 双站小交会角 5°（细长区域）", G,
         [(G[0] + 1200.0, G[1]),
          (G[0] + 1300.0 * math.cos(math.radians(5.0)),
           G[1] + 1300.0 * math.sin(math.radians(5.0)))])

    make("case4", "case4 双站同方位（方向锥重叠 → 区域无界）", G,
         [(G[0] + 1200.0, G[1]), (G[0] + 1300.0, G[1])])

    cases["case4b"] = dict(cases["case4"])
    cases["case4b"]["name"] = "case4b 同方位 + 1800 m 圆域约束（变为有界）"
    cases["case4b"]["region_disk_radius"] = ARENA_RADIUS

    make("case5", "case5 四站交会（比较对踵点对数量）", G,
         [(G[0] + 1200.0, G[1]),
          (G[0], G[1] + 1200.0),
          (G[0] + 1300.0 * math.cos(math.radians(210.0)),
           G[1] + 1300.0 * math.sin(math.radians(210.0))),
          (G[0] + 1400.0 * math.cos(math.radians(300.0)),
           G[1] + 1400.0 * math.sin(math.radians(300.0)))])

    # 多站算例：同一组 12 个检测点，按固定顺序取前 k 个。
    # 这样 k 增大时约束单调增多，定位区域**嵌套收缩**，直径单调不增，
    # 可直接看出"多站交会提高精度"的规律，并逼近 n = 2k 的顶点数上界。
    # 方向角经挑选：任意两方向的夹角 >= 15°（区域有界），且不存在恰好
    # 相差 180° 的一对（避免出现平行的边界直线）。
    multi = [
        (20.0, 1200.0), (120.0, 1350.0), (225.0, 1100.0), (310.0, 1450.0),
        (55.0, 1250.0), (140.0, 1300.0), (275.0, 1150.0), (35.0, 1400.0),
        (190.0, 1050.0), (105.0, 1420.0), (290.0, 1180.0), (205.0, 1330.0),
    ]

    def multi_case(key, name, k):
        pts = []
        for ang, rr in multi[:k]:
            a = math.radians(ang)
            pts.append((G[0] + rr * math.cos(a), G[1] + rr * math.sin(a)))
        bs = dfs(G, pts)
        cases[key] = dict(
            name=name,
            alpha_deg=1.0,
            detections=[[px, py, b] for (px, py), b in zip(pts, bs)],
            true_source=list(G),
        )

    for k, key in ((6, "case6"), (8, "case7"), (12, "case8")):
        multi_case(key, f"{key} {k}站交会", k)

    return cases


def nested_table():
    """同一组检测点按顺序取前 k 个，观察区域与直径随 k 的变化。"""
    G = (150.0, -250.0)
    multi = [
        (20.0, 1200.0), (120.0, 1350.0), (225.0, 1100.0), (310.0, 1450.0),
        (55.0, 1250.0), (140.0, 1300.0), (275.0, 1150.0), (35.0, 1400.0),
        (190.0, 1050.0), (105.0, 1420.0), (290.0, 1180.0), (205.0, 1330.0),
    ]
    print("=" * 82)
    print("多站交会：区域随检测点数 k 的嵌套收缩（真源 G=(150, -250)）")
    print("=" * 82)
    print(f"{'k':>3} {'n':>3} {'2k':>4} {'面积(m²)':>14} {'直径 D(m)':>14} "
          f"{'D 降幅':>10} {'对踵点对':>8} {'C(n,2)':>7}")
    prev = None
    for k in range(2, 13):
        pts = []
        for ang, rr in multi[:k]:
            a = math.radians(ang)
            pts.append((G[0] + rr * math.cos(a), G[1] + rr * math.sin(a)))
        bs = dfs(G, pts)
        dets = [(px, py, b) for (px, py), b in zip(pts, bs)]
        res = solve(dets, 1.0, name=f"k={k}")
        if res["status"] != "ok":
            print(f"{k:>3}  区域{res['status']}")
            continue
        drop = "" if prev is None else f"{(prev - res['diameter']) / prev:>9.2%}"
        n = res["n_vertices"]
        print(f"{k:>3} {n:>3} {2 * k:>4} {res['area']:>14.4f} "
              f"{res['diameter']:>14.6f} {drop:>10} {res['n_antipodal']:>8} "
              f"{n * (n - 1) // 2:>7}")
        prev = res["diameter"]
    print("\n注：n 为定位区域顶点数，上界为 2k；D 随 k 单调不增（区域嵌套收缩）。")


def run_case(case, save=True):
    dets = [tuple(d) for d in case["detections"]]
    kw = {}
    for k in ("region_disk_radius", "reception_max"):
        if case.get(k):
            kw[k] = case[k]
    res = solve(dets, case.get("alpha_deg", 1.0), name=case["name"], **kw)
    report(res)
    if save:
        os.makedirs(OUTDIR, exist_ok=True)
        key = case.get("key", "case")
        with open(os.path.join(OUTDIR, f"{key}.json"), "w", encoding="utf-8") as fh:
            json.dump(res, fh, ensure_ascii=False, indent=2)
        try:
            plot(res, os.path.join(FIGDIR, f"{key}.png"))
        except Exception as exc:  # noqa: BLE001
            print(f"  [warn] 绘图失败：{exc}")
    print()
    return res


# 绘图
def plot(res, path):
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

    fig, ax = plt.subplots(figsize=(7.6, 7.6))
    th = [2 * math.pi * i / 400 for i in range(401)]
    ax.plot([ARENA_RADIUS * math.cos(t) for t in th],
            [ARENA_RADIUS * math.sin(t) for t in th],
            color="0.75", ls="--", lw=1.0, label=f"目标区域 r={ARENA_RADIUS:.0f} m")

    reach = ARENA_RADIUS * 2
    for idx, (x, y, b) in enumerate(res["detections"], 1):
        for sgn in (-1, 1):
            a = math.radians(b + sgn * res["alpha_deg"])
            ax.plot([x, x + reach * math.cos(a)], [y, y + reach * math.sin(a)],
                    color="tab:blue", lw=0.7, alpha=0.5)
        a = math.radians(b)
        ax.annotate("", xy=(x + 650 * math.cos(a), y + 650 * math.sin(a)),
                    xytext=(x, y),
                    arrowprops=dict(arrowstyle="->", color="tab:blue", lw=1.1))
        ax.plot(x, y, "o", color="tab:blue", ms=6)
        ax.annotate(f"S{idx}", (x, y), textcoords="offset points",
                    xytext=(7, 6), color="tab:blue")

    if res["status"] == "ok" and res["n_vertices"] >= 3:
        reg = res["region"]
        ax.fill([p[0] for p in reg], [p[1] for p in reg],
                color="tab:red", alpha=0.32, zorder=3,
                label=f"定位区域（D={res['diameter']:.3f} m, n={res['n_vertices']}）")
        ax.plot([p[0] for p in reg] + [reg[0][0]],
                [p[1] for p in reg] + [reg[0][1]],
                color="tab:red", lw=2.0, zorder=4)
        i, j = res["diameter_pair"]
        A, B = reg[i], reg[j]
        ax.plot([A[0], B[0]], [A[1], B[1]], color="tab:green", lw=1.8, zorder=6,
                label=f"直径 V{i}--V{j}")
        for P in (A, B):
            ax.plot(P[0], P[1], "s", color="tab:green", ms=6, zorder=7)
        # 标出全部对踵点对的连线（浅色）
        for (p, q) in res["antipodal_pairs"]:
            if (p, q) == (i, j):
                continue
            ax.plot([reg[p][0], reg[q][0]], [reg[p][1], reg[q][1]],
                    color="tab:orange", lw=0.6, alpha=0.6, zorder=5)
        cx = sum(p[0] for p in reg) / len(reg)
        cy = sum(p[1] for p in reg) / len(reg)
        span = max(res["diameter"] * 2.6, 150.0)
        ax.set_xlim(cx - span, cx + span)
        ax.set_ylim(cy - span, cy + span)
    else:
        ax.set_xlim(-ARENA_RADIUS * 1.1, ARENA_RADIUS * 1.1)
        ax.set_ylim(-ARENA_RADIUS * 1.1, ARENA_RADIUS * 1.1)

    ax.set_aspect("equal")
    ax.grid(alpha=0.25)
    ax.set_xlabel("x / m（东）")
    ax.set_ylabel("y / m（北）")
    ax.set_title(f"旋转卡壳求定位区域直径 —— {res['name']}")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=170)
    plt.close(fig)
    print(f"  [图] {path}")


# 规模对比：卡壳 O(n) vs 暴力 O(n^2)


def benchmark(seed: int = 20260910):
    """比较旋转卡壳与暴力枚举：对踵点对数量、耗时随顶点数 n 的增长。

    用圆周附近的点构造凸包，以便把 n 做到几百，看清 O(n) 与 O(n^2) 的差别。
    """
    import random
    import time

    rnd = random.Random(seed)
    print("=" * 76)
    print("旋转卡壳 O(n) vs 暴力枚举 O(n^2)")
    print("=" * 76)
    print(f"{'n':>6} {'C(n,2)':>10} {'对踵点对':>10} {'卡壳(ms)':>10} "
          f"{'暴力(ms)':>10} {'比值':>8} {'最大偏差':>12}")
    for N in (16, 32, 64, 128, 256, 512, 1024, 2048):
        circ = [
            (3000.0 * math.cos(2 * math.pi * k / N) * (1 + rnd.uniform(-1e-3, 1e-3)),
             3000.0 * math.sin(2 * math.pi * k / N) * (1 + rnd.uniform(-1e-3, 1e-3)))
            for k in range(N)
        ]
        poly = rc.convex_hull(circ)
        nv = len(poly)
        if nv < 3:
            continue
        t0 = time.perf_counter()
        d1, _, _ = rc.rotating_calipers_diameter(poly)
        t1 = time.perf_counter()
        d2, _, _ = rc.brute_force_diameter(poly)
        t2 = time.perf_counter()
        ap = rc.antipodal_pairs(poly)
        tc, tb = (t1 - t0) * 1e3, (t2 - t1) * 1e3
        print(f"{nv:>6} {nv * (nv - 1) // 2:>10} {len(ap):>10} "
              f"{tc:>10.3f} {tb:>10.3f} {tb / tc if tc > 0 else 0:>8.1f} "
              f"{abs(d1 - d2):>12.3e}")
    print("\n注：本题定位区域的顶点数 n <= 2k（k 为检测点数），通常是个位数，")
    print("    因此实际计算用暴力枚举即可；旋转卡壳的价值在定理层面"
          "（保证只需考察 O(n) 个对踵点对），而非本题的运行速度。")


# CLI


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="问题 1 前半问：旋转卡壳求定位区域直径")
    ap.add_argument("--demo", default=None)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--json", default=None)
    ap.add_argument("--benchmark", action="store_true")
    ap.add_argument("--nested", action="store_true",
                    help="多站交会：区域随 k 的嵌套收缩表")
    args = ap.parse_args(argv)

    cases = demo_cases()
    if args.list:
        for k, v in cases.items():
            print(f"{k}\t{v['name']}")
        return 0
    if args.nested:
        nested_table()
        return 0
    if args.benchmark:
        benchmark()
        return 0
    if args.json:
        with open(args.json, encoding="utf-8") as fh:
            case = json.load(fh)
        case.setdefault("key", os.path.splitext(os.path.basename(args.json))[0])
        run_case(case)
        return 0
    if args.demo:
        if args.demo == "all":
            sel = list(cases)
        elif args.demo in cases:
            sel = [args.demo]
        else:
            sel = [k for k in cases if args.demo in k]
        if not sel:
            print(f"未找到算例 {args.demo!r}")
            return 2
        for k in sel:
            c = dict(cases[k])
            c["key"] = k
            run_case(c)
        return 0

    ap.print_help()
    return 0


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    raise SystemExit(main())
