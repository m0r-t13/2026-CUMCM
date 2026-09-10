"""B1: bearing-wedge intersection, polygon diameter and circle coverage.

Dependencies: numpy, scipy. Angles are degrees counterclockwise from +x.
Run: python b1_coverage.py --self-test
     python b1_coverage.py --demo --trials 400
     python b1_coverage.py --input observations.json --output result.json
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
import unittest

import numpy as np
from scipy.optimize import linprog, minimize

TOL = 1e-7  # absolute distance tolerance, metres (problem coordinates)


def cross(a, b):
    return float(a[0] * b[1] - a[1] * b[0])


def bearing_halfplanes(observations, error_deg=1.0):
    """Return unit-normal inequalities A @ point <= b for forward wedges."""
    if not observations or not 0 < error_deg < 90:
        raise ValueError('Require observations and 0 < error_deg < 90.')
    normals, offsets = [], []
    for obs in observations:
        s = np.asarray([obs['x'], obs['y']], dtype=float)
        angle = float(obs['bearing_deg'])
        if not np.isfinite(s).all() or not math.isfinite(angle):
            raise ValueError('All coordinates and bearings must be finite.')
        lo, hi = np.radians([angle - error_deg, angle + error_deg])
        # cross(u_lo, p-s) >= 0; cross(u_hi, p-s) <= 0.
        for normal in ([math.sin(lo), -math.cos(lo)],
                       [-math.sin(hi), math.cos(hi)]):
            normal = np.asarray(normal)
            normals.append(normal)
            offsets.append(float(normal @ s))
    return np.asarray(normals), np.asarray(offsets)


def convex_hull(points):
    """Andrew monotone chain, CCW order, excluding collinear interior points."""
    unique = []
    for p in sorted(map(tuple, points)):
        p = np.asarray(p, dtype=float)
        if not unique or np.linalg.norm(p - unique[-1]) > TOL:
            unique.append(p)
    if len(unique) < 3:
        return np.asarray(unique).reshape((-1, 2))
    chains = []
    for seq in (unique, unique[::-1]):
        chain = []
        for p in seq:
            while len(chain) >= 2 and cross(chain[-1]-chain[-2], p-chain[-1]) <= 0:
                chain.pop()
            chain.append(p)
        chains.append(chain[:-1])
    return np.asarray(chains[0] + chains[1])


def intersect_halfplanes(A, b):
    """Classify before enumerating vertices: never use an artificial box.

    For m constraints, intersection enumeration costs O(m^3); sorting costs
    O(m^2 log m). This transparent method suits a small number of stations.
    """
    lp_options = {'primal_feasibility_tolerance': 1e-9,
                  'dual_feasibility_tolerance': 1e-9}
    def solve(c):
        return linprog(c, A_ub=A, b_ub=b, bounds=[(None, None)]*2,
                       method='highs', options=lp_options)
    feasible = solve([0, 0])
    if feasible.status == 2:
        return 'empty', np.empty((0, 2))
    if not feasible.success:
        raise RuntimeError('Feasibility LP failed: ' + feasible.message)
    # A nonempty subset of R^2 is bounded iff both coordinates have finite bounds.
    for objective in ([1, 0], [-1, 0], [0, 1], [0, -1]):
        result = solve(objective)
        if result.status == 3:
            return 'unbounded', np.empty((0, 2))
        if not result.success:
            raise RuntimeError('Boundedness LP failed: ' + result.message)
    vertices = []
    for i, j in itertools.combinations(range(len(b)), 2):
        if abs(cross(A[i], A[j])) < 1e-13:
            continue
        point = np.linalg.solve(A[[i, j]], b[[i, j]])
        if np.max(A @ point - b) <= TOL:
            vertices.append(point)
    if not vertices:
        raise ArithmeticError('Bounded intersection but no stable vertices; rescale input.')
    return 'bounded', convex_hull(vertices)


def polygon_diameter(vertices):
    """Rotating calipers: O(h) for h convex vertices in CCW order."""
    p = np.asarray(vertices)
    n = len(p)
    if n == 0:
        raise ValueError('An empty set has no diameter pair here.')
    if n < 3:
        return float(np.linalg.norm(p[-1]-p[0])), (0, n-1)
    best, pair, j = -1.0, (0, 0), 1
    def consider(i, k):
        nonlocal best, pair
        d2 = float((p[i]-p[k]) @ (p[i]-p[k]))
        if d2 > best:
            best, pair = d2, (i, k)
    for i in range(n):
        ni = (i+1) % n
        edge = p[ni] - p[i]
        def area(k):
            return cross(edge, p[k]-p[i])
        for _ in range(n):
            nj = (j+1) % n
            if area(nj) <= area(j) + 1e-12:
                break
            j = nj
        consider(i, j)
        consider(ni, j)
        nj = (j+1) % n
        if abs(area(nj)-area(j)) <= 1e-10:
            consider(i, nj)
            consider(ni, nj)
    return math.sqrt(best), pair


def minimum_enclosing_circle(vertices):
    """Enumerate support sets of 1, 2 or 3 vertices; deterministic O(h^4).

    This is an independent small-polygon check, not needed for the yes/no test.
    A large-scale implementation can replace this with randomized incremental MEC.
    """
    points = np.asarray(vertices)
    best_radius, best_center, best_support = math.inf, None, None
    def consider(center, radius, support):
        nonlocal best_radius, best_center, best_support
        if radius >= best_radius:
            return
        actual = float(np.max(np.linalg.norm(points-center, axis=1)))
        if actual <= radius + TOL:
            best_radius, best_center, best_support = actual, center, support
    for i, p in enumerate(points):
        consider(p.copy(), 0.0, [i])
    for i, j in itertools.combinations(range(len(points)), 2):
        center = (points[i]+points[j])/2
        consider(center, float(np.linalg.norm(points[i]-center)), [i, j])
    for i, j, k in itertools.combinations(range(len(points)), 3):
        u, v = points[j]-points[i], points[k]-points[i]
        if abs(cross(u, v)) <= 1e-13 * max(1.0, np.linalg.norm(u)*np.linalg.norm(v)):
            continue
        delta = np.linalg.solve(2*np.array([u, v]), np.array([u@u, v@v]))
        consider(points[i]+delta, float(np.linalg.norm(delta)), [i, j, k])
    if best_center is None:
        raise ArithmeticError('No numerically stable enclosing circle found.')
    return best_center, best_radius, best_support


def analyze(observations, error_deg=1.0):
    A, b = bearing_halfplanes(observations, error_deg)
    status, vertices = intersect_halfplanes(A, b)
    out = {'status': status, 'observations': observations, 'error_deg': error_deg,
           'model': 'intersection of forward bearing wedges only', 'tolerance_m': TOL}
    if status != 'bounded':
        out['explanation'] = ('No common localization point.' if status == 'empty'
                              else 'Wedge intersection unbounded; no finite diameter.')
        return out
    diameter, pair = polygon_diameter(vertices)
    center = (vertices[pair[0]]+vertices[pair[1]])/2
    radius = diameter/2
    distances = np.linalg.norm(vertices-center, axis=1)
    worst = int(np.argmax(distances))
    mec_center, mec_radius, support = minimum_enclosing_circle(vertices)
    out.update(vertices=vertices.tolist(), vertex_count=len(vertices),
               diameter_m=diameter, diameter_pair_indices=list(pair),
               diameter_endpoints=vertices[list(pair)].tolist(),
               diameter_circle_center=center.tolist(), diameter_circle_radius_m=radius,
               covers=bool(np.max(distances) <= radius + TOL),
               max_distance_from_diameter_center_m=float(distances[worst]),
               maximum_radial_excess_m=float(distances[worst]-radius),
               witness=vertices[worst].tolist(),
               vertex_distances_m=distances.tolist(),
               mec_center=mec_center.tolist(), mec_radius_m=mec_radius,
               mec_support_indices=support, minimum_covering_diameter_m=2*mec_radius,
               radius_ratio=(mec_radius/radius if radius > TOL else 1.0))
    return out


def observations_for_source(stations, source=(0, 0), errors=None):
    if errors is None:
        errors = np.zeros(len(stations))
    return [{'x': float(x), 'y': float(y),
             'bearing_deg': (math.degrees(math.atan2(source[1]-y, source[0]-x))+float(e)) % 360}
            for (x, y), e in zip(stations, errors)]


def equilateral_example():
    # Each lower boundary is one triangle edge; each upper boundary is redundant.
    h = 20/math.sqrt(3)
    triangle = np.array([[0, h], [-10, -h/2], [10, -h/2]])
    observations = []
    for i in range(3):
        direction = (triangle[(i+1)%3]-triangle[i])/20
        s = triangle[i]-1000*direction
        observations.append({'x': float(s[0]), 'y': float(s[1]),
                             'bearing_deg': (math.degrees(math.atan2(direction[1], direction[0]))+1)%360})
    return observations


def verify_result(result):
    """Audit calipers against exhaustive pairs and both coverage criteria."""
    if result['status'] != 'bounded':
        return
    p = np.array(result['vertices'])
    brute = float(np.max(np.linalg.norm(p[:, None]-p[None, :], axis=2)))
    assert abs(brute-result['diameter_m']) <= 1e-6
    assert result['covers'] == (result['mec_radius_m'] <= brute/2+TOL)
    assert 1-1e-7 <= result['radius_ratio'] <= 2/math.sqrt(3)+1e-7
    A, b = bearing_halfplanes(result['observations'], result['error_deg'])
    assert np.max(A@p.T-b[:, None]) <= 2*TOL


def circle_svg(result, filename):
    """Standalone geometric diagram; all shapes use one isotropic scale."""
    p = np.array(result['vertices'])
    c, r = np.array(result['diameter_circle_center']), result['diameter_circle_radius_m']
    mc, mr = np.array(result['mec_center']), result['mec_radius_m']
    lo = np.min(np.vstack([p, c-r, mc-mr]), axis=0)
    hi = np.max(np.vstack([p, c+r, mc+mr]), axis=0)
    span = max(float(np.max(hi-lo)), 1.0)*1.25
    mid = (lo+hi)/2
    scale = 520/span
    def xy(pt):
        return 340+(pt[0]-mid[0])*scale, 330-(pt[1]-mid[1])*scale
    def circle(center, radius, color, dash=''):
        x, y = xy(center)
        return f'<circle cx="{x}" cy="{y}" r="{radius*scale}" fill="none" stroke="{color}" stroke-width="2.5" stroke-dasharray="{dash}"/>'
    shapes = [f'<svg xmlns="http://www.w3.org/2000/svg" width="920" height="650" viewBox="0 0 920 650"><rect width="920" height="650" fill="#fafbfe"/>',
              '<g font-family="Arial, sans-serif" fill="#182337">',
              '<text x="30" y="38" font-size="23">B1: Does the diameter circle cover the region?</text>',
              circle(mc, mr, '#13866b'), circle(c, r, '#c64449', '7 5')]
    polygon = ' '.join(f'{x:.8f},{y:.8f}' for x, y in map(xy, p))
    shapes.append(f'<polygon points="{polygon}" fill="#658fd0" fill-opacity="0.18" stroke="#336ca9" stroke-width="2"/>')
    a, b = [xy(pt) for pt in result['diameter_endpoints']]
    shapes.append(f'<line x1="{a[0]}" y1="{a[1]}" x2="{b[0]}" y2="{b[1]}" stroke="#c64449" stroke-width="2"/>')
    for i, pt in enumerate(p):
        x, y = xy(pt)
        shapes.append(f'<circle cx="{x}" cy="{y}" r="4" fill="#336ca9"/><text x="{x+7}" y="{y-8}" font-size="16">V{i}</text>')
    for center, color in [(c, '#c64449'), (mc, '#13866b')]:
        x, y = xy(center)
        shapes.append(f'<circle cx="{x}" cy="{y}" r="4" fill="{color}"/>')
    lines = [('Diameter circle', '#c64449'), (f'r = D/2 = {r:.6f} m', '#182337'),
             ('Minimum enclosing circle', '#13866b'), (f'R = {mr:.6f} m', '#182337'),
             (f'Covers region: {result["covers"]}', '#182337'),
             (f'Excess: {result["maximum_radial_excess_m"]:.6f} m', '#182337')]
    for i, (label, color) in enumerate(lines):
        shapes.append(f'<text x="645" y="{150+36*i}" fill="{color}" font-size="16">{label}</text>')
    shapes.append('</g></svg>')
    filename.write_text('\n'.join(shapes), encoding='utf-8')


def demo(outdir, trials):
    outdir.mkdir(parents=True, exist_ok=True)
    examples = {'equilateral_counterexample': equilateral_example(),
                'symmetric_cover': observations_for_source([(-1000, 0), (0, -1000)]),
                'unbounded': [{'x': 0, 'y': 0, 'bearing_deg': 0}],
                'inconsistent': [{'x': 0, 'y': 0, 'bearing_deg': 180},
                                 {'x': 100, 'y': 0, 'bearing_deg': 0}]}
    results = {name: analyze(obs) for name, obs in examples.items()}
    rng = np.random.default_rng(20260910)
    counts = {'bounded_covered': 0, 'bounded_not_covered': 0, 'unbounded': 0, 'empty': 0}
    by_n = {}
    optimization_checks = 0
    for t in range(trials):
        n = 2+t % 5
        angle = rng.uniform(0, 2*np.pi, n)
        distance = rng.uniform(200, 1000, n)
        stations = np.column_stack([distance*np.cos(angle), distance*np.sin(angle)])
        observations = observations_for_source(stations, errors=rng.uniform(-1, 1, n))
        result = analyze(observations)
        verify_result(result)
        A, b = bearing_halfplanes(observations)
        assert np.max(-b) <= TOL  # The simulated true source is feasible.
        label = ('bounded_covered' if result['covers'] else 'bounded_not_covered') if result['status']=='bounded' else result['status']
        counts[label] += 1
        by_n.setdefault(str(n), {k: 0 for k in counts})[label] += 1
        if n == 2 and label == 'bounded_not_covered' and 'two_station_counterexample' not in results:
            results['two_station_counterexample'] = result
        # Independent convex optimization audit of 20 minimum enclosing circles.
        if result['status'] == 'bounded' and optimization_checks < 20:
            p = np.array(result['vertices'])
            origin = p.mean(axis=0)
            normalization = max(result['diameter_m'], 1)
            q = (p-origin)/normalization
            opt = minimize(lambda z: z[2], [0, 0, 1.0], method='SLSQP',
                           bounds=[(None, None), (None, None), (0, None)],
                           constraints=[{'type': 'ineq', 'fun': lambda z: z[2]-np.linalg.norm(q-z[:2], axis=1)}],
                           options={'ftol': 1e-10, 'maxiter': 500})
            assert opt.success, opt.message
            assert abs(opt.x[2]*normalization-result['mec_radius_m']) < 1e-5
            optimization_checks += 1
    for name, result in results.items():
        verify_result(result)
        (outdir/f'{name}.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
        (outdir/f'{name}_input.json').write_text(json.dumps({'observations': result['observations']}, indent=2), encoding='utf-8')
        if result['status'] == 'bounded':
            circle_svg(result, outdir/f'{name}.svg')
    summary = {'seed': 20260910, 'trials': trials, 'counts': counts, 'by_station_count': by_n,
               'mec_independent_optimization_checks': optimization_checks,
               'sampling': 'Source (0,0); 2..6 stations; uniform azimuth, radius 200..1000 m; uniform error -1..1 deg.',
               'caution': 'Synthetic reproducible experiments, not simulator logs or a universal failure probability.'}
    (outdir/'test_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    for name, result in results.items():
        print(name, {k: result[k] for k in ['status', 'diameter_m', 'covers', 'mec_radius_m', 'maximum_radial_excess_m'] if k in result})


class GeometryTests(unittest.TestCase):
    def test_equilateral(self):
        r = analyze(equilateral_example())
        verify_result(r)
        self.assertFalse(r['covers'])
        self.assertEqual(r['vertex_count'], 3)
        self.assertAlmostEqual(r['diameter_m'], 20, places=7)
        self.assertAlmostEqual(r['mec_radius_m'], 20/math.sqrt(3), places=7)
        self.assertAlmostEqual(r['maximum_radial_excess_m'], 10*(math.sqrt(3)-1), places=7)

    def test_symmetric_cover(self):
        r = analyze(observations_for_source([(-1000, 0), (0, -1000)]))
        verify_result(r)
        self.assertTrue(r['covers'])

    def test_empty_and_unbounded(self):
        self.assertEqual(analyze([{'x': 0, 'y': 0, 'bearing_deg': 0}])['status'], 'unbounded')
        self.assertEqual(analyze([{'x': 0, 'y': 0, 'bearing_deg': 180},
                                 {'x': 100, 'y': 0, 'bearing_deg': 0}])['status'], 'empty')

    def test_degenerate(self):
        A = np.array([[1,0], [-1,0], [0,1], [0,-1]])
        for b, size, diameter in [([2,0,0,0], 2, 2), ([0,0,0,0], 1, 0)]:
            status, p = intersect_halfplanes(A, np.array(b))
            self.assertEqual(status, 'bounded')
            self.assertEqual(len(p), size)
            self.assertAlmostEqual(polygon_diameter(p)[0], diameter)
            self.assertAlmostEqual(minimum_enclosing_circle(p)[1], diameter/2)

    def test_angles_and_forward_direction(self):
        observations = [{'x': 0, 'y': 0, 'bearing_deg': 0}]
        A, b = bearing_halfplanes(observations)
        self.assertTrue(np.all(A@np.array([100, 0]) <= b+TOL))
        self.assertFalse(np.all(A@np.array([-100, 0]) <= b+TOL))
        C, d = bearing_halfplanes([{'x': 0, 'y': 0, 'bearing_deg': 360}])
        np.testing.assert_allclose(A, C, atol=1e-14)
        np.testing.assert_allclose(b, d)

    def test_calipers_against_brute_force(self):
        rng = np.random.default_rng(713)
        for _ in range(150):
            p = convex_hull(rng.normal(size=(60, 2)))
            brute = np.max(np.linalg.norm(p[:,None]-p[None,:], axis=2))
            self.assertAlmostEqual(polygon_diameter(p)[0], brute, places=10)
        square = np.array([[0,0], [1,0], [1,1], [0,1]])
        self.assertAlmostEqual(polygon_diameter(square)[0], math.sqrt(2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path)
    parser.add_argument('--output', type=Path, default=Path('result.json'))
    parser.add_argument('--demo', action='store_true')
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--trials', type=int, default=400)
    parser.add_argument('--outdir', type=Path, default=Path(__file__).parent/'results')
    args = parser.parse_args()
    if args.self_test:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(GeometryTests)
        if not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful():
            raise SystemExit(1)
    if args.demo:
        if args.trials < 0:
            parser.error('--trials must be nonnegative')
        demo(args.outdir, args.trials)
    if args.input:
        data = json.loads(args.input.read_text(encoding='utf-8-sig'))
        result = analyze(data['observations'], data.get('error_deg', 1.0))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
        print(json.dumps(result, ensure_ascii=False, indent=2))
    if not (args.self_test or args.demo or args.input):
        parser.print_help()


if __name__ == '__main__':
    main()
