"""Quality indicators and statistics."""
import math
from scipy.optimize import linprog
from scipy.stats import wilcoxon

from .core import dominates, nondominated


def bounds(point_sets):
    allp = [p for ps in point_sets for p in ps]
    ideal = tuple(min(p[j] for p in allp) for j in range(3))
    nadir = tuple(max(p[j] for p in allp) for j in range(3))
    return ideal, nadir


def normalise(points, ideal, nadir):
    out = []
    for p in points:
        out.append(tuple((p[j] - ideal[j]) / (nadir[j] - ideal[j]) if nadir[j] > ideal[j] else 0.0
                         for j in range(3)))
    return out


def _hv2(pts, rx, ry):
    pts = sorted(pts)
    area = 0.0
    best_y = ry
    for x, y in pts:
        if y < best_y:
            area += (rx - x) * (best_y - y)
            best_y = y
    return area


def hv3(points, ref=(1.1, 1.1, 1.1)):
    pts = [p for p in nondominated(points) if p[0] < ref[0] and p[1] < ref[1] and p[2] < ref[2]]
    if not pts:
        return 0.0
    pts.sort(key=lambda p: p[2])
    hv = 0.0
    for i, p in enumerate(pts):
        z_next = pts[i + 1][2] if i + 1 < len(pts) else ref[2]
        if z_next <= p[2]:
            continue
        hv += _hv2([(q[0], q[1]) for q in pts[:i + 1]], ref[0], ref[1]) * (z_next - p[2])
    return hv / ((ref[0]) * (ref[1]) * (ref[2]))


def coverage(A, B):
    """Share of points of B weakly dominated by some point of A."""
    if not B:
        return 0.0
    c = 0
    for b in B:
        if any(all(a[j] <= b[j] for j in range(3)) for a in A):
            c += 1
    return c / len(B)


def supported(points, eps=1e-6, tol=1e-9):
    """For each point of a non-dominated set: is it optimal for some weighted sum
    with strictly positive weights (w_j >= eps, sum w = 1)? Max-margin LP:
    maximise t subject to w.(q - p) >= t for all q != p; supported iff t >= -tol.
    The margin formulation makes the verdict independent of constraint order."""
    res = []
    for p in points:
        others = [q for q in points if q is not p]
        if not others:
            res.append(True)
            continue
        # variables w1, w2, w3, t ; minimise -t
        A_ub = [[-(q[j] - p[j]) for j in range(3)] + [1.0] for q in others]
        b_ub = [0.0] * len(others)
        r = linprog(c=[0, 0, 0, -1], A_ub=A_ub, b_ub=b_ub, A_eq=[[1, 1, 1, 0]], b_eq=[1],
                    bounds=[(eps, 1)] * 3 + [(None, None)], method="highs")
        res.append(r.status == 0 and -r.fun >= -tol)
    return res


def wil(x, y):
    d = [a - b for a, b in zip(x, y)]
    if all(abs(v) < 1e-12 for v in d):
        return 1.0
    return float(wilcoxon(x, y, alternative="two-sided", zero_method="wilcox").pvalue)


def holm(pvals):
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    prev = 0.0
    for r, i in enumerate(order):
        v = min(1.0, (m - r) * pvals[i])
        prev = max(prev, v)
        adj[i] = prev
    return adj
