"""Instance readers and generators.

* read_solomon_vrptw: classical VRPTW special case (one type, distance objective,
  no multiple trips), Euclidean distances truncated to one decimal (DIMACS rule).
* derive_mt_hfvrptw: deterministic MT-HFVRPTW instance derived from a Solomon file
  (rules of paper [9], subsection 5.7).
* generate: own instance generator (classes R, C, RC; sizes 25..200) with the
  parameters of paper [9], subsection 5.1.
"""
import math
import random
import os

from .core import Instance

DATA = os.path.join(os.path.dirname(__file__), "..", "data", "pyvrp_instances", "VRPTW", "Solomon")

# Vehicle types of paper [9], table 1
TYPES = [
    dict(name="truck", Q=220.0, W=48.0, F=190.0, c=1.45, speed=0.85, e=1.00, sig=1.40),
    dict(name="van", Q=110.0, W=30.0, F=105.0, c=1.00, speed=1.00, e=0.55, sig=1.00),
    dict(name="bike", Q=32.0, W=11.0, F=38.0, c=0.45, speed=0.70, e=0.00, sig=0.70),
]
ZONE_TRUCK = 1.55   # arcs inside the restricted zone are longer for trucks
BIKE_FACTOR = 0.88  # cargo bikes use links closed to motor traffic


def _parse_vrplib(path):
    sec = None
    meta = {}
    coords, dem, tw = {}, {}, {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line == "EOF":
                continue
            if line.endswith("SECTION"):
                sec = line
                continue
            if ":" in line and sec is None:
                k, v = line.split(":", 1)
                meta[k.strip()] = v.strip()
                continue
            p = line.split()
            if sec == "NODE_COORD_SECTION":
                coords[int(p[0])] = (float(p[1]), float(p[2]))
            elif sec == "DEMAND_SECTION":
                dem[int(p[0])] = float(p[1])
            elif sec == "TIME_WINDOW_SECTION":
                tw[int(p[0])] = (float(p[1]), float(p[2]))
    ids = sorted(coords)
    xy = [coords[i] for i in ids]
    q = [dem[i] for i in ids]
    a = [tw[i][0] for i in ids]
    b = [tw[i][1] for i in ids]
    st = float(meta.get("SERVICE_TIME", 0))
    s = [0.0] + [st] * (len(ids) - 1)
    return meta, xy, q, a, b, s


def solomon_bks(name):
    with open(os.path.join(DATA, name + ".sol")) as f:
        for line in f:
            if line.startswith("Cost"):
                return float(line.split()[1])


def solomon_names():
    return sorted(f[:-4] for f in os.listdir(DATA) if f.endswith(".vrp"))


def read_solomon_vrptw(name, n=None):
    meta, xy, q, a, b, s = _parse_vrplib(os.path.join(DATA, name + ".vrp"))
    if n is not None:
        xy, q, a, b, s = xy[:n + 1], q[:n + 1], a[:n + 1], b[:n + 1], s[:n + 1]
    Q = float(meta["CAPACITY"])
    types = [dict(name="vehicle", Q=Q, W=float("inf"), F=0.0, c=1.0, speed=1.0, e=0.0, sig=1.0)]
    vol = [0.0] * len(xy)
    return Instance("%s_%d" % (name, len(xy) - 1), xy, q, vol, a, b, s, types, H=b[0],
                    rho=0.0, multitrip=False, truncate=True)


def _profiles(xy, center, radius):
    """Distance matrices per type: trucks pay x1.55 on arcs whose midpoint lies in
    the restricted zone; bikes x0.88 on all arcs; vans the plain Euclidean metric."""
    N = len(xy)
    E = [[math.hypot(xy[i][0] - xy[j][0], xy[i][1] - xy[j][1]) for j in range(N)] for i in range(N)]
    Dt = [[0.0] * N for _ in range(N)]
    Db = [[0.0] * N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            mx = (xy[i][0] + xy[j][0]) / 2.0
            my = (xy[i][1] + xy[j][1]) / 2.0
            inside = math.hypot(mx - center[0], my - center[1]) <= radius
            Dt[i][j] = E[i][j] * (ZONE_TRUCK if inside else 1.0)
            Db[i][j] = E[i][j] * BIKE_FACTOR
    return [Dt, E, Db]


def derive_mt_hfvrptw(name, n, P=True, O=True, M=True):
    """Deterministic MT-HFVRPTW instance from a Solomon file (paper [9], 5.7).
    P/O/M switch on profiles, type-dependent service time, separate volume."""
    meta, xy, q, a, b, s = _parse_vrplib(os.path.join(DATA, name + ".vrp"))
    xy, q, a, b, s = xy[:n + 1], q[:n + 1], a[:n + 1], b[:n + 1], s[:n + 1]
    Qf = float(meta["CAPACITY"])
    vol = [0.0] + [q[i] * (0.15 + 0.10 * ((7 * i) % 5)) for i in range(1, n + 1)]
    mass = [1.00, 0.50, 0.15]
    volf = [0.24, 0.15, 0.055]
    types = []
    for k, t in enumerate(TYPES):
        d = dict(t)
        d["Q"] = mass[k] * Qf
        d["W"] = volf[k] * Qf
        types.append(d)
    xs = [p[0] for p in xy]
    ys = [p[1] for p in xy]
    center = ((min(xs) + max(xs)) / 2.0, (min(ys) + max(ys)) / 2.0)
    radius = 0.25 * min(max(xs) - min(xs), max(ys) - min(ys))
    return _build(name + "_%d" % n, xy, q, vol, a, b, s, types, b[0], center, radius, P, O, M)


def _build(name, xy, q, vol, a, b, s, types, H, center, radius, P=True, O=True, M=True):
    types = [dict(t) for t in types]
    if P:
        D = _profiles(xy, center, radius)
    else:
        D = None
    if not O:
        for t in types:
            t["sig"] = 1.0
    if not M:
        for t in types:
            t["W"] = float("inf")
    inst = Instance(name, xy, q, vol, a, b, s, types, H=H, rho=15.0, multitrip=True, D=D)
    inst.meta = dict(center=center, radius=radius, P=P, O=O, M=M)
    return inst


def generate(cls, n, idx, P=True, O=True, M=True):
    """Own generator (paper [9], 5.1): coordinates in [0,100]^2, depot in the centre,
    72 % small orders, four cargo densities, time windows around a reference route
    (width 120 for R/RC, 60 for C), horizon 720, reload time 15."""
    seed = {"R": 1, "C": 2, "RC": 3}[cls] * 100000 + n * 10 + idx
    rng = random.Random(seed)
    depot = (50.0, 50.0)
    pts = []
    if cls in ("C", "RC"):
        nc = 5 if n <= 50 else 8
        centers = [(rng.uniform(10, 90), rng.uniform(10, 90)) for _ in range(nc)]
    ncl = n if cls == "C" else (n // 2 if cls == "RC" else 0)
    for i in range(n):
        if i < ncl:
            cx, cy = centers[rng.randrange(len(centers))]
            pts.append((min(100, max(0, rng.gauss(cx, 5.0))), min(100, max(0, rng.gauss(cy, 5.0)))))
        else:
            pts.append((rng.uniform(0, 100), rng.uniform(0, 100)))
    xy = [depot] + pts
    q = [0.0]
    vol = [0.0]
    dens = [0.15, 0.25, 0.35, 0.50]
    for i in range(n):
        if rng.random() < 0.72:
            m = float(rng.randint(2, 9))
        else:
            m = float(rng.randint(18, 45))
        q.append(m)
        vol.append(m * dens[rng.randrange(4)])
    s = [0.0] + [10.0] * n
    H = 720.0
    width = 60.0 if cls == "C" else 120.0
    # reference route: nearest-neighbour tour at van speed, scaled to the horizon
    unv = set(range(1, n + 1))
    cur, t = 0, 0.0
    arr = [0.0] * (n + 1)
    while unv:
        nxt = min(unv, key=lambda j: math.hypot(xy[cur][0] - xy[j][0], xy[cur][1] - xy[j][1]))
        t += math.hypot(xy[cur][0] - xy[nxt][0], xy[cur][1] - xy[nxt][1]) + 10.0
        arr[nxt] = t
        unv.discard(nxt)
        cur = nxt
    scale = (H * 0.8) / t
    a = [0.0]
    b = [H]
    for i in range(1, n + 1):
        d0 = math.hypot(xy[i][0] - 50, xy[i][1] - 50)
        c = arr[i] * scale + rng.uniform(-60, 60)
        lo = max(d0 / 0.70 + 1.0, c - width / 2.0)
        hi = lo + width
        latest = H - 15.0 - 10.0 * 1.4 - d0 / 0.70 * 1.6
        if hi > latest:
            hi = latest
            lo = max(d0 / 0.70 + 1.0, hi - width)
        a.append(round(lo, 1))
        b.append(round(hi, 1))
    types = [dict(t) for t in TYPES]
    return _build("%s%d_%d" % (cls, n, idx), xy, q, vol, a, b, s, types, H, (50.0, 50.0), 25.0, P, O, M)


def own_set(sizes=(25, 50, 100)):
    return [(c, n, i) for n in sizes for c in ("R", "C", "RC") for i in (0, 1)]


def load_json(path):
    """Read an instance exported by export_instances.py (dataset/*/*.json)."""
    import json
    with open(path) as f:
        d = json.load(f)
    nodes = sorted(d["nodes"], key=lambda r: r["id"])
    xy = [(r["x"], r["y"]) for r in nodes]
    q = [r["mass"] for r in nodes]
    vol = [r["volume"] for r in nodes]
    a = [r["ready"] for r in nodes]
    b = [r["due"] for r in nodes]
    s = [r["service"] for r in nodes]
    types = [dict(name=t["name"], Q=t["mass_cap"], W=t["volume_cap"], F=t["fixed_cost"],
                  c=t["cost_per_dist"], speed=t["speed"], e=t["emission_per_dist"],
                  sig=t["service_mult"]) for t in d["vehicle_types"]]
    z = d["zone"]
    return _build(d["name"], xy, q, vol, a, b, s, types, d["horizon"],
                  tuple(z["center"]), z["radius"])
