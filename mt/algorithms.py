"""Search methods: MO-ILS (archive, cycled scalarisation directions) with ablation
variants, sequential weighted sum (WS), memetic NSGA-II, and single-objective ILS
for the classical VRPTW special case. All methods run under a wall-clock budget."""
import random
import time

from .core import (Sol, Archive, Counter, construct, local_search, destroy_repair,
                   gval, weight_vectors, dominates, crowding, eval_vehicle, veh_obj)

W15 = weight_vectors(4)
BAL = (1 / 3, 1 / 3, 1 / 3)


def scales(inst):
    """Characteristic objective scales s_j of an instance: objectives of a greedy
    cost-oriented construction (fixed seed, independent of the run seed)."""
    rng = random.Random(12345)
    s = construct(inst, rng, (1.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    o = s.objs()
    return (max(o[0], 1e-9), max(o[1], 1e-6 * max(o[0], 1.0)), max(o[2], 1e-9))


def _rand_probs(rng, K):
    x = [rng.random() + 0.05 for _ in range(K)]
    t = sum(x)
    return [v / t for v in x]


def mo_ils(inst, budget, seed, dirs="cycle", perturb="dr", strength="random",
           probe=False, sc=None):
    """Multi-objective ILS. Variants: dirs in {cycle, single}; perturb in {dr, restart};
    strength in {random, naive, cautious}. Returns dict with archive and statistics."""
    rng = random.Random(seed)
    t0 = time.perf_counter()
    now = time.perf_counter
    deadline = t0 + budget
    sc = sc or scales(inst)
    arch = Archive(50)
    cnt = Counter()
    W = W15 if dirs == "cycle" else [BAL]
    # initial solutions
    for j in range(3):
        w = W[(j * 5) % len(W)]
        x = construct(inst, rng, w, sc, _rand_probs(rng, inst.K))
        local_search(x, inst, w, sc, rng, deadline, now, cnt=cnt)
        arch.add(x)
    it = 0
    s_cur = 0.10        # naive rule
    ub = 0.35           # cautious rule upper bound
    stag = 0
    switches = improving = 0
    while now() < deadline:
        w = W[it % len(W)]
        if perturb == "restart":
            x = construct(inst, rng, w, sc, _rand_probs(rng, inst.K))
        else:
            items = arch.items
            if rng.random() < 0.6:
                base = min(items, key=lambda it_: gval(it_[0], w, sc))[1]
            else:
                base = items[rng.randrange(len(items))][1]
            if strength == "random":
                st = rng.uniform(0.10, 0.35)
            elif strength == "naive":
                st = s_cur
            else:
                st = rng.uniform(0.10, ub)
            x = destroy_repair(base, inst, w, sc, rng, st)
        local_search(x, inst, w, sc, rng, deadline, now, cnt=cnt)
        if now() > deadline and it > 0:
            # descent cut by the deadline: keep only if feasible improvement, still valid
            pass
        acc = arch.add(x)
        if probe and len(W) > 1 and now() < deadline:
            wn = W[(it + 1) % len(W)]
            switches += 1
            if local_search(x.copy(), inst, wn, sc, rng, deadline + 60, now, first_only=True):
                improving += 1
        if acc:
            stag = 0
            s_cur = 0.10
            ub = max(0.10, ub * 0.9)
        else:
            stag += 1
            if stag >= 8:
                stag = 0
                s_cur = min(0.60, s_cur * 1.3)
                ub = min(0.60, ub * 1.3)
        it += 1
    return dict(points=arch.points(), sols=[s for _, s in arch.items], iters=it,
                evals=cnt.evals, switches=switches, improving=improving, sc=sc)


def ws_seq(inst, budget, seed, sc=None):
    """Sequential weighted sum: the same local search and perturbation, budget split
    between the 15 fixed directions, each started from the previous direction's best."""
    rng = random.Random(seed)
    now = time.perf_counter
    t0 = now()
    sc = sc or scales(inst)
    cnt = Counter()
    slot = budget / len(W15)
    x = construct(inst, rng, W15[0], sc, _rand_probs(rng, inst.K))
    bests = []
    it = 0
    for d, w in enumerate(W15):
        dl = t0 + slot * (d + 1)
        local_search(x, inst, w, sc, rng, dl, now, cnt=cnt)
        best = x.copy()
        gb = gval(best.objs(), w, sc)
        while now() < dl:
            y = destroy_repair(best, inst, w, sc, rng, rng.uniform(0.10, 0.35))
            local_search(y, inst, w, sc, rng, dl, now, cnt=cnt)
            gy = gval(y.objs(), w, sc)
            it += 1
            if gy < gb - 1e-12:
                best, gb = y, gy
        bests.append(best)
        x = best.copy()
    arch = Archive(50)
    for s in bests:
        arch.add(s)
    return dict(points=arch.points(), sols=[s for _, s in arch.items], iters=it,
                evals=cnt.evals, sc=sc)


def _nd_sort(pts):
    n = len(pts)
    S = [[] for _ in range(n)]
    nd = [0] * n
    fronts = [[]]
    for p in range(n):
        for q in range(n):
            if p == q:
                continue
            if dominates(pts[p], pts[q]):
                S[p].append(q)
            elif dominates(pts[q], pts[p]):
                nd[p] += 1
        if nd[p] == 0:
            fronts[0].append(p)
    i = 0
    while fronts[i]:
        nxt = []
        for p in fronts[i]:
            for q in S[p]:
                nd[q] -= 1
                if nd[q] == 0:
                    nxt.append(q)
        i += 1
        fronts.append(nxt)
    return fronts[:-1]


def _route_crossover(p1, p2, inst, rng, sc, share=0.4):
    """Copy 40 % of p2's vehicles into p1, removing their customers elsewhere.
    Profile matrices need not satisfy the triangle inequality, so a vehicle that
    loses customers may become infeasible; its remaining customers are reinserted."""
    child = p1.copy()
    take = [v for v in p2.veh if rng.random() < share]
    if not take:
        take = [p2.veh[rng.randrange(len(p2.veh))]]
    moved = set(i for v in take for t in v[1] for i in t)
    nv, nc, orphans = [], [], []
    for k, trips in child.veh:
        nt = [[i for i in t if i not in moved] for t in trips]
        nt = [t for t in nt if t]
        if not nt:
            continue
        ev = eval_vehicle(inst, k, nt)
        if ev is None:
            orphans += [i for t in nt for i in t]
            continue
        nv.append([k, nt])
        nc.append(veh_obj(inst, k, ev))
    for k, trips in take:
        tr = [list(t) for t in trips]
        nv.append([k, tr])
        nc.append(veh_obj(inst, k, eval_vehicle(inst, k, tr)))
    child.veh, child.comp = nv, nc
    if orphans:
        from .core import insert_all
        rng.shuffle(orphans)
        insert_all(child, inst, orphans, W15[rng.randrange(len(W15))], sc, rng, 0.0)
    return child


def nsga2(inst, budget, seed, pop=20, sc=None):
    """Memetic NSGA-II with route crossover and limited local search on offspring."""
    rng = random.Random(seed)
    now = time.perf_counter
    t0 = now()
    deadline = t0 + budget
    sc = sc or scales(inst)
    cnt = Counter()
    P = []
    for j in range(pop):
        if now() > deadline and len(P) >= 2:
            break
        w = W15[rng.randrange(len(W15))]
        x = construct(inst, rng, w, sc, _rand_probs(rng, inst.K))
        local_search(x, inst, w, sc, rng, deadline, now, cnt=cnt)
        P.append(x)
    gens = 0

    def rank_crowd(pop_):
        pts = [s.objs() for s in pop_]
        fr = _nd_sort(pts)
        rank = [0] * len(pop_)
        cd = [0.0] * len(pop_)
        for r, f in enumerate(fr):
            c = crowding([pts[i] for i in f])
            for i, ci in zip(f, c):
                rank[i] = r
                cd[i] = ci
        return rank, cd

    while now() < deadline:
        rank, cd = rank_crowd(P)

        def tour():
            i, j = rng.randrange(len(P)), rng.randrange(len(P))
            if (rank[i], -cd[i]) < (rank[j], -cd[j]):
                return P[i]
            return P[j]

        off = []
        for _ in range(pop):
            if now() > deadline:
                break
            c = _route_crossover(tour(), tour(), inst, rng, sc)
            w = W15[rng.randrange(len(W15))]
            local_search(c, inst, w, sc, rng, deadline, now, cnt=cnt)
            off.append(c)
        if now() > deadline and len(off) < pop:
            # an incomplete generation is discarded
            break
        R = P + off
        rank, cd = rank_crowd(R)
        order = sorted(range(len(R)), key=lambda i: (rank[i], -cd[i]))
        P = [R[i] for i in order[:pop]]
        gens += 1
    rank, _ = rank_crowd(P)
    arch = Archive(10 ** 6)
    for s, r in zip(P, rank):
        if r == 0:
            arch.add(s)
    return dict(points=arch.points(), sols=[s for _, s in arch.items], iters=gens,
                evals=cnt.evals, sc=sc)


def ils_single(inst, budget, seed):
    """Single-objective ILS (distance) for the VRPTW special case."""
    rng = random.Random(seed)
    now = time.perf_counter
    deadline = now() + budget
    w = (1.0, 0.0, 0.0)
    sc = (1.0, 1.0, 1.0)
    cur = construct(inst, rng, w, sc)
    local_search(cur, inst, w, sc, rng, deadline, now)
    best = cur.copy()
    gc = gb = cur.objs()[0]
    it = 0
    while now() < deadline:
        y = destroy_repair(cur, inst, w, sc, rng, rng.uniform(0.10, 0.35))
        local_search(y, inst, w, sc, rng, deadline, now)
        gy = y.objs()[0]
        it += 1
        if gy < gc - 1e-9 or (gy <= gc * 1.02 and rng.random() < 0.05):
            cur, gc = y, gy
        if gy < gb - 1e-9:
            best, gb = y.copy(), gy
    return dict(best=best, cost=gb, iters=it)
