"""Core model, evaluation and search operators for the multi-trip heterogeneous
fleet vehicle routing problem with time windows (MT-HFVRPTW).

Solution representation: list of vehicles; vehicle = [k, trips], trips = list of
lists of customer indices (depot = 0 is implicit at both ends of every trip).
"""
import math
import random

INF = float("inf")


class Instance:
    """Problem data. Distances/travel times may differ per vehicle type (profiles)."""

    def __init__(self, name, xy, q, vol, a, b, s, types, H, rho=15.0,
                 multitrip=True, D=None, truncate=False, n_nbr=14):
        self.name = name
        self.xy = xy
        self.n = len(xy) - 1
        self.q, self.vol, self.a, self.b, self.s = q, vol, a, b, s
        self.types = types
        self.K = len(types)
        self.H = H
        self.rho = rho
        self.multitrip = multitrip
        N = self.n + 1
        E = [[0.0] * N for _ in range(N)]
        for i in range(N):
            xi, yi = xy[i]
            for j in range(N):
                d = math.hypot(xi - xy[j][0], yi - xy[j][1])
                if truncate:
                    d = math.floor(d * 10.0) / 10.0
                E[i][j] = d
        self.E = E
        self.D = D if D is not None else [E for _ in types]
        self.T = []
        self.S = []
        for k, t in enumerate(types):
            sp = t["speed"]
            Dk = self.D[k]
            self.T.append([[Dk[i][j] / sp for j in range(N)] for i in range(N)])
            self.S.append([s[i] * t["sig"] for i in range(N)])
        self.nbr = [[]]
        for i in range(1, N):
            order = sorted(range(1, N), key=lambda j: E[i][j])
            self.nbr.append([j for j in order if j != i][:n_nbr])


def eval_vehicle(inst, k, trips):
    """Return (distance, end_time) or None if infeasible."""
    T = inst.T[k]
    D = inst.D[k]
    S = inst.S[k]
    a, b, q, vol = inst.a, inst.b, inst.q, inst.vol
    tp = inst.types[k]
    Q, W = tp["Q"], tp["W"]
    t = a[0]
    L = 0.0
    first = True
    for trip in trips:
        if not first:
            t += inst.rho
        first = False
        m = 0.0
        v = 0.0
        for i in trip:
            m += q[i]
            v += vol[i]
        if m > Q + 1e-9 or v > W + 1e-9:
            return None
        prev = 0
        for i in trip:
            L += D[prev][i]
            t += T[prev][i]
            if t > b[i] + 1e-9:
                return None
            if t < a[i]:
                t = a[i]
            t += S[i]
            prev = i
        L += D[prev][0]
        t += T[prev][0]
        if t > inst.H + 1e-9:
            return None
    return L, t


def veh_obj(inst, k, ev):
    tp = inst.types[k]
    L, end = ev
    return (tp["F"] + tp["c"] * L, tp["e"] * L, end)


class Sol:
    """Solution with cached per-vehicle objective components."""

    __slots__ = ("veh", "comp")

    def __init__(self, veh=None, comp=None):
        self.veh = veh if veh is not None else []
        self.comp = comp if comp is not None else []

    def copy(self):
        return Sol([[v[0], [list(t) for t in v[1]]] for v in self.veh], list(self.comp))

    def objs(self):
        c = e = 0.0
        m = 0.0
        for cc in self.comp:
            c += cc[0]
            e += cc[1]
            if cc[2] > m:
                m = cc[2]
        return (c, e, m)

    def recompute(self, inst):
        self.comp = []
        for k, trips in self.veh:
            ev = eval_vehicle(inst, k, trips)
            if ev is None:
                raise ValueError("infeasible vehicle in recompute")
            self.comp.append(veh_obj(inst, k, ev))

    def clean(self):
        nv, nc = [], []
        for v, c in zip(self.veh, self.comp):
            v[1] = [t for t in v[1] if t]
            if v[1]:
                nv.append(v)
                nc.append(c)
        self.veh, self.comp = nv, nc


def gval(o, w, sc):
    return w[0] * o[0] / sc[0] + w[1] * o[1] / sc[1] + w[2] * o[2] / sc[2]


def dominates(a, b):
    return a[0] <= b[0] and a[1] <= b[1] and a[2] <= b[2] and (a[0] < b[0] or a[1] < b[1] or a[2] < b[2])


# ---------------------------------------------------------------- insertion

def _delta_with(sol, inst, w, sc, changes, newveh=None):
    """g of solution after replacing vehicle comps in `changes` {idx: comp}
    and appending `newveh` comp (or None). Returns g."""
    c = e = 0.0
    m = 0.0
    for i, cc in enumerate(sol.comp):
        if i in changes:
            cc = changes[i]
            if cc is None:
                continue
        c += cc[0]
        e += cc[1]
        if cc[2] > m:
            m = cc[2]
    if newveh is not None:
        c += newveh[0]
        e += newveh[1]
        if newveh[2] > m:
            m = newveh[2]
    return w[0] * c / sc[0] + w[1] * e / sc[1] + w[2] * m / sc[2]


class Totals:
    """Fast incremental g evaluation: sums plus sorted ends for max."""

    def __init__(self, sol):
        self.C = sum(c[0] for c in sol.comp)
        self.Em = sum(c[1] for c in sol.comp)
        self.ends = [c[2] for c in sol.comp]

    def max_excl(self, excl):
        m = 0.0
        for i, x in enumerate(self.ends):
            if i not in excl and x > m:
                m = x
        return m


def g_after(sol, tot, w, sc, old_idx, new_comps):
    """g after replacing vehicles old_idx (tuple) by new_comps (list of comps)."""
    C = tot.C
    Em = tot.Em
    for i in old_idx:
        C -= sol.comp[i][0]
        Em -= sol.comp[i][1]
    m = tot.max_excl(old_idx)
    for cc in new_comps:
        C += cc[0]
        Em += cc[1]
        if cc[2] > m:
            m = cc[2]
    return w[0] * C / sc[0] + w[1] * Em / sc[1] + w[2] * m / sc[2]


def best_insertion(sol, inst, u, w, sc, rng, noise=0.0, full=False, type_probs=None):
    """Find (noisy) cheapest feasible insertion of customer u. Returns move tuple."""
    tot = Totals(sol)
    g0 = gval(sol.objs(), w, sc)
    best = None
    bestv = INF
    if full:
        cand_trips = [(vi, ti) for vi, v in enumerate(sol.veh) for ti in range(len(v[1]))]
    else:
        loc = {}
        for vi, v in enumerate(sol.veh):
            for ti, trip in enumerate(v[1]):
                for i in trip:
                    loc[i] = (vi, ti)
        cand_trips = []
        seen = set()
        for x in inst.nbr[u]:
            p = loc.get(x)
            if p is not None and p not in seen:
                seen.add(p)
                cand_trips.append(p)
    for vi, ti in cand_trips:
        k, trips = sol.veh[vi]
        trip = trips[ti]
        for pos in range(len(trip) + 1):
            nt = trip[:pos] + [u] + trip[pos:]
            ntr = trips[:ti] + [nt] + trips[ti + 1:]
            ev = eval_vehicle(inst, k, ntr)
            if ev is None:
                continue
            gv = g_after(sol, tot, w, sc, (vi,), [veh_obj(inst, k, ev)]) - g0
            if noise:
                gv += noise * abs(gv) * (2 * rng.random() - 1)
            if gv < bestv:
                bestv = gv
                best = ("trip", vi, ti, pos, ev)
    if inst.multitrip:
        vis = range(len(sol.veh)) if full else sorted({p[0] for p in cand_trips})
        for vi in vis:
            k, trips = sol.veh[vi]
            ntr = trips + [[u]]
            ev = eval_vehicle(inst, k, ntr)
            if ev is None:
                continue
            gv = g_after(sol, tot, w, sc, (vi,), [veh_obj(inst, k, ev)]) - g0
            if noise:
                gv += noise * abs(gv) * (2 * rng.random() - 1)
            if gv < bestv:
                bestv = gv
                best = ("newtrip", vi, None, None, ev)
    # new vehicle
    ks = range(inst.K)
    if type_probs is not None:
        ks = [pick_type(rng, type_probs)]
    for k in ks:
        ev = eval_vehicle(inst, k, [[u]])
        if ev is None:
            continue
        gv = g_after(sol, tot, w, sc, (), [veh_obj(inst, k, ev)]) - g0
        if noise:
            gv += noise * abs(gv) * (2 * rng.random() - 1)
        if gv < bestv:
            bestv = gv
            best = ("newveh", k, None, None, ev)
    if best is None and type_probs is not None:
        for k in range(inst.K):
            ev = eval_vehicle(inst, k, [[u]])
            if ev is not None:
                return ("newveh", k, None, None, ev)
    return best


def pick_type(rng, probs):
    r = rng.random()
    acc = 0.0
    for k, p in enumerate(probs):
        acc += p
        if r <= acc:
            return k
    return len(probs) - 1


def apply_insertion(sol, inst, u, mv):
    kind, x, ti, pos, ev = mv
    if kind == "trip":
        sol.veh[x][1][ti].insert(pos, u)
        sol.comp[x] = veh_obj(inst, sol.veh[x][0], ev)
    elif kind == "newtrip":
        sol.veh[x][1].append([u])
        sol.comp[x] = veh_obj(inst, sol.veh[x][0], ev)
    else:
        sol.veh.append([x, [[u]]])
        sol.comp.append(veh_obj(inst, x, ev))


def insert_all(sol, inst, custs, w, sc, rng, noise=0.0, type_probs=None):
    for u in custs:
        mv = best_insertion(sol, inst, u, w, sc, rng, noise, False, type_probs)
        if mv is None or mv[0] == "newveh":
            mv2 = best_insertion(sol, inst, u, w, sc, rng, noise, True, type_probs)
            if mv2 is not None:
                mv = mv2
        if mv is None:
            mv = best_insertion(sol, inst, u, w, sc, rng, 0.0, True, None)
        if mv is None:
            raise RuntimeError("customer %d cannot be served" % u)
        apply_insertion(sol, inst, u, mv)
    return sol


def construct(inst, rng, w, sc, type_probs=None, order="tw"):
    custs = list(range(1, inst.n + 1))
    if order == "tw":
        custs.sort(key=lambda i: (inst.b[i], rng.random()))
    else:
        rng.shuffle(custs)
    sol = Sol()
    for u in custs:
        mv = best_insertion(sol, inst, u, w, sc, rng, 0.0, True, type_probs)
        if mv is None:
            mv = best_insertion(sol, inst, u, w, sc, rng, 0.0, True, None)
        apply_insertion(sol, inst, u, mv)
    return sol


# ---------------------------------------------------------------- local search

def locate(sol):
    loc = {}
    for vi, v in enumerate(sol.veh):
        for ti, trip in enumerate(v[1]):
            for p, i in enumerate(trip):
                loc[i] = (vi, ti, p)
    return loc


class Counter:
    __slots__ = ("evals", "type_evals")

    def __init__(self):
        self.evals = 0
        self.type_evals = 0


def local_search(sol, inst, w, sc, rng, deadline, now, use_type=True, cnt=None, first_only=False):
    """First-improvement descent on g_w. Moves: relocate, intra-trip 2-opt,
    vehicle type change, trip transfer. If first_only, return True as soon as an
    improving move exists (used as probe; solution is not modified)."""
    improved_any = False
    while True:
        if now() > deadline:
            return improved_any
        improved = False
        tot = Totals(sol)
        g0 = gval(sol.objs(), w, sc)
        loc = locate(sol)
        custs = list(range(1, inst.n + 1))
        rng.shuffle(custs)
        for u in custs:
            moved = False
            vu, tu, pu = loc[u]
            ku, tripsu = sol.veh[vu]
            for v in inst.nbr[u]:
                vv, tv, pv = loc[v]
                kv, tripsv = sol.veh[vv]
                # relocate u next to v (before / after)
                for ins_after in (0, 1):
                    if vu == vv:
                        trips = [list(t) for t in tripsu]
                        trips[tu].pop(pu)
                        tpv = tv
                        tripv = trips[tpv]
                        try:
                            pvv = tripv.index(v)
                        except ValueError:
                            continue
                        tripv.insert(pvv + ins_after, u)
                        if trips == tripsu:
                            continue
                        trips = [t for t in trips if t]
                        ev = eval_vehicle(inst, ku, trips)
                        if cnt:
                            cnt.evals += 1
                        if ev is None:
                            continue
                        gn = g_after(sol, tot, w, sc, (vu,), [veh_obj(inst, ku, ev)])
                        if gn < g0 - 1e-9:
                            if first_only:
                                return True
                            sol.veh[vu][1] = trips
                            sol.comp[vu] = veh_obj(inst, ku, ev)
                            improved = moved = True
                            break
                    else:
                        tripsA = [list(t) for t in tripsu]
                        tripsA[tu].pop(pu)
                        tripsA = [t for t in tripsA if t]
                        tripsB = [list(t) for t in tripsv]
                        tripsB[tv].insert(pv + ins_after, u)
                        evB = eval_vehicle(inst, kv, tripsB)
                        if cnt:
                            cnt.evals += 1
                        if evB is None:
                            continue
                        if tripsA:
                            evA = eval_vehicle(inst, ku, tripsA)
                            if cnt:
                                cnt.evals += 1
                            if evA is None:
                                continue
                            nc = [veh_obj(inst, ku, evA), veh_obj(inst, kv, evB)]
                        else:
                            nc = [veh_obj(inst, kv, evB)]
                        gn = g_after(sol, tot, w, sc, (vu, vv), nc)
                        if gn < g0 - 1e-9:
                            if first_only:
                                return True
                            sol.veh[vv][1] = tripsB
                            sol.comp[vv] = veh_obj(inst, kv, evB)
                            if tripsA:
                                sol.veh[vu][1] = tripsA
                                sol.comp[vu] = veh_obj(inst, ku, evA)
                            else:
                                sol.veh[vu][1] = []
                                sol.clean()
                            improved = moved = True
                            break
                if moved:
                    break
                # intra-trip 2-opt (reverse segment between u and v)
                if vu == vv and tu == tv and pu != pv:
                    i, j = (pu, pv) if pu < pv else (pv, pu)
                    trip = tripsu[tu]
                    nt = trip[:i] + trip[i:j + 1][::-1] + trip[j + 1:]
                    trips = tripsu[:tu] + [nt] + tripsu[tu + 1:]
                    ev = eval_vehicle(inst, ku, trips)
                    if cnt:
                        cnt.evals += 1
                    if ev is not None:
                        gn = g_after(sol, tot, w, sc, (vu,), [veh_obj(inst, ku, ev)])
                        if gn < g0 - 1e-9:
                            if first_only:
                                return True
                            sol.veh[vu][1] = trips
                            sol.comp[vu] = veh_obj(inst, ku, ev)
                            improved = moved = True
                            break
            if moved:
                tot = Totals(sol)
                g0 = gval(sol.objs(), w, sc)
                loc = locate(sol)
            if now() > deadline:
                return improved_any or improved
        if not improved and use_type:
            # vehicle type change
            for vi in range(len(sol.veh)):
                k, trips = sol.veh[vi]
                for k2 in range(inst.K):
                    if k2 == k:
                        continue
                    ev = eval_vehicle(inst, k2, trips)
                    if cnt:
                        cnt.evals += 1
                        cnt.type_evals += 1
                    if ev is None:
                        continue
                    gn = g_after(sol, tot, w, sc, (vi,), [veh_obj(inst, k2, ev)])
                    if gn < g0 - 1e-9:
                        if first_only:
                            return True
                        sol.veh[vi][0] = k2
                        sol.comp[vi] = veh_obj(inst, k2, ev)
                        improved = True
                        break
                if improved:
                    break
        if not improved and inst.multitrip:
            # trip transfer to another vehicle or to a new vehicle
            nveh = len(sol.veh)
            for vi in range(nveh):
                k, trips = sol.veh[vi]
                for ti in range(len(trips)):
                    rest = trips[:ti] + trips[ti + 1:]
                    if rest:
                        evA = eval_vehicle(inst, k, rest)
                        if evA is None:
                            continue
                        cA = [veh_obj(inst, k, evA)]
                    else:
                        cA = []
                    targets = [(vj, sol.veh[vj][0]) for vj in range(nveh) if vj != vi]
                    if rest:
                        targets += [(-1, k2) for k2 in range(inst.K)]
                    for vj, k2 in targets:
                        if vj >= 0:
                            newtr = sol.veh[vj][1] + [trips[ti]]
                        else:
                            newtr = [trips[ti]]
                        ev = eval_vehicle(inst, k2, newtr)
                        if cnt:
                            cnt.evals += 1
                        if ev is None:
                            continue
                        idx = (vi, vj) if vj >= 0 else (vi,)
                        gn = g_after(sol, tot, w, sc, idx, cA + [veh_obj(inst, k2, ev)])
                        if gn < g0 - 1e-9:
                            if first_only:
                                return True
                            mtrip = trips[ti]
                            if vj >= 0:
                                sol.veh[vj][1] = newtr
                                sol.comp[vj] = veh_obj(inst, k2, ev)
                            else:
                                sol.veh.append([k2, [mtrip]])
                                sol.comp.append(veh_obj(inst, k2, ev))
                            sol.veh[vi][1] = rest
                            if rest:
                                sol.comp[vi] = cA[0]
                            sol.clean()
                            improved = True
                            break
                    if improved:
                        break
                if improved:
                    break
        if not improved:
            return improved_any
        improved_any = True


# ---------------------------------------------------------------- perturbation

def destroy_repair(sol, inst, w, sc, rng, strength, noise=0.2):
    """Remove a share of customers geographically close to a random seed and
    reinsert them greedily with noise."""
    n = inst.n
    k = max(2, int(round(strength * n)))
    seed = rng.randint(1, n)
    E = inst.E
    order = sorted(range(1, n + 1), key=lambda j: E[seed][j])
    rem = set(order[:k])
    s2 = Sol()
    for (kk, trips) in sol.veh:
        nt = [[i for i in t if i not in rem] for t in trips]
        nt = [t for t in nt if t]
        if nt:
            ev = eval_vehicle(inst, kk, nt)
            if ev is None:  # removing customers can only shorten a trip, but keep safe
                for t in nt:
                    for i in t:
                        rem.add(i)
                continue
            s2.veh.append([kk, nt])
            s2.comp.append(veh_obj(inst, kk, ev))
    lst = list(rem)
    rng.shuffle(lst)
    insert_all(s2, inst, lst, w, sc, rng, noise)
    return s2


# ---------------------------------------------------------------- archive

def crowding(points):
    n = len(points)
    if n <= 2:
        return [INF] * n
    d = [0.0] * n
    for m in range(3):
        idx = sorted(range(n), key=lambda i: points[i][m])
        lo, hi = points[idx[0]][m], points[idx[-1]][m]
        d[idx[0]] = d[idx[-1]] = INF
        if hi - lo <= 0:
            continue
        for r in range(1, n - 1):
            d[idx[r]] += (points[idx[r + 1]][m] - points[idx[r - 1]][m]) / (hi - lo)
    return d


class Archive:
    def __init__(self, cap=50):
        self.cap = cap
        self.items = []  # (objs, sol)

    def add(self, sol):
        o = sol.objs()
        o = (round(o[0], 6), round(o[1], 6), round(o[2], 6))
        for (p, _) in self.items:
            if dominates(p, o) or p == o:
                return False
        self.items = [(p, s) for (p, s) in self.items if not dominates(o, p)]
        self.items.append((o, sol.copy()))
        if len(self.items) > self.cap:
            cd = crowding([p for p, _ in self.items])
            worst = min(range(len(cd)), key=lambda i: cd[i])
            self.items.pop(worst)
        return True

    def points(self):
        return [p for p, _ in self.items]


def weight_vectors(H=4):
    ws = []
    for i in range(H + 1):
        for j in range(H + 1 - i):
            ws.append((i / H, j / H, (H - i - j) / H))
    return ws  # 15 vectors for H = 4


def nondominated(points):
    pts = list(set(points))
    return [p for p in pts if not any(dominates(q, p) for q in pts if q != p)]
