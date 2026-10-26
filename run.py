"""Experiment runner. Usage: python3 run.py <experiment> [workers]
Results are appended to results/<experiment>.jsonl (one JSON record per run), so an
interrupted experiment resumes where it stopped."""
import json
import os
import sys
import time
from functools import lru_cache
from multiprocessing import Pool

from mt.instances import (read_solomon_vrptw, solomon_bks, solomon_names, derive_mt_hfvrptw,
                          generate, own_set)
from mt.algorithms import ils_single, mo_ils, ws_seq, nsga2, scales
from mt.checker import check

BUD = {25: 2.0, 50: 4.0, 100: 8.0, 200: 20.0}
SOL6 = ["R101", "R201", "C101", "C201", "RC101", "RC201"]


@lru_cache(maxsize=64)
def get_inst(spec, P=True, O=True, M=True):
    if spec[0] == "own":
        return generate(spec[1], spec[2], spec[3], P, O, M)
    if spec[0] == "der":
        return derive_mt_hfvrptw(spec[1], spec[2], P, O, M)
    raise ValueError(spec)


@lru_cache(maxsize=64)
def get_scales(spec):
    return scales(get_inst(spec))


def fleet(inst, sols):
    cnt = [0] * inst.K
    trips = 0
    nveh = 0
    for s in sols:
        for k, tr in s.veh:
            cnt[k] += 1
            trips += len(tr)
            nveh += 1
    tot = max(1, sum(cnt))
    return [c / tot for c in cnt], trips / max(1, len(sols)), nveh / max(1, len(sols))


def task(t):
    kind = t["kind"]
    t0 = time.time()
    if kind == "sol":
        inst = read_solomon_vrptw(t["name"])
        r = ils_single(inst, t["budget"], t["seed"])
        ok, viol, rec = check(inst, r["best"])
        out = dict(cost=r["cost"], recost=rec[0], bks=solomon_bks(t["name"]), iters=r["iters"],
                   ok=ok, nveh=len(r["best"].veh))
    elif kind == "mo":
        spec = tuple(t["spec"])
        inst = get_inst(spec)
        sc = get_scales(spec)
        n = inst.n
        b = BUD[n]
        m = t["method"]
        if m == "MO":
            r = mo_ils(inst, b, t["seed"], sc=sc, **t.get("var", {}))
        elif m == "WS":
            r = ws_seq(inst, b, t["seed"], sc=sc)
        else:
            r = nsga2(inst, b, t["seed"], sc=sc)
        bad = 0
        for s in r["sols"]:
            ok, viol, rec = check(inst, s, s.objs())
            bad += (not ok)
        shares, trips, nveh = fleet(inst, r["sols"])
        out = dict(points=r["points"], iters=r["iters"], evals=r["evals"], nsol=len(r["sols"]),
                   bad=bad, shares=shares, trips=trips, nveh=nveh,
                   switches=r.get("switches", 0), improving=r.get("improving", 0))
    elif kind == "fact":
        spec = tuple(t["spec"])
        P, O, M = t["P"], t["O"], t["M"]
        full = get_inst(spec)
        plan = get_inst(spec, P, O, M)
        sc = get_scales(spec)
        r = mo_ils(plan, BUD[full.n], t["seed"], sc=sc)
        feas = 0
        nviol = 0
        cdiff = []
        for s in r["sols"]:
            ok, viol, rec = check(full, s)
            feas += ok
            nviol += len(viol)
            if ok:
                cdiff.append(abs(s.objs()[0] - rec[0]) / rec[0])
        shares, trips, nveh = fleet(full, r["sols"])
        out = dict(nsol=len(r["sols"]), feas=feas, nviol=nviol, cdiff=cdiff, shares=shares,
                   trips=trips, nveh=nveh, iters=r["iters"])
    else:
        raise ValueError(kind)
    out.update(t)
    out["wall"] = time.time() - t0
    return out


def tasks_for(exp):
    T = []
    if exp == "solomon":
        for name in solomon_names():
            for seed in (1, 2):
                T.append(dict(kind="sol", name=name, seed=seed, budget=20.0))
    elif exp in ("abl_own", "abl_der"):
        specs = ([("own", c, n, i) for (c, n, i) in own_set()] if exp == "abl_own"
                 else [("der", nm, n) for n in (25, 50, 100) for nm in SOL6])
        V = {"A1": {}, "A2": {"dirs": "single"}, "A3": {"perturb": "restart"},
             "A0": {"strength": "naive"}, "A0b": {"strength": "cautious"}}
        for seed in range(1, 11):
            for sp in specs:
                for cfg, var in V.items():
                    T.append(dict(kind="mo", spec=sp, method="MO", var=var, cfg=cfg, seed=seed))
        if exp == "abl_own":
            for seed in (101, 102):
                for sp in specs:
                    T.append(dict(kind="mo", spec=sp, method="MO", var={"probe": True}, cfg="A1probe", seed=seed))
    elif exp in ("main_own", "main_der"):
        specs = ([("own", c, n, i) for (c, n, i) in own_set()] if exp == "main_own"
                 else [("der", nm, n) for n in (25, 50, 100) for nm in SOL6])
        for seed in range(1, 11):
            for sp in specs:
                for m in ("MO", "WS", "NS"):
                    T.append(dict(kind="mo", spec=sp, method=m, cfg=m, seed=seed))
    elif exp == "scale200":
        for c in ("R", "C", "RC"):
            for seed in (1, 2, 3):
                for m in ("MO", "NS"):
                    T.append(dict(kind="mo", spec=("own", c, 200, 0), method=m, cfg=m, seed=seed))
    elif exp == "fact":
        for seed in (1, 2, 3):
            for (c, n, i) in own_set():
                for P in (True, False):
                    for O in (True, False):
                        for M in (True, False):
                            T.append(dict(kind="fact", spec=("own", c, n, i), P=P, O=O, M=M, seed=seed))
    return T


def key(t):
    return json.dumps({k: t[k] for k in sorted(t) if k in
                       ("kind", "name", "spec", "method", "cfg", "seed", "P", "O", "M")}, sort_keys=True)


if __name__ == "__main__":
    exp = sys.argv[1]
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    os.makedirs("results", exist_ok=True)
    path = "results/%s.jsonl" % exp
    done = set()
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                done.add(key(json.loads(line)))
    T = [t for t in tasks_for(exp) if key(t) not in done]
    print(exp, "tasks to run:", len(T), flush=True)
    with Pool(workers) as pool, open(path, "a") as f:
        for i, r in enumerate(pool.imap_unordered(task, T)):
            f.write(json.dumps(r) + "\n")
            f.flush()
            if (i + 1) % 20 == 0:
                print(exp, i + 1, "/", len(T), flush=True)
    print(exp, "done", flush=True)
