"""Aggregate results/*.jsonl into summary tables (results/summary.json + stdout)."""
import json
import statistics as st
import sys
from collections import defaultdict

from mt.metrics import bounds, normalise, hv3, coverage, supported, wil, holm
from mt.core import nondominated


def load(exp):
    try:
        with open("results/%s.jsonl" % exp) as f:
            return [json.loads(l) for l in f]
    except FileNotFoundError:
        return []


def inst_key(r):
    return tuple(r["spec"])


def size(r):
    return r["spec"][2]


def mean(x):
    return sum(x) / len(x) if x else float("nan")


def solomon():
    R = load("solomon")
    if not R:
        return None
    by = defaultdict(list)
    for r in R:
        by[r["name"]].append(r)
    gaps, opt, ok = {}, 0, all(r["ok"] for r in R)
    cons = all(abs(r["cost"] - r["recost"]) < 1e-6 for r in R)
    for nm, rs in by.items():
        best = min(x["cost"] for x in rs)
        bks = rs[0]["bks"]
        g = 100 * (best / bks - 1)
        gaps[nm] = g
        if g <= 0.005:
            opt += 1
    allg = list(gaps.values())
    cls = defaultdict(list)
    for nm, g in gaps.items():
        c = nm.rstrip("0123456789")
        cls[c].append(g)
    meanrun = mean([100 * (r["cost"] / r["bks"] - 1) for r in R])
    return dict(n_inst=len(gaps), runs=len(R), mean_gap_best=mean(allg), median_gap=st.median(allg),
                max_gap=max(allg), reached=opt, mean_gap_run=meanrun, checker_ok=ok, cost_consistent=cons,
                by_class={c: mean(v) for c, v in sorted(cls.items())},
                mean_iters=mean([r["iters"] for r in R]), gaps=gaps)


def per_instance_hv(R, groupkey):
    """Normalise per instance over the union of all approximations in R."""
    byi = defaultdict(list)
    for r in R:
        byi[inst_key(r)].append(r)
    for ik, rs in byi.items():
        ideal, nadir = bounds([r["points"] for r in rs if r["points"]])
        for r in rs:
            r["hv"] = hv3(normalise(r["points"], ideal, nadir))
    return byi


def means_by(R, cfgs, field, filt=lambda r: True):
    """per-instance mean of field for each cfg -> {cfg: {inst: mean}}"""
    d = {c: defaultdict(list) for c in cfgs}
    for r in R:
        if r["cfg"] in d and filt(r):
            d[r["cfg"]][inst_key(r)].append(field(r))
    return {c: {i: mean(v) for i, v in d[c].items()} for c in cfgs}


def paired(m, a, b):
    ks = sorted(set(m[a]) & set(m[b]))
    return [m[a][k] for k in ks], [m[b][k] for k in ks]


def ablation(exp):
    R = [r for r in load(exp) if r["cfg"] != "A1probe"]
    if not R:
        return None
    cfgs = ["A1", "A2", "A3", "A0", "A0b"]
    per_instance_hv(R, "cfg")
    hv = means_by(R, cfgs, lambda r: r["hv"])
    pts = means_by(R, cfgs, lambda r: len(r["points"]))
    its = means_by(R, cfgs, lambda r: r["iters"])
    out = {}
    for c in cfgs:
        out[c] = dict(hv=mean(list(hv[c].values())), points=mean(list(pts[c].values())),
                      iters=mean(list(its[c].values())), runs=sum(1 for r in R if r["cfg"] == c))
    tests = {}
    for a, b in [("A1", "A2"), ("A1", "A3"), ("A1", "A0"), ("A1", "A0b"), ("A0b", "A0")]:
        x, y = paired(hv, a, b)
        tests["hv %s-%s" % (a, b)] = wil(x, y)
    for a, b in [("A1", "A2"), ("A1", "A3"), ("A1", "A0"), ("A1", "A0b"), ("A0b", "A0")]:
        x, y = paired(pts, a, b)
        tests["pts %s-%s" % (a, b)] = wil(x, y)
    names = list(tests)
    adj = holm([tests[n] for n in names])
    tests_holm = dict(zip(names, adj))
    # by size
    bysize = {}
    for n in (25, 50, 100):
        f = lambda r, n=n: size(r) == n
        hvn = means_by(R, cfgs, lambda r: r["hv"], f)
        ptn = means_by(R, cfgs, lambda r: len(r["points"]), f)
        bysize[n] = {c: dict(hv=mean(list(hvn[c].values())), points=mean(list(ptn[c].values()))) for c in cfgs}
    # contribution of each configuration to the joint approximation, per size
    contrib = {}
    byi = defaultdict(list)
    for r in R:
        byi[inst_key(r)].append(r)
    for n in (25, 50, 100):
        cc = defaultdict(int)
        for ik, rs in byi.items():
            if ik[2] != n:
                continue
            allp = []
            for r in rs:
                allp += [(tuple(p), r["cfg"]) for p in r["points"]]
            nd = set(nondominated([p for p, _ in allp]))
            for p, c in set(allp):
                if p in nd:
                    cc[c] += 1
        contrib[n] = dict(cc)
    bad = sum(r["bad"] for r in R)
    nsol = sum(r["nsol"] for r in R)
    res = dict(cfg=out, tests=tests, tests_holm=tests_holm, by_size=bysize, contrib=contrib,
               checked=nsol, bad=bad, runs=len(R))
    P = [r for r in load(exp) if r["cfg"] == "A1probe"]
    if P:
        sw = sum(r["switches"] for r in P)
        im = sum(r["improving"] for r in P)
        per = [r["improving"] / r["switches"] for r in P if r["switches"]]
        bys = {n: mean([r["improving"] / r["switches"] for r in P if r["switches"] and size(r) == n]) for n in (25, 50, 100)}
        res["probe"] = dict(runs=len(P), switches=sw, improving=im, pooled=im / sw if sw else None,
                            per_run_mean=mean(per), min=min(per), max=max(per), by_size=bys,
                            checked=sum(r["nsol"] for r in P), bad=sum(r["bad"] for r in P))
    return res


def main_exp(exp):
    R = load(exp)
    if not R:
        return None
    cfgs = ["MO", "WS", "NS"]
    byi = per_instance_hv(R, "cfg")
    hv = means_by(R, cfgs, lambda r: r["hv"])
    pts = means_by(R, cfgs, lambda r: len(r["points"]))
    cost = means_by(R, cfgs, lambda r: min(p[0] for p in r["points"]))
    its = means_by(R, cfgs, lambda r: r["iters"])
    nveh = means_by(R, cfgs, lambda r: r["nveh"])
    trips = means_by(R, cfgs, lambda r: r["trips"])
    out = {c: dict(hv=mean(list(hv[c].values())), points=mean(list(pts[c].values())),
                   best_cost=mean(list(cost[c].values())), iters=mean(list(its[c].values())),
                   nveh=mean(list(nveh[c].values())), trips=mean(list(trips[c].values())),
                   trips_per_veh=mean(list(trips[c].values())) / mean(list(nveh[c].values())))
           for c in cfgs}
    tests = {}
    for a, b in [("MO", "WS"), ("MO", "NS"), ("WS", "NS")]:
        for nm, m in (("hv", hv), ("pts", pts), ("cost", cost)):
            x, y = paired(m, a, b)
            tests["%s %s-%s" % (nm, a, b)] = wil(x, y)
    bysize = {}
    for n in (25, 50, 100):
        f = lambda r, n=n: size(r) == n
        d = {}
        for nm, fld in (("hv", lambda r: r["hv"]), ("points", lambda r: len(r["points"])), ("iters", lambda r: r["iters"])):
            m = means_by(R, cfgs, fld, f)
            d[nm] = {c: mean(list(m[c].values())) for c in cfgs}
            if nm in ("hv",):
                for a, b in [("MO", "WS"), ("MO", "NS")]:
                    x, y = paired(m, a, b)
                    d["p_%s_%s_%s" % (nm, a, b)] = wil(x, y)
        bysize[n] = d
    # coverage MO vs NS (pooled over runs, per instance union of each method)
    cov = {}
    for a, b in [("MO", "NS"), ("NS", "MO"), ("MO", "WS"), ("WS", "MO")]:
        vals = []
        for ik, rs in byi.items():
            A = nondominated([tuple(p) for r in rs if r["cfg"] == a for p in r["points"]])
            B = nondominated([tuple(p) for r in rs if r["cfg"] == b for p in r["points"]])
            vals.append(coverage(A, B))
        cov["%s>%s" % (a, b)] = mean(vals)
    # unsupported share on the best-known approximation per instance
    uns = {}
    tot = defaultdict(lambda: [0, 0])
    per_inst = []
    sens = {}
    for eps in (1e-9, 1e-6, 1e-4):
        u = n_all = 0
        for ik, rs in byi.items():
            ideal, nadir = bounds([r["points"] for r in rs])
            nd = nondominated([tuple(p) for r in rs for p in r["points"]])
            npts = normalise(nd, ideal, nadir)
            s = supported(npts, eps)
            u += s.count(False)
            n_all += len(s)
            if eps == 1e-6:
                tot[ik[2]][0] += s.count(False)
                tot[ik[2]][1] += len(s)
                per_inst.append(s.count(False) / len(s))
        sens[eps] = u / n_all
    uns = dict(by_size={n: dict(unsupported=v[0], points=v[1], share=v[0] / v[1]) for n, v in sorted(tot.items())},
               total_share=sens[1e-6], sensitivity=sens, per_inst_min=min(per_inst), per_inst_max=max(per_inst))
    bad = sum(r["bad"] for r in R)
    nsol = sum(r["nsol"] for r in R)
    return dict(cfg=out, tests=tests, by_size=bysize, coverage=cov, unsupported=uns, checked=nsol, bad=bad, runs=len(R))


def fact():
    R = load("fact")
    if not R:
        return None
    by = defaultdict(list)
    for r in R:
        by[(r["P"], r["O"], r["M"])].append(r)
    out = {}
    full_inst = defaultdict(list)
    for r in by[(True, True, True)]:
        full_inst[inst_key(r)].append(r["feas"] / r["nsol"])
    full_m = {k: mean(v) for k, v in full_inst.items()}
    for cfg, rs in sorted(by.items(), reverse=True):
        feas = sum(r["feas"] for r in rs) / sum(r["nsol"] for r in rs)
        viol = sum(r["nviol"] for r in rs) / sum(r["nsol"] for r in rs)
        sh = [mean([r["shares"][k] for r in rs]) for k in range(3)]
        cd = [x for r in rs for x in r["cdiff"]]
        pi = defaultdict(list)
        for r in rs:
            pi[inst_key(r)].append(r["feas"] / r["nsol"])
        m = {k: mean(v) for k, v in pi.items()}
        ks = sorted(m)
        p = wil([full_m[k] for k in ks], [m[k] for k in ks]) if cfg != (True, True, True) else None
        out["%s%s%s" % ("P" if cfg[0] else "-", "O" if cfg[1] else "-", "M" if cfg[2] else "-")] = dict(
            feasible=feas, viol_per_plan=viol, shares=sh, trips=mean([r["trips"] for r in rs]),
            max_cost_dev=max(cd) if cd else None, p_vs_full=p, runs=len(rs))
    return dict(cfg=out, runs=len(R))


def scale200():
    R = load("scale200")
    if not R:
        return None
    per_instance_hv(R, "cfg")
    out = {}
    for c in ("R", "C", "RC"):
        d = {}
        for m in ("MO", "NS"):
            rs = [r for r in R if r["spec"][1] == c and r["cfg"] == m]
            d[m] = dict(hv=mean([r["hv"] for r in rs]), points=mean([len(r["points"]) for r in rs]),
                        iters=mean([r["iters"] for r in rs]))
        out[c] = d
    return out


if __name__ == "__main__":
    S = dict(solomon=solomon(), abl_own=ablation("abl_own"), abl_der=ablation("abl_der"),
             main_own=main_exp("main_own"), main_der=main_exp("main_der"), fact=fact(), scale200=scale200())
    with open("results/summary.json", "w") as f:
        json.dump(S, f, indent=1, default=str)
    which = sys.argv[1:] or list(S)
    for k in which:
        print("=====", k)
        print(json.dumps(S[k], indent=1, default=str)[:6000])
