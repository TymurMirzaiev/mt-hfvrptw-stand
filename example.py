"""Solve one MT-HFVRPTW instance with the multi-objective ILS and print the trade-offs.

    python3 example.py dataset/generated/R25_0.json --time 2 --seed 1 [--out plans.json]
"""
import argparse
import json

from mt.instances import load_json
from mt.algorithms import mo_ils
from mt.checker import check


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("instance", help="path to an instance .json from dataset/")
    ap.add_argument("--time", type=float, default=2.0, help="time budget in seconds (default 2)")
    ap.add_argument("--seed", type=int, default=1, help="random seed (default 1)")
    ap.add_argument("--out", help="write the plans to this JSON file")
    args = ap.parse_args()

    inst = load_json(args.instance)
    res = mo_ils(inst, args.time, args.seed)
    plans = []
    for s in sorted(res["sols"], key=lambda s: s.objs()[0]):
        ok, viol, _ = check(inst, s, s.objs())
        cost, em, mk = s.objs()
        use = [0] * inst.K
        for k, _tr in s.veh:
            use[k] += 1
        plans.append(dict(cost=round(cost, 2), emissions=round(em, 2), makespan=round(mk, 2),
                          vehicles={inst.types[k]["name"]: use[k] for k in range(inst.K)},
                          feasible=ok,
                          routes=[dict(type=inst.types[k]["name"], trips=[list(t) for t in tr]) for k, tr in s.veh]))
    print("%s: %d customers, %d iterations, %d trade-offs\n" % (inst.name, inst.n, res["iters"], len(plans)))
    names = [t["name"] for t in inst.types]
    print("%10s %10s %10s   %s   %s" % ("cost", "emissions", "makespan", " ".join("%6s" % n for n in names), "check"))
    for p in plans:
        print("%10.1f %10.1f %10.1f   %s   %s" % (p["cost"], p["emissions"], p["makespan"],
              " ".join("%6d" % p["vehicles"][n] for n in names), "ok" if p["feasible"] else "VIOLATION"))
    if args.out:
        with open(args.out, "w") as f:
            json.dump(dict(instance=inst.name, time=args.time, seed=args.seed, plans=plans), f, indent=1)
        print("\nplans written to", args.out)


if __name__ == "__main__":
    main()
