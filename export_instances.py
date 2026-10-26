"""Export the 36 MT-HFVRPTW test instances of the MPZIS-2026 ablation study
(18 generated + 18 Solomon-derived) to plain-text and JSON files.

    python3 export_instances.py <out_dir>

The files hold raw data only (coordinates, demands, windows, vehicle types, zone);
distance matrices follow from the rules written in each file header and README.
"""
import json
import os
import shutil
import sys

from mt.instances import generate, own_set, derive_mt_hfvrptw, DATA, ZONE_TRUCK, BIKE_FACTOR

SOL6 = ["R101", "R201", "C101", "C201", "RC101", "RC201"]


def fmt(x):
    return ("%.4f" % x).rstrip("0").rstrip(".") if isinstance(x, float) else str(x)


def write(inst, folder, source):
    os.makedirs(folder, exist_ok=True)
    m = inst.meta
    types = inst.types
    head = [
        ("NAME", inst.name), ("TYPE", "MT-HFVRPTW"), ("SOURCE", source),
        ("CUSTOMERS", inst.n), ("HORIZON", fmt(inst.H)), ("RELOAD_TIME", fmt(inst.rho)),
        ("FLEET", "unlimited number of vehicles of each type"),
        ("ZONE_CENTER", "%s %s" % (fmt(m["center"][0]), fmt(m["center"][1]))),
        ("ZONE_RADIUS", fmt(m["radius"])),
        ("TRUCK_ZONE_FACTOR", fmt(ZONE_TRUCK)), ("BIKE_FACTOR", fmt(BIKE_FACTOR)),
        ("DISTANCE", "Euclidean, not rounded; truck x TRUCK_ZONE_FACTOR on arcs whose midpoint lies in the zone; bike x BIKE_FACTOR on all arcs; van plain"),
        ("TRAVEL_TIME", "distance / speed of the type"),
        ("SERVICE_TIME", "node service x service_mult of the type"),
    ]
    lines = ["%s : %s" % kv for kv in head]
    lines.append("VEHICLE_TYPES_SECTION")
    lines.append("# id name mass_cap volume_cap fixed_cost cost_per_dist speed emission_per_dist service_mult")
    for k, t in enumerate(types, 1):
        lines.append(" ".join([str(k), t["name"]] + [fmt(float(t[x])) for x in ("Q", "W", "F", "c", "speed", "e", "sig")]))
    lines.append("NODE_SECTION")
    lines.append("# id x y mass volume ready due service   (node 0 = depot)")
    for i in range(inst.n + 1):
        lines.append(" ".join([str(i)] + [fmt(float(v)) for v in (inst.xy[i][0], inst.xy[i][1], inst.q[i], inst.vol[i], inst.a[i], inst.b[i], inst.s[i])]))
    lines.append("EOF")
    with open(os.path.join(folder, inst.name + ".txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    js = dict(name=inst.name, type="MT-HFVRPTW", source=source, customers=inst.n, horizon=inst.H,
              reload_time=inst.rho, fleet="unlimited",
              zone=dict(center=list(m["center"]), radius=m["radius"], truck_factor=ZONE_TRUCK, bike_factor=BIKE_FACTOR),
              vehicle_types=[dict(name=t["name"], mass_cap=t["Q"], volume_cap=t["W"], fixed_cost=t["F"],
                                  cost_per_dist=t["c"], speed=t["speed"], emission_per_dist=t["e"],
                                  service_mult=t["sig"]) for t in types],
              nodes=[dict(id=i, x=inst.xy[i][0], y=inst.xy[i][1], mass=inst.q[i], volume=inst.vol[i],
                          ready=inst.a[i], due=inst.b[i], service=inst.s[i]) for i in range(inst.n + 1)])
    with open(os.path.join(folder, inst.name + ".json"), "w") as f:
        json.dump(js, f, indent=1)


def main(out):
    for (c, n, i) in own_set():
        seed = {"R": 1, "C": 2, "RC": 3}[c] * 100000 + n * 10 + i
        write(generate(c, n, i), os.path.join(out, "generated"),
              "own generator, class %s, %d customers, index %d, seed %d" % (c, n, i, seed))
    for n in (25, 50, 100):
        for nm in SOL6:
            write(derive_mt_hfvrptw(nm, n), os.path.join(out, "solomon_derived"),
                  "derived from Solomon %s, first %d customers" % (nm, n))
    dst = os.path.join(out, "solomon_original")
    os.makedirs(dst, exist_ok=True)
    for f in sorted(os.listdir(DATA)):
        if f.endswith((".vrp", ".sol")):
            shutil.copy(os.path.join(DATA, f), dst)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "dataset")
