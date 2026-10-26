"""Independent feasibility checker. Recomputes everything from the raw instance
fields (coordinates, profile matrices, type parameters) without using the cached
travel-time and service-time tables of the search code."""


def check(inst, sol, objs=None, tol=1e-6):
    """Return (ok, violations list, recomputed objectives)."""
    viol = []
    seen = {}
    cost = em = 0.0
    mk = 0.0
    for vi, (k, trips) in enumerate(sol.veh):
        tp = inst.types[k]
        Dk = inst.D[k]
        t = inst.a[0]
        L = 0.0
        for ti, trip in enumerate(trips):
            if ti > 0:
                if not inst.multitrip:
                    viol.append("vehicle %d: multiple trips not allowed" % vi)
                t += inst.rho
            mass = sum(inst.q[i] for i in trip)
            volume = sum(inst.vol[i] for i in trip)
            if mass > tp["Q"] + tol:
                viol.append("vehicle %d trip %d: mass %.2f > %.2f" % (vi, ti, mass, tp["Q"]))
            if volume > tp["W"] + tol:
                viol.append("vehicle %d trip %d: volume %.2f > %.2f" % (vi, ti, volume, tp["W"]))
            prev = 0
            for i in trip:
                seen[i] = seen.get(i, 0) + 1
                d = Dk[prev][i]
                L += d
                t += d / tp["speed"]
                if t > inst.b[i] + tol:
                    viol.append("customer %d: arrival %.2f after window end %.2f" % (i, t, inst.b[i]))
                t = max(t, inst.a[i])
                t += inst.s[i] * tp["sig"]
                prev = i
            L += Dk[prev][0]
            t += Dk[prev][0] / tp["speed"]
            if t > inst.H + tol:
                viol.append("vehicle %d: return %.2f after horizon %.2f" % (vi, t, inst.H))
        cost += tp["F"] + tp["c"] * L
        em += tp["e"] * L
        mk = max(mk, t)
    for i in range(1, inst.n + 1):
        c = seen.get(i, 0)
        if c != 1:
            viol.append("customer %d served %d times" % (i, c))
    rec = (cost, em, mk)
    if objs is not None:
        for j in range(3):
            if abs(rec[j] - objs[j]) > 1e-4 * max(1.0, abs(rec[j])):
                viol.append("objective %d mismatch: %.6f vs %.6f" % (j, objs[j], rec[j]))
    return (len(viol) == 0), viol, rec
