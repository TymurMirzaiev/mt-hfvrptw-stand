import json
from collections import defaultdict
from mt.metrics import bounds, normalise, hv3, wil
from mt.core import nondominated
R=[json.loads(l) for l in open('results/main_own.jsonl')]
byi=defaultdict(list)
for r in R: byi[tuple(r['spec'])].append(r)
for ref in (1.0,1.1,1.5):
    m={c:{} for c in ('MO','WS','NS')}
    for ik,rs in byi.items():
        ideal,nadir=bounds([r['points'] for r in rs])
        for c in m:
            v=[hv3(normalise(r['points'],ideal,nadir),(ref,)*3) for r in rs if r['cfg']==c]
            m[c][ik]=sum(v)/len(v)
    ks=sorted(m['MO'])
    print(ref, {c:round(sum(m[c].values())/len(ks),3) for c in m}, 'p MO-WS %.2g'%wil([m['MO'][k] for k in ks],[m['WS'][k] for k in ks]), 'p MO-NS %.2g'%wil([m['MO'][k] for k in ks],[m['NS'][k] for k in ks]))
