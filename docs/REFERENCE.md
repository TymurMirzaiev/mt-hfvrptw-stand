# Reference

Details moved out of the README: experiments, ablation configurations, result format, instances and the feasibility checker.

## Reproducing the experiments

```bash
python3 run.py <experiment> [workers]     # appends one JSON line per run to results/<experiment>.jsonl
python3 analyze.py                        # writes results/summary.json and prints the tables
```

| Experiment | What it does | Runs | Time on 2 cores |
|---|---|---|---|
| `solomon` | classical VRPTW special case on all 56 Solomon instances with 100 customers, 2 runs × 20 s | 112 | ~20 min |
| `abl_own` | ablation of the escape mechanisms (A1, A2, A3, A0, A0b) on 18 generated instances, 10 seeds; plus 36 probe runs | 936 | ~45 min |
| `abl_der` | the same on 18 Solomon-derived instances | 900 | ~45 min |
| `main_own` | MO-ILS vs weighted sum vs NSGA-II on 18 generated instances, 10 seeds | 540 | ~25 min |
| `main_der` | the same on 18 Solomon-derived instances | 540 | ~25 min |
| `fact` | 2×2×2 factorial study of model simplifications | 432 | ~20 min |
| `scale200` | three instances with 200 customers | 18 | ~3 min |

Time budgets per run: 2, 4, 8 and 20 s for 25, 50, 100 and 200 customers. The runner skips runs already present in the `.jsonl` file, so an interrupted experiment resumes where it stopped. To recompute from scratch, move the old file away first.

The budgets are wall-clock seconds, so absolute values depend on the machine. Rankings and significance should not.

### Ablation configurations (MPZIS-2026)

| Code | Change against the base method |
|---|---|
| A1 | base: destroy-and-repair with random strength 0.10–0.35, 15 scalarisation directions in cyclic order |
| A2 | one balanced direction (1/3, 1/3, 1/3) instead of 15 |
| A3 | no destroy-and-repair: a fresh construction every iteration (multistart) |
| A0 | naive adaptive strength: reset to 0.10 after an archive acceptance, ×1.3 after 8 idle iterations |
| A0b | cautious adaptive upper bound: ×1.3 after 8 idle iterations, ×0.9 after an acceptance, within 0.10–0.60 |

### Format of the raw results

One JSON object per line. Multi-objective runs contain `points` (objective vectors of the final archive: cost, emissions, makespan), `iters`, `nsol`, `bad` (checker violations, always 0), `shares` (fleet mix truck/van/bike), `spec` (instance), `cfg` (configuration) and `seed`. Solomon runs contain `cost`, `recost` (cost recomputed by the checker), `bks` and `ok`.

## The test instances

`dataset/` contains the 36 instances of the ablation study in a readable text format and in JSON, plus the 56 original Solomon instances. The format and all model rules (distances per vehicle type, service times, constraints, objectives) are described in [../dataset/README.md](../dataset/README.md). They are regenerated deterministically by `python3 export_instances.py dataset`.

## Feasibility checking

`mt/checker.py` is written separately from the search code. It recomputes loads, arrival times, time windows, reload times and shift length from the raw data for every solution and compares the objectives. In the published experiments it checked more than 39,000 solutions and found no violation.

