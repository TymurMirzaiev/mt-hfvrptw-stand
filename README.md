# MT-HFVRPTW stand

Multi-objective routing of a **heterogeneous fleet** (truck, van, cargo bike) with **multiple trips per vehicle** and **time windows**, where each vehicle type has its own road network, service speed and mass/volume capacity. Objectives: total cost, emissions and makespan.

This repository holds the code, test instances and raw results behind:

- T. R. Mirzaiev, H. Y. Mykhalchuk. *Escaping local optima in multi-objective iterated local search for heterogeneous fleet vehicle routing: an ablation study.* MPZIS-2026, Dnipro, 2026.
- T. R. Mirzaiev, H. Y. Mykhalchuk. *Multi-objective multi-trip heterogeneous fleet vehicle routing with time windows under vehicle-type-dependent parameters.* Manuscript, 2026.

Pure Python 3. No external solver; NumPy and SciPy are used only for statistics and one linear-programming test.

Українська інструкція — [нижче](#інструкція-українською).

---

## Contents

```
mt/                    the solver
  core.py              model, route evaluation, insertion, local search, destroy-and-repair, archive
  algorithms.py        MO-ILS (with ablation variants), weighted sum, memetic NSGA-II, single-objective ILS
  instances.py         instance generator, Solomon reader, Solomon-derived instances, JSON loader
  checker.py           independent feasibility checker (recomputes everything from raw data)
  metrics.py           hypervolume (exact 3-D), coverage, supportedness LP, Wilcoxon, Holm
example.py             solve one instance and print the trade-offs
run.py                 experiment runner
analyze.py             aggregate raw results into results/summary.json
figs.py                figures of the manuscript (needs matplotlib)
extra.py               hypervolume reference-point robustness check
export_instances.py    write the test instances to dataset/
dataset/               36 test instances (.txt and .json) + 56 original Solomon instances; see dataset/README.md
data/                  Solomon instances and best-known solutions (PyVRP instance collection, MIT licence)
results/               raw results of every run (*.jsonl) and the aggregated summary.json
figs/                  figures (PNG and EPS)
```

## Requirements

- Python 3.10 or newer (tested on 3.13)
- `pip install -r requirements.txt` (NumPy, SciPy; matplotlib only for `figs.py`)

## Quick start

Solve one instance for 2 seconds and print the set of trade-off plans:

```bash
python3 example.py dataset/generated/R25_0.json --time 2 --seed 1
```

```
R25_0: 25 customers, 459 iterations, 12 trade-offs

      cost  emissions   makespan    truck    van   bike   check
     813.0      174.2      538.9        0      1      3   ok
     820.4      166.7      538.9        0      1      3   ok
     ...
```

Each row is one plan; none is better than another in all three objectives at once. `check` is the verdict of the independent checker. Add `--out plans.json` to save the routes.

Use any file from `dataset/generated/` or `dataset/solomon_derived/`. Your own instance can be written in the same JSON format (see `dataset/README.md`).

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

`dataset/` contains the 36 instances of the ablation study in a readable text format and in JSON, plus the 56 original Solomon instances. The format and all model rules (distances per vehicle type, service times, constraints, objectives) are described in [dataset/README.md](dataset/README.md). They are regenerated deterministically by `python3 export_instances.py dataset`.

## Feasibility checking

`mt/checker.py` is written separately from the search code. It recomputes loads, arrival times, time windows, reload times and shift length from the raw data for every solution and compares the objectives. In the published experiments it checked more than 39,000 solutions and found no violation.

## Citation

If you use the code or the instances, please cite the papers above. Citation metadata is in `CITATION.cff`.

## Contact

Tymur Mirzaiev, Oles Honchar Dnipro National University — timur.mirzaiev@gmail.com · ORCID [0009-0006-0158-4329](https://orcid.org/0009-0006-0158-4329)

---

## Інструкція українською

**Що це.** Програма, на якій пораховано результати тез МПЗІС-2026 і статті про багатокритеріальну маршрутизацію різнорідного парку з повторними рейсами. Метод шукає набір компромісних планів доставки між вартістю, викидами й часом завершення роботи.

**Встановлення.** Потрібен Python 3.10 або новіший.

```bash
pip install -r requirements.txt
```

**Розв’язати одну задачу** (2 секунди, вивести всі компромісні плани):

```bash
python3 example.py dataset/generated/R25_0.json --time 2
```

Кожен рядок — окремий план: вартість, викиди, час завершення, скільки машин кожного типу і результат незалежної перевірки (`ok`). Щоб зберегти маршрути у файл, додайте `--out plans.json`.

**Повторити експерименти тез:**

```bash
python3 run.py abl_own 2    # абляція на 18 згенерованих задачах, 2 процеси (~45 хв)
python3 run.py abl_der 2    # абляція на 18 задачах, похідних від Соломона (~45 хв)
python3 run.py solomon 2    # перевірка на 56 задачах Соломона (~20 хв)
python3 analyze.py          # зведені таблиці в results/summary.json
```

Готові сирі результати всіх прогонів уже лежать у `results/`. Якщо файл результатів існує, `run.py` пропускає вже пораховані прогони; щоб порахувати заново, спершу перемістіть старий файл.

**Задачі.** Папка `dataset/` містить 18 згенерованих задач, 18 задач, похідних від набору Соломона, і 56 оригінальних задач Соломона. Формат і правила моделі описано в `dataset/README.md`.
