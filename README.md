# 🚚 🚐 🚲 MT-HFVRPTW stand

**Plan city deliveries with trucks, vans and cargo bikes, and get a set of good plans instead of just one.**

You give it a depot, orders (weight, volume, time window) and three vehicle types. It returns plans that trade off **cost**, **emissions** and **finish time**, so you can pick the one that suits your day.

Pure Python, no external solver. Українською — [нижче](#-українською).

## ⚡ Try it in a minute

```bash
pip install -r requirements.txt
python3 example.py dataset/generated/R25_0.json --time 2
```

```
      cost  emissions   makespan    truck    van   bike   check
     813.0      174.2      538.9        0      1      3   ok
     862.5      174.2      531.7        0      1      4   ok
     949.6      209.3      523.7        0      2      3   ok
     ...
```

Each row is one plan. `check: ok` means an independent checker confirmed that every rule holds. Add `--out plans.json` to save the routes.

## 📦 What's inside

| Folder | What you'll find |
|---|---|
| `mt/` | the solver: model, search method, feasibility checker, quality metrics |
| `dataset/` | 36 ready-to-use test instances + 56 classic Solomon instances ([format](dataset/README.md)) |
| `results/` | raw results of all 3,478 runs from the papers |
| `figs/` | figures (PNG, EPS) |
| `docs/REFERENCE.md` | how to rerun every experiment, configurations, result format |

## 🔁 Rerun the paper experiments

```bash
python3 run.py abl_own 2   # ablation study from the MPZIS-2026 paper (~45 min on 2 cores)
python3 analyze.py         # summary tables → results/summary.json
```

All experiments are listed in [docs/REFERENCE.md](docs/REFERENCE.md).

## 📄 Papers

- T. R. Mirzaiev, H. Y. Mykhalchuk. *Escaping local optima in multi-objective iterated local search for heterogeneous fleet vehicle routing: an ablation study.* MPZIS-2026, Dnipro.
- T. R. Mirzaiev, H. Y. Mykhalchuk. *Multi-objective multi-trip heterogeneous fleet vehicle routing with time windows under vehicle-type-dependent parameters.* Manuscript, 2026.

Using the code or data? Please cite it (see `CITATION.cff`).

## 👋 Contact

Tymur Mirzaiev · Oles Honchar Dnipro National University · timur.mirzaiev@gmail.com · [ORCID](https://orcid.org/0009-0006-0158-4329)

---

## 🇺🇦 Українською

**Що це.** Програма, яка планує доставку вантажівками, фургонами й вантажними велосипедами. Вона видає не один план, а кілька компромісних між вартістю, викидами й часом завершення, а ви обираєте потрібний.

**Як спробувати:**

```bash
pip install -r requirements.txt
python3 example.py dataset/generated/R25_0.json --time 2
```

**Що де лежить:** `mt/` — сам метод, `dataset/` — тестові задачі, `results/` — результати з наших статей, `docs/REFERENCE.md` — як повторити всі експерименти.
