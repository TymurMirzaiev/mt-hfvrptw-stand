# 🚚 🚐 🚲 MT-HFVRPTW stand

**Plan city deliveries with trucks, vans and cargo bikes, and get a set of good plans instead of just one.**

You give it a depot, orders (weight, volume, time window) and three vehicle types. It returns plans that trade off **cost**, **emissions** and **finish time**, so you can pick the one that suits your day.

Pure Python, no external solver. Українською — [нижче](#-українською).

## 🧩 The problem

Mixed fleets are great for cities, but hard to plan:

- **Vehicles differ in more than size.** A truck carries 220 units but detours around the city centre; a cargo bike carries 32, is slower and emits nothing, but can take car-free shortcuts. Each type also unloads at its own pace.
- **Weight and volume are separate limits.** A van can be full by volume long before it is full by weight.
- **Each vehicle can make several trips** a shift, reloading at the depot in between.
- **There is no single best plan.** Cheaper usually means more emissions or a later finish. A dispatcher needs the trade-offs, not one number.

Even for 100 orders, trying every plan is impossible, so we search smartly.

## 🧠 How the algorithm works

A **multi-objective iterated local search** that keeps an archive of the best trade-offs found so far (up to 50). Each round it:

1. **picks a viewpoint**, i.e. a set of weights for cost / emissions / time (15 viewpoints in turn);
2. **takes a plan** from the archive;
3. **breaks and rebuilds** part of it: removes 10–35 % of nearby customers and inserts them back;
4. **polishes** the plan with small moves: move a customer, reverse part of a route, swap vehicle type, move a trip;
5. **keeps it** if no plan in the archive beats it on all three objectives.

Local search alone gets stuck in plans that no small move improves (*local optima*). Two things get it out: the **break-and-rebuild** step, and **changing the viewpoint**, which turned a stuck plan into an improvable one in 56.5 % of switches. Our [ablation study](docs/REFERENCE.md#ablation-configurations-mpzis-2026) shows both are essential, while auto-tuning the break size gives nothing.

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

**У чому складність.** Машини різняться не лише розміром: вантажівка об’їжджає центр, велосипед їде пішохідними зв’язками; вага й об’єм — окремі обмеження; кожна машина може зробити кілька рейсів. І найкращого плану немає: дешевше зазвичай означає більше викидів або пізніше завершення.

**Як працює алгоритм.** Ітерований локальний пошук з архівом найкращих компромісів. У кожному циклі він обирає «точку зору» (ваги критеріїв), бере план з архіву, ламає й відбудовує його частину, доводить дрібними ходами і зберігає, якщо план не гірший за інші. Злам і зміна точки зору допомагають не застрягати в локальних оптимумах.

**Як спробувати:**

```bash
pip install -r requirements.txt
python3 example.py dataset/generated/R25_0.json --time 2
```

**Що де лежить:** `mt/` — сам метод, `dataset/` — тестові задачі, `results/` — результати з наших статей, `docs/REFERENCE.md` — як повторити всі експерименти.
