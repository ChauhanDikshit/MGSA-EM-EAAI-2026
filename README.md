# MGSA-EM-EAAI-2026

## MGSA-EM (Restart mechanism-based Multilevel GSA) — MATLAB

This repository contains a MATLAB implementation of **MGSA-EM**: a *restart mechanism-based multilevel gravitational search algorithm* (GSA) with competition–collaboration learning and an enhanced restart strategy.

It also includes an **IDE-EDA** (Improved Differential Evolution + Estimation of Distribution Algorithm) module used as a restart operator when stagnation is detected.

---

## Reference Paper (MGSA-EM)

**D. Chauhan**, *Restart mechanism-based multilevel gravitational search algorithm for global optimization and image segmentation*, **Engineering Applications of Artificial Intelligence**, 163 (2026) 112904.  
DOI: `10.1016/j.engappai.2025.112904`

---

## What’s inside

- `MGSA_EM.m`  
  Main MGSA-EM optimizer (multilevel structure + competition/collaboration updates + stagnation handling).
- `IDE_EDA_pop.m` (or `IDE_EDA.m`)  
  Population-level restart operator (DE/current-to-pbest + archive + optional EDA sampling).
- Benchmark wrappers (e.g., CEC2017) depending on your setup.

---

**Example layout**
- `figures/overview.png`
- `figures/restart_flow.png`
- `figures/convergence.png`

**Embedded figures**
<p align="center">
  <img src="figures/overview.png" width="800">
</p>
<p align="center">
  <em>Figure 1. MGSA-EM overview (multilevel learning + restart).</em>
</p>

---

## Core idea (high level)

MGSA-EM:
1. **Sorts the population** by fitness and partitions it into **multiple layers** (top → bottom).
2. Uses **competitive learning** inside each layer: individuals are paired, producing **winners** and **losers**.
3. Updates:
   - **Losers** learn from winners (+ acceleration term / pbest term).
   - **Winners** learn from upper-layer individuals (+ pbest term and cross-layer guidance).
4. Tracks **stagnation** using a counter per individual (increment if no improvement, reset when improved).
5. When stagnation is severe, a **restart mechanism** is triggered using:
   - **Differential mutation** (DE/current-to-pbest/1), optionally using an **archive**,
   - plus an optional **EDA sampling** step to inject diversity.

---

## Restart mechanism: how it’s triggered (important)

**Stagnation counters** are maintained per individual. When an individual’s counter exceeds the stagnation threshold `Sg`, the restart mechanism is activated (Algorithm 5–6 in the paper).

In code, you may use either:

### Option A — Individual-level restart (recommended for “only stagnated particles”)
Call the restart operator **only for individuals** with `counter(i) > Sg`.

### Option B — Population-level restart (what you implemented)
Trigger IDE-EDA for the whole population only when a large fraction is stagnated, e.g.:
```matlab
stagRate = mean(counter > Sg);
if stagRate > 0.5
    [fitness,p,...] = IDE_EDA_pop(...);
    counter(:) = 0;
end


<p align="center">
  <img src="figures/layers.png" width="800">
</p>
<p align="center">
  <em>Figure 2. Heat map of ranks on benchmark problems of the selected topological structures at 50 dimensions.</em>
</p>

<p align="center">
  <img src="figures/omega.png" width="800">
</p>
<p align="center">
  <em>Figure 3. Heat map of ranks on benchmark problems for each value of omega at 50 dimensions.</em>
</p>

<p align="center">
  <img src="figures/stagnation.png" width="800">
</p>
<p align="center">
  <em>Figure 4. Heat map of ranks on benchmark problems of each value of 𝑆𝑔 at 50 dimensions.</em>
</p>

<p align="center">
  <img src="figures/algorithms.png" width="800">
</p>
<p align="center">
  <em>Figure 5. The comparison of MGSA-EM’s dimension-wise average ranks with GSA variants and other state-of-the-art algorithms.</em>
</p>

---

### MGSA-EM
@article{chauhan2026restart,
  title={Restart mechanism-based multilevel gravitational search algorithm for global optimization and image segmentation},
  author={Chauhan, Dikshit},
  journal={Engineering Applications of Artificial Intelligence},
  volume={163},
  pages={112904},
  year={2026},
  doi = {10.1016/j.engappai.2025.112904},
  publisher={Elsevier}
}

### IDE-EDA (Restart Operator)
@article{li2023improved,
  title={An improved differential evolution by hybridizing with estimation-of-distribution algorithm},
  author={Li, Yintong and Han, Tong and Tang, Shangqin and Huang, Changqiang and Zhou, Huan and Wang, Yuan},
  journal={Information sciences},
  volume={619},
  pages={439--456},
  year={2023},
  doi = {10.1016/j.ins.2022.11.029},
  publisher={Elsevier}
}

