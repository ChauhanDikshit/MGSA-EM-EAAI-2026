# MGSA-EM-EAAI-2026

## MGSA-EM (Restart mechanism-based Multilevel GSA) — MATLAB

This repository contains a MATLAB implementation of **MGSA-EM**: a *restart mechanism-based multilevel gravitational search algorithm* (GSA) with competition–collaboration learning and an enhanced restart strategy [1].

It also includes an **IDE-EDA** (**I**mproved **D**ifferential **E**volution + **E**stimation of **D**istribution **A**lgorithm) module used as a restart operator when stagnation is detected [2].

---

## Reference Paper (MGSA-EM)

**D. Chauhan**, *Restart mechanism-based multilevel gravitational search algorithm for global optimization and image segmentation*, **Engineering Applications of Artificial Intelligence**, 163 (2026) 112904 [1].  
DOI: `10.1016/j.engappai.2025.112904`

---

## What’s inside

- `MGSA_EM.m`  
  Main MGSA-EM optimizer (multilevel structure + competition/collaboration updates + stagnation handling) [1].
- `IDE_EDA_pop.m` (or `IDE_EDA.m`)  
  Population-level restart operator (**DE/current-to-pbest/1** + archive + optional EDA sampling) [2].
- Benchmark wrappers (e.g., CEC2017), depending on your setup.

---

## Requirements

- MATLAB R2019b or later (64-bit recommended)
- (Optional) Statistics and Machine Learning Toolbox for `mvnrnd` (only needed if EDA sampling is enabled) [2]

---

## Figures (README preview)

<p align="center">
  <img src="figures/overview.png" width="800">
</p>
<p align="center">
  <em>Figure 1. MGSA-EM overview (multilevel learning + restart).</em>
</p>

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
  <em>Figure 4. Heat map of ranks on benchmark problems for each value of <code>Sg</code> at 50 dimensions.</em>
</p>

<p align="center">
  <img src="figures/algorithms.png" width="800">
</p>
<p align="center">
  <em>Figure 5. Comparison of MGSA-EM’s dimension-wise average ranks with GSA variants and other state-of-the-art algorithms.</em>
</p>

---

## Core idea (high level)

MGSA-EM [1]:
1. **Sorts the population** by fitness and partitions it into **multiple layers** (top → bottom).
2. Uses **competitive learning** inside each layer: individuals are paired, producing **winners** and **losers**.
3. Updates:
   - **Losers** learn from winners (plus acceleration / guidance terms).
   - **Winners** learn from upper-layer individuals (cross-layer guidance + exploitation).
4. Tracks **stagnation** using a counter per individual (increment if no improvement, reset when improved).
5. When stagnation is severe, a **restart mechanism** is triggered using:
   - **Differential mutation** (DE/current-to-pbest/1), optionally using an **archive** [2],
   - plus an optional **EDA sampling** step to inject diversity [2].

---

## Restart mechanism: how it’s triggered (important)

**Stagnation counters** are maintained per individual. When an individual’s counter exceeds the stagnation threshold `Sg`, the restart mechanism is activated (Algorithm 5–6 in the paper) [1].

In code, you may use either:

### Option A — Individual-level restart (only stagnated particles)
Call the restart operator **only for individuals** with `counter(i) > Sg` [1,2].

### Option B — Population-level restart (your current implementation)
Trigger IDE-EDA for the whole population only when a large fraction is stagnated, e.g. [1,2]:

```matlab
stagRate = mean(counter > Sg);
if stagRate > 0.5
    [fitness, p, ...] = IDE_EDA_pop(...);
    counter(:) = 0;
end
