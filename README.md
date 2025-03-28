# Scalable Global Methods for Dynare

This repository accompanies the paper [**Scalable Global Solution Techniques for High-Dimensional Models in Dynare**](https://arxiv.org/abs/2503.11464)

It provides Julia code to reproduce the figures and tables from the paper, using sparse grid and DDSG methods implemented in `Dynare.jl`.

## 📁 Project Structure

```
├── AnalyticalRBC.jl       # Analytical benchmark for the RBC model
├── AnalyticalDD.jl        # High-dimensional function benchmark (sin^c)
├── IRBC_SG.jl             # Sparse Grid solution for the IRBC model
├── IRBC_DD.jl             # DDSG solution for the IRBC model
├── rbc.mod                # Dynare mod file for RBC model
├── irbc_small.mod         # Dynare mod file for IRBC model
├── Project.toml           # Julia Project.toml file
├── Manifest.toml          # Julia Manifest.toml file
```

## 🛠 Setup

```julia
import Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## ✅ To-Do List

- [ ] In `SGapproximation` and `DDSGapproximation`, allow a warm-start from pre-existing
    - `TasmanianSG`/`DDSG` instance
    - `SparsegridsWs` instance
    - `Context` instance
- [ ] Allow user-defined time-iteration trajectories for `SGapproximation`. An example of such a trajectory is 20 iterations with `l=3`, 10 with `l=4`, 1 with `l=5`.
