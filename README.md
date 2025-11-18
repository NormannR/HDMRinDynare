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

## 📈 Replication

This repository reproduces all figures and tables from the paper *"Scalable Global Methods for Dynare"*. Below is a list of the replication scripts and their corresponding outputs:

- **`AnalyticalRBC.jl`**
  - Compares analytical solution with SG and DDSG approximations on the standard RBC model.

- **`AnalyticalDD.jl`**
  - Reproduces **Figure 6**
  - Benchmarks SG and DDSG on the high-dimensional function \( \left(\sum_i \sin(x_i)\right)^c \).

- **`IRBC_SG.jl`**
  - Reproduces **Tables 1–3**
  - Evaluates sparse grid performance on the IRBC model: error metrics, scalability, and iteration behavior.

- **`IRBC_DD.jl`**
  - Reproduces **Tables 4–5**
  - Benchmarks DDSG on the IRBC model and compares different initialization strategies.

- **`IRBC_SO_COMPARISON.jl`** and **so_ee_err.m**
  - Compare the Euler equation errors (average and 99.9% quantile) of the second-order-perturbation solution and the sparse-grid solution using simulation.
  - Uses model files `irbc_small_inc`, `irbc_small_so.mod` and `irbc_small_gl.mod`.
  - `IRBC_SO_COMPARISON.jl` compute the Euler equation errors (average and 99.9% quantile) of the sparse-grid solution using simulation. Importantly, it also saves the interpolation weights and nodes in `quadrature.mat`. The provided `quadrature.mat` file is for `N=2` and depth level 3. Note that if you change the number of dimensions of the IRBC model (`N` in `irbc_small_inc`) or the depth level of the sparse-grid solution (`l` in `IRBC_SO_COMPARISON.jl`), you must execute `IRBC_SO_COMPARISON.jl` first to get proper results in `so_ee_err.m`. 
  - `so_ee_err.m` compute the Euler equation errors (average and 99.9% quantile) of the second-order-perturbation solution using simulation. Importantly, it relies on the MATLAB/Octave version of Dynare (v6.4). 

- **`RBC.jl`** and **rbc_pol.m**
  - Compare the Euler equation errors (average and 99.9% quantile) of the second-order-perturbation solution and the sparse-grid solution.
  - Uses model files `rbc_inc.mod`, `rbc_so.mod` and `rbc_gl.mod`.
  - `RBC.jl` compute the Euler equation errors (average and 99.9% quantile) of the sparse-grid solution (i) on the state space using low-discrepancy sequences and (ii) sequentially using simulation. You must execute `rbc_pol.m` first.
  - `rbc_pol.m` computes the first-order and second-order perturbation solutions using Dynare and modifies the policy functions to take $(k_{t-1},z_t)$ as input (as in Schmitt-Grohé and Uribe (2004) JEDC) instead of $(k_{t-1},z_{t-1},e_t)$. The results are stored in `rbc_pol.mat`.

All results are computed from scratch using the `.mod` files and the provided code.

## 🛠 Setup

To replicate the results, clone the repository and navigate to its location. Julia v1.11 or more recent is necessary for the code to run smoothly. Launch Julia in a terminal and run:
```julia
julia>import Pkg
julia>Pkg.activate(".")
julia>Pkg.instantiate()
```
This setup step only needs to be done once. You can then run any replication script using:
```julia
julia>include("filename.jl")
```
or directly from a terminal:
```
$ julia --project=. filename.jl
```

## ✅ To-Do List

- [ ] In `SGapproximation` and `DDSGapproximation`, allow a warm-start from pre-existing
    - `TasmanianSG`/`DDSG` instance
    - `SparsegridsWs` instance
    - `Context` instance
- [ ] Allow user-defined time-iteration trajectories for `SGapproximation`. An example of such a trajectory is 20 iterations with `l=3`, 10 with `l=4`, 1 with `l=5`.
- [ ] Parallelize `simulation_approximation_error`
