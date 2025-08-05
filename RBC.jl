# %%
import Pkg; Pkg.activate(".")
using Revise
using Dynare
using Plots, LaTeXStrings
using QuasiMonteCarlo, Distributions
using MAT
using LinearAlgebra
using StaticArrays
using Statistics, StatsBase
using DataFrames
using CategoricalArrays
using VectorizedStatistics
# %%
# RBC Model
context = dynare("rbc_gl.mod", "stoponerror");
# %%
params = context.work.params;
symbol_table = context.symboltable;
# %%
# --- Parameters ---
const α = params[symbol_table["alpha"].orderintype]
const β = params[symbol_table["beta"].orderintype]
const ρ = params[symbol_table["rho"].orderintype]
const σ = params[symbol_table["sig"].orderintype]
const ss_l = params[symbol_table["ss_l"].orderintype]
const ss_k = params[symbol_table["ss_k"].orderintype]
# %%
# --- True Policy Functions ---
c_pol(K, Z) = (1 - α * β) * exp(Z) * K^α * ss_l^(1 - α)
k_pol(K, Z) = α * β * exp(Z) * K^α * ss_l^(1 - α)
y_pol(K, Z) = c_pol(K, Z) + k_pol(K, Z)
rk(K, Z) = α*y_pol(K, Z)/K
# %%
# --- Perturbation Policy Functions ---
vars = matread("rbc_pol.mat")
# %%
row_c = Integer(vars["row_c"]) # index of consumption in the Dynare policy function vector
row_rk = Integer(vars["row_rk"]) # index of consumption in the Dynare policy function vector
row_k = Integer(vars["row_k"])
g_0 = vars["g_0"]
g_1 = vars["g_1"]
g_2 = vars["g_2"]
ys_order_var = vars["ys_order_var"]
ys_order_var_k2 = vars["ys_order_var_k2"]
pol_1st(x) = ys_order_var+g_1*(x-ys_order_var_k2);
pol_2nd(x) = g_0+g_1*(x-ys_order_var_k2)+g_2*kron(x-ys_order_var_k2,x-ys_order_var_k2);
# %%
k_min = context.work.limits[:k].min
k_max = context.work.limits[:k].max
z_min = context.work.limits[:z].min
z_max = context.work.limits[:z].max
# %%
l = 3
sg_options = SGOptions(gridDepth=l, tol_ti=1e-7, ftol=1e-8, maxiter=2000, polUpdateWeight=1.)
SG_L3, sgws_L3 = SGapproximation(sg_options; context=context);
# %%
l = 5
sg_options = SGOptions(gridDepth=l, tol_ti=1e-7, ftol=1e-8, maxiter=2000, polUpdateWeight=1.)
SG_L5, sgws_L5 = SGapproximation(sg_options; context=context);
# %%
l = 7
sg_options = SGOptions(gridDepth=l, tol_ti=1e-7, ftol=1e-8, maxiter=2000, polUpdateWeight=1.)
SG_L7, sgws_L7 = SGapproximation(sg_options; context=context);
# %%
l = 3
ddsg_opts = DDSGOptions(gridDepth=l, tol_ti=1e-7, ftol=1e-8, polUpdateWeight=1.);
DDSG_L3, sgws_L3 = DDSGapproximation(ddsg_opts; context=context);
# %%
l = 5
ddsg_opts = DDSGOptions(gridDepth=l, tol_ti=1e-7, ftol=1e-8, polUpdateWeight=1.);
DDSG_L5, sgws_L5 = DDSGapproximation(ddsg_opts; context=context);
# %%
l = 7
ddsg_opts = DDSGOptions(gridDepth=l, tol_ti=1e-7, ftol=1e-8, polUpdateWeight=1.);
DDSG_L7, sgws_L7 = DDSGapproximation(ddsg_opts; context=context);
# %%
# Index of consumption in the (DD)SG policy function vector
row_c_sg = findfirst(context.models[1].i_dyn .== context.symboltable["c"].orderintype)
row_rk_sg = findfirst(context.models[1].i_dyn .== context.symboltable["rk"].orderintype)
row_k_sg = findfirst(context.models[1].i_dyn .== context.symboltable["k"].orderintype)
# %%
dynare_pol = Dict(
   "1st"        => pol_1st,
   "2nd"        => pol_2nd,
)
sg_pol = Dict(
   "SG-L3"      => x->Dynare.interpolate(SG_L3, x),
   "SG-L5"      => x->Dynare.interpolate(SG_L5, x),
   "SG-L7"      => x->Dynare.interpolate(SG_L7, x),
   "DDSG-L3"    => x->Dynare.interpolate(DDSG_L3, x),
   "DDSG-L5"    => x->Dynare.interpolate(DDSG_L5, x),
   "DDSG-L7"    => x->Dynare.interpolate(DDSG_L7, x)
)
# %%
# Plot the consumption policy functions
col_true        = :black                     # analytic benchmark
linestyle_sg    = [:dash, :dash, :dash, :dot, :dot, :dot]
palette_dynare = [:black, :red]
palette_sg = [:blue, :green, :orange, :purple, :cyan, :magenta]
# %%
σz = σ/sqrt(1 - ρ^2)
z_slices = [0.0 -2*σz  +2*σz -3*σz +3*σz]
nk = 200
k_grid = range(k_min, k_max; length=nk)
# %%
nk = 200
k_grid = range(k_min, k_max; length=nk)
p = plot(layout = (1, length(z_slices)),
         legend = :bottomright, size = (1000, 350))
for (jz, z_val) in enumerate(z_slices)
    c_vals = [c_pol(k, z_val) for k in k_grid]
    plot!(p, k_grid, c_vals; label = "True", lw = 1,
            color = :black, subplot = jz)
    for ((name, f), col) in zip(collect(dynare_pol), palette_dynare)
        c_vals = [f([k, z_val])[row_c] for k in k_grid]
        plot!(p, k_grid, c_vals; label = name, lw = 2,
              color = col, subplot = jz)
    end
    for ((name, f), col) in zip(collect(sg_pol), palette_sg)
        c_vals = [f([k, z_val])[row_c_sg] for k in k_grid]
        plot!(p, k_grid, c_vals; label = name, lw = 2,
              color = col, linestyle = :dash, subplot = jz)
    end
    xlabel!(p[jz], "K")
    ylabel!(p[jz], "C")
    title!(p[jz], "Z = $(round(z_val, digits = 3))")
end
display(p)      # show the composite figure
# %%
# Plot the consumption policy functions at Z=0.0
c_vals = [c_pol(k, 0.) for k in k_grid]
p = plot(k_grid, c_vals; label = "True",
        color = :black, lw=2)
for ((name, f), col) in zip(collect(dynare_pol), palette_dynare)
    c_vals = [f([k, 0.])[row_c] for k in k_grid]
    plot!(p, k_grid, c_vals; label = name,
            color = col)
end
for ((name, f), col, lin) in zip(collect(sg_pol), palette_sg, linestyle_sg)
    c_vals = [f([k, 0.])[row_c_sg] for k in k_grid]
    plot!(p, k_grid, c_vals; label = name, lw = 2,
            color = col, linestyle = lin)
end
vline!(p, [0.8435349564205504*ss_k, 1.1950389406183832*ss_k], lw=1, linestyle=:dash, color=:black, label="")
xlabel!(L"k_{t-1}")
ylabel!(L"c_t")
# %%
# Plot the consumption policy functions at Z=0.0
c_vals = [c_pol(k, 0.) for k in k_grid]
p = plot(k_grid, c_vals; label = "True",
        color = :black, lw=2)
# %%
for ((name, f), col) in zip(collect(dynare_pol), palette_dynare)
    c_vals = [f([k, 0.])[row_c] for k in k_grid]
    plot!(p, k_grid, c_vals; label = name,
            color = col)
end
# %%
f = x->Dynare.interpolate(SG_L3, x)
c_vals = [f([k, 0.])[row_c_sg] for k in k_grid]
plot!(p, k_grid, c_vals; label = "SG-L3",
        color = :blue, linestyle=:dash)
# %%
f = x->Dynare.interpolate(DDSG_L3, x)
c_vals = [f([k, 0.])[row_c_sg] for k in k_grid]
plot!(p, k_grid, c_vals; label = "DDSG-L3",
        color = :purple, linestyle=:dot)
# %%
vline!(p, [0.8435349564205504*ss_k, 1.1950389406183832*ss_k], lw=1, linestyle=:dash, color=:black, label="")
xlabel!(L"k_{t-1}")
ylabel!(L"c_t")
# %%
# Simulations:
# k_min = 0.8435349564205504*ss_k
# k_max = 1.1950389406183832*ss_k
# z_max = 3.8623283809251028*σz
# z_min = -3.401333795720321*σz
savefig("pol_c.pdf")
# %%
# Euler equation error on the domain
# %%
# Sobol in the (k,z) state-space
N  = 10_000
lb = [k_min, z_min]
ub = [k_max, z_max]
KZ = QuasiMonteCarlo.sample(N, lb, ub, SobolSample())
# %%
# Integration nodes
M = 1024
qmcs = QuasiMonteCarlo.sample(M, 1, SobolSample())
ϵ = quantile.(Normal(), qmcs)
# %%
make_wrapper(f, iC, iK, iRK) =
    x -> begin
        y = f(x)
        (y[iC],  y[iK],  y[iRK])
    end
# %%
function euler_residuals(policy)::Vector{Float64}
    # --- evaluate (c_t, k_{t+1}) on all N states
    X = [Vector(s) for s in eachslice(KZ,dims=2)]
    ck_rk = policy.(X)
    c_t   = first.(ck_rk)
    k_tp1 = getindex.(ck_rk, 2)

    # --- broadcast shocks → matrices (M × N)
    z_tp1 = ρ .* reshape(KZ[2,:], 1, :) .+ σ .* reshape(ϵ, :, 1)
    k_tp1_mat = zeros(size(z_tp1)) .+ reshape(k_tp1,1,:)
    X_tp1 = vcat(reshape(k_tp1_mat,1,:), reshape(z_tp1,1,:))
    X_tp1 = [Vector(s) for s in eachslice(X_tp1,dims=2)]

    ck_rk_next = policy.(X_tp1)
    c_tp1   = first.(ck_rk_next)
    rk_tp1  = getindex.(ck_rk_next, 3)

    c_tp1_mat  = reshape(c_tp1,  M, :)
    rk_tp1_mat = reshape(rk_tp1, M, :)

    # --- Euler equation errors
    rhs = β .* mean(rk_tp1_mat ./ c_tp1_mat; dims=1) |> vec
    lhs = 1 ./ c_t
    return abs.(lhs .- rhs)
end
# %%
stats = DataFrame(Method = String[], Mean = Float64[], P999 = Float64[])
euler_state_space = Dict()
# %%
# Check with the true policy function
function true_pol(x)
   K = x[1]
   Z = x[2]
   return (c_pol(K, Z), k_pol(K, Z), rk(K, Z))
end
# %%
e = euler_residuals(true_pol)
push!(stats, ("True",
              log10(mean(e)),
              log10(quantile(e, 0.999))))
# %%
for (name, pol) in dynare_pol
    e = euler_residuals(make_wrapper(pol, row_c, row_k, row_rk))
    euler_state_space[name] = e
    push!(stats, (name,
                  log10(mean(e)),
                  log10(quantile(e, 0.999))))
end
# %%
for (name, pol) in sg_pol
    euler_state_space[name] = e
    e = euler_residuals(make_wrapper(pol, row_c_sg, row_k_sg, row_rk_sg))
    push!(stats, (name,
                  log10(mean(e)),
                  log10(quantile(e, 0.999))))
end
# %%
# Scatter plot of errors
scatter(KZ[1, :], KZ[2, :];
        marker_z   = log10.(e),     # colour = residual
        ms         = 3,             # marker size
        colorbar   = true,
        colorbar_title = "log₁₀ |ε|",
        xlabel     = "capital k",
        ylabel     = "productivity z",
        title      = "Euler-equation residuals on Sobol nodes",
        legend     = false)
# %%
for (name,ee) in euler_state_space
   scatter(KZ[1, :], KZ[2, :];
         marker_z   = log10.(ee),
         ms         = 3,
         colorbar   = true,
         colorbar_title = "log₁₀ |ε|",
         xlabel     = "K",
         ylabel     = "Z",
         title      = "Euler-equation residuals: $(name)",
         legend     = false)
   savefig("euler_state_space_$(name).pdf")
end
# %%
# Heatmap of errors
# Build the bin edges
nk, nz = 45, 45
k_edges = range(k_min, k_max; length = nk + 1)
z_edges = range(z_min, z_max; length = nz + 1)
# Assign every Sobol node to a (k,z) cell
k_bins = cut(KZ[1, :], k_edges)
z_bins = cut(KZ[2, :], z_edges)
k_bins = levelcode.(k_bins)
z_bins = levelcode.(z_bins)
# Accumulate residuals per cell
e_sum    = zeros(nz, nk)
e_counts = zeros(Int, nz, nk)
for j in eachindex(e)
    ik = k_bins[j];  iz = z_bins[j]
    if !(ismissing(ik) || ismissing(iz))
        e_sum[iz, ik]    += e[j]
        e_counts[iz, ik] += 1
    end
end

e_avg = e_sum ./ max.(e_counts, 1)

# Plot log10 residuals as a heat-map
k_centers = (k_edges[1:end-1] .+ k_edges[2:end]) ./ 2
z_centers = (z_edges[1:end-1] .+ z_edges[2:end]) ./ 2

heatmap(k_centers,
        z_centers,
        log10.(e_avg);       # colour scale
        xlabel = "capital k",
        ylabel = "productivity z",
        colorbar_title = "log₁₀ |ε|",
        title = "Euler-equation residuals")
# %%
# Simulated Euler equation residuals
# %%
function simulated_euler_residuals(policy,e) 
   rep = size(e,1)
   T = size(e,2)
   ee_sim = zeros(rep,T)
   X = Matrix{Vector{Float64}}(undef,rep,T)
   for r = 1:rep
      X[r,1] = [ss_k,0]
   end
   for t=2:T
      # --- evaluate (c_t, k_{t+1})
      k = getindex.(X[:,t-1],1)
      z = ρ .* getindex.(X[:,t-1],2) .+ σ*e[:,t-1]
      ck_rk = policy.([collect(x) for x in zip(k,z)])
      c_t   = getindex.(ck_rk, 1)
      k_tp1 = getindex.(ck_rk, 2)
      X[:,t] = [collect(x) for x in zip(k_tp1,z)]

      # --- broadcast shocks → Vector (M × 1)
      z_tp1 = ρ * z .+ σ .* reshape(ϵ,1,:)
      k_tp1_mat = zeros(size(z_tp1)) .+ k_tp1
      X_tp1 = vcat(reshape(k_tp1_mat,1,:), reshape(z_tp1,1,:))
      X_tp1 = [Vector(s) for s in eachslice(X_tp1,dims=2)]

      ck_rk_next = policy.(X_tp1)
      c_tp1   = first.(ck_rk_next)
      rk_tp1  = getindex.(ck_rk_next, 3)
      # %%
      c_tp1_mat  = reshape(c_tp1,  rep, :)
      rk_tp1_mat = reshape(rk_tp1, rep, :)
      # %%
      # --- Euler equation errors
      rhs = β * mean(rk_tp1_mat ./ c_tp1_mat, dims=2)
      lhs = 1 ./ c_t
      ee_sim[:,t] = abs.(lhs .- rhs)
   end
   return (ee_sim, X)
end
# %%
euler_simulated = Dict()
stats_sim = DataFrame(Method = String[], Mean = Float64[], P999 = Float64[])
rep = 1
burn = 1000
T = 10000
e = randn(rep,T+burn)
# %%
eer,X = simulated_euler_residuals(true_pol,e)
euler_simulated["True"] = (eer,X)
q_eer = vquantile!(copy(eer[:,burn+1:T]), 0.999, dims=(2))
push!(stats_sim, ("True",
            log10(mean(eer[:,burn+1:T])),
            log10(mean(q_eer))))
# %%
for (name, pol) in dynare_pol
    eer,X = simulated_euler_residuals(make_wrapper(pol, row_c, row_k, row_rk),e)
    euler_simulated[name] = (eer,X)
    q_eer = vquantile!(copy(eer[:,burn+1:T]), 0.999, dims=(2))
    push!(stats_sim, (name,
                  log10(mean(eer[:,burn+1:T])),
                  log10(mean(q_eer))))
end
# %%
for (name, pol) in sg_pol
    eer,X = simulated_euler_residuals(make_wrapper(pol, row_c_sg, row_k_sg, row_rk_sg),e)
    euler_simulated[name] = (eer,X)
    q_eer = vquantile!(copy(eer[:,burn+1:T]), 0.999, dims=(2))
    push!(stats_sim, (name,
                  log10(mean(eer[:,burn+1:T])),
                  log10(mean(q_eer))))
end
# %%
# Checking state bounds for the simulations
k_max = maximum([maximum([x[1] for x in v[2]]) for v in values(euler_simulated)])
k_min = minimum([minimum([x[1] for x in v[2]]) for v in values(euler_simulated)])
z_max = maximum([maximum([x[2] for x in v[2]]) for v in values(euler_simulated)])
z_min = minimum([minimum([x[2] for x in v[2]]) for v in values(euler_simulated)])
# %%
k_max/ss_k
# %%
k_min/ss_k
# %%
z_max/σz
# %%
z_min/σz
# %%
