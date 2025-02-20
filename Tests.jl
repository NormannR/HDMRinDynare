# %%
import Pkg; Pkg.activate("./venvs/sparsegrids")
using Revise, Dynare
using Tasmanian
using Plots, LaTeXStrings
using Statistics
# %%
# RBC Model
context = dynare("rbc.mod", "stoponerror");
# %%
(SG_grid, sgws) = sparsegridapproximation(tol_ti=1e-6,gridDepth=7);
# %%
α     = 1/3;
β     = 0.99;
ρ     = 0.95;
σ     = 0.01;
ss_l  = 1/3;
ss_k  = ss_l*(β*α)^(1/(1-α));
ss_y  = ss_k^α*ss_l^(1-α);
ss_c  = ss_y-ss_k;
ss_rk = α*ss_y/ss_k;
theta = (1-α)*(1-ss_l)*ss_y/(ss_l*ss_c);
# %%
c_pol(K,Z) = (1-α*β)*exp(Z)*K^α*ss_l^(1-α);
k_pol(K,Z) = α*β*exp(Z)*K^α*ss_l^(1-α);
# %%
Z_vals = range(-2*σ/sqrt(1-ρ^2), 2*σ/sqrt(1-ρ^2), length=100)
K_vals = range(ss_k * 0.1, ss_k * 1.9, length=100) 
# %%
c_theoretical = [c_pol(K, Z) for Z in Z_vals, K in K_vals]
k_theoretical = [k_pol(K, Z) for Z in Z_vals, K in K_vals]
# %%
pol = x->evaluate(SG_grid, x)
# %%
# Declaration order in MATLAB
i_c = context.symboltable["c"].orderintype
i_k = context.symboltable["k"].orderintype
# Retrieve the order in the policy approximation
i_c = findall(i_c .== context.models[1].i_dyn)[1]
i_k = findall(i_k .== context.models[1].i_dyn)[1]
# %%
c_approx = [pol([K, Z])[i_c] for Z in Z_vals, K in K_vals]
k_approx = [pol([K, Z])[i_k] for Z in Z_vals, K in K_vals]
# %%
# Plot Consumption Policy Function for a mid-range value of Z
p1 = plot(K_vals, c_theoretical[50, :], label="Theoretical", title="Consumption Policy")
plot!(K_vals, c_approx[50, :], label="Approximated", linestyle=:dash)
ylabel!(L"C_t \mid Z_t = Z")

# Plot Capital Policy Function
p2 = plot(K_vals, k_theoretical[50, :], label="Theoretical", title="Capital Policy")
plot!(K_vals, k_approx[50, :], label="Approximated", linestyle=:dash)
ylabel!(L"K_t \mid Z_t = Z")

plot(p1, p2, layout=(1, 2))
xlabel!(L"K_{t-1}")
# %%
# Plot residuals as heatmaps
c_residuals = abs.(c_theoretical .- c_approx)
k_residuals = abs.(k_theoretical .- k_approx)
p1 = heatmap(K_vals, Z_vals, c_residuals', xlabel=L"K_{t-1}", ylabel=L"Z_t", title="Consumption Residuals")
p2 = heatmap(K_vals, Z_vals, k_residuals', xlabel=L"K_{t-1}", ylabel=L"Z_t", title="Capital Residuals")
plot(p1, p2, layout=(1, 2))
# %%
# IRBC model
context = dynare("irbc_small", "stoponerror");
# %%
# Reproduce Figure 8 in "USING ADAPTIVE SPARSE GRIDS TO SOLVE HIGH-DIMENSIONAL
# DYNAMIC MODELS", ECTA 2017, LHS panel
results_lhs = Vector{Any}()
errors_lhs =  Vector{Any}()
for l in [3,5,7,9]
   println(l)
   SG_grid, sgws = sparsegridapproximation(scaleCorrExclude=["lambda"], gridDepth=l)
   errors = simulation_approximation_error!(context=context,grid=SG_grid,sgws=sgws)
   push!(results_lhs, (l, SG_grid, sgws))  # Store (refinement level, grid, workspace)
   push!(errors_lhs, (l, errors))
end
# %%
T = 11000
burnin = 1000
nb_points = [getNumPoints(result[2]) for result in results_lhs]
mean_abs_error = [mean(abs.(error[2][:,burnin+1:end])) for error in errors_lhs]
max_error = [maximum(abs.(error[2][:,burnin+1:end])) for error in errors_lhs]
p1 = plot(nb_points, max_error, label="Max. Error", xaxis=:log10, yaxis=:log10)
plot!(nb_points, mean_abs_error, label="Avg. Error", xaxis=:log10, yaxis=:log10)
xlabel!("# Points")
ylabel!("Error")
# %%
# Reproduce Figure 8 in "USING ADAPTIVE SPARSE GRIDS TO SOLVE HIGH-DIMENSIONAL
# DYNAMIC MODELS", ECTA 2017, RHS panel
results_rhs = Vector{Any}()
errors_rhs =  Vector{Any}()
N_vec = [2,10,25,50]
for N in N_vec
   println(N)
   context = dynare("irbc_small", "-DN=$N", "stoponerror");
   SG_grid, sgws = sparsegridapproximation(scaleCorrExclude=["lambda"], gridDepth=3)
   errors = simulation_approximation_error!(context=context,grid=SG_grid,sgws=sgws)
   push!(results_rhs, (2*N, SG_grid, sgws))  # Store (refinement level, grid, workspace)
   push!(errors_rhs, (2*N, errors))
end
# %%
mean_abs_error = [mean(abs.(error[2][:,burnin+1:end])) for error in errors_lhs]
max_error = [maximum(abs.(error[2][:,burnin+1:end])) for error in errors_lhs]
p2 = plot(2*N_vec, max_error, label="Max. Error", xaxis=:log10, yaxis=:log10)
plot!(2*N_vec, mean_abs_error, label="Avg. Error", xaxis=:log10, yaxis=:log10)
xlabel!("# Dimensions")
ylabel!("Error")
# %%
plot(p1, p2, layout=(1, 2))