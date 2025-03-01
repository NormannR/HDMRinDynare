# %%
import Pkg; Pkg.activate("./venvs/sparsegrids")
using Dynare
using Tasmanian
using Plots, LaTeXStrings
using Statistics
# %%
# RBC Model
context = dynare("rbc.mod", "stoponerror");
# %%
(SG_grid, sgws) = sparsegridapproximation(tol_ti=1e-6,gridDepth=3,maxRef=0);
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
θ = (1-α)*(1-ss_l)*ss_y/(ss_l*ss_c);
# %%
c_pol(K,Z) = (1-α*β)*exp(Z)*K^α*ss_l^(1-α)
k_pol(K,Z) = α*β*exp(Z)*K^α*ss_l^(1-α)
y_pol(K,Z) = c_pol(K,Z)+k_pol(K,Z)
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
# Test of the user-provided initial policy guess
initialPolGuess = UserPolicyGuess(
    function (x)
        Z = x[1]
        K = x[2]
        return [
            c_pol(K,Z),
            k_pol(K,Z),
            y_pol(K,Z),
            α*y_pol(K,Z)/K
        ]
    end,
    ["z", "k"],
    ["c","k","y","rk"]
)
# %%
(SG_grid, sgws) = sparsegridapproximation(tol_ti=1e-6,gridDepth=3,maxRef=0,initialPolGuess=initialPolGuess,show_trace=false);
# %%
# IRBC model
# It's better to initialize the context variable out of a loop
context = dynare("irbc_small", "-DN=2","stoponerror");
# %%
# Reproduce Figure 8 in "USING ADAPTIVE SPARSE GRIDS TO SOLVE HIGH-DIMENSIONAL
# DYNAMIC MODELS", ECTA 2017, LHS panel
nb_points = Vector{Any}()
errors = Vector{Any}()
for l in [3,5,7]
   println(l)
   SG_grid, sgws = sparsegridapproximation(context=context, scaleCorrExclude=["lambda"], surplThreshold=0., gridDepth=l, maxRef=0, tol_ti=1e-7, ftol=1e-8, maxiter=400)
   errorMat = simulation_approximation_error!(context=context,grid=SG_grid,sgws=sgws)
   push!(nb_points, getNumPoints(SG_grid))
   push!(errors, errorMat)
end
# %%
# Check Euler Equations residuals only!
avg_errors = [ mean(abs.(sim[1:end-1, 1000:end,:])) for sim in errors]
q_errors = [ quantile(abs.(vec(sim[1:end-1, 1000:end,:])), 0.999) for sim in errors]
# %%
nb_points
# %%
log10.(avg_errors)
# %%
log10.(q_errors)
# %%
# Reproduce Figure 8 in "USING ADAPTIVE SPARSE GRIDS TO SOLVE HIGH-DIMENSIONAL
# DYNAMIC MODELS", ECTA 2017, RHS panel
errors =  Vector{Any}()
nb_points = Vector{Any}()
avg_time = Vector{Any}()
for N in [2,4,8]
   println(N)
   context = dynare("irbc_small", "-DN=$N", "stoponerror");
   SG_grid, sgws = sparsegridapproximation(context=context, scaleCorrExclude=["lambda"], surplThreshold=0., gridDepth=3, maxRef=0, tol_ti=1e-7, ftol=1e-8, maxiter=100)
   errorMat = simulation_approximation_error!(context=context,grid=SG_grid,sgws=sgws)
   results = context.results.model_results[1].sparsegrids
   push!(nb_points, getNumPoints(SG_grid))
   push!(errors, errorMat)
   push!(avg_time, results.average_iteration_time)
end
# %%
# Check Euler Equations residuals only!
avg_errors = [ mean(abs.(sim[1:end-1, 1000:end,:])) for sim in errors]
q_errors = [ quantile(abs.(vec(sim[1:end-1, 1000:end,:])), 0.999) for sim in errors]
# %%
nb_points
# %%
log10.(avg_errors)
# %%
log10.(q_errors)
# %%
# 1) Graph of the number of needed iterations to reach a tol_ti in function of the
# dimension
# 2) Graph of average iteration time for a single thread w.r.t # Dimensions
# 3) N=2, single thread, compute the number of needed iterations: compare first-order linear guess vs the full-zero policy guess
# DDSG: show the time spent is higher, but allows to reach higher dimensions than ASG or SG