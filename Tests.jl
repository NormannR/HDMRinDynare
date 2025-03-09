# %%
import Pkg; Pkg.activate("./venvs/sparsegrids")
using Revise, Dynare
# using Dynare
using Tasmanian
using Plots, LaTeXStrings
using Statistics
# %%
# RBC Model
context = dynare("rbc.mod", "stoponerror");
# %%
(SG_grid, sgws) = sparsegridapproximation(tol_ti=1e-6,gridDepth=3,maxRef=0, polUpdateWeight=1.);
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
initialPolGuess = UserPolicyGuess(
    function (x)
        Z = x[1]
        K = x[2]
        return [
            0.15,
            0.05,
            0.20,
            1.
        ]
    end,
    ["z", "k"],
    ["c","k","y","rk"]
)
# %%
(SG_grid, sgws) = sparsegridapproximation(tol_ti=1e-6,gridDepth=3,maxRef=0,initialPolGuess=initialPolGuess,show_trace=false,polUpdateWeight=1.);
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
   SG_grid, sgws = sparsegridapproximation(context=context, scaleCorrExclude=["lambda"], surplThreshold=0., gridDepth=l, maxRef=0, tol_ti=1e-7, ftol=1e-8, maxiter=400, polUpdateWeight=1)
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
# - Reproduce Figure 8 in "USING ADAPTIVE SPARSE GRIDS TO SOLVE HIGH-DIMENSIONAL
#   DYNAMIC MODELS", ECTA 2017, RHS panel
# - Count the number of iterations necessary to reach a time-iteration step
#   lower than 1e-4 with a native initial guess for policy functions
errors =  Vector{Any}()
nb_points = Vector{Any}()
avg_time = Vector{Any}()
for N in [2,4,8]
   println(N)
   context = dynare("irbc_small", "-DN=$N", "stoponerror");
   SG_grid, sgws = sparsegridapproximation(context=context, scaleCorrExclude=["lambda"], surplThreshold=0., gridDepth=3, maxRef=0, tol_ti=1e-7, ftol=1e-8, maxiter=100, polUpdateWeight=1  )
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
# Count the number of iterations necessary to reach a time-iteration step lower
# than 1e-4 with a Dynare-provided first-order initial guess for policy
# functions
iteration_numbers =  Vector{Any}()
for N in [2,4,8]
    println(N)
    context = dynare("irbc_small", "-DN=$N", "stoponerror");
    SG_grid, sgws = sparsegridapproximation(context=context, scaleCorrExclude=["lambda"], gridDepth=3, maxRef=0, tol_ti=1e-4, ftol=1e-8, polUpdateWeight=1)
    results = context.results.model_results[1].sparsegrids
    push!(iteration_numbers, results.iteration_number)
end
# %%
# Count the number of iterations necessary to reach a time-iteration step lower
# than 1e-4 with a native initial guess for policy functions
iteration_numbers =  Vector{Any}()
for N in [2,4,8]
    println(N)
    context = dynare("irbc_small", "-DN=$N", "stoponerror");
    initialPolGuess = UserPolicyGuess(
        x -> ones(N+1),
        vcat(["k_$i" for i=1:N], ["a_$i" for i=1:N]),
        vcat(["lambda"],["k_$i" for i=1:N])
    )
    SG_grid, sgws = sparsegridapproximation(context=context, scaleCorrExclude=["lambda"], gridDepth=3, maxRef=0, tol_ti=1e-4, ftol=1e-8,
    initialPolGuess=initialPolGuess, polUpdateWeight=1)
    results = context.results.model_results[1].sparsegrids
    push!(iteration_numbers, results.iteration_number)
end
# %%
iteration_numbers
# %%
# DDSG
# Figure 6 LHS Graph: Plotting 
nb_ddsg = Vector{Int}()
errors_ddsg = Vector{Float64}()
nb_sg = Vector{Int}()
errors_sg = Vector{Float64}()
# %%
c = 1
num_points = 1000
F_poly(X::AbstractVector{Float64}; c)::Vector{Float64} = [sum(sin, X)^c]
F_poly(X::AbstractMatrix{Float64}; c)::Matrix{Float64} = reduce(hcat, [F_poly(col;c=c) for col in eachcol(X)])
F(X::Vector{Float64}) = F_poly(X, c=c)
F(X::Matrix{Float64}) = F_poly(X, c=c)
dim = 20
k_max = 1
dof = 1
X_sample = rand(dim, num_points)
Y_exact = F(X_sample)
for l in 1:5
    ddsg = Dynare.DDSG(dim, dof, l, l, k_max)
    Dynare.DDSG_init!(ddsg)
    Dynare.DDSG_build!(
        ddsg,
        F,
        ddsg.centroid
    )
    # %%
    Y_ddsg = Dynare.DDSG_evaluate(ddsg, X_sample)
    push!(errors_ddsg, mean(@. abs((Y_ddsg - Y_exact)/Y_exact)))
    push!(nb_ddsg, sum(ddsg.grid_points))
    # %%
    sg = Tasmanian.TasmanianSG(dim, dof, l)
    Tasmanian.makeLocalPolynomialGrid!(sg)
    Tasmanian.setDomainTransform!(sg, ddsg.domain)  # Match grid to domain dimensions
    # Refinement based on surplus coefficients
    X = Tasmanian.getNeededPoints(sg)
    N = size(X,2)
    Y_val = F(X)
    Tasmanian.loadNeededPoints!(sg, Y_val)
    Y_sg = Tasmanian.evaluateBatch(sg,X_sample)
    # %%
    push!(errors_sg, mean(@. abs((Y_sg - Y_exact)/Y_exact)))
    push!(nb_sg, Tasmanian.getNumPoints(sg))
end
# %%
lhs = plot(nb_sg, errors_sg, color = :black, label=L"SG",
     markershape=:star6, markersize=5, linewidth = 4)
plot!(nb_ddsg, errors_ddsg, color = :blue, label=L"DDSG~(\mathcal{K} = 1)",
      markershape=:circle, markersize=5, linewidth = 4)
plot!(xscale=:log10, yscale=:log10, minorgrid=true)
xlabel!("Number of Grid Points")
ylabel!("Error")
# %%
# Figure 6 RHS Graph: Plotting 
nb_ddsg = Matrix{Int}(undef, 5, 3)
errors_ddsg = Matrix{Float64}(undef, 5, 3)
nb_sg = Vector{Int}()
errors_sg = Vector{Float64}()
# %%
c = 3
num_points = 1000
F_poly(X::AbstractVector{Float64}; c)::Vector{Float64} = [sum(sin, X)^c]
F_poly(X::AbstractMatrix{Float64}; c)::Matrix{Float64} = reduce(hcat, [F_poly(col;c=c) for col in eachcol(X)])
F(X::Vector{Float64}) = F_poly(X, c=c)
F(X::Matrix{Float64}) = F_poly(X, c=c)
dim = 20
dof = 1
X_sample = rand(dim, num_points)
Y_exact = F(X_sample)
# %%
for l in 1:5
    for k_max in 1:3
        ddsg = Dynare.DDSG(dim, dof, l, l, k_max)
        Dynare.DDSG_init!(ddsg)
        Dynare.DDSG_build!(
            ddsg,
            F,
            ddsg.centroid
        )
        # %%
        Y_ddsg = Dynare.DDSG_evaluate(ddsg, X_sample)
        errors_ddsg[l,k_max] = mean(@. abs((Y_ddsg - Y_exact)/Y_exact))
        nb_ddsg[l,k_max] =  sum(ddsg.grid_points)
    end
    # %%
    sg = Tasmanian.TasmanianSG(dim, dof, l)
    Tasmanian.makeLocalPolynomialGrid!(sg)
    domain = zeros(dim,2)
    domain[:,2] .= 1.
    Tasmanian.setDomainTransform!(sg, domain)  # Match grid to domain dimensions
    # Refinement based on surplus coefficients
    X = Tasmanian.getNeededPoints(sg)
    N = size(X,2)
    Y_val = F(X)
    Tasmanian.loadNeededPoints!(sg, Y_val)
    Y_sg = Tasmanian.evaluateBatch(sg,X_sample)
    # %%
    push!(errors_sg, mean(@. abs((Y_sg - Y_exact)/Y_exact)))
    push!(nb_sg, Tasmanian.getNumPoints(sg))
end
# %%
rhs = plot(nb_sg, errors_sg, color = :black, label=L"SG",
     markershape=:star6, markersize=5, linewidth = 4)
plot!(nb_ddsg[:,1], errors_ddsg[:,1], color = :blue, label=L"DDSG~(\mathcal{K} = 1)",
      markershape=:circle, markersize=5, linewidth = 4)
plot!(nb_ddsg[:,2], errors_ddsg[:,2], color = :red, label=L"DDSG~(\mathcal{K} = 2)",
      markershape=:xcross, markersize=5, linewidth = 4, linestyle = :dash)
plot!(nb_ddsg[:,3], errors_ddsg[:,3], color = :green, label=L"DDSG~(\mathcal{K} = 3)",
      markershape=:rect, markersize=5, linewidth = 4, linestyle = :dot)
plot!(xscale=:log10, yscale=:log10, minorgrid=true)
xlabel!("Number of Grid Points")
# %%
plot(lhs, rhs, layout=(1, 2), link=:y)
# %%
# 1) Graph of the number of needed iterations to reach a tol_ti in function of the
# dimension
# 2) Graph of average iteration time for a single thread w.r.t # Dimensions
# 3) N=2, single thread, compute the number of needed iterations: compare first-order linear guess vs the full-zero policy guess
# DDSG: show the time spent is higher, but allows to reach higher dimensions than ASG or SG