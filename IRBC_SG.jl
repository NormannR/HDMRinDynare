import Pkg; Pkg.activate(".")
using Dynare
using Tasmanian
using Statistics
# %%
function get_EE_errors(errorMat; burnin=0.0909, q=0.999)
   start=round(Int, burnin*size(errorMat,2))
   return (
      log10(mean(abs.(errorMat[end-1,start:end,:]))),
      log10(mean([quantile(abs.(errorMat[end-1,start:end,r]),q) for r in axes(errorMat)[3]]))
   )
end
# %%
# Table 1
context = dynare("irbc_small", "-DN=2","stoponerror");
nb_points = Vector{Int}()
avg_errors = Vector{Float64}()
q_errors = Vector{Float64}()
for l in [3,5,7]
   println(l)
   sg_options = SGOptions(scaleCorrExclude=["lambda"], gridDepth=l, tol_ti=1e-7, ftol=1e-8, maxiter=400, polUpdateWeight=1.)
   SG_grid, sgws = SGapproximation(sg_options; context=context)
   errorMat = simulation_approximation_error!(context=context,grid=SG_grid,sgws=sgws)
   avg_error, q_error = get_EE_errors(errorMat)
   push!(nb_points, getNumPoints(SG_grid))
   push!(avg_errors, avg_error)
   push!(q_errors, q_error)
end
# %%
nb_points
# %%
avg_errors
# %%
q_errors
# %%
N = 4
println(2*N)
context = dynare("irbc_small", "-DN=$N", "stoponerror");
sg_options = SGOptions(scaleCorrExclude=["lambda"], gridDepth=3, tol_ti=1e-7, ftol=1e-8, maxiter=100, polUpdateWeight=1.)
@profview SG_grid, sgws = SGapproximation(sg_options; context=context)
# %%
# Table 2
nb_points = Vector{Int}()
avg_errors = Vector{Float64}()
q_errors = Vector{Float64}()
avg_time = Vector{Float64}()
for N in [2,4,8]
   println(2*N)
   context = dynare("irbc_small", "-DN=$N", "stoponerror");
   sg_options = SGOptions(scaleCorrExclude=["lambda"], gridDepth=3, tol_ti=1e-7, ftol=1e-8, maxiter=100, polUpdateWeight=1.)
   SG_grid, sgws = SGapproximation(sg_options; context=context)
   errorMat = simulation_approximation_error!(context=context,grid=SG_grid,sgws=sgws)
   avg_error, q_error = get_EE_errors(errorMat)
   push!(nb_points, getNumPoints(SG_grid))
   push!(avg_errors, avg_error)
   push!(q_errors, q_error)
   push!(avg_time, context.results.model_results[1].sparsegrids.average_iteration_time/1e3)
end
# %%
nb_points
# %%
avg_errors
# %%
q_errors
# %%
avg_time
# %%
# Table 3
iteration_number_naive =  Vector{Int}()
iteration_number_fo =  Vector{Int}()
# Naive initial guess
for N in [2,4,8]
   println(2*N)
   context = dynare("irbc_small", "-DN=$N", "stoponerror");
   initialPolGuess = UserPolicyGuess(
      x -> ones(N+1),
      vcat(["k_$i" for i=1:N], ["a_$i" for i=1:N]),
      vcat(["lambda"],["k_$i" for i=1:N])
   )
   sg_options = SGOptions(scaleCorrExclude=["lambda"], gridDepth=3, tol_ti=1e-4, ftol=1e-8, initialPolGuess=initialPolGuess, polUpdateWeight=1.)
   SG_grid, sgws = SGapproximation(sg_options; context=context)
   push!(iteration_number_naive, context.results.model_results[1].sparsegrids.iteration_number)
end
# %%
# First-order Dynare-provided initial guess
for N in [2,4,8]
   println(2*N)
   context = dynare("irbc_small", "-DN=$N", "stoponerror");
   sg_options = SGOptions(scaleCorrExclude=["lambda"], gridDepth=3, tol_ti=1e-4, ftol=1e-8, polUpdateWeight=1.)
   SG_grid, sgws = SGapproximation(sg_options; context=context)
   push!(iteration_number_fo, context.results.model_results[1].sparsegrids.iteration_number)
end
# %%
iteration_number_naive
# %%
iteration_number_fo
# %%
