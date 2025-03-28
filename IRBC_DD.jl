# %%
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
# Table 4
nb_points = Vector{Int}()
avg_errors = Vector{Float64}()
q_errors = Vector{Float64}()
avg_time = Vector{Float64}()
# %%
for N in [2,4,8,25,50]
   context = dynare("irbc_small", "-DN=$N", "stoponerror");
   ddsg_opts = DDSGOptions(scaleCorrExclude=["lambda"], gridDepth=3, tol_ti=1e-7, ftol=1e-8, polUpdateWeight=1.);
   DDSG_grid, sgws = DDSGapproximation(ddsg_opts; context=context);
   errorMat = simulation_approximation_error!(context=context,grid=DDSG_grid,sgws=sgws, replications=100)
   avg_error, q_error = get_EE_errors(errorMat)
   push!(nb_points, Dynare.countPoints(DDSG_grid))
   push!(avg_errors, avg_error)
   push!(q_errors, q_error)
   push!(avg_time, context.results.model_results[1].sparsegrids.average_iteration_time/1e3)
end
# %%
# Table 5
context = dynare("irbc_small", "-DN=2","stoponerror");
# Count the number of iterations necessary to reach a time-iteration step lower
# than 1e-4 with a Dynare-provided first-order initial guess for policy
# functions
iteration_number_fo =  Vector{Int}()
nb_points = Vector{Int}()
for N in [2,4,8]
   println(N)
   context = dynare("irbc_small", "-DN=$N", "stoponerror");
   ddsg_opts = DDSGOptions(scaleCorrExclude=["lambda"], gridDepth=3, maxRef=0, tol_ti=1e-4, ftol=1e-8, polUpdateWeight=1.)
   DDSG_grid, sgws = DDSGapproximation(ddsg_opts; context=context)
   results = context.results.model_results[1].sparsegrids
   push!(iteration_number_fo, context.results.model_results[1].sparsegrids.iteration_number)
   push!(nb_points, Dynare.countPoints(DDSG_grid))
end
# %%
# Count the number of iterations necessary to reach a time-iteration step lower
# than 1e-4 with a native initial guess for policy functions
iteration_number_naive =  Vector{Int}()
for N in [2,4,8]
   println(N)
   context = dynare("irbc_small", "-DN=$N", "stoponerror");
   initialPolGuess = UserPolicyGuess(
      x -> ones(N+1),
      vcat(["k_$i" for i=1:N], ["a_$i" for i=1:N]),
      vcat(["lambda"],["k_$i" for i=1:N])
   )
   ddsg_opts = DDSGOptions(scaleCorrExclude=["lambda"], gridDepth=3, maxRef=0, tol_ti=1e-4, ftol=1e-8, polUpdateWeight=1., initialPolGuess=initialPolGuess)
   DDSG_grid, sgws = DDSGapproximation(ddsg_opts; context=context)
   push!(iteration_number_naive, context.results.model_results[1].sparsegrids.iteration_number)
end
# %%
iteration_number_fo
# %%
iteration_number_naive
# %%
nb_points