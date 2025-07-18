import Pkg; Pkg.activate(".")
using Dynare
using Tasmanian
using Statistics
using MAT
using VectorizedStatistics
# using Revise
# %%
function get_EE_errors(errorMat; burnin=0.0909, q=0.999)
    start=round(Int, burnin*size(errorMat,2))
    EEmat = abs.(errorMat[1:end-1, start:end, :])
    EEmat_copy = copy(EEmat)
    q_err = vquantile!(EEmat_copy, q, dims=(2))
    return (
        log10(mean(EEmat)),
        log10(mean(q_err))
    )
end
# %%
# Table 1
context = dynare("irbc_small_gl", "stoponerror");
# %%
l = 3
sg_options = SGOptions(scaleCorrExclude=["lambda"], gridDepth=l, tol_ti=1e-7, ftol=1e-8, maxiter=2000, polUpdateWeight=1.)
SG_grid, sgws = SGapproximation(sg_options; context=context);
# %%
errorMat = simulation_approximation_error!(context=context,grid=SG_grid,sgws=sgws)
# %%
avg_error, q_error = get_EE_errors(errorMat)
# %%
# Store the integration nodes and weights for use in the MATLAB file so_ee_err.m
int_nodes = sgws.monomial.nodes
int_nodes = reduce(hcat, int_nodes)
int_weights = sgws.monomial.weights
file = matopen("quadrature.mat","w")
write(file, "int_nodes", int_nodes)
write(file, "int_weights", int_weights)
close(file)
