# %%
import Pkg; Pkg.activate(".")
using Dynare
using Tasmanian
using Plots, LaTeXStrings
using Statistics
# %%
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
savefig("analytical_ddsg.pdf")