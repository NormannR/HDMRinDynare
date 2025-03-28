# %%
import Pkg; Pkg.activate(".")
using Dynare
using Plots, LaTeXStrings
# %%
# RBC Model
context = dynare("rbc.mod", "stoponerror");
# %%
# --- Parameters ---
const α = 1 / 3
const β = 0.99
const ρ = 0.95
const σ = 0.01
const ss_l = 1 / 3
const ss_k = ss_l * (β * α)^(1 / (1 - α))

const GRID_SIZE = 100
const Z_SLICE_IDX = 50

# --- Policy Functions ---
c_pol(K, Z) = (1 - α * β) * exp(Z) * K^α * ss_l^(1 - α)
k_pol(K, Z) = α * β * exp(Z) * K^α * ss_l^(1 - α)
y_pol(K, Z) = c_pol(K, Z) + k_pol(K, Z)

# --- Initial Policy Guess ---
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

# --- Evaluate Policy on Grid ---
function evaluate_policy(pol, idx, Z_vals, K_vals)
    return [pol([K, Z])[idx] for Z in Z_vals, K in K_vals]
end

# --- Plotting Function ---
function plot_policies(K_vals, theoretical, sg, ddsg1, ddsg2, label)
    p = plot(K_vals, theoretical[Z_SLICE_IDX, :], label = "Theoretical", title = "$label Policy")
    plot!(K_vals, sg[Z_SLICE_IDX, :],     label = "SG",         linestyle = :dash)
    plot!(K_vals, ddsg1[Z_SLICE_IDX, :], label = "DDSG (K = 1)", linestyle = :dash)
    plot!(K_vals, ddsg2[Z_SLICE_IDX, :], label = "DDSG (K = 2)", linestyle = :dash)
    return p
end
# %%
context = dynare("rbc.mod", "stoponerror")

sg_options = SGOptions(tol_ti = 1e-6, gridDepth = 3, maxRef = 1, polUpdateWeight = 1.0,
                        maxIterEarlyStopping = 10, initialPolGuess = initialPolGuess)
ddsg_options_1 = DDSGOptions(tol_ti = 1e-6, gridDepth = 3, polUpdateWeight = 0.1, initialPolGuess = initialPolGuess)
ddsg_options_2 = DDSGOptions(tol_ti = 1e-6, gridDepth = 3, polUpdateWeight = 1.0, k_max = 2, initialPolGuess = initialPolGuess)

(SG_grid, _) = SGapproximation(sg_options)
(DDSG_grid_1, _) = DDSGapproximation(ddsg_options_1)
(DDSG_grid_2, _) = DDSGapproximation(ddsg_options_2)

Z_vals = range(-2 * σ / sqrt(1 - ρ^2), 2 * σ / sqrt(1 - ρ^2), length = GRID_SIZE)
K_vals = range(ss_k * 0.1, ss_k * 1.9, length = GRID_SIZE)

c_theoretical = [c_pol(K, Z) for Z in Z_vals, K in K_vals]
k_theoretical = [k_pol(K, Z) for Z in Z_vals, K in K_vals]

pol_sg      = x -> Dynare.interpolate(SG_grid, x)
pol_ddsg_1  = x -> Dynare.interpolate(DDSG_grid_1, x)
pol_ddsg_2  = x -> Dynare.interpolate(DDSG_grid_2, x)

i_c = context.symboltable["c"].orderintype
i_k = context.symboltable["k"].orderintype
i_c = findall(i_c .== context.models[1].i_dyn)[1]
i_k = findall(i_k .== context.models[1].i_dyn)[1]

c_sg      = evaluate_policy(pol_sg,     i_c, Z_vals, K_vals)
c_ddsg_1  = evaluate_policy(pol_ddsg_1, i_c, Z_vals, K_vals)
c_ddsg_2  = evaluate_policy(pol_ddsg_2, i_c, Z_vals, K_vals)
k_sg      = evaluate_policy(pol_sg,     i_k, Z_vals, K_vals)
k_ddsg_1  = evaluate_policy(pol_ddsg_1, i_k, Z_vals, K_vals)
k_ddsg_2  = evaluate_policy(pol_ddsg_2, i_k, Z_vals, K_vals)

p1 = plot_policies(K_vals, c_theoretical, c_sg, c_ddsg_1, c_ddsg_2, "Consumption")
ylabel!(L"C_t \mid Z_t = Z")

p2 = plot_policies(K_vals, k_theoretical, k_sg, k_ddsg_1, k_ddsg_2, "Capital")
ylabel!(L"K_t \mid Z_t = Z")

plot(p1, p2, layout = (1, 2))
xlabel!(L"K_{t-1}")
savefig("rbc_comparison.pdf")