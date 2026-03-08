# https://juliafewbody.github.io/FiniteDifferenceMatrices.jl/stable/examples/#Discrete-Approximation-of-Hamiltonian

using FiniteDifferenceMatrices
using SparseArrays

function fdm(;
    n::Int = 2^10,
    x_min::Float64 = -50.0,
    x_max::Float64 = 50.0,
)
    Δx = (x_max - x_min) / (n - 1)
    X = range(x_min, x_max; length=n)
    H = -1/2 * fdmatrix(length(X), n=2, m=2, d=:c, h=Δx) + spdiagm([1/2*x^2 for x in X])
    ψ = [π^(-1//4) * exp(-1/2*x^2) for x in X]
    E = (ψ' * H * ψ) / (ψ' * ψ)
    return E
end