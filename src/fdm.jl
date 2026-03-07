using FiniteDifferenceMatrices
using SparseArrays

function fdm(; Δx = 0.1, X = -10:Δx:10)
    H = -1/2 * fdmatrix(length(X), n=2, m=2, d=:c, h=Δx) + spdiagm([1/2*x^2 for x in X])
    ψ = [π^(-1//4)*exp(-1/2*x^2) for x in X]
    E = (ψ' * H * ψ) / (ψ' * ψ)
    return E
end