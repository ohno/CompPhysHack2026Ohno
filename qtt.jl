#!/usr/bin/env julia

using LinearAlgebra
using Printf

struct QTTVector{T}
    cores::Vector{Array{T,3}}  # (rL, 2, rR)
end

struct QTTMPO{T}
    cores::Vector{Array{T,4}}  # (rL, 2, 2, rR)
end

order(tt::Union{QTTVector,QTTMPO}) = length(tt.cores)

function kinetic_fd_coeffs(c::Int, interval::Tuple{Float64,Float64})
    a, b = interval
    c >= 1 || throw(ArgumentError("c must be >= 1"))
    b > a || throw(ArgumentError("interval must satisfy a < b"))
    N = 1 << c
    h = (b - a) / (N - 1)
    lower = -0.5 / (h * h)
    main = 1.0 / (h * h)
    upper = -0.5 / (h * h)
    return lower, main, upper, h
end

function build_tridiag_toep(c::Int, lower::Float64, main::Float64, upper::Float64)
    c >= 2 || throw(ArgumentError("c must be >= 2"))

    I2 = Matrix{Float64}(I, 2, 2)
    J = [0.0 1.0; 0.0 0.0]
    JT = transpose(J)

    f_core = zeros(Float64, 1, 2, 2, 3)
    f_core[1, :, :, 1] = I2
    f_core[1, :, :, 2] = JT
    f_core[1, :, :, 3] = J

    m_core = zeros(Float64, 3, 2, 2, 3)
    m_core[1, :, :, 1] = I2
    m_core[1, :, :, 2] = JT
    m_core[1, :, :, 3] = J
    m_core[2, :, :, 2] = J
    m_core[3, :, :, 3] = JT

    l_core = zeros(Float64, 3, 2, 2, 1)
    l_core[1, :, :, 1] = main * I2 + upper * J + lower * JT
    l_core[2, :, :, 1] = lower * J
    l_core[3, :, :, 1] = upper * JT

    cores = Vector{Array{Float64,4}}(undef, c)
    cores[1] = f_core
    @inbounds for k in 2:(c - 1)
        cores[k] = m_core
    end
    cores[c] = l_core
    return QTTMPO(cores)
end

function qtt_quadratic_1d(c::Int; interval::Tuple{Float64,Float64}=(-1.0, 1.0))
    a, b = interval
    c >= 1 || throw(ArgumentError("c must be >= 1"))
    b > a || throw(ArgumentError("interval must satisfy a < b"))

    if c == 1
        core = zeros(Float64, 1, 2, 1)
        core[1, 1, 1] = 0.5 * a * a
        core[1, 2, 1] = 0.5 * b * b
        return QTTVector([core])
    end

    N = 1 << c
    h = (b - a) / (N - 1)
    right_vec = (0.5 * a * a, a, 0.5)

    cores = Vector{Array{Float64,3}}(undef, c)
    pos = 1
    for i in (c - 1):-1:0
        t = h * (1 << i)
        if i == c - 1
            core = zeros(Float64, 1, 2, 3)
            core[1, 1, 1] = 1.0
            core[1, 2, 1] = 1.0
            core[1, 2, 2] = t
            core[1, 2, 3] = t * t
            cores[pos] = core
        elseif i == 0
            core = zeros(Float64, 3, 2, 1)
            core[1, 1, 1] = right_vec[1]
            core[2, 1, 1] = right_vec[2]
            core[3, 1, 1] = right_vec[3]
            core[1, 2, 1] = right_vec[1] + t * right_vec[2] + (t * t) * right_vec[3]
            core[2, 2, 1] = right_vec[2] + 2.0 * t * right_vec[3]
            core[3, 2, 1] = right_vec[3]
            cores[pos] = core
        else
            core = zeros(Float64, 3, 2, 3)
            @inbounds for d in 1:3
                core[d, 1, d] = 1.0
            end
            core[1, 2, 1] = 1.0
            core[1, 2, 2] = t
            core[1, 2, 3] = t * t
            core[2, 2, 2] = 1.0
            core[2, 2, 3] = 2.0 * t
            core[3, 2, 3] = 1.0
            cores[pos] = core
        end
        pos += 1
    end
    return QTTVector(cores)
end

function build_diagonal_qtt_from_vec_qtt(qtt::QTTVector{T}) where {T}
    c = order(qtt)
    cores = Vector{Array{T,4}}(undef, c)
    @inbounds for k in 1:c
        G = qtt.cores[k]
        rL, _, rR = size(G)
        W = zeros(T, rL, 2, 2, rR)
        for a in 1:rL, s in 1:2, b in 1:rR
            W[a, s, s, b] = G[a, s, b]
        end
        cores[k] = W
    end
    return QTTMPO(cores)
end

function build_qtt_vec_from_dense(values::AbstractVector{<:Real}, c::Int; reltol::Float64=0.0, maxrank::Int=typemax(Int))
    N = 1 << c
    length(values) == N || throw(ArgumentError("values length must be 2^c"))

    dims = ntuple(_ -> 2, c)
    # Convert to MSB-first tensor order, matching QTT core construction above.
    Tmsb = permutedims(reshape(Float64.(values), dims), reverse(1:c))

    cores = Vector{Array{Float64,3}}(undef, c)
    r_prev = 1
    X = reshape(Tmsb, r_prev * 2, :)

    @inbounds for k in 1:(c - 1)
        F = svd(X; full=false)
        U, S, Vt = F.U, F.S, F.Vt
        r = min(length(S), maxrank)
        if reltol > 0
            thresh = reltol * S[1]
            rr = 0
            for i in 1:r
                S[i] >= thresh || break
                rr = i
            end
            r = max(rr, 1)
        end

        Utrunc = @view U[:, 1:r]
        cores[k] = reshape(Array(Utrunc), r_prev, 2, r)

        Vtrunc = @view Vt[1:r, :]
        Xnext = Array(Vtrunc)
        for i in 1:r
            @views Xnext[i, :] .*= S[i]
        end
        X = Xnext
        r_prev = r
        if k < c - 1
            X = reshape(X, r_prev * 2, :)
        end
    end

    cores[c] = reshape(X, r_prev, 2, 1)
    return QTTVector(cores)
end

function qtt_inner(u::QTTVector{Tu}, v::QTTVector{Tv}) where {Tu,Tv}
    order(u) == order(v) || throw(ArgumentError("u and v must have the same order"))
    T = promote_type(Tu, Tv)
    E = ones(T, 1, 1)

    @inbounds for k in 1:order(u)
        Gu = u.cores[k]
        Gv = v.cores[k]
        ruL, _, ruR = size(Gu)
        rvL, _, rvR = size(Gv)
        size(E, 1) == ruL || throw(ArgumentError("rank mismatch in u at core $k"))
        size(E, 2) == rvL || throw(ArgumentError("rank mismatch in v at core $k"))

        Enew = zeros(T, ruR, rvR)
        for a in 1:ruL, b in 1:rvL
            e = E[a, b]
            e == 0 && continue
            for p in 1:ruR, q in 1:rvR
                Enew[p, q] += e * (conj(Gu[a, 1, p]) * Gv[b, 1, q] + conj(Gu[a, 2, p]) * Gv[b, 2, q])
            end
        end
        E = Enew
    end

    return E[1, 1]
end

function qtt_expectation_mpo(psi::QTTVector{Tpsi}, mpo::QTTMPO{Tmpo}) where {Tpsi,Tmpo}
    order(psi) == order(mpo) || throw(ArgumentError("psi and mpo must have the same order"))
    T = promote_type(Tpsi, Tmpo)
    E = ones(T, 1, 1, 1)

    @inbounds for k in 1:order(psi)
        G = psi.cores[k]
        W = mpo.cores[k]
        rpsiL, _, rpsiR = size(G)
        rWL, _, _, rWR = size(W)
        size(E, 1) == rpsiL || throw(ArgumentError("psi left-rank mismatch at core $k"))
        size(E, 2) == rWL || throw(ArgumentError("mpo left-rank mismatch at core $k"))
        size(E, 3) == rpsiL || throw(ArgumentError("bra/ket rank mismatch at core $k"))

        Enew = zeros(T, rpsiR, rWR, rpsiR)
        for x in 1:rpsiL, y in 1:rWL, z in 1:rpsiL
            e = E[x, y, z]
            e == 0 && continue
            for p in 1:rpsiR, q in 1:rpsiR, w in 1:rWR
                b1 = conj(G[x, 1, p])
                b2 = conj(G[x, 2, p])
                k1 = G[z, 1, q]
                k2 = G[z, 2, q]
                Enew[p, w, q] += e * (b1 * (W[y, 1, 1, w] * k1 + W[y, 1, 2, w] * k2) +
                                      b2 * (W[y, 2, 1, w] * k1 + W[y, 2, 2, w] * k2))
            end
        end
        E = Enew
    end

    return E[1, 1, 1]
end

function add_mpo(A::QTTMPO{Ta}, B::QTTMPO{Tb}) where {Ta,Tb}
    c = order(A)
    c == order(B) || throw(ArgumentError("MPO orders must match"))
    T = promote_type(Ta, Tb)
    C = Vector{Array{T,4}}(undef, c)

    A1, B1 = A.cores[1], B.cores[1]
    rA1, rB1 = size(A1, 4), size(B1, 4)
    C1 = zeros(T, 1, 2, 2, rA1 + rB1)
    C1[1, :, :, 1:rA1] = A1
    C1[1, :, :, (rA1 + 1):(rA1 + rB1)] = B1
    C[1] = C1

    @inbounds for k in 2:(c - 1)
        Ak, Bk = A.cores[k], B.cores[k]
        rAL, rAR = size(Ak, 1), size(Ak, 4)
        rBL, rBR = size(Bk, 1), size(Bk, 4)
        Ck = zeros(T, rAL + rBL, 2, 2, rAR + rBR)
        Ck[1:rAL, :, :, 1:rAR] = Ak
        Ck[(rAL + 1):(rAL + rBL), :, :, (rAR + 1):(rAR + rBR)] = Bk
        C[k] = Ck
    end

    Al, Bl = A.cores[c], B.cores[c]
    rAL = size(Al, 1)
    rBL = size(Bl, 1)
    Cl = zeros(T, rAL + rBL, 2, 2, 1)
    Cl[1:rAL, :, :, 1] = Al[:, :, :, 1]
    Cl[(rAL + 1):(rAL + rBL), :, :, 1] = Bl[:, :, :, 1]
    C[c] = Cl

    return QTTMPO(C)
end

Base.:+(A::QTTMPO, B::QTTMPO) = add_mpo(A, B)

function print_energy_table(; c_start::Int=8, c_end::Int=12, interval::Tuple{Float64,Float64}=(-50.0, 50.0))
    c_end >= c_start || throw(ArgumentError("c_end must be >= c_start"))

    for c in c_start:c_end
        t0 = time_ns()

        lower, main, upper, h = kinetic_fd_coeffs(c, interval)
        a1 = build_tridiag_toep(c, lower, main, upper)
        a2 = build_diagonal_qtt_from_vec_qtt(qtt_quadratic_1d(c; interval=interval))
        a = a1 + a2

        # Dense reference state -> exact QTT via TT-SVD.
        x = range(interval[1], interval[2], length=(1 << c))
        phi_dense = exp.(-0.5 .* (x .^ 2))
        phi_qtt = build_qtt_vec_from_dense(phi_dense, c)

        num = qtt_expectation_mpo(phi_qtt, a)
        den = qtt_inner(phi_qtt, phi_qtt)
        energy = real(num / den)
        elapsed = (time_ns() - t0) * 1e-9

        @printf("N=%d, h=%.6e, energy=%.12f t=%.5f secs\n", 1 << c, h, energy, elapsed)
    end
end

# if abspath(PROGRAM_FILE) == @__FILE__
#     print_energy_table()
# end

print_energy_table()
