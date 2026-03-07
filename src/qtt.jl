using LinearAlgebra
import QuanticsTCI
import TensorCrossInterpolation

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

function qtt_value_at_index(qtt::QTTVector{T}, idx::Int) where {T}
    c = order(qtt)
    N = 1 << c
    (1 <= idx <= N) || throw(ArgumentError("idx must be in 1:$N"))

    idx0 = idx - 1
    v = ones(T, 1)
    @inbounds for k in 1:c
        G = qtt.cores[k]
        rL, _, rR = size(G)
        length(v) == rL || throw(ArgumentError("rank mismatch at site $k"))
        s = Int((idx0 >> (c - k)) & 0x1) + 1
        nxt = zeros(T, rR)
        for a in 1:rL, b in 1:rR
            nxt[b] += v[a] * G[a, s, b]
        end
        v = nxt
    end
    return v[1]
end

function qttvector_from_tci_1d(tci_obj, c::Int)
    tt = hasproperty(tci_obj, :tci) ? TensorCrossInterpolation.TensorTrain(getproperty(tci_obj, :tci)) :
                                      TensorCrossInterpolation.TensorTrain(tci_obj)
    # TCI may return TT with first core having rL>1; reverse ensures (1,2,r) structure
    tt = TensorCrossInterpolation.reverse(tt)
    tensors = collect(TensorCrossInterpolation.sitetensors(tt))
    length(tensors) == c || throw(ArgumentError("unexpected TCI length: got $(length(tensors)), expected $c"))

    cores = Vector{Array{Float64,3}}(undef, c)
    @inbounds for k in 1:c
        G = Array(tensors[k])
        ndims(G) == 3 || throw(ArgumentError("site tensor $k is not rank-3"))
        size(G, 2) == 2 || throw(ArgumentError("site tensor $k has local dimension $(size(G, 2)); expected 2"))
        cores[k] = Float64.(G)
    end
    return QTTVector(cores)
end

function orient_qtt_msb_first(qtt::QTTVector{Float64}, xvals::AbstractVector{<:Real}, f::Function)
    c = order(qtt)
    N = 1 << c
    length(xvals) == N || throw(ArgumentError("xvals length must be 2^c"))

    probe = unique(clamp.(Int[1, 2, 3, max(1, div(N, 2)), max(1, N - 2), N - 1, N], 1, N))
    # Reverse cores with permutedims so bond structure stays (1,2,r) for first core
    rev_cores = [permutedims(G, (3, 2, 1)) for G in reverse(qtt.cores)]
    rev_qtt = QTTVector(rev_cores)
    err_normal = 0.0
    err_reversed = 0.0
    @inbounds for i in probe
        target = f(xvals[i])
        err_normal += abs(qtt_value_at_index(qtt, i) - target)
        err_reversed += abs(qtt_value_at_index(rev_qtt, i) - target)
    end
    return err_reversed < err_normal ? rev_qtt : qtt
end

function build_qtt_vec_exp_tci(
    c::Int;
    interval::Tuple{Float64,Float64}=(-50.0, 50.0),
    tolerance::Float64=1e-10,
    maxbonddim::Int=64,
    maxiter::Int=100,
)
    a, b = interval
    c >= 1 || throw(ArgumentError("c must be >= 1"))
    b > a || throw(ArgumentError("interval must satisfy a < b"))
    N = 1 << c
    xvals = range(a, b; length=N)
    f(x) = exp(-0.5 * x * x)

    ci, _, _ = QuanticsTCI.quanticscrossinterpolate(
        Float64,
        f,
        xvals;
        tolerance=tolerance,
        maxbonddim=maxbonddim,
        maxiter=maxiter,
    )

    qtt = qttvector_from_tci_1d(ci, c)
    return orient_qtt_msb_first(qtt, xvals, f)
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

function qtt(;
    c::Int = 10,
    x_min::Float64 = -50.0,
    x_max::Float64 = 50.0,
    tci_tolerance::Float64 = 1e-10,
    tci_maxbonddim::Int = 64,
    tci_maxiter::Int = 100,
)
    lower, main, upper, _ = kinetic_fd_coeffs(c, (x_min, x_max))
    a1 = build_tridiag_toep(c, lower, main, upper)
    a2 = build_diagonal_qtt_from_vec_qtt(qtt_quadratic_1d(c; interval=(x_min, x_max)))
    a = a1 + a2

    phi_qtt = build_qtt_vec_exp_tci(
        c;
        interval=(x_min, x_max),
        tolerance=tci_tolerance,
        maxbonddim=tci_maxbonddim,
        maxiter=tci_maxiter,
    )

    num = qtt_expectation_mpo(phi_qtt, a)
    den = qtt_inner(phi_qtt, phi_qtt)
    return real(num / den)
end
