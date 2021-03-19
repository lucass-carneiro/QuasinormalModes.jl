"""
Cache of coefficient arrays for the AIM.
To each AIM problem corresponds a cache. As long as the problem doesn't change, the cache can be reused.

# Members
- `icda::Array{T,1}`: Hold the initial c data, i.e., c^i_0.
- `ccda::Array{T,1}`: Hold the coefficients for the current aim step, c^i_n.
- `pcda::Array{T,1}`: Hold the coefficients for the previous aim step, c^i_{n-1}.
- `bcda::Array{T,1}`: The work buffer used to actually compute the c coefficients in parallel.
- `idda::Array{T,1}`: Hold the initial d data, i.e., c^i_0.
- `cdda::Array{T,1}`: Hold the coefficients for the current aim step, d^i_n.
- `pdda::Array{T,1}`: Hold the coefficients for the previous aim step, d^i_{n-1}.
- `bdda::Array{T,1}`: The work buffer used to actually compute the d coefficients in parallel.
- `size::N`: The size of the arrays in the cache.
"""
struct AIMCache{N <: Unsigned, T <: Any}
    icda::Array{T,1}
    ccda::Array{T,1}
    pcda::Array{T,1}
    bcda::Array{T,1}

    idda::Array{T,1}
    cdda::Array{T,1}
    pdda::Array{T,1}
    bdda::Array{T,1}

    size::N
end

"""
    AIMCache(p::QuadraticEigenvalueProblem{N,T}) where {N <: Unsigned, T <: Number}

Create an AIMCache object suitable for Quadratic Eigenvalue Problems.

# Input
- `p::QuadraticEigenvalueProblem`: The problem data.

# Output
An `AIMCache{N,Polynomial{T}}` object.
"""
function AIMCache(p::QuadraticEigenvalueProblem{N,T}) where {N <: Unsigned, T <: Number}
    P = Polynomial{T}
    size = get_niter(p) + one(N)
    
    # The initial data arrays are initialized using the Taylor coefficient formula f^(k)(x0)/k!
    # f^(k)(x0) is a polynomial constructed by createPoly
    icda = [createPoly(p, i, λ0) for i in zero(N):get_niter(p)]
    ccda = zeros(P, size)
    pcda = zeros(P, size)
    bcda = zeros(P, size)

    idda = [createPoly(p, i, S0) for i in zero(N):get_niter(p)]
    cdda = zeros(P, size)
    pdda = zeros(P, size)
    bdda = zeros(P, size)

    AIMCache{N,P}(icda, ccda, pcda, bcda, idda, cdda, pdda, bdda, size)
end

"""
    AIMCache(p::NumericAIMProblem{N,T}) where {N <: Unsigned, T <: Number}

Create an AIMCache object suitable for Numeric Eigenvalue Problems.

# Input
- `p::NumericAIMProblem`: The problem data.

# Output
An `AIMCache{N,T}` object.
"""
function AIMCache(p::NumericAIMProblem{N,T}) where {N <: Unsigned, T <: Number}
    size = get_niter(p) + one(N)

    icda = zeros(T, size)
    ccda = zeros(T, size)
    pcda = zeros(T, size)
    bcda = zeros(T, size)

    idda = zeros(T, size)
    cdda = zeros(T, size)
    pdda = zeros(T, size)
    bdda = zeros(T, size)

    AIMCache{N,T}(icda, ccda, pcda, bcda, idda, cdda, pdda, bdda, size)
end

"""
    recomputeInitials(p::NumericAIMProblem{N,T}, c::AIMCache{N,T}, ω::T) where {N <: Unsigned, T <: Number}

Reevaluate (in-place) the initial data arrays. The initial data array elements are the Taylor expansion
coefficients of λ0 and S0 in the ODE variable x around x0 of order get_niter(p) at a point ω.

# Input
- `p::NumericAIMProblem`: The problem data.
- `c::AIMCache`: The problem data associated cache.
- `ω::T`: The value of the eigenvalue to evaluate the arrays in.

# Output
    nothing
"""
function recomputeInitials!(p::NumericAIMProblem{N,T}, c::AIMCache{N,T}, ω::T) where {N <: Unsigned, T <: Number}
    t = get_x0(p) + Taylor1(T, convert(Int64, get_niter(p)))
    c1 = λ0(p)(t,ω).coeffs
    c2 = S0(p)(t,ω).coeffs

    copy!(c.icda, c1)
    copy!(c.idda, c2)

    return nothing
end
