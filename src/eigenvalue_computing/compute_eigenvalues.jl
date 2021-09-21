"""
    computeEigenvalues(
        m::AIMSteppingMethod,
        p::QuadraticEigenvalueProblem{N,T},
        c::AIMCache{N,Polynomial{T}};
        plr_polish::Bool = true, 
        plr_epsilon::Real = convert(T, 1.0e-10)
        ) where {N <: Unsigned, T <: Number}

Compute the eigenvalues for the problem `p` with corresponding cache `c`.

# Input
- `m::AIMSteppingMethod`: The stepping method to use.
- `p::QuadraticEigenvalueProblem`: The previously defined problem data.
- `c::AIMCache`: The cache constructed from p.
- `plr_polish::Bool`: Tell PolynomialRoots to divide the original polynomial by each root found and polish the results using the full polynomial.
- `plr_epsilon::Real`: The stopping criterion described in Skowron & Gould paper. This is not the precision with which the roots will be calculated.

# Output
An object of type `Array{T,1}` containing the computed eigenvalues.
"""
function computeEigenvalues(
    m::AIMSteppingMethod,
    p::QuadraticEigenvalueProblem{N,T},
    c::AIMCache{N,Polynomial{T}};
    plr_polish::Bool=true, 
    plr_epsilon::Real=1.0e-10
    ) where {N <: Unsigned,T <: Number}
    
    # Compute the AIM "quantization condition"
    δ = computeDelta!(m, p, c)

    # Solve the quantization condition to obtain the eigenvalues
    eigenvalues = PolynomialRoots.roots(coeffs(δ), polish=plr_polish, epsilon=plr_epsilon)
    
    if isempty(eigenvalues)
        println("Warning: The computed mode array is empty. This means that no roots of the polynomial equation in ω were found.")
    end

    return eigenvalues
end 

"""
    computeEigenvalues(
        m::AIMSteppingMethod,
        p::NumericAIMProblem{N,T},
        c::AIMCache{N,T},
        guess::T;
        nls_xtol::Real = convert(T, 1.0e-10),
        nls_ftol::Real = convert(T, 1.0e-10),
        nls_iterations::Int = 1000
        ) where {N <: Unsigned, T <: Complex}

Compute a single eigenvalue for the problem `p` with corresponding cache `c`.

# Input
- `m::AIMSteppingMethod`: The stepping method to use.
- `p::NumericAIMProblem`: The previously defined problem data.
- `c::AIMCache`: The cache constructed from p.
- `nls_xtol::Real`: Norm difference in x between two successive iterates under which convergence is declared.
- `nls_ftol::Real`: Infinite norm of residuals under which convergence is declared.
- `nls_iterations::Int`: Maximum number of iterations performed by NLsolve.

# Output
An object of type `SolverResults` returned by `nlsolve`. See [NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl) for further details.
"""
function computeEigenvalues(
    m::AIMSteppingMethod,
    p::NumericAIMProblem{N,T},
    c::AIMCache{N,T},
    guess::T;
    nls_xtol::Real=1.0e-10,
    nls_ftol::Real=1.0e-10,
    nls_iterations::Int=1000
    ) where {N <: Unsigned,T <: Complex}

    # This function is passed to NLsolve to find the roots of δ
    function f!(F, x)
        y = computeDelta!(m, p, c, Complex(x[1], x[2]))
        F[1] = real(y)
        F[2] = imag(y)
    end
    
    # compute the roots of δ using NLsolve
    return nlsolve(f!, [real(guess), imag(guess)], ftol=nls_ftol, xtol=nls_xtol, iterations=nls_iterations)
end

"""
    computeEigenvalues(
        m::AIMSteppingMethod,
        p::NumericAIMProblem{N,T},
        c::AIMCache{N,T},
        guess::T;
        roots_atol::Real = 1.0e-10,
        roots_rtol::Real = 1.0e-10,
        roots_xatol::Real = 1.0e-10,
        roots_xrtol::Real = 1.0e-10,
        roots_maxevals::Int = 100
        ) where {N <: Unsigned, T <: Real}

Compute a single eigenvalue for the problem `p` with corresponding cache `c`.
For details on convergence settings see [Roots.jl](https://juliahub.com/docs/Roots/o0Xsi/1.0.7/reference/#Convergence).

# Input
- `m::AIMSteppingMethod`: The stepping method to use.
- `p::NumericAIMProblem{N,T}`: The previously defined problem data.
- `c::AIMCache{N,T}`: The cache constructed from p.
- `grid::Tuple{T,T}`: A tuple consisting of (start point, end point).
- `roots_atol::Real`: Absolute tolerance.
- `roots_rtol::Real`: Relative tolerance.
- `roots_xatol::Real`: Floating point comparison absolute tolerance.
- `roots_xrtol::Real`: Floating point comparison relative tolerance.
- `roots_maxevals::Int`: Number of algorithm iterations performed.

# Output
An object of type T containing the found eigenvalue.
"""
function computeEigenvalues(
    m::AIMSteppingMethod,
    p::NumericAIMProblem{N,T},
    c::AIMCache{N,T},
    guess::T;
    roots_atol::Real=1.0e-10,
    roots_rtol::Real=1.0e-10,
    roots_xatol::Real=1.0e-10,
    roots_xrtol::Real=1.0e-10,
    roots_maxevals::Int=100
    ) where {N <: Unsigned,T <: Real}

    try
        find_zero(
            x -> computeDelta!(m, p, c, x),
            guess,
            atol=roots_atol,
            rtol=roots_rtol,
            xatol=roots_xatol,
            xrtol=roots_xrtol,
            maxevals=roots_maxevals
            )
    catch
        println("find_zeros was unable to converge to an eigenvalue")
    end
end
