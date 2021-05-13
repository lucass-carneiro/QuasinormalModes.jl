"""
    eigenvaluesInGrid(
        m::AIMSteppingMethod,
        p::NumericAIMProblem{N,T},
        c::AIMCache{N,T},
        grid::Tuple{T,T,Int64,Int64};
        xtol::Real = 1.0e-10,
        ftol::Real = 1.0e-10,
        iterations::Int = 1000
        ) where {N <: Unsigned, T <: Complex}

Attempts to find eigenvalues using a grid of complex plane data points as initial guesses passed to nlsolve.

# Input
- `m::AIMSteppingMethod`: The stepping method to use.
- `p::NumericAIMProblem{N,T}`: The previously defined problem data.
- `c::AIMCache{N,T}`: The cache constructed from p.
- `grid::Tuple{T,T,Int64,Int64}`: A tuple consisting of (start point, end point, num. of real pts., num. of imag. pots.).
- `xtol::Real`: Norm difference in x between two successive iterates under which convergence is declared.
- `ftol::Real`: Infinite norm of residuals under which convergence is declared.
- `iterations::Int`: Maximum number of iterations performed by NLsolve.

# Output
An object of type `Array{T,1}` containing the modes found within the grid.
"""
function eigenvaluesInGrid(
    m::AIMSteppingMethod,
    p::NumericAIMProblem{N,T},
    c::AIMCache{N,T},
    grid::Tuple{T,T,Int64,Int64};
    xtol::Real = 1.0e-10,
    ftol::Real = 1.0e-10,
    iterations::Int = 1000
    ) where {N <: Unsigned, T <: Complex}

    if grid[3] < 0 || grid[4] < 0
        error("The number of points in the grid must be positive")
    end

    # Extract grid quantities
    re_start = real(grid[1])
    re_end = real(grid[2])
    re_size = grid[3]
    re_step = (re_end-re_start)/(re_size-1)
    re_range = re_start:re_step:re_end

    im_start = imag(grid[1])
    im_end = imag(grid[2])
    im_size = grid[4]
    im_step = (im_end-im_start)/(im_size-1)
    im_range = im_start:im_step:im_end

    # Check that the grid is  valid
    if (re_end <= re_start || re_step < 0 ) || (im_end <= im_start || im_step < 0 )
        error("The grid specified must be tuple in the form (start, end, real points, imag. points) where end > start.")
    end
    
    # Build the output data array with it's expected size
    eigenvalues::Array{T,1} = []

    for realPart in re_range
        for imagPart in im_range
            # Compute the solution, if any
            sol = computeEigenvalues(
                m,
                p,
                c,
                T(realPart, imagPart),
                nls_xtol = xtol,
                nls_ftol = ftol,
                nls_iterations = iterations
                )

            # If nlsolve converge the result is stored. If not take the appropriate action
            if sol.f_converged || sol.x_converged
                push!(eigenvalues, T(sol.zero[1], sol.zero[2]))
            end
        end
    end

    return eigenvalues
end

"""
    eigenvaluesInGrid(
        m::AIMSteppingMethod,
        p::NumericAIMProblem{N,T},
        c::AIMCache{N,T},
        grid::Tuple{T,T};
        roots_atol::Real = 1.0e-10,
        roots_rtol::Real = 1.0e-10,
        roots_xatol::Real = 1.0e-10,
        roots_xrtol::Real = 1.0e-10,
        roots_maxevals::Int = 100,
        roots_maxfnevals::Int = 100,
        ) where {N <: Unsigned, T <: Real}

Attempts to find eigenvalues using a range of real data points as a search region to find_zeros.
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
- `roots_maxfnevals::Int`:Number of function evaluations performed.

# Output
An object of type `Array{T,1}` containing the eigenvalues found within the grid.
"""
function eigenvaluesInGrid(
    m::AIMSteppingMethod,
    p::NumericAIMProblem{N,T},
    c::AIMCache{N,T},
    grid::Tuple{T,T};
    roots_atol::Real = 1.0e-10,
    roots_rtol::Real = 1.0e-10,
    roots_xatol::Real = 1.0e-10,
    roots_xrtol::Real = 1.0e-10,
    roots_maxevals::Int = 100,
    roots_maxfnevals::Int = 100,
    ) where {N <: Unsigned, T <: Real}

    if grid[1] > grid[2]
        error("The interval for searching eigenvalues must be ordered.")
    end

    try
        find_zeros(
            x -> computeDelta!(m, p, c, x),
            grid[1],
            grid[2],
            atol = roots_atol,
            rtol = roots_rtol,
            xatol = roots_xatol,
            xrtol = roots_xrtol,
            maxevals = roots_maxevals,
            maxfnevals = roots_maxfnevals
            )
    catch
        println("find_zeros was unable to converge to an eigenvalue")
    end
end
