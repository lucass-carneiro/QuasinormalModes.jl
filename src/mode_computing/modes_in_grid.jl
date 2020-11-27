"""
    function modesInGrid(data::GenFreqData,
                         grid::Tuple{Complex{Float64},Complex{Float64},Int64,Int64};
                         xtol::BigFloat = BigFloat("1.0e-15"),
                         ftol::BigFloat = BigFloat("1.0e-15"),
                         iterations::Int64 = 1000
                         )::Array{Complex{BigFloat},1}

Attempts to find QNMs using a grid of complex plane data points as initial guesses to nlsolve.

# Input
- `data::GenFreqData`: The space-time data to use.
- `grid::Tuple{Complex{Float64},Complex{Float64},Int64,Int64}`: A tuple consisting of (start point, end point, num. of real pts., num. of imag. pots.).
- `xtol::BigFloat`: Norm difference in x between two successive iterates under which convergence is declared.
- `ftol::BigFloat`: Infinite norm of residuals under which convergence is declared.
- `iterations::Int64`: Maximum number of iterations performed by NLsolve.

# Output
An object of type `Array{Complex{BigFloat},1}` containing the modes found within the grid.
"""
function modesInGrid(data::GenFreqData,
                     grid::Tuple{Complex{BigFloat},Complex{BigFloat},Int64,Int64};
                     xtol::BigFloat = BigFloat("1.0e-15"),
                     ftol::BigFloat = BigFloat("1.0e-15"),
                     iterations::Int64 = 1000
                     )::Array{Complex{BigFloat},1}

    # Extract grid quantities
    re_start = real(grid[1])
    re_end = real(grid[2])
    re_size = UInt64(grid[3])
    re_step = (re_end-re_start)/(re_size-1)
    re_range = re_start:re_step:re_end

    im_start = imag(grid[1])
    im_end = imag(grid[2])
    im_size = UInt64(grid[4])
    im_step = (im_end-im_start)/(im_size-1)
    im_range = im_start:im_step:im_end

    # Check that the grid is  valid
    if (re_end <= re_start || re_step < 0 ) || (im_end <= im_start || im_step < 0 )
        error("The grid specified in modesInGrid must be tuple in the form (start, end, real points, imag. points) where end > start.")
        throw(InvalidGridException)
    end
    
    # Build the output data array with it's expected size
    modes = zeros(Complex{BigFloat}, re_size*im_size)

    Threads.@threads for i in eachindex(re_range)
        for j in eachindex(im_range)
            # Compute the solution, if any
            sol = computeQNMs(data,Complex(re_range[i], im_range[j]),nls_xtol=xtol,nls_ftol=ftol,nls_iterations=iterations)

            # If nlsolve converge the result is stored. If not take the appropriate action
            if sol.f_converged || sol.x_converged
                @inbounds modes[(i*re_size + j) - re_size] = Complex(sol.zero[1],sol.zero[2])
            else
                println("nlsolve was unable to converge to a mode using using ($re_start,$re_end) as a guess. Returning 0.0 + 0.0*im to the mode list.")
                flush(stdout)
                
                @inbounds modes[(i*re_size + j) - re_size] = Complex(BigFloat(0.0), BigFloat(0.0))
            end
        end
    end

    sort!(modes, by = x -> imag(x))

    return modes
end
