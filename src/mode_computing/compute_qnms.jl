"""
    computeQNMs(data::QuadFreqData;
                plr_polish::Bool = true, 
                plr_epsilon::BigFloat = BigFloat(1.0e-20)
                )

Compute the quasinormal modes from the previouslly created `QuadFreqData` space-time data.

# Input
- `data::QuadFreqData`: The previously defined space-time data. 
- `plr_polish::Bool`: Tell PolynomialRoots to divide the original polynomial by each root found and polish the results using the full polynomial.
- `plr_epsilon::BigFloat`: The stopping criterion described in Skowron & Gould paper. This is not the precision with which the roots will be calculated.

# Output
An object of type `Array{Complex{BigFloat},1}` containing the computed modes.
"""
function computeQNMs(data::QuadFreqData;
                     plr_polish::Bool = true, 
                     plr_epsilon::BigFloat = BigFloat("1.0e-20")
                     ) 
    
                     # Compute the AIM "quantization condition"
    δ = computeDelta(data)

    # Solve the quantization condition to obtain Quasinormal modes and sort by their imaginary value
    qnms = PolynomialRoots.roots(coeffs(δ), polish = plr_polish, epsilon = plr_epsilon)
    
    # Empty qnm array warning
    if isempty(qnms)
        println("Warning: The computed mode array is empty. This means that no roots of the polynomial equation in ω were found.")
    end
    
    sort!(qnms, by = x -> imag(x))

    return qnms
end 

"""
    computeQNMs(data::GenFreqData,
                guess::Complex{BigFloat};

                nls_xtol::BigFloat = BigFloat("1.0e-15"),
                nls_ftol::BigFloat = BigFloat("1.0e-15"),
                nls_iterations::Int64 = 1000
                )
Compute a single quasinormal mode from the previously created `GenFreqData` space-time data.

# Input
- `data::GenFreqData`: The previously defined numeric space-time data.
- `guess::Complex{BigFloat}`: The initial guess used for computing the mode.
- `nls_xtol::BigFloat`: Norm difference in x between two successive iterates under which convergence is declared.
- `nls_ftol::BigFloat`: Infinite norm of residuals under which convergence is declared.
- `nls_iterations::Int64`: Maximum number of iterations performed by NLsolve.

# Output
An object of type `SolverResults` returned by `nlsolve`. See [NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl for further details)
"""
function computeQNMs(data::GenFreqData,
                     guess::Complex{BigFloat};

                     nls_xtol::BigFloat = BigFloat("1.0e-15"),
                     nls_ftol::BigFloat = BigFloat("1.0e-15"),
                     nls_iterations::Int64 = 1000
                     )

    # This function is passed to NLsolve to find the roots of δ
    function f!(F::Array{BigFloat,1}, x::Array{BigFloat,1})
        y = computeDelta(data, Complex(x[1], x[2]))
        F[1] = real(y)
        F[2] = imag(y)
    end
    
    # compute the roots of δ using NLsolve
    return nlsolve(f!, [real(guess), imag(guess)], ftol=nls_ftol, xtol=nls_xtol, iterations=nls_iterations)
end
