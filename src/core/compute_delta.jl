"""
    computeDelta!(p::QuadraticEigenvalueProblem{N,T}, c::AIMCache{N,Polynomial{T}}) where {N <: Unsigned, T <: Number}
    
Compute and return the AIM "quantization condition".

# Input
- `p::QuadraticEigenvalueProblem`: A quadratic frequency problem.
- `c::AIMCache`: An AIM cache object created from p.

# Output
An object of type `Polynomial{T}` whose roots are the the problem's eigenvalues.
"""
function computeDelta!(p::QuadraticEigenvalueProblem{N,T}, c::AIMCache{N,Polynomial{T}}) where {N <: Unsigned, T <: Number}
    
    if (get_niter(p) + one(N)) != c.size
        error("The provided cache cannot hold data for $(get_niter(p)) AIM iterations.
        Make sure that the passed AIMCache corresponds to the passed AIM problem")
    end

    # Initialize the current data arrays with the initial data for the first step
    copy!(c.ccda, c.icda)
    copy!(c.cdda, c.idda)

    # Perform the aim steps
    for i in one(N):get_niter(p)
        AIMStep!(p, c)
    end

    # Compute and return the AIM "quantization condition"
    return c.cdda[1]*c.pcda[1] - c.pdda[1]*c.ccda[1]
end

"""
    computeDelta!(p::NumericAIMProblem{N,T}, c::AIMCache{N,T}, ω::T) where {N <: Unsigned, T <: Number}

Compute and return the AIM "quantization condition".

# Input
- `p::QuadraticEigenvalueProblem`: A quadratic frequency problem.
- `c::AIMCache`: An AIM cache object created from p.
- `ω::T`: Point to evaluate the quantization condition.
    
# Output
An object of type `T` which represents the AIM quantization condition at point ω.
"""
function computeDelta!(p::NumericAIMProblem{N,T}, c::AIMCache{N,T}, ω::T) where {N <: Unsigned, T <: Number}
    
    if (get_niter(p) + one(N)) != c.size
        error("The provided cache cannot hold data for $(get_niter(p)) AIM iterations.
        Make sure that the passed AIMCache corresponds to the passed AIM problem")
    end

    # Fill iitial arrays with data corresponding to the point ω
    recomputeInitials!(p,c,ω)

    # Initialize the current data arryas with the initial data for the first step
    copy!(c.ccda, c.icda)
    copy!(c.cdda, c.idda)

    # Perform the aim steps
    for i in one(N):get_niter(p)
        AIMStep!(p, c)
    end

    # Compute and return the AIM "quantization condition" at point ω
    return c.cdda[1]*c.pcda[1] - c.pdda[1]*c.ccda[1]
end
