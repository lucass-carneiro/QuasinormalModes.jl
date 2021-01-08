"""
    dnx(n::T, f::Basic, v::Basic) where {T <: Unsigned}

Computes the n-th derivative of the the function f with respect to the variable v.
This function re implements SymEngine's own `diff` function using an early quitting strategy.

# Input
- `n::T`: The order of the derivative.
- `f::Basic`: The expression to derivate.
- `v::Basic`: The variable to derivate with.

# Output
A `SymEngine.Basic` object with the derived expression.
"""
function dnx(n::T, f::Basic, v::Basic) where {T <: Unsigned}
    b0 = Basic(0)
    zt = zero(T)
    ot = one(T)
    
    f == b0 && return b0
    n == zt && return f
    n == ot && return SymEngine.diff(f, v)
    n > 1 && return dnx(n - ot, SymEngine.diff(f, v), v)
end

"""
    dnx(p::AnalyticAIMProblem{N,T}, n::N, f::Function, v::Function) where {N <: Unsigned, T <: Number}

Computes the n-th derivative of the AIM expressions with respect to ODE's variable.
This function is only a thin wrapper around SymEngine's own `diff` function.
It works as a barrier function that produces a type stable `Basic` result.

# Input
- `p::AnalyticAIMProblem`: The problem data with the expressions to derivate.
- `n::Unsigned`: The order of the derivative.
- `f::Function`: The actual expression to derivate. Either λ0 or S0.
- `v::Function`: The variable to derivate with. Either get_ODEvar or get_ODEeigen.

# Output
A `SymEngine.Basic` object with the derived expression.
"""
function dnx(p::AnalyticAIMProblem{N,T}, n::N, f::Function, v::Function) where {N <: Unsigned, T <: Number}
    return convert(SymEngine.Basic, dnx(n, f(p), v(p)))
end

"""
    computePolynomialFactors(p::QuadraticEigenvalueProblem{N,T}, n::N, f::Function) where {N <: Unsigned, T <: Number}

Create a second order `Polynomial` object in the variable `ω` by computing derivatives of `λ0` or `S0`.

# Input
- `p::QuadraticEigenvalueProblem`: The problem data with the expressions to derivate.
- `n::Unsigned`: The order of the derivative.
- `f::Function`: the function to extract the polynomial from. Either λ0 or S0.

# Output
An object of type `Polynomial{T}` containing the polynomial resulting from the derivation of the expression.
"""
function computePolynomialFactors(p::QuadraticEigenvalueProblem{N,T}, n::N, f::Function) where {N <: Unsigned, T <: Number}
    SYB = SymEngine.Basic

    # Compute the derivative at x0 and expand it to obtain a polynomial in ω
    step1::SYB = convert(SYB, subs(dnx(p, n, f, get_ODEvar), get_ODEvar(p) => get_x0(p)))
    step2::SYB = convert(SYB, expand(step1))
    
    # Extract the factors of the polynomial in ω
    ω = get_ODEeigen(p)

    p0::T = convert(T, SymEngine.N(coeff(step2, ω, Basic(0))))
    p1::T = convert(T, SymEngine.N(coeff(step2, ω, Basic(1))))
    p2::T = convert(T, SymEngine.N(coeff(step2, ω, Basic(2))))

    poly = Polynomial{T}([p0, p1, p2])

    return poly
end

"""
    createPoly(p::QuadraticEigenvalueProblem{N,T}, n::N, f::Function)

Compute the n-th coefficient of the Taylor expansion around x0 for the functions `λ0` or `S0`

# Input
- `p::QuadraticEigenvalueProblem`: The problem data with the expressions to derivate.
- `n::Unsigned`: The order of the derivative.
- `f::Function`: the function to extract the polynomial from. Either λ0 or S0.

# Output
An object of type `Polynomial{T}` containing the polynomial Taylor coefficient.
"""
function createPoly(p::QuadraticEigenvalueProblem{N,T}, n::N, f::Function) where {N <: Unsigned, T <: Number}
    computePolynomialFactors(p, n, f)/factorial(n)
end

function createPoly(p::QuadraticEigenvalueProblem{N,T}, n::N, f::Function) where {N <: Unsigned, T <: BigFloat}
    computePolynomialFactors(p, n, f)/factorial(big(n))
end

function createPoly(p::QuadraticEigenvalueProblem{N,T}, n::N, f::Function) where {N <: Unsigned, T <: Complex{BigFloat}}
    computePolynomialFactors(p, n, f)/factorial(big(n))
end
