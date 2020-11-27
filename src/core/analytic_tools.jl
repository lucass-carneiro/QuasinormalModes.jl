"""
    function dnx(data::QuadFreqData, n::UInt32, idx::UInt32)

Computes n-th derivative of the AIM expressions with respect to ODE's variable, usually called `x`.

# Input
- `data::QuadFreqData`: The spacetime data with the expressions to derivate.
- `n::UInt32`: The order of the derivative.
- `idx::UInt32`: The actual expression to derivate. If `idx = 0x00001` it computes a derivative of `λ0`. If `idx = 0x00002` it computes a derivative of `S0`.

# Output
A `SymEngine.Basic` object with the derived expression.
"""
function dnx(data::QuadFreqData, n::UInt32, idx::UInt32)
    ret::Basic = data.exprs[idx]
    var::Basic = data.vars[1]
    if n != 0x00000
        for i in 1:n
            ret = SymEngine.diff(ret, var)
        end
    end
    return ret
end

"""
    createPoly(data::QuadFreqData, n::UInt32, idx::UInt32)::Polynomial{Complex{BigFloat}}

Create a 2nd order `Polynomial` object in the variable `ω` by computing derivatives of `λ0` or `S0`.

# Input
- `data::QuadFreqData`: The spacetime data with the expressions to derivate.
- `n::UInt32`: The order of the derivative.
- `idx::UInt32`: The actual expression to derivate. If `idx = 0x00001` it computes a derivative of λ0. If `idx = 0x00002` it computes a derivative iof S0.

# Output
An object of type `Polynomial{Complex{BigFloat}}` containing the polynomial resulting from the derivation of the expression.
"""
function createPoly(data::QuadFreqData, n::UInt32, idx::UInt32)
    # Compute the derivative at x0 and expand it to obtain a polynomial in ω
    der::Basic = expand(subs(dnx(data, n, idx), data.vars[1] => data.x0))

    # Extract the factors of the polynomial in ω as BigFloats
    p0::Complex{BigFloat} = Complex{BigFloat}(N(coeff(der, data.vars[2], Basic(0))))
    p1::Complex{BigFloat} = Complex{BigFloat}(N(coeff(der, data.vars[2], Basic(1))))
    p2::Complex{BigFloat} = Complex{BigFloat}(N(coeff(der, data.vars[2], Basic(2))))

    poly = Polynomial{Complex{BigFloat}}([p0, p1, p2])

    return poly
end 
