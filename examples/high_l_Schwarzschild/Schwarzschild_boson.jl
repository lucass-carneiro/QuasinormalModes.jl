# ------------------------------------------------------------------
# Spin 0, 1 and 2 field on a Schwarzschild background
# ------------------------------------------------------------------

struct Schwarzschild_boson{N,T} <: NumericAIMProblem{N,T}
    nIter::N
    x0::T

    M::T
    l::N
    s::N
end

QuasinormalModes.λ0(d::Schwarzschild_boson{N,T}) where {N,T} = (x,ω) -> (4 * d.M * im * ω * (2 * x^2 - 4 * x + 1) - (1 - 3 * x) * (1 - x))/(x * (1 - x)^2)
QuasinormalModes.S0(d::Schwarzschild_boson{N,T}) where {N,T} = (x,ω) -> (16 * d.M^2 * ω^2 * (x - 2) - 8 * d.M * im * ω * (1 - x) + d.l * (d.l + 1) + (1 - d.s^2)*(1 - x))/(x * (1-x)^2)
QuasinormalModes.get_niter(d::Schwarzschild_boson{N,T}) where {N,T} = d.nIter
QuasinormalModes.get_x0(d::Schwarzschild_boson{N,T}) where {N,T} = d.x0
