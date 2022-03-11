using QuasinormalModes
using SymEngine

# ------------------------------------------------------------------
# 1. Analytic Schwarzschild Black Hole
# ------------------------------------------------------------------

struct SchwarzschildData{N,T} <: QuadraticEigenvalueProblem{N,T}
    nIter::N
    x0::T

    vars::Tuple{Basic, Basic}
    exprs::Tuple{Basic, Basic}
end

function SchwarzschildData(nIter::N, x0::T, l::N, s::N) where {N,T}
    vars = @vars x ω

    λ0 = (-1 + (2*im)*ω + x*(4 - 3*x + (4*im)*(-2 + x)*ω))/((-1 + x)^2*x)
    S0 = (l + l^2 + (-1 + s^2)*(-1 + x) + (4*im)*(-1 + x)*ω + 4*(-2 + x)*ω^2)/((-1 + x)^2*x)

    return SchwarzschildData{N,T}(nIter, x0, vars, (λ0, S0))
end

QuasinormalModes.λ0(d::SchwarzschildData{N,T}) where {N,T} = d.exprs[1]
QuasinormalModes.S0(d::SchwarzschildData{N,T}) where {N,T}  = d.exprs[2]

QuasinormalModes.get_niter(d::SchwarzschildData{N,T}) where {N,T} = d.nIter
QuasinormalModes.get_x0(d::SchwarzschildData{N,T}) where {N,T} = d.x0

QuasinormalModes.get_ODEvar(d::SchwarzschildData{N,T}) where {N,T} = d.vars[1]
QuasinormalModes.get_ODEeigen(d::SchwarzschildData{N,T}) where {N,T} = d.vars[2]

# ------------------------------------------------------------------
# 2. Numeric Schwarzschild Black Hole
# ------------------------------------------------------------------

struct NSchwarzschildData{N,T} <: NumericAIMProblem{N,T}
    nIter::N
    x0::T
    l::N
    s::N
end

function NSchwarzschildData(nIter::N, x0::T, l::N, s::N) where {N,T}
    return NSchwarzschildData{N,T}(nIter, x0, l, s)
end

QuasinormalModes.λ0(::NSchwarzschildData{N,T}) where {N,T} = (x,ω) -> (-1 + (2*im)*ω + x*(4 - 3*x + (4*im)*(-2 + x)*ω))/((-1 + x)^2*x)
QuasinormalModes.S0(d::NSchwarzschildData{N,T}) where {N,T} = (x,ω) -> (d.l + d.l^2 + (-1 + d.s^2)*(-1 + x) + (4*im)*(-1 + x)*ω + 4*(-2 + x)*ω^2)/((-1 + x)^2*x)

QuasinormalModes.get_niter(d::NSchwarzschildData{N,T}) where {N,T} = d.nIter
QuasinormalModes.get_x0(d::NSchwarzschildData{N,T}) where {N,T} = d.x0

# ------------------------------------------------------------------
# 3. Constructing problems and caches
# ------------------------------------------------------------------

p_ana = SchwarzschildData(0x00030, Complex(BigFloat("0.43"), BigFloat("0.0")), 0x00000, 0x00000);
p_num = NSchwarzschildData(0x00030, Complex(0.43, 0.0), 0x00000, 0x00000);

c_ana = AIMCache(p_ana)
c_num = AIMCache(p_num)

# ------------------------------------------------------------------
# 4. Computing quasinormal modes
# ------------------------------------------------------------------

m_ana = computeEigenvalues(Serial(), p_ana, c_ana)

function printQNMs(qnms, cutoff, instab)
    println("-"^165)
    println("|", " "^36, "Re(omega)", " "^36, " ", " "^36, "Im(omega)", " "^36, "|")
    println("-"^165)

    for qnm in qnms
        if real(qnm) > cutoff && ( instab ? true : imag(qnm) < big"0.0" )
        println(real(qnm), "    ", imag(qnm))
        end
    end
    
    println("-"^165)

    return nothing
end

sort!(m_ana, by = x -> imag(x))

println("Analytic results")
printQNMs(m_ana[1:5], 1.0e-10, false)

ev = computeEigenvalues(Serial(), p_num, c_num, Complex(0.22, -0.20), nls_xtol = 1.0e-10, nls_ftol = 1.0e-10)

println("Numeric results:")
println(ev)

