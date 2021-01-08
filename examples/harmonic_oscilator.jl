using QuasinormalModes
using SymEngine

# ------------------------------------------------------------------
# 1. Analytic harmonic oscilator
# ------------------------------------------------------------------

struct HarmonicOscilatorData{N,T} <: QuadraticEigenvalueProblem{N,T}
    nIter::N
    x0::T

    vars::Tuple{Basic, Basic}
    exprs::Tuple{Basic, Basic}
end

function HarmonicOscilatorData(nIter::N, x0::T) where {N,T}
	
    vars = @vars x ω

    λ0 = 2*x
    S0 = 1 - ω

    return HarmonicOscilatorData{N,T}(nIter, x0, vars, (λ0, S0))
end

QuasinormalModes.λ0(d::HarmonicOscilatorData{N,T}) where {N,T} = d.exprs[1]
QuasinormalModes.S0(d::HarmonicOscilatorData{N,T}) where {N,T}  = d.exprs[2]

QuasinormalModes.get_niter(d::HarmonicOscilatorData{N,T}) where {N,T} = d.nIter
QuasinormalModes.get_x0(d::HarmonicOscilatorData{N,T}) where {N,T} = d.x0

QuasinormalModes.get_ODEvar(d::HarmonicOscilatorData{N,T}) where {N,T} = d.vars[1]
QuasinormalModes.get_ODEeigen(d::HarmonicOscilatorData{N,T}) where {N,T} = d.vars[2]

# ------------------------------------------------------------------
# 2. Numeric harmonic oscilator
# ------------------------------------------------------------------

struct NHarmonicOscilatorData{N,T} <: NumericAIMProblem{N,T}
    nIter::N
    x0::T
end

function NHarmonicOscilatorData(nIter::N, x0::T) where {N,T}
    return NHarmonicOscilatorData{N,T}(nIter, x0)
end

QuasinormalModes.λ0(::NHarmonicOscilatorData{N,T}) where {N,T} = (x,ω) -> 2*x
QuasinormalModes.S0(::NHarmonicOscilatorData{N,T}) where {N,T} = (x,ω) -> 1 - ω + x - x

QuasinormalModes.get_niter(d::NHarmonicOscilatorData{N,T}) where {N,T} = d.nIter
QuasinormalModes.get_x0(d::NHarmonicOscilatorData{N,T}) where {N,T} = d.x0

# ------------------------------------------------------------------
# 3. Constructing problems and caches
# ------------------------------------------------------------------

p_ana = HarmonicOscilatorData(0x0000A, 0.5);
p_num = NHarmonicOscilatorData(0x0000A, 0.5);

c_ana = AIMCache(p_ana)
c_num = AIMCache(p_num)

# ------------------------------------------------------------------
# 4. Computing quasinormal modes
# ------------------------------------------------------------------

ev_ana = computeEigenvalues(p_ana, c_ana)
ev_num = eigenvaluesInGrid(p_num, c_num, (0.0, 21.0))

function printEigen(eigenvalues)
    println("--------------------------------------")

    for i in eachindex(eigenvalues)
        println("n = $i, ω = $(eigenvalues[i])")
    end
    
    println("--------------------------------------")

    return nothing
end

println("Analytic results")
printEigen(reverse!(ev_ana))

println("Numeric results")
printEigen(ev_num)
