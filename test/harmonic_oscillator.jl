using QuasinormalModes
using SymEngine
using Test

# ------------------------------------------------------------------
# 1. Analytic harmonic oscillator
# ------------------------------------------------------------------

struct HarmonicOscilatorData{N,T} <: QuadraticEigenvalueProblem{N,T}
    nIter::N
    x0::T

    vars::Tuple{Basic,Basic}
    exprs::Tuple{Basic,Basic}
end

function HarmonicOscilatorData(nIter::N, x0::T) where {N,T}

    vars = @vars x ω

    λ0 = 2 * x
    S0 = 1 - ω

    return HarmonicOscilatorData{N,T}(nIter, x0, vars, (λ0, S0))
end

QuasinormalModes.λ0(d::HarmonicOscilatorData{N,T}) where {N,T} = d.exprs[1]
QuasinormalModes.S0(d::HarmonicOscilatorData{N,T}) where {N,T} = d.exprs[2]

QuasinormalModes.get_niter(d::HarmonicOscilatorData{N,T}) where {N,T} = d.nIter
QuasinormalModes.get_x0(d::HarmonicOscilatorData{N,T}) where {N,T} = d.x0

QuasinormalModes.get_ODEvar(d::HarmonicOscilatorData{N,T}) where {N,T} = d.vars[1]
QuasinormalModes.get_ODEeigen(d::HarmonicOscilatorData{N,T}) where {N,T} = d.vars[2]

# ------------------------------------------------------------------
# 2. Numeric harmonic oscillator
# ------------------------------------------------------------------

struct NHarmonicOscilatorData{N,T} <: NumericAIMProblem{N,T}
    nIter::N
    x0::T
end

function NHarmonicOscilatorData(nIter::N, x0::T) where {N,T}
    return NHarmonicOscilatorData{N,T}(nIter, x0)
end

QuasinormalModes.λ0(::NHarmonicOscilatorData{N,T}) where {N,T} = (x, ω) -> 2 * x
QuasinormalModes.S0(::NHarmonicOscilatorData{N,T}) where {N,T} = (x, ω) -> 1 - ω + x - x

QuasinormalModes.get_niter(d::NHarmonicOscilatorData{N,T}) where {N,T} = d.nIter
QuasinormalModes.get_x0(d::NHarmonicOscilatorData{N,T}) where {N,T} = d.x0

# ------------------------------------------------------------------
# 3. Serial Numeric Tests
# ------------------------------------------------------------------

@testset "Serial Numeric correctness - 10 iterations" begin
    p = NHarmonicOscilatorData(0x0000A, 0.5)
    c = AIMCache(p)
    ev = eigenvaluesInGrid(Serial(), p, c, (0.0, 20.0))
    ce = [computeEigenvalues(Serial(), p, c, guess) for guess in 1.0:2.0:19.0]

    @test length(ev) == 10
    @test round(ev[1]) == 1
    @test round(ev[2]) == 3
    @test round(ev[3]) == 5
    @test round(ev[4]) == 7
    @test round(ev[5]) == 9
    @test round(ev[6]) == 11
    @test round(ev[7]) == 13
    @test round(ev[8]) == 15
    @test round(ev[9]) == 17
    @test round(ev[10]) == 19

    @test length(ce) == 10
    @test round(ce[1]) == 1
    @test round(ce[2]) == 3
    @test round(ce[3]) == 5
    @test round(ce[4]) == 7
    @test round(ce[5]) == 9
    @test round(ce[6]) == 11
    @test round(ce[7]) == 13
    @test round(ce[8]) == 15
    @test round(ce[9]) == 17
    @test round(ce[10]) == 19
end

@testset "Serial Numeric correctness - 20 iterations" begin
    p = NHarmonicOscilatorData(0x00014, 0.5)
    c = AIMCache(p)
    ev = eigenvaluesInGrid(Serial(), p, c, (0.0, 20.0))
    ce = [computeEigenvalues(Serial(), p, c, guess) for guess in 1.0:2.0:19.0]

    @test length(ev) == 10
    @test round(ev[1]) == 1
    @test round(ev[2]) == 3
    @test round(ev[3]) == 5
    @test round(ev[4]) == 7
    @test round(ev[5]) == 9
    @test round(ev[6]) == 11
    @test round(ev[7]) == 13
    @test round(ev[8]) == 15
    @test round(ev[9]) == 17
    @test round(ev[10]) == 19

    @test length(ce) == 10
    @test round(ce[1]) == 1
    @test round(ce[2]) == 3
    @test round(ce[3]) == 5
    @test round(ce[4]) == 7
    @test round(ce[5]) == 9
    @test round(ce[6]) == 11
    @test round(ce[7]) == 13
    @test round(ce[8]) == 15
    @test round(ce[9]) == 17
    @test round(ce[10]) == 19
end

@testset "Serial Numeric correctness - 50 iterations" begin
    p = NHarmonicOscilatorData(0x00032, 0.5)
    c = AIMCache(p)
    ev = eigenvaluesInGrid(Serial(), p, c, (0.0, 20.0))
    ce = [computeEigenvalues(Serial(), p, c, guess) for guess in 1.0:2.0:19.0]

    @test length(ev) == 10
    @test round(ev[1]) == 1
    @test round(ev[2]) == 3
    @test round(ev[3]) == 5
    @test round(ev[4]) == 7
    @test round(ev[5]) == 9
    @test round(ev[6]) == 11
    @test round(ev[7]) == 13
    @test round(ev[8]) == 15
    @test round(ev[9]) == 17
    @test round(ev[10]) == 19

    @test length(ce) == 10
    @test round(ce[1]) == 1
    @test round(ce[2]) == 3
    @test round(ce[3]) == 5
    @test round(ce[4]) == 7
    @test round(ce[5]) == 9
    @test round(ce[6]) == 11
    @test round(ce[7]) == 13
    @test round(ce[8]) == 15
    @test round(ce[9]) == 17
    @test round(ce[10]) == 19
end

@testset "Serial Numeric correctness - 100 iterations" begin
    p = NHarmonicOscilatorData(0x00064, 0.5)
    c = AIMCache(p)
    ev = eigenvaluesInGrid(Serial(), p, c, (0.0, 20.0))
    ce = [computeEigenvalues(Serial(), p, c, guess) for guess in 1.0:2.0:19.0]

    @test length(ev) == 10
    @test round(ev[1]) == 1
    @test round(ev[2]) == 3
    @test round(ev[3]) == 5
    @test round(ev[4]) == 7
    @test round(ev[5]) == 9
    @test round(ev[6]) == 11
    @test round(ev[7]) == 13
    @test round(ev[8]) == 15
    @test round(ev[9]) == 17
    @test round(ev[10]) == 19

    @test length(ce) == 10
    @test round(ce[1]) == 1
    @test round(ce[2]) == 3
    @test round(ce[3]) == 5
    @test round(ce[4]) == 7
    @test round(ce[5]) == 9
    @test round(ce[6]) == 11
    @test round(ce[7]) == 13
    @test round(ce[8]) == 15
    @test round(ce[9]) == 17
    @test round(ce[10]) == 19
end

@testset "Serial Numeric correctness - 25 iterations with BigFloat" begin
    p = NHarmonicOscilatorData(0x00019, BigFloat("0.5"))
    c = AIMCache(p)
    ev = eigenvaluesInGrid(Serial(), p, c, (BigFloat("0.0"), BigFloat("20.0")))
    ce = [computeEigenvalues(Serial(), p, c, guess) for guess in BigFloat("1.0"):BigFloat("2.0"):BigFloat("19.0")]

    @test length(ev) == 10
    @test round(ev[1]) == 1
    @test round(ev[2]) == 3
    @test round(ev[3]) == 5
    @test round(ev[4]) == 7
    @test round(ev[5]) == 9
    @test round(ev[6]) == 11
    @test round(ev[7]) == 13
    @test round(ev[8]) == 15
    @test round(ev[9]) == 17
    @test round(ev[10]) == 19

    @test length(ce) == 10
    @test round(ce[1]) == 1
    @test round(ce[2]) == 3
    @test round(ce[3]) == 5
    @test round(ce[4]) == 7
    @test round(ce[5]) == 9
    @test round(ce[6]) == 11
    @test round(ce[7]) == 13
    @test round(ce[8]) == 15
    @test round(ce[9]) == 17
    @test round(ce[10]) == 19

end

# ------------------------------------------------------------------
# 4. Serial Analytic Tests
# ------------------------------------------------------------------

@testset "Serial Analytic correctness - 10 iterations" begin
    p = HarmonicOscilatorData(0x0000A, 0.5)
    c = AIMCache(p)
    ev = computeEigenvalues(Serial(), p, c)

    @test length(ev) >= 10
    @test round(real(ev[end])) == 1
    @test round(real(ev[end-1])) == 3
    @test round(real(ev[end-2])) == 5
    @test round(real(ev[end-3])) == 7
    @test round(real(ev[end-4])) == 9
    @test round(real(ev[end-5])) == 11
    @test round(real(ev[end-6])) == 13
    @test round(real(ev[end-7])) == 15
    @test round(real(ev[end-8])) == 17
    @test round(real(ev[end-9])) == 19
end

@testset "Serial Analytic correctness - 20 iterations" begin
    p = HarmonicOscilatorData(0x00014, 0.5)
    c = AIMCache(p)
    ev = computeEigenvalues(Serial(), p, c)

    @test length(ev) >= 10
    @test round(real(ev[end])) == 1
    @test round(real(ev[end-1])) == 3
    @test round(real(ev[end-2])) == 5
    @test round(real(ev[end-3])) == 7
    @test round(real(ev[end-4])) == 9
    @test round(real(ev[end-5])) == 11
    @test round(real(ev[end-6])) == 13
    @test round(real(ev[end-7])) == 15
    @test round(real(ev[end-8])) == 17
    @test round(real(ev[end-9])) == 19
end

@testset "Serial Analytic correctness - 25 iterations with BigFloat" begin
    p = HarmonicOscilatorData(0x00019, BigFloat("0.5"))
    c = AIMCache(p)
    ev = computeEigenvalues(Serial(), p, c)

    @test length(ev) >= 10
    @test round(real(ev[end])) == 1
    @test round(real(ev[end-1])) == 3
    @test round(real(ev[end-2])) == 5
    @test round(real(ev[end-3])) == 7
    @test round(real(ev[end-4])) == 9
    @test round(real(ev[end-5])) == 11
    @test round(real(ev[end-6])) == 13
    @test round(real(ev[end-7])) == 15
end

# ------------------------------------------------------------------
# 5. Threaded Numeric Tests
# ------------------------------------------------------------------

@testset "Threaded Numeric correctness - 10 iterations" begin
    p = NHarmonicOscilatorData(0x0000A, 0.5)
    c = AIMCache(p)
    ev = eigenvaluesInGrid(Threaded(), p, c, (0.0, 20.0))
    ce = [computeEigenvalues(Threaded(), p, c, guess) for guess in 1.0:2.0:19.0]

    @test length(ev) == 10
    @test round(ev[1]) == 1
    @test round(ev[2]) == 3
    @test round(ev[3]) == 5
    @test round(ev[4]) == 7
    @test round(ev[5]) == 9
    @test round(ev[6]) == 11
    @test round(ev[7]) == 13
    @test round(ev[8]) == 15
    @test round(ev[9]) == 17
    @test round(ev[10]) == 19

    @test length(ce) == 10
    @test round(ce[1]) == 1
    @test round(ce[2]) == 3
    @test round(ce[3]) == 5
    @test round(ce[4]) == 7
    @test round(ce[5]) == 9
    @test round(ce[6]) == 11
    @test round(ce[7]) == 13
    @test round(ce[8]) == 15
    @test round(ce[9]) == 17
    @test round(ce[10]) == 19
end

@testset "Threaded Numeric correctness - 20 iterations" begin
    p = NHarmonicOscilatorData(0x00014, 0.5)
    c = AIMCache(p)
    ev = eigenvaluesInGrid(Threaded(), p, c, (0.0, 20.0))
    ce = [computeEigenvalues(Threaded(), p, c, guess) for guess in 1.0:2.0:19.0]

    @test length(ev) == 10
    @test round(ev[1]) == 1
    @test round(ev[2]) == 3
    @test round(ev[3]) == 5
    @test round(ev[4]) == 7
    @test round(ev[5]) == 9
    @test round(ev[6]) == 11
    @test round(ev[7]) == 13
    @test round(ev[8]) == 15
    @test round(ev[9]) == 17
    @test round(ev[10]) == 19

    @test length(ce) == 10
    @test round(ce[1]) == 1
    @test round(ce[2]) == 3
    @test round(ce[3]) == 5
    @test round(ce[4]) == 7
    @test round(ce[5]) == 9
    @test round(ce[6]) == 11
    @test round(ce[7]) == 13
    @test round(ce[8]) == 15
    @test round(ce[9]) == 17
    @test round(ce[10]) == 19
end

@testset "Threaded Numeric correctness - 50 iterations" begin
    p = NHarmonicOscilatorData(0x00032, 0.5)
    c = AIMCache(p)
    ev = eigenvaluesInGrid(Threaded(), p, c, (0.0, 20.0))
    ce = [computeEigenvalues(Threaded(), p, c, guess) for guess in 1.0:2.0:19.0]

    @test length(ev) == 10
    @test round(ev[1]) == 1
    @test round(ev[2]) == 3
    @test round(ev[3]) == 5
    @test round(ev[4]) == 7
    @test round(ev[5]) == 9
    @test round(ev[6]) == 11
    @test round(ev[7]) == 13
    @test round(ev[8]) == 15
    @test round(ev[9]) == 17
    @test round(ev[10]) == 19

    @test length(ce) == 10
    @test round(ce[1]) == 1
    @test round(ce[2]) == 3
    @test round(ce[3]) == 5
    @test round(ce[4]) == 7
    @test round(ce[5]) == 9
    @test round(ce[6]) == 11
    @test round(ce[7]) == 13
    @test round(ce[8]) == 15
    @test round(ce[9]) == 17
    @test round(ce[10]) == 19
end

@testset "Threaded Numeric correctness - 100 iterations" begin
    p = NHarmonicOscilatorData(0x00064, 0.5)
    c = AIMCache(p)
    ev = eigenvaluesInGrid(Threaded(), p, c, (0.0, 20.0))
    ce = [computeEigenvalues(Threaded(), p, c, guess) for guess in 1.0:2.0:19.0]

    @test length(ev) == 10
    @test round(ev[1]) == 1
    @test round(ev[2]) == 3
    @test round(ev[3]) == 5
    @test round(ev[4]) == 7
    @test round(ev[5]) == 9
    @test round(ev[6]) == 11
    @test round(ev[7]) == 13
    @test round(ev[8]) == 15
    @test round(ev[9]) == 17
    @test round(ev[10]) == 19

    @test length(ce) == 10
    @test round(ce[1]) == 1
    @test round(ce[2]) == 3
    @test round(ce[3]) == 5
    @test round(ce[4]) == 7
    @test round(ce[5]) == 9
    @test round(ce[6]) == 11
    @test round(ce[7]) == 13
    @test round(ce[8]) == 15
    @test round(ce[9]) == 17
    @test round(ce[10]) == 19
end

@testset "Threaded Numeric correctness - 25 iterations with BigFloat" begin
    p = NHarmonicOscilatorData(0x00019, BigFloat("0.5"))
    c = AIMCache(p)
    ev = eigenvaluesInGrid(Threaded(), p, c, (big"0.0", big"20.0"))
    ce = [computeEigenvalues(Threaded(), p, c, guess) for guess in BigFloat("1.0"):BigFloat("2.0"):BigFloat("19.0")]

    @test length(ev) == 10
    @test round(ev[1]) == 1
    @test round(ev[2]) == 3
    @test round(ev[3]) == 5
    @test round(ev[4]) == 7
    @test round(ev[5]) == 9
    @test round(ev[6]) == 11
    @test round(ev[7]) == 13
    @test round(ev[8]) == 15
    @test round(ev[9]) == 17
    @test round(ev[10]) == 19

    @test length(ce) == 10
    @test round(ce[1]) == 1
    @test round(ce[2]) == 3
    @test round(ce[3]) == 5
    @test round(ce[4]) == 7
    @test round(ce[5]) == 9
    @test round(ce[6]) == 11
    @test round(ce[7]) == 13
    @test round(ce[8]) == 15
    @test round(ce[9]) == 17
    @test round(ce[10]) == 19
end

# ------------------------------------------------------------------
# 5. Threaded Analytic Tests
# ------------------------------------------------------------------

@testset "Threaded Analytic correctness - 10 iterations" begin
    p = HarmonicOscilatorData(0x0000A, 0.5)
    c = AIMCache(p)
    ev = computeEigenvalues(Threaded(), p, c)

    @test length(ev) >= 10
    @test round(real(ev[end])) == 1
    @test round(real(ev[end-1])) == 3
    @test round(real(ev[end-2])) == 5
    @test round(real(ev[end-3])) == 7
    @test round(real(ev[end-4])) == 9
    @test round(real(ev[end-5])) == 11
    @test round(real(ev[end-6])) == 13
    @test round(real(ev[end-7])) == 15
    @test round(real(ev[end-8])) == 17
    @test round(real(ev[end-9])) == 19
end

@testset "Threaded Analytic correctness - 20 iterations" begin
    p = HarmonicOscilatorData(0x00014, 0.5)
    c = AIMCache(p)
    ev = computeEigenvalues(Threaded(), p, c)

    @test length(ev) >= 10
    @test round(real(ev[end])) == 1
    @test round(real(ev[end-1])) == 3
    @test round(real(ev[end-2])) == 5
    @test round(real(ev[end-3])) == 7
    @test round(real(ev[end-4])) == 9
    @test round(real(ev[end-5])) == 11
    @test round(real(ev[end-6])) == 13
    @test round(real(ev[end-7])) == 15
    @test round(real(ev[end-8])) == 17
    @test round(real(ev[end-9])) == 19
end

@testset "Threaded Analytic correctness - 25 iterations with BigFloat" begin
    p = HarmonicOscilatorData(0x00019, BigFloat("0.5"))
    c = AIMCache(p)
    ev = computeEigenvalues(Threaded(), p, c)

    @test length(ev) >= 10
    @test round(real(ev[end])) == 1
    @test round(real(ev[end-1])) == 3
    @test round(real(ev[end-2])) == 5
    @test round(real(ev[end-3])) == 7
    @test round(real(ev[end-4])) == 9
    @test round(real(ev[end-5])) == 11
    @test round(real(ev[end-6])) == 13
    @test round(real(ev[end-7])) == 15
end
