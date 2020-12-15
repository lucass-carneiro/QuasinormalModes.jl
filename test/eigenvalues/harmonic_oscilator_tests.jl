using QuasinormalModes
using SymEngine
using TaylorSeries
using Test

#============================================================
== Numeric and semi-analytic Schwarzschild data structures ==
============================================================#
struct HarmonicOscilatorData <: QuadFreqData # Indicate that we will operate semi-analytically
    nIter::UInt32                            # The number of iteration the AIM will perform
    x0::Float64                              # The evaluation point for the AIM functions

    vars::Tuple{Basic, Basic}                # The variables x and ω as SymEngine expressions.
    exprs::Tuple{Basic, Basic}               # The functions λ0 and S0 as SymEngine expressions.
    
    function HarmonicOscilatorData(nIter::UInt32, x0::Float64 = 0.0)
	
        vars = @vars x ω
	
        λ0 = 2*x
        S0 = 1 - ω

        return new(nIter, x0, vars, (λ0, S0))
    end
end 

struct NHarmonicOscilatorData <: GenFreqData # Indicate that we will operate semi-analytically
    nIter::UInt32                            # The number of iteration the AIM will perform
    x0::Float64                              # The evaluation point for the AIM functions
    
    function NHarmonicOscilatorData(nIter::UInt32, x0::Float64 = 0.0)
        return new(nIter, x0)
    end
end 

function (data::NHarmonicOscilatorData)(idx::UInt32,
                                        x::TaylorSeries.Taylor1{Complex{BigFloat}},
                                        ω::Complex{BigFloat}
                                        )::Array{Complex{BigFloat},1}
    
    if idx == 0x00001
        ts = 2*x # The expression for λ0
        return ts.coeffs
    elseif idx == 0x00002
        ts = 1 - ω + x - x # The expression for S0. We add and subtract x in order for the taylor expansion to work
        return ts.coeffs
    end
end

#=====================================
== Tolerances and iteration numbers ==
=====================================#
const iter_low = UInt32(20)
const iter_high = UInt32(50)

const tol_low = BigFloat("1.0e-20")
const tol_high = BigFloat("1.0e-50")

#===================
== Test functions ==
===================#

function numeric_correctnes(iter, expected, tol, eps)
    data = NHarmonicOscilatorData(iter)

    mode = computeQNMs(data, expected, nls_xtol=tol, nls_ftol=tol)

    re_correct = abs(mode.zero[1] - real(expected)) < eps
    im_correct = abs(mode.zero[2] - imag(expected)) < eps

    if mode.x_converged || mode.f_converged
        return re_correct && im_correct
    else
        return false
    end
end

function analytic_correctnes(iter, expected, idx, eps)
    data = HarmonicOscilatorData(iter)

    modes = computeQNMs(data, plr_epsilon=BigFloat("1.0e-200"))
    sort!(modes, by = x -> real(x))

    if length(modes) != 0
        
        re_correct = abs(real(modes[idx]) - real(expected)) < eps
        im_correct = abs(imag(modes[idx]) - imag(expected)) < eps

        return re_correct && im_correct
    else
        return false
    end
end

#=====================
== Expected results ==
=====================#
const expected_0 = Complex{BigFloat}(BigFloat("1.0"), BigFloat("0.0"))
const expected_1 = Complex{BigFloat}(BigFloat("3.0"), BigFloat("0.0"))
const expected_2 = Complex{BigFloat}(BigFloat("5.0"), BigFloat("0.0"))

#==========
== Tests ==
==========#
@testset "Numeric correctness - low accuracy with n = 0" begin
    @test numeric_correctnes(iter_low, expected_0, tol_low, BigFloat("1.0e-10"))
    @test numeric_correctnes(iter_low, expected_1, tol_low, BigFloat("1.0e-10"))
    @test numeric_correctnes(iter_low, expected_2, tol_low, BigFloat("1.0e-10"))
end

@testset "Numeric correctness - high accuracy with n = 0" begin
    @test numeric_correctnes(iter_high, expected_0, tol_high, BigFloat("1.0e-10"))
    @test numeric_correctnes(iter_high, expected_1, tol_high, BigFloat("1.0e-10"))
    @test numeric_correctnes(iter_high, expected_2, tol_high, BigFloat("1.0e-10"))
end

@testset "Analytic correctness - low accuracy with n = 0" begin
    @test analytic_correctnes(iter_low, expected_0, 1, BigFloat("1.0e-10"))
    @test analytic_correctnes(iter_low, expected_1, 2, BigFloat("1.0e-10"))
    @test analytic_correctnes(iter_low, expected_2, 3, BigFloat("1.0e-10"))
end

@testset "Analytic correctness - high accuracy with n = 0" begin
    @test analytic_correctnes(iter_high, expected_0, 1, BigFloat("1.0e-10"))
    @test analytic_correctnes(iter_high, expected_1, 2, BigFloat("1.0e-10"))
    @test analytic_correctnes(iter_high, expected_2, 3, BigFloat("1.0e-10"))
end