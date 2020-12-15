using QuasinormalModes
using SymEngine
using TaylorSeries
using Test

#============================================================
== Numeric and semi-analytic Schwarzschild data structures ==
============================================================#
struct SchwarzschildData <: QuadFreqData
    nIter::UInt32
    x0::Float64

    l::UInt32
    s::UInt32
    vars::Tuple{Basic, Basic}
    exprs::Tuple{Basic, Basic}
    
    function SchwarzschildData(nIter::UInt32, l::UInt32, s::UInt32, x0::Float64 = 0.5)
	
        if x0 >= 1.0 || x0 <= 0.0
            error("x0 must be a number in the open interval 0 < x < 1.")
            throw(AIMEvalPointException)
        end

        vars = @vars x ω
	
        λ0 = (-1 + (2*im)*ω + x*(4 - 3*x + (4*im)*(-2 + x)*ω))/((-1 + x)^2*x)
        S0 = (l + l^2 + (-1 + s^2)*(-1 + x) + (4*im)*(-1 + x)*ω + 4*(-2 + x)*ω^2)/((-1 + x)^2*x)

        return new(nIter, x0, l, s, vars, (λ0, S0))
    end
end 

struct NSchwarzschildData <: GenFreqData
    nIter::UInt32
    x0::Float64

    l::UInt32
    s::UInt32
    
    function NSchwarzschildData(nIter::UInt32, l::UInt32, s::UInt32, x0::Float64 = 0.5)

        if x0 >= 1.0 || x0 <= 0.0
            error("x0 must be a number in the open interval 0 < x < 1.")
            throw(AIMEvalPointException)
        end

        return new(nIter, x0, l, s)
    end
end 

function (data::NSchwarzschildData)(idx::UInt32,
                                    x::TaylorSeries.Taylor1{Complex{BigFloat}},
                                    ω::Complex{BigFloat}
                                    )::Array{Complex{BigFloat},1}
    
    if idx == 0x00001
        ts = (-1 + (2*im)*ω + x*(4 - 3*x + (4*im)*(-2 + x)*ω))/((-1 + x)^2*x)
        return ts.coeffs
    elseif idx == 0x00002
        ts =(data.l + data.l^2 + (-1 + data.s^2)*(-1 + x) + (4*im)*(-1 + x)*ω + 4*(-2 + x)*ω^2)/((-1 + x)^2*x)
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

function numeric_correctnes(iter, l, s, expected, tol, eps)
    data = NSchwarzschildData(iter, UInt32(l), UInt32(s))

    mode = computeQNMs(data, expected, nls_xtol=tol, nls_ftol=tol)

    re_correct = abs(mode.zero[1] - real(expected)) < eps
    im_correct = abs(mode.zero[2] - imag(expected)) < eps

    if mode.x_converged || mode.f_converged
        return re_correct && im_correct
    else
        return false
    end
end

function analytic_correctnes(iter, l, s, expected, eps)
    data = SchwarzschildData(iter, UInt32(l), UInt32(s))

    modes = computeQNMs(data, plr_epsilon=BigFloat("1.0e-200"))
    filterQNMs!(modes)

    if length(modes) != 0
        
        re_correct = abs(real(modes[end]) - real(expected)) < eps
        im_correct = abs(imag(modes[end]) - imag(expected)) < eps

        return re_correct && im_correct
    else
        return false
    end
end

#=====================
== Expected results ==
=====================#
const expected_000 = Complex{BigFloat}(BigFloat("0.2209098781608393"), BigFloat("-0.2097914341737619"))
const expected_011 = Complex{BigFloat}(BigFloat("0.4965265283562174"), BigFloat("-0.1849754359058844"))
const expected_022 = Complex{BigFloat}(BigFloat("0.7473433688360838"), BigFloat("-0.1779246313778714"))

#==========
== Tests ==
==========#
@testset "Numeric correctness - low accuracy with n = 0" begin
    @test numeric_correctnes(iter_low, 0, 0, expected_000, tol_low, BigFloat("1.0e-2"))
    @test numeric_correctnes(iter_low, 1, 1, expected_011, tol_low, BigFloat("1.0e-2"))
    @test numeric_correctnes(iter_low, 2, 2, expected_022, tol_low, BigFloat("1.0e-2"))
end

@testset "Numeric correctness - high accuracy with n = 0" begin
    @test numeric_correctnes(iter_high, 0, 0, expected_000, tol_high, BigFloat("1.0e-3"))
    @test numeric_correctnes(iter_high, 1, 1, expected_011, tol_high, BigFloat("1.0e-3"))
    @test numeric_correctnes(iter_high, 2, 2, expected_022, tol_high, BigFloat("1.0e-3"))
end

@testset "Analytic correctness - low accuracy with n = 0" begin
    @test analytic_correctnes(iter_low, 0, 0, expected_000, BigFloat("1.0e-2"))
    @test analytic_correctnes(iter_low, 1, 1, expected_011, BigFloat("1.0e-2"))
    @test analytic_correctnes(iter_low, 2, 2, expected_022, BigFloat("1.0e-2"))
end

@testset "Analytic correctness - high accuracy with n = 0" begin
    @test analytic_correctnes(iter_high, 0, 0, expected_000, BigFloat("1.0e-3"))
    @test analytic_correctnes(iter_high, 1, 1, expected_011, BigFloat("1.0e-3"))
    @test analytic_correctnes(iter_high, 2, 2, expected_022, BigFloat("1.0e-3"))
end