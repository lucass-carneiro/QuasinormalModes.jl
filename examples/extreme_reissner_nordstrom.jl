using QuasinormalModes
using SymEngine
using TaylorSeries

#=======================================================
== Analytic Extreme Reissner-Nordstrom data structure ==
=======================================================#

struct ExtremeRNData <: QuadFreqData
    nIter::UInt32
    x0::Float64

    M1::Float64
    l::UInt32

    vars::Tuple{Basic, Basic}
    exprs::Tuple{Basic, Basic}
    
    function ExtremeRNData(nIter::UInt32, M1::Float64, l::UInt32, x0::Float64 = 0.5)
        if x0 >= 1.0 || x0 <= 0.0
            error("x0 must be a number in the open interval 0 < x < 1.")
            throw(AIMEvalPointException)
        end
        
        vars = @vars x ω

        λ0 = (2*(-1 + 2*x)*((1 - x)*x + im*M1*(-1 + 2*(-1 + x)*x)*ω))/((-1 + x)^2*x^2)
        S0 = (l + l^2 - 2*(-1 + x)*x + (12*im)*M1*(-1 + x)*x*ω + 4*M1^2*(-3 + 2*x)*(1 + 2*x)ω^2)/((-1 + x)^2*x^2)

        return new(nIter, x0, M1, l, vars, (λ0, S0))
    end
end

#=====================================================
== Numerc Extreme Reissner-Nordstrom data structure ==
=====================================================#

struct NExtremeRNData <: GenFreqData
    nIter::UInt32
    x0::Float64

    M1::Float64
    l::UInt32
    
    function NExtremeRNData(nIter::UInt32, M1::Float64, l::UInt32, x0::Float64 = 0.5)

        if x0 >= 1.0 || x0 <= 0.0
            error("x0 must be a number in the open interval 0 < x < 1.")
            throw(AIMEvalPointException)
        end

        return new(nIter, x0, M1, l)
    end
end 

function (data::NExtremeRNData)(idx::UInt32,
                                x::TaylorSeries.Taylor1{Complex{BigFloat}},
                                ω::Complex{BigFloat}
                                )::Array{Complex{BigFloat},1}
    
    if idx == 0x00001
        ts = (2*(-1 + 2*x)*((1 - x)*x + im*data.M1*(-1 + 2*(-1 + x)*x)*ω))/((-1 + x)^2*x^2)
        return ts.coeffs
    elseif idx == 0x00002
        ts = (data.l + data.l^2 - 2*(-1 + x)*x + (12*im)*data.M1*(-1 + x)*x*ω + 4*data.M1^2*(-3 + 2*x)*(1 + 2*x)ω^2)/((-1 + x)^2*x^2)
        return ts.coeffs
    end
end

#==================================================
== Computing modes for the semi-analytic version ==
==================================================#
function semi_analytic()
    data = ExtremeRNData(UInt32(140), 1.0, UInt32(0))
    modes = computeQNMs(data, plr_epsilon=BigFloat("1.0e-200"))
    printQNMs(modes)
end


#============================================
== Computing modes for the numeric version ==
============================================#
function numeric()
    data = NExtremeRNData(UInt32(50), 1.0, UInt32(0))
    mode = computeQNMs(data, Complex{BigFloat}(BigFloat("0.13"), BigFloat("-0.09")), nls_xtol=BigFloat("1.0e-50"), nls_ftol=BigFloat("1.0e-50"))

    if mode.x_converged || mode.f_converged
        println(mode.zero[1], "    ", mode.zero[2])
    else
        println("Did not converge to any modes :-(")
    end
end


#============================
== Finding modes in a grid ==
============================#
function grid()
    data = NExtremeRNData(UInt32(50), 1.0, UInt32(0))

    grid_start = Complex{BigFloat}(BigFloat("0.01"), BigFloat("-1.0"))
    grid_end = Complex{BigFloat}(BigFloat("1.0"), BigFloat("-0.01"))

    real_pts = 5
    imag_pts = 5

    grid = (grid_start, grid_end, real_pts, imag_pts)

    modes = modesInGrid(data, grid, xtol=BigFloat("1.0e-50"), ftol=BigFloat("1.0e-50"))

    printQNMs(modes)
end
