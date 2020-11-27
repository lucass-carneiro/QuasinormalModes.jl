using QuasinormalModes
using SymEngine
using TaylorSeries

#==========================================
== Analytic Poschl Teller data structure ==
==========================================#

struct PoschlTellerData <: QuadFreqData # Indicate that we will operate semi-analytically
    nIter::UInt32                       # The number of iteration the AIM will perform
    x0::Float64                         # The evaluation point for the AIM functions

    vars::Tuple{Basic, Basic}           # The variables x and ω as SymEngine expressions.
    exprs::Tuple{Basic, Basic}          # The functions λ0 and S0 as SymEngine expressions.
    
    function PoschlTellerData(nIter::UInt32, x0::Float64 = 0.5)
	
        vars = @vars x ω
	
        λ0 = (2*x*(1-im*ω))/(1-x^2)
        S0 = (1 - 2*im*ω - 2*ω^2)/(2*(1-x^2))

        return new(nIter, x0, vars, (λ0, S0))
    end
end 

#========================================
== Numerc Schwarzschild data structure ==
========================================#

struct NPoschlTellerData <: GenFreqData # Indicate that we will operate semi-analytically
    nIter::UInt32                       # The number of iteration the AIM will perform
    x0::Float64                         # The evaluation point for the AIM functions
    
    function NPoschlTellerData(nIter::UInt32, x0::Float64 = 0.5)
        return new(nIter, x0)
    end
end 

function (data::NPoschlTellerData)(idx::UInt32,
                                   x::TaylorSeries.Taylor1{Complex{BigFloat}},
                                   ω::Complex{BigFloat}
                                   )::Array{Complex{BigFloat},1}
    
    if idx == 0x00001
        ts = (2*x*(1-im*ω))/(1-x^2) # The expression for λ0
        return ts.coeffs
    elseif idx == 0x00002
        ts = (1 - 2*im*ω - 2*ω^2)/(2*(1-x^2)) # The expression for S0. We add and subtract x in order for the taylor expansion to work
        return ts.coeffs
    end
end

#==================================================
== Computing modes for the semi-analytic version ==
==================================================#
function semi_analytic()
    data = PoschlTellerData(UInt32(50), 0.0)
    modes = computeQNMs(data, plr_epsilon=BigFloat("1.0e-200"))
    printQNMs(modes)
end


#============================================
== Computing modes for the numeric version ==
============================================#
function numeric()
    data = NPoschlTellerData(UInt32(50), 0.0)
    mode = computeQNMs(data, Complex{BigFloat}(BigFloat("0.5"), BigFloat("2.5")), nls_xtol=BigFloat("1.0e-50"), nls_ftol=BigFloat("1.0e-50"))

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
    data = NPoschlTellerData(UInt32(50), 0.0)

    grid_start = Complex{BigFloat}(BigFloat("0.45"), BigFloat("-10.5"))
    grid_end = Complex{BigFloat}(BigFloat("0.55"), BigFloat("-0.5"))

    real_pts = 2
    imag_pts = 10

    grid = (grid_start, grid_end, real_pts, imag_pts)

    modes = modesInGrid(data, grid, xtol=BigFloat("1.0e-50"), ftol=BigFloat("1.0e-50"))

    printQNMs(modes)
end
