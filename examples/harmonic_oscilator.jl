using QuasinormalModes
using SymEngine
using TaylorSeries

#===============================================
== Analytic Harmonic Oscilator data structure ==
===============================================#

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

#========================================
== Numerc Schwarzschild data structure ==
========================================#

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

#==================================================
== Computing modes for the semi-analytic version ==
==================================================#
function semi_analytic()
    data = HarmonicOscilatorData(UInt32(50))
    modes = computeQNMs(data, plr_epsilon=BigFloat("1.0e-20"))

    sort!(modes, by = x -> real(x))
    
    for mode in modes
        println(round(Int64,real(mode)))
    end
end


#============================================
== Computing modes for the numeric version ==
============================================#
function numeric()
    data = NHarmonicOscilatorData(UInt32(50))
    mode = computeQNMs(data, Complex{BigFloat}(BigFloat("15.0"), BigFloat("1.0")), nls_xtol=BigFloat("1.0e-20"), nls_ftol=BigFloat("1.0e-20"))

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
    data = NHarmonicOscilatorData(UInt32(50))

    grid_start = Complex{BigFloat}(BigFloat("1.0"), BigFloat("-1.0"))
    grid_end = Complex{BigFloat}(BigFloat("20.0"), BigFloat("-0.01"))

    real_pts = 5
    imag_pts = 5

    grid = (grid_start, grid_end, real_pts, imag_pts)

    modes = modesInGrid(data, grid, xtol=BigFloat("1.0e-20"), ftol=BigFloat("1.0e-20"))

    sort!(modes, by = x -> real(x))
    
    for mode in modes
        println(round(Int64,real(mode)))
    end
end
