using QuasinormalModes
using SymEngine
using TaylorSeries

#==========================================
== Analytic Schwarzschild data structure ==
==========================================#

struct SchwarzschildData <: QuadFreqData # Indicate that we will operate semi-analytically
    nIter::UInt32                        # The number of iteration the AIM will perform
    x0::Float64                          # The evaluation point for the AIM functions

    l::UInt32                            # The angualar number.
    s::UInt32                            # The perturbation spin.
    vars::Tuple{Basic, Basic}            # The variables x and ω as SymEngine expressions.
    exprs::Tuple{Basic, Basic}           # The functions λ0 and S0 as SymEngine expressions.
    
    function SchwarzschildData(nIter::UInt32, l::UInt32, s::UInt32, x0::Float64 = 0.5)
	
        # Make sure that teh evaluation point is inside the ODE interval.
        # The ODE is singular at the endpoints so they are not allowed.
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

#========================================
== Numerc Schwarzschild data structure ==
========================================#

struct NSchwarzschildData <: GenFreqData # Indicate that we will operate semi-analytically
    nIter::UInt32                        # The number of iteration the AIM will perform
    x0::Float64                          # The evaluation point for the AIM functions

    l::UInt32                            # The angualar number.
    s::UInt32                            # The perturbation spin.
    
    function NSchwarzschildData(nIter::UInt32, l::UInt32, s::UInt32, x0::Float64 = 0.5)
        
        # Make sure that teh evaluation point is inside the ODE interval.
        # The ODE is singular at the endpoints so they are not allowed.
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
        ts = (-1 + (2*im)*ω + x*(4 - 3*x + (4*im)*(-2 + x)*ω))/((-1 + x)^2*x) # The expression for λ0
        return ts.coeffs
    elseif idx == 0x00002
        ts =(data.l + data.l^2 + (-1 + data.s^2)*(-1 + x) + (4*im)*(-1 + x)*ω + 4*(-2 + x)*ω^2)/((-1 + x)^2*x) # The expression for S0
        return ts.coeffs
    end
end

#==================================================
== Computing modes for the semi-analytic version ==
==================================================#
function semi_analytic()
    data = SchwarzschildData(UInt32(50), UInt32(0), UInt32(0))
    modes = computeQNMs(data, plr_epsilon=BigFloat("1.0e-20"))
    printQNMs(modes)
end


#============================================
== Computing modes for the numeric version ==
============================================#
function numeric()
    data = NSchwarzschildData(UInt32(50), UInt32(0), UInt32(0))
    mode = computeQNMs(data, Complex{BigFloat}(BigFloat("0.220"), BigFloat("-0.209")), nls_xtol=BigFloat("1.0e-20"), nls_ftol=BigFloat("1.0e-20"))

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
    data = NSchwarzschildData(UInt32(50), UInt32(0), UInt32(0))

    grid_start = Complex{BigFloat}(BigFloat("0.01"), BigFloat("-1.0"))
    grid_end = Complex{BigFloat}(BigFloat("1.0"), BigFloat("-0.01"))

    real_pts = 5
    imag_pts = 5

    grid = (grid_start, grid_end, real_pts, imag_pts)

    modes = modesInGrid(data, grid, xtol=BigFloat("1.0e-20"), ftol=BigFloat("1.0e-20"))

    printQNMs(modes)
end
