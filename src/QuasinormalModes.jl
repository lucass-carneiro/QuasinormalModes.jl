#= QuasinormalModes.jl
 = This package contains routines for computing the qusinormal modes (QNMs) of black holes in General Ralativity
 = using the "Assymptotic Iteration Method" (https://arxiv.org/abs/math-ph/0309066v1). using the implementation
 = based on the "improved" version of the AIM, described in (https://arxiv.org/abs/1111.5024).
 = In general this package can be used to compute the eigenvalues of any 2nd order ordinary differential equation
 = if they are written in the propepr format.
=#
__precompile__()

module QuasinormalModes

using SymEngine
using Polynomials
using PolynomialRoots
using TaylorSeries
using NLsolve
using CSV
using DataFrames

# Export QNM data abstract types
export QNMData
export GenFreqData
export QuadFreqData

# Export exceptions that can be used to construct user code
export AIMEvalPointException

# Mode-computing exports
export computeQNMs, modesInGrid

# Output exports
export filterQNMs!, printQNMs, saveQNMs

"""
Parent type of all QNM data types.
"""
abstract type QNMData end

"""
Parent type of all space-times that have QNM equations that are generic functions of ω and x.
"""
abstract type GenFreqData <: QNMData end

"""
Parent type of all space-times that have QNM equations that are quadratic in the frequency ω.
"""
abstract type QuadFreqData <: QNMData end

"""
Exception that can be thrown by the user when the evaluation point of the QNM data is out of the allowed range.
"""
abstract type AIMEvalPointException <: Exception end

"""
Exception that is thrown when the grid passed to modesInGrid is not properly formatted
"""
abstract type InvalidGridException <: Exception end

# Include core functionality
include("./core/aim_step.jl")
include("./core/analytic_tools.jl")
include("./core/compute_delta.jl")

# Include mode-computng functionality
include("./mode_computing/compute_qnms.jl")
include("./mode_computing/modes_in_grid.jl")

# Include output functionality
include("./output/filter_qnms.jl")
include("./output/print_qnms.jl")
include("./output/save_qnms.jl")

end # QuasinormalModes
