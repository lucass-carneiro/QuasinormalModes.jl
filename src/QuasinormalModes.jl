__precompile__()

"""
This package contains routines for computing eigenvalues of second
order ordinary differential equations and in particular the
quasinormal modes (QNMs) of black holes in General Relativity using
the "Asymptotic Iteration Method" [1] using the implementation
based on the "improved" version of the AIM, described in [2].

References:

[1] (https://arxiv.org/abs/math-ph/0309066v1)
[2] (https://arxiv.org/abs/1111.5024)
"""
module QuasinormalModes

# ------------------------------------------------------------------
# 1. Imports
# ------------------------------------------------------------------

using SymEngine
using Polynomials
using PolynomialRoots
using TaylorSeries
using NLsolve
using Roots

# ------------------------------------------------------------------
# 2. Public API (Exports)
# ------------------------------------------------------------------

# --- Type hierarchy ---
export AIMProblem
export AnalyticAIMProblem
export NumericAIMProblem
export QuadraticEigenvalueProblem

# --- Extra (concrete) types ---
export AIMCache

# --- Mandatory methods ---
export λ0, S0, get_niter, get_x0
export get_ODEvar, get_ODEeigen

# --- AIM methods ---
export computeDelta!
export computeEigenvalues
export eigenvaluesInGrid

# ------------------------------------------------------------------
# 3. Type hierarchy
# ------------------------------------------------------------------

"""
Parent super-type of all problems that can be solved using the AIM.
"""
abstract type AIMProblem{N <: Unsigned, T <: Number} end

"""
Parent super-type of all problems that can be solved using the AIM semi-analytically.
"""
abstract type AnalyticAIMProblem{N <: Unsigned, T <: Number} <: AIMProblem{N,T} end

"""
Parent super-type of all problems that can be solved using the AIM numerically.
"""
abstract type NumericAIMProblem{N <: Unsigned, T <: Number} <: AIMProblem{N,T} end

"""
Parent super-type of all problems whose eigenvalue is a quadratic polynomial.
"""
abstract type QuadraticEigenvalueProblem{N <: Unsigned, T <: Number} <: AnalyticAIMProblem{N,T} end

# ------------------------------------------------------------------
# 4. Traits
# ------------------------------------------------------------------

"""
Super-type of traits describing the analyticity of eigenvalue problems.
"""
abstract type AnalyticityTrait end

"""
All problems with eigenvalues that *can* be described by analytic functions have this trait.
"""
struct IsAnalytic <: AnalyticityTrait end

"""
All problems with eigenvalues that *can't* be described by analytic functions have this trait.
"""
struct IsNumeric <: AnalyticityTrait end

"""
The default trait of AIMProblem(s).
"""
AnalyticityTrait(::Type{<:AIMProblem}) = IsNumeric()

"""
The trait of AnalyticAIMProblem(s).
"""
AnalyticityTrait(::Type{<:AnalyticAIMProblem}) = IsAnalytic()

"""
The trait of NumericAIMProblem(s).
"""
AnalyticityTrait(::Type{<:NumericAIMProblem}) = IsNumeric()

"""
All problem types must implement a λ0 function.
This behaviour is enforced by the default implementations.
"""
λ0(x::T) where {T} = λ0(AnalyticityTrait(T), x)
λ0(::IsAnalytic, x) = error("Please implement a λ0 function for ", typeof(x))
λ0(::IsNumeric, x) = error("Please implement a λ0 function for ", typeof(x))

"""
All problem types must implement a S0 function.
This behaviour is enforced by the default implementations.
"""
S0(x::T) where {T} = S0(AnalyticityTrait(T), x)
S0(::IsAnalytic, x) = error("Please implement a S0 function for ", typeof(x))
S0(::IsNumeric, x) = error("Please implement a S0 function for ", typeof(x))

"""
All problem types must implement get_niter to return the number of iterations to perform.
"""
get_niter(x::T) where {T} = get_niter(AnalyticityTrait(T), x)
get_niter(::IsAnalytic, x) = error("Please implement a get_niter function for ", typeof(x))
get_niter(::IsNumeric, x) = error("Please implement a get_niter function for ", typeof(x))

"""
All problem types must implement get_x0 to return AIM's point of evaluation.
"""
get_x0(x::T) where {T} = get_x0(AnalyticityTrait(T), x)
get_x0(::IsAnalytic, x) = error("Please implement a get_x0 function for ", typeof(x))
get_x0(::IsNumeric, x) = error("Please implement a get_x0 function for ", typeof(x))

"""
Analytic problems must implement an acessor to the variable of the ODE.
"""
get_ODEvar(x::T) where {T} = get_ODEvar(AnalyticityTrait(T), x)
get_ODEvar(::IsAnalytic, x) = error("Please implement a get_ODEvar function for ", typeof(x))
get_ODEvar(::IsNumeric, x) = error("Numeric problem", typeof(x), "cannot implement get_ODEvar")

"""
Analytic problems must implement an acessor to the eigenvalue of the ODE.
"""
get_ODEeigen(x::T) where {T} = get_ODEeigen(AnalyticityTrait(T), x)
get_ODEeigen(::IsAnalytic, x) = error("Please implement a get_ODEeigen function for ", typeof(x))
get_ODEeigen(::IsNumeric, x) = error("Numeric problem", typeof(x), "cannot implement get_ODEeigen")

# ------------------------------------------------------------------
# 5. Includes
# ------------------------------------------------------------------

# --- Include core functionality ---
include("./core/aim_cache.jl")
include("./core/analytic_tools.jl")
include("./core/aim_step.jl")
include("./core/compute_delta.jl")

# --- Include eigenvalue computing functionality ---
include("./eigenvalue_computing/compute_eigenvalues.jl")
include("./eigenvalue_computing/eigenvalues_in_grid.jl")

end # QuasinormalModes
