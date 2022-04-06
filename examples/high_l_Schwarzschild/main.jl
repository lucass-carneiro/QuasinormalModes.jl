# ------------------------------------------------------------------
# Compute Schwarzschild QNMs perturbed by a scalar field
# ------------------------------------------------------------------

using QuasinormalModes
using ProgressMeter
using Dates

include("Schwarzschild_boson.jl")
include("guesses.jl")
include("compute_boson.jl")

# Creates a list of Schwarzschild modes with spin 0
compute_boson(convert(UInt32, 100), 0)
