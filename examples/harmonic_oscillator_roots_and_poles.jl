using QuasinormalModes
using RootsAndPoles

# ------------------------------------------------------------------
# 1. Numeric harmonic oscillator
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
# 2. Constructing problems and caches
# ------------------------------------------------------------------

const m = Serial()
const p = NHarmonicOscilatorData(0x0000A, Complex(0.5, 0.0))
const c = AIMCache(p)

# ------------------------------------------------------------------
# 3. Creatting a wrapper function to pass to RootsAndPoles.jl
# ------------------------------------------------------------------

function δ_wrapper(z)
    computeDelta!(m, p, c, z)
end

# ------------------------------------------------------------------
# 3. RootsAndPoles.jl search domain and mesh construction
# ------------------------------------------------------------------

const xb = 0.0   # real part begin
const xe = 21.0  # real part end
const yb = 0.0   # imag part begin
const ye = 21.0  # imag part end
const r = 1.0    # initial mesh step

const origcoords = rectangulardomain(Complex(xb, yb), Complex(xe, ye), r)

# ------------------------------------------------------------------
# 4. RootsAndPoles.jl search settings
# ------------------------------------------------------------------

# For details, see https://github.com/fgasdia/RootsAndPoles.jl
params = GRPFParams(
    100,     # the maximum number of refinement iterations before `grpf` returns.
    50000,   # the maximum number of Delaunay tessalation nodes before `grpf` returns.
    3,       # maximum ratio of the longest to shortest side length of Delaunay triangles before they are split during `grpf` refinement iterations.
    5000,    # provide a size hint to the total number of expected nodes in the Delaunay tesselation. Setting this number approximately correct can improve performance
    1.0e-12, # maximum allowed edge length of the tesselation defined in the `origcoords` domain before returning
    false    # use `Threads.@threads` to run the user-provided function `fcn` across the `DelaunayTriangulation`
)

# ------------------------------------------------------------------
# 5. Root finding and printing
# ------------------------------------------------------------------

roots, poles = grpf(δ_wrapper, origcoords, params)
sort!(roots, by = z -> real(z))
sort!(poles, by = z -> real(z))

println("Roots:")

for root in roots
    println(root)
end

println("-------------------------------------------")
println("Poles:")

for pole in poles
    println(pole)
end

