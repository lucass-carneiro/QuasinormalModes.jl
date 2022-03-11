using QuasinormalModes
using RootsAndPoles

# ------------------------------------------------------------------
# 1. Numeric Schwarzschild Black Hole
# ------------------------------------------------------------------

struct NSchwarzschildData{N,T} <: NumericAIMProblem{N,T}
    nIter::N
    x0::T
    l::N
    s::N
end

function NSchwarzschildData(nIter::N, x0::T, l::N, s::N) where {N,T}
    return NSchwarzschildData{N,T}(nIter, x0, l, s)
end

QuasinormalModes.λ0(::NSchwarzschildData{N,T}) where {N,T} = (x,ω) -> (-1 + (2*im)*ω + x*(4 - 3*x + (4*im)*(-2 + x)*ω))/((-1 + x)^2*x)
QuasinormalModes.S0(d::NSchwarzschildData{N,T}) where {N,T} = (x,ω) -> (d.l + d.l^2 + (-1 + d.s^2)*(-1 + x) + (4*im)*(-1 + x)*ω + 4*(-2 + x)*ω^2)/((-1 + x)^2*x)

QuasinormalModes.get_niter(d::NSchwarzschildData{N,T}) where {N,T} = d.nIter
QuasinormalModes.get_x0(d::NSchwarzschildData{N,T}) where {N,T} = d.x0

# ------------------------------------------------------------------
# 2. Constructing problems and caches
# ------------------------------------------------------------------

const m = Serial()
const p = NSchwarzschildData(0x00030, Complex(0.5, 0.0), 0x00000, 0x00000);
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

const xb = 0.1   # real part begin
const xe = 1.0   # real part end
const yb = -1.0  # imag part begin
const ye = 0.1   # imag part end
const r = 6.0e-3 # initial mesh step

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

println("Roots:")

for root in roots
    println(root)
end

println("-------------------------------------------")
println("Poles:")

for pole in poles
    println(pole)
end

