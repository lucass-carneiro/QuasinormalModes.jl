# Complete Example: The Harmonic Oscillator

We will now turn away from general relativity and use `QuasinormalModes.jl` to compute to compute the energy eigenvalues of the quantum harmonic oscillator following [this paper](https://arxiv.org/abs/1111.5024).

# Mathematical preliminaries

If we measure the energy of the system in units of ``\hbar\omega`` and distance in units of ``\sqrt{\hbar/(m\omega)}`` the time independent Schrödinger equation for the quantum harmonic oscillator is written as

```math
-\psi^{\prime\prime}(x) + x^2\psi(x) = \epsilon\psi(x),
```
where we defined ``\epsilon \equiv 2 E`` and ``E`` is the quantum state's energy. Imposing that ``\psi(x)`` decays like a Gaussian distribution asymptotically, we apply the ansatz

```math
\psi(x) = e^{-x^2/2}f(x)
```
which, substituting in the original equation, yields

```math
f^{\prime\prime}(x) = 2 x f^\prime(x) + (1-\epsilon)f(x).
```

This allows us to easily identify ``\lambda_0 = 2x`` and ``s_0 = 1 - \epsilon``. In all our implementations we shall refer the sought eigenvalue ``\epsilon`` using the variable `ω` in order to maintain consistency with the previous example.

# Implementing the master equation as an analytic problem

The first step is to load the required packages to run this example: `QuasinormalModes` and `SymEngine`:

```julia
using QuasinormalModes
using SymEngine
```

Next, we create a parametric type that sub-types `AnalyticAIMProblem`. As the eigenvalue in the master equation is a quadratic polynomial, we will sub-type `QuadraticEigenvalueProblem` with the following structure:

```julia
struct HarmonicOscilatorData{N,T} <: QuadraticEigenvalueProblem{N,T}
    nIter::N
    x0::T

    vars::Tuple{Basic, Basic}
    exprs::Tuple{Basic, Basic}
end
```
Now we implement the constructor and extend the default implementations:

```julia
function HarmonicOscilatorData(nIter::N, x0::T) where {N,T}
	
    vars = @vars x ω

    λ0 = 2*x
    S0 = 1 - ω

    return HarmonicOscilatorData{N,T}(nIter, x0, vars, (λ0, S0))
end

QuasinormalModes.λ0(d::HarmonicOscilatorData{N,T}) where {N,T} = d.exprs[1]
QuasinormalModes.S0(d::HarmonicOscilatorData{N,T}) where {N,T}  = d.exprs[2]

QuasinormalModes.get_niter(d::HarmonicOscilatorData{N,T}) where {N,T} = d.nIter
QuasinormalModes.get_x0(d::HarmonicOscilatorData{N,T}) where {N,T} = d.x0

QuasinormalModes.get_ODEvar(d::HarmonicOscilatorData{N,T}) where {N,T} = d.vars[1]
QuasinormalModes.get_ODEeigen(d::HarmonicOscilatorData{N,T}) where {N,T} = d.vars[2]
```

# Implementing the master equation as a numeric problem

The structure, constructor and extensions are

```julia
struct NHarmonicOscilatorData{N,T} <: NumericAIMProblem{N,T}
    nIter::N
    x0::T
end

function NHarmonicOscilatorData(nIter::N, x0::T) where {N,T}
    return NHarmonicOscilatorData{N,T}(nIter, x0)
end

QuasinormalModes.λ0(::NHarmonicOscilatorData{N,T}) where {N,T} = (x,ω) -> 2*x
QuasinormalModes.S0(::NHarmonicOscilatorData{N,T}) where {N,T} = (x,ω) -> 1 - ω + x - x

QuasinormalModes.get_niter(d::NHarmonicOscilatorData{N,T}) where {N,T} = d.nIter
QuasinormalModes.get_x0(d::NHarmonicOscilatorData{N,T}) where {N,T} = d.x0
```

# Constructing problems and initializing the cache

Once again, we create our problems and cache objects by calling the constructors:

```julia
p_ana = HarmonicOscilatorData(0x0000A, 0.5);
p_num = NHarmonicOscilatorData(0x0000A, 0.5);

c_ana = AIMCache(p_ana)
c_num = AIMCache(p_num)
```

Here we are setting up problems to be solved using 10 iterations with `x0 = 0.5`

# Computing the eigenvalues

Once again we compute the eigenvalues by calling

```julia
ev_ana = computeEigenvalues(Serial(), p_ana, c_ana)
ev_num = eigenvaluesInGrid(Serial(), p_num, c_num, (0.0, 21.0))
```

The results are two arrays, containing the eigenvalues. As before, we define a function to print the results to `stdout`

```julia
function printEigen(eigenvalues)
    println("--------------------------------------")

    for i in eachindex(eigenvalues)
        println("n = $i, ω = $(eigenvalues[i])")
    end
    
    println("--------------------------------------")

    return nothing
end

println("Analytic results")
printEigen(reverse!(ev_ana))

println("Numeric results")
printEigen(ev_num)
```

The complete source file for this example can be found in [harmonic_oscillator.jl](https://github.com/lucass-carneiro/QuasinormalModes.jl/blob/master/examples/harmonic_oscillator.jl). The output is agreement with the expected result for the eigenenergies of the harmonic oscillator, that is, ``E_n = n + 1/2``

```
Analytic results
--------------------------------------
n = 1, ω = 0.9999999999999999 + 0.0im
n = 2, ω = 2.9999999999999964 + 0.0im
n = 3, ω = 4.999999999999426 + 0.0im
n = 4, ω = 7.000000000006788 + 0.0im
n = 5, ω = 8.999999999980533 + 0.0im
n = 6, ω = 10.999999804542819 + 0.0im
n = 7, ω = 13.000000959453153 + 0.0im
n = 8, ω = 14.999998108295404 + 0.0im
n = 9, ω = 17.00000187312756 + 0.0im
n = 10, ω = 18.999999068409203 - 0.0im
n = 11, ω = 21.000000186185098 + 0.0im
--------------------------------------
Numeric results
--------------------------------------
n = 1, ω = 1.0
n = 2, ω = 3.000000000000006
n = 3, ω = 5.000000000000002
n = 4, ω = 7.000000000000006
n = 5, ω = 8.999999999999988
n = 6, ω = 11.0
n = 7, ω = 12.999999999999977
n = 8, ω = 15.000000000000014
n = 9, ω = 16.999999999999908
n = 10, ω = 19.000000000000025
n = 11, ω = 21.0
--------------------------------------
```
