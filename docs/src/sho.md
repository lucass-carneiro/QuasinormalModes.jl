# Complete Example: The Harmonic Oscilator

We can use `QuasinormalModes.jl` to compute the eigenvalues of any 2nd order differential equation. To illustrate this we will show how to compute the energy eigenvalues of the quantum harmonic oscillator following Ref. [[3]](https://arxiv.org/abs/1111.5024).

# Mathematical preliminaries

If we measure the energy of the system in units of ``\hbar\omega`` and distance in units of ``\sqrt{\hbar/(m\omega)}`` the time independent Schrödinger equation for the system is written as

```math
-\psi^{\prime\prime}(x) + x^2\psi(x) = \epsilon\psi(x)
```
where we defined ``\epsilon \equiv 2 E`` and ``E`` is the quantum state's energy. Imposing that ``\psi(x)`` decays like a Gaussian distribution asymptotically, we apply the ansatz

```math
\psi(x) = e^{-x^2/2}f(x)
```
which substituting in the original equation yeilds

```math
f^{\prime\prime}(x) = 2 x f^\prime(x) + (1-\epsilon)f(x)
```

This allows us to easily identify ``\lambda_0 = 2x`` and ``s_0 = 1 - \epsilon``. In all our implementations we shall refer the sough eigenvalue ``\epsilon`` using the variable `ω` in order to maintain consistency with the previous example.

# Implementing Using `QuadFreqData`

```julia
struct HarmonicOscilatorData <: QuadFreqData # Indicate that we will operate semi-analytically
    nIter::UInt32                            # The number of iteration the AIM will perform
    x0::Float64                              # The evaluation point for the AIM functions

    vars::Tuple{Basic, Basic}                # The variables x and ω as SymEngine expressions.
    exprs::Tuple{Basic, Basic}               # The functions λ0 and S0 as SymEngine expressions.
    
    function HarmonicOscilatorData(nIter::UInt32, x0::Float64 = 0.5)
	
        vars = @vars x ω
	
        λ0 = 2*x
        S0 = 1 - ω

        return new(nIter, x0, vars, (λ0, S0))
    end
end 
```

# Implementing Using `GenFreqData`

```julia
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
```

!!! note "Add and subtract `x`"
    In the expression for `S0` we have added and subtracted the variable `x`. This is done because as `S0` is independent of `x` the variable `ts` would be of type `Complex{BigFloat}` and the `ts.coeffs` instruction would fail. By adding and subtracting s we transform the expression for `S0` in a `TaylorSeries` object that can be expanded correctly.

# Computing the modes

To compute the modes we proceed exactly as in the previous Schwarzschild example. The eigenenergies of the harmonic oscillator are well know

```math
\epsilon_n = 2 n + 1
```
thus we will use custom printing routines to display only the real part of the computed energies and sort them by ascending order. In the semi-analytic case this can be done with

```
data = HarmonicOscilatorData(UInt32(50))
modes = computeQNMs(data, plr_epsilon=BigFloat("1.0e-20"))

sort!(modes, by = x -> real(x))
    
for mode in modes
    println(round(Int64,real(mode)))
end	
```

which will produce a list of the 50 first even numbers which corresponds to the expected result. We can also verify that the purely numeric case also produces correct results by computing a single energy with

```julia
data = NHarmonicOscilatorData(UInt32(50))
mode = computeQNMs(data, Complex{BigFloat}(BigFloat("15.0"), BigFloat("1.0")), nls_xtol=BigFloat("1.0e-20"), nls_ftol=BigFloat("1.0e-20"))

if mode.x_converged || mode.f_converged
    println(mode.zero[1], "    ", mode.zero[2])
else
    println("Did not converge to any modes :-(")
end
```
or by sweeping the complex grid with

```julia
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
```
