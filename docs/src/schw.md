# Complete Example: Schwarzschild Quasinormal Modes

To illustrate how to use `QuasinormalModes.jl` we will show from start to finish how to compute the quasinormal modes of a Schwarzschild black hole perturbed by a field. This section will follow closely Emanuele Berti's lectures on black hole perturbation theory, which can be found [here](https://www.dropbox.com/sh/9th1um175m8gco9/AACCIkvNa3h-zdHBMZkQ31Baa?dl=0) and in Ref. [[3]](https://arxiv.org/abs/1111.5024)

# Mathematical preliminaries

Let's say that our Schwarzschild black hole is being perturbed by an external field ``\psi_{ls}`` where ``s`` is the spin of the field (``s = 0, 1, 2`` for scalar, electromagnetic and gravitational perturbations, respectively) and ``l`` is the angular index of the perturbation. Using mass units such that ``2M=1`` the "master" radial equation governing the perturbation is 

```math
r(r-1)\psi_{ls}^{\prime\prime}(r) + \psi_{ls}^{\prime}(r) - \left[ l(l+1) - \frac{s^2-1}{r} - \frac{\omega^2 r^3}{r-1} \right]\psi_{ls}(r) = 0
```
where primes denote derivatives with respect to the radial coordinate ``r`` and ``\omega`` are the quasinormal frequencies. Since we are solving for quasinormal modes, we need to enforce the proper boundary conditions in the master equation: Classically no wave can escape from the BH's event horizon and at spatial infinity waves can only "leave" the space-time. It's thus said that our field is purely *ingoing* in the event horizon (when ``r\rightarrow 1``) and purely *outgoing* at spatial infinity (when ``r\rightarrow\infty``). Mathematically, this means that the solution to the master equation must be of the form

```math
\psi_{ls}(r) = (r-1)^{-i \omega} r^{2 i \omega} e^{i \omega (r-1)}f_{ls}(r)
```

By substituting this solution *ansatz* in the master equation, we obtain a new 2nd order ODE, now for the function ``f_{ls}(r)``. This new ODE is enforcing the correct boundary quasinormal mode boundary conditions. This process usually referred to as incorporating the boundary conditions into the differential equation. The resulting equation reads

```math
r \left((r-1) r f^{\prime\prime}(r)+\left(1+2 i \left(r^2-2\right) \omega \right) f^\prime(r)\right)+f(r) \left(-r \left(l^2+l-4 \omega ^2\right)+s^2+(2 \omega +i)^2\right) = 0
```

The last step, although not strictly required, facilitates the numerical handling of the equation. Because the radial coordinates extends from the event horizon to infinity, that is, ``r\in [1,\infty]`` and computers can't handle infinities, we re-scale the ODE's domain to a finite one. This can be easily done with the change of variables

```math
x = 1 - \frac{1}{r}
```
which implies that when ``r=1`` we have ``x=0`` and when ``r\rightarrow\infty`` we have ``x = 1``. Thus the solution domain has been successfully compactifyied in the interval ``x\in[0,1]``. By making this change of variables we get to the final form of the master equation which we will actually feed to `QuasinormalModes.jl`

```math
-x (x-1)^2 f^{\prime\prime}(x) + (x (4 i (x-2) \omega -3 x+4)+2 i \omega -1) f^\prime(x)+f(x) \left(l^2+l+\left(s^2-1\right) (x-1)+4 (x-2) \omega ^2+4 i (x-1) \omega \right) = 0
```

# Implementing the master equation

`QuasinormalModes.jl` uses Julia's type system to implement structures that can used to solve the eigenvalue problem via the AIM. All ODEs are implemented as structures that are subtypes of one of the following abstract types. Each abstract type listed bellow is a subtype of the parent type `QNMData`:

1. `QuadFreqData` - Aimed at ODEs in which the eigenvalue is a quadratic polynomial.
2. `GenFreqData` - Aimed at ODEs in which the eigenvalue a general function.

With structures that are a subtype of `QuadFreqData`, `QuasinormalModes.jl` takes advantage of the polynomial nature of the eigenvalues in the ODE and operates in a semi-analytic way. Modes are computed by finding roots of the large polynomials in the eigenvalue produced by the iteration steps of the AIM using [PolynomialRoots.jl](https://github.com/giordano/PolynomialRoots.jl). Since all roots are found, this approach can generate extensive lists of quasinormal modes without the need of an initial "guess".

With structures that are a subtype of `GenFreqData`, `QuasinormalModes.jl` makes no assumptions about the eigenvalue and operates in a purely numeric way. In fact, the ODE might even contain numeric functions of the eigenvalue. Since it's hard to determine all the roots and poles of generic complex functions, when operating in this mode `QuasinormalModes.jl` requires in initial "guess" to be given for root finding which is carried out by [NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl).

!!! note "Semi-analytic VS numeric approach"
    Because of the semi-analytic nature of the operation performed when a structure is a subtype of `QuadFreqData`, `QuasinormalModes.jl` is naturally slower to compute modes in this case. One may also find that for a large number of iterations the AIM might even fail to find modes. A general good approach would be to use the semi-analytic mode to generate lists of eigenvalues for a number of iterations that runs reasonably fast and then use these results as initial guesses for the numeric mode with a high number of iterations. 

The AIM can be used to solve generic 2nd order ODES of the form

```math
y^{\prime\prime}(x) = \lambda_0(x)y^\prime(x) + s_0(x)y(x)
```

The structures we construct must contain the ``\lambda_0(x)`` and ``s_0(x)`` functions, which we will derive from or master equation. If we are sub-typing `QuadFreqData` we must implement these functions as well as ``x`` and ``\omega`` as symbolic variables and include [SymEngine.jl](https://github.com/symengine/SymEngine.jl) as a dependency. If we are sub-typing `GenFreqData` we must include [TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl) in the project as this package is responsible for computing the high order derivatives necessary to the AIM.

## Using `QuadFreqData`

Here is how to implement the Schwarzchild master equation we derived previously taking advantage of the fact that ``\omega`` is a quadratic polynomial

```julia
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
```

Any structure that subtypes `QuadFreqData` **must** contain the following fields:
1. `nIter::UInt32`
2. `x0::Float64`
3. `vars::Tuple{Basic, Basic}`
4. `exprs::Tuple{Basic, Basic}`

This is because `QuasinormalModes.jl` internally expects all concrete types that are sub-types of `QuadFreqData` to have these fields in order to function correctly. We recommend that users use the above example as a template to implement their own space-times.

## Using `GenFreqData`

Here is how to implement the Schwarzchild master equation we derived previously in a purely numeric way

```@julia
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
```

Any structure that subtypes `GenFreqData` **must** contain the following fields:
1. `nIter::UInt32`
2. `x0::Float64`

Additionally the operator `()` must be overloaded to work with the newly defined type and it's signature must be **exactly** that of the example. In the expressions for `λ0` and `S0` the space-time specific parameters `l` and `s` have been replaced by `data.l` and `data.s` so that they can be located during runtime. Again we recommend that users use the above example as a template to implement their own space-times.

!!! note "This is too much code!"
    At this point the user might be scared with the amount of boilerplate code required to implement a spacetime so that `QuasinormalMods.jl` can do it's job. The authors are well aware of this fact and intend to improve the package so that it eventually becomes easier to implement these structures with the correct signatures. In the meantime, we recommend that the user uses the provided `.jl` example files as starting points to their projects. In the following section we will show how to use these structures to compute QNMs. Rest assured that once the structures are constructed, computing the modes is much easier.

# Computing the modes

In order to compute the modes it is only necessary to invoke the function `computeQNMs` with the space-time structure and in the case of a structure that sub-types `GenFreqData` the initial guess as a second argument.

For example, if our structure is a syb-type of `QuadFreqData` we would do

```julia
data = SchwarzschildData(UInt32(50), UInt32(0), UInt32(0))
modes = computeQNMs(data, plr_epsilon=BigFloat("1.0e-20"))
```
and `modes` would be an object of type `Array{Complex{BigFloat},1}` containing the computed modes with 20 iterations of the AIM.

If our structure is a syb-type of `GenFreqData` we would do

```julia
data = NSchwarzschildData(UInt32(50), UInt32(0), UInt32(0))
mode = computeQNMs(data, Complex{BigFloat}( BigFloat("0.220"), BigFloat("-0.209")))
```
and `mode` would be an object of type `SolverResults` returned by `NLsolve`. We encourage the reader to read `NLsolve`'s documentation, but what users usually will want to know is if the solver converged to a mode and the actual value of the mode. Convergence can be tested with

```julia
converged = mode.x_converged || mode.f_converged
```

and the real and imaginary part of the solution are obtained with

```julia
real_part = mode.zero[1]
imag_part = mode.zero[2]
```

We chose to return `NLsolve`'s solution object instead of a more "neatly" formatted complex number to allow the users more freedom in their applications.

Aside from `computeQNMs`, if our space-time structure is a sub-type of `GenFreqData`, we can explore the existence of modes inside a grid of points in the complex plane. This is accomplished by the function `modesInGrid` In code we wold do

```julia
data = NSchwarzschildData(UInt32(50), UInt32(0), UInt32(0))

grid_start = Complex{BigFloat}(BigFloat("0.01"), BigFloat("-1.0"))
grid_end = Complex{BigFloat}(BigFloat("1.0"), BigFloat("-0.01"))

real_pts = 5
imag_pts = 5

grid = (grid_start, grid_end, real_pts, imag_pts)

modes = modesInGrid(data, grid)
```

The variable `modes` would be an array of type `Array{Complex{BigFloat},1}` containing the modes found in the grid. We specify a grid by start and end points in the complex plane and choosing how many points to compute in the real an imaginary part of the grid. In the example we will try to find modes for 25 (5x5) initial conditions inside the chosen rectangle in the complex plane. We might see that some modes repeat inside the `modes` array. this is of course to be expected as multiple initial conditions will lead `NLsolve` to converge to a certain mode.

# Printing the modes

Having computed a single mode or an array of modes we can easily operate on these results in any way we want. For convenience we provide functions to save print and filter modes. To save modes (an object of type `Array{Complex{BigFloat},1}` under the name `"qnms.dat"` we would simply do

```julia
saveQNMs(modes, filename = "qnms.dat")
```

There are also arguments controlling the cutoff of the real part of the modes (the value from which modes are no longer saved in the output file) and whether to save instabilities (modes with positive imaginary part).

The function `printQNMs` has the exact same signature as `saveQNMs` but instead of outputting to a file it writes the modes to `stdout`. Finally the function `filterQNMs!` has also the same signature and filters an array of QNMs in-place according to the `cutoff` and `instab` parameters.

!!! note "Complete code"
    The complete code of this example can be found under **TODOl: LINK TO EXAMPLE FILE**
