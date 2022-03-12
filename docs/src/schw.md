# Complete Example: Schwarzschild Quasinormal Modes

To illustrate how to use `QuasinormalModes.jl` we will show from start to finish how to compute the quasinormal modes of a Schwarzschild black hole perturbed by an external field. This section will follow closely Emanuele Berti's lectures on black hole perturbation theory, which can be found [here](https://www.dropbox.com/sh/9th1um175m8gco9/AACCIkvNa3h-zdHBMZkQ31Baa?dl=0) and in Ref. [[3]](https://arxiv.org/abs/1111.5024)

# Mathematical preliminaries

Let's say that our Schwarzschild black hole is being perturbed by an external field ``\psi_{ls}`` where ``s`` is the spin of the field (``s = 0, 1, 2`` for scalar, electromagnetic and gravitational perturbations, respectively) and ``l`` is the angular index of the perturbation. Using mass units such that ``2M=1`` the "master" radial equation governing the perturbation is 

```math
r(r-1)\psi_{ls}^{\prime\prime}(r) + \psi_{ls}^{\prime}(r) - \left[ l(l+1) - \frac{s^2-1}{r} - \frac{\omega^2 r^3}{r-1} \right]\psi_{ls}(r) = 0
```
where primes denote derivatives with respect to the radial coordinate ``r`` and ``\omega`` are the quasinormal frequencies. Since we are solving for quasinormal modes, we need to enforce the proper boundary conditions in the master equation: classically no wave can escape from the BH's event horizon and at spatial infinity waves can only "leave" the space-time. It's thus said that our field is purely *ingoing* in the event horizon (when ``r\rightarrow 1``) and purely *outgoing* at spatial infinity (when ``r\rightarrow\infty``). Mathematically, this means that the solution to the master equation must be of the form

```math
\psi_{ls}(r) = (r-1)^{-i \omega} r^{2 i \omega} e^{i \omega (r-1)}f_{ls}(r).
```

By substituting this solution *ansatz* in the master equation, we obtain a new 2nd order ODE, now for the function ``f_{ls}(r)``. This new ODE is enforcing the correct quasinormal mode boundary conditions. This process usually referred to as incorporating the boundary conditions into the differential equation. The resulting equation reads

```math
r \left((r-1) r f^{\prime\prime}(r)+\left(1+2 i \left(r^2-2\right) \omega \right) f^\prime(r)\right)+f(r) \left(-r \left(l^2+l-4 \omega ^2\right)+s^2+(2 \omega +i)^2\right) = 0.
```

The last step, although not strictly required, facilitates the numerical handling of the equation. Because the radial coordinates extends from the event horizon to infinity, that is, ``r\in [1,\infty]`` and computers can't handle infinities, we re-scale the ODE's domain to a finite one. This can be easily done with the change of variables

```math
x = 1 - \frac{1}{r}
```
which implies that when ``r=1`` we have ``x=0`` and when ``r\rightarrow\infty`` we have ``x = 1``. Thus the solution domain has been successfully compactifyied in the interval ``x\in[0,1]``. By making this change of variables we get to the final form of the master equation which we will actually feed to `QuasinormalModes.jl`

```math
-x (x-1)^2 f^{\prime\prime}(x) + (x (4 i (x-2) \omega -3 x+4)+2 i \omega -1) f^\prime(x)+f(x) \left(l^2+l+\left(s^2-1\right) (x-1)+4 (x-2) \omega ^2+4 i (x-1) \omega \right) = 0.
```

# Implementing the master equation as an analytic problem

The first step is to load the required packages to run this example: `QuasinormalModes` and `SymEngine`:

```julia
using QuasinormalModes
using SymEngine
```

Next, we create a parametric type that sub-types `AnalyticAIMProblem`. As the eigenvalue in the master equation is a quadratic polynomial, we will sub-type `QuadraticEigenvalueProblem` with the following structure:

```julia
struct SchwarzschildData{N,T} <: QuadraticEigenvalueProblem{N,T}
    nIter::N
    x0::T

    vars::Tuple{Basic, Basic}
    exprs::Tuple{Basic, Basic}
end
```

As the reader might notice the structure is quite simple. The variables `nIter` and `x0` store the AIM's number of iterations and expansion point, respectively while `vars` will be responsible for storing the `SymEngine` variables representing the ODE's variable and eigenvalue, respectively, as a tuple. Finally `exprs` will store the `SymEngine` expressions for the `λ0` and `S0` parts of the ODE.

Next we create a parametric constructor for `SchwarzschildData` that will initializes the fields:

```julia
function SchwarzschildData(nIter::N, x0::T, l::N, s::N) where {N,T}
    vars = @vars x ω

    λ0 = (-1 + (2*im)*ω + x*(4 - 3*x + (4*im)*(-2 + x)*ω))/((-1 + x)^2*x)
    S0 = (l + l^2 + (-1 + s^2)*(-1 + x) + (4*im)*(-1 + x)*ω + 4*(-2 + x)*ω^2)/((-1 + x)^2*x)

    return SchwarzschildData{N,T}(nIter, x0, vars, (λ0, S0))
end
```

This constructor can be used by passing the values directly instead of explicitly declaring type parameters. The final step is to extend the default accessors functions to operate on `SchwarzschildData`:

```julia
QuasinormalModes.λ0(d::SchwarzschildData{N,T}) where {N,T} = d.exprs[1]
QuasinormalModes.S0(d::SchwarzschildData{N,T}) where {N,T}  = d.exprs[2]

QuasinormalModes.get_niter(d::SchwarzschildData{N,T}) where {N,T} = d.nIter
QuasinormalModes.get_x0(d::SchwarzschildData{N,T}) where {N,T} = d.x0

QuasinormalModes.get_ODEvar(d::SchwarzschildData{N,T}) where {N,T} = d.vars[1]
QuasinormalModes.get_ODEeigen(d::SchwarzschildData{N,T}) where {N,T} = d.vars[2]
```
These functions are fairly straightforward accessors and require no additional comment.

# Implementing the master equation as a numeric problem

Again we start by defining a structure but this time around we sub-type `NumericAIMProblem`

```julia
struct NSchwarzschildData{N,T} <: NumericAIMProblem{N,T}
    nIter::N
    x0::T
    l::N
    s::N
end
```

Here `nIter` and `x0` have the same meaning as before, but now instead of storing symbolic variables and expressions we store two additional unsigned integers, `l` and `s`. These are the angular and spin parameters of the master equation. Here we must store them in the struct as they can't be "embedded" into the expressions for `λ0` and `S0` as in the analytic case.

We proceed once again by creating a more convenient constructor. This time no intermediate computation is required upon construction:

```julia
function NSchwarzschildData(nIter::N, x0::T, l::N, s::N) where {N,T}
    return NSchwarzschildData{N,T}(nIter, x0, l, s)
end
```

Finally we extend the default implementations

```julia
QuasinormalModes.λ0(::NSchwarzschildData{N,T}) where {N,T} = (x,ω) -> (-1 + (2*im)*ω + x*(4 - 3*x + (4*im)*(-2 + x)*ω))/((-1 + x)^2*x)
QuasinormalModes.S0(d::NSchwarzschildData{N,T}) where {N,T} = (x,ω) -> (d.l + d.l^2 + (-1 + d.s^2)*(-1 + x) + (4*im)*(-1 + x)*ω + 4*(-2 + x)*ω^2)/((-1 + x)^2*x)

QuasinormalModes.get_niter(d::NSchwarzschildData{N,T}) where {N,T} = d.nIter
QuasinormalModes.get_x0(d::NSchwarzschildData{N,T}) where {N,T} = d.x0
```

This time `λ0` and `S0` return two parameters lambda functions that will be called multiple times during the evaluation of the AIM. As we've previously mentioned, the first parameters is assumed to be the ODE's variables while the second the ODE's eigenvalue. The body of each lambda is the expression for their respective parts on the ODE.

# Constructing problems and initializing the cache

We create our problems and cache objects by calling the constructors:

```julia
p_ana = SchwarzschildData(0x00030, Complex(BigFloat("0.43"), BigFloat("0.0")), 0x00000, 0x00000);
p_num = NSchwarzschildData(0x00030, Complex(0.43, 0.0), 0x00000, 0x00000);

c_ana = AIMCache(p_ana)
c_num = AIMCache(p_num)
```

Here we are setting up problems to be solved using 48 iterations with `x0 = 0.43 + 0.0*im` and `l = s = 0`.

# Computing the eigenvalues

Finally, to compute the quasinormal frequencies we will call `computeEigenvalues(Serial(), p_ana, c_ana)` (or `computeEigenvalues(Threaded(), p_ana, c_ana)` if you wish, but don't forget to start julia with the `--threads` option). This returns an array with all the roots of the quantization condition. We will sort the array by descending order in the imaginary part and after that we will filter the array to remove entries whose real part is too small or with a positive imaginary part and print the result to `stdout`:

```julia
m_ana = computeEigenvalues(Serial(), p_ana, c_ana)

function printQNMs(qnms, cutoff, instab)
    println("-"^165)
    println("|", " "^36, "Re(omega)", " "^36, " ", " "^36, "Im(omega)", " "^36, "|")
    println("-"^165)

    for qnm in qnms
        if real(qnm) > cutoff && ( instab ? true : imag(qnm) < big"0.0" )
        println(real(qnm), "    ", imag(qnm))
        end
    end
    
    println("-"^165)

    return nothing
end

sort!(m_ana, by = x -> imag(x))
printQNMs(m_ana, 1.0e-10, false)
```

Remember that (as was discussed in [here](org.md#Computing-eigenvalues-and-general-workflow-guidelines)) not all values are actually eigenvalues of the ODE (that is, quasinormal modes). Next we will call

```julia
ev = computeEigenvalues(Serial(), p_num, c_num, Complex(0.22, -0.20), nls_xtol = 1.0e-10, nls_ftol = 1.0e-10)
```

The variable `ev` now contains a `SolverResults` object from the [NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl) package. The first solution element represents the real part of the computed mode while the second represents the imaginary part. The object also contains information about the convergence of the method. Note that with a numerical problem we can only find one mode at a time using a certain initial guess. This can be somewhat remedied by using `eigenvaluesInGrid`, which uses multiple initial conditions as a guess and collects the converged results. The complete source code of this example can be found [here](https://github.com/lucass-carneiro/QuasinormalModes.jl/blob/master/examples/schwarzschild.jl).

# Interpreting `SolverResults` output

If we print `ev` to `stdout`, we will see something like

```
Results of Nonlinear Solver Algorithm
 * Algorithm: Trust-region with dogleg and autoscaling
 * Starting Point: [0.22, -0.2]
 * Zero: [0.22090807949439797, -0.20979076729905097]
 * Inf-norm of residuals: (a large number)
 * Iterations: 5
 * Convergence: true
   * |x - x'| < 1.0e-10: true
   * |f(x)| < 1.0e-10: false
 * Function Calls (f): 6
 * Jacobian Calls (df/dx): 6
```

This is precisely the information returned by the `SolverResults` object in human readable format. Most of these fields are self explanatory, but we must pay close attention to the `Convergence` field where we see that for convergence to be declared, either of two tests must pass:

1. The `|x - x'| < 1.0e-10: true` test indicates that the difference between two solution candidates obtained by two consecutive iterations of the trust region algorithm are differing by an amount smaller than `1.0e-10`
2. The `|f(x)| < 1.0e-10: false` test indicates that the Inf-norm of the residuals (the value of the function we are trying to find the root for yields after we substitute a solution candidate back into it) is not smaller than `1.0e-10`. In fact, in this example it is a large number.

The fact that the first test passes and the second one does not, suggests that the root finding method is converging, since at each iteration the difference between candidate solutions gets smaller, but not to a true root of the AIM quantization condition (and therefore a quasinormal mode) since that we get a number far from zero when we substitute this value back at the quantization condition. To resolve this kind of ambiguity, we can employ a different root finding method that does not "polish" roots like Newton's method but "brackets" them, like the bisection method, that is guaranteed to converge to a root if it exists within an interval. Such a method for finding all roots and poles of a complex function within an interval exists and is know as the [GRPF](https://github.com/PioKow/GRPF) method, implemented in Julia by the [RootsAndPoles.jl](https://github.com/fgasdia/RootsAndPoles.jl) package. This algorithm works by subdividing a complex plane domain in triangles and applying a discretized version of Cauchy's argument principle, which counts the roots and poles of a complex function within a region. Despite relying internally on `NLsolve.jl` in the `computeEigenvalues` or `eigenvaluesInGrid` family of functions, `QuasinormalModes.jl`'s core responsibility is to compute the AIM quantization condition at any point in the complex plane for a given problem. Such computation is provided by the `computeDelta!` family of functions. By exposing this functionality, the user has complete freedom to choose what root finding method will be employed. In fact, internally, `computeEigenvalues` simply wraps `computeDelta!` into another function that is on the format accepted by `NLsolve.jl`. To see a concrete example of finding Schwarzschild modes with `RootsAndPoles.jl` look at [this](https://github.com/lucass-carneiro/QuasinormalModes.jl/blob/master/examples/schwarzschild_roots_and_poles.jl) example, which when asked to search for roots in the ``[0.1 - 1.0 i, 1.0 + 1 .0 i]`` range reports

```
Roots:
0.2082018398054518 - 0.7014133118267079im
0.2209861849099763 - 0.2095673605939246im
0.2082018397886055 - 0.7014133118222261im
-------------------------------------------
Poles:
0.2082018397894562 - 0.7014133118211809mi
```

Within this list, we can find the fundamental and first excited mode. By comparing these two results we can be sure that a found root is in fact a root of the quantization condition. Weather or not this root is an actual eigenvalue of the original second order ODE is a different story and was [previously](https://lucass-carneiro.github.io/QuasinormalModes.jl/dev/org/#Computing-eigenvalues-and-general-workflow-guidelines) discussed.