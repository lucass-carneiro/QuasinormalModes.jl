# Package organization

## A Brief description of the AIM

`QuasinormalModes.jl` is in its core an implementation of the Asymptotic Iteration Method. For a complete description of the general method the reader is encouraged to read [this paper](https://arxiv.org/abs/math-ph/0309066v1). Our implementation is based on the variation of the method described in [this paper](https://arxiv.org/abs/1111.5024). The method requires 3 basic steps:

1. Incorporate the asymptotic boundary conditions into the ODE.
2. Compactify the domain of the problem (if it isn't already compact).
3. Write the ODE in the form ``y^{\prime\prime}(x) = \lambda_0(x)y^{\prime}(x) + s_0(x)y(x)``

From the ODE coefficients ``\lambda_0(x)`` and ``s_0(x)`` the AIM computes the eigenvalues by requiring that the "quantization condition"

```math
\delta_n = s_n\lambda_{n-1} - s_{n-1}\lambda_{n} = 0
```

is satisfied, where

```math
\lambda_n = \lambda^\prime_{n-1} + s_{n-1} + \lambda_0 \lambda_{n-1}
```

and

```math
s_n = s^\prime_{n-1} + s_0 s_{n-1}.
```

`QuasinormalModes.jl` expects as input the ``\lambda_0(x)`` and ``s_0(x)`` coefficients computed with the 3 steps described above. This documentation contains two practical examples of how to obtain and feed such coefficients to the package.

## The type hierarchy

`QuasinormalModes.jl` employs two main strategies in order to find eigenvalues using the AIM: problems can be solved in a semi-analytic or purely numeric fashion. We make use of Julia's type system in order to implement structures that reflect these operation modes. All of the package's exported functionality is designed to operate on sub-types of abstract types that reflect the desired solution strategy (semi-analytic or numeric). The user is responsible for constructing concrete types that are sub-types of the exported abstract types with the actual problem specific information. It's thus useful to start by inspecting the package's exported type hierarchy:

```@raw html
<table border="0"><tr>
<td>
	<figure>
		<img src='../assets/types.svg' alt='missing'><br>
		<figcaption><em>QuasinormalModes.jl type hierarchy</em></figcaption>
	</figure>
</td>
</tr></table>
```

1. `AIMProblem` is the parent type of all problems that can be solved with this package. All problems must be sub-type it and a user can use it to construct functions that operate on all AIM solvable problems.
2. `NumericAIMProblem` is the parent type of all problems that can be solved using a numeric approach.
3. `AnalyticAIMProlem` is the parent type of all problems that can be solved using a semi-analytic approach.
4. `QuadraticEigenvaluePoblem` is a specific type of analytic problem whose eigenvalues appear in the ODE as a (possibly incomplete) quadratic polynomial.

All types are parameterized by two parameters: `N <: Unsigned` and `T <: Number` which represent respectively, the type used to represent the number of iterations the AIM will perform and the type used in the numeric computations of the method.

## Type traits

Type traits are non-exported abstract types that help the user to ensure that their sub-types implement the correct functions. Currently there is only one defined trait, called `AnalyticityTrait`. This trait can have two possible "values": `IsAnalytic` and `IsNumeric`, that are represented by concrete types. The default trait of an `AIMProblem` is `IsNumeric`, while any sub-type of `AnalyticAIMProblem` has the `IsAnalytic` and `NumericAIMProblem` have the `IsNumeric` trait

With these traits, we enforce that the user must implement for all problem types, the following functions:

1. `λ0`: Return the λ0 component of the ODE. The actual implementation depends heavily on the problem type.
2. `S0`: Return the S0 component of the ODE. The actual implementation depends heavily on the problem type.
3. `get_niter`: Return the number of iterations that the AIM will perform.
4. `get_x0`: Return the expansion point of the AIM.

For problems with the `IsAnalytic` trait, the user must implement the following functions function:
1. `get_ODEvar` which returns an object that represents the ODE's variable.
2. `get_ODEeigen` which returns an object that represents the ODE's eigenvalue.

Failure to implement these functions returns an error with the appropriate message. Note that these traits only check that such functions are implemented for a certain problem type and not that they follow a particular implementation pattern. The contract on the functions implementations is *soft* and will be clarified further on. Failure to abide by these soft contracts results in undefined behavior.

## Extending the default functionality

The following assumes that the package `SymEngine` is installed. If a problem type `P{N,T}` is a sub-type of `AnalyticAIMProblem{N,T}`, the user must extend the default implementations abiding by the following rules
1. `QuasinormalModes.λ0(p::P{N,T}) where {N,T}` must return a `SymEngine.Basic` object representing the symbolic expression for the `λ0` part of the ODE.
2. `QuasinormalModes.S0(p::P{N,T}) where {N,T}` must return a `SymEngine.Basic` object representing the symbolic expression for the `S0` part of the ODE.
3. `QuasinormalModes.get_ODEvar(p::P{N,T}) where {N,T}` must return a `SymEngine.Basic` objects representing the `SymEngine` variable associated with the ODE's variable.
4. `QuasinormalModes.get_ODEeigen(p::P{N,T}) where {N,T}` must return a `SymEngine.Basic` objects representing the `SymEngine` variable associated with the ODE's eigenvalue.

If a problem type `P{N,T}` is a sub-type of `NumericAIMProblem{N,T}`, the user must extend the default implementations abiding by the following rules
1. `QuasinormalModes.λ0(p::P{N,T}) where {N,T}` must return a lambda function of two parameters, the first representing the ODE's variable and the second representing the ODE's eigenvalue where the body represents the expression for the `λ0` part of the ODE.
2. `QuasinormalModes.S0(p::P{N,T}) where {N,T}` must return a lambda function of two parameters, the first representing the ODE's variable and the second representing the ODE's eigenvalue where the body represents the expression for the `S0` part of the ODE.

All problems `P{N,T}` that are a sub-type of `AIMProblem{N,T}` must extend the default implementations abiding by the following rules
1. `QuasinormalModes.get_niter(p::P{N,T}) where {N,T}` must return an unsigned number of type `N` representing the number of iterations for the AIM to perform.
2. `QuasinormalModes.get_x0(p::P{N,T}) where {N,T}` must return a number of type `T` representing the evaluation point of the AIM.

In the following sections, concrete examples of problems will be illustrated in order to better acquaint the user with the package and hopefully clear out any remaining misunderstandings.

!!! note "Semi-analytic VS numeric approach"
    Because of the semi-analytic nature of the operation performed when a structure is a subtype of `AnalyticAIMProblem`, `QuasinormalModes.jl` is naturally slower to compute modes in this case. One may also find that for a large number of iterations the AIM might even fail to find modes. A general good approach would be to use the semi-analytic mode to generate lists of eigenvalues for a number of iterations that runs reasonably fast and then use these results as initial guesses for the numeric mode with a high number of iterations.

## The memory cache

In order to minimize memory allocations, all functions that actually compute eigenvalues require a `AIMCache` object. Given a certain a problem `P{N,T}` it initializes memory for 8 arrays of size `get_niter(p) + one(N)` elements of type `T`. These arrays are used to store intermediate and final computation results. By using a cache object, we guarantee that memory for the computation data is allocated only once and not at each step of the AIM.

## Stepping methods

In order to compute ``\delta_n``, `QuasinormalModes.jl` evolves ``\lambda_0`` and ``s_0`` according to the previously stated equations. The evolution from step 1 to step ``n`` must happen sequentially but the step itself, that is, the computation of new values of ``\lambda`` and ``s`` from old ones can be performed in parallel. We've provided singleton types that allow the user to control this behavior by passing instances of those types to the eigenvalue computing functions. The user can currently choose the following stepping methods:

1. `Serial`: Each instruction in a single AIM step is executed sequentially.
2. `Threaded`: Instruction in a single AIM step is executed in parallel using Julia's built-in `Threads.@threads` macro.

## Computing eigenvalues and general workflow guidelines 

To compute eigenvalues, 3 functions are provided:
1. `computeDelta!`: Computes the AIM "quantization condition" ``\delta_n``.
2. `computeEigenvalues`: Computes a single, or a list of eigenvalues.
3. `eigenvaluesInGrid`: Find all eigenvalues in a certain numerical grid.

Depending on the problem type, these functions return and behave differently. In a `QuadraticEigenvalueProblem`.
1. `computeDelta!`: Returns a polynomial whose roots are the eigenvalues of the ODE.
2. `computeEigenvalues`: Computes the complete list of eigenvalues given by the roots of the computed polynomial.

In a `NumericAIMProblem`,
1. `computeDelta!`: Returns a value of the quantization condition at a given point in the complex plane.
2. `computeEigenvalues`: Computes a single eigenvalues from an initial trial frequency.
3. `eigenvaluesInGrid`: Attempts to find eigenvalues using a grid of real or complex data points as initial trial frequencies passed to `NLSolve`.

For more detail on these functions and their behaviors with each problem type refer to the the [API Reference](api_ref.md) where specific descriptions can be found. 

The AIM provides the user with "two degrees of freedom" when computing eigenvalues: The number of iterations to perform (which we refer by ``n``) and the point around which the ODE functions will be expanded (which we refer by``x_0``). Additionally, our implementation asks for an initial guess in `NumericAIMProblem`s to find the roots of ``\delta_n``, adding yet another degree of freedom to the method. So far, the literature around the AIM cannot provide a general prescription for choosing optimal values for ``n`` or ``x_0``, however, it is know that ``x_0`` can affect the speed at which the method converges to a correct solution and if ``n`` is chosen to be too small, no eigenvalues will be found. Furthermore, because `computeEigenvalues` employ a Newton-like rootfinding method (provided by `NLSolve`) that is based on an initial guess for the root, choosing this guess "too far" from the correct solution might not converge to a root or it might be that the root is unstable and any small perturbation around an initial guess produces wildly different results. That being said, we can still outline a general procedure that works empirically when finding quasinormal modes based on the different problem types.

First, when the optimal values of ``n`` and ``x_0`` are unknown, start with ``x_0`` in the midpoint of the compactified domain and ``n`` around 20 or 30. This number of iterations will not yield the the most accurate results but it will be enough to determine if we are on the right track while also not being too computationally expensive. From here, we can take one of two different paths:

1. If we have a `QuadraticEigenvalueProblem`, we will have a list of several eigenvalue candidates that are roots of the ``\delta_n`` polynomial but are not necessarily eigenvalues of the ODE. To determine the true eigenvalues, we need to call `computeEigenvalues` repeatedly with an increasing number of AIM iterations. Eigenvalues that persist or change slowly when the number of iterations changes are very likely to be true eigenvalues of the ODE. Other values are likely to be spurious numerical results. This procedure is similar to the one employed when computing eigenvalues using pseudospectral methods: Various spurious results are produced and the true ones are found by repeatedly refining and comparing results. Once true eigenvalues start to emerge, we can start to play around with ``x_0`` to see if more eigenvalues emerge in the list.

2. If we have a `NumericAIMProblem`, a call to `computeEigenvalues` can only produce a single eigenvalue based on an initial guess. Assuming that the `NLSolve` actually converges to a solution this mode is also under the peril of returning spurious results. Here, the wisdom of the `QuadraticEigenvalueProblem`s remains: True results must be refined when the number of iterations increase (indicating numerical convergence). If the returned eigenvalue changes wildly for a fixed initial guess this might indicate that the result is spurious. Once an eigenvalue is found, fine tuning to ``x_0`` can be made. A good value for ``x_0`` will make `NLSolve` converge to a root faster (with less iterations) than a bad one. Also, note that the optimal ``x_0`` value for a certain eigenvalue might not be optimal for all eigenvalues in the spectrum of the ODE (this has been observed empirically). This means that if we are sure that there is an eigenvalue in the vicinity of an initial guess (because we have obtained it with another method, for instance) and `computeEigenvalues` cannot find it even when the number of AIM iterations is high, tuning ``x_0`` might make these modes emerge. Furthermore, in `NumericAIMProblem`s the function `computeDelta!` is a point-wise function that returns the value of ``\delta_n`` anywhere in the complex plane. Using this function, the user can employ a different root finding method than `NLSolve`. This flexibility allows one to eliminate the additional degree of freedom imposed by the initial guess. We can, for instance, use [RootsAndPoles.jl](https://github.com/fgasdia/RootsAndPoles.jl) to find all roots of ``\delta_n`` or any other root finding method desired. An implementation of this idea is presented [here](https://github.com/lucass-carneiro/QuasinormalModes.jl/blob/master/examples/schwarzschild_roots_and_poles.jl) and [here](https://github.com/lucass-carneiro/QuasinormalModes.jl/blob/master/examples/harmonic_oscillator_roots_and_poles.jl).

!!! note "Determining initial guesses"
	Finding initial guesses to supply to `NumericAIMProblem`s can be difficult when solving a new physics problem. If possible, one could firs try and find eigenvalue candidates implementing the problem of interest as a `QuadraticEigenvalueProblem` or extending the code to work semi-analytically with other function types. Furthermore, one could guess a reasonable region where modes would be and use a root bracketing scheme, as described [here](https://lucass-carneiro.github.io/QuasinormalModes.jl/dev/schw/#Interpreting-SolverResults-output) and exemplified [here](https://github.com/lucass-carneiro/QuasinormalModes.jl/blob/master/examples/harmonic_oscillator_roots_and_poles.jl) and [here](https://github.com/lucass-carneiro/QuasinormalModes.jl/blob/master/examples/schwarzschild_roots_and_poles.jl)