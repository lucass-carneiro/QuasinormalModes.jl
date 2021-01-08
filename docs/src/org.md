# Package organization

`QuasinormalModes.jl` employs two main strategies in order to find eigenvalues using the AIM: problems can be solved in a semi-analytic or purely numeric fashion. We make use of Julia's type system in order to implement structures that reflect these operation modes. All of the package's exported functionality is designed to operate on sub-types of abstract types that reflect the desired solution strategy (semi-analytic or numeric). The user is responsible for constructing concrete types that are sub-types of the exported abstract types with the actual problem specific information. It's thus usefully to start by inspecting the package's exported type hierarchy:

## The type hierarchy

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
2. `AnalyticAIMProlem` is the parent type of all problems that can be solved using a semi-analytic approach.
3. `NumericAIMProblem` is the parent type of all problems that can be solved using a numeric approach.
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

Failure to implement these functions returns an error with the appropriate message. Note that these traits only check that such functions are implemented for a certain problem type and not that they follow a particular implementation pattern. The contract on the functions implementations is *soft* and will be clarified further on. Failure to abide by these soft contracts results in undefined behaviour.

## Extending the default functionality

The following assumes that the package `SymEngine` is installed. If a problem type `P{N,T}` is a sub-type of `AnalyticAIMProblem{N,T}`, the user must extend the default implementations abiding by the following rules
1. `QuasinormalModes.λ0(p::P{N,T}) where {N,T}` must return a `SymEngine.Basic` objects representing the symbolic expression for the `λ0` part of the ODE.
2. `QuasinormalModes.S0(p::P{N,T}) where {N,T}` must return a `SymEngine.Basic` objects representing the symbolic expression for the `S0` part of the ODE.
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
