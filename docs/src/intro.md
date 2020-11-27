# Introduction

This package's primary objective is to compute qusinormal modes (QNMs) black holes in General Ralativity fast an accurately. QNMs are the characteristic oscillations produced by black holes when perturbed. These oscillations decay exponentially in time and thus it's said that QNMs contain a real ``\omega_R`` oscillation frequency and an imaginary ``\omega_I`` frequency that represents the mode's decay rate. For a comprehensive review see [[1]](https://arxiv.org/abs/0905.2975).

To compute qusinormal frequencies this package uses the Asymptotic Iteration Method (AIM) [[2]](https://arxiv.org/abs/math-ph/0309066v1), more specifically the "improved" version of the AIM as described in [[3]](https://arxiv.org/abs/1111.5024). The AIM can be used to find the eigenvectors and eigenvalues of any second order differential equation (the class of problems with which the quasi normal modes belong) and thus this package can be used not only in the context of General Relativity but can also be used in to find the eigenvalues of any 2nd order ODE such as finding eigenenergies  of a quantum system described by the time independent Schrödinger equation.

In the following sections we will describe QuasinormalModes.jl API and how to use it in a series of (hopefully sufficient) examples so that new users can quickly set up their equations within the infrastructure and get results.

# Installing

This package can be installed using the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run

```julia
pkg> add QuasinormalModes
```
and then type backspece to exit back to the REPL.

# Arbtrary precision
By default, all of the computations performed by this package are carried out using arbitrary precision arithmetic. This is desirable since in root finding catastrophic numerical cancellations can occur if only floating point accuracy is used. To control the number of significant digits globally use Julia's [setprecision](https://docs.julialang.org/en/v1/base/numbers/#Base.MPFR.setprecision) and [setrounding](https://docs.julialang.org/en/v1/base/numbers/#Base.Rounding.setrounding-Tuple{Type,Any}) methods. We also recommend the reading of the [Arbitrary Precision Arithmetic](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/#Arbitrary-Precision-Arithmetic) section of the official documentation.
