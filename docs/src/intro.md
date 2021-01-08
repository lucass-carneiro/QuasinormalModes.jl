# Introduction

This package's primary objective is to compute the discrete eigenvalues of second order ordinary differential equations. It was written with the intent to be used for computing qusinormal modes (QNMs) of black holes in General Relativity efficiently and accurately. QNMs are the discrete spectrum of characteristic oscillations produced by black holes when perturbed. These oscillations decay exponentially in time and thus it's said that QNMs contain a real ``\omega_R`` oscillation frequency and an imaginary ``\omega_I`` frequency that represents the mode's decay rate. These frequencies are often described by a discrete eigenvalue in a second order ODE. For a comprehensive review see [[1]](https://arxiv.org/abs/0905.2975).

To compute eigenvalues (and thus qusinormal frequencies) this package uses the Asymptotic Iteration Method (AIM) [[2]](https://arxiv.org/abs/math-ph/0309066v1), more specifically the "improved" version of the AIM as described in [[3]](https://arxiv.org/abs/1111.5024). The AIM can be used to find the eigenvectors and eigenvalues of any second order differential equation (the class of problems with which the quasi normal modes belong) and thus this package can be used not only in the context of General Relativity but can also be used to find the discrete eigenvalues of other systems such as the eigenenergies of a quantum system described by the time independent Schrödinger equation.

In the following sections you will find the QuasinormalModes.jl API and instructions on how to use it in a series of (hopefully sufficient) examples.

# Installing

This package can be installed using the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run

```julia
pkg> add QuasinormalModes
```
and then type backspace to exit back to the REPL.
