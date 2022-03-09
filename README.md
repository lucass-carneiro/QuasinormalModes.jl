# QuasinormalModes.jl

This is a [Julia](http://julialang.org) package whose primary objective is to compute the discrete eigenvalues of second order ordinary differential equations. It was written with the intent to be used for computing quasinormal modes (QNMs) of black holes in General Relativity efficiently and accurately. QNMs are the discrete spectrum of characteristic oscillations produced by black holes when perturbed. These oscillations decay exponentially in time and thus it's said that QNMs contain a real ``\omega_R`` oscillation frequency and an imaginary ``\omega_I`` frequency that represents the mode's decay rate. These frequencies are often described by a discrete eigenvalue in a second order ODE. For a comprehensive review see [[1]](https://arxiv.org/abs/0905.2975).

To compute eigenvalues (and thus quasinormal frequencies) this package uses the Asymptotic Iteration Method (AIM) [[2]](https://arxiv.org/abs/math-ph/0309066v1), more specifically the "improved" version of the AIM as described in [[3]](https://arxiv.org/abs/1111.5024). The AIM can be used to find the eigenvectors and eigenvalues of any second order differential equation (the class of problems with which the quasinormal modes belong) and thus this package can be used not only in the context of General Relativity but can also be used to find the discrete eigenvalues of other systems such as the eigenenergies of a quantum system described by the time independent Schrödinger equation.

[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://lucass-carneiro.github.io/QuasinormalModes.jl/stable)
![Build Status](https://github.com/lucass-carneiro/QuasinormalModes.jl/actions/workflows/CI.yml/badge.svg)
[![codecov](https://codecov.io/gh/lucass-carneiro/QuasinormalModes.jl/branch/master/graph/badge.svg?token=GK9052NQK2)](https://codecov.io/gh/lucass-carneiro/QuasinormalModes.jl)

# Author
[Lucas T. Sanches](lucas.t@ufabc.edu.br), Centro de Ciências Naturais e Humanas, Universidade Federal do ABC (UFABC).

# License

`QuasinormalModes` is licensed under the [MIT license](./LICENSE.md).

# Installation

This package can be installed using the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run

```julia
pkg> add QuasinormalModes
```
and then type backspace to exit back to the REPL.

# Using

For detailed usage instructions please read the [documentation](https://lucass-carneiro.github.io/QuasinormalModes.jl/).

# Contributing

There are many ways to contribute to this package:

- Report an issue if you encounter some odd behavior, or if you have suggestions to improve the package.
- Contribute with code addressing some open issues, that add new functionality or that improve the performance.
- When contributing with code, add docstrings and comments, so others may understand the methods implemented.
- Contribute by updating and improving the documentation.
