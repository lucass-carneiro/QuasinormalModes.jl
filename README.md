# QuasinormalModes.jl

A [Julia](http://julialang.org) package for computing Quasinormal Modes or other eigenvalues of second order ordinary differential equations.

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://lucass-carneiro.github.io/QuasinormalModes.jl/)

# Author
[Lucas T. Sanches](lucas.t@ufabc.edu.br), Centro de Ciências Naturais e Humanas, Universidade Federal do ABC (UFABC).

# License

`QuasinormalModes` is licensed under the [MIT license](./LICENSE.md).

# Installation

Currently this package is not registered within Julia's ecosystem. To use it, clone the git repository and it's parent directory to Julia's load path:

```julia
julia> ;
shell> git clone https://github.com/lucass-carneiro/QuasinormalModes.jl
julia> push!(LOAD_PATH, "./")
julia> using QuasinormalModes
```

# Using

For usage instruction please read the [documentation](https://lucass-carneiro.github.io/QuasinormalModes.jl/).
Exemple scripts can be found in the [examples](./examples/) folder

# Contributing

There are many ways to contribute to this package:

- Report an issue if you encounter some odd behavior, or if you have suggestions to improve the package.
- Contribute with code addressing some open issues, that add new functionality or that improve the performance.
- When contributing with code, add docstrings and comments, so others may understand the methods implemented.
- Contribute by updating and improving the documentation.
