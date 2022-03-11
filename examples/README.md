# QuasinormalModes.jl Examples

In this folder you can find examples of `QuasinormalModes.jl` in action. Here is a list and brief description of each example:

1. `schwarzschild.jl`: Shows how to compute the quasinormal modes of a Schwarzschild black hole both numerically and semi-analytically. A complete description of how to assemble this example can be found [here](https://lucass-carneiro.github.io/QuasinormalModes.jl/dev/schw/).
2. `harmonic_oscillator.jl`: Shows how to obtain the eigenenergies from a quantum harmonic oscillator. A complete description of how to assemble this example can be found [here](https://lucass-carneiro.github.io/QuasinormalModes.jl/dev/sho/).
3. `schwarzschild_roots_and_poles.jl` Shows how one can use the `computeDelta` function to use any root finding scheme, or package, in order to find eigenvalues. This example makes use of the `RootsAndPoles.jl` package, that finds all roots of the AIM quantization condition in a given region of the complex plane.
4. `harmonic_oscillator_roots_and_poles.jl` Shows how one can use the `computeDelta` and `RootsAndPoles.jl` to find the first few eigenenergies of the quantum harmonic oscillator.

Note that this folder has it's own `Project.toml` and `Manifest.toml`. This allows tracking example dependencies separately from the main package code.