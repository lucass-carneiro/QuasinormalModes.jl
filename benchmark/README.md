# Benchmarks
This folder contains all the necessary resources to reproduce the benchmark and error convergence results presented in the JOSS paper:

1. `benchmark.jl`: The Julia script that produces the benchmark results. The benchmark consists in a user configurable number of "runs". Each run computes a quasinormal frequency with a varying (and also configurable) number of AIM iterations. The function that performs the benchmark is called twice, so that Julia's compilation time is not erroneously included at the first iteration.
2. `perf.py`: This Python script takes as an input a number of runs and produces a plot of the average runtime versus the number of iterations from a set of benchmark result files.
3. `err.py`: This Python script produces a plot of the error of the computed frequencies versus the number of AIM iterations from a set of benchmark result files. Since the convergence rate does not change with the number of runs, it computes the plot using the data for the first run.

Additionally, it contains the raw data files with the benchmark results that were used in the production of the JOSS paper compressed inside the `runs.tar.gz` file, alongside the original `.pdf` images included in the paper.

# Hardware and commands

The included data files were produced using 16 threads on a Intel(R) Core(TM) i9-7900X @ 3.30GHz CPU with 256 bit precision floating point numbers. The command issued to produce the results was

```
julia -O3 --startup-file=no --threads=16 benchmark.jl
```

Once the data files were produced with the command above, the `perf.pdf` plot was created by issuing

```
python perf.py 20
```

and finally, the convergence rate plot was produced with

```
python err.py
```
