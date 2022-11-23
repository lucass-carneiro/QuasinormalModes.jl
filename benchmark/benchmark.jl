using QuasinormalModes

# ------------------------------------------------------------------
# 1. Setting up Schwarzschild black hole data structure for
#    computing the QNMs numerically.
# ------------------------------------------------------------------

struct NSchwarzschildData{N,T} <: NumericAIMProblem{N,T}
    nIter::N
    x0::T
    l::N
    s::N
end

function NSchwarzschildData(nIter::N, x0::T, l::N, s::N) where {N,T}
    return NSchwarzschildData{N,T}(nIter, x0, l, s)
end

QuasinormalModes.λ0(::NSchwarzschildData{N,T}) where {N,T} = (x,ω) -> (-1 + (2*im)*ω + x*(4 - 3*x + (4*im)*(-2 + x)*ω))/((-1 + x)^2*x)
QuasinormalModes.S0(d::NSchwarzschildData{N,T}) where {N,T} = (x,ω) -> (d.l + d.l^2 + (-1 + d.s^2)*(-1 + x) + (4*im)*(-1 + x)*ω + 4*(-2 + x)*ω^2)/((-1 + x)^2*x)

QuasinormalModes.get_niter(d::NSchwarzschildData{N,T}) where {N,T} = d.nIter
QuasinormalModes.get_x0(d::NSchwarzschildData{N,T}) where {N,T} = d.x0

# ------------------------------------------------------------------
# 2. Benchmark
# ------------------------------------------------------------------

"""
    bench_iter(iter_start, iter_end, x0, l, s, xt, ft, nls_iter, reference_ωr, reference_ωi, repeat)
Measures the time required to compute a pair of quasinormal frequencies as a function of the
number of AIM iterations performed. The result is also compared with a know result in order
to provide a error measurement.
# Input
- `iter_start`: First number of iterations to perform.
- `iter_end`: Last number of iterations.
- `x0`: The point around which AIM functions will be Taylor expanded.
- `l`: The l value of the sought mode.
- `s`: The s value of the sought mode.
- `xt`: The `xtol` value passed to NLSolve.
- `ft`: The `ftol` value passed to NLSolve.
- `nls_iter`: The number of iterations that NLSolve will perform.
- `reference_ωr`: The reference value for the real QNM frequency.
- `reference_ωi`: The reference value for the imaginary QNM frequency.
- `repeat`: The amount of times to repeat the operation.
# Output
    nothing
"""
function bench_iter(iter_start, iter_end, x0, l, s, xt, ft, nls_iter, reference_ωr, reference_ωi, repeat)

    for i in 1:repeat
        println("Benchmark iteration ", i, ":")
        file = open("run_$(i).dat", "w")
        println(file, "# 1:time(ns) 2:iter 3:x0 4:l 5:s 6:xt 7:ft 8:w_r 9:w_i 10:error in w_r 11:error in w_i")
        
        for iter in iter_start:iter_end
            println("    AIM iteration ", iter)
            p_num = NSchwarzschildData(iter, x0, l, s);
            c_num = AIMCache(p_num)

            t0 = time_ns();
            ev = computeEigenvalues(Threaded(), p_num, c_num, typeof(x0)(reference_ωr, reference_ωi), nls_xtol = xt, nls_ftol = ft, nls_iterations = nls_iter)
            elapsed_time = time_ns() - t0

            if(ev.x_converged || ev.f_converged)
                println(file,
                        elapsed_time,          "    ",
                        iter,                  "    ",
                        convert(Float64, x0),  "    ",
                        l,                     "    ",
                        s,                     "    ",
                        convert(Float64, xt),  "    ",
                        convert(Float64, ft),  "    ",
                        ev.zero[1],            "    ",
                        ev.zero[2],            "    ",
                        ev.zero[1] - reference_ωr, "    ",
                        ev.zero[2] - reference_ωi
                        )
                
                flush(file)
            end
        end
        close(file)
    end
end

# We call the function once and discard the results in order so that the compilation time does not get included in the benchmark.
# Reference values are obtained from https://pages.jh.edu/eberti2/ringdown/
bench_iter(convert(UInt32, 1), convert(UInt32, 100), Complex(big"0.39", big"0.0"), 0x00000, 0x00000, big"1.0e-55", big"1.0e-55", 5000000, big"0.2209098781608393", big"-0.2097914341737619", 1)
bench_iter(convert(UInt32, 1), convert(UInt32, 100), Complex(big"0.39", big"0.0"), 0x00000, 0x00000, big"1.0e-55", big"1.0e-55", 5000000, big"0.2209098781608393", big"-0.2097914341737619", 20)
