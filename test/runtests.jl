using SafeTestsets

@time begin
    @time @safetestset "Schwarzschild tests" begin include("eigenvalues/schwarzschild_tests.jl") end
    @time @safetestset "Quantum Harmonic oscilator tests" begin include("eigenvalues/harmonic_oscilator_tests.jl") end
end