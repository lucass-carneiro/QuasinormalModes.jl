using SafeTestsets

@time begin
    @time @safetestset "Quantum Harmonic oscillator tests" begin
        include("harmonic_oscillator.jl")
    end
    @time @safetestset "Schwarzschild quasinormal modes tests" begin
        include("schwarzschild_qnm.jl")
    end
end
