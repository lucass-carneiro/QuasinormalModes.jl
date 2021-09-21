using SafeTestsets

@time begin
    @time @safetestset "Quantum Harmonic oscilator tests" begin include("harmonic_oscilator.jl") end
    @time @safetestset "Schwarzschild quasinormal modes tests" begin include("schwarzschild_qnm.jl") end
end
