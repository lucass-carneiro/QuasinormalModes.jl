using QuasinormalModes
using SymEngine
using Test

# ------------------------------------------------------------------
# 1. Analytic Schwarzschild Black Hole
# ------------------------------------------------------------------

struct SchwarzschildData{N,T} <: QuadraticEigenvalueProblem{N,T}
    nIter::N
    x0::T

    vars::Tuple{Basic,Basic}
    exprs::Tuple{Basic,Basic}
end

function SchwarzschildData(nIter::N, x0::T, l::N, s::N) where {N,T}
    vars = @vars x ω

    λ0 = (-1 + (2 * im) * ω + x * (4 - 3 * x + (4 * im) * (-2 + x) * ω)) / ((-1 + x)^2 * x)
    S0 = (l + l^2 + (-1 + s^2) * (-1 + x) + (4 * im) * (-1 + x) * ω + 4 * (-2 + x) * ω^2) / ((-1 + x)^2 * x)

    return SchwarzschildData{N,T}(nIter, x0, vars, (λ0, S0))
end

QuasinormalModes.λ0(d::SchwarzschildData{N,T}) where {N,T} = d.exprs[1]
QuasinormalModes.S0(d::SchwarzschildData{N,T}) where {N,T}  = d.exprs[2]

QuasinormalModes.get_niter(d::SchwarzschildData{N,T}) where {N,T} = d.nIter
QuasinormalModes.get_x0(d::SchwarzschildData{N,T}) where {N,T} = d.x0

QuasinormalModes.get_ODEvar(d::SchwarzschildData{N,T}) where {N,T} = d.vars[1]
QuasinormalModes.get_ODEeigen(d::SchwarzschildData{N,T}) where {N,T} = d.vars[2]

# ------------------------------------------------------------------
# 2. Numeric Schwarzschild Black Hole
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

QuasinormalModes.λ0(::NSchwarzschildData{N,T}) where {N,T} = (x, ω) -> (-1 + (2 * im) * ω + x * (4 - 3 * x + (4 * im) * (-2 + x) * ω)) / ((-1 + x)^2 * x)
QuasinormalModes.S0(d::NSchwarzschildData{N,T}) where {N,T} = (x, ω) -> (d.l + d.l^2 + (-1 + d.s^2) * (-1 + x) + (4 * im) * (-1 + x) * ω + 4 * (-2 + x) * ω^2) / ((-1 + x)^2 * x)

QuasinormalModes.get_niter(d::NSchwarzschildData{N,T}) where {N,T} = d.nIter
QuasinormalModes.get_x0(d::NSchwarzschildData{N,T}) where {N,T} = d.x0

# ------------------------------------------------------------------
# 3. Nmeric Tests
# ------------------------------------------------------------------

@testset "Serial Numeric correctness - BigFloat, 100 iterations, s = 0, l = 0" begin
    p = NSchwarzschildData(0x00064, Complex(BigFloat("0.43"), BigFloat("0.0")), 0x00000, 0x00000);
    c = AIMCache(p)
    
    # --- n = 0 ---
    ev = computeEigenvalues(
        Serial(),
        p,
        c,
        Complex(BigFloat("0.22"), BigFloat("-0.20")),
        nls_xtol=BigFloat("1.0e-50"),
        nls_ftol=BigFloat("1.0e-50")
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - BigFloat("0.2209098781608393")) < BigFloat("1.0e-7")
    @test abs(ev.zero[2] + BigFloat("0.2097914341737619")) < BigFloat("1.0e-7")
end

@testset "Serial Numeric correctness - BigFloat, 100 iterations, s = 1, l = 1" begin
    p = NSchwarzschildData(0x00064, Complex(BigFloat("0.43"), BigFloat("0.0")), 0x00001, 0x00001);
    c = AIMCache(p)
    
    # --- n = 0 ---
    ev = computeEigenvalues(
        Serial(),
        p,
        c,
        Complex(BigFloat("0.49"), BigFloat("-0.18")),
        nls_xtol=BigFloat("1.0e-50"),
        nls_ftol=BigFloat("1.0e-50")
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - BigFloat("0.4965265283562174")) < BigFloat("1.0e-7")
    @test abs(ev.zero[2] + BigFloat("0.1849754359058844")) < BigFloat("1.0e-7")
end

@testset "Serial Numeric correctness - BigFloat, 100 iterations, s = 1, l = 1" begin
    p = NSchwarzschildData(0x00064, Complex(BigFloat("0.43"), BigFloat("0.0")), 0x00002, 0x00002);
    c = AIMCache(p)
    
    # --- n = 0 ---
    ev = computeEigenvalues(
        Serial(),
        p,
        c,
        Complex(BigFloat("0.74"), BigFloat("-0.17")),
        nls_xtol=BigFloat("1.0e-50"),
        nls_ftol=BigFloat("1.0e-50")
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - BigFloat("0.7473433688360838")) < BigFloat("1.0e-7")
    @test abs(ev.zero[2] + BigFloat("0.1779246313778714")) < BigFloat("1.0e-7")
end
    
@testset "Serial Numeric correctness - Float64, 48 iterations, s = 0, l = 0" begin
    p = NSchwarzschildData(0x00030, Complex(0.43, 0.0), 0x00000, 0x00000);
    c = AIMCache(p)
    
    # --- n = 0 ---
    ev = computeEigenvalues(
        Serial(),
        p,
        c,
        Complex(0.22, -0.20),
        nls_xtol=1.0e-10,
        nls_ftol=1.0e-10
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - 0.2209098781608393) < 1.0e-4
    @test abs(ev.zero[2] + 0.2097914341737619) < 1.0e-4

    # --- n = 1 ---
    ev = computeEigenvalues(
        Serial(),
        p,
        c,
        Complex(0.17, -0.69),
        nls_xtol=1.0e-10,
        nls_ftol=1.0e-10
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - 0.1722338366727985) < 1.0e-3
    @test abs(ev.zero[2] + 0.6961048936129209) < 1.0e-3

    # --- n = 2 ---
    ev = computeEigenvalues(
        Serial(),
        p,
        c,
        Complex(0.15, -1.2),
        nls_xtol=1.0e-10,
        nls_ftol=1.0e-10
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - 0.1514838710703517) < 1.0e-2
    @test abs(ev.zero[2] + 0.1202157180071607e1) < 1.0e-2
end

@testset "Serial Numeric correctness - Float64, 47 iterations, s = 1, l = 1" begin
    p = NSchwarzschildData(0x0002F, Complex(0.43, 0.0), 0x00001, 0x00001);
    c = AIMCache(p)
    
    # --- n = 0 ---
    ev = computeEigenvalues(
        Serial(),
        p,
        c,
        Complex(0.49, -0.18),
        nls_xtol=1.0e-10,
        nls_ftol=1.0e-10
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - 0.4965265283562174) < 1.0e-9
    @test abs(ev.zero[2] + 0.1849754359058844) < 1.0e-9

    # --- n = 1 ---
    ev = computeEigenvalues(
        Serial(),
        p,
        c,
        Complex(0.42, -0.58),
        nls_xtol=1.0e-10,
        nls_ftol=1.0e-10
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - 0.4290308391272117) < 1.0e-6
    @test abs(ev.zero[2] + 0.5873352910914573) < 1.0e-6

    # --- n = 2 ---
    ev = computeEigenvalues(
        Serial(),
        p,
        c,
        Complex(0.34, -1.05),
        nls_xtol=1.0e-10,
        nls_ftol=1.0e-10
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - 0.3495471352140215) < 1.0e-4
    @test abs(ev.zero[2] + 0.1050375198717648e1) < 1.0e-4
end

@testset "Serial Numeric correctness - Float64, 46 iterations, s = 2, l = 2" begin
    p = NSchwarzschildData(0x0002D, Complex(0.43, 0.0), 0x00002, 0x00002);
    c = AIMCache(p)
    
    # --- n = 0 ---
    ev = computeEigenvalues(
        Serial(),
        p,
        c,
        Complex(0.74, -0.17),
        nls_xtol=1.0e-10,
        nls_ftol=1.0e-10
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - 0.7473433688360838) < 1.0e-9
    @test abs(ev.zero[2] + 0.1779246313778714) < 1.0e-9

    # --- n = 1 ---
    ev = computeEigenvalues(
        Serial(),
        p,
        c,
        Complex(0.69, -0.54),
        nls_xtol=1.0e-10,
        nls_ftol=1.0e-10
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - 0.6934219937583268) < 1.0e-8
    @test abs(ev.zero[2] + 0.5478297505824697) < 1.0e-8

    # --- n = 1 ---
    ev = computeEigenvalues(
        Serial(),
        p,
        c,
        Complex(0.60, -0.95),
        nls_xtol=1.0e-10,
        nls_ftol=1.0e-10
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - 0.6021069092247328) < 1.0e-6
    @test abs(ev.zero[2] + 0.9565539664461437) < 1.0e-6
end

# ------------------------------------------------------------------
# 4. Threaded Nmeric Tests
# ------------------------------------------------------------------

@testset "Threaded Numeric correctness - BigFloat, 100 iterations, s = 0, l = 0" begin
    p = NSchwarzschildData(0x00064, Complex(BigFloat("0.43"), BigFloat("0.0")), 0x00000, 0x00000);
    c = AIMCache(p)
    
    # --- n = 0 ---
    ev = computeEigenvalues(
        Threaded(),
        p,
        c,
        Complex(BigFloat("0.22"), BigFloat("-0.20")),
        nls_xtol=BigFloat("1.0e-50"),
        nls_ftol=BigFloat("1.0e-50")
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - BigFloat("0.2209098781608393")) < BigFloat("1.0e-7")
    @test abs(ev.zero[2] + BigFloat("0.2097914341737619")) < BigFloat("1.0e-7")
end

@testset "Threaded Numeric correctness - BigFloat, 100 iterations, s = 1, l = 1" begin
    p = NSchwarzschildData(0x00064, Complex(BigFloat("0.43"), BigFloat("0.0")), 0x00001, 0x00001);
    c = AIMCache(p)
    
    # --- n = 0 ---
    ev = computeEigenvalues(
        Threaded(),
        p,
        c,
        Complex(BigFloat("0.49"), BigFloat("-0.18")),
        nls_xtol=BigFloat("1.0e-50"),
        nls_ftol=BigFloat("1.0e-50")
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - BigFloat("0.4965265283562174")) < BigFloat("1.0e-7")
    @test abs(ev.zero[2] + BigFloat("0.1849754359058844")) < BigFloat("1.0e-7")
end

@testset "Threaded Numeric correctness - BigFloat, 100 iterations, s = 1, l = 1" begin
    p = NSchwarzschildData(0x00064, Complex(BigFloat("0.43"), BigFloat("0.0")), 0x00002, 0x00002);
    c = AIMCache(p)
    
    # --- n = 0 ---
    ev = computeEigenvalues(
        Threaded(),
        p,
        c,
        Complex(BigFloat("0.74"), BigFloat("-0.17")),
        nls_xtol=BigFloat("1.0e-50"),
        nls_ftol=BigFloat("1.0e-50")
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - BigFloat("0.7473433688360838")) < BigFloat("1.0e-7")
    @test abs(ev.zero[2] + BigFloat("0.1779246313778714")) < BigFloat("1.0e-7")
end

@testset "Threaded Numeric correctness - Float64, 48 iterations, s = 0, l = 0" begin
    p = NSchwarzschildData(0x00030, Complex(0.43, 0.0), 0x00000, 0x00000);
    c = AIMCache(p)
    
    # --- n = 0 ---
    ev = computeEigenvalues(
        Threaded(),
        p,
        c,
        Complex(0.22, -0.20),
        nls_xtol=1.0e-10,
        nls_ftol=1.0e-10
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - 0.2209098781608393) < 1.0e-4
    @test abs(ev.zero[2] + 0.2097914341737619) < 1.0e-4

    # --- n = 1 ---
    ev = computeEigenvalues(
        Threaded(),
        p,
        c,
        Complex(0.17, -0.69),
        nls_xtol=1.0e-10,
        nls_ftol=1.0e-10
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - 0.1722338366727985) < 1.0e-3
    @test abs(ev.zero[2] + 0.6961048936129209) < 1.0e-3

    # --- n = 2 ---
    ev = computeEigenvalues(
        Threaded(),
        p,
        c,
        Complex(0.15, -1.2),
        nls_xtol=1.0e-10,
        nls_ftol=1.0e-10
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - 0.1514838710703517) < 1.0e-2
    @test abs(ev.zero[2] + 0.1202157180071607e1) < 1.0e-2
end

@testset "Threaded Numeric correctness - Float64, 47 iterations, s = 1, l = 1" begin
    p = NSchwarzschildData(0x0002F, Complex(0.43, 0.0), 0x00001, 0x00001);
    c = AIMCache(p)
    
# --- n = 0 ---
    ev = computeEigenvalues(
        Threaded(),
        p,
        c,
        Complex(0.49, -0.18),
        nls_xtol=1.0e-10,
        nls_ftol=1.0e-10
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - 0.4965265283562174) < 1.0e-9
    @test abs(ev.zero[2] + 0.1849754359058844) < 1.0e-9

    # --- n = 1 ---
    ev = computeEigenvalues(
        Threaded(),
        p,
        c,
        Complex(0.42, -0.58),
        nls_xtol=1.0e-10,
        nls_ftol=1.0e-10
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - 0.4290308391272117) < 1.0e-6
    @test abs(ev.zero[2] + 0.5873352910914573) < 1.0e-6

    # --- n = 2 ---
    ev = computeEigenvalues(
        Threaded(),
        p,
        c,
        Complex(0.34, -1.05),
        nls_xtol=1.0e-10,
        nls_ftol=1.0e-10
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - 0.3495471352140215) < 1.0e-4
    @test abs(ev.zero[2] + 0.1050375198717648e1) < 1.0e-4
end

@testset "Threaded Numeric correctness - Float64, 46 iterations, s = 2, l = 2" begin
    p = NSchwarzschildData(0x0002D, Complex(0.43, 0.0), 0x00002, 0x00002);
    c = AIMCache(p)
    
    # --- n = 0 ---
    ev = computeEigenvalues(
        Threaded(),
        p,
        c,
        Complex(0.74, -0.17),
        nls_xtol=1.0e-10,
        nls_ftol=1.0e-10
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - 0.7473433688360838) < 1.0e-9
    @test abs(ev.zero[2] + 0.1779246313778714) < 1.0e-9

    # --- n = 1 ---
    ev = computeEigenvalues(
        Threaded(),
        p,
        c,
        Complex(0.69, -0.54),
        nls_xtol=1.0e-10,
        nls_ftol=1.0e-10
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - 0.6934219937583268) < 1.0e-8
    @test abs(ev.zero[2] + 0.5478297505824697) < 1.0e-8

    # --- n = 1 ---
    ev = computeEigenvalues(
        Threaded(),
        p,
        c,
        Complex(0.60, -0.95),
        nls_xtol=1.0e-10,
        nls_ftol=1.0e-10
        )
    
    @test (ev.x_converged || ev.f_converged) == true
    @test abs(ev.zero[1] - 0.6021069092247328) < 1.0e-6
    @test abs(ev.zero[2] + 0.9565539664461437) < 1.0e-6
end

# ------------------------------------------------------------------
# 5. Grid search tests
# ------------------------------------------------------------------

@testset "Serial grid search - Float64, 46 iterations, s = 2, l = 2" begin
    p = NSchwarzschildData(0x0002D, Complex(0.43, 0.0), 0x00002, 0x00002);
    c = AIMCache(p)
    ev = eigenvaluesInGrid(Serial(), p, c, (Complex(0.60, -0.95), Complex(0.74, -0.17), 3, 3))
    
    cutoff = 1.0e-5

    filter!(x -> real(x) > cutoff && imag(x) < 0.0, ev)

    @test abs(real(ev[1]) - 0.6021069092247328) < 1.0e-6 && abs(imag(ev[1]) + 0.9565539664461437) < 1.0e-6
    @test abs(real(ev[5]) - 0.6934219937583268) < 1.0e-6 && abs(imag(ev[5]) + 0.5478297505824697) < 1.0e-6
    @test abs(real(ev[6]) - 0.7473433688360838) < 1.0e-6 && abs(imag(ev[6]) + 0.1779246313778714) < 1.0e-6
end

@testset "Serial grid search - BigFloat, 46 iterations, s = 2, l = 2" begin
    p = NSchwarzschildData(0x0002D, Complex(big"0.43", big"0.0"), 0x00002, 0x00002);
    c = AIMCache(p)
    ev = eigenvaluesInGrid(Threaded(), p, c, (Complex(big"0.60", big"-0.95"), Complex(big"0.74", big"-0.17"), 3, 3))
    
    cutoff = big"1.0e-5"

    filter!(x -> real(x) > cutoff && imag(x) < big"0.0", ev)

    @test abs(real(ev[1]) - big"0.6021069092247328") < 1.0e-6 && abs(imag(ev[1]) + big"0.9565539664461437") < 1.0e-6
    @test abs(real(ev[5]) - big"0.6934219937583268") < 1.0e-6 && abs(imag(ev[5]) + big"0.5478297505824697") < 1.0e-6
    @test abs(real(ev[6]) - big"0.7473433688360838") < 1.0e-6 && abs(imag(ev[6]) + big"0.1779246313778714") < 1.0e-6
end
