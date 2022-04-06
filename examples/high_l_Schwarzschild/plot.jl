using Plots
using CSV
using DataFrames
using LaTeXStrings

file = CSV.File("high_l_Schwarzschild/qnm_2022-04-05 17:49:10_s_0.dat", header=["iter", "l", "n", "real_guess", "imag_guess", "real", "imag"], delim=' ', skipto=2, ignorerepeated=true)
df = DataFrame(file)

function plot_ns(save=true)
    dfs = [filter(row -> row.n == n, df) for n in 0:3]
    
    p1 = plot(dfs[1].real,
              -dfs[1].imag,

              seriestype=:scatter,

              xlabel=L"\Re(\omega)",
              ylabel=L"-\Im(\omega)",

              label=L"n = 0",

              frame=:box,
              grid=false,

              color=:red,

              fontfamily="Computer Modern",

              legend=:outertopright
              )
    
    p2 = plot!(p1, dfs[2].real, -dfs[2].imag, seriestype=:scatter, label=L"n = 1", color=:black)
    p3 = plot!(p2, dfs[3].real, -dfs[3].imag, seriestype=:scatter, label=L"n = 2", color=:blue)
    p4 = plot!(p3, dfs[4].real, -dfs[4].imag, seriestype=:scatter, label=L"n = 3", color=:green)

    if save
        savefig(p4, "high_l_Schwarzschild/l_plot.pdf")
    end
    
    return p4
end

plot_ns()
