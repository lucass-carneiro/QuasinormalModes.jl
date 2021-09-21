# ------------------------------------------------------------------
# Compute the code coverage locally.
# See https://github.com/JuliaCI/Coverage.jl
# ------------------------------------------------------------------
using Coverage

coverage = process_folder()
coverage = merge_coverage_counts(coverage, filter!(
    let prefixes = joinpath(pwd(), "src", "")
    c -> any(p -> startswith(c.filename, p), prefixes)
end,
    LCOV.readfolder("test")))

covered_lines, total_lines = get_summary(coverage)

println("Total lines: ", total_lines)
println("Covered lines: ", covered_lines)
println("Missing lines: ", total_lines - covered_lines)
println("Code convergence: ", covered_lines / total_lines)
