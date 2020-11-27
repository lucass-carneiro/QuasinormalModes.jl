if !("../src/" in LOAD_PATH)
    push!(LOAD_PATH, "../src/")
end

using Documenter
using QuasinormalModes

makedocs(sitename="QuasinormalModes.jl",
         modules = [QuasinormalModes],
         pages = [
             "index.md",
             "intro.md",
             "schw.md",
             "sho.md",
             "api_ref.md"
             ]
         )

deploydocs(repo = "https://github.com/lucass-carneiro/QuasinormalModes.jl.git")
