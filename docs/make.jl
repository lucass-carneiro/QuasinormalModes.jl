if !("./" in LOAD_PATH)
    push!(LOAD_PATH, "./")
end

using Documenter
using QuasinormalModes

makedocs(sitename="QuasinormalModes.jl",
         modules = [QuasinormalModes],
         pages = [
             "index.md",
             "intro.md",
             "org.md",
             "schw.md",
             "sho.md",
             "api_ref.md"
             ]
         )

# To view in browser start a server in the build dir:
# python -m http.server --bind localhost
