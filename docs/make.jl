using Documenter, Mathieu, PyPlot

DocMeta.setdocmeta!(Mathieu, :DocTestSetup, :(using Mathieu); recursive=true)

makedocs(
    sitename = "Mathieu.jl",
    authors = "Jérémy Béjanin",
    modules = [Mathieu],
    format = Documenter.HTML(prettyurls = false),
    linkcheck = true,
    clean = false,
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "examples/mathieu_funs.md",
            "examples/mathieu_surf.md",
        ],
        "API" => [
            "api/values.md",
            "api/functions.md"
        ],
    ]
)

deploydocs(
    repo = "github.com/jebej/Mathieu.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
