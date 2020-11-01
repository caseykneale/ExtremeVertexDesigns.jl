using ExtremeVertexDesigns
using Documenter

makedocs(;
    modules=[ExtremeVertexDesigns],
    authors="Casey Kneale",
    repo="https://github.com/caseykneale/ExtremeVertexDesigns.jl/blob/{commit}{path}#L{line}",
    sitename="ExtremeVertexDesigns.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://caseykneale.github.io/ExtremeVertexDesigns.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/caseykneale/ExtremeVertexDesigns.jl",
)
