using Landscapes
using Documenter

DocMeta.setdocmeta!(Landscapes, :DocTestSetup, :(using Landscapes); recursive=true)

makedocs(;
    modules=[Landscapes],
    authors="Rafael Schouten <rafaelschouten@gmail.com>",
    repo="https://github.com/rafaqz/Landscapes.jl/blob/{commit}{path}#{line}",
    sitename="Landscapes.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://rafaqz.github.io/Landscapes.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/rafaqz/Landscapes.jl",
    devbranch="main",
)
