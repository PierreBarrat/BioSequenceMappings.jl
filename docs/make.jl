using BioSequenceMappings
using Documenter

DocMeta.setdocmeta!(BioSequenceMappings, :DocTestSetup, :(using BioSequenceMappings); recursive=true)

makedocs(;
    modules=[BioSequenceMappings],
    authors="Pierre Barrat-Charlaix",
    sitename="BioSequenceMappings.jl",
    format=Documenter.HTML(;
        canonical="https://PierreBarrat.github.io/BioSequenceMappings.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/PierreBarrat/BioSequenceMappings.jl",
    devbranch="master",
)
