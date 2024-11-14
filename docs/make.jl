using BioSequenceMappings
using Documenter

DocMeta.setdocmeta!(
    BioSequenceMappings, :DocTestSetup, :(using BioSequenceMappings)
    ; recursive=true
)

makedocs(;
    modules=[BioSequenceMappings],
    authors="Pierre Barrat-Charlaix",
    sitename="BioSequenceMappings.jl",
    format=Documenter.HTML(;
        canonical="https://pierrebarrat.github.io/BioSequenceMappings.jl",
        edit_link="master",
        assets=String[],
    ),
    pages = [
        "Quickstart" => "index.md",
        "Manual" => [
            "Alphabets" => "alphabets.md",
            "Alignments" => "alignments.md",
            "Utilities" => "utilities.md"
        ],
        "Reference" => "reference.md"
    ],
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/PierreBarrat/BioSequenceMappings.jl.git",
    devbranch="master",
)
