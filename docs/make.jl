using ADAPT
using Documenter

DocMeta.setdocmeta!(ADAPT, :DocTestSetup, :(using ADAPT); recursive=true)

makedocs(;
    modules=[ADAPT],
    authors="Kyle Sherbert <kyle.sherbert@vt.edu> and contributors",
    repo="https://github.com/kmsherbertvt/ADAPT.jl/blob/{commit}{path}#{line}",
    sitename="ADAPT.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kmsherbertvt.github.io/ADAPT.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kmsherbertvt/ADAPT.jl",
    devbranch="main",
)
