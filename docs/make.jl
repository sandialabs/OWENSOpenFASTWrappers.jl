using Documenter, Literate, OWENSOpenFASTWrappers

# Build documentation
makedocs(;
    modules = [OWENSOpenFASTWrappers],
    pages = [
        "Home" => "index.md",
        "API Reference" => joinpath("reference", "reference.md"),
    ],
    sitename = "OWENSOpenFASTWrappers.jl",
    authors = "Kevin R. Moore <kevmoor@sandia.gov>",
)

deploydocs(
    repo = "github.com/sandialabs/OWENSOpenFASTWrappers.jl.git",
)