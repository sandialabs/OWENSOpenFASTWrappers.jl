using Documenter, Literate, OWENSOpenFASTWrappers

# Build documentation
makedocs(;
    modules = [OWENSOpenFASTWrappers],
    pages = [
        "Home" => "index.md",
        "Quickstart" => "quickstart.md",
        "OpenFAST Artifacts" => "openfast_artifacts.md",
        "Wrapper Lifecycle" => "lifecycle.md",
        "AeroDyn and InflowWind" => "aerodyn_inflowwind.md",
        "HydroDyn and MoorDyn" => "hydrodyn_moordyn.md",
        "Frames, Units, and Validation" => joinpath("theory", "frames_units_validation.md"),
        "API Reference" => joinpath("reference", "reference.md"),
    ],
    sitename = "OWENSOpenFASTWrappers.jl",
    authors = "Kevin R. Moore <kevmoor@sandia.gov>",
    remotes = nothing,
    format = Documenter.HTML(
        repolink = "https://github.com/sandialabs/OWENSOpenFASTWrappers.jl",
        edit_link = "master",
    ),
)

if get(ENV, "CI", "false") == "true"
    deploydocs(
        repo = "github.com/sandialabs/OWENSOpenFASTWrappers.jl.git",
        devbranch = "master",
    )
end
