# Quickstart

This package wraps native OpenFAST shared libraries through `OWENSOpenFAST_jll`. Use the wrappers from Julia when a workflow needs direct library calls instead of launching a complete OpenFAST executable.

## Install

```julia
using Pkg
Pkg.add(PackageSpec(url = "https://github.com/sandialabs/OWENSOpenFASTWrappers.jl.git"))
```

For toolkit development from sibling checkouts:

```julia
using Pkg
Pkg.develop(path = "../OWENSOpenFASTWrappers.jl")
```

The docs project uses local package sources and should be built with Julia 1.11
or newer.

## Minimal Lifecycle

All library families follow the same broad lifecycle:

```julia
using OWENSOpenFASTWrappers

# Example shape only: choose the concrete init function for the library being used.
# init(...)
# values = calc_output(position_or_mesh, time)
# update_states(...)
# end_function()
```

For InflowWind test fixtures, the concrete form is:

```julia
using OWENSOpenFASTWrappers

turbsim_filename = joinpath(pkgdir(OWENSOpenFASTWrappers), "test", "data", "ifw", "test.bts")
OWENSOpenFASTWrappers.ifwinit(; turbsim_filename)
velocity = OWENSOpenFASTWrappers.ifwcalcoutput([0.0, 0.0, 100.0], 0.1)
OWENSOpenFASTWrappers.ifwend()
```

Always end the initialized library. `OWENSOpenFASTWrappers.endAll()` is available as a cleanup guard when a test or coupled run may have more than one active wrapper.

## First Checks

Run:

```julia
using Pkg
Pkg.test("OWENSOpenFASTWrappers")
```

The full test suite touches native libraries and may take longer than pure Julia unit tests. Unit tests under `test/pure_helpers_unit.jl` cover helper text generation and state-independent behavior.
