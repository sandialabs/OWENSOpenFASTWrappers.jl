# OpenFAST Artifacts

OWENSOpenFASTWrappers loads native OpenFAST libraries through `OWENSOpenFAST_jll`. Julia resolves the matching artifact for the current platform during package installation or test setup, then the wrapper functions pass those library paths to `ccall`.

Use `openfastLibraryPaths()` to inspect the resolved shared libraries:

```julia
using OWENSOpenFASTWrappers

paths = openfastLibraryPaths()
paths.aerodyn_inflow
paths.hydrodyn
paths.inflowwind
paths.moordyn
```

For install debugging and platform smoke tests, use `openfastLibraryArtifactStatus()`:

```julia
status = openfastLibraryArtifactStatus()
status.inflowwind.exists
status.inflowwind.can_load
```

Each status entry contains:

| Field | Meaning |
| --- | --- |
| `path` | Absolute path resolved from `OWENSOpenFAST_jll`. |
| `exists` | Whether that shared-library file exists locally. |
| `can_load` | Whether `Libdl` can load the library on the current platform. |

If `exists` is false, re-run `Pkg.instantiate()` or reinstall `OWENSOpenFAST_jll`. If `can_load` is false while the file exists, the platform artifact is present but the dynamic loader rejected it; check the operating system, architecture, and any dependent shared-library errors from a direct wrapper test.

Most users should keep the artifact-provided libraries. For local OpenFAST development, the wrapper entry points still accept explicit library paths, such as `ifwinit(inflowlib_filename=...)`, `HD_Init(hdlib_filename=...)`, `MD_Init(mdlib_filename=...)`, and the `adi_lib` argument to `setupTurb`.
