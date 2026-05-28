# Developer Guide

This guide is for changes to wrappers, tests, examples, and documentation. The native calls are thin enough that small-looking edits can change physics behavior, memory layout, or cleanup state.

## Before Editing

Read the wrapper source and the corresponding test case together:

| Source | Test or fixture |
|:---|:---|
| `src/inflowwind.jl` | `test/ifw_run.jl`, `test/data/ifw` |
| `src/aerodyn.jl` | `test/aerodyn_run.jl`, `test/AeroDynInputs`, `test/HAWT_verification` |
| `src/hydrodyn.jl` | `test/hydrodyn_run.jl`, `test/data/*HydroDyn.dat`, `test/data/*SeaState.dat`, `test/data/potflow` |
| `src/moordyn.jl` | `test/moordyn_run.jl`, `test/data/moordyn_unit.h5` |

Check the native OpenFAST function signature before changing a `ccall` argument type. The wrapper often converts Julia arrays to `Cfloat`, `Cdouble`, or `Cint` immediately at the call boundary.

## Native State

The wrappers keep library handles, native function symbols, error objects, and active flags in module globals. This keeps user calls simple, but it means wrapper changes must preserve cleanup behavior:

1. Open the native library during initialization.
2. Set the active flag only after the library and symbols are ready.
3. Reset and check native error status after each call.
4. On abort-level errors, end all active wrappers before throwing.
5. Close the native library during the matching end function.

When adding tests, prefer `try`/`finally` so cleanup runs even if a native module fails midway through a time loop.

## Input Files and Generated Text

The wrappers pass OpenFAST inputs either as file-path strings or as NUL-separated input text, depending on the module and function path. Keep this distinction explicit in docs and tests.

Do not assume that relative paths inside OpenFAST files are resolved relative to the file itself. The native module may resolve them from the process working directory. Tests that rely on relative paths should set the working directory deliberately or use paths that are valid from the test process.

Generated AeroDyn, OLAF, blade, and InflowWind files are useful for examples but are sensitive to OpenFAST input-format changes. For regression tests and physics comparisons, prefer reviewed fixture files checked into `test`.

## Physics Review Checklist

Before accepting a wrapper behavior change, check:

| Topic | Questions |
|:---|:---|
| Units | Are time, length, speed, force, moment, density, gravity, and angular quantities still in the documented units? |
| Frames | Did any `+X`, `+Y`, `+Z`, hub, blade, platform, or global-frame convention change? |
| Signs | Do forces, moments, azimuth, roll, pitch, yaw, and rotor speed signs match the reference output? |
| Ordering | Are mesh nodes, blade roots, struts, output channels, and six-DOF arrays in the native order expected by OpenFAST? |
| State | Does the wrapper reinitialize when inputs, mesh topology, time step, or output channels change? |
| Cleanup | Can the same test run twice in one Julia process without stale native state? |

## Documentation

Keep docs source in `docs/src`. Update `docs/make.jl` when adding or moving pages. The reference page uses Documenter autodocs, so public docstrings should stay accurate and concise.

For source-level docstrings, document the unit, frame, mutation behavior, and cleanup expectation. Avoid documenting planned APIs as if they exist in the current checkout.

## Local Build

Build documentation from the package root with:

```sh
julia --project=OWENSOpenFASTWrappers.jl/docs OWENSOpenFASTWrappers.jl/docs/make.jl
```

If dependencies are already instantiated locally, this should not need network access. If a build tries to resolve missing packages, instantiate the docs environment outside a no-network workflow and rerun the build.
