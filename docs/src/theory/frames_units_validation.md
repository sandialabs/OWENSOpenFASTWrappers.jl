# Frames, Units, and Validation

The wrappers inherit OpenFAST module conventions. Most quantities are SI, but each module has its own coordinate, mesh, and channel-order assumptions.

## Units and State

| Quantity | Convention |
| --- | --- |
| time | seconds |
| position and displacement | meters |
| velocity | meters per second |
| angular velocity | radians per second unless inherited input-file text states otherwise |
| force | newtons |
| moment | newton meters |
| density | kg/m^3 |
| gravity | m/s^2 |

Input files may contain additional unit annotations. Treat the OpenFAST input file as the source of truth when a field name is ambiguous.

## Frame Checks

Wrapper tests should verify:

- point and mesh coordinates in the expected OpenFAST module frame;
- blade element ordering and blade index mapping;
- AeroDyn rotation direction, including whether the blade root starts at the top or bottom of a VAWT blade and whether local span points upward or downward;
- force and moment signs after converting to OWENS frames;
- platform position and moment reference for HydroDyn and MoorDyn;
- InflowWind point locations relative to reference height and reference length.

## Current Test Evidence

| Area | Evidence |
| --- | --- |
| InflowWind | `test/ifw_run.jl` initializes IFW and samples a velocity. |
| AeroDyn | `test/aerodyn_run.jl` and `test/HAWT_verification` exercise the AeroDyn wrapper and HAWT mapping path. `test/hawt_output_parity.jl` parses tracked selected-channel HAWT AeroDyn-driver fixtures and pins wrapper-vs-standalone RMSE/max-error metrics for rotor coefficients, power, hub loads, and a representative blade-node force triplet. |
| HydroDyn | `test/hydrodyn_run.jl` checks HydroDyn initialization and output calls. |
| MoorDyn | `test/moordyn_run.jl` checks MoorDyn initialization and output calls. |
| Pure helpers | `test/pure_helpers_unit.jl` pins state-independent input text and helper behavior. |

`readOpenFASTOutputTable(path)` and
`openfastOutputChannelMetrics(reference, candidate, channels; rows=...)` are
the reusable pieces behind those checks. They are intentionally table-based so
future validation work can compare native OpenFAST outputs, wrapper outputs,
and reduced examples without plotting or GUI dependencies. When a wrapper
initialization row is known to differ from the standalone driver, pass an
explicit row range and document why the excluded row is not part of the physics
comparison.

## Hardening Rules

- Pin concrete channel values, vector lengths, and element types.
- Keep generated files isolated by test case.
- Call the relevant end function in `try`/`finally` blocks when adding tests.
- Add frame-transform tests around OWENS-coupled loads, not only native-library smoke tests.
- Document whether a fixture is a unit test, wrapper regression, or physics validation case.
