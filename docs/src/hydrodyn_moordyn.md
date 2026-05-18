# HydroDyn and MoorDyn

HydroDyn and MoorDyn wrappers support offshore platform and mooring workflows used by OWENS coupled simulations.

## HydroDyn

The HydroDyn lifecycle is:

```julia
HD_Init(...)
HD_CalcOutput(...)
HD_UpdateStates(...)
HD_End()
```

Typical inputs include the HydroDyn input file, SeaState input file, potential-flow file, water density, water depth, mean sea level offset, time step, and maximum simulation time. The tests under `test/hydrodyn_run.jl` show the expected fixture layout.

## MoorDyn

The MoorDyn lifecycle is:

```julia
MD_Init(...)
MD_CalcOutput(...)
MD_UpdateStates(...)
MD_End()
```

Typical inputs include the MoorDyn input file, water density, water depth, initial platform position, gravity, time step, and interpolation order. The tests under `test/moordyn_run.jl` show the current calling pattern.

## Coupling Notes

For offshore OWENS work, hydrodynamic and mooring forces should be tested at the wrapper boundary and after structural-frame mapping. Platform-position ordering, moment reference points, interpolation order, and state update timing can all change loads without necessarily causing a native-library error.
