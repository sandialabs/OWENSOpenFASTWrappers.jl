# AeroDyn and InflowWind

AeroDyn and InflowWind are the aerodynamic OpenFAST-facing paths in this package. They can be driven directly or through higher-level OWENS setup functions.

## InflowWind

`ifwinit` prepares an InflowWind instance from either generated steady-input content or a TurbSim `.bts` file. `ifwcalcoutput(position, time)` returns velocity at a point and time. `ifwend` releases the native instance.

Important inputs include:

| Input | Meaning |
| --- | --- |
| `HWindSpeed` | steady horizontal wind speed used by generated inputs. |
| `RefHt` | reference height. |
| `RefLength` | reference length used in generated InflowWind files. The current writer scales this by 10 in the file text, matching existing tests. |
| `turbsim_filename` | path to turbulent inflow `.bts` data. |

## AeroDyn

`adiInit`, `adiCalcOutput`, `adiUpdateStates`, and `adiEnd` expose the AeroDyn-Inflow lifecycle. The higher-level `setupTurb` path builds the mesh and orientation inputs used by OWENS HAWT-style verification cases.

Both `adiInit` and `setupTurb` can initialize AeroDyn and InflowWind from filenames or generated input text. The default `*_input_file_passed=0` passes the filename through to the native library. Set `*_input_file_passed=1` with `*_input_source=:file` to read a file in Julia and pass its NUL-separated text to OpenFAST, or set `*_input_source=:text` when the caller already has generated input text or a vector of input lines. This keeps the existing file workflow intact while allowing direct data transfer from OWENS-side input generators.

`setupTurb` accepts a `rotation_direction` keyword so the requested AeroDyn frame convention is explicit. The currently validated combinations are `rotation_direction = :ccw` for VAWT/cross-flow turbines and `rotation_direction = :cw` for HAWT turbines. Unsupported combinations throw before native AeroDyn initialization because the opposite root-order, span-direction, velocity, and load-frame mappings still need dedicated validation.

When preparing an AeroDyn run, keep these inputs together in the test fixture:

- AeroDyn input file;
- InflowWind input file;
- output root name;
- blade shape arrays;
- blade and element index mappings;
- mesh and orientation structures;
- output-channel count and names.

## Coupled Runs

OWENS-coupled AeroDyn runs should pin both OpenFAST output-channel values and OWENS-mapped blade loads. That catches regressions in native wrapper calls, channel ordering, mesh indexing, and force/moment sign conventions.
