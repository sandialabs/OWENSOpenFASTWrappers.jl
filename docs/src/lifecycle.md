# Wrapper Lifecycle

OpenFAST modules keep internal native state. The Julia wrappers expose that state through explicit init, calculate-output, update-state, and end calls.

## Standard Pattern

| Step | InflowWind | AeroDyn | HydroDyn | MoorDyn |
| --- | --- | --- | --- | --- |
| initialize | `ifwinit` | `adiInit` | `HD_Init` | `MD_Init` |
| calculate outputs | `ifwcalcoutput` | `adiCalcOutput` | `HD_CalcOutput` | `MD_CalcOutput` |
| update states | not applicable for simple IFW calls | `adiUpdateStates` | `HD_UpdateStates` | `MD_UpdateStates` |
| cleanup | `ifwend` | `adiEnd` | `HD_End` | `MD_End` |

Call the matching cleanup function after every successful initialization. If a run errors after several modules are initialized, call `OWENSOpenFASTWrappers.endAll()` before retrying in the same Julia process.

## State Rules

- Do not call a `CalcOutput` or `UpdateStates` function before its corresponding init function.
- Re-initialize when input files, mesh topology, output channels, library filenames, or platform state assumptions change.
- Avoid running independent OpenFAST wrapper cases concurrently in the same Julia process unless each wrapper has been audited for independent native handles.
- Keep generated OpenFAST input and output files in test or run-specific directories so coupled cases do not overwrite each other.

## Error Handling

The wrappers call native libraries and then convert outputs back to Julia arrays and structures. When debugging, capture both Julia exceptions and OpenFAST-generated log/output files. A wrapper error after initialization should be treated as requiring cleanup before another case is attempted.
