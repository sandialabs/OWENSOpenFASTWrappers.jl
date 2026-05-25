module OWENSOpenFASTWrappers
import Libdl
using LinearAlgebra     # for aerodyn cross-product

using OWENSOpenFAST_jll: libaerodyn_inflow_c_binding, libhydrodyn_c_binding,
    libifw_c_binding, libmoordyn_c_binding, turbsim, aerodyn_driver, hydrodyn_driver, inflowwind_driver, moordyn_driver
const path = splitdir(@__FILE__)[1]
OFWpath = path

# InflowWind
export ifwinit, ifwcalcoutput, ifwend, OFWpath, turbsim, inflowwind_driver
export openfastLibraryPaths, openfastLibraryArtifactStatus

# AD15
export Turbine, Environment, Structure, adiInit, adiCalcOutput, adiUpdateStates, adiEnd, aerodyn_driver
export normalizeADIRotationDirection, rotationDirectionSign, validateADIRotationDirection
export normalizeOpenFASTInputSource, openfastInputString
export readOpenFASTOutputTable, openfastOutputChannelMetrics

# HydroDyn routines
export HD_Init, HD_CalcOutput, HD_UpdateStates, HD_End, hydrodyn_driver

# MoorDyn routines
export MD_Init, MD_CalcOutput, MD_UpdateStates, MD_End, moordyn_driver

include("./aerodyn.jl")
include("./hydrodyn.jl")
include("./inflowwind.jl")
include("./moordyn.jl")
include("./output_tables.jl")

"""
    openfastLibraryPaths() -> NamedTuple

Return the native OpenFAST shared-library paths resolved through
`OWENSOpenFAST_jll`.
"""
function openfastLibraryPaths()
    return (
        aerodyn_inflow = String(libaerodyn_inflow_c_binding),
        hydrodyn = String(libhydrodyn_c_binding),
        inflowwind = String(libifw_c_binding),
        moordyn = String(libmoordyn_c_binding),
    )
end

function _isOpenFASTArtifactPath(path)
    path_string = String(path)
    return any(artifact_path -> path_string == artifact_path, openfastLibraryPaths())
end

function _shouldCheckOpenFASTLibraryLoad(path)
    return !(Sys.iswindows() && _isOpenFASTArtifactPath(path))
end

function _canLoadOpenFASTLibrary(path)
    handle = Libdl.dlopen_e(path)
    can_load = handle != C_NULL
    can_load && Libdl.dlclose(handle)
    return can_load
end

function _openfast_artifact_status(path)
    path_string = String(path)
    exists = isfile(path_string)
    can_load = exists ?
        (_shouldCheckOpenFASTLibraryLoad(path_string) ?
            _canLoadOpenFASTLibrary(path_string) :
            missing) :
        false
    return (path = path_string, exists = exists, can_load = can_load)
end

function _checkedOpenFASTLibraryPath(label, path)
    path_string = String(path)
    if !isfile(path_string)
        throw(ArgumentError("$label library path does not exist: $path_string"))
    end

    if _shouldCheckOpenFASTLibraryLoad(path_string) &&
        !_canLoadOpenFASTLibrary(path_string)
        throw(ArgumentError("$label library path exists but cannot be loaded: $path_string"))
    end
    return path_string
end

"""
    openfastLibraryArtifactStatus() -> NamedTuple

Return a smoke-test status for each native OpenFAST shared library resolved
through `OWENSOpenFAST_jll`. Each entry includes the absolute artifact path,
whether the file exists, and whether `Libdl` can load it on this platform.
Windows artifact load checks are skipped and reported as `missing` because
preflight `dlopen`/`dlclose` checks can hang on hosted Windows runners; wrapper
initialization still loads the library normally when a module is used.
"""
function openfastLibraryArtifactStatus()
    paths = openfastLibraryPaths()
    return (
        aerodyn_inflow = _openfast_artifact_status(paths.aerodyn_inflow),
        hydrodyn = _openfast_artifact_status(paths.hydrodyn),
        inflowwind = _openfast_artifact_status(paths.inflowwind),
        moordyn = _openfast_artifact_status(paths.moordyn),
    )
end

function endAll()
    if adi_active
        adiEnd()
    end
    if hd_active
        HD_End()
    end
    if ifw_active
        ifwend()
    end
    if md_active
        MD_End()
    end
end

end
