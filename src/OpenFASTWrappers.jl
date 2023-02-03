module OpenFASTWrappers
import Libdl

export ifwinit, ifwcalcoutput, ifwend

# HydroDyn routines
export HD_Init
export HD_CalcOutput
export HD_UpdateStates
export HD_End

# MoorDyn routines
export MD_Init
export MD_CalcOutput
export MD_UpdateStates
export MD_End

# Platform Point Mesh routines
export SolvePtfmAccels
export SolvePtfmLoads

include("./inflowwind.jl")
include("./hydrodyn.jl")
include("./moordyn.jl")

end
