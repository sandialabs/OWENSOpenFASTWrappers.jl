module OWENSOpenFASTWrappers
import Libdl
using LinearAlgebra     # for aerodyn cross-product

const path = splitdir(@__FILE__)[1]
OFWpath = path

# InflowWind
export ifwinit, ifwcalcoutput, ifwend, OFWpath

# AD15
export Turbine, Environment, Structure, adiInit, adiCalcOutput, adiUpdateStates, adiEnd

# HydroDyn routines
export HD_Init, HD_CalcOutput, HD_UpdateStates, HD_End

# MoorDyn routines
export MD_Init, MD_CalcOutput, MD_UpdateStates, MD_End

# Platform Point Mesh routines
export SolvePtfmAccels, SolvePtfmLoads

include("./aerodyn.jl")
include("./hydrodyn.jl")
include("./inflowwind.jl")
include("./moordyn.jl")

end
