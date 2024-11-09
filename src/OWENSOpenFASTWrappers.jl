module OWENSOpenFASTWrappers
import Libdl
using LinearAlgebra     # for aerodyn cross-product

using OWENSOpenFAST_jll: libaerodyn_inflow_c_binding, libhydrodyn_c_binding,
    libifw_c_binding, libmoordyn_c_binding, turbsim, aerodyn_driver, hydrodyn_driver, inflowwind_driver, moordyn_driver
const path = splitdir(@__FILE__)[1]
OFWpath = path

# InflowWind
export ifwinit, ifwcalcoutput, ifwend, OFWpath, turbsim, inflowwind_driver

# AD15
export Turbine, Environment, Structure, adiInit, adiCalcOutput, adiUpdateStates, adiEnd, aerodyn_driver

# HydroDyn routines
export HD_Init, HD_CalcOutput, HD_UpdateStates, HD_End, hydrodyn_driver

# MoorDyn routines
export MD_Init, MD_CalcOutput, MD_UpdateStates, MD_End, moordyn_driver

# Platform Point Mesh routines
export SolvePtfmAccels, SolvePtfmLoads

include("./aerodyn.jl")
include("./hydrodyn.jl")
include("./inflowwind.jl")
include("./moordyn.jl")

end
