module OpenFASTWrappers
import Libdl

export ifwinit, ifwcalcoutput, ifwend, hdInit, hdCalcOutput, hdUpdateStates, hdEnd, mdInit, mdCalcOutput, mdUpdateStates, mdEnd

include("./inflowwind.jl")
include("./hydrodyn.jl")
include("./moordyn.jl")

end
