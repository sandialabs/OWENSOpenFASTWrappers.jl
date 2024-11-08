import OWENSOpenFASTWrappers
using Test
# Usage Example
# cd openfast
# mkdir build
# cd build
# cmake -DBUILD_SHARED_LIBS=ON ..
# make ifw_c_binding
# openfast/build/modules/inflowwind is where the dynamic library will be
# Move the dynamic library to VAWTHydro.jl/bin or a custom location if you specify the path in the call
turbsim_filename = "$path/data/ifw/test.bts"
# inflowlib_filename = "$path/../deps/openfast/build/modules/inflowwind/libifw_c_binding"
OWENSOpenFASTWrappers.ifwinit(; turbsim_filename)
velocity = OWENSOpenFASTWrappers.ifwcalcoutput([0.0,0.0,100.0],0.1)
println(velocity)
OWENSOpenFASTWrappers.ifwend()
@test isapprox(velocity,[15.636677, -0.3423329, 0.24170564];atol=1e-4)
