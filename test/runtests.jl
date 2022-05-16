import OpenFASTWrappers
using Test
path,_ = splitdir(@__FILE__)

@testset "OpenFASTWrappers.jl" begin
    # Usage Example
    # mkdir build
    # cd build
    # cmake ..
    # make ifw_c_binding
    # the the library is located in: deps/bin
    turbsim_filename = "$path/data/ifw/test.bts"
    OpenFASTWrappers.ifwinit(;turbsim_filename)
    velocity = OpenFASTWrappers.ifwcalcoutput([0.0,0.0,100.0],0.1)
    println(velocity)
    OpenFASTWrappers.ifwend()
    @test isapprox(velocity,[15.636677, -0.3423329, 0.24170564];atol=1e-4)
end
