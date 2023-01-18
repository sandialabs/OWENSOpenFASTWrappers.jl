import OpenFASTWrappers
using Test
path,_ = splitdir(@__FILE__)

@testset "OpenFASTWrappers.jl inflowwind" begin
    # Usage Example
    # cd openfast
    # mkdir build
    # cd build
    # cmake -DBUILD_SHARED_LIBS=ON ..
    # make ifw_c_binding
    # openfast/build/modules/inflowwind is where the dynamic library will be
    # Move the dynamic library to VAWTHydro.jl/bin or a custom location if you specify the path in the call
    turbsim_filename = "$path/data/ifw/test.bts"
    inflowlib_filename = "$path/../../openfast/build/modules/inflowwind/libifw_c_binding"
    OpenFASTWrappers.ifwinit(;inflowlib_filename,turbsim_filename)
    velocity = OpenFASTWrappers.ifwcalcoutput([0.0,0.0,100.0],0.1)
    println(velocity)
    OpenFASTWrappers.ifwend()
    @test isapprox(velocity,[15.636677, -0.3423329, 0.24170564];atol=1e-4)
end

@testset "HydroDyn" begin
    include("hydrodyn_run.jl")
end

@testset "MoorDyn" begin
    include("moordyn_run.jl")
end
