using Test

path,_ = splitdir(@__FILE__)

@testset "InflowWind" begin
    include("ifw_run.jl")
end

# @testset "HydroDyn" begin
#     include("hydrodyn_run.jl")
# end

@testset "MoorDyn" begin
    include("moordyn_run.jl")
end

# @testset "AeroDyn" begin
#     include("aerodyn_run.jl")
# end
