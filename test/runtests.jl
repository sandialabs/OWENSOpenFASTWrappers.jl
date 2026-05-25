using Test

path,_ = splitdir(@__FILE__)

include("pure_helpers_unit.jl")
include("root_motion_unit.jl")
include("hawt_output_parity.jl")
include("vawt_output_parity.jl")

if Sys.iswindows()
    @info "Skipping native OpenFAST integration tests on Windows; helper and parser tests still run."
else
    @testset "InflowWind" begin
        include("ifw_run.jl")
    end

    @testset "HydroDyn" begin
        include("hydrodyn_run.jl")
    end

    @testset "MoorDyn" begin
        include("moordyn_run.jl")
    end

    @testset "AeroDyn" begin
        include("aerodyn_run.jl")
        include("aerodyn_run.jl") #Test running twice
    end
end
