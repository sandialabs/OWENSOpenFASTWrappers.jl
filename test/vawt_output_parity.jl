using Test
using OWENSOpenFASTWrappers

@testset "VAWT AeroDyn wrapper output parity" begin
    data_dir = joinpath(@__DIR__, "data")
    standalone = OWENSOpenFASTWrappers.readOpenFASTOutputTable(
        joinpath(data_dir, "vawt_standalone_selected.dat"),
    )
    wrapped = OWENSOpenFASTWrappers.readOpenFASTOutputTable(
        joinpath(data_dir, "vawt_wrapper_selected.dat"),
    )

    @test size(standalone.data, 1) == 49
    @test size(wrapped.data, 1) == 49
    @test standalone.channels == wrapped.channels
    @test standalone.channels == [
        "Time",
        "RtAeroCp",
        "RtAeroCq",
        "RtAeroCt",
        "RtAeroPwr",
        "RtAeroFxh",
        "RtAeroFyh",
        "RtAeroFzh",
        "RtAeroMxh",
        "RtAeroMyh",
        "RtAeroMzh",
        "AB4N005STVx",
        "AB4N005STVy",
        "AB4N005STVz",
        "AB4N005Vx",
        "AB4N005Vy",
        "AB4N005Alpha",
        "AB4N005Fx",
        "AB4N005Fy",
    ]
    @test standalone.data[:, 1] == wrapped.data[:, 1]
    @test first(standalone.data[:, 1]) == 0.0
    @test last(standalone.data[:, 1]) == 0.48

    channels = standalone.channels[2:end]
    metrics = OWENSOpenFASTWrappers.openfastOutputChannelMetrics(
        standalone,
        wrapped,
        channels,
    )

    @test metrics["RtAeroCp"].rmse == 0.0
    @test metrics["RtAeroCq"].rmse == 0.0
    @test metrics["RtAeroCt"].rmse == 0.0
    @test metrics["RtAeroPwr"].rmse ≈ 3905.182583712624 atol=1e-9
    @test metrics["RtAeroFxh"].rmse ≈ 17.703266412836808 atol=1e-12
    @test metrics["RtAeroFyh"].rmse ≈ 300.7481398090607 atol=1e-10
    @test metrics["RtAeroFzh"].rmse ≈ 521.4980006673471 atol=1e-10
    @test metrics["RtAeroMxh"].rmse ≈ 745.8349568796287 atol=1e-10
    @test metrics["RtAeroMyh"].rmse ≈ 70.79885377960191 atol=1e-11
    @test metrics["RtAeroMzh"].rmse ≈ 45.98426174479197 atol=1e-11
    @test metrics["AB4N005STVx"].rmse ≈ 1.2008326033201033e-7 atol=1e-18
    @test metrics["AB4N005STVy"].rmse ≈ 0.010686340826312493 atol=1e-15
    @test metrics["AB4N005STVz"].rmse ≈ 2.0704010388394978e-7 atol=1e-18
    @test metrics["AB4N005Vx"].rmse ≈ 0.40232367528038904 atol=1e-14
    @test metrics["AB4N005Vy"].rmse ≈ 0.665462095145674 atol=1e-14
    @test metrics["AB4N005Alpha"].rmse ≈ 1.1384450273890505 atol=1e-14
    @test metrics["AB4N005Fx"].rmse ≈ 11.395820793687411 atol=1e-12
    @test metrics["AB4N005Fy"].rmse ≈ 1.7278121550093446 atol=1e-13

    @test metrics["RtAeroPwr"].max_abs_error ≈ 7363.158199999998 atol=1e-9
    @test metrics["RtAeroFyh"].max_abs_error ≈ 719.3679700000002 atol=1e-10
    @test metrics["RtAeroFzh"].max_abs_error ≈ 833.212751 atol=1e-10
    @test metrics["RtAeroMxh"].max_abs_error ≈ 1406.25972 atol=1e-10
    @test metrics["AB4N005Vx"].max_abs_error ≈ 0.5244381819999984 atol=1e-15
    @test metrics["AB4N005Alpha"].max_abs_error ≈ 1.881326099999999 atol=1e-14
    @test metrics["AB4N005Fx"].max_abs_error ≈ 17.514237099999995 atol=1e-12

    @test metrics["RtAeroPwr"].final_error ≈ 4579.824699999997 atol=1e-9
    @test metrics["RtAeroFxh"].final_error ≈ -21.245927000000023 atol=1e-12
    @test metrics["RtAeroFyh"].final_error ≈ -148.24448000000007 atol=1e-10
    @test metrics["AB4N005Alpha"].final_error ≈ -1.4196324200000001 atol=1e-14
    @test metrics["AB4N005Fy"].final_error ≈ 3.5362967899999997 atol=1e-13
end
