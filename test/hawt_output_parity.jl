using Test
using OWENSOpenFASTWrappers

@testset "OpenFAST output table parser" begin
    mktempdir() do dir
        output_file = joinpath(dir, "driver.out")
        write(
            output_file,
            "preamble\0text\n\n Time  RtAeroCp RtAeroCt\n (s) (-) (-)\n 0.0 1.0 2.0\n bad row\n 0.1 1.5 2.5\n",
        )

        table = OWENSOpenFASTWrappers.readOpenFASTOutputTable(output_file)
        @test table.channels == ["Time", "RtAeroCp", "RtAeroCt"]
        @test table.units == ["(s)", "(-)", "(-)"]
        @test table.data == [0.0 1.0 2.0; 0.1 1.5 2.5]

        candidate_file = joinpath(dir, "candidate.out")
        write(
            candidate_file,
            " Time  RtAeroCp RtAeroCt\n (s) (-) (-)\n 0.0 1.1 1.8\n 0.1 1.4 2.7\n",
        )
        metrics = OWENSOpenFASTWrappers.openfastOutputChannelMetrics(
            output_file,
            candidate_file,
            ["RtAeroCp", "RtAeroCt"],
        )
        @test metrics["RtAeroCp"].n == 2
        @test metrics["RtAeroCp"].rmse ≈ 0.1 atol=1e-15
        @test metrics["RtAeroCp"].mean_bias == 0.0
        @test metrics["RtAeroCp"].mean_abs_error ≈ 0.1 atol=1e-15
        @test metrics["RtAeroCp"].max_abs_error ≈ 0.1 atol=1e-15
        @test metrics["RtAeroCp"].reference_final == 1.5
        @test metrics["RtAeroCp"].candidate_final == 1.4
        @test metrics["RtAeroCp"].final_error ≈ -0.1 atol=1e-15
        @test metrics["RtAeroCt"].mean_bias ≈ 0.0 atol=1e-15
        @test_throws ArgumentError OWENSOpenFASTWrappers.openfastOutputChannelMetrics(
            output_file,
            candidate_file,
            ["missing"],
        )
        @test_throws ArgumentError OWENSOpenFASTWrappers.openfastOutputChannelMetrics(
            output_file,
            candidate_file,
            ["RtAeroCp"];
            rows = 3:3,
        )
    end
end

@testset "HAWT AeroDyn wrapper output parity" begin
    hawt_dir = joinpath(@__DIR__, "data")
    standalone = OWENSOpenFASTWrappers.readOpenFASTOutputTable(
        joinpath(hawt_dir, "hawt_standalone_selected.dat"),
    )
    wrapped = OWENSOpenFASTWrappers.readOpenFASTOutputTable(
        joinpath(hawt_dir, "hawt_wrapper_selected.dat"),
    )
    @test size(standalone.data, 1) == 50
    @test size(wrapped.data, 1) == 51
    @test "RtAeroCp" in standalone.channels
    @test "AB2N005Fxi" in standalone.channels

    channels = [
        "RtAeroCp",
        "RtAeroCq",
        "RtAeroCt",
        "RtAeroPwr",
        "RtAeroFxh",
        "RtAeroMxh",
        "AB2N005Fxi",
        "AB2N005Fyi",
        "AB2N005Fzi",
    ]
    common_rows_after_initialization = 2:min(size(standalone.data, 1), size(wrapped.data, 1))
    metrics = OWENSOpenFASTWrappers.openfastOutputChannelMetrics(
        standalone,
        wrapped,
        channels;
        rows = common_rows_after_initialization,
    )

    @test metrics["RtAeroCp"].rmse ≈ 1.2317235314925243e-9 atol=1e-18
    @test metrics["RtAeroCq"].rmse ≈ 4.134474278213488e-10 atol=1e-18
    @test metrics["RtAeroCt"].rmse ≈ 2.9721496376610603e-9 atol=1e-18
    @test metrics["RtAeroPwr"].rmse ≈ 0.002344076338800445 atol=1e-15
    @test metrics["RtAeroFxh"].rmse ≈ 0.0013037136908755046 atol=1e-15
    @test metrics["RtAeroMxh"].rmse ≈ 0.0034123425429810343 atol=1e-15
    @test metrics["AB2N005Fxi"].rmse ≈ 3.621803586801804e-5 atol=1e-18
    @test metrics["AB2N005Fyi"].rmse ≈ 3.0783158659956865e-6 atol=1e-18
    @test metrics["AB2N005Fzi"].rmse ≈ 7.621960752142232e-6 atol=1e-18

    @test metrics["RtAeroCp"].max_abs_error ≈ 3.800000000428461e-9 atol=1e-18
    @test metrics["RtAeroCt"].max_abs_error ≈ 7.89999999965818e-9 atol=1e-18
    @test metrics["RtAeroPwr"].max_abs_error ≈ 0.007400000002235174 atol=1e-15
    @test metrics["RtAeroFxh"].max_abs_error ≈ 0.004149999999754073 atol=1e-15
    @test metrics["RtAeroMxh"].max_abs_error ≈ 0.00999999999476131 atol=1e-15
    @test metrics["AB2N005Fxi"].max_abs_error ≈ 6.369999999833453e-5 atol=1e-18
    @test metrics["AB2N005Fyi"].max_abs_error ≈ 6.700000000137152e-6 atol=1e-18
    @test metrics["AB2N005Fzi"].max_abs_error ≈ 1.3299999999105694e-5 atol=1e-18
    @test metrics["RtAeroCp"].final_error ≈ 1.6000000005456094e-9 atol=1e-18

    initial_metrics = OWENSOpenFASTWrappers.openfastOutputChannelMetrics(
        standalone,
        wrapped,
        ["RtAeroCt"];
        rows = 1:1,
    )
    @test initial_metrics["RtAeroCt"].max_abs_error ≈ 0.4756504058 atol=1e-12
end
