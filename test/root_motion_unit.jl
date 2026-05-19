using OWENSOpenFASTWrappers

@testset "root hub rotation kinematics" begin
    turbine = OWENSOpenFASTWrappers.Turbine(
        1.0,
        [2.0],
        zeros(3),
        1,
        1,
        1,
        [1 1],
        [1 1],
        nothing,
        nothing,
        false,
    )

    root_pos = reshape(Float32[1.0, 0.0, 0.0], 3, 1)
    root_vel, root_acc = OWENSOpenFASTWrappers.getRootVelAcc(
        turbine,
        root_pos,
        zeros(Float32, 6),
        zeros(Float32, 6),
        0.0,
        2.0,
        3.0,
        zeros(Float32, 3),
        zeros(Float32, 3),
        zeros(Float32, 3),
        zeros(Float32, 6),
        zeros(Float32, 6),
    )

    @test root_vel isa Matrix{Float32}
    @test root_acc isa Matrix{Float32}
    @test size(root_vel) == (6, 1)
    @test size(root_acc) == (6, 1)
    @test root_vel[:, 1] == Float32[0.0, 2.0, 0.0, 0.0, 0.0, 2.0]
    @test root_acc[:, 1] == Float32[0.0, 3.0, 0.0, 0.0, 0.0, 3.0]
end

@testset "HAWT root x-axis rotation kinematics" begin
    turbine = OWENSOpenFASTWrappers.Turbine(
        1.0,
        [2.0],
        zeros(3),
        1,
        1,
        1,
        [1 1],
        [1 1],
        nothing,
        nothing,
        true,
    )

    root_pos = reshape(Float32[0.0, 1.0, 0.0], 3, 1)
    root_vel, root_acc = OWENSOpenFASTWrappers.getRootVelAcc(
        turbine,
        root_pos,
        zeros(Float32, 6),
        zeros(Float32, 6),
        0.0,
        2.0,
        3.0,
        zeros(Float32, 3),
        zeros(Float32, 3),
        zeros(Float32, 3),
        zeros(Float32, 6),
        zeros(Float32, 6),
    )

    @test root_vel isa Matrix{Float32}
    @test root_acc isa Matrix{Float32}
    @test size(root_vel) == (6, 1)
    @test size(root_acc) == (6, 1)
    @test root_vel[:, 1] == Float32[0.0, 0.0, 2.0, 2.0, 0.0, 0.0]
    @test root_acc[:, 1] == Float32[0.0, 0.0, 3.0, 3.0, 0.0, 0.0]
end
