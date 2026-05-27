using OWENSOpenFASTWrappers
using LinearAlgebra: I

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

@testset "HAWT mesh position, velocity, and orientation kinematics" begin
    mesh = (
        x = [1.0, 2.0, 3.0, 4.0],
        y = [10.0, 20.0, 30.0, 40.0],
        z = [100.0, 200.0, 300.0, 400.0],
    )
    ort = (
        Psi_d = zeros(4),
        Theta_d = zeros(4),
        Twist_d = zeros(4),
    )
    blade_idx = [1 2; 4 3]
    turbine = OWENSOpenFASTWrappers.Turbine(
        50.0,
        [2.0],
        zeros(3),
        1,
        2,
        4,
        blade_idx,
        blade_idx,
        mesh,
        ort,
        true,
    )

    u_j = zeros(24)
    udot_j = zeros(24)
    uddot_j = zeros(24)
    for inode in 1:4
        offset = (inode - 1) * 6
        u_j[offset+1:offset+3] = [0.1inode, 0.2inode, -0.3inode]
        udot_j[offset+1:offset+6] = [inode, 10inode, 100inode, 0.1inode, 0.2inode, 0.3inode]
        uddot_j[offset+1:offset+6] = [-inode, -10inode, -100inode, -0.1inode, -0.2inode, -0.3inode]
    end

    nac_pos = [5.0, -2.0, 1.0]
    hub_pos = zeros(3)
    hub_angle = zeros(3)

    root_pos = OWENSOpenFASTWrappers.getRootPos(
        turbine,
        u_j,
        0.0,
        nac_pos,
        hub_pos,
        hub_angle,
    )
    @test root_pos isa Matrix{Float32}
    @test root_pos ≈ Float32[
        6.1 9.4
        8.2 38.8
        100.7 399.8
    ] atol = 1f-5

    mesh_pos = OWENSOpenFASTWrappers.getAD15MeshPos(
        turbine,
        u_j,
        0.0,
        nac_pos,
        hub_pos,
        hub_angle,
    )
    @test mesh_pos isa Matrix{Float32}
    @test mesh_pos ≈ Float32[
        104.7 204.4 403.8 304.1
        8.2 18.4 38.8 28.6
        -0.1 -1.2 -3.4 -2.3
    ] atol = 1f-5

    blade_mesh_pos = Float32[
        0.0 0.0 0.0 0.0
        1.0 2.0 4.0 3.0
        0.0 0.0 -1.0 1.0
    ]
    hub_vel = [0.5, -0.25, 1.0, 9.0, 9.0, 9.0]
    hub_acc = [-0.5, 0.25, -1.0, 8.0, 8.0, 8.0]
    mesh_vel, mesh_acc = OWENSOpenFASTWrappers.getAD15MeshVelAcc(
        turbine,
        blade_mesh_pos,
        udot_j,
        uddot_j,
        0.0,
        2.0,
        3.0,
        zeros(3),
        hub_pos,
        hub_angle,
        hub_vel,
        hub_acc,
    )
    @test mesh_vel isa Matrix{Float32}
    @test mesh_acc isa Matrix{Float32}
    @test mesh_vel ≈ Float32[
        1.5 2.5 4.5 3.5
        9.75 19.75 41.75 27.75
        103.0 205.0 409.0 307.0
        2.1 2.2 2.4 2.3
        0.2 0.4 0.8 0.6
        0.3 0.6 1.2 0.9
    ] atol = 1f-6
    @test mesh_acc ≈ Float32[
        -1.5 -2.5 -4.5 -3.5
        -9.75 -19.75 -36.75 -32.75
        -98.0 -195.0 -389.0 -292.0
        2.9 2.8 2.6 2.7
        -0.2 -0.4 -0.8 -0.6
        -0.3 -0.6 -1.2 -0.9
    ] atol = 1f-6

    expected_hawt_dcm = Float32[
        -1.0,
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        0.0,
        -1.0,
    ]
    root_dcm = OWENSOpenFASTWrappers.getRootDCM(turbine, u_j, 0.0, hub_angle)
    @test root_dcm isa Matrix
    @test root_dcm ≈ repeat(expected_hawt_dcm, 1, 2) atol = 1f-6

    dcm_turbine = OWENSOpenFASTWrappers.Turbine(
        50.0,
        [2.0],
        zeros(3),
        1,
        2,
        6,
        blade_idx,
        blade_idx,
        mesh,
        ort,
        true,
    )
    mesh_dcm = OWENSOpenFASTWrappers.getAD15MeshDCM(dcm_turbine, u_j, 0.0, hub_angle)
    @test mesh_dcm isa Matrix
    @test mesh_dcm ≈ repeat(expected_hawt_dcm, 1, 6) atol = 1f-6
end

@testset "deformAD15 updates HAWT motion state before inactive-library guard" begin
    mesh = (
        x = [0.0, 0.0, 0.0],
        y = [1.0, 2.0, 3.0],
        z = [0.0, 0.0, 0.0],
        numNodes = 3,
    )
    ort = (
        Psi_d = zeros(3),
        Theta_d = zeros(3),
        Twist_d = zeros(3),
    )
    turbine = OWENSOpenFASTWrappers.Turbine(
        2.0,
        [0.0],
        zeros(3),
        1,
        1,
        3,
        [1 3],
        [1 2],
        mesh,
        ort,
        true,
    )
    structure = OWENSOpenFASTWrappers.Structure(
        zeros(Float32, 3),
        vec(Matrix{Float64}(I, 3, 3)),
        zeros(Float32, 6),
        zeros(Float32, 6),
        zeros(Float32, 3),
        Matrix{Float64}(I, 3, 3),
        zeros(Float32, 6),
        zeros(Float32, 6),
        zeros(Float32, 3, 1),
        zeros(9, 1),
        zeros(Float32, 6, 1),
        zeros(Float32, 6, 1),
        zeros(Float32, 3, 3),
        zeros(9, 3),
        zeros(Float32, 6, 3),
        zeros(Float32, 6, 3),
    )
    Core.eval(OWENSOpenFASTWrappers, :(adi_active = false))
    Core.eval(OWENSOpenFASTWrappers, :(turbine = [$turbine]))
    Core.eval(OWENSOpenFASTWrappers, :(turbstruct = [$structure]))

    hub_vel = [zeros(Float32, 6)]
    hub_acc = [zeros(Float32, 6)]
    deform_error = try
        OWENSOpenFASTWrappers.deformAD15(
            [zeros(18)],
            [zeros(18)],
            [zeros(18)],
            [0.0],
            [7.0],
            [0.7],
            [zeros(Float32, 3)],
            [zeros(3)],
            hub_vel,
            hub_acc,
        )
        nothing
    catch err
        err
    end
    @test deform_error isa ErrorException
    @test occursin("Could not set rotor motion turbine 1", sprint(showerror, deform_error))
    @test hub_vel[1] == Float32[0.0, 0.0, 0.0, 7.0, 0.0, 0.0]
    @test hub_acc[1] == Float32[0.0, 0.0, 0.0, 0.7, 0.0, 0.0]
    @test OWENSOpenFASTWrappers.turbstruct[1].rootPos ≈ Float32[0.0; 1.0; 0.0;;] atol = 1f-6
    @test OWENSOpenFASTWrappers.turbstruct[1].rootVel[:, 1] ≈
          Float32[0.0, 0.0, 7.0, 7.0, 0.0, 0.0] atol = 1f-6
    @test OWENSOpenFASTWrappers.turbstruct[1].rootAcc[:, 1] ≈
          Float32[0.0, 0.0, 0.7, 0.7, 0.0, 0.0] atol = 1f-6
end
