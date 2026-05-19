using LinearAlgebra
using OWENSOpenFASTWrappers

@testset "exported API is defined" begin
    exported_names = names(OWENSOpenFASTWrappers)
    @test all(name -> isdefined(OWENSOpenFASTWrappers, name), exported_names)
    @test :SolvePtfmAccels ∉ exported_names
    @test :SolvePtfmLoads ∉ exported_names
end

@testset "pure AeroDyn helper functions" begin
    @test OWENSOpenFASTWrappers.getAD15numMeshNodes([1 3; 5 2]) == 7

    rot_z_90 = [0.0 1.0 0.0; -1.0 0.0 0.0; 0.0 0.0 1.0]
    rot_x_90 = [1.0 0.0 0.0; 0.0 0.0 1.0; 0.0 -1.0 0.0]
    @test OWENSOpenFASTWrappers.createSingleRotationDCM(90.0, 3) ≈ rot_z_90 atol=eps()
    @test OWENSOpenFASTWrappers.createSingleRotationDCM(90.0, 1) ≈ rot_x_90 atol=eps()
    @test OWENSOpenFASTWrappers.createGeneralTransformationMatrix([90.0, 90.0], [3, 1]) ≈ rot_z_90 * rot_x_90 atol=eps()
    @test OWENSOpenFASTWrappers.calcHubRotMat(zeros(3), pi / 2; rot_axis=3) ≈ rot_z_90 atol=eps()
    @test OWENSOpenFASTWrappers.calcHubRotMat(zeros(3), pi / 2; rot_axis=1) ≈ rot_x_90 atol=eps()

    @test OWENSOpenFASTWrappers.normalizeADIRotationDirection(:ccw) === :ccw
    @test OWENSOpenFASTWrappers.normalizeADIRotationDirection("counter-clockwise") === :ccw
    @test OWENSOpenFASTWrappers.normalizeADIRotationDirection(:clockwise) === :cw
    @test OWENSOpenFASTWrappers.normalizeADIRotationDirection(1) === :ccw
    @test OWENSOpenFASTWrappers.normalizeADIRotationDirection(-1.0) === :cw
    @test OWENSOpenFASTWrappers.rotationDirectionSign(:ccw) == 1
    @test OWENSOpenFASTWrappers.rotationDirectionSign(:cw) == -1
    @test OWENSOpenFASTWrappers.validateADIRotationDirection(false, :ccw) === :ccw
    @test OWENSOpenFASTWrappers.validateADIRotationDirection(true, :cw) === :cw

    @test OWENSOpenFASTWrappers.normalizeOpenFASTInputSource(:file) === :file
    @test OWENSOpenFASTWrappers.normalizeOpenFASTInputSource("path") === :file
    @test OWENSOpenFASTWrappers.normalizeOpenFASTInputSource(:direct) === :text
    @test OWENSOpenFASTWrappers.normalizeOpenFASTInputSource("lines") === :text
    @test OWENSOpenFASTWrappers.openfastInputString(["A", "", "B"]; source=:text) == "A\0\0B"
    @test OWENSOpenFASTWrappers.openfastInputString("A\n\nB\n"; source=:text) == "A\0\0B"
    @test OWENSOpenFASTWrappers.openfastInputString("A\r\nB\r\n"; source=:text) == "A\0B"
    @test OWENSOpenFASTWrappers.openfastInputString("A\0B"; source=:text) == "A\0B"
    mktempdir() do dir
        input_file = joinpath(dir, "openfast_input.dat")
        write(input_file, "A\n\nB\n")
        @test OWENSOpenFASTWrappers.openfastInputString(input_file; source=:file, label="Unit") == "A\0\0B"
    end

    @test_throws ErrorException OWENSOpenFASTWrappers.createSingleRotationDCM(90.0, 4)
    @test_throws ArgumentError OWENSOpenFASTWrappers.normalizeADIRotationDirection(:upwind)
    @test_throws ArgumentError OWENSOpenFASTWrappers.normalizeADIRotationDirection(0)
    @test_throws ArgumentError OWENSOpenFASTWrappers.validateADIRotationDirection(false, :cw)
    @test_throws ArgumentError OWENSOpenFASTWrappers.validateADIRotationDirection(true, :ccw)
    @test_throws ErrorException OWENSOpenFASTWrappers.openfastInputString("none"; source=:file)
    @test_throws ArgumentError OWENSOpenFASTWrappers.normalizeOpenFASTInputSource(:invalid)
    @test_throws ArgumentError OWENSOpenFASTWrappers.openfastInputString(1; source=:text)
end

@testset "pure input-file writers" begin
    mktempdir() do dir
        blade_file = joinpath(dir, "blade.dat")
        blade_text = OWENSOpenFASTWrappers.writeADbladeFile(
            blade_file;
            NumBlNds=2,
            BlSpn=[0.0, 2.0],
            BlChord=[0.25, 0.5],
            BlAFID=[1, 3],
        )
        @test read(blade_file, String) == blade_text
        @test occursin("2   NumBlNds", blade_text)
        @test occursin("2.0 0.0 0.0 0.0 0.0 0.5 3", blade_text)

        olaf_file = joinpath(dir, "OLAF.dat")
        olaf_text = OWENSOpenFASTWrappers.writeOLAFfile(olaf_file; nNWPanel=11, nFWPanels=4)
        @test read(olaf_file, String) == olaf_text
        @test occursin("11           nNWPanel", olaf_text)
        @test occursin("4      nFWPanels", olaf_text)

        ifw_file = joinpath(dir, "IW.dat")
        ifw_text = OWENSOpenFASTWrappers.writeIWfile(
            12.3456789,
            ifw_file;
            RefHt=42.0,
            RefLength=3.0,
            WindType=2,
            windINPfilename="wind.dat",
        )
        @test read(ifw_file, String) == ifw_text
        @test occursin("12.345679   HWindSpeed", ifw_text)
        @test occursin("42.0   RefHt", ifw_text)
        @test occursin("30.0   RefLength", ifw_text)
        @test occursin("\"wind.dat\"", ifw_text)
    end
end
