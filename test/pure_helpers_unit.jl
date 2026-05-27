using LinearAlgebra
using OWENSOpenFASTWrappers

@testset "exported API is defined" begin
    exported_names = names(OWENSOpenFASTWrappers)
    @test all(name -> isdefined(OWENSOpenFASTWrappers, name), exported_names)
    @test :SolvePtfmAccels ∉ exported_names
    @test :SolvePtfmLoads ∉ exported_names
end

@testset "OpenFAST artifact smoke checks" begin
    paths = OWENSOpenFASTWrappers.openfastLibraryPaths()
    status = OWENSOpenFASTWrappers.openfastLibraryArtifactStatus()
    expected_keys = (:aerodyn_inflow, :hydrodyn, :inflowwind, :moordyn)

    @test keys(paths) == expected_keys
    @test keys(status) == expected_keys
    for key in expected_keys
        @test getproperty(paths, key) == getproperty(status, key).path
        @test isabspath(getproperty(status, key).path)
        @test getproperty(status, key).exists === true
        if Sys.iswindows()
            @test ismissing(getproperty(status, key).can_load)
        else
            @test getproperty(status, key).can_load === true
            @test OWENSOpenFASTWrappers._checkedOpenFASTLibraryPath(
                String(key),
                getproperty(paths, key),
            ) == getproperty(paths, key)
        end
    end

    missing_path = joinpath(tempdir(), "missing_openfast_override_library")
    @test !isfile(missing_path)
    missing_error = try
        OWENSOpenFASTWrappers._checkedOpenFASTLibraryPath("UnitTest", missing_path)
        nothing
    catch err
        err
    end
    @test missing_error isa ArgumentError
    @test occursin("UnitTest", sprint(showerror, missing_error))
    @test occursin(missing_path, sprint(showerror, missing_error))

    mktemp() do path, io
        write(io, "not a shared library")
        close(io)
        load_error = try
            OWENSOpenFASTWrappers._checkedOpenFASTLibraryPath("UnitTest", path)
            nothing
        catch err
            err
        end
        if Sys.iswindows()
            @test load_error === nothing
        else
            @test load_error isa ArgumentError
            @test occursin("cannot be loaded", sprint(showerror, load_error))
            @test occursin(path, sprint(showerror, load_error))
        end
    end
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

@testset "native wrapper failure guards and error-state resets" begin
    Core.eval(OWENSOpenFASTWrappers, :(adi_active = false))
    Core.eval(OWENSOpenFASTWrappers, :(hd_active = false))
    Core.eval(OWENSOpenFASTWrappers, :(ifw_active = false))
    Core.eval(OWENSOpenFASTWrappers, :(md_active = false))

    Core.eval(OWENSOpenFASTWrappers, :(backup_Vx = 11.25))
    Core.eval(OWENSOpenFASTWrappers, :(ifw_err = IFW_Error([0], "inflow-ok")))
    @test OWENSOpenFASTWrappers.ifwcalcoutput([0.0, 1.0, 2.0], 3.0) ==
          [11.25, 0.0, 0.0]
    OWENSOpenFASTWrappers.ifwend()

    @test_throws ErrorException OWENSOpenFASTWrappers.HD_CalcOutput(
        0.0,
        zeros(6),
        zeros(6),
        zeros(6),
        zeros(Float32, 6),
        zeros(Float32, 1),
    )
    @test_throws ErrorException OWENSOpenFASTWrappers.HD_UpdateStates(
        0.0,
        0.1,
        zeros(6),
        zeros(6),
        zeros(6),
    )
    @test_throws ErrorException OWENSOpenFASTWrappers.MD_CalcOutput(
        0.0,
        zeros(6),
        zeros(6),
        zeros(6),
        zeros(Float32, 6),
        zeros(Float32, 1),
    )
    @test_throws ErrorException OWENSOpenFASTWrappers.MD_UpdateStates(
        0.0,
        0.1,
        zeros(6),
        zeros(6),
        zeros(6),
    )
    OWENSOpenFASTWrappers.HD_End()
    OWENSOpenFASTWrappers.MD_End()
    OWENSOpenFASTWrappers.endAll()

    @test_throws ErrorException OWENSOpenFASTWrappers.adiSetRotorMotion(
        1,
        zeros(3),
        zeros(9),
        zeros(6),
        zeros(6),
        zeros(3),
        zeros(9),
        zeros(6),
        zeros(6),
        zeros(3, 1),
        zeros(9, 1),
        zeros(6, 1),
        zeros(6, 1),
        1,
        zeros(3, 1),
        zeros(9, 1),
        zeros(6, 1),
        zeros(6, 1),
    )
    @test_throws ErrorException OWENSOpenFASTWrappers.adiCalcOutput(0.0, 3)
    Core.eval(OWENSOpenFASTWrappers, :(turbine = [Turbine(1.0, [0.0], zeros(3), 1, 1, 1, [1 1], [1 1], nothing, nothing, true)]))
    @test_throws ErrorException OWENSOpenFASTWrappers.adiGetRotorLoads(1)
    @test_throws ErrorException OWENSOpenFASTWrappers.adiUpdateStates(0.0, 0.1)
    @test_throws ErrorException OWENSOpenFASTWrappers.setRotorMotion(1)

    Core.eval(OWENSOpenFASTWrappers, :(hd_abort_error_level = 4))
    Core.eval(OWENSOpenFASTWrappers, :(hd_err = HD_Error([0], "hydro-ok")))
    OWENSOpenFASTWrappers.hd_check_error()
    @test OWENSOpenFASTWrappers.hd_err.error_status == [0]
    @test length(OWENSOpenFASTWrappers.hd_err.error_message) == 1025

    Core.eval(OWENSOpenFASTWrappers, :(hd_err = HD_Error([2], "hydro-warning")))
    @test_logs (:warn, r"Error status 2: hydro-warning") OWENSOpenFASTWrappers.hd_check_error()
    @test OWENSOpenFASTWrappers.hd_err.error_status == [0]

    hydro_fatal = Ref{Any}(nothing)
    Core.eval(OWENSOpenFASTWrappers, :(hd_err = HD_Error([4], "hydro-fatal")))
    @test_logs (:warn, r"Error status 4: hydro-fatal") begin
        hydro_fatal[] = try
            OWENSOpenFASTWrappers.hd_check_error()
            nothing
        catch err
            err
        end
    end
    @test hydro_fatal[] isa ErrorException
    @test occursin("HydroDyn terminated prematurely", sprint(showerror, hydro_fatal[]))

    Core.eval(OWENSOpenFASTWrappers, :(ifw_abort_error_level = 4))
    Core.eval(OWENSOpenFASTWrappers, :(ifw_err = IFW_Error([0], "inflow-ok")))
    OWENSOpenFASTWrappers.ifw_check_error()
    @test OWENSOpenFASTWrappers.ifw_err.error_status == [0]
    @test length(OWENSOpenFASTWrappers.ifw_err.error_message) == 1025

    Core.eval(OWENSOpenFASTWrappers, :(ifw_err = IFW_Error([2], "inflow-warning")))
    @test_logs (:warn, r"Error status 2: inflow-warning") OWENSOpenFASTWrappers.ifw_check_error()
    @test OWENSOpenFASTWrappers.ifw_err.error_status == [0]

    inflow_fatal = Ref{Any}(nothing)
    Core.eval(OWENSOpenFASTWrappers, :(ifw_err = IFW_Error([4], "inflow-fatal")))
    @test_logs (:warn, r"Error status 4: inflow-fatal") begin
        inflow_fatal[] = try
            OWENSOpenFASTWrappers.ifw_check_error()
            nothing
        catch err
            err
        end
    end
    @test inflow_fatal[] isa ErrorException
    @test occursin("InflowWind terminated prematurely", sprint(showerror, inflow_fatal[]))

    Core.eval(OWENSOpenFASTWrappers, :(md_abort_error_level = 4))
    Core.eval(OWENSOpenFASTWrappers, :(md_err = MD_Error([0], "moor-ok")))
    OWENSOpenFASTWrappers.md_check_error()
    @test OWENSOpenFASTWrappers.md_err.error_status == [0]
    @test length(OWENSOpenFASTWrappers.md_err.error_message) == 1025

    Core.eval(OWENSOpenFASTWrappers, :(md_err = MD_Error([2], "moor-warning")))
    @test_logs (:warn, r"Error status 2: moor-warning") OWENSOpenFASTWrappers.md_check_error()
    @test OWENSOpenFASTWrappers.md_err.error_status == [0]

    moor_fatal = Ref{Any}(nothing)
    Core.eval(OWENSOpenFASTWrappers, :(md_err = MD_Error([4], "moor-fatal")))
    @test_logs (:warn, r"Error status 4: moor-fatal") begin
        moor_fatal[] = try
            OWENSOpenFASTWrappers.md_check_error()
            nothing
        catch err
            err
        end
    end
    @test moor_fatal[] isa ErrorException
    @test occursin("MoorDyn terminated prematurely", sprint(showerror, moor_fatal[]))

    Core.eval(OWENSOpenFASTWrappers, :(adi_abort_error_level = 4))
    Core.eval(OWENSOpenFASTWrappers, :(adi_err = adiError([0], "aero-ok")))
    OWENSOpenFASTWrappers.adi_check_error()
    @test OWENSOpenFASTWrappers.adi_err.error_status == [0]
    @test length(OWENSOpenFASTWrappers.adi_err.error_message) == 1025

    Core.eval(OWENSOpenFASTWrappers, :(adi_err = adiError([2], "aero-warning")))
    @test_logs (:warn, r"Error status 2: aero-warning") OWENSOpenFASTWrappers.adi_check_error()
    @test OWENSOpenFASTWrappers.adi_err.error_status == [0]

    aero_fatal = Ref{Any}(nothing)
    Core.eval(OWENSOpenFASTWrappers, :(adi_err = adiError([4], "aero-fatal")))
    @test_logs (:error, r"Error status 4: aero-fatal") begin
        aero_fatal[] = try
            OWENSOpenFASTWrappers.adi_check_error()
            nothing
        catch err
            err
        end
    end
    @test aero_fatal[] isa ErrorException
    @test occursin("AeroDyn-Inflow terminated prematurely", sprint(showerror, aero_fatal[]))

    mktempdir() do dir
        hydro_file = joinpath(dir, "hydrodyn.dat")
        write(hydro_file, "header\nline\n")
        @test_throws SystemError OWENSOpenFASTWrappers.HD_Init(
            output_root_name=joinpath(dir, "HD"),
            hd_input_file=hydro_file,
            ss_input_file=joinpath(dir, "missing_seastate.dat"),
        )

        hydro_direct_error = Ref{Any}(nothing)
        @test_logs (:error, r"Hydrodyn file") (:error, r"Hydrodyn file") begin
            hydro_direct_error[] = try
                OWENSOpenFASTWrappers.HD_Init(output_root_name=joinpath(dir, "HD-direct"))
                nothing
            catch err
                err
            end
        end
        @test hydro_direct_error[] isa UndefVarError

        uniform_file = joinpath(dir, "uniform.wnd")
        write(uniform_file, "! header\n\n")
        @test_throws BoundsError OWENSOpenFASTWrappers.ifwinit(
            turbsim_filename=uniform_file,
        )
        @test_throws BoundsError OWENSOpenFASTWrappers.ifwinit(
            inflowlib_filename=nothing,
            turbsim_filename=uniform_file,
        )

        @test_throws SystemError OWENSOpenFASTWrappers.MD_Init(
            md_input_file=joinpath(dir, "missing_moordyn.dat"),
        )
    end
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

        ad_input_file = joinpath(dir, "AD input.dat")
        spaced_blade_file = joinpath(dir, "blade with space.dat")
        spaced_airfoil_file = joinpath(dir, "airfoil with space.dat")
        spaced_olaf_file = joinpath(dir, "OLAF with space.dat")
        ad_input_text = OWENSOpenFASTWrappers.writeADinputFile(
            ad_input_file,
            [spaced_blade_file],
            [spaced_airfoil_file],
            spaced_olaf_file,
        )
        @test read(ad_input_file, String) == ad_input_text
        @test occursin("\"$spaced_blade_file\" ADBlFile(1)", ad_input_text)
        @test occursin("\"$spaced_airfoil_file\"    AFNames", ad_input_text)
        @test occursin("\"$spaced_olaf_file\" OLAFInputFileName", ad_input_text)

        @test OWENSOpenFASTWrappers.quotedOpenFASTString("plain.dat", "unit") ==
              "\"plain.dat\""
        quote_error = try
            OWENSOpenFASTWrappers.quotedOpenFASTString("bad\"name.dat", "unit filename")
            nothing
        catch err
            err
        end
        @test quote_error isa ArgumentError
        @test occursin(
            "unit filename cannot contain double quotes: bad\"name.dat",
            sprint(showerror, quote_error),
        )

        single_airfoil_file = joinpath(dir, "single airfoil.dat")
        two_blade_input_file = joinpath(dir, "AD two blades.dat")
        two_blade_text = OWENSOpenFASTWrappers.writeADinputFile(
            two_blade_input_file,
            [spaced_blade_file, "second blade.dat"],
            single_airfoil_file,
            spaced_olaf_file,
        )
        @test read(two_blade_input_file, String) == two_blade_text
        @test occursin("1                      NumAFfiles", two_blade_text)
        @test occursin("\"$single_airfoil_file\"    AFNames", two_blade_text)
        @test occursin("\"$spaced_blade_file\" ADBlFile(1)", two_blade_text)
        @test occursin("\"second blade.dat\" ADBlFile(2)", two_blade_text)

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
