global adilib
global sym_calcoutput
global sym_updatestates
global sym_end

path,_ = splitdir(@__FILE__)

mutable struct adiError
error_status
error_message
end



"""
    adiInit(adilib_filename output_root_name ; )

calls aerodyn_inflow_init to initialize AeroDyn and InflowWind together

# Inputs:
* `adilib_filename::string`: path and name of AeroDyn-Inflow dynamic library

* `ad_input_file_passed::bool`: flag to indicate the AD15 input file is passed as a string
                                (set to false if passing input file name instead, NOT SUPPORTED YET)
* `ad_input_file::string`: name of input file for AD15 -- this is read by julia and passed to AD15
* `ifw_input_file_passed::bool`: flag to indicate the InflowWind input file is passed as a string
                                (set to false if passing input file name instead, NOT SUPPORTED YET)
* `ifw_input_file::string`: name of input file for InflowWind -- this is read by julia and passed to InflowWind

* `gravity::float`:     optional, gravity value (default: 9.80665 m/s^2)
* `defFldDens::float`:  optional, air density (default: 1.225 kg/m^3)
* `defKinVisc::float`:  optional, kinematic viscosity of working fluid (default: 1.464E-05 m^2/s)
* `defSpdSound::float`: optional, speed of sound in working fluid (default: 335.0 m/s)
* `defPatm::float`:     optional, atmospheric pressure (default:  103500.0 Pa) [used only for an MHK turbine cavitation check]
* `defPvap::float`:     optional, vapour pressure of working fluid (default: 1700.0 Pa) [used only for an MHK turbine cavitation check]
* `WtrDpth::float`:     optional, water depth (default: 0.0 m) [used only for an MHK turbine]
* `MSL2SWL::float`:     optional, offset between still-water level and mean sea level (default: 0.0 m) [positive upward, used only for an MHK turbine]
* `AeroProjMod::int`:   optional, aero projection mode:
*                               1.      APM_BEM_NoSweepPitchTwist  "Original AeroDyn model where momentum balance is done in the WithoutSweepPitchTwist system"
*                               2.      APM_BEM_Polar              "Use staggered polar grid for momentum balance in each annulus"
*                               3.      APM_LiftingLine            "Use the blade lifting line (i.e. the structural) orientation (currently for OLAF with VAWT)"

* `storeHHVel::bool`:   optional, internal parameter for adi_library.  Exposed for convenience, but not needed.
* `transposeDCM::bool`: optional, transpose DCM internally in ADI to match calling code convention for direction cosine matrices (default: true)
* `WrVTK::int`:         optional, write VTK output files at all timesteps to visualize AeroDyn 15 meshes [0 none (default), 1 ref, 2 motion]
* `WrVTK_Type::int`:    optional, write VTK output files as [1 surfaces (default), 2 lines, 3 both]
* `VTKNacDim::Array(float*6)`   optional, Nacelle Dimension for VTK visualization x0,y0,z0,Lx,Ly,Lz (m)
* `VTKHubRad::float`:   optional, HubRadius for VTK visualization (m)

* `initHubPos::Array(float)`: required, (x,y,z) position of hub
* `initHubOrient::Array(float)`: required, orientation of hub as 9 element vector of flattened DCM

* `initNacellePos::Array(float)`: required, (x,y,z) position of nacelle
* `initNacelleOrient::Array(float)`: required, orientation of nacelle as 9 element vector of flattened DCM

* `numBlades::int`: required, number of blades
* `initRootPos::Array(float)`: required, size (numBlades,3) position vectors of roots
* `initRootOrient::Array(float)`: required, size (numBlades,9) orientation DCMs flattened to array of 9 element vectors

* `numMeshNodes::int`: required, number of structural mesh points (total across all blades)
* `initMeshPos::Array(float)`: required, size (numMeshNodes,3) position vectors of mesh points
* `initMeshOrient::Array(float)`: required, size (numMeshNodes,9) orientation DCMs flattened to array of 9 element vectors

* `interp_order::int`: optional, interpolation order used internally [1 first order (default), 2 second order]

* `t_initial::float`:   optional, initial time of simulation (default: 0.0)
* `dt::float64`:        required, timestep for AD15 (needed for setting internal constants)
* `t_max::float`:       required, total expected simulation time -- used only for setting VTK counter width

# Outputs:
* `num_channels::int`: number of output channels
* `channel_names::string`: string of output channel names from ADI
* `channel_units::string`: string of output channel units from ADI

"""
function adiInit(adilib_filename, output_root_name;
    ad_input_file_passed= true,
    ad_input_file="none",
    ifw_input_file_passed= true,
    ifw_input_file="none",
    gravity     =   9.80665,  # Gravitational acceleration (m/s^2)
    defFldDens  =     1.225,  # Air density (kg/m^3)
    defKinVisc  = 1.464E-05,  # Kinematic viscosity of working fluid (m^2/s)
    defSpdSound =     335.0,  # Speed of sound in working fluid (m/s)
    defPatm     =  103500.0,  # Atmospheric pressure (Pa) [used only for an MHK turbine cavitation check]
    defPvap     =    1700.0,  # Vapour pressure of working fluid (Pa) [used only for an MHK turbine cavitation check]
    WtrDpth     =       0.0,  # Water depth (m)
    MSL2SWL     =       0.0,  # Offset between still-water level and mean sea level (m) [positive upward]
    AeroProjMod =         1,  # see note
    storeHHVel  = false, 
    transposeDCM= true,       # transpose DCM internally for calculations
    WrVTK       = 0,          # write VTK files from adi [0 none, 1 ref, 2 motion]
    WrVTK_Type  = 1,          # write VTK files from adi [1 surfaces, 2 lines, 3 both]
    VTKNacDim   = [-1 ,-1 ,-1 ,2 ,2 ,2],        # Nacelle Dimension for VTK visualization x0,y0,z0,Lx,Ly,Lz (m)
    VTKHubRad   = 0.1,                          # HubRadius for VTK visualization (m)
    initHubPos         = zeros(3),  # initial position vector of hub
    initHubOrient      = zeros(9),  # initial orientation of hub (flattened 3x3 DCM)
    initNacellePos     = zeros(3),  # initial position vector of nacelle 
    initNacelleOrient  = zeros(9),  # initial orientation of nacelle (flattened 3x3 DCM)
    numBlades          = 3,         # number of blades in system
    initRootPos        = zeros(numBlades,3),    # initial root position vectors
    initRootOrient     = zeros(numBlades,9),    # initial root orientation DCMs
    numMeshNodes       = 1,         # number of mesh points representing structural mesh of rotor
    initMeshPos        = zeros(numMeshNodes,3),   # initial position vectors of all mesh points
    initMeshOrient     = zeros(numMeshNodes,9),   # initial orientations of all mesh points
    interp_order=1,
    t_initial=0.0,
    dt=0.01,
    t_max=60.0)

    global adi_abort_error_level = 4

    # AeroDyn 15 input file
    if ad_input_file_passed == false
        ad_input_string = ad_input_file
    else
        if ad_input_file == "none"
            ad_input_string_array = [
            ]
            error("Default AeroDyn input file not setup yet.")
        else
            println("Reading AeroDyn data from $ad_input_file.")
            fid = open(ad_input_file, "r") 
            ad_input_string_array = readlines(fid)
            ad_input_string        = join(ad_input_string_array, "\0")
            close(fid)
        end
    end
    ad_input_string_length = length(ad_input_string)


    # InflowWind input file
    if ifw_input_file_passed == false
        ifw_input_string = ifw_input_file
    else
        if ifw_input_file == "none"
            ifw_input_string_array = [
            ]
            error("Default InflowWind input file not setup yet.")
        else
            println("Reading InfloWind data from $ifw_input_file.")
            fid = open(ifw_input_file, "r") 
            ifw_input_string_array = readlines(fid)
            ifw_input_string        = join(ifw_input_string_array, "\0")
            close(fid)
        end
    end
    ifw_input_string_length = length(ifw_input_string)


    # Allocate Outputs
    num_channels = [0]
    channel_names = string(repeat(" ", 20 * 8000))      # This must match value for MaxADIOutputs in the library
    channel_units = string(repeat(" ", 20 * 8000))      # This must match value for MaxADIOutputs in the library

    try
        println("Attempting to access AeroDyn-Inflow at: $adilib_filename")
        global adilib = Libdl.dlopen(adilib_filename) # Open the library explicitly.
        global adi_active = true

        adi_sym_init = Libdl.dlsym(adilib, :AeroDyn_Inflow_C_Init)   # Get a symbol for the function to call.
        global adi_sym_calcoutput = Libdl.dlsym(adilib, :AeroDyn_Inflow_C_CalcOutput)   # Get a symbol for the function to call.
        global adi_sym_updatestates = Libdl.dlsym(adilib, :AeroDyn_Inflow_C_UpdateStates)
        global adi_sym_end = Libdl.dlsym(adilib, :AeroDyn_Inflow_C_End) # !!! "c" is capitalized in library, change if errors
        global adi_err = adiError([0], string(repeat(" ", 1025)))
        
        ccall(adi_sym_init,Cint,
            (Ref{Cint},         # IN: ad input file passed as string (c_bool)
            Ptr{Ptr{Cchar}},    # IN: ad_input_string_array
            Ref{Cint},          # IN: ad_input_string_length
            Ref{Cint},          # IN: ifw input file passed as string (c_bool)
            Ptr{Ptr{Cchar}},    # IN: ifw_input_string_array
            Ref{Cint},          # IN: ifw_input_string_length
            Cstring,            # IN: output_root_name
            Ref{Cfloat},        # IN: gravity
            Ref{Cfloat},        # IN: defFldDens
            Ref{Cfloat},        # IN: defKinVisc
            Ref{Cfloat},        # IN: defSpdSound
            Ref{Cfloat},        # IN: defPatm
            Ref{Cfloat},        # IN: defPvap
            Ref{Cfloat},        # IN: WtrDpth
            Ref{Cfloat},        # IN: MSL2SWL
            Ref{Cint},          # IN: AeroProjMod
            Ref{Cint},          # IN: interp_order
            Ref{Cdouble},       # IN: t_initial
            Ref{Cdouble},       # IN: dt
            Ref{Cdouble},       # IN: t_max
            Ref{Cint},          # IN: storeHHVel    (c_bool)
            Ref{Cint},          # IN: transposeDCM  (c_bool)
            Ref{Cint},          # IN: WrVTK
            Ref{Cint},          # IN: WrVTK_Type
            Ref{Cfloat},        # IN: VTKNacDim
            Ref{Cfloat},        # IN: VTKHubRad
            Ref{Cfloat},        # IN: initHubPos
            Ref{Cdouble},       # IN: initHubOrient (do we need to flatten this, or just do fortran index order???)
            Ref{Cfloat},        # IN: initNacellePos
            Ref{Cdouble},       # IN: initNacelleOrient (do we need to flatten this, or just do fortran index order???)
            Ref{Cint},          # IN: numBlades
            Ref{Cfloat},        # IN: initRootPos
            Ref{Cdouble},       # IN: initRootOrient (do we need to flatten this, or just do fortran index order???)
            Ref{Cint},          # IN: numMeshNodes
            Ref{Cfloat},        # IN: initMeshPos
            Ref{Cdouble},       # IN: initMeshOrient (do we need to flatten this, or just do fortran index order???)
            Ptr{Cint},          # OUT: num_channels
            Cstring,            # OUT: channel_names
            Cstring,            # OUT: channel_units
            Ptr{Cint},          # OUT: error_status
            Cstring),           # OUT: error_message
            Cint.(ad_input_file_passed),
            [ad_input_string],
            ad_input_string_length,
            Cint.(ifw_input_file_passed),
            [ifw_input_string],
            ifw_input_string_length,
            output_root_name,
            gravity,
            defFldDens,
            defKinVisc,
            defSpdSound,
            defPatm,
            defPvap,
            WtrDpth,
            MSL2SWL,
            AeroProjMod,
            interp_order,
            t_initial,
            dt,
            t_max,
            Cint.(storeHHVel),
            Cint.(transposeDCM),
            WrVTK,
            WrVTK_Type,
            Cfloat.(VTKNacDim),
            VTKHubRad,
            Cfloat.(initHubPos),
            Cdouble.(initHubOrient),
            Cfloat.(initNacellePos),
            Cdouble.(initNacelleOrient),
            numBlades,
            Cfloat.(initRootPos),
            Cdouble.(initRootOrient),
            numMeshNodes,
            Cfloat.(initMeshPos),
            Cdouble.(initMeshOrient),
            num_channels,
            channel_names,
            channel_units,
            adi_err.error_status,
            adi_err.error_message)
 
        adi_check_error()
    catch
        error("AeroDyn-InflowWind library could not initialize")
        global adi_active = false
    end

    return num_channels, channel_names, channel_units

end

function adiCalcOutput(time, 
                 hubPos, hubOrient, hubVel, hubAcc,
                 nacPos, nacOrient, nacVel, nacAcc,
                 rootPos, rootOrient, rootVel, rootAcc,
                 numMeshNodes,
                 meshPos, meshOrient, meshVel, meshAcc,
                 num_channels)
    meshFrcMom = zeros(Cfloat,6,numMeshNodes);
    out_channel_vals = zeros(Cfloat,1,num_channels)

    if adi_active
        ccall(adi_sym_calcoutput,Cint,
            (Ptr{Cdouble},      # IN: time
            Ref{Cfloat},        # IN: hubPos 
            Ref{Cdouble},       # IN: hubOrient
            Ref{Cfloat},        # IN: hubVel 
            Ref{Cfloat},        # IN: hubAcc 
            Ref{Cfloat},        # IN: nacPos 
            Ref{Cdouble},       # IN: nacOrient
            Ref{Cfloat},        # IN: nacVel 
            Ref{Cfloat},        # IN: nacAcc 
            Ref{Cfloat},        # IN: rootPos 
            Ref{Cdouble},       # IN: rootOrient
            Ref{Cfloat},        # IN: rootVel 
            Ref{Cfloat},        # IN: rootAcc 
            Ref{Cint},          # IN: numMeshNodes
            Ref{Cfloat},        # IN: meshPos 
            Ref{Cdouble},       # IN: meshOrient
            Ref{Cfloat},        # IN: meshVel 
            Ref{Cfloat},        # IN: meshAcc 
            Ptr{Cfloat},        # OUT: node_force
            Ptr{Cfloat},        # OUT: out_channel_vals
            Ptr{Cint},          # OUT: error_status
            Cstring),           # OUT: error_message 
            [time],
            Cfloat.(hubPos),
            Cdouble.(hubOrient),
            Cfloat.(hubVel),
            Cfloat.(hubAcc),
            Cfloat.(nacPos),
            Cdouble.(nacOrient),
            Cfloat.(nacVel),
            Cfloat.(nacAcc),
            Cfloat.(rootPos),
            Cdouble.(rootOrient),
            Cfloat.(rootVel),
            Cfloat.(rootAcc),
            numMeshNodes,
            Cfloat.(meshPos),
            Cdouble.(meshOrient),
            Cfloat.(meshVel),
            Cfloat.(meshAcc),
            meshFrcMom,
            out_channel_vals,
            adi_err.error_status,
            adi_err.error_message) 

        adi_check_error()

    else
        error("AerooDyn-Inflow instance has not been initialized. Use adiInit() function.")
    end

    return meshFrcMom, out_channel_vals
end

function adiUpdateStates(time, next_time,
                 hubPos, hubOrient, hubVel, hubAcc,
                 nacPos, nacOrient, nacVel, nacAcc,
                 rootPos, rootOrient, rootVel, rootAcc,
                 numMeshNodes,
                 meshPos, meshOrient, meshVel, meshAcc)

    if adi_active
        ccall(adi_sym_updatestates,Cint,
            (Ptr{Cdouble},      # IN: time
            Ptr{Cdouble},       # IN: next_time
            Ref{Cfloat},        # IN: hubPos 
            Ref{Cdouble},       # IN: hubOrient
            Ref{Cfloat},        # IN: hubVel 
            Ref{Cfloat},        # IN: hubAcc 
            Ref{Cfloat},        # IN: nacPos 
            Ref{Cdouble},       # IN: nacOrient
            Ref{Cfloat},        # IN: nacVel 
            Ref{Cfloat},        # IN: nacAcc 
            Ref{Cfloat},        # IN: rootPos 
            Ref{Cdouble},       # IN: rootOrient
            Ref{Cfloat},        # IN: rootVel 
            Ref{Cfloat},        # IN: rootAcc 
            Ref{Cint},          # IN: numMeshNodes
            Ref{Cfloat},        # IN: meshPos 
            Ref{Cdouble},       # IN: meshOrient
            Ref{Cfloat},        # IN: meshVel 
            Ptr{Cint},          # OUT: error_status
            Cstring),           # OUT: error_message 
            [time],
            [next_time],
            Cfloat.(hubPos),
            Cdouble.(hubOrient),
            Cfloat.(hubVel),
            Cfloat.(hubAcc),
            Cfloat.(nacPos),
            Cdouble.(nacOrient),
            Cfloat.(nacVel),
            Cfloat.(nacAcc),
            Cfloat.(rootPos),
            Cdouble.(rootOrient),
            Cfloat.(rootVel),
            Cfloat.(rootAcc),
            numMeshNodes,
            Cfloat.(meshPos),
            Cdouble.(meshOrient),
            Cfloat.(meshVel),
            adi_err.error_status,
            adi_err.error_message) 

        #println("    return from adi_updatestates: $(adi_err.error_status)")
        adi_check_error()

    else
        error("AeroDyn-Inflow instance has not been initialized. Use adiInit() function.")
    end
end

function adiEnd()
    if adi_active
        global adi_active = false
        try
            ccall(adi_sym_end,Cint,
                (Ptr{Cint},         # OUT: ErrStat_C
                Cstring),           # OUT: ErrMsg_C
                adi_err.error_status,
                adi_err.error_message)
        catch
            @error "AeroDyn-Inflow did not end properly"
        end

        Libdl.dlclose(adilib) # Close the library explicitly.
    end
end

function adi_check_error()
    if adi_err.error_status[1] == 0
        adi_err.error_status = [0] # reset error status/message
        adi_err.error_message = string(repeat(" ", 1025))
    elseif adi_err.error_status[1] < adi_abort_error_level
        @warn("Error status " * string(adi_err.error_status[1]) * ": " * string(adi_err.error_message))
        adi_err.error_status = [0] # reset error status/message
        adi_err.error_message = string(repeat(" ", 1025))
    else
        @error("Error status " * string(adi_err.error_status[1]) * ": " * string(adi_err.error_message))
        adi_err.error_status = [0] # reset error status/message
        adi_err.error_message = string(repeat(" ", 1025))
        adiEnd()
        error("AeroDyn-Inflow terminated prematurely.")
    end
end
