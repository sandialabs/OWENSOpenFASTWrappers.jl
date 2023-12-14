global adilib
global sym_calcoutput
global sym_updatestates
global sym_end

mutable struct adiError
    error_status
    error_message
end

"""
    adiPreInit(adilib_filename numTurbines transposeDCM)

Does some pre-initializing of the ADI library to setup arrays for each turbine

# Inputs:
* `adilib_filename::string`: path and name of AeroDyn-Inflow dynamic library
* `numTurbines::int`: required, number of turbines to setup ADI to handle
* `transposeDCM::int`: required, transpose DCM internally in ADI to match calling code convention for direction cosine matrices (default: 1==true)
"""
function adiPreInit(adilib_filename, numTurbines,transposeDCM,;adi_debug=0)

    # Set the error level
    global adi_abort_error_level = 4

    try
        println("Opening AeroDyn-Inflow library at: $adilib_filename")
        global adilib = Libdl.dlopen(adilib_filename) # Open the library explicitly.
        global adi_active = true

        # Get symbols for function calls.
        global adi_sym_preinit          = Libdl.dlsym(adilib, :ADI_C_PreInit)        # Setup turbine data storage internally
        global adi_sym_setuprotor       = Libdl.dlsym(adilib, :ADI_C_SetupRotor)            # Setup for one rotor (initial root positions etc)
        global adi_sym_setrotormotion   = Libdl.dlsym(adilib, :ADI_C_SetRotorMotion)        # set motions on one rotor
        global adi_sym_getrotorloads    = Libdl.dlsym(adilib, :ADI_C_GetRotorLoads)         # get loads from one rotor
        global adi_sym_init             = Libdl.dlsym(adilib, :ADI_C_Init)           # Initialize AeroDyn + InflowWind after rotors setup
        global adi_sym_calcoutput       = Libdl.dlsym(adilib, :ADI_C_CalcOutput)     # Calculate outputs for given outputs and states
        global adi_sym_updatestates     = Libdl.dlsym(adilib, :ADI_C_UpdateStates)   # Advance to next timestep
        global adi_sym_end              = Libdl.dlsym(adilib, :ADI_C_End)            # !!! "c" is capitalized in library, change if errors
        global adi_err = adiError([0], string(repeat(" ", 1025)))
    catch
        error("AeroDyn-InflowWind library could not initialize")
        global adi_active = false
    end

    global adi_abort_error_level = 4

    try 
        ccall(adi_sym_preinit,Cint,
            (Ref{Cint},         # IN: number of turbines to setup ADI for
            Ref{Cint},         # IN: transposeDCM flag
            Ref{Cint},         # IN: adi debug
            Ptr{Cint},          # OUT: error_status
            Cstring),           # OUT: error_message 
            numTurbines,
            transposeDCM,
            adi_debug,
            adi_err.error_status,
            adi_err.error_message)

        adi_check_error()
    catch
        error("AeroDyn-InflowWind library could not initialize turbines")
        global adi_active = false
    end
end


"""
    adiSetupRotor(iturb; )

Sets the initial locations of a single rotor (root orientations/positions etc)

# Inputs:
* `iturb::int`: required, current turbine number
* `isHAWT::bool`: required, false: VAWT or cross-flow turbine, true: HAWT
* `intTurbPos::Array(float)`: required, (x,y,z) position of turbine
* `initHubPos::Array(float)`: required, (x,y,z) position of hub
* `initHubOrient::Array(float)`: required, orientation of hub as 9 element vector of flattened DCM

* `initNacellePos::Array(float)`: required, (x,y,z) position of nacelle
* `initNacelleOrient::Array(float)`: required, orientation of nacelle as 9 element vector of flattened DCM

* `numBlades::int`: required, number of blades
* `initRootPos::Array(float)`: required, size (numBlades,3) position vectors of roots
* `initRootOrient::Array(float)`: required, size (numBlades,9) orientation DCMs flattened to array of 9 element vectors

* `numMeshNodes::Array(int)`: required, number of structural mesh points (total across all blades)
* `initMeshPos::Array(float)`: required, size (sum(numMeshNodes),3) position vectors of mesh points
* `initMeshOrient::Array(float)`: required, size (sum(numMeshNodes),9) orientation DCMs flattened to array of 9 element vectors

"""
function adiSetupRotor(iturb;
    isHAWT             = false,     # Assume VAWT if not specified
    initTurbPos        = zeros(3),  # initial turbine position
    initHubPos         = zeros(3),  # initial position vector of hub
    initHubOrient      = zeros(9),  # initial orientation of hub (flattened 3x3 DCM)
    initNacellePos     = zeros(3),  # initial position vector of nacelle 
    initNacelleOrient  = zeros(9),  # initial orientation of nacelle (flattened 3x3 DCM)
    numBlades          = 3,         # number of blades in system
    initRootPos        = zeros(3),  # initial root position vectors
    initRootOrient     = zeros(9),  # initial root orientation DCMs
    numMeshNodes       = 1,         # number of mesh points representing structural mesh of rotor
    initMeshPos        = zeros(3),  # initial position vectors of all mesh points
    initMeshOrient     = zeros(9)   # initial orientations of all mesh points
    )

    try 
        ccall(adi_sym_setuprotor,Cint,
            (Ref{Cint},         # IN: turbine number
            Ref{Cint},          # IN: isHAWT
            Ref{Cfloat},        # IN: initTurbPos
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
            Ptr{Cint},          # OUT: error_status
            Cstring),           # OUT: error_message
            iturb,
            Cint.(isHAWT),
            Cfloat.(initTurbPos), 
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
            adi_err.error_status,
            adi_err.error_message)

        adi_check_error()
    catch
        error("AeroDyn-InflowWind library could not initialize turbines")
        global adi_active = false
    end

end


"""
    adiSetRotorMotion(iturb; )

Sets the motions for a single turbine rotor

# Inputs:
* `iturb::int`:         required, current turbine number

* `HubPos::Array(float)`: required, (x,y,z) position of hub
* `HubOrient::Array(float)`: required, orientation of hub as 9 element vector of flattened DCM
* `HubVel::Array(float)`: required, (TVx,TVy,TVz,RVx,RVy,RVz) velocity of hub, does not include rotational velocity, so this is extra like from a platform
* `HubAcc::Array(float)`: required, (TAx,TAy,TAz,RAx,RAy,RAz) acceleration of hub, does not include rotational accel, so this is extra like from a platform

* `NacPos::Array(float)`: required, (x,y,z) position of nacelle
* `NacOrient::Array(float)`: required, orientation of nacelle as 9 element vector of flattened DCM
* `NacVel::Array(float)`: required, (TVx,TVy,TVz,RVx,RVy,RVz) velocity of nacelle
* `NacAcc::Array(float)`: required, (TAx,TAy,TAz,RAx,RAy,RAz) acceleration of nacelle

* `RootPos::Array(float)`: required, size (numBlades,3) position vectors of roots
* `RootOrient::Array(float)`: required, size (numBlades,9) orientation DCMs flattened to array of 9 element vectors
* `RootVel::Array(float)`: required, size (numBlades,6) velocity vectors of roots
* `RootAcc::Array(float)`: required, size (numBlades,6) acceleration vectors of roots

* `numMeshNodes::Array(int)`: required, number of structural mesh points (total across all blades)
* `MeshPos::Array(float)`: required, size (sum(numMeshNodes),3) position vectors of mesh points
* `MeshOrient::Array(float)`: required, size (sum(numMeshNodes),9) orientation DCMs flattened to array of 9 element vectors
* `MeshVel::Array(float)`: required, size (sum(numMeshNodes),6) velocity vectors of mesh points
* `MeshAcc::Array(float)`: required, size (sum(numMeshNodes),6) acceleration vectors of mesh points

# Outputs:
* `none`:
"""
function adiSetRotorMotion(iturb,
    hubPos, hubOrient, hubVel, hubAcc,
    nacPos, nacOrient, nacVel, nacAcc,
    rootPos, rootOrient, rootVel, rootAcc,
    numMeshNodes,
    meshPos, meshOrient, meshVel, meshAcc
    )

    if adi_active
        try
            ccall(adi_sym_setrotormotion,Cint,
                (Ref{Cint},         # IN: iturb
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
                Ptr{Cint},          # OUT: error_status
                Cstring),           # OUT: error_message 
                iturb,
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
                adi_err.error_status,
                adi_err.error_message) 
 
            adi_check_error()
        catch
            error("AeroDyn-InflowWind SetRotorMotion could not be called")
            global adi_active = false
        end
    else
        error("AerooDyn-Inflow instance has not been initialized. Use adiInit() function.")
    end
end

"""
    adiInit( output_root_name ; )

calls aerodyn_inflow_init to initialize AeroDyn and InflowWind together

# Inputs:
* `ad_input_file_passed::int`: flag to indicate the AD15 input file is passed as a string [0=false, 1=true]
                                (set to false if passing input file name instead, NOT SUPPORTED YET)
* `ad_input_file::string`: name of input file for AD15 -- this is read by julia and passed to AD15
* `ifw_input_file_passed::int`: flag to indicate the InflowWind input file is passed as a string [0=false, 1=true]
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

* `storeHHVel::int`:   optional, internal parameter for adi_library.  Exposed for convenience, but not needed. [0=false, 1=true]
* `WrVTK::int`:         optional, write VTK output files at all timesteps to visualize AeroDyn 15 meshes [0 none (default), 1 ref, 2 motion]
* `WrVTK_Type::int`:    optional, write VTK output files as [1 surfaces (default), 2 lines, 3 both]
* `VTKNacDim::Array(float*6)`   optional, Nacelle Dimension for VTK visualization x0,y0,z0,Lx,Ly,Lz (m)
* `VTKHubRad::float`:   optional, HubRadius for VTK visualization (m)

* `wrOuts::int`:        optional, file format for writing outputs [0 none (default), 1 txt, 2 binary, 3 both]
* `DT_Outs::float64`:   optional, timestep for outputs to file [0.0 (default) for every timestep]

* `interp_order::int`: optional, interpolation order used internally [1 first order (default), 2 second order]

* `dt::float64`:        required, timestep for AD15 (needed for setting internal constants)
* `t_max::float`:       required, total expected simulation time -- used only for setting VTK counter width

# Outputs:
* `num_channels::int`: number of output channels
* `channel_names::string`: string of output channel names from ADI
* `channel_units::string`: string of output channel units from ADI

"""
function adiInit(output_root_name;
    ad_input_file_passed= 1,
    ad_input_file="none",
    ifw_input_file_passed= 0,
    ifw_input_file="none",
    gravity     =   9.80665,  # Gravitational acceleration (m/s^2)
    defFldDens  =     1.225,  # Air density (kg/m^3)
    defKinVisc  = 1.464E-05,  # Kinematic viscosity of working fluid (m^2/s)
    defSpdSound =     335.0,  # Speed of sound in working fluid (m/s)
    defPatm     =  103500.0,  # Atmospheric pressure (Pa) [used only for an MHK turbine cavitation check]
    defPvap     =    1700.0,  # Vapour pressure of working fluid (Pa) [used only for an MHK turbine cavitation check]
    WtrDpth     =       0.0,  # Water depth (m)
    MSL2SWL     =       0.0,  # Offset between still-water level and mean sea level (m) [positive upward]
    storeHHVel  = 0,          # some internal library stuff we probably don't need to expose [0=false, 1=true]
    WrVTK       = 0,          # write VTK files from adi [0 none, 1 ref, 2 motion]
    WrVTK_Type  = 1,          # write VTK files from adi [1 surfaces, 2 lines, 3 both]
    VTKNacDim   = [-1 ,-1 ,-1 ,2 ,2 ,2],        # Nacelle Dimension for VTK visualization x0,y0,z0,Lx,Ly,Lz (m)
    VTKHubRad   = 0.1,                          # HubRadius for VTK visualization (m)
    wrOuts      = 0,
    DT_Outs     = 0.0,
    interp_order=1,
    dt=0.01,
    t_max=60.0)

    # AeroDyn 15 input file
    if ad_input_file_passed == 0
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
    if ifw_input_file_passed == 0 
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
            Ref{Cint},          # IN: interp_order
            Ref{Cdouble},       # IN: dt
            Ref{Cdouble},       # IN: t_max
            Ref{Cint},          # IN: storeHHVel
            Ref{Cint},          # IN: WrVTK
            Ref{Cint},          # IN: WrVTK_Type
            Ref{Cfloat},        # IN: VTKNacDim
            Ref{Cfloat},        # IN: VTKHubRad
            Ref{Cint},          # IN: wrOuts
            Ref{Cdouble},       # IN: DT_Outs
            Ptr{Cint},          # OUT: num_channels
            Cstring,            # OUT: channel_names
            Cstring,            # OUT: channel_units
            Ptr{Cint},          # OUT: error_status
            Cstring),           # OUT: error_message
            ad_input_file_passed,
            [ad_input_string],
            ad_input_string_length,
            ifw_input_file_passed,
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
            interp_order,
            Cdouble(dt),
            Cdouble(t_max),
            storeHHVel,
            WrVTK,
            WrVTK_Type,
            Cfloat.(VTKNacDim),
            VTKHubRad,
            wrOuts,
            Cdouble(DT_Outs),
            num_channels,
            channel_names,
            channel_units,
            adi_err.error_status,
            adi_err.error_message)

        adi_check_error()
    catch
        error("AeroDyn_Inflow_C_Init failed")
        global adi_active = false
    end

    return num_channels[1], channel_names, channel_units

end


"""
    adiCalcOutput( )

calls AeroDyn_Inflow_C_CalcOutput to calculate resulting loads.  Call this only after SetRotorMotion on all rotors/turbines.

# Inputs:
# `time::c_double`: current timestep
* `num_channels::int`: number of output channels
"""
function adiCalcOutput(time, num_channels)
    out_channel_vals = zeros(Cfloat,1,num_channels)
    if adi_active
        try
            ccall(adi_sym_calcoutput,Cint,
                (Ptr{Cdouble},      # IN: time
                Ptr{Cfloat},        # OUT: out_channel_vals
                Ptr{Cint},          # OUT: error_status
                Cstring),           # OUT: error_message 
                [time],
                out_channel_vals,
                adi_err.error_status,
                adi_err.error_message) 
       
            adi_check_error()
        catch
            error("AeroDyn_Inflow_C_CalcOutput failed")
            global adi_active = false
        end

        return out_channel_vals

    else
        error("AeroDyn-Inflow instance has not been initialized. Use adiInit() function.")
    end
end


"""
    adiGetRotorLoads(iturb; )

Gets the loads for a single turbine rotor

# Inputs:
* `iturb::int`:         required, current turbine number

# Outputs:
* `meshFrcMom::Array(float)`: loads from ADI at mesh nodes
"""
function adiGetRotorLoads(iturb)
    global turbine

    meshFrcMom = zeros(Cfloat,6*turbine[iturb].numMeshNodes)

    if adi_active
        try
            ccall(adi_sym_getrotorloads,Cint,
                (Ref{Cint},         # IN: iturb current turbine number
                Ref{Cint},          # IN: numMeshNodes -- number of attachment points expected (where motions are transferred into HD)
                Ptr{Cfloat},        # OUT: meshFrcMom resulting forces/moments array
                Ptr{Cint},          # OUT: error_status
                Cstring),           # OUT: error_message 
                iturb,
                turbine[iturb].numMeshNodes,
                meshFrcMom,
                adi_err.error_status,
                adi_err.error_message) 
 
            adi_check_error()

            return meshFrcMom

        catch
            error("AeroDyn-InflowWind adiGetRotorLoads could not be called")
            global adi_active = false
        end
    else
        error("AerooDyn-Inflow instance has not been initialized. Use adiInit() function.")
    end
end

"""
    adiUpdateStates( )

calls AeroDyn_Inflow_C_UpdateStates to step ADI forward to the next timestep  Call this only after SetRotorMotion on all rotors/turbines.

# Inputs:
# `time::c_double`: current timestep
# `time_next::c_double`: next timestep
"""
function adiUpdateStates(time, next_time)
    if adi_active
        try
            ccall(adi_sym_updatestates,Cint,
                (Ptr{Cdouble},      # IN: time
                Ptr{Cdouble},       # IN: next_time
                Ptr{Cint},          # OUT: error_status
                Cstring),           # OUT: error_message 
                [time],
                [next_time],
                adi_err.error_status,
                adi_err.error_message) 
         
            adi_check_error()
        catch
            error("AeroDyn_Inflow_C_UpdateStates failed")
            global adi_active = false
        end

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


# Common Structs

"""
    Turbine(R::TF1,omega::TAF1,B::TI1,adi_numbl::TI2,numMeshNodes::TI3,bladeIdx::TAI1,bladeElem::TAI2)
    Turbine(R,omega,B,adi_numbl,numMeshNodes,bladeIdx,bladeElem) = Turbine(R,r,twist,B,adi_numbl,numMeshNodes,bladeIdx,bladeElem)

Contains specications for turbine

NOTE: struts are modeled as blades in AD15 with their root at the hub radius.  In OWENS the strut connects to the tower, so we
  find the node that is closest to the hub radius along the OWENS strut (move it to exactly the hub radius if necessary) and map
  that node as the blade root.  This means we don't pass all the OWENS strut nodes/orientations to AD15, but only a subset.  This
  necessetated adding the bladeIdx and bladeElem indexing arrays.

# Inputs
* `R::TF`: Nominal turbine radius (m)
* `omega::TAF1`: Array of rotational rate corresponding to each azimuthal position, allows for active blade deformation (rad/s)
* `refPos::TAF2`: turbine reference position
* `B::TI`: Number of blades
* `adi_numbl::TI2`: number of adi blades (includes struts)
* `numMeshNodes::TI3`: number of mesh nodes we are passing motions/loads to/from.  Derived from bladeIdx
* `bladeIdx::TAI1`: an array of start and end nodes for each AD15 blade. Used to index into mesh
* `bladeElem::TAI2`: an array of start and end indices into the orientation angles stored in myort.Psi_d, myort.Theta_d, and myort.Twist_d.
* `Mesh::GyricFEA:Mesh`: mesh without deformation. Storing a copy here to clean up the code a bit.
* `Ort::GyricFEA:Ort`: orientations of all elements without deformation.  Storing here since we don't have a good way of passing this info during simulation
* `isHAWT::bool`: flag to switch between VAWT and HAWT conventions 
 

# Outputs:
* `none`:

"""
#FIXME: change refPos to baseOriginInit -- follow AD15 driver naming here.
struct Turbine{TF1,TAF1,TAF2,TI1,TI2,TI3,TAI1,TAI2,TM,TO,TB}
    R::TF1
    omega::TAF1
    refPos::TAF2
    B::TI1
    adi_numbl::TI2
    numMeshNodes::TI3
    bladeIdx::TAI1
    bladeElem::TAI2
    Mesh::TM
    Ort::TO
    isHAWT::TB
end
Turbine(R,omega,refPos,B,adi_numbl,numMeshNodes,bladeIdx,bladeElem,mymesh,myort) = Turbine(R,omega,refPos,B,adi_numbl,numMeshNodes,bladeIdx,bladeElem,mymesh,myort,true)

"""
Structure

Contains the position and orientation info passed to AD15.  These positions are all in global and include the turbine RefPos.

# Inputs:
* `hubPos::TAF1`:
* `hubOrient::TAF2`:
* `hubVel::TAF3`:
* `hubAcc::TAF4`:
* `nacPos::TAF5`:
* `nacOrient::TAF6`:
* `nacVel::TAF7`:
* `nacAcc::TAF8`:
* `rootPos::TAF9`:
* `rootOrient::TAF10`:
* `rootVel::TAF11`:
* `rootAcc::TAF12`:
* `meshPos::TAF13`:
* `meshOrient::TAF14`:
* `meshVel::TAF15`:
* `meshAcc::TAF16`:

# Outputs:
* `none`:
"""
mutable struct Structure{TAF1,TAF2,TAF3,TAF4,TAF5,TAF6,TAF7,TAF8,TAF9,TAF10,TAF11,TAF12,TAF13,TAF14,TAF15,TAF16}
    hubPos::TAF1
    hubOrient::TAF2
    hubVel::TAF3
    hubAcc::TAF4
    nacPos::TAF5
    nacOrient::TAF6
    nacVel::TAF7
    nacAcc::TAF8
    rootPos::TAF9
    rootOrient::TAF10
    rootVel::TAF11
    rootAcc::TAF12
    meshPos::TAF13
    meshOrient::TAF14
    meshVel::TAF15
    meshAcc::TAF16
end

"""
Environment(rho::TF1, dt::TF2, steplast::TAI)
Environment(rho) = Environment(rho)

Contains specications for turbine environment/operating conditions as well as some backend memory

# Inputs
* `rho::TF1`: Working fluid density (kg/m^3)
* `dt::TF2`: timestep for ADI
* `num_channels::TI1`: number of output channels from AD15
* `steplast::TAI`: prior simulation step index, used for unsteady wake propogation

# Outputs:
* `none`:

"""
struct Environment{TF1,TF2,TI1,TAI}
    rho::TF1
    dt::TF2
    num_channels::TI1
    steplast::TAI
end

Environment(rho,dt,num_channels) = Environment(rho,dt,num_channels,zeros(Int,1))

global adi_initialized = false
global turbine
global turbenv
global turbstruct

"""
setupTurb(adi_lib,ad_input_file,ifw_input_file,adi_rootname,bld_x,bld_z,B;
    rho = 1.225,
    gravity     =   9.80665,
    defKinVisc  = 1.464E-05,
    defSpdSound =     335.0,
    defPatm     =  103500.0,
    defPvap     =    1700.0,
    WtrDpth     =       0.0,
    MSL2SWL     =       0.0,
    storeHHVel  = 0,    #false
    transposeDCM= 1,    #true
    WrVTK       = 2,
    WrVTK_Type  = 3,
    VTKNacDim   = [-.10 ,-.10 ,-.10 ,.2 ,.2 ,.2],
    VTKHubRad   = 1.0,
    adi_nstrut  = 2,
    adi_dt      = 0.05,
    adi_tmax    = 10,
    adi_wrOuts  = 0,
    adi_DT_Outs = 0.0,
    isHAWT      = false,                          # false: VAWT or cross-flow turbine, true: HAWT
    numTurbines = 1)

Initializes aerodynamic models and sets up backend persistent memory to simplify intermittent calling within coupled solver loops

# Inputs
* `adi_lib`: path to adi library (.so, .dylib, .dll)
* `ad_input_file`: input file for aerodyn15
* `ifw_input_file`: input file for inflow wind
* `adi_rootname`: rootname for vtk outputs
* `bld_x`: Blade x shape
* `bld_z`: Blade z shape
* `B`: Number of blades
* `rho`: working fluid density (kg/m^3)
* `gravity`: Gravitational acceleration (m/s^2)
* `defKinVisc`: Kinematic viscosity of working fluid (m^2/s)
* `defSpdSound`: Speed of sound in working fluid (m/s)
* `defPatm`: Atmospheric pressure (Pa) [used only for an MHK turbine cavitation check]
* `defPvap`: Vapour pressure of working fluid (Pa) [used only for an MHK turbine cavitation check]
* `WtrDpth`: Water depth (m)
* `MSL2SWL`: Offset between still-water level and mean sea level (m) [positive upward]
* `storeHHVel`: unused here
* `transposeDCM`: 0=false, 1=true transpose DCM internally for calculations
* `WrVTK`: write VTK files from adi to directory adi-vtk [0 none, 1 ref, 2 motion]
* `WrVTK_Type`: write VTK files from adi to directory adi-vtk [1 surface, 2 lines, 3 both]
* `VTKNacDim`: Nacelle Dimension for VTK visualization x0,y0,z0,Lx,Ly,Lz (m)
* `VTKHubRad`: HubRadius for VTK visualization (m)
* `adi_wrOuts`: file format to write to
* `adi_DT_Outs`: output timestep to write at
* `adi_nstrut`: create_mesh_struts is hard coded for 2 struts per blade
* `adi_dt`: timestep
* `adi_tmax`: maximum time
* `hubPos`: hub position in global coordinates, 3-vector (m). NOTE: AD15 assumes a different hub location than OWENS
* `hubAngle`: hub axis angle, 3-vector (deg), no rotation from spinning
* `nacPos`: nacelle position in global coordinates, 3-vector (m). NOTE: AD15 assumes a different hub location than OWENS
* `nacAngle`: nacelle axis angle, 3-vector (deg)
* `numTurbines`: number of turbines
* `isHAWT`: # false: VAWT or cross-flow turbine, true: HAWT


# Outputs:
* `none`:

"""
function setupTurb(adi_lib,ad_input_file,ifw_input_file,adi_rootname,bld_x,bld_z,B,Ht,mymesh,myort,bladeIdx,bladeElem;
    rho         =     1.225,  # fluid density (kg/m^3)
    gravity     =   9.80665,  # Gravitational acceleration (m/s^2)
    defKinVisc  = 1.464E-05,  # Kinematic viscosity of working fluid (m^2/s)
    defSpdSound =     335.0,  # Speed of sound in working fluid (m/s)
    defPatm     =  103500.0,  # Atmospheric pressure (Pa) [used only for an MHK turbine cavitation check]
    defPvap     =    1700.0,  # Vapour pressure of working fluid (Pa) [used only for an MHK turbine cavitation check]
    WtrDpth     =       0.0,  # Water depth (m)
    MSL2SWL     =       0.0,  # Offset between still-water level and mean sea level (m) [positive upward]
    transposeDCM= 1,          # 0=false, 1=true transpose DCM internally for calculations
    WrVTK       = 2,          # write VTK files from adi to directory adi-vtk [0 none, 1 ref, 2 motion]
    WrVTK_Type  = 3,          # write VTK files from adi to directory adi-vtk [1 surface, 2 lines, 3 both]
    VTKHubRad   = 1.0,                          # HubRadius for VTK visualization (m)
    adi_wrOuts  = 0,                            # output file format [0 none, 1 txt, 2 binary, 3 both]
    adi_DT_Outs = 0.0,                          # output frequency (resets to dt if smaller)
    adi_dt      = 0.05,                         # random default
    adi_tmax    = 10,                           # end time
    isHAWT      = false,                          # false: VAWT or cross-flow turbine, true: HAWT
    adi_debug   = 0,                              #0 is no debug outputs
    VTKNacDim   = [-.10 ,-.10 ,-.10 ,.2 ,.2 ,.2],        # Nacelle Dimension for VTK visualization x0,y0,z0,Lx,Ly,Lz (m)
    numTurbines = 1,
    adi_nstrut  = [2 for i=1:numTurbines],                            # create_mesh_struts is hard coded for 2 struts per blade
    omega       = [0 for i=1:numTurbines],                          # rad/s
    hubPos      = [zeros(3) for i=1:numTurbines],                      # m
    hubAngle    = [zeros(3) for i=1:numTurbines],                      # deg , no rotation from spinning
    refPos      = [zeros(3) for i=1:numTurbines],                      # turbine location
    nacPos      = [zeros(3) for i=1:numTurbines],                      # m
    nacAngle    = [zeros(3) for i=1:numTurbines],                      # m
    )

    # load library and set number of turbines
    try
        adiPreInit(adi_lib,numTurbines,transposeDCM;adi_debug)
    catch
        error("AeroDyn-InflowWind library could not initialize")
        global adi_active = false
    end


    global turbine = Array{Turbine}(undef, numTurbines)
    global turbstruct = Array{Structure}(undef, numTurbines)

    for iturb = 1:numTurbines

        # Set up structs for the entire turbine
        if isHAWT
            adi_numbl = B[iturb]
        else #TODO N struts
            adi_numbl = B[iturb] + B[iturb]*adi_nstrut[iturb]    # Count struts as blades (strut each side of tower counted separately)
        end
        Radius = maximum(bld_x[iturb])
        
        numMeshNodes = getAD15numMeshNodes(bladeIdx[iturb])

        turbine[iturb] = Turbine(Radius,omega[iturb],refPos[iturb],B[iturb],adi_numbl,
            numMeshNodes,bladeIdx[iturb],bladeElem[iturb],mymesh[iturb],myort[iturb],isHAWT)

        # Mesh info for ADI
        # set the origin for AD15 at the top of the "tower" (Ht in this setup)
        refPosTmp=Float32.(refPos[iturb][1:3])
        # nacelle -- not actually used here since we don't consider loads.
        nacPosTmp = Float32.(nacPos[iturb][1:3])
        nacOrient = createGeneralTransformationMatrix(nacAngle[iturb],[1,2,3])
        # nacOrient = Float64.([1,0,0,0,1,0,0,0,1])
        # hub -- align to origin for now
        hubPosTmp = Float32.(hubPos[iturb][1:3])
        # hubOrient = Float64.([1,0,0,0,1,0,0,0,1])

        # set initial motion to 0
        u_j     = zeros(mymesh[iturb].numNodes*6)
        udot_j  = zeros(mymesh[iturb].numNodes*6)
        uddot_j = zeros(mymesh[iturb].numNodes*6)
        azi     = 0.0

        # blade roots (2nd is rotated 180 degrees about z)
        rootPos     = getRootPos(turbine[iturb],u_j,azi,nacPosTmp,hubPosTmp,hubAngle[iturb])       # get root positions of all AD15 blades (blades + struts in OWENS)
        rootOrient  = getRootDCM(turbine[iturb],u_j,azi,hubAngle[iturb])              # get orientations of all AD15 blades   (blades + struts in OWENS)

        # Multiple mesh points along all blades for full structural mesh representation in ADI
        meshPos      = getAD15MeshPos(turbine[iturb],u_j,azi,nacPosTmp,hubPosTmp,hubAngle[iturb])  # get positions of all AD15 nodes (blades + struts in OWENS)
        meshOrient   = getAD15MeshDCM(turbine[iturb],u_j,azi,hubAngle[iturb])         # get orientations of all AD15 blades   (blades + struts in OWENS)


        # AD15 node velocities/accelerations
        #TODO: If OWENS sets these at initialization, need to transfer values here. 
        hubVel    = zeros(Float32,2*size(hubPosTmp,1))
        hubAcc    = zeros(Float32,2*size(hubPosTmp,1))
        nacVel    = zeros(Float32,2*size(nacPosTmp,1))
        nacAcc    = zeros(Float32,2*size(nacPosTmp,1))
        # rootVel   = zeros(Float32,2*size(rootPos,1),size(rootPos,2))
        # rootAcc   = zeros(Float32,2*size(rootPos,1),size(rootPos,2))
        # meshVel   = zeros(Float32,2*size(meshPos,1),size(meshPos,2))
        # meshAcc   = zeros(Float32,2*size(meshPos,1),size(meshPos,2))

              
        rootVel,rootAcc = getRootVelAcc(turbine[iturb],rootPos,udot_j,uddot_j,azi,omega[iturb],zero(omega[iturb]),nacPosTmp,hubPosTmp,hubAngle[iturb],hubVel,hubAcc)       # get root vel/acc of all AD15 blades   (blades + struts in OWENS)

        meshVel,meshAcc = getAD15MeshVelAcc(turbine[iturb],meshPos,udot_j,uddot_j,azi,omega[iturb],zero(omega[iturb]),nacPosTmp,hubPosTmp,hubAngle[iturb],hubVel,hubAcc)   # get mesh vel/acc of all AD15 blades   (blades + struts in OWENS)

        # hub
        CG2H = calcHubRotMat(hubAngle[iturb], azi;rot_axis = 1) 
   
        if !turbine[iturb].isHAWT
            CN2P = createGeneralTransformationMatrix([-90,180],[2,3])  
            CG2H = CG2H*CN2P
        end
        
        hubOrient    = vec(CG2H')

        # Initialize outputs and resulting mesh forces
        try
            adiSetupRotor(iturb;
                initTurbPos        = refPos[iturb],     #
                isHAWT             = isHAWT,            # 0: VAWT or cross-flow turbine, 1: HAWT
                initHubPos         = hubPosTmp,         # 3
                initHubOrient      = hubOrient,         # 9
                initNacellePos     = nacPosTmp,            # 3 -- not actually used
                initNacelleOrient  = nacOrient,         # 9 -- not actually used
                numBlades          = adi_numbl,
                initRootPos        = rootPos,           # 3*adi_numbl
                initRootOrient     = rootOrient,        # 9*adi_numbl
                numMeshNodes       = numMeshNodes,
                initMeshPos        = meshPos,           # 3*numMeshNodes
                initMeshOrient     = meshOrient)        # 9*numMeshNodes
        catch
            error("Could not initialize rotor")
            global adi_active = false
        end

        turbstruct[iturb]=Structure(
                hubPosTmp,  hubOrient,  hubVel,  hubAcc,
                nacPosTmp,  nacOrient,  nacVel,  nacAcc,
                rootPos,   rootOrient, rootVel, rootAcc,
                meshPos,   meshOrient, meshVel, meshAcc,
            )

    end

    num_channels, channel_names, channel_units = adiInit(adi_rootname;
        ad_input_file_passed= 0,
        ad_input_file=ad_input_file,
        ifw_input_file_passed= 0,
        ifw_input_file=ifw_input_file,
        gravity     = gravity,
        defFldDens  = rho,
        defKinVisc  = defKinVisc,
        defSpdSound = defSpdSound,
        defPatm     = defPatm,
        defPvap     = defPvap,
        WtrDpth     = WtrDpth,
        MSL2SWL     = MSL2SWL,
        storeHHVel  = false,    # unused
        WrVTK       = WrVTK,
        WrVTK_Type  = WrVTK_Type,
        VTKNacDim   = VTKNacDim,
        VTKHubRad   = VTKHubRad,
        wrOuts      = adi_wrOuts,
        DT_Outs     = adi_DT_Outs,
        interp_order=1,
        dt=adi_dt,
        t_max=adi_tmax,
        )

    # can add some additional things here if needed
    adi_initialized = true

    # store channel info
    global turbenv = Environment(rho,adi_dt,num_channels)

    # Set the motions for each turbine
    for iturb = 1:numTurbines
        setRotorMotion(iturb)
    end
    t_initial = 0.0

    # calculate initial values
    out_channel_vals = adiCalcOutput(t_initial,
            turbenv.num_channels);

    # get the resulting forces and moments for each turbine rotor
    for iturb = 1:numTurbines
        meshFrcMom = adiGetRotorLoads(iturb)
    end

end


"""
    setRotorMotion(iturb,turbstruct

Set the rotor motion for one turbine

# Inputs
* `iturb::int`:         required, current turbine number
"""
function setRotorMotion(iturb)
    global turbine
    global turbstruct
    # Initialize outputs and resulting mesh forces
    try
        adiSetRotorMotion(iturb,
            turbstruct[iturb].hubPos,  turbstruct[iturb].hubOrient,  turbstruct[iturb].hubVel,  turbstruct[iturb].hubAcc,
            turbstruct[iturb].nacPos,  turbstruct[iturb].nacOrient,  turbstruct[iturb].nacVel,  turbstruct[iturb].nacAcc,
            turbstruct[iturb].rootPos, turbstruct[iturb].rootOrient, turbstruct[iturb].rootVel, turbstruct[iturb].rootAcc,
            turbine[iturb].numMeshNodes,
            turbstruct[iturb].meshPos, turbstruct[iturb].meshOrient, turbstruct[iturb].meshVel, turbstruct[iturb].meshAcc);
    catch
        error("Could not set rotor motion turbine $iturb")
        global adi_active = false
    end
end



"""
    deformAD15(u_j,udot_j,uddot_j,azi,Omega_rad,OmegaDot_rad,hubPos,hubAngle,hubVel,hubAcc)

Sets the inputs for AD15.

# Inputs
* `u_j`:            mesh displacements -- in hub coordinates, (m,rad)
* `udot_j`:         mesh velocity      -- in hub coordinates, (m/s,rad/s)
* `uddot_j`:        mesh velocity      -- in hub coordinates, (m/s,rad/s)
* `azi`:            current azimuth (rad)
* `Omega_rad`:      angular velocity of hub about hub axis (rad/s)
* `OmegaDot_rad`:   angular acceleration of hub about hub axis (rad/s^2)
* `hubPos`:         current global hubPos (x,y,z) vector (m)
* `hubAngle`:       3 angle set for hub orientation (rad), no rotation from spinning
* `hubVel`:         hub velocity in global coords, 6-vector (m/s,rad/s), does not include rotational velocity, so this is extra like from a platform
* `hubAcc`:         hub acceleration in global coords, 6-vector (m/s^2,rad/s^2), does not include rotational accel, so this is extra like from a platform
"""
function deformAD15(u_j,udot_j,uddot_j,azi,Omega_rad,OmegaDot_rad,hubPos,hubAngle,hubVel,hubAcc)

    global turbine
    global turbenv
    global turbstruct
    
    numTurbines = length(turbine)
    for iturb = 1:numTurbines
        if turbine[iturb].isHAWT #Add in the rotational components here 
            hubVel[iturb][4] = Omega_rad[iturb]
            hubAcc[iturb][4] = OmegaDot_rad[iturb]
        else
            hubVel[iturb][6] = Omega_rad[iturb]
            hubAcc[iturb][6] = OmegaDot_rad[iturb]
        end

        nacPos = turbstruct[iturb].nacPos #TODO: nacelle position deformation
        # Root
        turbstruct[iturb].rootPos                    = getRootPos(turbine[iturb],u_j[iturb],azi[iturb],nacPos,hubPos[iturb],hubAngle[iturb])                                                                         # get root positions of all AD15 blades (blades + struts in OWENS)
        turbstruct[iturb].rootOrient                 = getRootDCM(turbine[iturb],u_j[iturb],azi[iturb],hubAngle[iturb])                                                                                # get orientations of all AD15 blades   (blades + struts in OWENS)
        turbstruct[iturb].rootVel,turbstruct[iturb].rootAcc = getRootVelAcc(turbine[iturb],turbstruct[iturb].rootPos,udot_j[iturb],uddot_j[iturb],azi[iturb],Omega_rad[iturb],OmegaDot_rad[iturb],nacPos,hubPos[iturb],hubAngle[iturb],hubVel[iturb],hubAcc[iturb])       # get root vel/acc of all AD15 blades   (blades + struts in OWENS)

        # Mesh
        turbstruct[iturb].meshPos                    = getAD15MeshPos(turbine[iturb],u_j[iturb],azi[iturb],nacPos,hubPos[iturb],hubAngle[iturb])                                                                     # get mesh positions of all AD15 blades (blades + struts in OWENS)
        turbstruct[iturb].meshOrient                 = getAD15MeshDCM(turbine[iturb],u_j[iturb],azi[iturb],hubAngle[iturb])                                                                            # get orientations of all AD15 blades   (blades + struts in OWENS)
        turbstruct[iturb].meshVel,turbstruct[iturb].meshAcc = getAD15MeshVelAcc(turbine[iturb],turbstruct[iturb].meshPos,udot_j[iturb],uddot_j[iturb],azi[iturb],Omega_rad[iturb],OmegaDot_rad[iturb],nacPos,hubPos[iturb],hubAngle[iturb],hubVel[iturb],hubAcc[iturb])   # get mesh vel/acc of all AD15 blades   (blades + struts in OWENS)

        # hub
        #FIXME: this is not complete.  The hubVel is probably not correctly set.
        CG2H = calcHubRotMat(hubAngle[iturb], azi[iturb];rot_axis = 1)  
        if !turbine[iturb].isHAWT #orient hub to expected rotatation about x, so align x with global z axis
            CN2P = createGeneralTransformationMatrix([-90,180],[2,3])  
            CG2H = CG2H*CN2P
        end
        # end 
        turbstruct[iturb].hubPos       = hubPos[iturb]
        turbstruct[iturb].hubOrient    = vec(CG2H') # this should have the rotation about x, so for a vawt, x is up.
        turbstruct[iturb].hubVel       = hubVel[iturb] # this should be in global
        turbstruct[iturb].hubAcc       = hubAcc[iturb]

        # Nacelle   FIXME: add this later

        #TODO:  Transfer mesh structural deformation from OWENS, convert to global coords (if needed), and then apply to turbstruct
        #            turbstruct.nacPos,  turbstruct.nacOrient,  turbstruct.nacVel,  turbstruct.nacAcc,

        # Set the motions for each turbine (passed in global values)
        setRotorMotion(iturb)

    end
end


"""
    advanceAD15(t_new;ts=2*pi/(turbine.omega[1]*turbine.ntheta))

Runs a previously initialized aero model (see ?setupTurb) in the unsteady mode (can be repeateadly called, or called for a specific time, or repeatedly called for sections of time)

# Inputs
* `t_new::float`: new time (s); will run from last time specified from the last call, to the current time specified, or from t=ts if the first time called
* `mesh::`:     GyricFEA mesh for the turbine structure
* `azi::`:      hub azimuth (radians)
* `dt::float`:  optional timestep

# Outputs:
* `n_steps`: number timesteps taken
* `Fx`: Array(sum(numMeshNodes),ntheta) Turbine Fx (N)
* `Fy`: Array(sum(numMeshNodes),ntheta) Turbine Fy (N)
* `Fz`: Array(sum(numMeshNodes),ntheta) Turbine Fz (N)
* `Mx`: Array(sum(numMeshNodes),ntheta) Turbine Mx (N-m)
* `My`: Array(sum(numMeshNodes),ntheta) Turbine My (N-m)
* `Mz`: Array(sum(numMeshNodes),ntheta) Turbine Mz (N-m)

"""
function advanceAD15(t_new,mesh,azi;dt=turbenv.dt)
    global turbine
    global turbenv
    global turbstruct
    global dt
    dt = turbenv.dt

    n_steps=1   # hard code for now
    numTurbines = length(turbine)

    # Loads -- set the output array to the size of the OWENS mesh
    Fx = Array{Matrix{Float64}}(undef, numTurbines)
    Fy = Array{Matrix{Float64}}(undef, numTurbines)
    Fz = Array{Matrix{Float64}}(undef, numTurbines)
    Mx = Array{Matrix{Float64}}(undef, numTurbines)
    My = Array{Matrix{Float64}}(undef, numTurbines)
    Mz = Array{Matrix{Float64}}(undef, numTurbines)
    for iturb = 1:numTurbines
        Fx[iturb] = zeros(mesh[iturb].numNodes,n_steps) #zeros(mesh[iturb].numNodes,n_steps)
        Fy[iturb] = zeros(mesh[iturb].numNodes,n_steps) #zeros(mesh[iturb].numNodes,n_steps)
        Fz[iturb] = zeros(mesh[iturb].numNodes,n_steps) #zeros(mesh[iturb].numNodes,n_steps)
        Mx[iturb] = zeros(mesh[iturb].numNodes,n_steps) #zeros(mesh[iturb].numNodes,n_steps)
        My[iturb] = zeros(mesh[iturb].numNodes,n_steps) #zeros(mesh[iturb].numNodes,n_steps)
        Mz[iturb] = zeros(mesh[iturb].numNodes,n_steps) #zeros(mesh[iturb].numNodes,n_steps)
    end

    # conversion to hub coordinates (rotating)
    step1 = 0 #initialize scope
    for istep = 1:n_steps

        t = t_new - dt      # time ADI is starting from (used to copy states)
        # Call update states to go from T to T+dt       (NOTE: T and T+dt must be different)
        adiUpdateStates(t, t_new)

        out_channel_vals = adiCalcOutput(t_new,
                turbenv.num_channels);

        for iturb = 1:numTurbines
            # get the loads for the turbine
            meshFrcMom = adiGetRotorLoads(iturb)

            # copy mesh forces from meshFrcMom
            #   the AD15 mesh is a subset of points from the OWENS mesh. The bladeIdx array indicates the start and end nodes in the
            #   OWENS mesh for the current AD15 blade (struts are blades in AD15).
            #   This unpacking could be made much more compact by using only a single returned force array, but that can be done some
            #   other time. For now I want something that I can test -- ADP
            iBl = 0                                     # index before next blade in meshForceMom from AD15
            for ibld = 1:size(turbine[iturb].bladeIdx,1)       # blade or strut index
                idx1 = turbine[iturb].bladeIdx[ibld,1]         # start index of this blade in OWENS mesh
                idx2 = turbine[iturb].bladeIdx[ibld,2]         # end   index of this blade in OWENS mesh
                # AD15 blade may be upside down compared to OWENS mesh, so invert if needed
                sgn = sign(idx2-idx1)
                nNodes = abs(idx2-idx1)+1               # number of nodes on this blade/strut
                # step through all nodes on this blade/strut
                for iNode=0:nNodes-1                    # shift indices for simplicity
                    Fx[iturb][idx1 + sgn*iNode,istep] = meshFrcMom[iBl + iNode*6 + 1]      # Fx
                    Fy[iturb][idx1 + sgn*iNode,istep] = meshFrcMom[iBl + iNode*6 + 2]      # Fy
                    Fz[iturb][idx1 + sgn*iNode,istep] = meshFrcMom[iBl + iNode*6 + 3]      # Fz
                    Mx[iturb][idx1 + sgn*iNode,istep] = meshFrcMom[iBl + iNode*6 + 4]      # Mx
                    My[iturb][idx1 + sgn*iNode,istep] = meshFrcMom[iBl + iNode*6 + 5]      # My
                    Mz[iturb][idx1 + sgn*iNode,istep] = meshFrcMom[iBl + iNode*6 + 6]      # Mz
                end
                iBl = iBl + nNodes*6        # node index before next blade start from AD15
            end

            # convert to hub rotating coordinates with azimuth angle
            #   There is undoubtedly a much more elegant way to do this
            # TODO: need to check that this is actually correct
            # conversion to hub coordinates (rotating)
            # TODO: for HAWT ensure loads are properly mapped from the aerodyn Rx FOR to the owens Rz FOR
            if turbine[iturb].isHAWT
                CG2H = calcHubRotMat(zeros(3), azi[iturb];rot_axis = 1)
            else
                CG2H = calcHubRotMat(zeros(3), azi[iturb];rot_axis = 3)
            end
            for iNode=1:mesh[iturb].numNodes
                FMg = [Fx[iturb][iNode] Fy[iturb][iNode] Fz[iturb][iNode] Mx[iturb][iNode] My[iturb][iNode] Mz[iturb][iNode]]
                FM = frame_convert(FMg, CG2H)
                if turbine[iturb].isHAWT #TODO incorporate into the CG2H matrix
                    Fx[iturb][iNode,istep] = FM[3]
                    Fy[iturb][iNode,istep] = FM[2]
                    Fz[iturb][iNode,istep] = FM[1]
                    Mx[iturb][iNode,istep] = FM[6]
                    My[iturb][iNode,istep] = FM[5]
                    Mz[iturb][iNode,istep] = FM[4]
                else
                    Fx[iturb][iNode,istep] = FM[1]
                    Fy[iturb][iNode,istep] = FM[2]
                    Fz[iturb][iNode,istep] = FM[3]
                    Mx[iturb][iNode,istep] = FM[4]
                    My[iturb][iNode,istep] = FM[5]
                    Mz[iturb][iNode,istep] = FM[6]
                end
            end
        end
    end

    return n_steps,Fx,Fy,Fz,Mx,My,Mz
end


"""
End ADI and clear data
"""
function endTurb()
    adiEnd()
    adi_initialized = false
end

# Outstanding questions
#   1. why is azi negative in finding positions, but positive in orientations?
#   2. For calculating tangential velocity terms, we use the Z axis of the hub orientation DCM.
#       Q:  Should this be taken as CG2H[:,3] or CG2H[3,:]?  See the getRootVelAcc and getMeshVelAcc routines 

# NOTE:
#   - Omega_j and OmegaDot_j in OWENS are in rev/s and rev/s^2.
#   - Omega_rad and OmegaDot_rad in this file are in rad/s and rad/s^2.  So convert before passing in.



"""
getRootPos(turbine,u_j,azi,nacPos,hubPos,hubAngle)

Extract the root positions for all ADI blades

# Inputs
* `turbine`:    turbine data storage
* `u_j`:        mesh displacements -- in hub coordinates, (m,rad)
* `azi`:        current azimuth (rad)
* `nacPos`:     current global nacPos (x,y,z) vector (m)
* `hubPos`:     current global hubPos (x,y,z) vector (m)
* `hubAngle`:   3 angle set for hub orientation (rad), no rotation from spinning
"""
function getRootPos(turbine,u_j,azi,nacPos,hubPos,hubAngle)
    RootPos     = zeros(Float32,3,turbine.adi_numbl)
    # conversion from hub coordinates to global
    if turbine.isHAWT
        CG2H = calcHubRotMat(hubAngle, azi;rot_axis = 1)
    else
        CG2H = calcHubRotMat(hubAngle, azi;rot_axis = 3)
    end
    CH2G = CG2H
    # blades
    for ibld = 1:turbine.adi_numbl
        idx=turbine.bladeIdx[ibld,1]
        x=turbine.Mesh.x[idx] + u_j[(idx-1)*6+1]
        y=turbine.Mesh.y[idx] + u_j[(idx-1)*6+2]
        z=turbine.Mesh.z[idx] + u_j[(idx-1)*6+3]
        #println("Blade $ibld bottom at [$x,$y,$z] at index $idx")
        RootPos[:,ibld] = [x y z] * CH2G + nacPos' #+ hubPos'
    end
    return RootPos
end

"""
getRootVelAcc(turbine,rootPos,udot_j,uddot_j,azi,Omega_rad,OmegaDot_rad,nacPos,hubPos,hubAngle,hubVel,hubAcc)

Extract the root velocities and accelerations for all ADI blades

# Inputs
* `turbine`:        turbine data storage
* `rootPos`:        root positions from call to getRootPos
* `azi`:            current azimuth (rad)
* `Omega_rad`:      angular velocity of hub about hub axis (rad/s)
* `OmegaDot_rad`:   angular acceleration of hub about hub axis (rad/s^2)
* `nacPos`:         current global nacPos (x,y,z) vector (m)
* `hubPos`:         current global hubPos (x,y,z) vector (m)
* `hubAngle`:       3 angle set for hub orientation (rad) , no rotation from spinning
* `hubVel`:         hub velocity in global coords, 6-vector (m/s,rad/s), does not include rotational velocity, so this is extra like from a platform
* `hubAcc`:         hub acceleration in global coords, 6-vector (m/s^2,rad/s^2), does not include rotational accel, so this is extra like from a platform
"""
function getRootVelAcc(turbine,rootPos,udot_j,uddot_j,azi,Omega_rad,OmegaDot_rad,nacPos,hubPos,hubAngle,hubVel,hubAcc)
    RootVel     = zeros(Float32,6,turbine.adi_numbl)
    RootAcc     = zeros(Float32,6,turbine.adi_numbl)
    ### 1. calculate relative velocity from mesh distortions
    # conversion from hub coordinates to global
    if turbine.isHAWT
        CG2H = calcHubRotMat(hubAngle, azi;rot_axis = 1)
    else
        CG2H = calcHubRotMat(hubAngle, azi;rot_axis = 3)
    end
    # end
    CH2G = CG2H
    # blades #TODO: check the uddot, may be all 0s and shouldn't be, will be an issue for MHK turbines added mass
    for ibld = 1:turbine.adi_numbl
        tmp=turbine.bladeIdx[ibld,1]
        idx=(tmp-1)*6   # just before the node of interest
        RootVel[1:3,ibld] =  udot_j[idx+1:idx+3]' * CH2G         # translation Vel (m/2)
        RootVel[4:6,ibld] =  udot_j[idx+4:idx+6]' * CH2G         # rotation    Vel (rad/s)
        RootAcc[1:3,ibld] = uddot_j[idx+1:idx+3]' * CH2G         # translation Acc (m/s^2)
        RootAcc[4:6,ibld] = uddot_j[idx+4:idx+6]' * CH2G         # rotation    Acc (rad/s^2)
    end
  
    ### 2. Tangential velocity due to hub rotation
    # calculate distance of point from hub axis, multiply by Omega_rad for tangential velocity component
    # hub axis vector in global coordinates
    if turbine.isHAWT
        hubAxis = CG2H[:,1]
    else
        hubAxis = CG2H[:,3]
    end

    for ibld = 1:size(RootVel,2)
        # Global coordinates
        # tangential velocity and acceleration, based on distance to hub axis 
        TanVel = cross(    Omega_rad*hubAxis, (rootPos[1:3,ibld]-nacPos)) / norm(hubAxis)
        TanAcc = cross( OmegaDot_rad*hubAxis, (rootPos[1:3,ibld]-nacPos)) / norm(hubAxis)
        RootVel[1:3,ibld] = RootVel[1:3,ibld] + TanVel
        RootVel[4:6,ibld] = RootVel[4:6,ibld] + Omega_rad*hubAxis
        RootVel[1:3,ibld] = RootVel[1:3,ibld] + TanAcc
        RootVel[4:6,ibld] = RootVel[4:6,ibld] + OmegaDot_rad*hubAxis
    end

    ### 3. add in contributions from hub motion in global coordinates #TODO: nac velocity
    for ibld = 1:size(RootVel,2)
        RootVel[1:3,ibld] = RootVel[1:3,ibld] + hubVel[1:3]
        RootAcc[1:3,ibld] = RootAcc[1:3,ibld] + hubAcc[1:3]
    end
    return RootVel,RootAcc
end


"""
getAD15numMeshNodes(bladeIdx)

Find the number of mesh points we will pass
"""
function getAD15numMeshNodes(bladeIdx)
    # Find number of nodes -- note that we skip some OWENS mesh nodes along struts.  Also note that this method allows for a
    # different number of nodes in each blade
    numMeshNodes = 0
    for i=1:size(bladeIdx,1)
        numMeshNodes = numMeshNodes + abs(bladeIdx[i,2] - bladeIdx[i,1]) + 1     # include all nodes in range 
    end
    return numMeshNodes
end

"""
getAD15MeshPos(turbine,u_j,azi,nacPos,hubPos,hubAngle)

Extract the mesh points for all the AD15 blades
  ordering here is important
      1. root to tip of blades,        in blade order
      2. root to tip of bottom struts, in blade order
      3. root to tip of top    struts, in blade order

# Inputs
* `turbine`:    turbine data storage
* `u_j`:        mesh displacements -- in hub coordinates, (m,rad)
* `azi`:        current azimuth (rad)
* `hubPos`:     current global hubPos (x,y,z) vector (m)
* `nacPos`:     current global nacPos (x,y,z) vector (m)
* `hubAngle`:   3 angle set for hub orientation (rad) , no rotation from spinning
"""
function getAD15MeshPos(turbine,u_j,azi,nacPos,hubPos,hubAngle)
    #display(turbine.bladeIdx)
    MeshPos     = zeros(Float32,3,turbine.numMeshNodes)
    iNode = 1
    # conversion from hub coordinates to global
    if turbine.isHAWT
        CG2H = calcHubRotMat(hubAngle, azi;rot_axis = 1)
        CG2H = createGeneralTransformationMatrix([90],[2]) * CG2H
    else
        CG2H = calcHubRotMat(hubAngle, azi;rot_axis = 3)
    end
    CH2G = CG2H
    # blades, bottom struts, top struts
    for ibld = 1:size(turbine.bladeIdx,1)
        npts = turbine.bladeIdx[ibld,2] - turbine.bladeIdx[ibld,1]
        sgn = 1*sign(npts)
        for idx=turbine.bladeIdx[ibld,1]:sgn:turbine.bladeIdx[ibld,2]
            x=turbine.Mesh.x[idx] + u_j[(idx-1)*6+1]
            y=turbine.Mesh.y[idx] + u_j[(idx-1)*6+2]
            z=turbine.Mesh.z[idx] + u_j[(idx-1)*6+3]

            MeshPos[:,iNode] = [x y z] * CH2G + nacPos' #+ hubPos'
            iNode += 1
        end
    end
    return MeshPos
end

"""
getAD15MeshVelAcc(turbine,meshPos,udot_j,uddot_j,azi,Omega_rad,OmegaDot_rad,nacPos,hubPos,hubAngle,hubVel,hubAcc)

Extract the mesh velocities and accelerations for all the AD15 blades
  ordering here is important
      1. root to tip of blades,        in blade order
      2. root to tip of bottom struts, in blade order
      3. root to tip of top    struts, in blade order

# Inputs
* `turbine`:        turbine data storage
* `rootPos`:        root positions from call to getAD15MeshPos
* `u_j`:            mesh displacements -- in hub coordinates, (m,rad)
* `udot_j`:         mesh velocity      -- in hub coordinates, (m/s,rad/s)
* `uddot_j`:        mesh velocity      -- in hub coordinates, (m/s,rad/s)
* `azi`:            current azimuth (rad)
* `Omega_rad`:      angular velocity of hub about hub axis (rad/s)
* `OmegaDot_rad`:   angular acceleration of hub about hub axis (rad/s^2)
* `hubPos`:         current global hubPos (x,y,z) vector (m)
* `nacPos`:         current global nacPos (x,y,z) vector (m)
* `hubAngle`:       3 angle set for hub orientation (rad), no rotation from spinning
* `hubVel`:         hub velocity in global coords (m/s,rad/s), does not include rotational velocity, so this is extra like from a platform
* `hubAcc`:         hub acceleration in global coords (m/s^2,rad/s^2), does not include rotational accel, so this is extra like from a platform
"""
function getAD15MeshVelAcc(turbine,meshPos,udot_j,uddot_j,azi,Omega_rad,OmegaDot_rad,nacPos,hubPos,hubAngle,hubVel,hubAcc)
    #display(turbine.bladeIdx)
    MeshVel     = zeros(Float32,6,turbine.numMeshNodes)
    MeshAcc     = zeros(Float32,6,turbine.numMeshNodes)
    # conversion from hub coordinates to global
    if turbine.isHAWT
        CG2H = calcHubRotMat(hubAngle, azi;rot_axis = 1)
    else
        CG2H = calcHubRotMat(hubAngle, azi;rot_axis = 3)
    end
    # end
    CH2G = CG2H 

    # hub axis vector in global coordinates
    if turbine.isHAWT
        hubAxis = CG2H[:,1]
    else
        hubAxis = CG2H[:,3]
    end
    # blades, bottom struts, top struts
    iNode = 1
    for ibld = 1:size(turbine.bladeIdx,1)
        npts = turbine.bladeIdx[ibld,2] - turbine.bladeIdx[ibld,1]
        sgn = 1*sign(npts)
        for tmpIdx=turbine.bladeIdx[ibld,1]:sgn:turbine.bladeIdx[ibld,2]
            idx=(tmpIdx-1)*6   # just before the node of interest
            ### 1. relative velocity from mesh distortions
            MeshVel[1:3,iNode] =  udot_j[idx+1:idx+3]' * CH2G     # translation Vel (m/2)
            MeshVel[4:6,iNode] =  udot_j[idx+4:idx+6]' * CH2G     # rotation    Vel (rad/s)
            MeshAcc[1:3,iNode] = uddot_j[idx+1:idx+3]' * CH2G     # translation Acc (m/s^2)
            MeshAcc[4:6,iNode] = uddot_j[idx+4:idx+6]' * CH2G     # rotation    Acc (rad/s^2)

            #FIXME: missing tangential velocity components due to hub rotational velocity not about hub axis
            ### 2. Tangential velocity due to hub rotation
            # tangential velocity and acceleration, based on distance to hub axis
            # TODO: meshPos-nacPos may be incorrect for the struts - velocity needs to be in the global FOR, and aerodyn converts via the orientations to go into blade orientation (STXv)
            Vel = cross(    Omega_rad*hubAxis, (meshPos[1:3,iNode]-nacPos)) / norm(hubAxis)
            Acc = cross( OmegaDot_rad*hubAxis, (meshPos[1:3,iNode]-nacPos)) / norm(hubAxis)
            MeshVel[1:3,iNode] = MeshVel[1:3,iNode] + Vel
            MeshVel[4:6,iNode] = MeshVel[4:6,iNode] + Omega_rad*hubAxis
            MeshAcc[1:3,iNode] = MeshAcc[1:3,iNode] + Acc
            MeshAcc[4:6,iNode] = MeshAcc[4:6,iNode] + OmegaDot_rad*hubAxis

            ### 3. add in contributions from hub motion in global coordinates
            MeshVel[1:3,iNode] = MeshVel[1:3,iNode] + hubVel[1:3]
            MeshAcc[1:3,iNode] = MeshAcc[1:3,iNode] + hubAcc[1:3]

            # increment counter for next point
            iNode += 1
        end
    end
    #println("           Max vel mesh:   ",maximum(MeshVel),"  at ",Omega_rad," rad/s")
    return MeshVel,MeshAcc
end



"""
getRootDCM(turbine,u_j,azi,hubAngle)

Note on angles
    The GyricFEA mesh in OWENS uses +X as the blade/strut long axis.  In AeroDyn, the blade axis is +Z.  So for transforming the
    blades from GyricFEA mesh coordinates, first rotate by -90 degrees about Y, then do the 3,2,1 coordinate transforms with
    (Twist,Theta,Psi) = (Rz,Ry,Rx) = (Yaw,Pitch,Roll)
        Psi_d   -- rotation about Z axis    -- Yaw (Rz)     -- degrees
        Theta_d -- rotation about Y axis    -- Pitch (Ry)   -- degrees
        Twist_d -- rotation about X axis    -- Roll (Rx)    -- degrees
    The rotation sequence is Roll --> Pitch --> Yaw.  In rotation matrix form, it is R = Rz*Ry*Rx (a [3,2,1] matrix order).

# Inputs
* `turbine`:    turbine data storage
* `u_j`:        mesh displacements -- in hub coordinates, (m,rad)
* `azi`:        current azimuth (rad)
* `hubAngle`:   3 angle set for hub orientation (rad), no rotation from spinning

#FIXME: add averaging of orientations to get nodes within blade/strut
"""
function getRootDCM(turbine,u_j,azi,hubAngle)
    RootOrient  = zeros(9,turbine.adi_numbl)
    # conversion from hub coordinates to global
    if turbine.isHAWT
        CG2H = calcHubRotMat(hubAngle, azi;rot_axis = 1)
    else
        CG2H = calcHubRotMat(hubAngle, azi;rot_axis = 3)
    end
    CH2G = CG2H'

    for i=1:size(turbine.bladeElem,1)
        idx=turbine.bladeElem[i]    #,1]
        Psi     = turbine.Ort.Psi_d[idx]    #- rad2deg(u_j[(idx-1)*6+4]) # we should not include deformation of the angles since those are in the mesh, this is just aligning the overall blade angles, there is deformation in the overall position so the blades follow the deformation of the tower.
        Theta   = turbine.Ort.Theta_d[idx]  #- rad2deg(u_j[(idx-1)*6+5])
        Twist   = turbine.Ort.Twist_d[idx]  #- rad2deg(u_j[(idx-1)*6+6])

        if turbine.isHAWT
            # CW is for the standard clockwise rotation of a HAWT (when standing in front of it looking towards the rotor along
            # +X global direction)
            # CCW has not been developed for the HAWT
            angle_axes = [2,3,1]
            #if turbine.rotateCCW
            #ang1 = [180+Theta,180+Twist,Psi]    # CCW
            #else
            ang1 = [180+Theta,Twist,Psi]        # CW
            #end
        else
            # blades
            #   Y is always towards trailing edge in both OWENS and AD15.
            #   X in OWENS is always outward
            #   AD15 CCW, AD15 blade root is at top with +Z pointing downwards along span
            #   AD15 CW,  AD15 blade root is at bottom with +Z upwards along span
            angle_axes = [2,1,2,3]
            #if turbine.rotateCCW
            ang1 = [-90,Twist*0.0,-90.0,Psi]      # CCW
            #else
#FIXME:CW for clockwise, the blade root will be at the bottom of the blade instead of at the top, so Z is upwards and Y is to
#trailing edge.  New logic is needed here to setup the blade roots correctly.  I don't have time right now to do that.
            #end

            # struts
            #   Y is always towards trailing edge in both OWENS and AD15
            #   OWENS always has Z point towards hub. AD15 always has Z point away from hub.
            #if turbine.rotateCCW
            ang2 = [90,Twist,Theta,Psi]     # CCW
            #else
            #ang2 = [90,Twist+180,Theta,Psi] # CW -- FIXME:CW this has not been tested!!!!
            #end
        end

        if i<=turbine.B
            #FIME: the following is for a CCW spinning rotor.  some things need changing for a CW spinning rotor.
            # flip +z towards X, then apply Twist (Roll, Rx) -> Theta (Pitch, Ry) -> Psi (Yaw, Rz) 
            DCM = CH2G * createGeneralTransformationMatrix(ang1,angle_axes)'
        else
            DCM = CH2G * createGeneralTransformationMatrix(ang2,angle_axes)'
        end

        # if turbine.isHAWT
        #     DCM = DCM * createGeneralTransformationMatrix([90.0],[2])
        # end

        RootOrient[:,i] = vec(DCM)
    end
    return RootOrient
end



"""
getAD15MeshDCM(turbine,u_j,azi,hubAngle)

Extract the mesh points orientations for all the AD15 blades
  ordering here is important
      1. root to tip of blades,        in blade order
      2. root to tip of bottom struts, in blade order
      3. root to tip of top    struts, in blade order

# Inputs
* `turbine`:    turbine data storage
* `u_j`:        mesh displacements -- in hub coordinates, (m,rad)
* `azi`:        current azimuth (rad)
* `hubAngle`:   3 angle set for hub orientation (rad), no rotation from spinning

#FIXME: add averaging of orientations to get nodes within blade/strut
"""
function getAD15MeshDCM(turbine,u_j,azi,hubAngle)
    #display(turbine.bladeElem)
    MeshOrient  = zeros(9,turbine.numMeshNodes)
    # conversion from hub coordinates to global
    if turbine.isHAWT
        CG2H = calcHubRotMat(hubAngle, azi;rot_axis = 1)
    else
        CG2H = calcHubRotMat(hubAngle, azi;rot_axis = 3)
    end
    CH2G = CG2H'
    iNode=0
    for ibld=1:size(turbine.bladeElem,1)
        lenbld = turbine.bladeElem[ibld,2] - turbine.bladeElem[ibld,1]
        sgn = 1*sign(lenbld)
        if sgn > 0      # normal ordering
            idxarry=collect(turbine.bladeElem[ibld,1]:turbine.bladeElem[ibld,2])
        else
            idxarry=reverse(collect(turbine.bladeElem[ibld,2]:turbine.bladeElem[ibld,1]))
        end
            #println("+ direction blade: $(turbine.bladeElem[ibld,1]):$(turbine.bladeElem[ibld,2])")
        for idx in idxarry

            Psi     = turbine.Ort.Psi_d[idx]    - rad2deg(u_j[(idx-1)*6+4])
            Theta   = turbine.Ort.Theta_d[idx]  - rad2deg(u_j[(idx-1)*6+5])
            Twist   = turbine.Ort.Twist_d[idx]  - rad2deg(u_j[(idx-1)*6+6])

            if turbine.isHAWT
                angle_axes = [2,1,2,3,2]
                ang = [90,Twist,Theta,Psi,90]
            else
                # The OWENS mesh for VAWT is always setup with X pointing away from axis of rotation, and Y towards the trailing edge
                angle_axes = [2,1,2,3]
                ang = [-90,Twist,Theta,Psi]
            end

            DCM =  CH2G * createGeneralTransformationMatrix(ang,angle_axes)'

            iNode += 1
            MeshOrient[:,iNode] = vec(DCM)
            # duplicate last node #TODO: map elements to nodes instead
            if idx==turbine.bladeElem[ibld,2]
                iNode += 1
                MeshOrient[:,iNode] = vec(DCM)
            end
        end     
    end
    return MeshOrient
end


### NOTE: the following are copied from OWENS.jl/src/meshing_utilities.jl
#          (I couldn't figure out how to use OWENS here without a mess)
"""
createGeneralTransformationMatrix(angleArray,axisArray)

Calculates the transformation matrix assocaited with a general Euler rotation sequence.

#Input
* `angleArray`:      = array of angles for Euler rotation sequence
* `axisArray`:       = array of axis of rotations for Euler rotation

#Output
* `dcmTotal`:        = transformation matrix of specified euler rotation sequence
"""
function createGeneralTransformationMatrix(angleArray,axisArray)

    numRotations = length(angleArray) #get number of rotation to perform
    dcmArray = zeros(3,3,numRotations) #initialize individual rotation direction cosine matrix arrays

    for i=1:numRotations #calculate individual single rotatio direction cosine matrices
        dcmArray[:,:,i] = createSingleRotationDCM(angleArray[i],axisArray[i])
    end

    dcmTotal = dcmArray[:,:,1] #initialize dcmTotal as first rotation

    #multiply consecutive rotation sequency direction cosine matrices to arrive at overall transformation matrix
    for i=2:1:numRotations
        dcmTotal = dcmTotal*dcmArray[:,:,i]
    end

    return dcmTotal
end

"""
createSingleRotationDCM(angleDeg,axisNum)

Creates a direction cosine matrix (dcm) associated with a rotation of angleDeg about axisNum.
"""
function createSingleRotationDCM(angleDeg,axisNum)

    angleRad = angleDeg*pi/180.0 #convert angle to radians

    if axisNum == 1 #direction cosine matrix about 1 axis
        dcm = [1.0 0.0 0.0
        0.0 cos(angleRad) sin(angleRad)
        0.0 -sin(angleRad) cos(angleRad)]
    elseif axisNum == 2 #direction cosine matrix about 2 axis
        dcm = [cos(angleRad) 0.0 -sin(angleRad)
        0.0 1.0 0.0
        sin(angleRad) 0.0 cos(angleRad)]
    elseif axisNum == 3 #direction cosine matrix about 3 axis
        dcm = [cos(angleRad) sin(angleRad) 0.0
        -sin(angleRad) cos(angleRad) 0.0
        0.0 0.0 1.0]
    else  #error catch
        error("Error: createSingleRotationDCM. Axis number must be 1, 2, or 3.")
    end

    return dcm

end



"""
transMat(theta1, theta2, theta3)

Internal, computes the 3x3 transformation matrix for given input rotations. The generated matrix
is the closest orthonormal matrix to the Bernoulli-Euler transformation matrix from beam theory,
which assumes small rotations. A full description of this matrix is found in the
"FASTCoordinateSystems.doc" document by Jason Jonkman.
"""
function transMat(theta1, theta2, theta3)
    theta11      = theta1 * theta1
    theta22      = theta2 * theta2
    theta33      = theta3 * theta3

    sqrdSum      = theta11 + theta22 + theta33
    sqrt1sqrdSum = sqrt( 1.0 + sqrdSum )
    comDenom     = sqrdSum * sqrt1sqrdSum

    theta12S     = theta1 * theta2 * ( sqrt1sqrdSum - 1.0 )
    theta13S     = theta1 * theta3 * ( sqrt1sqrdSum - 1.0 )
    theta23S     = theta2 * theta3 * ( sqrt1sqrdSum - 1.0 )


    # Define the transformation matrix:
    transMat = Array{Float32}(undef, 3,3)
    if comDenom == 0.0  # All angles are zero and matrix is ill-conditioned (the matrix is derived assuming that the angles are not zero); return identity

        transMat = [1 0 0; 0 1 0; 0 0 1]

    else  # At least one angle is nonzero

        transMat[1,1] = ( theta11*sqrt1sqrdSum + theta22              + theta33              ) / comDenom
        transMat[2,2] = ( theta11              + theta22*sqrt1sqrdSum + theta33              ) / comDenom
        transMat[3,3] = ( theta11              + theta22              + theta33*sqrt1sqrdSum ) / comDenom
        transMat[1,2] = (  theta3*sqrdSum + theta12S ) / comDenom
        transMat[2,1] = ( -theta3*sqrdSum + theta12S ) / comDenom
        transMat[1,3] = ( -theta2*sqrdSum + theta13S ) / comDenom
        transMat[3,1] = (  theta2*sqrdSum + theta13S ) / comDenom
        transMat[2,3] = (  theta1*sqrdSum + theta23S ) / comDenom
        transMat[3,2] = ( -theta1*sqrdSum + theta23S ) / comDenom

    end

    return transMat

end


"""
frame_convert(init_frame_vals, trans_mat)

Internal, transfers 6 DOFs element-wise to a new reference frame

# Input
* `init_frame_vals::Vector{<:float}`: Values in 6 degrees of freedom in the initial reference frame
* `trans_mat::Array{<:float}`: Transformation matrix to the output reference frame

# Output
* `out_frame_vals`: Values in 6 degrees of freedom in the output reference frame
"""
function frame_convert(init_frame_vals, trans_mat)
    out_frame_vals = copy(init_frame_vals).*0.0
    out_frame_vals[1:3] = trans_mat * init_frame_vals[1:3]
    out_frame_vals[4:6] = trans_mat * init_frame_vals[4:6]
    return out_frame_vals
end

function calcHubRotMat(ptfmRot, azi_j;rot_axis = 3)
    # CN2P = transMat(ptfmRot[1], ptfmRot[2], ptfmRot[3])
    CP2H = createGeneralTransformationMatrix([azi_j*180/pi],[rot_axis])
    CN2P = createGeneralTransformationMatrix(ptfmRot,[1,2,3])
    # Since we are modeling the rotor in the VAWT frame of reference, rotate about the azimuth, then rotate the hub 
    CN2H = CP2H*CN2P
    return CN2H
end


#ADP: I would not include a hard coded file like this.  This could be a maintenance nightmare as the AD15 input file changes very frequently (almost every OpenFAST release has a change in this file)
function writeADinputFile(filename="ADInput.dat",blade_filenames=nothing,airfoil_filenames=nothing,OLAF_filename=nothing)
    
    OLAF_str = "\"$OLAF_filename\" OLAFInputFileName - Input file for OLAF [used only when WakeMod=3]"

    #TODO: multiple airfoils, other relevant inputs
    AF_str ="\"$airfoil_filenames\"    AFNames            - Airfoil file names (NumAFfiles lines) (quoted strings)"

    blades_str = ""
    for ibld = 1:length(blade_filenames) #includes struts since they are modeled as b
        if ibld == length(blade_filenames)
            blades_str =blades_str*"$(blade_filenames[ibld]) ADBlFile($ibld) - Name of file containing distributed aerodynamic properties for Blade #1 (-)"
        else
            blades_str =blades_str*"$(blade_filenames[ibld]) ADBlFile($ibld) - Name of file containing distributed aerodynamic properties for Blade #1 (-) \n"
        end

    end

    input_string_array = 
"------- AERODYN v15.03.* INPUT FILE ------------------------------------------------
Generated with OWENS driver
======  General Options  ============================================================================
True                   Echo        - Echo the input to \"<rootname>.AD.ech\"?  (flag)
\"default\"              DTAero      - Time interval for aerodynamic calculations {or \"default\"} (s)
3                      WakeMod     - Type of wake/induction model (switch) {0=none, 1=BEMT, 2=DBEMT, 3=OLAF}
#2                      AFAeroMod     - Type of blade airfoil aerodynamics model (switch) {1=steady model, 2=Beddoes-Leishman unsteady model} [must be 1 when linearizing]
1                      AFAeroMod     - Type of blade airfoil aerodynamics model (switch) {1=steady model, 2=Beddoes-Leishman unsteady model} [must be 1 when linearizing]
0                      TwrPotent   - Type tower influence on wind based on potential flow around the tower (switch) {0=none, 1=baseline potential flow, 2=potential flow with Bak correction}
0                      TwrShadow   - Calculate tower influence on wind based on downstream tower shadow? (flag)
False                  TwrAero     - Calculate tower aerodynamic loads? (flag)
False                  FrozenWake  - Assume frozen wake during linearization? (flag) [used only when WakeMod=1 and when linearizing]
False                  CavitCheck  - Perform cavitation check? (flag) TRUE will turn off unsteady aerodynamics
False                  Buoyancy           - Include buoyancy effects? (flag)
False                  CompAA      - Flag to compute AeroAcoustics calculation [only used when WakeMod=1 or 2]
\"AeroAcousticsInput.dat\" AA_InputFile - Aeroacoustics input file
======  Environmental Conditions  ===================================================================
1.225                   AirDens     - Air density (kg/m^3)
1.4639E-05              KinVisc     - Kinematic air viscosity (m^2/s)
3.350000000000000e+02 SpdSound    - Speed of sound (m/s)
1.035000000000000e+05 Patm        - Atmospheric pressure (Pa) [used only when CavitCheck=True]
1.700000000000000e+03 Pvap        - Vapour pressure of fluid (Pa) [used only when CavitCheck=True]
======  Blade-Element/Momentum Theory Options  ====================================================== [used only when WakeMod=1]
2                      SkewMod     - Type of skewed-wake correction model (switch) {1=uncoupled, 2=Pitt/Peters, 3=coupled} [used only when WakeMod=1]
\"default\"              SkewModFactor - Constant used in Pitt/Peters skewed wake model {or \"default\" is 15/32*pi} (-) [used only when SkewMod=2; unused when WakeMod=0]
True                   TipLoss     - Use the Prandtl tip-loss model? (flag) [used only when WakeMod=1]
True                   HubLoss     - Use the Prandtl hub-loss model? (flag) [used only when WakeMod=1]
True                   TanInd      - Include tangential induction in BEMT calculations? (flag) [used only when WakeMod=1]
True                   AIDrag      - Include the drag term in the axial-induction calculation? (flag) [used only when WakeMod=1]
True                   TIDrag      - Include the drag term in the tangential-induction calculation? (flag) [used only when WakeMod=1 and TanInd=TRUE]
\"Default\"              IndToler    - Convergence tolerance for BEMT nonlinear solve residual equation {or \"default\"} (-) [used only when WakeMod=1]
100                    MaxIter     - Maximum number of iteration steps (-) [used only when WakeMod=1]
======  Dynamic Blade-Element/Momentum Theory Options  ====================================================== [used only when WakeMod=1]
1                      DBEMT_Mod   - Type of dynamic BEMT (DBEMT) model {1=constant tau1, 2=time-dependent tau1} (-) [used only when WakeMod=2]
20                     tau1_const  - Time constant for DBEMT (s) [used only when WakeMod=2 and DBEMT_Mod=1]
======  OLAF -- cOnvecting LAgrangian Filaments (Free Vortex Wake) Theory Options  ================== [used only when WakeMod=3]
$OLAF_str
#\"OLAF-prescr.dat\"               OLAFInputFileName - Input file for OLAF [used only when WakeMod=3]
======  Beddoes-Leishman Unsteady Airfoil Aerodynamics Options  ===================================== [used only when AFAeroMod=2]
2                      UAMod       - Unsteady Aero Model Switch (switch) {1=Baseline model (Original), 2=Gonzalez's variant (changes in Cn,Cc,Cm), 3=Minemma/Pierce variant (changes in Cc and Cm)} [used only when AFAeroMod=2]
True                   FLookup     - Flag to indicate whether a lookup for f' will be calculated (TRUE) or whether best-fit exponential equations will be used (FALSE); if FALSE S1-S4 must be provided in airfoil input files (flag) [used only when AFAeroMod=2]
======  Airfoil Information =========================================================================
2                      AFTabMod    - Interpolation method for multiple airfoil tables {1=1D interpolation on AoA (first table only); 2=2D interpolation on AoA and Re; 3=2D interpolation on AoA and UserProp} (-)
1                      InCol_Alfa  - The column in the airfoil tables that contains the angle of attack (-)
2                      InCol_Cl    - The column in the airfoil tables that contains the lift coefficient (-)
3                      InCol_Cd    - The column in the airfoil tables that contains the drag coefficient (-)
4                      InCol_Cm    - The column in the airfoil tables that contains the pitching-moment coefficient; use zero if there is no Cm column (-)
0                      InCol_Cpmin - The column in the airfoil tables that contains the Cpmin coefficient; use zero if there is no Cpmin column (-)
1                      NumAFfiles  - Number of airfoil files used (-)
$AF_str
======  Rotor/Blade Properties  =====================================================================
True                   UseBlCm     - Include aerodynamic pitching moment in calculations?  (flag)
$blades_str
======  Hub Properties ============================================================================== [used only when Buoyancy=True]
0.0   VolHub             - Hub volume (m^3)
0.0   HubCenBx           - Hub center of buoyancy x direction offset (m)
======  Nacelle Properties ========================================================================== [used only when Buoyancy=True]
0.0   VolNac             - Nacelle volume (m^3)
0,0,0 NacCenB            - Position of nacelle center of buoyancy from yaw bearing in nacelle coordinates (m)
======  Tail fin Aerodynamics ======================================================================== 
False         TFinAero           - Calculate tail fin aerodynamics model (flag)
\"unused\"      TFinFile           - Input file for tail fin aerodynamics [used only when TFinAero=True]
======  Tower Influence and Aerodynamics ============================================================ [used only when TwrPotent/=0, TwrShadow/=0, TwrAero=True, or Buoyancy=True]
        3    NumTwrNds         - Number of tower nodes used in the analysis  (-) [used only when TwrPotent/=0, TwrShadow/=0, TwrAero=True, or Buoyancy=True]
TwrElev        TwrDiam        TwrCd          TwrTI          TwrCb !TwrTI used only with TwrShadow=2, TwrCb used only with Buoyancy=True
(m)              (m)           (-)            (-)           (-)
0.00    0.60e+00  1.000e+00  0.1  0.0
5.00    0.60e+00  1.000e+00  0.1  0.0
10.00    0.60e+00  1.000e+00  0.1  0.0
======  Outputs  ====================================================================================
True                    SumPrint    - Generate a summary file listing input options and interpolated properties to \"<rootname>.AD.sum\"?  (flag)
0                      NBlOuts     - Number of blade node outputs [0 - 9] (-)
1, 4, 7, 13, 18, 23, 26, 28, 30 BlOutNd     - Blade nodes whose values will be output  (-)
0                      NTwOuts     - Number of tower node outputs [0 - 9]  (-)
1, 2, 6                TwOutNd     - Tower nodes whose values will be output  (-)
                OutList             - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)
\"RtAeroCp\"
\"RtAeroCq\"
\"RtAeroCt\"
\"RtAeroFxh\"
\"RtAeroFyh\"
\"RtAeroFzh\"
\"RtAeroMxh\"
\"RtAeroMyh\"
\"RtAeroMzh\"
\"RtAeroPwr\"
\"RtAeroFxi\"
\"RtAeroFyi\"
\"RtAeroFzi\"
\"RtAeroMxi\"
\"RtAeroMyi\"
\"RtAeroMzi\"
\"RtArea\"
\"RtSkew\"
\"RtSpeed\"
\"RtTSR\"
\"RtVAvgxh\"
\"RtVAvgyh\"
\"RtVAvgzh\"
\"B1Pitch\"
\"B2Pitch\"
\"B3Pitch\"
\"B1Azimuth\"
\"B2Azimuth\"
\"B3Azimuth\"
\"B1AeroFx\"
\"B1AeroFy\"
\"B1AeroFz\"
\"B1AeroMx\"
\"B1AeroMy\"
\"B1AeroMz\"
\"B1AeroPwr\"
\"B1AeroFxi\"
\"B1AeroFyi\"
\"B1AeroFzi\"
\"B1AeroMxi\"
\"B1AeroMyi\"
\"B1AeroMzi\"
\"B2AeroFx\"
\"B2AeroFy\"
\"B2AeroFz\"
\"B2AeroMx\"
\"B2AeroMy\"
\"B2AeroMz\"
\"B2AeroPwr\"
\"B2AeroFxi\"
\"B2AeroFyi\"
\"B2AeroFzi\"
\"B2AeroMxi\"
\"B2AeroMyi\"
\"B2AeroMzi\"
\"B3AeroFx\"
\"B3AeroFy\"
\"B3AeroFz\"
\"B3AeroMx\"
\"B3AeroMy\"
\"B3AeroMz\"
\"B3AeroPwr\"
\"B3AeroFxi\"
\"B3AeroFyi\"
\"B3AeroFzi\"
\"B3AeroMxi\"
\"B3AeroMyi\"
\"B3AeroMzi\"
\"B4AeroFx\"
\"B4AeroFy\"
\"B4AeroFz\"
\"B4AeroMx\"
\"B4AeroMy\"
\"B4AeroMz\"
\"B4AeroPwr\"
\"B4AeroFxi\"
\"B4AeroFyi\"
\"B4AeroFzi\"
\"B4AeroMxi\"
\"B4AeroMyi\"
\"B4AeroMzi\"
END of input file (the word \"END\" must appear in the first 3 columns of this last OutList line)
====== Outputs for all blade stations (same ending as above for B1N1.... =========================== [optional section]
9              BldNd_BladesOut     - Number of blades to output all node information at.  Up to number of blades on turbine. (-)
\"All\"          BldNd_BlOutNd       - Future feature will allow selecting a portion of the nodes to output.  Not implemented yet. (-)
                OutListAD             - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)
Alpha
Cl
Cd
Cm
Fx
Fy
Vx
Vy
VUndx
VUndy
VUndz
VDisx
VDisy
VDisz
STVx
STVy
STVz
Vrel
Vindx
Vindy
Gam
Fn
Ft
Cx
Cy
Cn
Ct
Fl
Fd
Mm
END of input file (the word \"END\" must appear in the first 3 columns of this last OutList line)
---------------------------------------------------------------------------------------
"
    
    # ifw_input_string = join(input_string_array, "\0")
    
    # Save the file
    open(filename, "w") do file
        write(file, input_string_array)
    end

    return input_string_array
end

function writeADbladeFile(filename="./blade.dat";NumBlNds=10,BlSpn=collect(LinRange(0,12,NumBlNds)),
    BlCrvAC=zeros(NumBlNds),BlSwpAC=zeros(NumBlNds),BlCrvAng=zeros(NumBlNds),
    BlTwist=zeros(NumBlNds),BlChord=ones(NumBlNds).*0.5,BlAFID=ones(Int,NumBlNds))

    data_str = ""
    for inode = 1:NumBlNds #includes struts since they are modeled as blades
        data_str =data_str*"$(BlSpn[inode]) $(BlCrvAC[inode]) $(BlSwpAC[inode]) $(BlCrvAng[inode]) $(BlTwist[inode]) $(BlChord[inode]) $(BlAFID[inode])	\n"
    end

    input_string_array = 
"------- AERODYN v15.00.* BLADE DEFINITION INPUT FILE -------------------------------------						
Generated with OWENS driver						
======  Blade Properties =================================================================						
$NumBlNds   NumBlNds    - Number of blade nodes used in the analysis (-)
    BlSpn   BlCrvAC  BlSwpAC   BlCrvAng  BlTwist   BlChord   BlAFID						
    (m)      (m)      (m)      (deg)      (deg)      (m)      (-)						
$data_str
"
    
    # ifw_input_string = join(input_string_array, "\0")
    
    # Save the file
    open(filename, "w") do file
        write(file, input_string_array)
    end

    return input_string_array
end

function writeOLAFfile(filename="./OLAF.dat";nNWPanel=200,nFWPanels=10)
    input_string_array = 
"--------------------------- FREE WAKE INPUT FILE ---------------------------------------------- 
Free wake input file for the turbine
--------------------------- GENERAL OPTIONS ---------------------------------------------------
5             IntMethod     - Integration method {1: RK4, 5: Forward Euler 1st order, default: 5} (switch)
default       DTfvw         - Time interval for wake propagation. {default: dtaero} (s)
0.0           FreeWakeStart - Time when wake is free. (-) value = always free. {default: 0.0} (s)
0.1           FullCircStart - Time at which full circulation is reached. {default: 0.0} (s)
--------------------------- CIRCULATION SPECIFICATIONS ----------------------------------------
1             CircSolvingMethod - Circulation solving method {1: Cl-Based, 2: No-Flow Through, 3: Prescribed, default: 1 }(switch)
default       CircSolvConvCrit - Convergence criteria {default: 0.001} [only if CircSolvingMethod=1] (-)
default       CircSolvRelaxation - Relaxation factor {default: 0.1} [only if CircSolvingMethod=1] (-)
default       CircSolvMaxIter - Maximum number of iterations for circulation solving {default: 30} (-)
\"GammaPrescr.csv\" PrescribedCircFile - File containing prescribed circulation [only if CircSolvingMethod=3] (quoted string)
===============================================================================================
--------------------------- WAKE OPTIONS ------------------------------------------------------
------------------- WAKE EXTENT AND DISCRETIZATION --------------------------------------------
$nNWPanel           nNWPanel      - Number of near-wake panels (-)
default nNWPanelsFree      - Number of free near-wake panels (-) {default: nNWPanels}
$nFWPanels      nFWPanels          - Number of far-wake panels (-) {default: 0}
default nFWPanelsFree      - Number of free far-wake panels (-) {default: nFWPanels}
default FWShedVorticity    - Include shed vorticity in the far wake {default: False}
------------------- WAKE REGULARIZATIONS AND DIFFUSION -----------------------------------------
0             DiffusionMethod - Diffusion method to account for viscous effects {0: None, 1: Core Spreading, \"default\": 0}
3             RegDeterMethod - Method to determine the regularization parameters {0: Manual, 1: Optimized, 2: Chord, 3: Span, default: 0 }
1             RegFunction   - Viscous diffusion function {0: None, 1: Rankine, 2: LambOseen, 3: Vatistas, 4: Denominator, \"default\": 3} (switch)
3             WakeRegMethod - Wake regularization method {1: Constant, 2: Stretching, 3: Age, default: 1} (switch)
1.0           WakeRegFactor - Wake regularization factor (m or -)
1.0           WingRegFactor - Wing regularization factor (m or -)
2000          CoreSpreadEddyVisc - Eddy viscosity in core spreading methods, typical values 1-1000
------------------- WAKE TREATMENT OPTIONS ---------------------------------------------------
True          TwrShadowOnWake - Include tower flow disturbance effects on wake convection {default:false} [only if TwrPotent or TwrShadow]
0             ShearModel    - Shear Model {0: No treatment, 1: Mirrored vorticity, default: 0}
------------------- SPEEDUP OPTIONS -----------------------------------------------------------
1             VelocityMethod - Method to determine the velocity {1:Biot-Savart Segment, 2:Particle tree, default: 1}
1.5           TreeBranchFactor - Branch radius fraction above which a multipole calculation is used {default: 2.0} [only if VelocityMethod=2]
1             PartPerSegment - Number of particles per segment [only if VelocityMethod=2]
===============================================================================================
--------------------------- OUTPUT OPTIONS  ---------------------------------------------------
1             WrVTk         - Outputs Visualization Toolkit (VTK) (independent of .fst option) {False: NoVTK, True: Write VTK at each time step} (flag)
3             nVTKBlades    - Number of blades for which VTK files are exported {0: No VTK per blade, n: VTK for blade 1 to n} (-)
1             VTKCoord      - Coordinate system used for VTK export. {1: Global, 2: Hub, \"default\": 1}
\"default\"     VTK_fps       - Frame rate for VTK output (frames per second) {\"all\" for all glue code timesteps, \"default\" for all FVW timesteps} [only if WrVTK=1]
0             nGridOut      - Number of grid outputs
GridName	DTOut         - XStart XEnd nX YStart YEnd nY ZStart ZEnd nZ
(-)  		(s)           - (m) (m) (-) (m) (m) (-) (m) (m) (-)
\"hori\"    default      -10       30.   400     -70.      70.    150   20   20  1
\"vert\"    default      -10       30.   400     -0.      0.     1     0    80  100 
\"vert2\"   default      -10       30.   400     5.      5.     1     5.136    15.136  100 
===============================================================================================
"
    
    # ifw_input_string = join(input_string_array, "\0")
    
    # Save the file
    open(filename, "w") do file
        write(file, input_string_array)
    end

    return input_string_array
end

function writeIWfile(Vinf,filename = "./AD15-input/IW_test.dat";WindType=1,windINPfilename=nothing)
    HWindSpeed_str = "$(round(Vinf,digits=6))   HWindSpeed     - Horizontal windspeed                            (m/s)"
    turbsim_str = "\"$windINPfilename\"      filename_bts   - name of the full field wind file to use (.bts)"
    uniformWind_str = "\"$windINPfilename\"      FileName_Uni   - Filename of time series data for uniform wind field.      (-)"

    input_string_array = 
"------- InflowWind INPUT FILE -------------------------------------------------------------------------
Input 
---------------------------------------------------------------------------------------------------------------
False         Echo           - Echo input data to <RootName>.ech (flag)
            $WindType   WindType       - switch for wind file type (1=steady; 2=uniform; 3=binary TurbSim FF; 4=binary Bladed-style FF; 5=HAWC format; 6=User defined; 7=native Bladed FF)
            0   PropagationDir - Direction of wind propagation (meteorological rotation from aligned with X (positive rotates towards -Y) -- degrees) (not used for native Bladed format WindType=7)
            0   VFlowAng       - Upflow angle (degrees) (not used for native Bladed format WindType=7)
        False   VelInterpCubic - Use cubic interpolation for velocity in time (false=linear, true=cubic) [Used with WindType=2,3,4,5,7]
            1   NWindVel       - Number of points to output the wind velocity    (0 to 9)
            0   WindVxiList    - List of coordinates in the inertial X direction (m)
            0   WindVyiList    - List of coordinates in the inertial Y direction (m)
            50   WindVziList    - List of coordinates in the inertial Z direction (m)
================== Parameters for Steady Wind Conditions [used only for WindType = 1] =========================
        $HWindSpeed_str
            50   RefHt          - Reference height for horizontal wind speed      (m)
            0   PLExp          - Power law exponent                              (-)
================== Parameters for Uniform wind file   [used only for WindType = 2] ============================
    $uniformWind_str    FileName_Uni   - Filename of time series data for uniform wind field.      (-)
        100   RefHt_Uni      - Reference height for horizontal wind speed                (m)
        125.88   RefLength      - Reference length for linear horizontal and vertical sheer (-)
================== Parameters for Binary TurbSim Full-Field files   [used only for WindType = 3] ==============
$turbsim_str
================== Parameters for Binary Bladed-style Full-Field files   [used only for WindType = 4 or WindType = 7] =========
\"unused\"      FileNameRoot   - WindType=4: Rootname of the full-field wind file to use (.wnd, .sum); WindType=7: name of the intermediate file with wind scaling values
False         TowerFile      - Have tower file (.twr) (flag) ignored when WindType = 7
================== Parameters for HAWC-format binary files  [Only used with WindType = 5] =====================
\"unused\"      FileName_u     - name of the file containing the u-component fluctuating wind (.bin)
\"unused\"      FileName_v     - name of the file containing the v-component fluctuating wind (.bin)
\"unused\"      FileName_w     - name of the file containing the w-component fluctuating wind (.bin)
            64   nx             - number of grids in the x direction (in the 3 files above) (-)
            32   ny             - number of grids in the y direction (in the 3 files above) (-)
            32   nz             - number of grids in the z direction (in the 3 files above) (-)
            16   dx             - distance (in meters) between points in the x direction    (m)
            3   dy             - distance (in meters) between points in the y direction    (m)
            3   dz             - distance (in meters) between points in the z direction    (m)
        15000   RefHt_Hawc     - reference height; the height (in meters) of the vertical center of the grid (m)
    -------------   Scaling parameters for turbulence   ---------------------------------------------------------
            2   ScaleMethod    - Turbulence scaling method   [0 = none, 1 = direct scaling, 2 = calculate scaling factor based on a desired standard deviation]
            1   SFx            - Turbulence scaling factor for the x direction (-)   [ScaleMethod=1]
            1   SFy            - Turbulence scaling factor for the y direction (-)   [ScaleMethod=1]
            1   SFz            - Turbulence scaling factor for the z direction (-)   [ScaleMethod=1]
        1.2   SigmaFx        - Turbulence standard deviation to calculate scaling from in x direction (m/s)    [ScaleMethod=2]
        0.8   SigmaFy        - Turbulence standard deviation to calculate scaling from in y direction (m/s)    [ScaleMethod=2]
        0.2   SigmaFz        - Turbulence standard deviation to calculate scaling from in z direction (m/s)    [ScaleMethod=2]
    -------------   Mean wind profile parameters (added to HAWC-format files)   ---------------------------------
            12   URef           - Mean u-component wind speed at the reference height (m/s)
            2   WindProfile    - Wind profile type (0=constant;1=logarithmic,2=power law)
        0.2   PLExp_Hawc     - Power law exponent (-) (used for PL wind profile type only)
        0.03   Z0             - Surface roughness length (m) (used for LG wind profile type only)
            0   XOffset         - Initial offset in +x direction (shift of wind box)
================== LIDAR Parameters ===========================================================================
            0   SensorType          - Switch for lidar configuration (0 = None, 1 = Single Point Beam(s), 2 = Continuous, 3 = Pulsed)
            0   NumPulseGate        - Number of lidar measurement gates (used when SensorType = 3)
            30   PulseSpacing        - Distance between range gates (m) (used when SensorType = 3)
            0   NumBeam             - Number of lidar measurement beams (0-5)(used when SensorType = 1)
        -200   FocalDistanceX      - Focal distance co-ordinates of the lidar beam in the x direction (relative to hub height) (only first coordinate used for SensorType 2 and 3) (m)
            0   FocalDistanceY      - Focal distance co-ordinates of the lidar beam in the y direction (relative to hub height) (only first coordinate used for SensorType 2 and 3) (m)
            0   FocalDistanceZ      - Focal distance co-ordinates of the lidar beam in the z direction (relative to hub height) (only first coordinate used for SensorType 2 and 3) (m)
0.0 0.0 0.0   RotorApexOffsetPos  - Offset of the lidar from hub height (m)
            17   URefLid             - Reference average wind speed for the lidar[m/s]
        0.25   MeasurementInterval - Time between each measurement [s]
        False   LidRadialVel        - TRUE => return radial component, FALSE => return 'x' direction estimate
            1   ConsiderHubMotion   - Flag whether to consider the hub motion's impact on Lidar measurements
====================== OUTPUT ==================================================
False         SumPrint     - Print summary data to <RootName>.IfW.sum (flag)
                OutList      - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)
\"Wind1VelX\"               X-direction wind velocity at point WindList(1)
\"Wind1VelY\"               Y-direction wind velocity at point WindList(1)
\"Wind1VelZ\"               Z-direction wind velocity at point WindList(1)
END of input file (the word \"END\" must appear in the first 3 columns of this last OutList line)
---------------------------------------------------------------------------------------"
       
    # ifw_input_string = join(input_string_array, "\0")
    
    # Save the file
    open(filename, "w") do file
        write(file, input_string_array)
    end

    return input_string_array
end
