global adilib
global sym_calcoutput
global sym_updatestates
global sym_end

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

* `wrOuts::int`:        optional, file format for writing outputs [0 none (default), 1 txt, 2 binary, 3 both]
* `DT_Outs::float64`:   optional, timestep for outputs to file [0.0 (default) for every timestep]

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
    wrOuts      = 0,
    DT_Outs     = 0.0,
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
        println("Opening AeroDyn-Inflow library at: $adilib_filename")
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
            Ref{Cdouble},       # IN: dt
            Ref{Cdouble},       # IN: t_max
            Ref{Cint},          # IN: storeHHVel    (c_bool)
            Ref{Cint},          # IN: transposeDCM  (c_bool)
            Ref{Cint},          # IN: WrVTK
            Ref{Cint},          # IN: WrVTK_Type
            Ref{Cfloat},        # IN: VTKNacDim
            Ref{Cfloat},        # IN: VTKHubRad
            Ref{Cint},          # IN: wrOuts
            Ref{Cdouble},       # IN: DT_Outs
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
            Cdouble(dt),
            Cdouble(t_max),
            Cint.(storeHHVel),
            Cint.(transposeDCM),
            WrVTK,
            WrVTK_Type,
            Cfloat.(VTKNacDim),
            VTKHubRad,
            wrOuts,
            Cdouble(DT_Outs),
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

    return num_channels[1], channel_names, channel_units

end

function adiCalcOutput(time, 
                 hubPos, hubOrient, hubVel, hubAcc,
                 nacPos, nacOrient, nacVel, nacAcc,
                 rootPos, rootOrient, rootVel, rootAcc,
                 numMeshNodes,
                 meshPos, meshOrient, meshVel, meshAcc,
                 num_channels)
    meshFrcMom = zeros(Cfloat,6*numMeshNodes);
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
            Ref{Cfloat},        # IN: meshAcc 
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
            Cfloat.(meshAcc),
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
* `B::TI`: Number of blades
* `adi_numbl::TI2`: number of adi blades (includes struts)
* `numMeshNodes::TI3`: number of mesh nodes we are passing motions/loads to/from.  Derived from bladeIdx
* `bladeIdx::TAI1`: an array of start and end nodes for each AD15 blade. Used to index into mesh
* `bladeElem::TAI2`: an array of start and end indices into the orientation angles stored in myort.Psi_d, myort.Theta_d, and myort.Twist_d.
* `Mesh::GyricFEA:Mesh`: mesh without deformation. Storing a copy here to clean up the code a bit.
* `Ort::GyricFEA:Ort`: orientations of all elements without deformation.  Storing here since we don't have a good way of passing this info during simulation
 

# Outputs:
* `none`:

"""
struct Turbine{TF1,TAF1,TI1,TI2,TI3,TAI1,TAI2,TM,TO}
    R::TF1
    omega::TAF1
    B::TI1
    adi_numbl::TI2
    numMeshNodes::TI3
    bladeIdx::TAI1
    bladeElem::TAI2
    Mesh::TM
    Ort::TO
end

#Turbine(R,omega,B,adi_numbl,numMeshNodes,bladeIdx,bladeElem,myort) = Turbine(R,omega,B,adi_numbl,numMeshNodes,bladeIdx,bladeElem,Mesh,Ort)

"""
Structure

Contains the position and orientation info passed to AD15

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
    storeHHVel  = false,
    transposeDCM= true,
    WrVTK       = 2,
    WrVTK_Type  = 3,
    VTKNacDim   = [-.10 ,-.10 ,-.10 ,.2 ,.2 ,.2],
    VTKHubRad   = 1.0,
    adi_nstrut  = 2,
    adi_dt      = 0.05,
    adi_tmax    = 10,
    adi_wrOuts  = 0,
    adi_DT_Outs = 0.0)

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
* `transposeDCM`: this is default -- how data passed across interface is handled
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
* `hubAngle`: hub axis angle, 3-vector (radians)


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
    transposeDCM= true,       # this is default -- how data passed across interface is handled
    WrVTK       = 2,          # write VTK files from adi to directory adi-vtk [0 none, 1 ref, 2 motion]
    WrVTK_Type  = 3,          # write VTK files from adi to directory adi-vtk [1 surface, 2 lines, 3 both]
    VTKNacDim   = [-.10 ,-.10 ,-.10 ,.2 ,.2 ,.2],        # Nacelle Dimension for VTK visualization x0,y0,z0,Lx,Ly,Lz (m)
    VTKHubRad   = 1.0,                          # HubRadius for VTK visualization (m)
    adi_wrOuts  = 0,                            # output file format [0 none, 1 txt, 2 binary, 3 both]
    adi_DT_Outs = 0.0,                          # output frequency (resets to dt if smaller)
    adi_nstrut  = 2,                            # create_mesh_struts is hard coded for 2 struts per blade
    adi_dt      = 0.05,                         # random default
    adi_tmax    = 10,                           # end time
    omega       = 0,                            # rad/s
    hubPos      = [0,0,0],                      # m
    hubAngle    = [0,0,0]                       # rad
    )

    # Set up structs for the entire turbine
    adi_numbl = B + B*adi_nstrut    # Count struts as blades (strut each side of tower counted separately)
    Radius = maximum(bld_x)
    global turbine
    global turbenv
    numMeshNodes = getAD15numMeshNodes(bladeIdx)
    turbine = Turbine(Radius,omega,B,adi_numbl,numMeshNodes,bladeIdx,bladeElem,mymesh,myort)

    # Mesh info for ADI
    # set the origin for AD15 at the top of the "tower" (Ht in this setup)
    # nacelle -- not actually used here since we don't consider loads.
    nacPos    = Float32.([0,0,Ht])
    nacOrient = Float64.([1,0,0,0,1,0,0,0,1])
    # hub -- align to origin for now
    ADhubPos    = Float32.([0,0,Ht])
    ADhubOrient = Float64.([1,0,0,0,1,0,0,0,1])

    # set initial motion to 0
    u_j     = zeros(mymesh.numNodes*6)
    udot_j  = zeros(mymesh.numNodes*6)
    uddot_j = zeros(mymesh.numNodes*6)
    azi     = 0.0

    # blade roots (2nd is rotated 180 degrees about z)
    rootPos     = getRootPos(turbine,u_j,azi,hubPos,hubAngle)       # get root positions of all AD15 blades (blades + struts in OWENS)
    rootOrient  = getRootDCM(turbine,u_j,azi,hubAngle)              # get orientations of all AD15 blades   (blades + struts in OWENS)

    # Multiple mesh points along all blades for full structural mesh representation in ADI
    meshPos      = getAD15MeshPos(turbine,u_j,azi,hubPos,hubAngle)  # get positions of all AD15 nodes (blades + struts in OWENS)
    meshOrient   = getAD15MeshDCM(turbine,u_j,azi,hubAngle)         # get orientations of all AD15 blades   (blades + struts in OWENS)

    num_channels, channel_names, channel_units = adiInit(adi_lib,adi_rootname;
        ad_input_file_passed= false,
        ad_input_file=ad_input_file,
        ifw_input_file_passed= false,
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
        transposeDCM= transposeDCM,
        WrVTK       = WrVTK,
        WrVTK_Type  = WrVTK_Type,
        VTKNacDim   = VTKNacDim,
        VTKHubRad   = VTKHubRad,
        wrOuts      = adi_wrOuts,
        DT_Outs     = adi_DT_Outs,
        initHubPos         = ADhubPos,          # 3
        initHubOrient      = ADhubOrient,       # 9
        initNacellePos     = nacPos,            # 3 -- not actually used
        initNacelleOrient  = nacOrient,         # 9 -- not actually used
        numBlades          = adi_numbl,
        initRootPos        = rootPos,           # 3*adi_numbl
        initRootOrient     = rootOrient,        # 9*adi_numbl
        numMeshNodes       = numMeshNodes,
        initMeshPos        = meshPos,           # 3*numMeshNodes
        initMeshOrient     = meshOrient,        # 9*numMeshNodes
        interp_order=1,
        dt=adi_dt,
        t_max=adi_tmax
        )

    # can add some additional things here if needed
    adi_initialized = true

    # store channel info
    turbenv = Environment(rho,adi_dt,num_channels)

    # AD15 node velocities/accelerations
    #TODO: If OWENS sets these at initialization, need to transfer values here. 
    hubVel         = zeros(Float32,2*size(ADhubPos,1))
    hubAcc         = zeros(Float32,2*size(ADhubPos,1))
    nacVel         = zeros(Float32,2*size(nacPos,1))
    nacAcc         = zeros(Float32,2*size(nacPos,1))
    rootVel        = zeros(Float32,2*size(rootPos,1),size(rootPos,2))
    rootAcc        = zeros(Float32,2*size(rootPos,1),size(rootPos,2))
    meshVel        = zeros(Float32,2*size(meshPos,1),size(meshPos,2))
    meshAcc        = zeros(Float32,2*size(meshPos,1),size(meshPos,2))

    # Store turbine structure info for mesh points here
    global turbstruct=Structure(
            ADhubPos,  ADhubOrient,  hubVel,  hubAcc,
            nacPos,  nacOrient,  nacVel,  nacAcc,
            rootPos, rootOrient, rootVel, rootAcc,
            meshPos, meshOrient, meshVel, meshAcc,
        )

    # Initialize outputs and resulting mesh forces
    meshFrcMom = zeros(Cfloat,6*numMeshNodes);
    out_channel_vals = zeros(Cfloat,1,turbenv.num_channels)
    t_initial = 0.0

    # calculate initial values
    meshFrcMom, out_channel_vals = adiCalcOutput(t_initial,
            turbstruct.hubPos,  turbstruct.hubOrient,  turbstruct.hubVel,  turbstruct.hubAcc,
            turbstruct.nacPos,  turbstruct.nacOrient,  turbstruct.nacVel,  turbstruct.nacAcc,
            turbstruct.rootPos, turbstruct.rootOrient, turbstruct.rootVel, turbstruct.rootAcc,
            turbine.numMeshNodes,
            turbstruct.meshPos, turbstruct.meshOrient, turbstruct.meshVel, turbstruct.meshAcc,
            turbenv.num_channels);


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
* `hubAngle`:       3 angle set for hub orientation (rad)
* `hubVel`:         hub velocity in global coords, 6-vector (m/s,rad/s)
* `hubAcc`:         hub acceleration in global coords, 6-vector (m/s^2,rad/s^2)
"""
function deformAD15(u_j,udot_j,uddot_j,azi,Omega_rad,OmegaDot_rad,hubPos,hubAngle,hubVel,hubAcc)

    global turbine
    global turbenv
    global turbstruct

    # Hub       FIXME: add this later
    # Nacelle   FIXME: add this later

    # Root
    turbstruct.rootPos                    = getRootPos(turbine,u_j,azi,hubPos,hubAngle)                                                                         # get root positions of all AD15 blades (blades + struts in OWENS)
    turbstruct.rootOrient                 = getRootDCM(turbine,u_j,azi,hubAngle)                                                                                # get orientations of all AD15 blades   (blades + struts in OWENS)
    turbstruct.rootVel,turbstruct.rootAcc = getRootVelAcc(turbine,turbstruct.rootPos,udot_j,uddot_j,azi,Omega_rad,OmegaDot_rad,hubPos,hubAngle,hubVel,hubAcc)       # get root vel/acc of all AD15 blades   (blades + struts in OWENS)

    # Mesh
    turbstruct.meshPos                    = getAD15MeshPos(turbine,u_j,azi,hubPos,hubAngle)                                                                     # get mesh positions of all AD15 blades (blades + struts in OWENS)
    turbstruct.meshOrient                 = getAD15MeshDCM(turbine,u_j,azi,hubAngle)                                                                            # get orientations of all AD15 blades   (blades + struts in OWENS)
    turbstruct.meshVel,turbstruct.meshAcc = getAD15MeshVelAcc(turbine,turbstruct.meshPos,udot_j,uddot_j,azi,Omega_rad,OmegaDot_rad,hubPos,hubAngle,hubVel,hubAcc)   # get mesh vel/acc of all AD15 blades   (blades + struts in OWENS)

    # hub
#FIXME: this is not complete.  The hubVel is probably not correctly set.
    CG2H = calcHubRotMat(hubAngle, -azi)
    turbstruct.hubPos       = hubPos
    turbstruct.hubOrient    = vec(CG2H)
    turbstruct.hubVel       = hubVel
    turbstruct.hubAcc       = hubAcc

#TODO:  Transfer mesh structural deformation from OWENS, convert to global coords (if needed), and then apply to turbstruct
#            turbstruct.nacPos,  turbstruct.nacOrient,  turbstruct.nacVel,  turbstruct.nacAcc,
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
* `Fx`: Array(numMeshNodes,ntheta) Turbine Fx (N)
* `Fy`: Array(numMeshNodes,ntheta) Turbine Fy (N)
* `Fz`: Array(numMeshNodes,ntheta) Turbine Fz (N)
* `Mx`: Array(numMeshNodes,ntheta) Turbine Mx (N-m)
* `My`: Array(numMeshNodes,ntheta) Turbine My (N-m)
* `Mz`: Array(numMeshNodes,ntheta) Turbine Mz (N-m)

"""
function advanceAD15(t_new,mesh,azi;dt=turbenv.dt)
    global turbine
    global turbenv
    global turbstruct
    global dt
    dt = turbenv.dt

    n_steps=1   # hard code for now

    # Loads -- set the output array to the size of the OWENS mesh
    Fx = zeros(mesh.numNodes,n_steps)
    Fy = zeros(mesh.numNodes,n_steps)
    Fz = zeros(mesh.numNodes,n_steps)
    Mx = zeros(mesh.numNodes,n_steps)
    My = zeros(mesh.numNodes,n_steps)
    Mz = zeros(mesh.numNodes,n_steps)

    # conversion to hub coordinates (rotating)
    CG2H = calcHubRotMat(zeros(3), -azi)
    step1 = 0 #initialize scope
    for istep = 1:n_steps

        meshFrcMom = zeros(Cfloat,6*turbine.numMeshNodes);
        out_channel_vals = zeros(Cfloat,1,turbenv.num_channels)

        t = t_new - dt      # time ADI is starting from (used to copy states)
        # Call update states to go from T to T+dt       (NOTE: T and T+dt must be different)
        adiUpdateStates(t, t_new,
                turbstruct.hubPos,  turbstruct.hubOrient,  turbstruct.hubVel,  turbstruct.hubAcc,
                turbstruct.nacPos,  turbstruct.nacOrient,  turbstruct.nacVel,  turbstruct.nacAcc,
                turbstruct.rootPos, turbstruct.rootOrient, turbstruct.rootVel, turbstruct.rootAcc,
                turbine.numMeshNodes,
                turbstruct.meshPos, turbstruct.meshOrient, turbstruct.meshVel, turbstruct.meshAcc);

        meshFrcMom, out_channel_vals = adiCalcOutput(t_new,
                turbstruct.hubPos,  turbstruct.hubOrient,  turbstruct.hubVel,  turbstruct.hubAcc,
                turbstruct.nacPos,  turbstruct.nacOrient,  turbstruct.nacVel,  turbstruct.nacAcc,
                turbstruct.rootPos, turbstruct.rootOrient, turbstruct.rootVel, turbstruct.rootAcc,
                turbine.numMeshNodes,
                turbstruct.meshPos, turbstruct.meshOrient, turbstruct.meshVel, turbstruct.meshAcc,
                turbenv.num_channels);

        # copy mesh forces from meshFrcMom
        #   the AD15 mesh is a subset of points from the OWENS mesh. The bladeIdx array indicates the start and end nodes in the
        #   OWENS mesh for the current AD15 blade (struts are blades in AD15).
        #   This unpacking could be made much more compact by using only a single returned force array, but that can be done some
        #   other time. For now I want something that I can test -- ADP
        iBl = 0                                     # index before next blade in meshForceMom from AD15
        for ibld = 1:size(turbine.bladeIdx,1)       # blade or strut index
            idx1 = turbine.bladeIdx[ibld,1]         # start index of this blade in OWENS mesh
            idx2 = turbine.bladeIdx[ibld,2]         # end   index of this blade in OWENS mesh
            # AD15 blade may be upside down compared to OWENS mesh, so invert if needed
            sgn = sign(idx2-idx1)
            nNodes = abs(idx2-idx1)+1               # number of nodes on this blade/strut
            # step through all nodes on this blade/strut
            for iNode=0:nNodes-1                    # shift indices for simplicity
                Fx[idx1 + sgn*iNode,istep] = meshFrcMom[iBl + iNode*6 + 1]      # Fx
                Fy[idx1 + sgn*iNode,istep] = meshFrcMom[iBl + iNode*6 + 2]      # Fy
                Fz[idx1 + sgn*iNode,istep] = meshFrcMom[iBl + iNode*6 + 3]      # Fz
                Mx[idx1 + sgn*iNode,istep] = meshFrcMom[iBl + iNode*6 + 4]      # Mx
                My[idx1 + sgn*iNode,istep] = meshFrcMom[iBl + iNode*6 + 5]      # My
                Mz[idx1 + sgn*iNode,istep] = meshFrcMom[iBl + iNode*6 + 6]      # Mz
            end
            iBl = iBl + nNodes*6        # node index before next blade start from AD15
        end

        # convert to hub rotating coordinates with azimuth angle
        #   There is undoubtedly a much more elegant way to do this
        # TODO: need to check that this is actually correct
        for iNode=1:mesh.numNodes
            FMg = [Fx[iNode] Fy[iNode] Fz[iNode] Mx[iNode] My[iNode] Mz[iNode]]
            FM = frame_convert(FMg, CG2H)
            Fx[iNode,istep] = FM[1]
            Fy[iNode,istep] = FM[2]
            Fz[iNode,istep] = FM[3]
            Mx[iNode,istep] = FM[4]
            My[iNode,istep] = FM[5]
            Mz[iNode,istep] = FM[6]
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
getRootPos(turbine,u_j,azi,hubPos,hubAngle)

Extract the root positions for all ADI blades

# Inputs
* `turbine`:    turbine data storage
* `u_j`:        mesh displacements -- in hub coordinates, (m,rad)
* `azi`:        current azimuth (rad)
* `hubPos`:     current global hubPos (x,y,z) vector (m)
* `hubAngle`:   3 angle set for hub orientation (rad)
"""
function getRootPos(turbine,u_j,azi,hubPos,hubAngle)
    RootPos     = zeros(Float32,3,turbine.adi_numbl)
    # conversion from hub coordinates to global
    CG2H = calcHubRotMat(hubAngle, -azi)
    CH2G = transpose(CG2H)
    # blades
    for ibld = 1:turbine.B
        idx=turbine.bladeIdx[ibld,1]
        x=turbine.Mesh.x[idx] + u_j[(idx-1)*6+1]
        y=turbine.Mesh.y[idx] + u_j[(idx-1)*6+2]
        z=turbine.Mesh.z[idx] + u_j[(idx-1)*6+3]
        #println("Blade $ibld bottom at [$x,$y,$z] at index $idx")
        RootPos[:,ibld] = CH2G * [x; y; z] + hubPos
    end
    # bottom strut
    for ibld = turbine.B+1:2*turbine.B
        idx=turbine.bladeIdx[ibld,1]
        x=turbine.Mesh.x[idx] + u_j[(idx-1)*6+1]
        y=turbine.Mesh.y[idx] + u_j[(idx-1)*6+2]
        z=turbine.Mesh.z[idx] + u_j[(idx-1)*6+3]
        #println("Blade strut $ibld bottom at [$x,$y,$z] at index $idx")
        RootPos[:,ibld] = CH2G * [x; y; z] + hubPos
    end
    # top strut
    for ibld = 2*turbine.B+1:3*turbine.B
        idx=turbine.bladeIdx[ibld,1]
        x=turbine.Mesh.x[idx] + u_j[(idx-1)*6+1]
        y=turbine.Mesh.y[idx] + u_j[(idx-1)*6+2]
        z=turbine.Mesh.z[idx] + u_j[(idx-1)*6+3]
        #println("Blade strut $ibld top at [$x,$y,$z] at index $idx")
        RootPos[:,ibld] = CH2G * [x; y; z] + hubPos
    end
    return RootPos
end

"""
getRootVelAcc(turbine,udot_j,uddot_j,azi,Omega_rad,OmegaDot_rad)

Extract the root velocities and accelerations for all ADI blades

# Inputs
* `turbine`:        turbine data storage
* `rootPos`:        root positions from call to getRootPos
* `azi`:            current azimuth (rad)
* `Omega_rad`:      angular velocity of hub about hub axis (rad/s)
* `OmegaDot_rad`:   angular acceleration of hub about hub axis (rad/s^2)
* `hubPos`:         current global hubPos (x,y,z) vector (m)
* `hubAngle`:       3 angle set for hub orientation (rad)
* `hubVel`:         hub velocity in global coords, 6-vector (m/s,rad/s)
* `hubAcc`:         hub acceleration in global coords, 6-vector (m/s^2,rad/s^2)
"""
function getRootVelAcc(turbine,rootPos,udot_j,uddot_j,azi,Omega_rad,OmegaDot_rad,hubPos,hubAngle,hubVel,hubAcc)
    RootVel     = zeros(Float32,6,turbine.adi_numbl)
    RootAcc     = zeros(Float32,6,turbine.adi_numbl)
    ### 1. calculate relative velocity from mesh distortions
    # conversion from hub coordinates to global
    CG2H = calcHubRotMat(hubAngle, azi)
    CH2G = transpose(CG2H)
    # blades
    for ibld = 1:turbine.B
        tmp=turbine.bladeIdx[ibld,1]
        idx=(tmp-1)*6   # just before the node of interest
        RootVel[1:3,ibld] = CH2G *  udot_j[idx+1:idx+3]         # translation Vel (m/2)
        RootVel[4:6,ibld] = CH2G *  udot_j[idx+4:idx+6]         # rotation    Vel (rad/s)
        RootAcc[1:3,ibld] = CH2G * uddot_j[idx+1:idx+3]         # translation Acc (m/s^2)
        RootAcc[4:6,ibld] = CH2G * uddot_j[idx+4:idx+6]         # rotation    Acc (rad/s^2)
    end
    # bottom strut
    for ibld = turbine.B+1:2*turbine.B
        tmp=turbine.bladeIdx[ibld,1]
        idx=(tmp-1)*6   # just before the node of interest
        RootVel[1:3,ibld] = CH2G *  udot_j[idx+1:idx+3]         # translation Vel (m/2)
        RootVel[4:6,ibld] = CH2G *  udot_j[idx+4:idx+6]         # rotation    Vel (rad/s)
        RootAcc[1:3,ibld] = CH2G * uddot_j[idx+1:idx+3]         # translation Acc (m/s^2)
        RootAcc[4:6,ibld] = CH2G * uddot_j[idx+4:idx+6]         # rotation    Acc (rad/s^2)
    end
    # top strut
    for ibld = 2*turbine.B+1:3*turbine.B
        tmp=turbine.bladeIdx[ibld,1]
        idx=(tmp-1)*6   # just before the node of interest
        RootVel[1:3,ibld] = CH2G *  udot_j[idx+1:idx+3]         # translation Vel (m/2)
        RootVel[4:6,ibld] = CH2G *  udot_j[idx+4:idx+6]         # rotation    Vel (rad/s)
        RootAcc[1:3,ibld] = CH2G * uddot_j[idx+1:idx+3]         # translation Acc (m/s^2)
        RootAcc[4:6,ibld] = CH2G * uddot_j[idx+4:idx+6]         # rotation    Acc (rad/s^2)
    end
    ### 2. Tangential velocity due to hub rotation
    # calculate distance of point from hub axis, multiply by Omega_rad for tangential velocity component
#FIXME: is this correct, or is it CG2H[3,:]???????????
    # hub axis vector in global coordinates
    hubAxis = CG2H[:,3]
    for ibld = 1:size(RootVel,2)
        # Global coordinates
        # tangential velocity and acceleration, based on distance to hub axis 
        TanVel = cross(    Omega_rad*hubAxis, (rootPos[1:3,ibld]-hubPos)) / norm(hubAxis)
        TanAcc = cross( OmegaDot_rad*hubAxis, (rootPos[1:3,ibld]-hubPos)) / norm(hubAxis)
        RootVel[1:3,ibld] = RootVel[1:3,ibld] + TanVel
        RootVel[4:6,ibld] = RootVel[4:6,ibld] + Omega_rad*hubAxis
        RootVel[1:3,ibld] = RootVel[1:3,ibld] + TanAcc
        RootVel[4:6,ibld] = RootVel[4:6,ibld] + OmegaDot_rad*hubAxis
    end

    ### 3. add in contributions from hub motion in global coordinates
    for ibld = 1:size(RootVel,2)
        RootVel[1:3,ibld] = RootVel[1:3,ibld] + hubVel[1:3]
        RootVel[4:6,ibld] = RootVel[4:6,ibld] + hubVel[4:6]     # rad/s
        RootAcc[1:3,ibld] = RootAcc[1:3,ibld] + hubAcc[1:3]
        RootAcc[4:6,ibld] = RootAcc[4:6,ibld] + hubAcc[4:6]     # rad/s^2
    end
    return RootVel,RootAcc
end


"""
getAD15numMeshNodes(bladeIdx)

Find the number of mesh points we will pass
"""
function getAD15numMeshNodes(bladeIdx)
    # Find number of nodes -- note that we skip some OWENS mesh nodes along struts.  Also note that this method allows for a
    # different number of nodes in each blade (OWENS does not allow this at present)
    numMeshNodes = 0
    for i=1:size(bladeIdx,1)
        numMeshNodes = numMeshNodes + abs(bladeIdx[i,2] - bladeIdx[i,1]) + 1     # include all nodes in range 
    end
    return numMeshNodes
end

"""
getAD15MeshPos(turbine,u_j,azi,hubPos,hubAangle)

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
* `hubAngle`:   3 angle set for hub orientation (rad)
"""
function getAD15MeshPos(turbine,u_j,azi,hubPos,hubAngle)
    #display(turbine.bladeIdx)
    MeshPos     = zeros(Float32,3,turbine.numMeshNodes)
    iNode = 1
    # conversion from hub coordinates to global
    CG2H = calcHubRotMat(hubAngle, azi)
    CH2G = transpose(CG2H)
    # blades, bottom struts, top struts
    for ibld = 1:size(turbine.bladeIdx,1)
        npts = turbine.bladeIdx[ibld,2] - turbine.bladeIdx[ibld,1]
        sgn = 1*sign(npts)
        #println("                   iBld: $ibld     npts: $npts     sgn: $sgn")
        for idx=turbine.bladeIdx[ibld,1]:sgn:turbine.bladeIdx[ibld,2]
            x=turbine.Mesh.x[idx] + u_j[(idx-1)*6+1]
            y=turbine.Mesh.y[idx] + u_j[(idx-1)*6+2]
            z=turbine.Mesh.z[idx] + u_j[(idx-1)*6+3]
            #println("Blade $ibld at [$x,$y,$z], node $idx")
            MeshPos[:,iNode] = CH2G * [x; y; z] + hubPos
            iNode += 1
        end
    end
    return MeshPos
end

"""
getAD15MeshVelAcc(turbine,udot_j,uddot_j,azi,Omega_rad,OmegaDot_rad)

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
* `hubAngle`:       3 angle set for hub orientation (rad)
* `hubVel`:         hub velocity in global coords (m/s,rad/s)
* `hubAcc`:         hub acceleration in global coords (m/s^2,rad/s^2)
"""
function getAD15MeshVelAcc(turbine,meshPos,udot_j,uddot_j,azi,Omega_rad,OmegaDot_rad,hubPos,hubAngle,hubVel,hubAcc)
    #display(turbine.bladeIdx)
    MeshVel     = zeros(Float32,6,turbine.numMeshNodes)
    MeshAcc     = zeros(Float32,6,turbine.numMeshNodes)
    # conversion from hub coordinates to global
    CG2H = calcHubRotMat(hubAngle, azi)
    CH2G = transpose(CG2H)
#FIXME: is this correct, or is it CG2H[3,:]???????????
    # hub axis vector in global coordinates
    hubAxis = CG2H[:,3]
    # blades, bottom struts, top struts
    iNode = 1
    for ibld = 1:size(turbine.bladeIdx,1)
        npts = turbine.bladeIdx[ibld,2] - turbine.bladeIdx[ibld,1]
        sgn = 1*sign(npts)
        for tmpIdx=turbine.bladeIdx[ibld,1]:sgn:turbine.bladeIdx[ibld,2]
            idx=(tmpIdx-1)*6   # just before the node of interest
            ### 1. relative velocity from mesh distortions
            MeshVel[1:3,iNode] = CH2G *  udot_j[idx+1:idx+3]     # translation Vel (m/2)
            MeshVel[4:6,iNode] = CH2G *  udot_j[idx+4:idx+6]     # rotation    Vel (rad/s)
            MeshAcc[1:3,iNode] = CH2G * uddot_j[idx+1:idx+3]     # translation Acc (m/s^2)
            MeshAcc[4:6,iNode] = CH2G * uddot_j[idx+4:idx+6]     # rotation    Acc (rad/s^2)

#FIXME: missing tangential velocity components due to hub rotational velocity not about hub axis
            ### 2. Tangential velocity due to hub rotation
            # tangential velocity and acceleration, based on distance to hub axis 
            Vel = cross(    Omega_rad*hubAxis, (meshPos[1:3,iNode]-hubPos)) / norm(hubAxis)
            Acc = cross( OmegaDot_rad*hubAxis, (meshPos[1:3,iNode]-hubPos)) / norm(hubAxis)
            MeshVel[1:3,iNode] = MeshVel[1:3,iNode] + Vel
            MeshVel[4:6,iNode] = MeshVel[4:6,iNode] + Omega_rad*hubAxis
            MeshAcc[1:3,iNode] = MeshAcc[1:3,iNode] + Acc
            MeshAcc[4:6,iNode] = MeshAcc[4:6,iNode] + OmegaDot_rad*hubAxis

            ### 3. add in contributions from hub motion in global coordinates
            MeshVel[1:3,iNode] = MeshVel[1:3,iNode] + hubVel[1:3]
            MeshVel[4:6,iNode] = MeshVel[4:6,iNode] + hubVel[4:6]       # rad/s
            MeshAcc[1:3,iNode] = MeshAcc[1:3,iNode] + hubAcc[1:3]
            MeshAcc[4:6,iNode] = MeshAcc[4:6,iNode] + hubAcc[4:6]       # rad/s

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
* `hubAngle`:   3 angle set for hub orientation (rad)

#FIXME: add averaging of orientations to get nodes within blade/strut
"""
function getRootDCM(turbine,u_j,azi,hubAngle)
    RootOrient  = zeros(9,turbine.adi_numbl)
    # conversion from hub coordinates to global
    CG2H = calcHubRotMat(hubAngle, azi)
    CH2G = transpose(CG2H)
    for i=1:size(turbine.bladeElem,1)
        idx=turbine.bladeElem[i]    #,1]
        Psi     = turbine.Ort.Psi_d[idx]    - rad2deg(u_j[(idx-1)*6+4])
        Theta   = turbine.Ort.Theta_d[idx]  - rad2deg(u_j[(idx-1)*6+5])
        Twist   = turbine.Ort.Twist_d[idx]  - rad2deg(u_j[(idx-1)*6+6])
        if i<=turbine.B
            #FIME: the following is for a CCW spinning rotor.  some things need changing for a CW spinning rotor.
            # flip +z towards X, then apply Twist (Roll, Rx) -> Theta (Pitch, Ry) -> Psi (Yaw, Rz) 
            ang = [Psi, Theta, Twist, -90];
            DCM = CH2G * createGeneralTransformationMatrix(ang,[3,2,1,2]);
        else
            ang = [Psi, Theta, Twist,  90];
            DCM = CH2G * transpose(createGeneralTransformationMatrix(ang,[3,2,1,2]));
        end
        RootOrient[:,i] = transpose(vec(DCM))
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
* `hubAngle`:   3 angle set for hub orientation (rad)

#FIXME: add averaging of orientations to get nodes within blade/strut
"""
function getAD15MeshDCM(turbine,u_j,azi,hubAngle)
    #display(turbine.bladeElem)
    MeshOrient  = zeros(9,turbine.numMeshNodes)
    # conversion from hub coordinates to global
    CG2H = calcHubRotMat(hubAngle, azi)
    CH2G = transpose(CG2H)
    iNode=0
    for ibld=1:size(turbine.bladeElem,1)
        sgn = 1*sign(turbine.bladeElem[ibld,2] - turbine.bladeElem[ibld,1])
        if sgn > 0      # normal ordering
            #println("+ direction blade: $(turbine.bladeElem[ibld,1]):$(turbine.bladeElem[ibld,2])")
            for idx=(turbine.bladeElem[ibld,1]):(turbine.bladeElem[ibld,2])
                Psi     = turbine.Ort.Psi_d[idx]    - rad2deg(u_j[(idx-1)*6+4])
                Theta   = turbine.Ort.Theta_d[idx]  - rad2deg(u_j[(idx-1)*6+5])
                Twist   = turbine.Ort.Twist_d[idx]  - rad2deg(u_j[(idx-1)*6+6])
                #if idx==turbine.bladeElem[ibld,1]
                #    @printf("            normal   %i      %7.3f + %7.3f     %7.3f + %7.3f     %7.3f + %7.3f\n", ibld, turbine.Ort.Psi_d[idx], rad2deg(u_j[(idx-1)*6+4]), turbine.Ort.Theta_d[idx], rad2deg(u_j[(idx-1)*6+5]), turbine.Ort.Twist_d[idx], rad2deg(u_j[(idx-1)*6+6]))
                #end
                if ibld<=turbine.B
                    ang = [Psi, Theta, Twist, -90];
                    DCM = CH2G * createGeneralTransformationMatrix(ang,[3,2,1,2]);
                else
                    ang = [Psi, Theta, Twist,  90];
                    DCM = CH2G * transpose(createGeneralTransformationMatrix(ang,[3,2,1,2]));
                end
                iNode += 1
                MeshOrient[:,iNode] = transpose(vec(DCM))
                #println("+Blade $ibld Node $iNode orient: $ang            from Elem $idx")
                # duplicate last node
                if idx==turbine.bladeElem[ibld,2]
                    iNode += 1
                    MeshOrient[:,iNode] = transpose(vec(DCM))
                    #println("+Blade $ibld Node $iNode orient: $ang            from Elem $idx")
                end
            end
        else
            #println("- direction blade: $(turbine.bladeElem[ibld,1]):$(turbine.bladeElem[ibld,2])")
            for idx=(turbine.bladeElem[ibld,1]):-1:(turbine.bladeElem[ibld,2])
                Psi     = turbine.Ort.Psi_d[idx]    - rad2deg(u_j[(idx-1)*6+4])
                Theta   = turbine.Ort.Theta_d[idx]  - rad2deg(u_j[(idx-1)*6+5])
                Twist   = turbine.Ort.Twist_d[idx]  - rad2deg(u_j[(idx-1)*6+6])
                #if idx==turbine.bladeElem[ibld,1]
                #    @printf("            reverse  %i      %7.3f + %7.3f     %7.3f + %7.3f     %7.3f + %7.3f\n", ibld, turbine.Ort.Psi_d[idx], rad2deg(u_j[(idx-1)*6+4]), turbine.Ort.Theta_d[idx], rad2deg(u_j[(idx-1)*6+5]), turbine.Ort.Twist_d[idx], rad2deg(u_j[(idx-1)*6+6]))
                #end
                if ibld<=turbine.B
                    ang = [Psi, Theta, Twist, -90];
                    DCM = CH2G * createGeneralTransformationMatrix(ang,[3,2,1,2]);
                else
                    ang = [Psi, Theta, Twist,  90];
                    DCM = CH2G * transpose(createGeneralTransformationMatrix(ang,[3,2,1,2]));
                end
                iNode += 1
                MeshOrient[:,iNode] = transpose(vec(DCM))
                #println("-Blade $ibld Node $iNode orient: $ang            from Elem $idx")
                # duplicate last node
                if idx==turbine.bladeElem[ibld,2]
                    iNode += 1
                    MeshOrient[:,iNode] = transpose(vec(DCM))
                    #println("-Blade $ibld Node $iNode orient: $ang            from Elem $idx")
                end
            end
        end
    end
    return MeshOrient
end


### NOTE: the following are copied from ModelGen.jl/src/meshing_utilities.jl
#          (I couldn't figure out how to use ModelGen here without a mess)
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
        dcmTotal = dcmArray[:,:,i]*dcmTotal
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

function calcHubRotMat(ptfmRot, azi_j)
    CN2P = transMat(ptfmRot[1], ptfmRot[2], ptfmRot[3])
    CP2H = [cos(azi_j) sin(azi_j) 0; -sin(azi_j) cos(azi_j) 0;0 0 1]
    CN2H = CN2P*CP2H
    return CN2H
end
