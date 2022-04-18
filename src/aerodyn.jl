global adilib
global sym_calcoutput
global sym_updatestates
global sym_end
global backup_Vx

path,_ = splitdir(@__FILE__)

mutable struct ADI_Error
error_status
error_message
end

function ADI_Init(adilib_filename, output_root_name;
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
    storeHHVel  = false, 
    WrVTK       = 2,          # write VTK files from adi [0 none, 1 ref, 2 motion]
    transposeDCM= true,
    initHubPos         = zeros(3),
    initHubOrient      = zeros(9),
    initNacellePos     = zeros(3),
    initNacelleOrient  = zeros(9),
    numBlades          = 3,
    initRootPos        = zeros(numBlades,3),
    initRootOrient     = zeros(numBlades,9),
    #FIXME: need to define meshpoints
    numMeshPts         = 3,     # this must be set on call.... somehow
    initMeshPos        = zeros(numMeshPts,3),
    initMeshOrient     = zeros(numMeshPts,9),
    interp_order=1,
    t_initial=0.0,
    dt=0.01,
    t_max=60.0)

    global adi_abort_error_level = 4


    if ad_input_file == "none"
        ad_input_string_array = [
        ]
        error("Default AeroDyn input file not setup yet.")
    else
        println("Reading AeroDyn data from $ad_input_file.")
        fid = open(ad_input_file, "r") 
        ad_input_string_array = readlines(fid)
        close(fid)
    end

    ad_input_string        = join(ad_input_string_array, "\0")
    ad_input_string_length = length(ad_input_string)


    if ifw_input_file == "none"
        ifw_input_string_array = [
        ]
        error("Default InflowWind input file not setup yet.")
    else
        println("Reading InfloWind data from $ifw_input_file.")
        fid = open(ifw_input_file, "r") 
        ifw_input_string_array = readlines(fid)
        close(fid)
    end
    ifw_input_string        = join(ifw_input_string_array, "\0")
    ifw_input_string_length = length(ad_input_string)


    # Allocate Outputs
    num_channels = [0]
    channel_names = string(repeat(" ", 20 * 8000))
    channel_units = string(repeat(" ", 20 * 8000))

    global adilib = Libdl.dlopen(adilib_filename) # Open the library explicitly.
    global adi_active = true

    adi_sym_init = Libdl.dlsym(adilib, :AeroDyn_Inflow_C_Init)   # Get a symbol for the function to call.
    global adi_sym_calcoutput = Libdl.dlsym(adilib, :AeroDyn_Inflow_C_CalcOutput)   # Get a symbol for the function to call.
    global adi_sym_updatestates = Libdl.dlsym(adilib, :AeroDyn_Inflow_C_UpdateStates)
    global adi_sym_end = Libdl.dlsym(adilib, :AeroDyn_Inflow_C_End) # !!! "c" is capitalized in library, change if errors
    global adi_err = ADI_Error([0], string(repeat(" ", 1025)))
    
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
        Ref{Cdouble},       # IN: t_initial
        Ref{Cdouble},       # IN: dt
        Ref{Cdouble},       # IN: t_max
        Ref{Cint},          # IN: storeHHVel    (c_bool)
        Ref{Cint},          # IN: WrVTK
        Ref{Cint},          # IN: transposeDCM  (c_bool)
        Ref{Cfloat},        # IN: initHubPos
        Ref{Cdouble},       # IN: initHubOrient (do we need to flatten this, or just do fortran index order???)
        Ref{Cfloat},        # IN: initNacellePos
        Ref{Cdouble},       # IN: initNacelleOrient (do we need to flatten this, or just do fortran index order???)
        Ref{Cint},          # IN: numBlades
        Ref{Cfloat},        # IN: initRootPos
        Ref{Cdouble},       # IN: initRootOrient (do we need to flatten this, or just do fortran index order???)
        Ref{Cint},          # IN: numMeshPts
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
        interp_order,
        t_initial,
        dt,
        t_max,
        Cint.(storeHHVel),
        WrVTK,
        Cint.(transposeDCM),
        Cfloat.(initHubPos),
        Cdouble.(initHubOrient),
        Cfloat.(initNacellePos),
        Cdouble.(initNacelleOrient),
        numBlades,
        Cfloat.(initRootPos),
        Cdouble.(initRootOrient),
        numMeshPts,
        Cfloat.(initMeshPos),
        Cdouble.(initMeshOrient),
        num_channels,
        channel_names,
        channel_units,
        adi_err.error_status,
        adi_err.error_message)

    adi_check_error()

    return num_channels, channel_names, channel_units

end

function ADI_CalcOutput(time, 
                 hubPos, hubOrient, hubVel, hubAcc,
                 nacPos, nacOrient, nacVel, nacAcc,
                 rootPos, rootOrient, rootVel, rootAcc,
                 meshPos, meshOrient, meshVel, meshAcc, meshFrcMom,
                 out_channel_vals; num_node_pts=numMeshPts)
    # error_message = string(repeat(" ", 1025))
    # error_status = [0]

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
            Ref{Cint},          # IN: numMeshPts
            Ref{Cfloat},        # IN: meshPos 
            Ref{Cdouble},       # IN: meshOrient
            Ref{Cfloat},        # IN: meshVel 
            Ref{Cfloat},        # IN: meshAcc 
            Ptr{Cfloat},        # OUT: node_force
            Ptr{Cfloat},        # OUT: out_channel_vals
            Ptr{Cint},          # OUT: error_status
            Cstring),      # OUT: error_message 
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
            num_node_pts,
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
        error("AerooDyn-Inflow instance has not been initialized. Use ADI_Init() function.")
    end

    return meshFrcMom, out_channel_vals
end

function ADI_UpdateStates(time, next_time,
                 hubPos, hubOrient, hubVel, hubAcc,
                 nacPos, nacOrient, nacVel, nacAcc,
                 rootPos, rootOrient, rootVel, rootAcc,
                 meshPos, meshOrient, meshVel, meshAcc; num_node_pts=numMeshPts)
    # error_message = string(repeat(" ", 1025))
    # error_status = [0]

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
            Ref{Cint},          # IN: numMeshPts
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
            num_node_pts,
            Cfloat.(meshPos),
            Cdouble.(meshOrient),
            Cfloat.(meshVel),
            adi_err.error_status,
            adi_err.error_message) 

        adi_check_error()

    else
        error("AeroDyn-Inflow instance has not been initialized. Use ADI_Init() function.")
    end
end

function ADI_End()

    if adi_active

        global adi_active = false
        ccall(adi_sym_end,Cint,
        (Ptr{Cint},         # OUT: ErrStat_C
        Cstring),           # OUT: ErrMsg_C
        adi_err.error_status,
        adi_err.error_message)

        Libdl.dlclose(adilib) # Close the library explicitly.
        adi_check_error()

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
        @warn("Error status " * string(adi_err.error_status[1]) * ": " * string(adi_err.error_message))
        adi_err.error_status = [0] # reset error status/message
        adi_err.error_message = string(repeat(" ", 1025))
        ADI_End()
        error("AeroDyn-Inflow terminated prematurely.")
    end
end
