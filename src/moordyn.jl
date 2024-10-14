global hdlib
global sym_calcoutput
global sym_updatestates
global sym_end
global backup_Vx

mutable struct MD_Error
    error_status
    error_message
end

def_anchor_pts = [[418.8, 725.383, -200.],
                    [-837.6, 0.0, -200.],
                    [418.8, -725.383, -200.]]
def_fairlead_pts = [[20.434, 35.393, -14.],
                    [-40.868, 0., -14.],
                    [20.434, -35.393, -14.]]

function MD_Init(;mdlib_filename=nothing, md_input_file="none", WtrDens=1025, WtrDpth=200, init_ptfm_pos=zeros(6), gravity=9.80665, dt=0.01, interp_order=1)

    if isnothing(mdlib_filename)
        mdlib_filename="$path/../deps/openfast/build/modules/moordyn/libmoordyn_c_binding"
    end

    global md_abort_error_level = 4

    if md_input_file == "none"
        # TODO specify pts and line type parameters when directly passing variables
        input_string_array = [
            "--------------------- MoorDyn Input File ------------------------------------",
            "Mooring system for OC4-DeepCwind Semi",
            "FALSE    Echo      - echo the input file data (flag)",
            "----------------------- LINE TYPES ------------------------------------------",
            "Name     Diam      MassDen       EA    BA/-zeta    EI    Cd      Ca     CdAx   CaAx",
            "(-)       (m)      (kg/m)        (N)    (N-s/-)    (-)   (-)     (-)    (-)    (-)",
            "main     0.0766    113.35     7.536E8     -1.0      0    2.0     0.8    0.4    0.25",
            "---------------------- POINTS --------------------------------",
            "ID     Attachment  X        Y         Z      M      V       CdA   CA",
            "(-)    (-)        (m)      (m)       (m)    (kg)   (m^3)   (m^2)  (-)",
            "1      Fixed    418.8    725.383  -200.0     0       0       0     0",
            "2      Fixed   -837.6      0.0    -200.0     0       0       0     0",
            "3      Fixed    418.8   -725.383  -200.0     0       0       0     0",
            "4      Vessel    20.434   35.393   -14.0     0       0       0     0",
            "5      Vessel   -40.868    0.0     -14.0     0       0       0     0",
            "6      Vessel    20.434  -35.393   -14.0     0       0       0     0",
            "---------------------- LINES --------------------------------------",
            "ID      LineType   AttachA   AttachB  UnstrLen  NumSegs   Outputs",
            "(-)       (-)        (-)       (-)       (m)      (-)       (-)",
            "1         main        1         4       835.35     20        -",
            "2         main        2         5       835.35     20        -",
            "3         main        3         6       835.35     20        -",
            "---------------------- SOLVER OPTIONS ---------------------------------------",
            "0.001    dtM       - time step to use in mooring integration (s)",
            "3.0e6    kbot      - bottom stiffness (Pa/m)",
            "3.0e5    cbot      - bottom damping (Pa-s/m)",
            "2.0      dtIC      - time interval for analyzing convergence during IC gen (s)",
            "80.0     TmaxIC    - max time for ic gen (s)",
            "4.0      CdScaleIC - factor by which to scale drag coefficients during dynamic relaxation (-)",
            "0.01    threshIC  - threshold for IC convergence (-)",
            "------------------------ OUTPUTS --------------------------------------------",
            "FairTen1",
            "FairTen2",
            "FairTen3",
            "AnchTen1",
            "AnchTen2",
            "AnchTen3",
            "END",
            "------------------------- need this line --------------------------------------",
        ]

    else
        println("Reading MoorDyn data from $md_input_file.")
        fid = open(md_input_file, "r")
        input_string_array = readlines(fid)
        close(fid)
    end

    input_string        = join(input_string_array, "\0")
    input_string_length = length(input_string)

    # Allocate Outputs
    num_channels = [0]
    channel_names = string(repeat(" ", 20 * 4000))
    channel_units = string(repeat(" ", 20 * 4000))

    global mdlib = Libdl.dlopen(mdlib_filename) # Open the library explicitly.
    global md_active = true
    global md_sym_init = Libdl.dlsym(mdlib, :MD_C_Init)   # Get a symbol for the function to call.
    global md_sym_calcoutput = Libdl.dlsym(mdlib, :MD_C_CalcOutput)   # Get a symbol for the function to call.
    global md_sym_updatestates = Libdl.dlsym(mdlib, :MD_C_UpdateStates)
    global md_sym_end = Libdl.dlsym(mdlib, :MD_C_End) # !!! "c" is capitalized in library, change if errors
    global md_err = MD_Error([0], string(repeat(" ", 1025)))

    ccall(md_sym_init,Cint,
        (Ptr{Ptr{Cchar}},   # IN: input_string
        Ref{Cint},          # IN: input_string_length
        Ref{Cdouble},       # IN: dt
        Ref{Cfloat},        # IN: g
        Ref{Cfloat},        # IN: rho
        Ref{Cfloat},        # IN: depth
        Ref{Cfloat},        # IN: init_node_pos
        Ref{Cint},          # IN: interp_order
        Ptr{Cint},          # OUT: num_channels
        Cstring,            # OUT: channel_names
        Cstring,            # OUT: channel_units
        Ptr{Cint},          # OUT: error_status
        Cstring),           # OUT: error_message,
        [input_string],
        input_string_length,
        dt,
        Cfloat.(gravity),
        Cfloat.(WtrDens),
        Cfloat.(WtrDpth),
        Cfloat.(init_ptfm_pos),
        interp_order,
        num_channels,
        channel_names,
        channel_units,
        md_err.error_status,
        md_err.error_message)

    md_check_error()

end

function MD_CalcOutput(time, positions, velocities, accelerations, forces, out_channel_vals)

    if md_active

        ccall(md_sym_calcoutput,Cint,
            (Ptr{Cdouble},      # IN: time
            Ref{Cfloat},        # IN: positions
            Ref{Cfloat},        # IN: velocities
            Ref{Cfloat},        # IN: accelerations
            Ptr{Cfloat},        # OUT: forces
            Ptr{Cfloat},        # OUT: out_channel_vals
            Ptr{Cint},          # OUT: error_status
            Cstring),           # OUT: error_message
            [time],
            Cfloat.(positions),
            Cfloat.(velocities),
            Cfloat.(accelerations),
            forces,
            out_channel_vals,
            md_err.error_status,
            md_err.error_message)

        md_check_error()

    else
        error("MoorDyn instance has not been initialized. Use MD_Init() function.")
    end

    return forces, out_channel_vals

end

function MD_UpdateStates(curr_time, next_time, positions, velocities, accelerations)

    if md_active

        ccall(md_sym_updatestates,Cint,
            (Ptr{Cdouble},       # IN: curr_time (t)
            Ptr{Cdouble},       # IN: next_time (t+dt)
            Ref{Cfloat},        # IN: positions
            Ref{Cfloat},        # IN: velocities
            Ref{Cfloat},        # IN: accelerations
            Ptr{Cint},          # OUT: error_status
            Cstring),           # OUT: error_message
            [curr_time],
            [next_time],
            Cfloat.(positions),
            Cfloat.(velocities),
            Cfloat.(accelerations),
            md_err.error_status,
            md_err.error_message)

        md_check_error()

    else
        error("MoorDyn instance has not been initialized. Use MD_Init() function.")
    end
end

function MD_End()

    if md_active
        global md_active  = false
        ccall(md_sym_end,Cint,
        (Ptr{Cint},         # OUT: ErrStat_C
        Cstring),           # OUT: ErrMsg_C
        md_err.error_status,
        md_err.error_message)

        Libdl.dlclose(mdlib) # Close the library explicitly.

    end
end

function md_check_error()
    if md_err.error_status[1] == 0
        md_err.error_status = [0] # reset error status/message
        md_err.error_message = string(repeat(" ", 1025))
    elseif md_err.error_status[1] < md_abort_error_level
        @warn("Error status " * string(md_err.error_status[1]) * ": " * string(md_err.error_message))
        md_err.error_status = [0] # reset error status/message
        md_err.error_message = string(repeat(" ", 1025))
    else
        @warn("Error status " * string(md_err.error_status[1]) * ": " * string(md_err.error_message))
        md_err.error_status = [0] # reset error status/message
        md_err.error_message = string(repeat(" ", 1025))
        MD_End()
        error("MoorDyn terminated prematurely.")
    end
end
