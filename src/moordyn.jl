global mdlib
global symCalcOutput
global symUpdateStates
global symEnd

path,_ = splitdir(@__FILE__)

mutable struct MDError
    error_status
    error_message
end

def_anchor_pts = [[418.8, 725.383, -200.],
                    [-837.6, 0.0, -200.],
                    [418.8, -725.383, -200.]]
def_fairlead_pts = [[20.434, 35.393, -14.],
                    [-40.868, 0., -14.],
                    [20.434, -35.393, -14.]]

function mdInit(;mdLib_filename="$path/../deps/bin/MoorDyn_c_lib_x64.dll", input_file="none", WtrDens=1025, WtrDpth=200, init_ptfm_pos=zeros(6), gravity=9.80665, dt=0.01, interp_order=1)

    global abort_error_level = 4
    
    if input_file == "none"
        # TODO specify pts and line type parameters when directly passing variables
        input_string_array = [
            "--------------------- MoorDyn v2.a8 Input File ------------------------------                                                                                                      ",
            "Mooring system for OC4-DeepCwind Semi                                                                                                                                              ",
            "---------------------- LINE TYPES -------------------------------------                                                                                                            ",
            "TypeName   Diam    Mass/m     EA     BA/-zeta   EI     Cd    Ca   CdAx  CaAx                                                                                                       ",
            "(-)        (m)     (kg/m)     (N)    (N-s/-)  (N-m^2)  (-)   (-)  (-)   (-)                                                                                                        ",
            "main     0.0766    113.35   7.536E8   -1.0      0      2.0   0.8  0.4   0.25                                                                                                       ",
            "------------------------ POINTS ---------------------------------------                                                                                                            ",
            "PointID Type     X        Y         Z     Mass  Volume  CdA    Ca                                                                                                                  ",
            "(-)     (-)     (m)      (m)       (m)    (kg)  (m?3)   (m^2)  (-)                                                                                                                 ",
            "1      Fixed   418.8    725.383  -200.0    0      0      0      0                                                                                                                  ",
            "2      Fixed  -837.6      0.0    -200.0    0      0      0      0                                                                                                                  ",
            "3      Fixed   418.8   -725.383  -200.0    0      0      0      0                                                                                                                  ",
            "4     Coupled   20.434   35.393   -14.0    0      0      0      0                                                                                                                  ",
            "5     Coupled  -40.868    0.0     -14.0    0      0      0      0                                                                                                                  ",
            "6     Coupled   20.434  -35.393   -14.0    0      0      0      0                                                                                                                  ",
            "------------------------ LINES ----------------------------------------                                                                                                            ",
            "LineID  LineType  UnstrLen   NumSegs  AttachA   AttachB  LineOutputs                                                                                                               ",
            "(-)       (-)       (m)        (-)    (point#)  (point#)    (-)                                                                                                                    ",
            "1         main     835.35      20        1         4         -                                                                                                                     ",
            "2         main     835.35      20        2         5         -                                                                                                                     ",
            "3         main     835.35      20        3         6         -                                                                                                                     ",
            "------------------------ OPTIONS ---------------------------------------                                                                                                           ",
            "0.001    dtM       - time step to use in mooring integration (s)                                                                                                                   ",
            "3.0e6    kbot      - bottom stiffness (Pa/m)                                                                                                                                       ",
            "3.0e5    cbot      - bottom damping (Pa-s/m)                                                                                                                                       ",
            "2.0      dtIC      - time interval for analyzing convergence during IC gen (s)                                                                                                     ",
            "150.0     TmaxIC    - max time for ic gen (s)                                                                                                                                       ",
            "4.0      CdScaleIC - factor by which to scale drag coefficients during dynamic relaxation (-)                                                                                      ",
            "0.01     threshIC  - threshold for IC convergence (-)                                                                                                                              ",
            "------------------------ OUTPUTS ---------------------------------------                                                                                                           ",
            "FairTen1                                                                                                                                                                           ",
            "FairTen2                                                                                                                                                                           ",
            "FairTen3                                                                                                                                                                           ",
            "AnchTen1                                                                                                                                                                           ",
            "AnchTen2                                                                                                                                                                           ",
            "AnchTen3                                                                                                                                                                           ",
            "END                                                                                                                                                                                ",
            "---------------------- need this line ----------------------------------                                                                                                           ",
            ""
            ]
    
    else
        println("Reading MoorDyn data from $input_file.")
        fid = open(input_file, "r") 
        input_string_array = readlines(fid)
        close(fid)
    end

    input_string        = join(input_string_array, "\0")
    input_string_length = length(input_string)

    # Allocate Outputs
    num_channels = [0]
    channel_names = string(repeat(" ", 20 * 4000))
    channel_units = string(repeat(" ", 20 * 4000))

    global mdLib = Libdl.dlopen(mdLib_filename) # Open the library explicitly.
    global mdActive = true
    global symInit = Libdl.dlsym(mdLib, :MD_INIT_C)   # Get a symbol for the function to call.
    global symCalcOutput = Libdl.dlsym(mdLib, :MD_CALCOUTPUT_C)
    global symUpdateStates = Libdl.dlsym(mdLib, :MD_UPDATESTATES_C)
    global symEnd = Libdl.dlsym(mdLib, :MD_END_C)
    global mdErr = MDError([0], string(repeat(" ", 1025)))

    ccall(symInit,Cint,
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
        Cstring),            # OUT: error_message
        [input_string],
        input_string_length,
        dt,
        gravity,
        WtrDens,
        WtrDpth,
        Cfloat.(init_ptfm_pos),
        interp_order,
        num_channels,
        channel_names,
        channel_units,
        mdErr.error_status,
        mdErr.error_message)

    mdCheckError() 

end

function mdCalcOutput(time, positions, velocities, accelerations, forces, out_channel_vals)

    if mdActive

        ccall(symCalcOutput,Cint,
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
            mdErr.error_status,
            mdErr.error_message) 

        mdCheckError()  

    else
        error("MoorDyn instance has not been initialized. Use mdInit() function.")
    end

    return forces, out_channel_vals

end

function mdUpdateStates(prev_time, curr_time, next_time, positions, velocities, accelerations)

    if mdActive

        ccall(symUpdateStates,Cint,
            (Ptr{Cdouble},      # IN: prev_time (t-dt)
            Ptr{Cdouble},       # IN: curr_time (t)
            Ptr{Cdouble},       # IN: next_time (t+dt)
            Ref{Cfloat},        # IN: positions
            Ref{Cfloat},        # IN: velocities
            Ref{Cfloat},        # IN: accelerations
            Ptr{Cint},          # OUT: error_status
            Cstring),           # OUT: error_message 
            [prev_time],
            [curr_time],
            [next_time],
            Cfloat.(positions),
            Cfloat.(velocities),
            Cfloat.(accelerations),
            mdErr.error_status,
            mdErr.error_message) 

        mdCheckError()

    else
        error("MoorDyn instance has not been initialized. Use mdInit() function.")
    end
end

function mdEnd()

    if mdActive
        global mdActive  = false
        ccall(symEnd,Cint,
        (Ptr{Cint},         # OUT: ErrStat_C
        Cstring),           # OUT: ErrMsg_C
        mdErr.error_status,
        mdErr.error_message)

        Libdl.dlclose(mdLib) # Close the library explicitly.

    end
end

function mdCheckError()
    if mdErr.error_status[1] == 0
        mdErr.error_status = [0] # reset error status/message
        mdErr.error_message = string(repeat(" ", 1025))
    elseif mdErr.error_status[1] < abort_error_level
        @warn("Error status " * string(mdErr.error_status[1]) * ": " * mdErr.error_message)
        mdErr.error_status = [0] # reset error status/message
        mdErr.error_message = string(repeat(" ", 1025))
    else
        @error("Error status " * string(mdErr.error_status[1]) * ": " * mdErr.error_message)
        mdErr.error_status = [0] # reset error status/message
        mdErr.error_message = string(repeat(" ", 1025))
        mdEnd()
        error("MoorDyn terminated prematurely.")
    end
end
