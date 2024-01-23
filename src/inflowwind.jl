global inflowlib
global sym_calcoutput
global sym_end
global backup_Vx

path,_ = splitdir(@__FILE__)

ifw_active = true


"""
    ifwinit(inflowlib_filename ;HWindSpeed=6.87,turbsim_filename="path/test.bts")

calls inflow wind init

# Inputs
* `inflowlib_filename::string`: path and name of inflow-wind dynamic library
* `HWindSpeed::float`: optional, backup steady windspeed (m/s)
* `turbsim_filename::string`: path and name of turbsim data e.g. "path/test.bts"

# Outputs:
* `none`:

"""
function ifwinit(;inflowlib_filename="$path/../deps/bin/libifw_c_binding",HWindSpeed=10.125,turbsim_filename="$path/test.bts")

    # Where the input is manipulated
    HWindSpeed_str = "       $(round(HWindSpeed,digits=1))   HWindSpeed     - Horizontal windspeed                            (m/s)"
    turbsim_str = "\"$turbsim_filename\"      filename_bts   - name of the full field wind file to use (.bts)"
    uniformWind_str = "\"$turbsim_filename\"      FileName_Uni   - Filename of time series data for uniform wind field.      (-)"

    if turbsim_filename[end-3:end] == ".bts"
        WindType = 3
    else
        WindType = 2
    end

    #\x00------ InflowWind v3.01.* INPUT FILE -------------------------------------------------------------------------
    input_string_array = [
        "------ InflowWind v3.01.* INPUT FILE -------------------------------------------------------------------------",
        "Steady 15 m/s winds with no shear for IEA 15 MW Offshore Reference Turbine",
        "--------------------------------------------------------------------------------------------------------------",
            "false  Echo           - Echo input data to <RootName>.ech (flag)",
                "$WindType   WindType       - switch for wind file type (1=steady; 2=uniform; 3=binary TurbSim FF; 4=binary Bladed-style FF; 5=HAWC format; 6=User defined; 7=native Bladed FF)",
                "0   PropagationDir - Direction of wind propagation (meteoroligical rotation from aligned with X (positive rotates towards -Y) -- degrees)",
                "0   VFlowAng       - Upflow angle (degrees) (not used for native Bladed format WindType=7)",
            "False   VelInterpCubic - Use cubic interpolation for velocity in time (false=linear, true=cubic) [Used with WindType=2,3,4,5,7]",
                "1   NWindVel       - Number of points to output the wind velocity    (0 to 9)",
                "0   WindVxiList    - List of coordinates in the inertial X direction (m)",
                "0   WindVyiList    - List of coordinates in the inertial Y direction (m)",
                "1   WindVziList    - List of coordinates in the inertial Z direction (m)",
        "================== Parameters for Steady Wind Conditions [used only for WindType = 1] =========================",
                "$HWindSpeed_str",
                "150   RefHt          - Reference height for horizontal wind speed      (m)",
                "0.0   PLexp          - Power law exponent                              (-)",
        "================== Parameters for Uniform wind file   [used only for WindType = 2] ============================",
   "$uniformWind_str    FileName_Uni   - Filename of time series data for uniform wind field.      (-)",
                "15   RefHt_Uni      - Reference height for horizontal wind speed                (m)",
            "10   RefLength      - Reference length for linear horizontal and vertical sheer (-)",
        "================== Parameters for Binary TurbSim Full-Field files   [used only for WindType = 3] ==============",
        "$turbsim_str",
        "================== Parameters for Binary Bladed-style Full-Field files   [used only for WindType = 4] =========",
        "\"unused\"      FilenameRoot   - Rootname of the full-field wind file to use (.wnd, .sum)",
        "False         TowerFile      - Have tower file (.twr) (flag)",
        "================== Parameters for HAWC-format binary files  [Only used with WindType = 5] =====================",
        "\"unused\"      FileName_u     - name of the file containing the u-component fluctuating wind (.bin)",
        "\"unused\"      FileName_v     - name of the file containing the v-component fluctuating wind (.bin)",
        "\"unused\"      FileName_w     - name of the file containing the w-component fluctuating wind (.bin)",
                "64   nx             - number of grids in the x direction (in the 3 files above) (-)",
                "32   ny             - number of grids in the y direction (in the 3 files above) (-)",
                "32   nz             - number of grids in the z direction (in the 3 files above) (-)",
                "16   dx             - distance (in meters) between points in the x direction    (m)",
                "3   dy             - distance (in meters) between points in the y direction    (m)",
                "3   dz             - distance (in meters) between points in the z direction    (m)",
                "150   RefHt_HAWC     - reference height; the height (in meters) of the vertical center of the grid (m)",
        "-------------   Scaling parameters for turbulence   ---------------------------------------------------------",
                "2   ScaleMethod    - Turbulence scaling method   [0 = none, 1 = direct scaling, 2 = calculate scaling factor based on a desired standard deviation]",
                "1   SFx            - Turbulence scaling factor for the x direction (-)   [ScaleMethod=1]",
                "1   SFy            - Turbulence scaling factor for the y direction (-)   [ScaleMethod=1]",
                "1   SFz            - Turbulence scaling factor for the z direction (-)   [ScaleMethod=1]",
                "1.2   SigmaFx        - Turbulence standard deviation to calculate scaling from in x direction (m/s)    [ScaleMethod=2]",
                "0.8   SigmaFy        - Turbulence standard deviation to calculate scaling from in y direction (m/s)    [ScaleMethod=2]",
                "0.2   SigmaFz        - Turbulence standard deviation to calculate scaling from in z direction (m/s)    [ScaleMethod=2]",
        "-------------   Mean wind profile parameters (added to HAWC-format files)   ---------------------------------",
                "12   URef           - Mean u-component wind speed at the reference height (m/s)",
                "2   WindProfile    - Wind profile type (0=constant;1=logarithmic,2=power law)",
                "0.2   PLExp_HAWC     - Power law exponent (-) (used for PL wind profile type only)",
            "0.03   Z0             - Surface roughness length (m) (used for LG wind profile type only)",
                "0   XOffset        - Initial offset in +x direction (shift of wind box)",
        "================== LIDAR Parameters ===========================================================================",
                "0   SensorType          - Switch for lidar configuration (0 = None, 1 = Single Point Beam(s), 2 = Continuous, 3 = Pulsed)",
                "0   NumPulseGate        - Number of lidar measurement gates (used when SensorType = 3)",
                "30   PulseSpacing        - Distance between range gates (m) (used when SensorType = 3)",
                "0   NumBeam             - Number of lidar measurement beams (0-5)(used when SensorType = 1)",
            "-200   FocalDistanceX      - Focal distance co-ordinates of the lidar beam in the x direction (relative to hub height) (only first coordinate used for SensorType 2 and 3) (m)",
                "0   FocalDistanceY      - Focal distance co-ordinates of the lidar beam in the y direction (relative to hub height) (only first coordinate used for SensorType 2 and 3) (m)",
                "0   FocalDistanceZ      - Focal distance co-ordinates of the lidar beam in the z direction (relative to hub height) (only first coordinate used for SensorType 2 and 3) (m)",
        "0.0 0.0 0.0   RotorApexOffsetPos  - Offset of the lidar from hub height (m)",
                "17   URefLid             - Reference average wind speed for the lidar[m/s]",
            "0.25   MeasurementInterval - Time between each measurement [s]",
            "False   LidRadialVel        - TRUE => return radial component, FALSE => return 'x' direction estimate",
                "1   ConsiderHubMotion   - Flag whether to consider the hub motion's impact on Lidar measurements",
        "====================== OUTPUT ==================================================",
        "False         SumPrint       - Print summary data to <RootName>.IfW.sum (flag)",
                    "OutList        - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)",
        "\"Wind1VelX,Wind1VelY,Wind1VelZ\"     - Wind velocity at point WindVxiList(1),WindVyiList(1),WindVziList(1).  X, Y, and Z direction components.",
        "END of input file (the word \"END\" must appear in the first 3 columns of this last OutList line)",
    ]

    input_string        = join(input_string_array, "\0")
    input_string_length = length(input_string)

    # input_string_length = [length(input_string_array[1])]
    # Only needed for WindType = 2, can leave it empty if not used, but still need as input
    # ifw_uniform_string_array = [""] # if not used
    # could be an arbitrary number of lines long
    # NOTE: that there is an issue when passing in a binary coded string to fortran that an extra line/character is prepended, so we have to removed the first line to compensate
    # ! OpenFAST InflowWind uniform wind input file for 15 m/s wind.\x00

    if WindType == 2
        lines = open(turbsim_filename) do fid   #open mesh file
            readlines(fid)
        end
        
        linestart = 1
        for iline = 1:length(lines)
            if lines[iline][1]=='!' && iline !=1
                linestart = iline
                break
            end
        end
        
        uniform_string_array = lines[linestart:end]
        
    else

        uniform_string_array = [
            "! Time Wind  Wind  Vert. Horiz. Vert. LinV  Gust   Upflow",
            "!      Speed Dir   Speed Shear  Shear Shear Speed  Angle",
            "! (sec) (m/s) (deg) (m/s) (-)    (-)   (-)  (m/s)  (deg)",
            "0.0  15.0  0.0   0.0   5.0    0.0   0.0   0.0    0.0",
            "0.1  15.0  0.0   0.0   0.0    0.0   0.0   0.0    0.0",
            "1.0  15.0  0.0   0.0   0.0    0.0   0.0   0.0    0.0",
        ]
    end

    numWindPts = [1] #[Int(length(positions)/3)]     # total number of wind points requesting velocities for at each time step. must be integer

    uniform_string        = join(uniform_string_array, "\0")
    uniform_string_length = length(uniform_string)

    # Allocate Outputs
    numChannels = [0]
    channel_names = string(repeat(" ", 20 * 4000))
    channel_units = string(repeat(" ", 20 * 4000))
    error_status = [0]
    error_message = string(repeat(" ", 1025))


    # try
        println("Attempting to access inflow wind at: $inflowlib_filename")
        global inflowlib = Libdl.dlopen(inflowlib_filename) # Open the library explicitly.
        sym_init = Libdl.dlsym(inflowlib, :IfW_C_Init)   # Get a symbol for the function to call.
        global sym_calcoutput = Libdl.dlsym(inflowlib, :IfW_C_CalcOutput)   # Get a symbol for the function to call.
        global sym_end = Libdl.dlsym(inflowlib, :IfW_C_End)

        ccall(sym_init,Cint,
            (Ptr{Ptr{Cchar}}, # IN: input_string_array
            Ref{Cint},     # IN: input file string length
            Ptr{Ptr{Cchar}}, # IN: uniform_string_array
            Ref{Cint},     # IN: input file string length
            Ptr{Cint}, # IN: numWindPts
            Ptr{Cdouble}, # IN: dt
            Ptr{Cint}, # OUT: numChannels
            Cstring, # OUT: channel_names
            Cstring, # OUT: channel_units
            Ptr{Cint}, # OUT: error_status
            Cstring), # OUT: error_message
            [input_string], #InputFileStrings_C
            input_string_length,
            [uniform_string], #InputUniformStrings_C
            uniform_string_length,
            numWindPts, #NumWindPts_C
            [0.0], #DT_C !!!NOT USED!!!
            numChannels,#NumChannels_C
            channel_names, #OutputChannelNames_C
            channel_units, #OutputChannelUnits_C
            error_status, #ErrStat_C
            error_message) #ErrMsg_C
            if error_status[1] != 0
                @warn "Inflow Wind Init Error"
                println(error_message)
            end
        # catch
        #     global backup_Vx = HWindSpeed
        #     @warn "inflow wind library is not configured properly, using a steady wind profile [$backup_Vx,0.0,0.0] for Vx, Vy, Vz"
        #     global ifw_active = false
        # end

end


"""
    ifwcalcoutput(position,time)

calls inflow wind clacoutput

# Inputs
* `position::Array(float)`: x, y, z sample position (m)
* `time::float`: sample time (s)

# Outputs:
* `velocities`: x, y, z velocity at sample position

"""
function ifwcalcoutput(position,time)

    outputChannelValues = zeros(3)

    #Reset error message
    error_message = string(repeat(" ", 1025))
    error_status = [0]
    velocities = zeros(Float32,3) # output velocities (N x 3)

    if ifw_active

        ccall( sym_calcoutput,Cint,
        (Ptr{Cdouble}, #Time_C
        Ptr{Cfloat}, #Positions_C
        Ptr{Cfloat}, #Velocities_C
        Ptr{Cfloat}, #OutputChannelValues_C
        Ptr{Cint}, #ErrStat_C
        Cstring), #ErrMsg_C
        [time],                 # IN: time at which to calculate velocities
        Float32.(position),                      # IN: positions - specified by user
        velocities,                     # OUT: velocities at desired positions
        outputChannelValues,                 # OUT: output channel values as described in input file
        error_status,              # ErrStat_C
        error_message)                     # ErrMsg_C

        if error_status[1] != 0
            @warn "Inflow Wind Calcoutput Error"
            println(error_message)
        end

        return velocities
    else
        # @warn "using 10m/s inflow in the x direction since inflow wind library is misconfigured"
        return [backup_Vx,0.0,0.0]
    end
end

"""
    ifwend()

calls inflow wind end function and cleanup
"""
function ifwend()
    #Reset error message
    error_message = string(repeat(" ", 1025))
    error_status = [0]

    if ifw_active
        ccall( sym_end,Cint,
        (Ptr{Cint}, #ErrStat_C
        Cstring), #ErrMsg_C
        error_status,
        error_message)

        if error_status[1] != 0
            @warn "Inflow Wind End Error"
            println(error_message)
        end

        Libdl.dlclose(inflowlib) # Close the library explicitly.

    end
end
