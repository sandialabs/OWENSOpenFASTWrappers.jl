# Example usage the Julia wrapper for the AeroDyn DLL

import FLOWMath
import GyricFEA
import OWENS
import VAWTAero
import OWENS

import DelimitedFiles
path = splitdir(@__FILE__)[1]
# import OpenFASTWrappers
using Test

import PyPlot
PyPlot.pygui(true)
PyPlot.rc("figure", figsize=(15, 15))
PyPlot.rc("font", size=10.0)
PyPlot.rc("lines", linewidth=1.5)
PyPlot.rc("lines", markersize=4.0)
PyPlot.rc("legend", frameon=true)
PyPlot.rc("axes.spines", right=false, top=false)
PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
PyPlot.rc("figure",max_open_warning=500)
# PyPlot.rc("axes", prop_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

# println("Calculate/Set up sectional properties")
# include("$path/setup_XFlow_25kW.jl")
# XFlow specific mesh generator
# include("$path/setupOWENShawt.jl")
include("$path/../../src/OpenFASTWrappers.jl")

# # Run the standalone aerodyn
# run(`$path/../../../openfastandy/build/modules/aerodyn/aerodyn_driver $path/HAWT_standalone_test.dvr`)

adi_lib = "$path/../../../openfastandy/build/modules/aerodyn/libaerodyn_inflow_c_binding" #change this to match your local path of the AeroDyn DLL
# adi_lib = "/builds/8921-VAWT-TOOLS/OpenFASTWrappers.jl/openfast/build/modules/AeroDyn/libAeroDyn_c_binding" #change this to match your local path of the AeroDyn DLL

# output files from ADI
adi_rootname = "$path/HAWT_OWENS_AD15"
adi_rootname_twoway = "$path/HAWT_OWENS_AD15_stiff_one_way"

num_corrections = 0


dt = 0.01
t_initial = 0.0
t_max = 0.1
ts = collect(t_initial:dt:t_max)
numTS = length(ts)

saveName = "$path/vtk/OWENS_HAWT_test"
rho = 1.225
mu = 1.4639E-05 
Nslices = 30
RPM = 6.0
Ht = 137.0
B = Nbld = 3
H = 10.54#-1.0   # m 
R = 61.5     # m
Area = H*R*2
Vinf = 10.0 #m/s
chord = 0.5*ones(Nslices)
omega = RPM / 60 * 2 * pi
blade_joint_angle_Degrees = 0.0

# These aren't used by aerodyn...
shapeX = LinRange(0,R,Nslices+1)
shapeY = zero(shapeX)

mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
    stiff_twr, stiff_bld,RefArea,bld_precompinput,
    bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
    twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,RefArea,
    mass_breakout_blds,mass_breakout_twr,bladeIdx,bladeElem,system,assembly,sections = OWENS.setupOWENShawt(VAWTAero,path;
    rho,
    Nslices,
    RPM,
    Vinf,
    eta = 0.5,
    B,
    H,
    R,
    hubR = 0.0,
    shapeY,
    shapeX,
    NuMad_geom_xlscsv_file_twr = nothing,
    NuMad_mat_xlscsv_file_twr = nothing,
    NuMad_geom_xlscsv_file_bld = nothing,
    NuMad_mat_xlscsv_file_bld = nothing,
    stack_layers_bld = nothing,
    Ht=0.0,
    ntelem = 10, #tower elements
    nbelem = 19, #blade elements
    ncelem = 10,
    joint_type = 0,
    RPI=true,
    biwing=false,
    hub_depth = 1.0, #Hub Beam Depth
    AD15_ccw = true,
    angularOffset = -pi/2,
    # R_root = 10.0, # m biwing radius
    # R_biwing = 30.0, # outer radius
    # R_tip = 54.014, # outer radius
    # nbelem_root = 30, #biwing elements, for each 
    # nbelem_biwing = 30, #tip elements
    # nbelem_tip = 30, #tip elements
    # bshapex_root = LinRange(0.0,R_root,nbelem_root+1), #Blade shape, magnitude is relevant
    # bshapez_root = zeros(nbelem_root+1), #Blade shape, magnitude is relevant
    # bshapex_biwing_U = LinRange(R_root,R_biwing,nbelem_biwing+1), #Blade shape, magnitude is relevant
    # bshapez_biwing_U = zeros(nbelem_biwing+1), #Blade shape, magnitude is relevant
    # bshapex_biwing_L = LinRange(R_root,R_biwing,nbelem_biwing+1), #Blade shape, magnitude is relevant
    # bshapez_biwing_L = zeros(nbelem_biwing+1), #Blade shape, magnitude is relevant
    # bshapex_tip = LinRange(R_biwing,R_tip,nbelem_tip+1), #Blade shape, magnitude is relevant
    # bshapez_tip = zeros(nbelem_tip+1), #Blade shape, magnitude is relevant
    )




    # PyPlot.figure()
    # PyPlot.scatter3D(mymesh.x,mymesh.y,mymesh.z,color=plot_cycle[1])
    # # PyPlot.zlim(([90,120]))
    
    # # Allocate the node angle arrays
    # Psi_d_Nodes = zeros(mymesh.numNodes)
    # Theta_d_Nodes = zeros(mymesh.numNodes)
    # Twist_d_Nodes = zeros(mymesh.numNodes)
    
    # function rotate_normal(i_el, mymesh, myort;vec=[0,1,0.0],normal_len=1)
    #     # * `Psi_d::Array{<:float}`: length NumEl, element rotation about 3 in global FOR (deg) These angles are used to transform from the global coordinate frame to the local element/joint frame via a 3-2 Euler rotation sequence.
    #     # * `Theta_d::Array{<:float}`: length NumEl, element rotation about 2 (deg)
    #     # * `Twist_d::Array{<:float}`: length NumEl, element twist (deg)
        
    #     # Map from element to node
    #     nodenum1 = Int(mymesh.conn[i_el,1])
    #     nodenum2 = Int(mymesh.conn[i_el,2])
    
    #     # Extract the element angles
    #     Psi_d_el = myort.Psi_d[i_el]
    #     Theta_d_el = myort.Theta_d[i_el]
    #     Twist_d_el = myort.Twist_d[i_el]
    
    #     # Map the element angles to the node angles
    #     Psi_d_Nodes[[nodenum1,nodenum2]] .= Psi_d_el
    #     Theta_d_Nodes[[nodenum1,nodenum2]] .= Theta_d_el
    #     Twist_d_Nodes[[nodenum1,nodenum2]] .= Twist_d_el
    
    #     # Use a line and rotate it about the angles, different starting vectors show different angles.
    #     myvec = vec'.*normal_len
    
    #     # apply the twist rotation, which is about the x (1) axis
    #     myvec = myvec*[1.0 0.0 0.0
    #         0.0 cosd(Twist_d_el) sind(Twist_d_el)
    #         0.0 -sind(Twist_d_el) cosd(Twist_d_el)]
    
    #     # apply theta rotation, which is the tilt angle, or about the y (2) axis in global
    #     myvec = myvec*[cosd(Theta_d_el) 0.0 -sind(Theta_d_el)
    #         0.0 1.0 0.0
    #         sind(Theta_d_el) 0.0 cosd(Theta_d_el)]
    
    #     # apply Psi rotation, which is about Z (3) axis in global
    #     myvec = myvec*[cosd(Psi_d_el) sind(Psi_d_el) 0.0
    #         -sind(Psi_d_el) cosd(Psi_d_el) 0.0
    #         0.0 0.0 1.0]
    
    #     # Get the location of the element
    #     x_el = (mymesh.x[nodenum1]+mymesh.x[nodenum2])/2
    #     y_el = (mymesh.y[nodenum1]+mymesh.y[nodenum2])/2
    #     z_el = (mymesh.z[nodenum1]+mymesh.z[nodenum2])/2
    
    #     # Offset the myvector by the location of the element
    #     myvec = myvec + [x_el,y_el,z_el]'
    #     x_el_plot = [x_el,myvec[1]]
    #     y_el_plot = [y_el,myvec[2]]
    #     z_el_plot = [z_el,myvec[3]]
    #     return x_el_plot, y_el_plot, z_el_plot
    # end
    
    # # Add the orientation vectors, ort is on elements
    # for i_el = 1:mymesh.numEl
    #     x_el_plot, y_el_plot, z_el_plot = rotate_normal(i_el, mymesh, myort;vec=[5,0,0.0])
    #     PyPlot.plot3D(x_el_plot,y_el_plot,z_el_plot,"-",color=plot_cycle[2])
    # end
    
    # for i_el = 1:mymesh.numEl
    #     x_el_plot, y_el_plot, z_el_plot = rotate_normal(i_el, mymesh, myort;vec=[0,5,0.0])
    #     PyPlot.plot3D(x_el_plot,y_el_plot,z_el_plot,"-",color=plot_cycle[3])
    # end
    
    # for i_el = 1:mymesh.numEl
    #     x_el_plot, y_el_plot, z_el_plot = rotate_normal(i_el, mymesh, myort;vec=[0,0,5.0])
    #     PyPlot.plot3D(x_el_plot,y_el_plot,z_el_plot,"-",color=plot_cycle[4])
    # end
    
    # PyPlot.plot3D([0.0],[0.0],[0.0],"-",color=plot_cycle[2],label="X-norm")
    # PyPlot.plot3D([0.0],[0.0],[0.0],"-",color=plot_cycle[3],label="Y-norm")
    # PyPlot.plot3D([0.0],[0.0],[0.0],"-",color=plot_cycle[4],label="Z-norm")
    # PyPlot.legend()
    # # for i_joint = 1:length(myjoint[:,1])
    # #     i_el = findall(x->x==myjoint[i_joint,2],mymesh.conn[:,1])
    # #     if length(i_el)==0 #Use the other element associated with the joint
    # #         i_el = findall(x->x==myjoint[i_joint,2],mymesh.conn[:,2])
    # #     end
    # #     if length(i_el)==0 #Use the other element associated with the joint
    # #         i_el = findall(x->x==myjoint[i_joint,3],mymesh.conn[:,2])
    # #     end
    # #     x_el_plot, y_el_plot, z_el_plot = rotate_normal(i_el[1], mymesh, myort;normal_len=3)
    # #     PyPlot.plot3D(x_el_plot,y_el_plot,z_el_plot,"-",color=plot_cycle[5])
    # # end
    
    # PyPlot.xlabel("x")
    # PyPlot.ylabel("y")
    # PyPlot.zlabel("z")
    # PyPlot.axis("equal")

# # Create AeroDyn Files

# # input files for ADI
ad_input_file="$path/ad_primary.dat"
ifw_input_file="$path/ifw_primary.dat"
# blade_filename="$path/AD15-input/blade2.dat"
# lower_strut_filename="$path/AD15-input/lower_longer_arm2.dat"
# upper_strut_filename="$path/AD15-input/upper_shorter_arm2.dat"
# airfoil_filenames = "$path/AD15-input/Airfoils/NACA_0018_AllRe.dat"
# OLAF_filename = "$path/AD15-input/OLAF2.dat"

# blade_filenames = [blade_filename,blade_filename,blade_filename,lower_strut_filename,lower_strut_filename,lower_strut_filename,upper_strut_filename,upper_strut_filename,upper_strut_filename]


# OpenFASTWrappers.writeADinputFile(ad_input_file,blade_filenames,airfoil_filenames,OLAF_filename)

# NumBlNds=10
# bld_len = [10.54,10.54,10.54,4.790113480486895,4.790113480486895,4.790113480486895,4.686747309213206,4.686747309213206,4.686747309213206]
# for (ifile,filename) in enumerate(blade_filenames)
#     BlSpn=collect(LinRange(0,bld_len[ifile],NumBlNds))
#     BlCrvAC=zeros(NumBlNds)
#     BlSwpAC=zeros(NumBlNds)
#     BlCrvAng=zeros(NumBlNds)
#     BlTwist=zeros(NumBlNds)
#     BlChord=ones(NumBlNds).*0.5
#     BlAFID=ones(Int,NumBlNds)
#     OpenFASTWrappers.writeADbladeFile(filename;NumBlNds,BlSpn,BlCrvAC,BlSwpAC,BlCrvAng,BlTwist,BlChord,BlAFID)
# end

# OpenFASTWrappers.writeOLAFfile(OLAF_filename;nNWPanel=200,nFWPanels=10)

# OpenFASTWrappers.writeIWfile(10.0,ifw_input_file;turbsim_filename=nothing)

# Run AeroDyn
OpenFASTWrappers.setupTurb(adi_lib,ad_input_file,ifw_input_file,adi_rootname,[shapeX],[shapeY],[B],[Ht],[mymesh],[myort],[bladeIdx],[bladeElem];
        rho     = rho,
        adi_dt  = dt,
        adi_nstrut  = 0,
        adi_tmax= t_max,
        omega   = [omega],
        adi_wrOuts = 1,     # write output file [0 none, 1 txt, 2 binary, 3 both]
        adi_DT_Outs = dt,    # output frequency
        hubAngle = [[0.0,0.0,0.0]], #deg
        hubPos      = [[0,0,Ht]],
        nacAngle = [[0.0,0.0,0.0]], #deg
        nacPos      = [[0,0,0.0]],  
        refPos = [[0,0.0,0.0]],
        isHAWT = true
        )

# # Time marching
# ForceValHist = zeros(Int(mymesh.numNodes*6),length(ts[1:end-1]))
# AziHist = zeros(length(ts[1:end-1]))
# Fzhist = zeros(mymesh.numNodes,length(ts[1:end-1]))

# for (tidx, t) in enumerate(ts[1:end-1])
#     # println("time $t")
#     for correction = 1:num_corrections+1
#         # println("correction $correction")

#         #mapAD15(t,azi_j,mesh,advanceAD15)
#         azi_j = omega*(t+dt) # rad/s * s
#         AziHist[tidx] = azi_j
#         u_j     = zeros(mymesh.numNodes*6)
#         udot_j  = zeros(mymesh.numNodes*6)
#         uddot_j = zeros(mymesh.numNodes*6)
#         hubPos      = [0,0,Ht]                      # m
#         hubAngle    = [0,90,0.0]                       # rad
#         hubVel = zeros(6)
#         hubAcc = zeros(6)
#         OpenFASTWrappers.deformAD15([u_j],[udot_j],[uddot_j],[azi_j],[omega],[zero(omega)],[hubPos],[hubAngle],[hubVel],[hubAcc])
#         n_steps,Fx,Fy,Fz,Mx,My,Mz = OpenFASTWrappers.advanceAD15(t,[mymesh],azi_j)
#         Fzhist[:,tidx] = Fz[1][:,1]
#         # NOTE on AD15 advanceTurb values (Fx,Fy,Fz,Mx,My,Mz)
#         #       - forces/moments are in hub coordinates (converted in advanceAD15)
#         #       - array length is the number of OWENS mesh points
#         #       - This includes the struts (and could include tower when we add that to the AD15 interface)
    
#         #     [~,~,timeLen] = size(aeroDistLoadsArrayTime)
        
    
#         # Map loads over from advanceTurb
#         for i=1:mymesh.numNodes
#             ForceValHist[(i-1)*6+1,tidx] = Fz[1][i,1]
#             ForceValHist[(i-1)*6+2,tidx] = Fy[1][i,1]
#             ForceValHist[(i-1)*6+3,tidx] = -Fx[1][i,1]
#             ForceValHist[(i-1)*6+4,tidx] = Mz[1][i,1]
#             ForceValHist[(i-1)*6+5,tidx] = My[1][i,1]
#             ForceValHist[(i-1)*6+6,tidx] = -Mx[1][i,1]
#         end
#         # DOFs are sequential through all nodes
#         ForceDof=collect(1:1:mymesh.numNodes*6)
#     end

# end
# println("End Aerodynamics")
# OpenFASTWrappers.endTurb()


#########################################
### Run STIFF one way coupling
#########################################

dt = 0.01
t_initial = t_initial = 0.0
t_max = 0.5
ts = collect(t_initial:dt:t_max)
numTS = length(ts)

mymesh.hubPos =  [0,0.0,137.0]#TODO: make this streamlined in the code initialization, consider putting it elsewhere.
mymesh.hubAngle = [0.0,-90,0.0]

#########################################
### Create Aero Functions
#########################################

# #FIXME: this is a placeholder call returning no loads
# OpenFASTWrappers.setupTurb(adi_lib,ad_input_file,ifw_input_file,adi_rootname_twoway,[shapeX],[shapeY],[B],[Ht],[mymesh],[myort],[bladeIdx],[bladeElem];
#         rho     = rho,
#         adi_dt  = dt,
#         adi_nstrut  = 0,
#         adi_tmax= t_max,
#         omega   = [omega],
#         adi_wrOuts = 1,     # write output file [0 none, 1 txt, 2 binary, 3 both]
#         adi_DT_Outs = dt,    # output frequency
#         hubAngle = [[0.0,0.0,0.0]], #deg
#         refPos = [[0,0.0,137.0]],
#         isHAWT = true
#         )

# Routine for getting aero forces from aD15

aeroForcesAD15(t,azi) = OWENS.mapAD15(t,azi,[mymesh],OpenFASTWrappers.advanceAD15;alwaysrecalc=true,verbosity=1)


######################################
#### Perform Stiff One Way Test
#######################################

pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0]

model = OWENS.Inputs(;analysisType = "GX",
    tocp = [0.0,100000.1],
    Omegaocp = [RPM,RPM] ./ 60,
    turbineStartup = 0,
    generatorOn = false,
    useGeneratorFunction = false,
    numTS = numTS,
    delta_t = dt,
    aeroLoadsOn = 2.0,        # 1: one way; 1.5: aero once in iteration in 2 way; 2: two-way
    AD15On = true,
    topsideOn = true,
    interpOrder = 2)
model.iteration_parameters.MAXITER=20   # temporary for testing

feamodel = GyricFEA.FEAModel(;analysisType = "GX",
    outFilename = "none",
    joint = myjoint,
    platformTurbineConnectionNodeNumber = 1,
    pBC,
    nlOn = false,
    gravityOn = false, # turn off since aero only doesn't include gravity, also turned off rotational effects in setupxxx.jl
    numNodes = mymesh.numNodes,
    RayleighAlpha = 0.05,
    RayleighBeta = 0.05,
    return_all_reaction_forces = true,
    iterationType = "DI")

println("Running Unsteady")
t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,
FReactionHist_stiff_one_way,FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,
uHist_prp,epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,
kappa_y_hist,kappa_z_hist,FPtfmHist,FHydroHist,FMooringHist,
topFexternal_hist,rbDataHist = OWENS.Unsteady_Land(model;topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForcesAD15,
    deformAero=OpenFASTWrappers.deformAD15,system=system,assembly=assembly) #,meshcontrolfunction=mymeshcontrolfunction2,userDefinedGenerator=userDefinedGenerator,



println("End Aerodynamics")
OpenFASTWrappers.endTurb()
####################################################################################################################################################################################################################################


# println("Saving VTK time domain files")
# OWENS.gyricFEA_VTK("$path/vtk/OWENS_test_HAWT1",t,uHist,system,assembly,sections;scaling=1)#,azi=aziHist)

# ######################################
# #### Plot Comparison
# #######################################

# # Load Standalone Data
# standalone_AD = DelimitedFiles.readdlm("$path/HAWT_standalone_test.out",skipstart=8)
# header_standalone = DelimitedFiles.readdlm("$path/HAWT_standalone_test.out",header=true,skipstart=6)[2]

# # Compare with Current Data on the AeroDyn Side
# library_ADside = DelimitedFiles.readdlm("$adi_rootname.out",skipstart=8)
# header_library = DelimitedFiles.readdlm("$adi_rootname.out",header=true,skipstart=6)[2]

# # library_ADside_owens_stiff_one_way = DelimitedFiles.readdlm("$path/$adi_rootname_twoway.out",skipstart=8)
# # header_library_owens_stiff_one_way = DelimitedFiles.readdlm("$path/$adi_rootname_twoway.out",header=true,skipstart=6)[2]

# # headerNames = ["AB4N006Vrel","AB1N002Alpha","AB4N006STVx","AB4N006STVy","AB4N006STVz","B1AeroFx","B1AeroFy","B1AeroFz"]
# bladenum= "1"
# node = "018"
# bladenum2= bladenum#"1"
# bladenum3 = 1
# # headerNames1 = ["B$(bladenum)AeroFxg","B$(bladenum)AeroFyg","B$(bladenum)AeroFzg","B$(bladenum)AeroMxg","B$(bladenum)AeroMyg","B$(bladenum)AeroMzg"] #["AB$(bladenum)N$(node)Vrel","AB$(bladenum)N$(node)Alpha","AB$(bladenum)N$(node)STVx","AB$(bladenum)N$(node)STVy","AB$(bladenum)N$(node)STVz","AB$(bladenum)N$(node)Fn","AB$(bladenum)N$(node)Ft","B$(bladenum)AeroFx","B$(bladenum)AeroFy","B$(bladenum)AeroFz"]
# # headerNames2 = ["B$(bladenum2)AeroFxg","B$(bladenum2)AeroFyg","B$(bladenum2)AeroFzg","B$(bladenum2)AeroMxg","B$(bladenum2)AeroMyg","B$(bladenum2)AeroMzg"] #["AB$(bladenum2)N$(node)Vrel","AB$(bladenum2)N$(node)Alpha","AB$(bladenum2)N$(node)STVx","AB$(bladenum2)N$(node)STVy","AB$(bladenum2)N$(node)STVz","AB$(bladenum2)N$(node)Fn","AB$(bladenum2)N$(node)Ft","B$(bladenum2)AeroFx","B$(bladenum2)AeroFy","B$(bladenum2)AeroFz"]

# headerNames1 = headerNames2 = ["AB$(bladenum)N$(node)Ft"]#B1Azimuth"]#AB1N001Alpha"]#["B$(bladenum)AeroFxg"]

# plotdata1_standalone = zeros(length(headerNames1),length(standalone_AD[:,1]))
# plotdata1_library = zeros(length(headerNames1),length(standalone_AD[:,1]))
# time_library = zeros(length(standalone_AD[:,1]))

# plotdata1_library_stiff_one_way = zeros(length(headerNames1),length(standalone_AD[:,1]))
# time_library_stiff_one_way = zeros(length(standalone_AD[:,1]))


# # Since we did correction stepping with the library, we need to filter out each unique timestep, preferrably the aeroDistLoadsArrayTime
# for (idx_t_standalone, time) in enumerate(standalone_AD[1:end-1,1])

#     idx_t_library = findlast(x->x==time,library_ADside[:,1])
#     time_library[idx_t_standalone] = library_ADside[idx_t_library,1]
    
#     # idx_t_library_stiff_one_way = findlast(x->x==time,library_ADside_owens_stiff_one_way[:,1])
#     # time_library_stiff_one_way[idx_t_standalone] = library_ADside_owens_stiff_one_way[idx_t_library_stiff_one_way,1]

#     # println(idx_t_library)
#     #AOA, structural velocity, and structural loads
#     for (ihead,header_name) in enumerate(headerNames1)
#         header_idx_standalone = findfirst(x->x==header_name,header_standalone)
        
#         header_idx_library = findfirst(x->x==headerNames2[ihead],header_library)
        
#         # println("$header_name $header_idx_standalone $header_idx_library")

#         # @test isapprox(standalone_AD[idx_t_standalone,header_idx[2]],library_ADside[idx_t_library,header_idx[2]])
#         plotdata1_standalone[ihead,idx_t_standalone] = standalone_AD[idx_t_standalone,header_idx_standalone[2]]
#         plotdata1_library[ihead,idx_t_standalone] = library_ADside[idx_t_library,header_idx_library[2]]

#         # header_idx_library_owens_stiff_one_way = findfirst(x->x==headerNames2[ihead],header_library_owens_stiff_one_way)
#         # plotdata1_library_stiff_one_way[ihead,idx_t_standalone] = library_ADside_owens_stiff_one_way[idx_t_library_stiff_one_way,header_idx_library_owens_stiff_one_way[2]]
        
#     end
# end

# import PyPlot
# PyPlot.pygui(true)
# PyPlot.rc("figure", figsize=(15, 15))
# PyPlot.close("all")

# for (ihead,header_name) in enumerate(headerNames1)

#     PyPlot.figure()
#     PyPlot.plot(standalone_AD[1:end-1,1],plotdata1_standalone[ihead,1:end-1],"k-",label="standalone")
#     PyPlot.plot(time_library[1:end-1].+dt*2,plotdata1_library[ihead,1:end-1],color=plot_cycle[1],":",label="aerodyn side library")

#     # PyPlot.plot(time_library[1:end-1].+dt,fx_owensside[1:end-1],color=plot_cycle[1],"--",label="xowens side library")

    
#     PyPlot.ylabel(header_name)
#     PyPlot.xlabel("Time (s)")
#     PyPlot.legend()
#     # PyPlot.ylim([-4000,4000])
# end
