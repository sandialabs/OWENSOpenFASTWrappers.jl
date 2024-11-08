# Example usage the Julia wrapper for the AeroDyn DLL

import FLOWMath
import OWENSAero
import OWENS

import DelimitedFiles
path = splitdir(@__FILE__)[1]
# import OWENSOpenFASTWrappers
using Test

import PyPlot
PyPlot.pygui(true)
PyPlot.rc("figure", figsize=(15, 15))
PyPlot.rc("font", size=10.0)
PyPlot.rc("lines", linewidth=1.5)
PyPlot.rc("lines", markersize=4.0)
PyPlot.rc("legend", frameon=true)
PyPlot.rc("axes.spines", right=false, top=false)
PyPlot.rc("figure.subplot", left=.21, bottom=.17, top=0.9, right=.9)
PyPlot.rc("figure",max_open_warning=500)
# PyPlot.rc("axes", prop_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

# println("Calculate/Set up sectional properties")
# include("$path/setup_XFlow_25kW.jl")
# XFlow specific mesh generator
# include("$path/setupOWENShawt.jl")
include("$path/../../src/OWENSOpenFASTWrappers.jl")

# # Run the standalone aerodyn
# run(`$path/../../../../openfast/build/modules/aerodyn/aerodyn_driver $path/HAWT_standalone_test.dvr`)

adi_lib = nothing#"$path/../../../../openfast/build/modules/aerodyn/libaerodyn_inflow_c_binding" #change this to match your local path of the AeroDyn DLL
# adi_lib = "$path/../..//openfast/build/modules/AeroDyn/libAeroDyn_c_binding" #change this to match your local path of the AeroDyn DLL

# output files from ADI
adi_rootname_direct = "$path/HAWT_OWENS_AD15"
adi_rootname_FSI = "$path/HAWT_OWENS_AD15FSI"

num_corrections = 0


dt = 0.01
t_initial = 0.0
t_max = 0.5
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
shapeZ = zero(shapeX)

mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
    FSI_twr, FSI_bld,bld_precompinput,
    bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
    twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,
    mass_breakout_blds,mass_breakout_twr,bladeIdx,bladeElem,system,assembly,sections = OWENS.setupOWENShawt(OWENSAero,path;
    rho,
    Nslices,
    RPM,
    Vinf,
    eta = 0.5,
    B,
    H,
    R,
    hubR = 0.0,
    shapeZ,
    shapeX,
    NuMad_geom_xlscsv_file_twr = nothing,
    NuMad_mat_xlscsv_file_twr = nothing,
    NuMad_geom_xlscsv_file_bld = nothing,
    NuMad_mat_xlscsv_file_bld = nothing,
    stack_layers_bld = nothing,
    Ht = 1e-6,
    ntelem = 10, #tower elements
    nbelem = 18, #blade elements
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


PyPlot.figure()
PyPlot.plot(mymesh.x,mymesh.y,"b-")
 for myi = 1:length(mymesh.y)
     PyPlot.text(mymesh.x[myi].+rand()/30,mymesh.y[myi].+rand()/30,"$myi",ha="center",va="center")
     PyPlot.draw()
     sleep(0.1)
 end
PyPlot.xlabel("x")
PyPlot.ylabel("y")
# PyPlot.axis("equal")

# Insert Sectional Properties

NREL5MWprops = DelimitedFiles.readdlm("$path/NRELOffshrBsline5MW_Blade_51.csv",',',Float64,skipstart=1)

NREL5MW_r = NREL5MWprops[:,1]
NREL5MW_BlFract = NREL5MWprops[:,2]
NREL5MW_PitchAxis = NREL5MWprops[:,3]
NREL5MW_StrcTwst = NREL5MWprops[:,4]
NREL5MW_BMassDen = NREL5MWprops[:,5]
NREL5MW_FlpStff = NREL5MWprops[:,6]
NREL5MW_EdgStff = NREL5MWprops[:,7]

for iel = 1:length(myel.props)
    # myel.props[iel].ac
    # myel.props[iel].twist
    myel.props[iel].rhoA .= NREL5MW_BMassDen[10]
    myel.props[iel].EIyy .= NREL5MW_FlpStff[10]
    myel.props[iel].EIzz .= NREL5MW_EdgStff[10]
    myel.props[iel].GJ .= NREL5MW_EdgStff[10]*5
    myel.props[iel].EA .= NREL5MW_EdgStff[10]*5
    # myel.props[iel].rhoIyy
    # myel.props[iel].rhoIzz
    # myel.props[iel].rhoJ
end

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


# OWENSOpenFASTWrappers.writeADinputFile(ad_input_file,blade_filenames,airfoil_filenames,OLAF_filename)

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
#     OWENSOpenFASTWrappers.writeADbladeFile(filename;NumBlNds,BlSpn,BlCrvAC,BlSwpAC,BlCrvAng,BlTwist,BlChord,BlAFID)
# end

# OWENSOpenFASTWrappers.writeOLAFfile(OLAF_filename;nNWPanel=200,nFWPanels=10)

# OWENSOpenFASTWrappers.writeIWfile(10.0,ifw_input_file;turbsim_filename=nothing)

# Run AeroDyn
OWENSOpenFASTWrappers.setupTurb(adi_lib,ad_input_file,ifw_input_file,adi_rootname_direct,[shapeX],[shapeZ],[B],[Ht],[mymesh],[myort],[bladeIdx],[bladeElem];
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
        nacPos      = [[0,0,Ht]],  
        refPos = [[0,0.0,0.0]],
        isHAWT = true
        )

# Time marching
ForceValHist = zeros(length(ts[1:end-1]),Int(mymesh.numNodes*6))
AziHist = zeros(length(ts[1:end-1]))
Fxhist = zeros(mymesh.numNodes,length(ts[1:end-1]))
Fyhist = zeros(mymesh.numNodes,length(ts[1:end-1]))
Fzhist = zeros(mymesh.numNodes,length(ts[1:end-1]))
Mxhist = zeros(mymesh.numNodes,length(ts[1:end-1]))
Myhist = zeros(mymesh.numNodes,length(ts[1:end-1]))
Mzhist = zeros(mymesh.numNodes,length(ts[1:end-1]))

for (tidx, t) in enumerate(ts[1:end-1])
    # println("time $t")
    for correction = 1:num_corrections+1
        # println("correction $correction")

        #mapAD15(t,azi_j,mesh,advanceAD15)
        azi_j = omega*(t+dt) # rad/s * s
        AziHist[tidx] = azi_j
        u_j     = zeros(mymesh.numNodes*6)
        udot_j  = zeros(mymesh.numNodes*6)
        uddot_j = zeros(mymesh.numNodes*6)
        hubPos      = [0,0,Ht]                      # m
        hubAngle    = [0,0,0]                       # rad
        rotvel = omega
        hubVel = [0,0,0,0.0,0,0]#zeros(6)
        hubAcc = zeros(6) #TODO: may eventually need this for MHK?

        OWENSOpenFASTWrappers.deformAD15([u_j],
        [udot_j],
        [uddot_j],
        [azi_j],
        [omega],
        [zero(omega)],
        [hubPos],
        [hubAngle],
        [hubVel],
        [hubAcc])

        n_steps,Fx,Fy,Fz,Mx,My,Mz = OWENSOpenFASTWrappers.advanceAD15(t,[mymesh],[azi_j])
        Fxhist[:,tidx] = Fx[1][:,1]
        Fyhist[:,tidx] = Fy[1][:,1]
        Fzhist[:,tidx] = Fz[1][:,1]
        Mxhist[:,tidx] = Mx[1][:,1]
        Myhist[:,tidx] = My[1][:,1]
        Mzhist[:,tidx] = Mz[1][:,1]

        # Map loads over from advanceTurb
        for i=1:mymesh.numNodes
            ForceValHist[tidx,(i-1)*6+1] = Fx[1][i,1]
            ForceValHist[tidx,(i-1)*6+2] = Fy[1][i,1]
            ForceValHist[tidx,(i-1)*6+3] = Fz[1][i,1]
            ForceValHist[tidx,(i-1)*6+4] = Mx[1][i,1]
            ForceValHist[tidx,(i-1)*6+5] = My[1][i,1]
            ForceValHist[tidx,(i-1)*6+6] = Mz[1][i,1]
        end
        # DOFs are sequential through all nodes
        ForceDof=collect(1:1:mymesh.numNodes*6)
    end

end
println("End Aerodynamics")
OWENSOpenFASTWrappers.endTurb()

######################################
#### Plot Comparison
#######################################

# Load Standalone Data
standalone_AD = DelimitedFiles.readdlm("$path/HAWT_standalone_test.out",skipstart=8)
header_standalone = DelimitedFiles.readdlm("$path/HAWT_standalone_test.out",header=true,skipstart=6)[2]

# Compare with Current Data on the AeroDyn Side
library_ADside_direct = DelimitedFiles.readdlm("$adi_rootname_direct.out",skipstart=8)
header_library_direct = DelimitedFiles.readdlm("$adi_rootname_direct.out",header=true,skipstart=6)[2]

bladenum = 2
nodeint = 5
node = "00$nodeint"
headerNames1 = ["RtAeroFxh","RtAeroFyh","RtAeroFzh","AB$(bladenum)N$(node)Fxi","AB$(bladenum)N$(node)Fyi","AB$(bladenum)N$(node)Fzi"]# ","AB$(bladenum)N$(node)STVx","AB$(bladenum)N$(node)STVy","AB$(bladenum)N$(node)STVz","AB$(bladenum)N$(node)Vx","AB$(bladenum)N$(node)Vy","AB$(bladenum)N$(node)Alpha","B$(bladenum)AeroFxi","B$(bladenum)AeroFyi"]

# Extract the blade and node reaction force
meshNodeidx_plt = Int(mymesh.structuralNodeNumbers[bladenum,nodeint])
locallen = sqrt((mymesh.x[meshNodeidx_plt+1]-mymesh.x[meshNodeidx_plt])^2 + (mymesh.y[meshNodeidx_plt+1]-mymesh.y[meshNodeidx_plt])^2 + (mymesh.z[meshNodeidx_plt+1]-mymesh.z[meshNodeidx_plt])^2)

strtspnidx = bladeIdx[bladenum,1]
endspnidx = bladeIdx[bladenum,2]

spanlen = sqrt((mymesh.x[endspnidx]-mymesh.x[strtspnidx])^2+(mymesh.y[endspnidx]-mymesh.y[strtspnidx])^2+(mymesh.z[endspnidx]-mymesh.z[strtspnidx])^2)


plotdata1_standalone = zeros(length(headerNames1),length(standalone_AD[:,1]))

plotdata1_direct = zeros(length(headerNames1),length(standalone_AD[:,1]))
time_library_direct = zeros(length(standalone_AD[:,1]))


# Since we did correction stepping with the library, we need to filter out each unique timestep, preferrably the aeroDistLoadsArrayTime
for (idx_t_standalone, time) in enumerate(standalone_AD[1:end-1,1])

    idx_t_direct = findlast(x->x==time,library_ADside_direct[:,1])
    time_library_direct[idx_t_standalone] = library_ADside_direct[idx_t_direct,1]

    for (ihead,header_name) in enumerate(headerNames1)
        header_idx_standalone = findfirst(x->x==header_name,header_standalone)
        
        header_idx_library = findfirst(x->x==headerNames1[ihead],header_library_direct)
        
        # println("$header_name $header_idx_standalone $header_idx_library")
        plotdata1_standalone[ihead,idx_t_standalone] = standalone_AD[idx_t_standalone,header_idx_standalone[2]]
        plotdata1_direct[ihead,idx_t_standalone] = library_ADside_direct[idx_t_direct,header_idx_library[2]]
    end
end

# Now get the mapped forces:
fm_direct = zeros(6,length(ForceValHist[:,1]))

for i_time = 1:length(ForceValHist[:,1])
    # Get values from the library, but at the owens loads side.
    for idof = 1:6
        dof_end = maximum(bladeIdx[bladenum,:])*6-6+idof
        dof_start = minimum(bladeIdx[bladenum,:])*6-6+idof
        fm_direct[idof,i_time] = sum(ForceValHist[i_time,dof_start:6:dof_end])
    end
end

# DCM = OWENSOpenFASTWrappers.createGeneralTransformationMatrix(ang1,angle_axes)

fy_directcall = Fyhist[strtspnidx+nodeint-1,:].*cos.(AziHist).-Fxhist[strtspnidx+nodeint-1,:].*sin.(AziHist)
fz_directcall = Fyhist[strtspnidx+nodeint-1,:].*sin.(AziHist).+Fxhist[strtspnidx+nodeint-1,:].*cos.(AziHist)
fx_directcall = Fzhist[strtspnidx+nodeint-1,:]
strtspnidx2 = strtspnidx+nodeint-1
endspnidx2 = strtspnidx2+1
localspanlen = sqrt((mymesh.x[endspnidx2]-mymesh.x[strtspnidx2])^2+(mymesh.y[endspnidx2]-mymesh.y[strtspnidx2])^2+(mymesh.z[endspnidx2]-mymesh.z[strtspnidx2])^2)

mx_directcall = fm_direct[4,:].*cos.(AziHist).-fm_direct[5,:].*sin.(AziHist)
my_directcall = fm_direct[4,:].*sin.(AziHist).+fm_direct[5,:].*cos.(AziHist)

totalFxh_direct = [sum(Fxhist[:,itime]) for itime = 1:length(Fxhist[1,:])]
totalFyh_direct = [sum(Fyhist[:,itime]) for itime = 1:length(Fxhist[1,:])]
totalFzh_direct = [sum(Fzhist[:,itime]) for itime = 1:length(Fxhist[1,:])]
# Mxhist
# Myhist
# Mzhist

PyPlot.rc("figure", figsize=(4.5, 3))
PyPlot.close("all")

for (ihead,header_name) in enumerate(headerNames1)

    PyPlot.figure(ihead)
    PyPlot.plot(standalone_AD[2:end-1,1],plotdata1_standalone[ihead,2:end-1],"k-",label="Standalone OLAF Output")
    PyPlot.plot(time_library_direct[2:end-1].+dt,plotdata1_direct[ihead,1:end-2],color=plot_cycle[1],":",label="Direct Library OLAF Side")
    PyPlot.ylabel("$header_name")
    PyPlot.xlabel("Time (s)")

    if contains(header_name,"Fxi")
        PyPlot.plot(ts[2:end],fx_directcall/localspanlen,color=plot_cycle[1],"-",label="Direct Library OWENS Side")
    end
    if contains(header_name,"Fyi")
        PyPlot.plot(ts[2:end],fy_directcall/localspanlen,color=plot_cycle[1],"-",label="Direct Library OWENS Side")
    end
    if contains(header_name,"Fzi")
        PyPlot.plot(ts[2:end],fz_directcall/localspanlen,color=plot_cycle[1],"-",label="Direct Library OWENS Side")
    end

    if contains(header_name,"Fxh")
        PyPlot.plot(ts[2:end],totalFzh_direct,color=plot_cycle[1],"-",label="Direct Library OWENS Side")
    end
    if contains(header_name,"Fyh")
        PyPlot.plot(ts[2:end],totalFyh_direct,color=plot_cycle[1],"-",label="Direct Library OWENS Side")
    end
    if contains(header_name,"Fzh")
        PyPlot.plot(ts[2:end],totalFxh_direct,color=plot_cycle[1],"-",label="Direct Library OWENS Side")
    end
    
    PyPlot.legend()
    # PyPlot.xlim([0.4,0.5])
    PyPlot.savefig("$(path)/$header_name.pdf",transparent = true)
end


##########################################
############## FSI Simulation ############
##########################################
mymesh.hubPos = [0,0,Ht] #TODO: make this more integrated or automatically set within
OWENSOpenFASTWrappers.setupTurb(adi_lib,ad_input_file,ifw_input_file,adi_rootname_FSI,[shapeX],[shapeZ],[B],[Ht],[mymesh],[myort],[bladeIdx],[bladeElem];
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
        nacPos      = [[0,0,Ht]],  
        refPos = [[0,0.0,0.0]],
        isHAWT = true
        )

# Routine for getting aero forces from aD15
aeroForcesAD15(t,azi) = OWENS.mapAD15(t,azi,[mymesh],OWENSOpenFASTWrappers.advanceAD15;alwaysrecalc=true,verbosity=1)

model = OWENS.Inputs(;analysisType = "TNB",
    tocp = [0.0,100000.1],
    Omegaocp = [RPM,RPM] ./ 60,
    turbineStartup = 0,
    generatorOn = false,
    useGeneratorFunction = false,
    numTS = numTS,
    delta_t = dt,
    aeroLoadsOn = 2.0,        # 1: one way; 1.5: aero once in iteration in 2 way; 2: two-way
    topsideOn = true,
    interpOrder = 2,
    AD15On=true)
    model.iteration_parameters.MAXITER=20   # temporary for testing

pBC = [ 1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0]

feamodel = OWENS.FEAModel(;analysisType = "TNB",
    dataOutputFilename = "none",
    joint = myjoint,
    platformTurbineConnectionNodeNumber = 1,
    pBC,
    nlOn = false,
    gravityOn = false, # for HAWT, gravity would be in the +x direction in the OWENS FOR
    numNodes = mymesh.numNodes,
    RayleighAlpha = 500.05,
    RayleighBeta = 500.05,
    iterationType = "DI")

println("Running Unsteady")
t_FSI, aziHist_FSI,OmegaHist_FSI,OmegaDotHist_FSI,gbHist_FSI,gbDotHist_FSI,gbDotDotHist_FSI,
ForceValHist_FSI,FTwrBsHist_FSI,genTorque_FSI,genPower_FSI,torqueDriveShaft_FSI,uHist_FSI,
uHist_prp_FSI,epsilon_x_hist_FSI,epsilon_y_hist_FSI,epsilon_z_hist_FSI,kappa_x_hist_FSI,
kappa_y_hist_FSI,kappa_z_hist_FSI,FPtfmHist_FSI,FHydroHist_FSI,FMooringHist_FSI,
topFexternal_hist_FSI,rbDataHist_FSI = OWENS.Unsteady_Land(model;topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForcesAD15,
deformAero=OWENSOpenFASTWrappers.deformAD15,system=system,assembly=assembly) #,meshcontrolfunction=mymeshcontrolfunction2,userDefinedGenerator=userDefinedGenerator,

println("Saving VTK time domain files")
OWENS.OWENSFEA_VTK("$path/vtk/fsi1",t_FSI,uHist_FSI,system,assembly,sections;scaling=1,azi=aziHist_FSI)


println("End Aerodynamics")
OWENSOpenFASTWrappers.endTurb()

######################################
#### Plot Comparison
#######################################

# Load Standalone Data
standalone_AD = DelimitedFiles.readdlm("$path/HAWT_standalone_test.out",skipstart=8)
header_standalone = DelimitedFiles.readdlm("$path/HAWT_standalone_test.out",header=true,skipstart=6)[2]

# Compare with Current Data on the AeroDyn Side
library_ADside_direct = DelimitedFiles.readdlm("$adi_rootname_direct.out",skipstart=8)
header_library_direct = DelimitedFiles.readdlm("$adi_rootname_direct.out",header=true,skipstart=6)[2]

# Compare with FSI Data on the AeroDyn Side
library_ADside_FSI = DelimitedFiles.readdlm("$adi_rootname_FSI.out",skipstart=8)
header_library_FSI = DelimitedFiles.readdlm("$adi_rootname_FSI.out",header=true,skipstart=6)[2]


# Extract the blade and node reaction force
meshNodeidx_plt = Int(mymesh.structuralNodeNumbers[bladenum,nodeint])
locallen = sqrt((mymesh.x[meshNodeidx_plt+1]-mymesh.x[meshNodeidx_plt])^2 + (mymesh.y[meshNodeidx_plt+1]-mymesh.y[meshNodeidx_plt])^2 + (mymesh.z[meshNodeidx_plt+1]-mymesh.z[meshNodeidx_plt])^2)

plotdata1_standalone = zeros(length(headerNames1),length(standalone_AD[:,1]))

plotdata1_direct = zeros(length(headerNames1),length(standalone_AD[:,1]))
time_library_direct = zeros(length(standalone_AD[:,1]))

plotdata1_FSI = zeros(length(headerNames1),length(standalone_AD[:,1]))
time_library_FSI = zeros(length(standalone_AD[:,1]))


# Since we did correction stepping with the library, we need to filter out each unique timestep, preferrably the aeroDistLoadsArrayTime
for (idx_t_standalone, time) in enumerate(standalone_AD[1:end-1,1])

    idx_t_direct = findlast(x->x==time,library_ADside_direct[:,1])
    time_library_direct[idx_t_standalone] = library_ADside_direct[idx_t_direct,1]

    idx_t_FSI = findlast(x->x==time,library_ADside_FSI[:,1])
    time_library_FSI[idx_t_standalone] = library_ADside_FSI[idx_t_FSI,1]

    for (ihead,header_name) in enumerate(headerNames1)
        header_idx_standalone = findfirst(x->x==header_name,header_standalone)
        
        header_idx_library = findfirst(x->x==headerNames1[ihead],header_library_direct)

        header_idx_FSI = findfirst(x->x==headerNames1[ihead],header_library_FSI)
        
        # println("$header_name $header_idx_standalone $header_idx_library")
        plotdata1_standalone[ihead,idx_t_standalone] = standalone_AD[idx_t_standalone,header_idx_standalone[2]]
        plotdata1_direct[ihead,idx_t_standalone] = library_ADside_direct[idx_t_direct,header_idx_library[2]]
        plotdata1_FSI[ihead,idx_t_standalone] = library_ADside_FSI[idx_t_FSI,header_idx_library[2]]
    end
end


# # Now get the mapped forces:
# fm_FSI = zeros(6,length(ForceValHist_FSI[:,1]))

# for i_time = 1:length(ForceValHist_FSI[:,1])
#     for idof = 1:6
#         dof_end = maximum(bladeIdx[bladenum,:])*6-6+idof
#         dof_start = minimum(bladeIdx[bladenum,:])*6-6+idof
#         # fm_FSI[idof,i_time] = sum(ForceValHist_FSI[i_time,dof_start:6:dof_end])
#     end
# end


# totalFxh_FSI = [sum(ForceValHist_FSI[itime,1:6:end]) for itime = 1:length(ForceValHist_FSI[:,1])]
# totalFyh_FSI = [sum(ForceValHist_FSI[itime,2:6:end]) for itime = 1:length(ForceValHist_FSI[:,1])]
# totalFzh_FSI = [sum(ForceValHist_FSI[itime,3:6:end]) for itime = 1:length(ForceValHist_FSI[:,1])]
# Mxhist
# Myhist
# Mzhist

# fx_FSIcall = fm_FSI[1,:].*cos.(aziHist_FSI).-fm_FSI[2,:].*sin.(aziHist_FSI)
# fy_FSIcall = fm_FSI[1,:].*sin.(aziHist_FSI).+fm_FSI[2,:].*cos.(aziHist_FSI)
# fz_FSIcall = fm_FSI[3,:]
# mx_FSIcall = fm_FSI[4,:].*cos.(aziHist_FSI).-fm_FSI[5,:].*sin.(aziHist_FSI)
# my_FSIcall = fm_FSI[4,:].*sin.(aziHist_FSI).+fm_FSI[5,:].*cos.(aziHist_FSI)

# fxi_node_FSI = (ForceValHist_FSI[:,(meshNodeidx_plt)*6-6+1]./locallen).*cos.(AziHist) .- (ForceValHist_FSI[:,(meshNodeidx_plt)*6-6+2]./locallen).*sin.(AziHist)
# fyi_node_FSI = (ForceValHist_FSI[:,(meshNodeidx_plt)*6-6+1]./locallen).*sin.(AziHist) .+ (ForceValHist_FSI[:,(meshNodeidx_plt)*6-6+2]./locallen).*cos.(AziHist)
# fzi_node_FSI = (ForceValHist_FSI[:,(meshNodeidx_plt)*6-6+3]./locallen)

Fxhist = ForceValHist_FSI[:,1:6:end]
Fyhist = ForceValHist_FSI[:,2:6:end]
Fzhist = ForceValHist_FSI[:,3:6:end]

fy_directcall = Fyhist[:,strtspnidx+nodeint-1].*cos.(AziHist).-Fxhist[:,strtspnidx+nodeint-1].*sin.(AziHist)
fz_directcall = Fyhist[:,strtspnidx+nodeint-1].*sin.(AziHist).+Fxhist[:,strtspnidx+nodeint-1].*cos.(AziHist)
fx_directcall = Fzhist[:,strtspnidx+nodeint-1]

fy_FSIcall = fy_directcall#Fyhist_FSI[:,strtspnidx+nodeint-1].*cos.(AziHist_FSI).-Fxhist_FSI[:,strtspnidx+nodeint-1].*sin.(AziHist_FSI)
fz_FSIcall = fz_directcall#Fyhist_FSI[:,strtspnidx+nodeint-1].*sin.(AziHist_FSI).+Fxhist_FSI[:,strtspnidx+nodeint-1].*cos.(AziHist_FSI)
fx_FSIcall = fx_directcall#Fzhist_FSI[:,strtspnidx+nodeint-1]

PyPlot.close("all")

for (ihead,header_name) in enumerate(headerNames1)

    PyPlot.figure(ihead)
    PyPlot.plot(time_library_FSI[2:end-1].+dt,plotdata1_FSI[ihead,1:end-2],color=plot_cycle[2],":",label="Stiff OLAF Side")
    PyPlot.ylabel("$header_name")
    PyPlot.xlabel("Time (s)")

    if bladenum>3
        plotname="Strut"
        localbladenum = bladenum-B
    else
        plotname="Blade"
        localbladenum = bladenum
    end

    if contains(header_name,"Fxi")
        PyPlot.plot(ts[2:end],fx_FSIcall,color=plot_cycle[2],"-",label="FSI OWENS Side")
        PyPlot.ylabel("$plotname $localbladenum Node $nodeint Global X Force (N)")
    end
    if contains(header_name,"Fyi")
        PyPlot.plot(ts[2:end],fy_FSIcall,color=plot_cycle[2],"-",label="FSI OWENS Side")
        PyPlot.ylabel("$plotname $localbladenum Node $nodeint Global Y Force (N)")
    end
    if contains(header_name,"Fzi")
        PyPlot.plot(ts[2:end],fz_FSIcall,color=plot_cycle[2],"-",label="FSI OWENS Side")
        PyPlot.ylabel("$plotname $localbladenum Node $nodeint Global Z Force (N)")
    end

    if contains(header_name,"Fxh")
        PyPlot.plot(ts[2:end],ForceValHist_FSI[:,3],color=plot_cycle[2],"-",label="Stiff OWENS Side")
        PyPlot.ylabel("Hub X (Windward) Force (N)")
    end
    if contains(header_name,"Fyh")
        PyPlot.plot(ts[2:end],ForceValHist_FSI[:,2],color=plot_cycle[2],"-",label="Stiff OWENS Side")
        PyPlot.ylabel("Hub Y (Crosswind at t=0) Force (N)")
    end
    if contains(header_name,"Fzh")
        PyPlot.plot(ts[2:end],-ForceValHist_FSI[:,1],color=plot_cycle[2],"-",label="Stiff OWENS Side")
        PyPlot.ylabel("Hub Z (Vertical at t=0) Force (N)")
    end
    
    PyPlot.legend()
    # PyPlot.xlim([0.4,0.5])
    PyPlot.savefig("$(path)/$header_name.pdf",transparent = true)
end