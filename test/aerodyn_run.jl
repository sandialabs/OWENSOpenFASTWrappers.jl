# import PyPlot
# PyPlot.pygui(true)
# PyPlot.rc("figure", figsize=(15, 15))
# PyPlot.rc("font", size=10.0)
# PyPlot.rc("lines", linewidth=1.5)
# PyPlot.rc("lines", markersize=4.0)
# PyPlot.rc("legend", frameon=true)
# PyPlot.rc("axes.spines", right=false, top=false)
# PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
# PyPlot.rc("figure",max_open_warning=500)
# # PyPlot.rc("axes", prop_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
# plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]
# PyPlot.close("all")

import FLOWMath
import DelimitedFiles
import LinearAlgebra
path = splitdir(@__FILE__)[1]
import OWENSOpenFASTWrappers
using Test

# Run the standalone aerodyn
run(`$(OWENSOpenFASTWrappers.turbsim()) $path/AeroDynInputs/DLC1_1Vinf10.0.inp`)
run(`$(OWENSOpenFASTWrappers.aerodyn_driver()) $path/AeroDynInputs/HVAWT_standalone_test.dvr`)

include("$path/AeroDynInputs/meshdeps.jl")

adi_lib = nothing #change this to match your local path of the AeroDyn DLL

# output files from ADI
adi_rootname_direct = "$path/ADI_OWENS_direct"
adi_rootname_stiff = "$path/ADI_OWENS_stiff"

num_corrections = 0

AD15On = true
dt = 0.01
t_initial = 0.0
t_max = 0.5
ts = collect(t_initial:dt:t_max)
numTS = length(ts)

saveName = "$path/vtk/OWENS_test"
rho = 1.225
mu = 1.4639E-05 
Nslices = 30
RPM = 50.0
Ht = 10.0
B = Nbld = 3
H = 10.0#-1.0   # m 
R = 10.0/2     # m
Area = H*R*2
Vinf = 10.0 #m/s
chord = 0.5*ones(Nslices)
omega = RPM / 60 * 2 * pi
blade_joint_angle_Degrees = 0.0
Hub_Height = Ht

shapeY = collect(LinRange(0,H,Nslices+1))
shapeX = R.*vec(ones(Nslices+1,1))

println("Create Mesh")

mymesh, myort, myjoint, bladeIdx, bladeElem = create_mesh_struts(;Htwr_base = Ht-H/2,
    Htwr_blds = H/2,
    Hbld = H, #blade height
    R = R, # m bade radius
    AD15hubR = 1.0,
    nblade = 3,
    ntelem = 10, #tower elements
    nbelem = 10, #blade elements
    nselem = 10,
    strut_twr_mountpoint = [0.46,0.54], # This puts struts at top and bottom, as a fraction of the blade position
    strut_bld_mountpoint = [0.17,0.83], # This puts struts at bottom 0, mid 0.5, and top 1.0 as a fraction of the blade position
    bshapex = shapeX, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    bshapez = shapeY,
    bshapey = zero(shapeX), # but magnitude for this is relevant
    angularOffset = pi/2.0+0*pi/180, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    AD15_ccw = true,
    verbosity = 0, # 0 nothing, 1 basic, 2 lots: amount of printed information
    connectBldTips2Twr = false)


# PyPlot.figure()
# for icon = 1:length(mymesh.conn[:,1])
# idx1 = mymesh.conn[icon,1]
# idx2 = mymesh.conn[icon,2]
# PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"b.-")
# PyPlot.text3D(mymesh.x[idx1].+rand()/30,mymesh.y[idx1].+rand()/30,mymesh.z[idx1].+rand()/30,"$idx1",ha="center",va="center")
# # sleep(0.1)
# end

# for ijoint = 1:length(myjoint[:,1])
# idx2 = Int(myjoint[ijoint,2])
# idx1 = Int(myjoint[ijoint,3])
# PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"y.-")
# sleep(0.1)
# end


    

# Offset mesh by hub position
# mymesh.z .-= Hub_Height

# Create AeroDyn Files

# input files for ADI
ad_input_file="$path/AeroDynInputs/ADInputFile_oneTurbine.dat"
ifw_input_file="$path/AeroDynInputs/IW.dat"
blade_filename="$path/AeroDynInputs/blade.dat"
lower_strut_filename="$path/AeroDynInputs/lower_arm.dat"
upper_strut_filename="$path/AeroDynInputs/upper_arm.dat"
airfoil_filenames = "$path/airfoils/NACA_0018_AllRe.dat"
OLAF_filename = "$path/AeroDynInputs/OLAF.dat"

# Run AeroDyn
OWENSOpenFASTWrappers.setupTurb(adi_lib,ad_input_file,ifw_input_file,adi_rootname_direct,[shapeX],[shapeY],[B],[Ht],[mymesh],[myort],[bladeIdx],[bladeElem];
        rho     = rho,
        adi_dt  = dt,
        adi_tmax= t_max,
        omega   = [omega],
        adi_wrOuts = 1,     # write output file [0 none, 1 txt, 2 binary, 3 both]
        adi_DT_Outs = dt,    # output frequency
        hubPos = [[0,0,Hub_Height]],
        hubAngle = [[0,0,0]],
        isHAWT = false,
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
        hubPos      = [0,0,Hub_Height]                      # m
        hubAngle    = [0,0,0]                       # rad
        rotvel = omega
        hubVel = [0,0,0,0,0,rotvel]#zeros(6)
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
standalone_AD = DelimitedFiles.readdlm("$path/AeroDynInputs/HVAWT_standalone_test.out",skipstart=8)
header_standalone = DelimitedFiles.readdlm("$path/AeroDynInputs/HVAWT_standalone_test.out",header=true,skipstart=6)[2]

# Compare with Current Data on the AeroDyn Side
library_ADside_direct = DelimitedFiles.readdlm("$adi_rootname_direct.out",skipstart=8)
header_library_direct = DelimitedFiles.readdlm("$adi_rootname_direct.out",header=true,skipstart=6)[2]

bladenum = 4
node = "005"
headerNames1 = ["RtAeroFxh","RtAeroFyh","RtAeroFzh","AB$(bladenum)N$(node)STVx","AB$(bladenum)N$(node)STVy","AB$(bladenum)N$(node)STVz","AB$(bladenum)N$(node)Vx","AB$(bladenum)N$(node)Vy","AB$(bladenum)N$(node)Alpha","AB$(bladenum)N$(node)Fx","AB$(bladenum)N$(node)Fy"] 

plotdata1_standalone = zeros(length(headerNames1),length(standalone_AD[:,1]))
# plotdata1_stiff = zeros(length(headerNames1),length(standalone_AD[:,1]))
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
Fxh_direct = zeros(length(ForceValHist[:,1]))
Fyh_direct = zeros(length(ForceValHist[:,1]))
Fzh_direct = zeros(length(ForceValHist[:,1]))

for i_time = 1:length(ForceValHist[:,1])
    # Get values from the library, but at the owens loads side.
    Fxh_direct[i_time] = sum(Fxhist[:,i_time])
    Fyh_direct[i_time] = sum(Fyhist[:,i_time])
    Fzh_direct[i_time] = sum(Fzhist[:,i_time])
    for idof = 1:6
        dof_end = maximum(bladeIdx[bladenum,:])*6-6+idof
        dof_start = minimum(bladeIdx[bladenum,:])*6-6+idof
        fm_direct[idof,i_time] = sum(ForceValHist[i_time,dof_start:6:dof_end])
    end
end

# DCM = OWENSOpenFASTWrappers.createGeneralTransformationMatrix(ang1,angle_axes)

fx_directcall = fm_direct[1,:].*cos.(AziHist).-fm_direct[2,:].*sin.(AziHist)
fy_directcall = fm_direct[1,:].*sin.(AziHist).+fm_direct[2,:].*cos.(AziHist)
mx_directcall = fm_direct[4,:].*cos.(AziHist).-fm_direct[5,:].*sin.(AziHist)
my_directcall = fm_direct[4,:].*sin.(AziHist).+fm_direct[5,:].*cos.(AziHist)

# PyPlot.rc("figure", figsize=(4.5, 3))
# PyPlot.close("all")
percentTolDirect_Standalone = 14.0
percentTolADside_OWENSside = 4.0
mintol = 1e-4
common_spline_time = LinRange(standalone_AD[2,1],standalone_AD[end-1,1],100)
for (ihead,header_name) in enumerate(headerNames1)
    # println(ihead)
    # PyPlot.figure(ihead)
    # PyPlot.plot(standalone_AD[2:end-1,1],plotdata1_standalone[ihead,2:end-1],"k-",label="Standalone OLAF Output")
    # PyPlot.plot(time_library_direct[2:end-1].+dt,plotdata1_direct[ihead,1:end-2],color=plot_cycle[3],"+",label="Direct Library OLAF Side")
    # PyPlot.ylabel("$header_name")
    # PyPlot.xlabel("Time (s)")

    standalone_spline = FLOWMath.akima(standalone_AD[2:end-1,1],plotdata1_standalone[ihead,2:end-1],common_spline_time)
    direct_spline = FLOWMath.akima(time_library_direct[2:end-1].+dt,plotdata1_direct[ihead,1:end-2],common_spline_time)    

    atol=max(mintol,maximum(abs.(standalone_spline))*percentTolDirect_Standalone/100)
    for iel = 1:length(standalone_spline)
        @test isapprox(standalone_spline[iel],direct_spline[iel];atol)
    end

    if contains(header_name,"Fxi")
        wrapper_spline = FLOWMath.akima(ts[2:end],fx_directcall,common_spline_time)
        @test isapprox(wrapper_spline,direct_spline;atol=max(mintol,maximum(abs.(wrapper_spline))*percentTolADside_OWENSside/100))
        # PyPlot.plot(ts[2:end],fx_directcall,color=plot_cycle[1],".",label="Direct Library OWENS Side")        
    end
    if contains(header_name,"Fyi")
        wrapper_spline = FLOWMath.akima(ts[2:end],fy_directcall,common_spline_time)
        @test isapprox(wrapper_spline,direct_spline;atol=max(mintol,maximum(abs.(wrapper_spline))*percentTolADside_OWENSside/100))
        # PyPlot.plot(ts[2:end],fy_directcall,color=plot_cycle[1],".",label="Direct Library OWENS Side")        
    end

    if contains(header_name,"Fxh")
        wrapper_spline = FLOWMath.akima(ts[2:end],Fzh_direct,common_spline_time)
        @test isapprox(wrapper_spline,direct_spline;atol=max(mintol,maximum(abs.(wrapper_spline))*percentTolADside_OWENSside/100))
        # PyPlot.plot(ts[2:end],Fzh_direct,color=plot_cycle[1],".",label="Direct Library OWENS Side")        
    end
    if contains(header_name,"Fyh")
        wrapper_spline = FLOWMath.akima(ts[2:end],-Fyh_direct,common_spline_time)
        @test isapprox(wrapper_spline,direct_spline;atol=max(mintol,maximum(abs.(wrapper_spline))*percentTolADside_OWENSside/100))
        # PyPlot.plot(ts[2:end],-Fyh_direct,color=plot_cycle[1],".",label="Direct Library OWENS Side")        
    end
    if contains(header_name,"Fzh")
        wrapper_spline = FLOWMath.akima(ts[2:end],Fxh_direct,common_spline_time)
        @test isapprox(wrapper_spline,direct_spline;atol=max(mintol,maximum(abs.(wrapper_spline))*percentTolADside_OWENSside/100))
        # PyPlot.plot(ts[2:end],Fxh_direct,color=plot_cycle[1],".",label="Direct Library OWENS Side")        
    end
    
    # PyPlot.legend()
    # PyPlot.savefig("$(path)/$header_name.pdf",transparent = true)
end

# PyPlot.figure()
# PyPlot.plot(standalone_AD[2:end-1,1],(plotdata1_standalone[1,3:end].-fx_directcall[2:end-1])./plotdata1_standalone[1,3:end].*100,color=plot_cycle[1],"-")
# PyPlot.ylabel("Percent Difference (%)")
# PyPlot.xlabel("Time (s)")
# PyPlot.legend()