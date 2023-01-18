# Example usage the Julia wrapper for the HydroDyn DLL
import Libdl
import DelimitedFiles
import HDF5
path = splitdir(@__FILE__)[1]
import OpenFASTWrappers
using Test

cd(path)

md_lib_filename = "$path/../openfast/build/modules/moordyn/libmoordyn_c_binding" #change this to match your local path of the MoorDyn DLL
ptfm_motions_filename = "$path/OpenFAST_DisplacementTimeseries.dat"
# md_input_file = "$path/NRELOffshrBsline5MW_OC4DeepCwindSemi_MoorDynv2.dat"
md_input_file = "$path/data/moordyn_test.dat"
t_initial = 0.0
dt = 0.0125 #0.0125
t_max = 60.0  #60.0
interp_order = 2
num_corrections = 1

# load motions test file
ptfm_ts = DelimitedFiles.readdlm(ptfm_motions_filename, ',')
ptfm_pos_ts = ptfm_ts[:, 1:6]
ptfm_vel_ts = ptfm_ts[:, 7:12]
ptfm_acc_ts = ptfm_ts[:, 13:18]

# Preallocate matrices
ts = collect(t_initial:dt:t_max)
ptfm_force_ts = Array{Float64,2}(undef, length(ts), 6)
line_tensions_ts = Array{Float64,2}(undef, length(ts), 6)

forces = Vector{Float32}(undef, 6)
line_tensions = Vector{Float32}(undef, 6)

## Run MoorDyn
OpenFASTWrappers.MD_Init(md_lib_filename; init_ptfm_pos=ptfm_pos_ts[1,:], interp_order=interp_order)

# Time step zero
forces[:], line_tensions[:] = OpenFASTWrappers.MD_CalcOutput(t_initial, ptfm_pos_ts[1,:], ptfm_vel_ts[1,:], ptfm_acc_ts[1,:], forces, line_tensions)
ptfm_force_ts[1, :] = forces
line_tensions_ts[1, :] = line_tensions

# Time marching
for (idx, t) in enumerate(ts[1:end-2])
    for correction in range(1, num_corrections+1)
        # If there are correction steps, the inputs would be updated using outputs
        # from the other modules.
        ptfm_pos = ptfm_pos_ts[idx+1, :]
        ptfm_vel = ptfm_vel_ts[idx+1, :]
        ptfm_acc = ptfm_acc_ts[idx+1, :]

        OpenFASTWrappers.MD_UpdateStates(t, t+dt, ptfm_pos, ptfm_vel, ptfm_acc)

        OpenFASTWrappers.MD_CalcOutput(t+dt, ptfm_pos, ptfm_vel, ptfm_acc, forces, line_tensions)
    
        ptfm_force_ts[idx+1, :] = forces[:]
        line_tensions_ts[idx+1, :] = line_tensions[:]

        # When coupled to a different code, this is where the Force/Moment info
        # would be passed to the aerodynamic solver.

    end

end
OpenFASTWrappers.MD_End()

filename = "$path/data/moordyn_unit.h5"
# HDF5.h5open(filename, "w") do file
#     HDF5.write(file,"ts",ts) #power rating
#     HDF5.write(file,"ptfm_force_ts",ptfm_force_ts) #speed RPM
#     HDF5.write(file,"ptfm_pos_ts",ptfm_pos_ts) #Torque Nm
# end

# ts = HDF5.h5read(filename,"ts")
ptfm_force_ts_unit = HDF5.h5read(filename,"ptfm_force_ts")
ptfm_pos_ts_unit = HDF5.h5read(filename,"ptfm_pos_ts")

for itest = 1:length(ptfm_force_ts[:,1])-1
    for jtest = 1:length(ptfm_force_ts[1,:])-1
        @test isapprox(ptfm_force_ts_unit[itest,jtest],ptfm_force_ts[itest,jtest],atol=1e-6)
        @test isapprox(ptfm_pos_ts_unit[itest,jtest],ptfm_pos_ts[itest,jtest],atol=1e-6)
    end
end


# import PyPlot
# fig,ax = PyPlot.subplots()
# ax2 = ax.twinx()
# ax.plot(ts, ptfm_force_ts[:,1], color="red")
# ax.set_xlabel("Time [s]")
# ax.set_ylabel("Ptfm Force, Surge [N]", color="red")
#
# ax2.plot(ts, ptfm_pos_ts[:,1], color="blue")
# ax2.set_ylabel("Ptfm Displacement, Surge [m]", color="blue")
#
# PyPlot.show()
