# Example usage the Julia wrapper for the HydroDyn DLL
import Libdl
import DelimitedFiles
path = splitdir(@__FILE__)[1]
import OpenFASTWrappers
using Test

cd(path)

# hd_lib_filename = "$path/../deps/bin/libhydrodyn_c_binding" #change this to match your local path of the HydroDyn DLL
hd_lib_filename = "/builds/8921-VAWT-TOOLS/OpenFASTWrappers.jl/openfast/build/modules/hydrodyn/libhydrodyn_c_binding" #change this to match your local path of the HydroDyn DLL
ptfm_motions_filename = "$path/OpenFAST_DisplacementTimeseries.dat"
output_root_name = "$path/hd_wrapper_test"
potmod_dir = "$path/data/potential_flow_data/marin_semi"
hd_input_file = "$path/NRELOffshrBsline5MW_OC4DeepCwindSemi_HydroDyn.dat"
num_corrections = 0

t_initial = 0.0
dt = 0.1
t_max = 1.0  #600.0

# load motions test file
ptfm_ts = DelimitedFiles.readdlm(ptfm_motions_filename, ',')
ptfm_pos_ts = ptfm_ts[:, 1:6]
ptfm_vel_ts = ptfm_ts[:, 7:12]
ptfm_acc_ts = ptfm_ts[:, 13:18]

# Preallocate matrices
ts = collect(t_initial:dt:t_max)
ptfm_force_ts = zeros(Cfloat,length(ts), 6)
out_vals_ts = zeros(Cfloat,length(ts), 1)
forces = zeros(Cfloat,6)
out_vals = zeros(Cfloat,43)
# Run HydroDyn
OpenFASTWrappers.HD_Init(hd_lib_filename, output_root_name, hd_input_file=hd_input_file, PotFile=potmod_dir, t_initial=t_initial, dt=dt, t_max=t_max)

# Time step zero
OpenFASTWrappers.HD_CalcOutput(t_initial, ptfm_pos_ts[1,:], ptfm_vel_ts[1,:], ptfm_acc_ts[1,:], forces, out_vals)

ptfm_force_ts[1, :] = forces
out_vals_ts[1, 1] = out_vals[1]

# Time marching
for (idx, t) in enumerate(ts[1:end-1])
    for correction = 1:num_corrections+1

        # If there are correction steps, the inputs would be updated using outputs
        # from the other modules.
        ptfm_pos = ptfm_pos_ts[idx+1, :]
        ptfm_vel = ptfm_vel_ts[idx+1, :]
        ptfm_acc = ptfm_acc_ts[idx+1, :]

        OpenFASTWrappers.HD_UpdateStates(t, t+dt, ptfm_pos, ptfm_vel, ptfm_acc)

        # New values may be available at this point from the structural solver, so update them here
        # (this is redundant for uncoupled code, but I am keeping this in so it is not forgotten later)
        ptfm_pos = ptfm_pos_ts[idx+1, :]
        ptfm_vel = ptfm_vel_ts[idx+1, :]
        ptfm_acc = ptfm_acc_ts[idx+1, :]

        OpenFASTWrappers.HD_CalcOutput(t+dt, ptfm_pos, ptfm_vel, ptfm_acc, forces, out_vals)
        ptfm_force_ts[idx+1, :] = forces
        out_vals_ts[idx+1, 1] = out_vals[1]

        # When coupled to a different code, this is where the Force/Moment info
        # would be passed to the aerodynamic solver.

    end

end
OpenFASTWrappers.HD_End()

out_vals_ts_unit = [-0.9369897842407227, -0.8816386461257935, -0.826287567615509, -0.7560029625892639, -0.6857184171676636, -0.6023154258728027, -0.5189124345779419, -0.4251416325569153, -0.3313708007335663, -0.2307586669921875, -0.13014645874500275]

for itest = 1:length(out_vals_ts_unit)
    @test isapprox(out_vals_ts_unit[itest],out_vals_ts[itest],atol = 1e-6)
end
