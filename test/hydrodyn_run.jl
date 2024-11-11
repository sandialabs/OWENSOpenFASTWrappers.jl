# Example usage the Julia wrapper for the HydroDyn DLL
import DelimitedFiles
path = splitdir(@__FILE__)[1]
import OWENSOpenFASTWrappers
using Test

cd(path)

# hdlib_filename = "$path/../deps/bin/libhydrodyn_c_binding" #change this to match your local path of the HydroDyn DLL
# hdlib_filename = "$path/../deps/openfast/build/modules/hydrodyn/libhydrodyn_c_binding" #change this to match your local path of the HydroDyn DLL
ptfm_motions_filename = "$path/data/OpenFAST_DisplacementTimeseries.dat"
output_root_name = "$path/data/hd_wrapper_test"
PotFile = "$path/data/potential_flow_data/marin_semi"
hd_input_file = "$path/data/NRELOffshrBsline5MW_OC4DeepCwindSemi_HydroDyn.dat"
ss_input_file = "$path/data/NRELOffshrBsline5MW_OC4DeepCwindSemi_SeaState.dat"
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
OWENSOpenFASTWrappers.HD_Init(; output_root_name, hd_input_file, ss_input_file, PotFile, t_initial, dt, t_max)

# Time step zero
OWENSOpenFASTWrappers.HD_CalcOutput(t_initial, ptfm_pos_ts[1,:], ptfm_vel_ts[1,:], ptfm_acc_ts[1,:], forces, out_vals)

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

        OWENSOpenFASTWrappers.HD_UpdateStates(t, t+dt, ptfm_pos, ptfm_vel, ptfm_acc)

        # New values may be available at this point from the structural solver, so update them here
        # (this is redundant for uncoupled code, but I am keeping this in so it is not forgotten later)
        ptfm_pos = ptfm_pos_ts[idx+1, :]
        ptfm_vel = ptfm_vel_ts[idx+1, :]
        ptfm_acc = ptfm_acc_ts[idx+1, :]

        OWENSOpenFASTWrappers.HD_CalcOutput(t+dt, ptfm_pos, ptfm_vel, ptfm_acc, forces, out_vals)
        ptfm_force_ts[idx+1, :] = forces
        out_vals_ts[idx+1, 1] = out_vals[1]

        # When coupled to a different code, this is where the Force/Moment info
        # would be passed to the aerodynamic solver.

    end

end
OWENSOpenFASTWrappers.HD_End()

out_vals_ts_unit = [0.11191913; 0.099896014; 0.0878729; 0.07078146; 0.053690024; 0.032114245; 0.010538463; -0.014760992; -0.040060446; -0.06819408; -0.09632771;;]

for itest = 1:length(out_vals_ts_unit)
    @test isapprox(out_vals_ts_unit[itest],out_vals_ts[itest],atol = 1e-6)
end
