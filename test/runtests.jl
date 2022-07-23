import OpenFASTWrappers
using Test
path,_ = splitdir(@__FILE__)

@testset "OpenFASTWrappers.jl" begin
    # Usage Example
    # mkdir build
    # cd build
    # cmake ..
    # make ifw_c_binding
    # the the library is located in: deps/bin
    turbsim_filename = "$path/data/ifw/test.bts"
    OpenFASTWrappers.ifwinit(;turbsim_filename)
    velocity = OpenFASTWrappers.ifwcalcoutput([0.0,0.0,100.0],0.1)
    println(velocity)
    OpenFASTWrappers.ifwend()
    @test isapprox(velocity,[15.636677, -0.3423329, 0.24170564];atol=1e-4)

    ptfm_pos = [5.0, 0.0, 0.0, 0.0, 0.033161256, 0.0]
    ptfm_vel = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    ptfm_acc = [-0.024818709, -9.5636366e-05, 0.024033478, 3.4662767e-06, -0.0066289646, -1.1669358e-07]

    hd_forces = Vector{Float32}(undef, 6)
    out_vals = Vector{Float32}(undef, 1)
    hd_potdir = "$path/data/hd/potflow/marin_semi"
    OpenFASTWrappers.hdInit(PotFile=hd_potdir, t_initial=0.0, dt=0.0125, t_max=0.0125)
    hd_forces[:], _ = OpenFASTWrappers.hdCalcOutput(0.0, ptfm_pos, ptfm_vel, ptfm_acc, hd_forces, out_vals)
    println(hd_forces)
    OpenFASTWrappers.hdEnd()
    exp_hd_forces = [658294.5625000, 732.4172974, 45545132.0000000, -24004.8085938, -5331063.5000000, 788.6155396]
    @test isapprox(hd_forces,exp_hd_forces; rtol=1e-3)

    md_forces = Vector{Float32}(undef, 6)
    line_tensions = Vector{Float32}(undef, 6)
    OpenFASTWrappers.mdInit(init_ptfm_pos=ptfm_pos, interp_order=2)
    md_forces[:], _ = OpenFASTWrappers.mdCalcOutput(0.0, ptfm_pos, ptfm_vel, ptfm_acc, md_forces, line_tensions)
    println(md_forces)
    OpenFASTWrappers.mdEnd()
    exp_md_forces = [-378879.84, 0.35762787, -1.9042124e6, -14.305115, -2.4444795e6, 1.9222498]
    @test isapprox(md_forces,exp_md_forces; rtol=1e-3)
end