localpath = dirname(@__FILE__)
if Sys.isunix()

    println("Clone the correct openfast branch")
    cd(localpath)
    run(`rm -rf openfast/`)
    run(`git clone --depth 1 https://github.com/andrew-platt/openfast.git`)
    cd(joinpath(localpath, "openfast"))
    run(`git remote set-branches origin '*'`)
    run(`git fetch --depth 1 origin f/ADI_c_binding_multiRotor`)
    run(`git checkout f/ADI_c_binding_multiRotor`)
    
    println("Set up the build directory")
    run(`mkdir build`)
    cd(joinpath(localpath, "openfast/build"))
    
    println("run cmake")
    run(`cmake -DBUILD_SHARED_LIBS=ON -DOPENMP=ON ..`)
    # run(`cmake -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DBUILD_SHARED_LIBS=ON -DOPENMP=ON ..`)
    
    println("and compile")
    run(`make ifw_c_binding`)
    run(`make moordyn_c_binding`)
    run(`make hydrodyn_c_binding`)
    run(`make aerodyn_inflow_c_binding`)
    run(`make aerodyn_driver`)
    run(`make turbsim`)

else
    error("windows currently unsupported, please install openfast libraries manually")
end
