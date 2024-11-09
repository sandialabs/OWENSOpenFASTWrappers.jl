
"""
mymesh, myort, myjoint, AD15bldNdIdxRng, AD15bldElIdxRng = create_mesh_struts(;Htwr_base = 15.0,
Htwr_blds = 147.148-15.0,
    Hbld = 147.148-15.0, #blade height
    R = 54.014, # m bade radius
    AD15hubR = 2.0,
    nblade = 3,
    ntelem = 30, #tower elements
    nbelem = 30, #blade elements
    nselem = 4,
    strut_twr_mountpoint = [0.01,0.5,0.9],
    strut_bld_mountpoint = [0.01,0.5,0.9],
    bshapex = zeros(nbelem+1), #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    bshapez = zeros(nbelem+1),
    bshapey = zeros(nbelem+1), # but magnitude for this is relevant
    angularOffset = 0.0, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    AD15_ccw = false,
    verbosity=0, # 0 nothing, 1 basic, 2 lots: amount of printed information
    connectBldTips2Twr =true)

Standard Mesh Matching 5MW, 34m configurations

#Inputs
* `Htwr_base::float`: height of tower before blades attach (m)
* `Htwr_blds::float`: height of the tower once the blades attach (m)
* `Hbld::float`: blade height (m)
* `R::float`: bade radius (m)
* `nblade::int`: number of blades
* `ntelem::int`: number of tower elements
* `nbelem::int`: number of blade elements
* `nselem::int`: number of strut elements
* `strut_twr_mountpoint::float` = [0.01,0.5,0.9], # factor of blade height where the bottom strut attaches on the tower # This puts struts at top and bottom, as a fraction of the blade position
* `strut_bld_mountpoint::float` = [0.01,0.5,0.9], # factor of blade height where the bottom strut attaches on the blade # This puts struts at bottom 0, mid 0.5, and top 1.0 as a fraction of the blade position
* `bshapex::Array{<:float}`: Blade shape, magnitude is irrelevant, scaled based on height and radius above
* `bshapez::Array{<:float}`: Blade shape, magnitude is irrelevant, scaled based on height and radius above
* `bshapey::Array{<:float}`: Blade shape, magnitude IS relevant #TODO: resolve this
* `angularOffset::float`: (rad) angular offset of mesh generation, typically used to match CACTUS input.  Value of 0 puts blade 1 at the "north" position and the others populate counterclockwise when looking down
* `AD15_ccw::boolean`: Use AD15 convention of VAWT counter-clockwise with blade root at top (blade points down)
* `AD15hubR::float`: AD15 has a hub radius, so the struts do not go all the way to the center of the axis of rotation, while the structural mesh does.
# `verbosity::int`: 0 nothing, 1 basic, 2 lots: amount of printed information 
* `connectBldTips2Twr::book`: True for Darrieus style, false for H-VAWT, but the blade shapes should be appropriate

#Outputs
* `mymesh::OWENSFEA.Mesh`: see ?OWENSFEA.Mesh
* `myort::OWENSFEA.Ort`: see ?OWENSFEA.Ort
* `myjoint:Array{<:float}`: see ?OWENSFEA.FEAModel.joint
* `AD15bldNdIdxRng`: indices for start and end of all blades for AD15 (includes struts).  Note that strut start nodes may be inside the strut (strut connects to tower, AD15 blade connects to hub wich is a few nodes away from tower)
* `AD15bldElIdxRng`: range of elements for start and end of all AD15 blades (includes struts)
"""
function create_mesh_struts(;Htwr_base = 15.0,
    Htwr_blds = 147.148-15.0,
    Hbld = 147.148-15.0, #blade height
    R = 54.014, # m bade radius
    AD15hubR = 2.0,
    nblade = 3,
    ntelem = 30, #tower elements
    nbelem = 30, #blade elements
    nselem = 4,
    strut_twr_mountpoint = [0.125,0.5,0.95], # This puts struts at top and bottom, as a fraction of the blade position
    strut_bld_mountpoint = [0.25,0.5,0.75], # This puts struts at bottom 0, mid 0.5, and top 1.0 as a fraction of the blade position
    bshapex = zeros(nbelem+1), #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    bshapez = zeros(nbelem+1),
    bshapey = zeros(nbelem+1), # but magnitude for this is relevant
    angularOffset = 0.0, #Blade shape, magnitude is irrelevant, scaled based on height and radius above
    AD15_ccw = false,
    verbosity = 0, # 0 nothing, 1 basic, 2 lots: amount of printed information
    connectBldTips2Twr = true)

    Nstrut = length(strut_bld_mountpoint)

    if length(strut_bld_mountpoint) != length(strut_twr_mountpoint)
        error("strut_twr_mountpoint must be the same length as strut_bld_mountpoint")
    end

    ##################################
    #             _
    #           /_|_\
    #          |  |  )
    #           \-|-/
    ####################################
    ###------------Tower--------------##
    ####################################
    mesh_z = collect(LinRange(0,Htwr_blds+Htwr_base,ntelem+1))

    # Insert mount point base
    mesh_z = sort([mesh_z;Htwr_base])
    t_botidx = findall(x->isapprox(x,Htwr_base,atol=1e-5),mesh_z)[1]#:nblade]

    t2s_idx = zeros(Int, Nstrut)
    # Insert strut mount points
    for istrut = 1:Nstrut
        if maximum(Hbld*strut_twr_mountpoint[istrut]+Htwr_base .== mesh_z) # if we are at exactly an existing node, then offset our mount point
            @warn "Mesh warning: two points are directly on top of one another, consider adjusting number of elements to space out the mesh"
            strut_twr_mountpoint[istrut] += 1e-6
        end
        mesh_z = sort([mesh_z;Hbld*strut_twr_mountpoint[istrut]+Htwr_base])

        # pick out the strut mounting indices
        t2s_idx[istrut] = findall(x->isapprox(x,Hbld*strut_twr_mountpoint[istrut]+Htwr_base,atol=1e-5*Hbld),mesh_z)[1]
    end

    # Create the x and y components of same size as mesh_z now that the strut mount points are inserted
    mesh_x = zero(mesh_z)
    mesh_y = zero(mesh_z)

    # pick out the tower top index
    t_topidx = length(mesh_z)

    # intra-tower connectivity
    conn = zeros(length(mesh_z)-1,2)
    conn[:,1] = collect(1:length(mesh_z)-1)
    conn[:,2] = collect(2:length(mesh_z))

    #####################################
    ###------------Blades--------------##
    #####################################

    #connection points on tower are simply the bottom of the tower offset connecting to the bottom of the blades, which is jointed if Darrieus in the joint matrix below
    bld_Z = collect(LinRange(0.0,Hbld,nbelem+1))

    # Insert bottom strut mount point
    for istrut = 1:Nstrut
        if maximum(Hbld*strut_bld_mountpoint[istrut] .== bld_Z) # if we are at exactly an existing node, then offset our mount point
            @warn "Mesh warning: two points are directly on top of one another, consider adjusting number of elements to space out the mesh"
            strut_bld_mountpoint[istrut] -= 1e-6
        end
        bld_Z = sort([bld_Z;Hbld*strut_bld_mountpoint[istrut]])
    end

    if bshapex == zeros(nbelem+1)
        bld_Y = R.*(1.0.-4.0.*(bld_Z/Hbld.-.5).^2)
    else
        # Ensure the blade shape conforms to the turbine height and radius specs
        bshapex = R .* bshapex./maximum(bshapex)
        bshapez = Hbld .* bshapez./maximum(bshapez)
        bld_Y = safeakima(bshapez,bshapex,bld_Z)
    end

    if bshapey == zeros(nbelem+1)
        bld_X = zero(bld_Y)
    else
        bld_X = safeakima(bshapez,bshapey,bld_Z)
    end

    # AeroDyn Compatability
    AD15bldNdIdxRng = zeros(Int64,0,2)

    bld_Z .+= Htwr_base

    b_Z = []
    b_X = []
    b_Y = []
    # Now using standard VAWT convention, blade 1 is zero degrees at top dead center, or North/Y+
    # and they are offset counter clockwise
    b_topidx = zeros(Int,nblade)
    b_botidx = zeros(Int,nblade) .+ length(mesh_z)
    conn_b = zeros(length(bld_Z)-1,2)
    for ibld = 1:nblade
        myangle = (ibld-1)*2.0*pi/nblade + angularOffset
        b_Z = [b_Z;bld_Z]
        b_X = [b_X;bld_Y.*sin(myangle).+bld_X.*cos(myangle)]
        b_Y = [b_Y;-bld_X.*sin(myangle).+bld_Y.*cos(myangle)]

        # Element joint indices
        b_botidx[ibld] = length(mesh_z)+1 + length(bld_Z)*(ibld-1)
        b_topidx[ibld] = length(mesh_z)+1 + length(bld_Z)*ibld-1

        # Intraconnectivity
        conn_b[:,1] = collect(b_botidx[ibld]:1:b_topidx[ibld]-1)
        conn_b[:,2] = collect(b_botidx[ibld]+1:1:b_topidx[ibld])
        conn = [conn;conn_b]

        if AD15_ccw #Clockwise, the blades roots are at the top, trailing edge is always positive y
            AD15bldNdIdxRng = [AD15bldNdIdxRng; b_topidx[ibld] b_botidx[ibld]]    # top of blade is root 
        elseif !(AD15_ccw) #Clockwise, the blades roots are at the bottom
            AD15bldNdIdxRng = [AD15bldNdIdxRng; b_botidx[ibld] b_topidx[ibld]]    # bottom of blade is root
        end
    end

    # pick out the strut mounting indices
    b2s_idx = zeros(Int,nblade,Nstrut)
    for istrut = 1:Nstrut
        b2s_idx[:,istrut] = findall(x->x==Hbld*strut_bld_mountpoint[istrut]+Htwr_base,b_Z)[1:nblade] .+ length(mesh_z)
    end

    # Add to the mesh
    mesh_z = [mesh_z;b_Z]
    mesh_x = [mesh_x;b_X]
    mesh_y = [mesh_y;b_Y]

    #####################################
    ###------------Struts--------------##
    #####################################

    function createstrut(sstartidx,sendidx,mesh_x,mesh_y,mesh_z,conn,AD15hubR)

        sxstart = mesh_x[sstartidx]
        systart = mesh_y[sstartidx]
        szstart = mesh_z[sstartidx]

        sxend = mesh_x[sendidx]
        syend = mesh_y[sendidx]
        szend = mesh_z[sendidx]

        # Now draw the lines
        s_x = collect(LinRange(sxstart,sxend,nselem+1))
        s_y = collect(LinRange(systart,syend,nselem+1))
        s_z = collect(LinRange(szstart,szend,nselem+1))

        hubIdx = 1
        if AD15hubR > 1e-6
            lenXY = sqrt((sxend - sxstart)^2 + (syend - systart)^2)   # strut length in XY
            minR2 = lenXY 
            for i = 1:nselem+1  # step through to find closest point to hub radius on x-y plane
                R2 = AD15hubR - sqrt((s_x[i] - sxstart)^2 + (s_y[i] - systart)^2)
                if abs(R2) < abs(minR2)
                    hubIdx = i
                    minR2 = R2
                end
            end
            R_temp = minR2

            s_x[hubIdx] = s_x[hubIdx] + R_temp/lenXY*(sxend-sxstart)
            s_y[hubIdx] = s_y[hubIdx] + R_temp/lenXY*(syend-systart)
            s_z[hubIdx] = s_z[hubIdx] + R_temp/lenXY*(szend-szstart)

            if verbosity>0
                println("Hub crossing at idx $hubIdx and radially at $R_temp with AD15 hub radius of $AD15hubR")
                println("Moving strut point from [$(s_x[hubIdx]),$(s_y[hubIdx]),$(s_z[hubIdx])] to [$(s_x[hubIdx]),$(s_y[hubIdx]),$(s_z[hubIdx])]")
            end

        end

        hubIdx += length(mesh_z)

        # joint connections
        s2b_idx_internal = length(mesh_z)+1
        s2t_idx_internal = s2b_idx_internal+length(s_z)-1

        # and add to the mesh
        mesh_x = [mesh_x;s_x]
        mesh_y = [mesh_y;s_y]
        mesh_z = [mesh_z;s_z]

        # Intraconnectivity
        conn_s = zeros(nselem,2)
        conn_s[:,1] = collect(s2b_idx_internal:1:s2t_idx_internal-1)
        conn_s[:,2] = collect(s2b_idx_internal+1:1:s2t_idx_internal)
        conn = [conn;conn_s]

        return s2b_idx_internal,s2t_idx_internal,mesh_x,mesh_y,mesh_z,conn,hubIdx
    end

    #Connect from the tower to the blades
    # For each blade, find the mounting location and draw a line
    s2b_idx = zeros(Int,nblade,Nstrut)
    s2t_idx = zeros(Int,nblade,Nstrut)

    # Bottom Struts
    for istrut = 1:Nstrut
        for ibld = 1:nblade
            s2t_idx[ibld,istrut],s2b_idx[ibld,istrut],mesh_x,mesh_y,mesh_z,conn,hubIsectIdx = createstrut(t2s_idx[istrut],b2s_idx[ibld,istrut],mesh_x,mesh_y,mesh_z,conn,AD15hubR)
            
            AD15bldNdIdxRng = [AD15bldNdIdxRng; hubIsectIdx  s2b_idx[ibld,istrut]] #AD15 struts always start at hub regardless of rotation, but watch out for airfoil orientation!
        end
    end

    #######################################
    ###   Cleanup/Derived parameters    ###
    #######################################

    numNodes = length(mesh_z)
    nodeNum = collect(LinRange(1,numNodes,numNodes))
    numEl = length(conn[:,1])
    elNum = collect(LinRange(1,numEl,numEl))

    # Define Mesh Types
    # Mesh Type: 0-blade 1-tower, treat struts like blades
    meshtype = zeros(Int,numEl)

    # Find elnum associated with t_topidx
    topel_idx = findall(x->x==t_topidx,conn[:,2])
    meshtype[1:topel_idx[1]] .= 1 #Tower

    #########################
    # .bld equivalent
    #########################

    meshSeg = zeros(Int,1+nblade+Nstrut*nblade) #tower, blades, and struts

    meshSeg[1] = topel_idx[1]
    meshSeg[2:nblade+1] .= nbelem+Nstrut #+Nstrut for strut mount points
    meshSeg[nblade+2:end] .= nselem

    # For each blade
    structuralSpanLocNorm = zeros(nblade,length(bld_Z))
    structuralNodeNumbers = zeros(nblade,length(bld_Z))
    structuralElNumbers = zeros(nblade,length(bld_Z))

    for iblade = 1:nblade

        # Normalized Span
        span_len = bld_Z.-Htwr_base#sqrt.(bld_X.^2.0.+bld_Y.^2.0.+(bld_Z.-Htwr_base).^2.0)
        structuralSpanLocNorm[iblade,:] = span_len./maximum(span_len)

        # Node Numbers
        structuralNodeNumbers[iblade,:] = collect(b_botidx[iblade]:b_topidx[iblade])

        # Element Numbers
        structuralElNumbers[iblade,:] = structuralNodeNumbers[iblade,:].-iblade
        structuralElNumbers[iblade,end] = -1 #TODO: figure out why this is in the original OWENS setup and if it is used
    end

    mymesh = Mesh(nodeNum,numEl,numNodes,mesh_x,mesh_y,mesh_z,elNum,Int.(conn),meshtype,meshSeg,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers)

    ######################################
    ####----------Joint Matrix----------##
    ######################################

    # Connect Tower Top to Blades bottom, then each cable to each blade bottom
    # Then each cable to each blade top, then the latter two blade tops to the first

    njoint = nblade*2+Nstrut*nblade*2 #create the full array and then pull out zeros below if H-VAWT where the blades aren't connected to the tower.
    jointconn = zeros(Int,njoint,2)
    jointtype = zeros(njoint)
    for ibld = 1:nblade
        if connectBldTips2Twr
            # connect tower to blades
            jointconn[ibld,:] = [t_botidx b_botidx[ibld]]
        end

        for istrut = 1:Nstrut
            # connect tower to each strut
            jointconn[ibld+nblade*istrut,:] = [t2s_idx[istrut] s2t_idx[ibld,istrut]]
        end

        if connectBldTips2Twr
            # connect tower to blades tops
            jointconn[ibld+nblade*(Nstrut+1),:] = [t_topidx b_topidx[ibld]]
        end

        for istrut = 1:Nstrut
            # connect strut to blade bottom
            jointconn[ibld+nblade*(Nstrut+1+istrut),:] = [s2b_idx[ibld,istrut] b2s_idx[ibld,istrut]]
        end
    end

    # Reduce the matrix based on if the blades got connected or not, throwing out all the zero rows
    bitlogic = jointconn[:,1] .!= 0.0
    jointconn = jointconn[bitlogic,:]
    jointtype = jointtype[bitlogic]

   

    njoint = length(jointconn[:,1]) # reset the number of joints
    myort = calculateElementOrientation(mymesh)
    # println("start")
    Psi_d_joint = zeros(njoint)
    Theta_d_joint = zeros(njoint)
    for jnt = 1:njoint
        elnum_of_joint = findall(x->x==jointconn[jnt,1],conn[:,1]) #gives index of the elNum vector which contains the point index we're after. (the elNum vector is a map between point index and element index)
        if length(elnum_of_joint)==0 #Use the other element associated with the joint
            elnum_of_joint = findall(x->x==jointconn[jnt,1],conn[:,2])
        end
        if length(elnum_of_joint)==0 #Use the other element associated with the joint
            elnum_of_joint = findall(x->x==jointconn[jnt,2],conn[:,2])
        end
        Psi_d_joint[jnt] = myort.Psi_d[elnum_of_joint[1]]
        Theta_d_joint[jnt] = myort.Theta_d[elnum_of_joint[1]]
    end
    #Joint Types: (0 = weld(fixed), 1=pinned, 2 = hinge joint with axis about slave node element’s e2 axis, 3 = hinge joint axis about slave node element’s e1 axis, 4 = hinge joint axis about slave node element’s e3 axis)

    #Joint Number,   Joint Connections, Joint Type, Joint Mass, Not Used, Psi_D, Theta_D
    myjoint = [Float64.(1:1:njoint) jointconn jointtype zeros(njoint) zeros(njoint) Psi_d_joint Theta_d_joint]

    # Blade and strut starting and ending node and element numbers
    AD15bldElIdxRng = zeros(Int64,0,2)
    for i = 1:size(AD15bldNdIdxRng,1)
        if AD15bldNdIdxRng[i,2] > AD15bldNdIdxRng[i,1]  # ascending order
            idx1 = findfirst(x->x==AD15bldNdIdxRng[i,1], mymesh.conn[:,1])
            idx2 = findfirst(x->x==AD15bldNdIdxRng[i,2], mymesh.conn[:,2])
        else    # upside down oriented blade
            idx1 = findlast(x->x==AD15bldNdIdxRng[i,1], mymesh.conn[:,2])
            idx2 = findlast(x->x==AD15bldNdIdxRng[i,2], mymesh.conn[:,1])
        end

        if isnothing(idx2)
            idx2 = findlast(x->x==AD15bldNdIdxRng[i,2], mymesh.conn[:,2])
        end
        AD15bldElIdxRng = [AD15bldElIdxRng; idx1 idx2]
    end
    return mymesh, myort, myjoint, AD15bldNdIdxRng, AD15bldElIdxRng
end

"""

    calculateElementOrientation(mesh)

Calculates the orientation of elements in a mesh.

#Input
* `mesh::OWENSFEA.Mesh`  see ?OWENSFEA.Mesh object containing mesh data

#Output
* `elOr::OWENSFEA.Ort`  see ?OWENSFEA.Ort object containing element orientation data
"""
function calculateElementOrientation(mesh)

    # Note on gimbal lock:
    #   when calculating a (roll -> pitch -> yaw) rotation sequence (twist -> theta -> psi) on a vertical element, it is ambiguous
    #   if the roll or yaw (twist or psi) should be set.  It is therefore up to the designer to pick which is used.  For gimbal lock
    #   the yaw (psi) is currently defined as zero by choice and the roll (twist) set to non-zero.  It could be equivalently
    #   calculated as the roll (twist) is zero, the yaw (psi) is non-zero.  The latter scenario may actually be simpler for
    #   calculating DCM's when coupling to other codes.

    numEl = mesh.numEl #get number of elements
    Psi_d=zeros(numEl) #initialize Psi, Theta, Twist, and Offset Arrays
    Theta_d=zeros(numEl)
    twist_d=zeros(numEl)
    Offset=zeros(3,numEl)    #offset is the hub frame coordinate of node 1 of the element
    elNum=zeros(numEl,2) #initialize element number array

    lenv = zeros(numEl)
    for i = 1:numEl #loop over elements

        n1 = Int(mesh.conn[i,1]) #n1 := node number for node 1 of element i
        n2 = Int(mesh.conn[i,2]) #n2 := node number for node 2 of element i

        p1 = [mesh.x[n1] mesh.y[n1] mesh.z[n1]] #nodal coordinates of n1
        p2 = [mesh.x[n2] mesh.y[n2] mesh.z[n2]] #nodal coordinates of n2
        Offset[:,i] = p1 #set offset as position of n1
        elpos = (p1.+p2)./2

        Psi_d[i] = -atand(elpos[1],elpos[2]).-90.0 #global yaw position, this calculation is agnostic to vertical elements, keep in mind that top dead center is 0 degrees yaw

        if mesh.type[i] == 4 # treat tangentially aligned mesh components different than radially aligned
            Psi_d[i] -= 90.0
            twist_d[i] += 90.0
        end

        # Now with the global yaw position know, get the node points in a consistent frame of reference to calculate the delta, or the slope of the element
        p1[1],p1[2],p1[3] = rigidBodyRotation(p1[1],p1[2],p1[3],[-Psi_d[i]],[3])
        p2[1],p2[2],p2[3] = rigidBodyRotation(p2[1],p2[2],p2[3],[-Psi_d[i]],[3])
        
        v=p2-p1
        v[abs.(v).<1e-7] .= 0.0 #zero out close to zero differences

        Theta_d[i]  = atand(v[1],v[3]).-90.0

        lenv[i] = LinearAlgebra.norm(v) #calculate element length

        elNum[i,:] = mesh.conn[i,:] #get node number map

    end

    #assign data to element orientation (Ort) object
    # yaw, blade slope, angle of attack
    return Ort(Psi_d,Theta_d,twist_d,lenv,elNum,Offset)
end


"""
    Ort(Psi_d,Theta_d,Twist_d,Length,elNum,Offset)

Struct with element orientation

# Inputs
* `Psi_d::Array{<:float}`: length NumEl, element rotation about 3 in global FOR (deg) These angles are used to transform from the global coordinate frame to the local element/joint frame via a 3-2 Euler rotation sequence.
* `Theta_d::Array{<:float}`: length NumEl, element rotation about 2 (deg)
* `Twist_d::Array{<:float}`: length NumEl, element twist (deg)
* `Length::Array{<:float}`: length NumEl, element length (m)
* `elNum::Array{<:float}`: Element number the other arrays are associated with
* `Offset::Array{<:float}`: hub frame coordinate of node 1 of the element

# Outputs:
* `none`:

"""
mutable struct Ort
    Psi_d
    Theta_d
    Twist_d
    Length
    elNum
    Offset
end


"""
    Mesh(nodeNum,numEl,numNodes,x,y,z,elNum,conn,type,meshSeg,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers)

Struct with mesh definition

# Inputs
* `nodeNum::Array{<:int}`: Number mapping of nodes (typically 1:Nnodes)
* `numEl::int`: total number of elements
* `numNodes::int`: total number of nodes
* `x::Array{<:float}`: Nodal x position
* `y::Array{<:float}`: Nodal y position
* `z::Array{<:float}`: Nodal z position
* `elNum::Array{<:int}`: Number mapping of elements (typically 1:Nelements)
* `conn::Array{<:int}`: Nelemx2 connectivity between nodes, gaps between joints (which are defined in the joints)
* `type::Array{<:int}`: 0-blade 1-tower 2-strut
* `meshSeg::Array{<:int}`: number of nodes within each segment, with segments consisting of tower, blade 1 2 etc, struts
* `structuralSpanLocNorm::Array{<:float}`: Should be named heigh loc norm - unitized position along the blade height, used for aeroload mapping
* `structuralNodeNumbers::Array{<:int}`: Node numbers associated with blades for aero loads mapping
* `structuralElNumbers::Array{<:int}`: Element numbers associated with blades for aero loads mapping
* `nonRotating::Array{<:int}`: size(Nsections,numNodes) where nsections are the number of sections of the mesh that are non-rotating, like if for some reason you had two towers, or if you had multiple guy wires
* `hubNodeNum::int`: Node number where the rotating part of the turbine starts, assumes meshing always starts with tower, then blades, etc.

# Outputs:
* `none`:

"""
mutable struct Mesh
    nodeNum
    numEl
    numNodes
    x
    y
    z
    elNum
    conn
    type
    meshSeg
    structuralSpanLocNorm
    structuralNodeNumbers
    structuralElNumbers
    nonRotating
    hubNodeNum
    hubPos
    hubAngle
end

Mesh(nodeNum,numEl,numNodes,x,y,z,elNum,conn,type,meshSeg,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers) = Mesh(nodeNum,numEl,numNodes,x,y,z,elNum,conn,type,meshSeg,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers,0,1,zeros(3),zeros(3))

function safeakima(x,y,xpt)
    if minimum(xpt)<(minimum(x)-abs(minimum(x))*0.1) || maximum(xpt)>(maximum(x)+abs(maximum(x))*0.1)
        msg="Extrapolating on akima spline results in undefined solutions minimum(xpt)<minimum(x) $(minimum(xpt))<$(minimum(x)) or maximum(xpt)<maximum(x) $(maximum(xpt))>$(maximum(x))"
        throw(OverflowError(msg))
    end
    return FLOWMath.akima(x,y,xpt)
end


"""

    rigidBodyRotation(B1,B2,B3,AngleArray,AxisArray)

Performs a coordinate transformation from a local
body "B"(element) frame to a common hub frame "H" via a 3-2-3 euler
rotation sequence

#Input
* `B1`:         array containing body frame 1 coordinates of points to be mapped to the hub frame
* `B2`:         array containing body frame 2 coordinates of points to be mapped to the hub frame
* `B3`:         array containing body frame 3 coordinates of points to be mapped to the hub frame
* `AngleArray`:  Array of angles for Euler rotation sequence
* `AxisArray`:   Array of axes for Euler rotation sequence

#Output
* `H1`:         array containg hub frame 1 coordinates of points mapped to the hub frame from body frame
* `H2`:         array containg hub frame 2 coordinates of points mapped to the hub frame from body frame
* `H3`:         array containg hub frame 3 coordinates of points mapped to the hub frame from body frame

That is CHtoB = [M3(SweepAngle)][M2(Theta)][M3(Psi)];
"""
function rigidBodyRotation(B1,B2,B3,AngleArray,AxisArray)

    #calculate coordinate transformation matrix from element frame to
    #hub frame (CBtoH)
    dcm = createGeneralTransformationMatrix(AngleArray,AxisArray)
    C = dcm'

    #transform body coordinatized vector to be coordinatized in the hub
    #frame
    H1 = C[1,1].*B1 + C[1,2].* B2 + C[1,3].*B3
    H2 = C[2,1].*B1 + C[2,2].* B2 + C[2,3].*B3
    H3 = C[3,1].*B1 + C[3,2].* B2 + C[3,3].*B3

    return H1,H2,H3
end


"""

    createGeneralTransformationMatrix(angleArray,axisArray)

Calculates the transformation matrix assocaited with a general Euler rotation sequence.

#Input
* `angleArray`:      = array of angles for Euler rotation sequence
* `axisArray`:       = array of axis of rotatoins for Euler rotation

#Output
* `dcmTotal`:        = transformation matrix of specified euler rotation sequence
"""
function createGeneralTransformationMatrix(angleArray,axisArray)

    numRotations = length(angleArray) #get number of rotation to perform
    dcmArray = zeros(3,3,numRotations) #initialize individual rotation direction cosine matrix arrays

    for i=1:numRotations #calculate individual single rotatio direction cosine matrices
        dcmArray[:,:,i] = createSingleRotationDCM(angleArray[i],axisArray[i])
    end

    dcmTotal = dcmArray[:,:,1] #initialize dcmTotal as first rotation

    #multiply consecutive rotation sequency direction cosine matrices to arrive at overall transformation matrix
    for i=2:1:numRotations
        dcmTotal = dcmArray[:,:,i]*dcmTotal
    end

    return dcmTotal
end


"""
Creates a direction cosine matrix (dcm) associated with a rotation of angleDeg about axisNum.
"""
function createSingleRotationDCM(angleDeg,axisNum)

    angleRad = angleDeg*pi/180.0 #convert angle to radians

    if axisNum == 1 #direction cosine matrix about 1 axis
        dcm = [1.0 0.0 0.0
        0.0 cos(angleRad) sin(angleRad)
        0.0 -sin(angleRad) cos(angleRad)]
    elseif axisNum == 2 #direction cosine matrix about 2 axis
        dcm = [cos(angleRad) 0.0 -sin(angleRad)
        0.0 1.0 0.0
        sin(angleRad) 0.0 cos(angleRad)]
    elseif axisNum == 3 #direction cosine matrix about 3 axis
        dcm = [cos(angleRad) sin(angleRad) 0.0
        -sin(angleRad) cos(angleRad) 0.0
        0.0 0.0 1.0]
    else  #error catch
        error("Error: createSingleRotationDCM. Axis number must be 1, 2, or 3.")
    end

    return dcm

end