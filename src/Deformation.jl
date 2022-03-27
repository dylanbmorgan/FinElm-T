
#################################
#--------- Deformation ---------#
#################################

using LinearAlgebra


function planestrain(E, ν)
    constant = E / ((1 + ν) * (1 - (2*ν)))
    mat = [[1 - ν, ν, 0] [ν, 1 - ν, 0] [0, 0, 0.5*(1 - (2 * ν))]]
    C = constant * mat

    return C
end

function straindisp(rc, ξ, η)
    # Natural coors of quadrilateral
    # rc stands for real coors
    corners = zeros(4,2)
    for j = 1:2
        for i = 1:4
            corners[i,j] = rc[i][j]
        end
    end

    natcoors = [[-1, 1, 1, -1] [-1, -1, 1, 1]]

    # Derivatives of shape functions wrt natural coords
    dNdnat = zeros((2, 4))
    dNdnat[1,:] = 0.25 .* natcoors[:,1] .* (1 .+ natcoors[:,2] .* η)
    dNdnat[2,:] = 0.25 .* natcoors[:,2] .* (1 .+ natcoors[:,1] .* ξ)

    # Elemental Jacobian matrix
    Jmat = dNdnat * corners
    JmatInv = inv(Jmat)
    dNdx = sum(JmatInv, dims=2) .* dNdnat

    dsB = zeros((3, 8))
    dsB[1, 1:2:8] = dNdx[1,:]
    dsB[2, 2:2:8] = dNdx[2,:]
    dsB[3, 1:2:8] = dNdx[2,:]
    dsB[3, 2:2:8] = dNdx[1,:]

    return dsB, Jmat
end

function stiffmatrix(elem, C)
    # Stiffness element
    Ke = zeros((8,8))

    # Location of Gauss points
    a = 1/sqrt(3)

    # Weights of functions
    w = 1

    # Matrix of Gauss points
    gauss = [[-a, a, a, -a] [-a, -a, a, a]]

    for i = 1:4
        # Natural coors
        ξ = gauss[i,1]
        η = gauss[i,2]

        # B matrix
        dsB, J = straindisp(elem, ξ, η)

        # Jacobian determinant
        detJ = det(J)

        # Calculate K
        dsBT = dsB'
        dot1 = dsBT * C
        dot2 = dot1 * dsB
        Ke = Ke + dot2 * detJ * w * 10
    end

    return Ke
end

function forcevec(xyz, top, right, forcearea)
    f = zeros(2 * length(xyz))

    Threads.@threads for i = 1:length(xyz)
        # Extract y positions of top nodes
        ycoor = filter(isequal(top), xyz[i])
        append!(ycoor, 0)
        ynode = ycoor[1]

        f[i * 2] = ynode
    end

    # Find node indices to apply force to
    nnodes = findall(!isequal(0), f)

    # Total force on face and force per node
    tforce = forcearea * right
    nforce = tforce / length(nnodes)

    Threads.@threads for i = 1:length(nnodes)
        f[nnodes[i]] = nforce
    end

    return f
end

function reducemat!(xyz, k, f)
    ### Set boundary conditions to bottom of the shape ###
    bci1 = findall(coor -> 0 in coor[2], xyz)
    bci2 = bci1 .+ bci1[end]
    bci = vcat(bci1, bci2)

    ### Reduce force vector ###
    deleteat!(f, bci)

    ### Reduce glob stiff matrix ###
    # Find values in k that aren't in bci
    freenodes = filter(val -> !(val in bci), 1:size(k, 1))
    newsize = 1:size(k, 1) - length(bci)  # Size of new array
    newk = zeros(length(newsize), length(newsize))

    Threads.@threads for (i1, j1) = collect(zip(freenodes, newsize))
        for (i2, j2) = zip(freenodes, newsize)
            newk[j2, j1] = k[i2, i1]
        end
    end

    return f, newk, bci
end

function deform(xyz, D, scalefac)
    # Deform nodes
    # Convert xyz to match dims of D
    defxyz = zeros(length(xyz) * 2)
    Threads.@threads for (i, j) in collect(zip(1:2:length(xyz)*2, 1:length(xyz)))
        defxyz[i] = xyz[j][1]
        defxyz[i + 1] = xyz[j][2]
    end

    # Deform nodes including scaling factor
    scaleD = copy(D)
    scaleD .*= scalefac
    defxyz .+= scaleD

    return defxyz
end
