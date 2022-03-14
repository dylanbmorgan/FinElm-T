#!/usr/bin/env julia

# Beginning of meshgrid module

"""
    pointsdiff(totpoints, xstart, eright, ystart, etop)

Calculate the number of meshpoints required for each side of the shape
"""
function pointsdiff(totpoints, xstart, eright, ystart, etop)
    # Generate a factor for difference of shape length
    xdir = abs(eright - xstart)
    ydir = abs(etop - ystart)
    factor1 = xdir / ydir
    factor2 = ydir - 1 / xdir

    # Work out number of grid points in x and y directions
    if factor1 > 1
        factor = factor1
        ypoints = round(sqrt(totpoints / factor))
        xpoints = round(factor * ypoints)
    else
        factor = factor2
        xpoints = round(sqrt(totpoints / factor))
    end

    if iseven(xpoints)
        xpoints += 1
    end

    return xpoints
end

"""
    elements(xpoints, ypoints, xstart, ystart, etop, eright)

Create arrays of meshpoints for each side of the shape
"""
function mesharrays(xpoints, ypoints, xstart, ystart, etop, eright)
    # Generate points in x-direction
    bot = zeros(xpoints)
    Threads.@threads for (i,num) in collect(enumerate(LinRange(xstart, eright, xpoints)))
        bot[i] = num
    end

    # Generate points in y-direction
    left = zeros(ypoints)
    Threads.@threads for (j,num) in collect(enumerate(LinRange(ystart, etop, ypoints)))
        left[j] = num
    end

    right = copy(left)
    top = copy(bot)

    return bot, top, left, right
end

"""
    mesh(bot, top, left, right)

Generate a finite element mesh grid for some 2D object
"""
function mesh(bot, top, left, right)
    # Number of nodes and elements in the domain 
    nnodesx = size(bot, 1)
    nnodesy = size(left, 1)
    nelx = nnodesx - 1
    nely = nnodesy - 1
    nnodes = nnodesx * nnodesy

    # Dimensions of the domain 
    lx = bot[nnodesx] - bot[1]
    ly = left[nnodesy] - left[1]

    # Generate coords of nodes 'xyz' with array of arrays
    xyz = fill(Float64[], nnodes)

    Threads.@threads for i in 1:nnodesy
        yl = left[i] - left[1]
        dy = right[i] - left[i]

        for j in 1:nnodesx
            xb = bot[j] - bot[1]
            dx = top[j] - bot[j]

            xcoor = (dx*yl + xb*ly) / (ly - dx*dy / lx)
            ycoor = dy / lx*xcoor + yl

            xyz[j + (i-1)*nnodesx] = Float64[xcoor + bot[1], ycoor + left[1]]
        end
    end

    # Elements in the solid
    nel = nelx * nely
    elem = fill(Vector[], nel)

    Threads.@threads for i in 1:nely
        for j in 1:nelx
            elem[j + (i-1)*nelx] = [xyz[j + (i-1)*nnodesx],
                                    xyz[j + (i-1)*nnodesx + 1],
                                    xyz[j + i*nnodesx + 1],
                                    xyz[j + i*nnodesx]]
        end
    end

    return xyz, elem
end

function condof(elements, xyz)
    # Obtain connective indices between nodes to form elements
    con = fill(Int[], length(elements))

    Threads.@threads for (i, element) in collect(enumerate(elements))
        con[i] = [j for j in findall(pos -> pos in element, xyz)]
    end

    # Global DOF for each element (4-node (linear) quadrilateral element)
    dof = fill(Int[], size(con, 1))

    Threads.@threads for (i, elem) in collect(enumerate(con))
        dof[i] = [elem[1] * 2, elem[2] * 2 - 1, elem[2] * 2, elem[2] * 2 + 1,
                  elem[3] * 2, elem[3] * 2 + 1, elem[4] * 2, elem[4] * 2 + 1]
    end

    return dof
end

# New module shapefunctions

function planestrain(E, ν)
    constant = E / ((1 + ν) * (1 - (2*ν)))
    mat = [[1 - ν, ν, 0] [ν, 1 - ν, 0] [0, 0, 0.5*(1 - (2 * ν))]]
    C = constant * mat

    return C
end

function cmat(ξ, η)  # Come back to this
    ξ += 1
    cmat = 0.125 * (ξ + 3*η - ξ*η) + 0.625

    return cmat
end

function shapefuncs(ξ, η)
    N = zeros(Float64, 4)
    N[1] = 0.25 * (1 - ξ) * (1 + η)
    N[2] = 0.25 * (1 - ξ) * (1 - η)
    N[3] = 0.25 * (1 + ξ) * (1 + η)
    N[4] = 0.25 * (1 + ξ) * (1 - η)

    return N
end

# New module Elemental

using LinearAlgebra

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

function stiffmatrix(realcoors, Ce)
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
        dsB, J = straindisp(realcoors, ξ, η)

        # Jacobian determinant
        detJ = det(J)

        # Calculate K
        dsBT = dsB'
        dot1 = dsBT * Ce
        dot2 = dot1 * dsB
        Ke = Ke + dot2 * detJ * w
    end

    return Ke
end

function deformation(E, ν, fext, coors, C, Ke)
    # Boundary condition
    ndof = size(coors, 1)
    BC = ones(Int, (ndof, 2))

    # Find first and last indices in xyz between 1.5 and 2.5 in x
    y1 = findfirst(i -> i == 0, coors[:,2])
    y2 = findlast(j -> j == 0, coors[:,2])
    lb = findfirst(k -> (k >= 1.5) && (k <= 2.5), coors[:,1])
    ub = findlast(l -> (l >= 1.5) && (l <= 2.5), coors[y1:y2, 1])

    # Set these elements = 0
    BC[lb:ub] .= 0
    BC = reshape(BC', (length(BC), 1)) # TODO Is there a better way to do this?
    BCid = findall(!iszero, BC)

    # Apply the loads in the x direction to the top of the solid
    ytop = findall(i -> i == 5, coors[:,2])
    rhs = zeros((size(coors, 1), 2))
    Threads.@threads for j in ytop
        rhs[j,1] = fext
    end
    rhs1 = reshape(rhs', (length(rhs), 1)) # TODO Is there a better way to do this?

    # Non zero elements
    rhs2 = zeros(length(BCid))
    Threads.@threads for i = 1:length(BCid)
        rhs2[i] = rhs1[BCid[i]]
    end

    # TODO Find a better way of doing this
    iBCid = zeros(Int, length(BCid))
    for i = 1:length(BCid)
        iBCid[i] = i
    end
    npBCid = NumPyArray(iBCid)
    npKe = NumPyArray(Ke)
    # display(npBCid)
    # display(npKe)
    # TODO HELP!!!
    npKe = npKe[np.ix_(npBCid,npBCid)]

    ue = gesvx!(npKe, BCid)

    de = zeros(2 * length(coors))
    # de[BCid] = ue
end

function stress(realcoors, Ce, de)
    # Location of Gauss points
    a = 1/sqrt(3)

    # Weights of functions
    w = 1

    # Matrix of Gauss points
    gauss = [[-a, a, a, -a] [-a, -a, a, a]]

    σ = zeros(3)

    vol = 0

    for i in 1:4
        ξ = gauss[0,1]

        # TODO

    end

end


using Plots; gr()
using LaTeXStrings

function main(points::Int=1000)
    ### Generate mesh grid ###
    if points < 10
        println(points, " grid point have be specified.")
        println("This can not run with less than 10 grid points!")
        println("It is also recommended to use at least 100 grid points.")
        throw(DomainError(points, "points argument must be ≥ 30"))
    elseif points < 100
        println("It is recommended to use at least 100 grid points.")
    end

    # Corners of both parts of solid
    xstart1 = 1.5
    ystart1 = 0
    etop1 = 4
    eright1 = 2.5

    xstart2 = 0
    ystart2 = 4
    etop2 = 5
    eright2 = 4

    # Assign points to each solid part
    xpoints1 = pointsdiff(points, xstart1, eright1, ystart1, etop1)
    ypoints1 = (4 * (xpoints1 - 1))
    ypoints2 = xpoints1
    xpoints2 = (4 * xpoints1) - 3

    println(Int(xpoints1 * ypoints1) + Int(xpoints2 + ypoints2),
            " grid points used."); println()

    print("Calculating meshpoint arrays...")

    # Obtain mesh arrays for solid part edges
    bot1, top1, left1, right1 = mesharrays(Int(xpoints1), Int(ypoints1), xstart1, ystart1,
                                     etop1, eright1)
    bot2, top2, left2, right2 = mesharrays(Int(xpoints2), Int(ypoints2), xstart2, ystart2,
                                     etop2, eright2)

    println("Done")
    print("Calculating meshgrid coordinates and solid elements...")

    # Return grid of nodes for each solid part
    xyz1, elem1 = mesh(bot1, top1, left1, right1)
    xyz2, elem2 = mesh(bot2, top2, left2, right2)

    println("Done")

    # Solid nodes concatenated with duplicates removed
    xyz = union(xyz1, xyz2)

    # Axes for plotting nodes
    X = [i[1] for i in xyz]
    Y = [j[2] for j in xyz]

    # Combine elements from 2 solid shapes
    elem = vcat(elem1, elem2)

    print("Calculating connecting grid and elemental degrees of freedom...")

    dof = condof(elem, xyz)

    println("Done")
    print("Plotting meshgrid...")

    ### Plotting ###
    meshgrid = plot(
        X, Y,
        seriestype = :scatter,
        aspect_ratio = :equal,
        xlim = (-0.2, 4.2),
        ylim = (-0.2, 5.2),
        legend = :none,
        markersize = 1, markershape = :circle, markercolor = :black,
        xlabel = L"$x_1$",
        ylabel = L"$x_2$",
        title = "Meshgrid"
    )

    display(meshgrid)

    println("Done")

    ### Deformation ###
    # Material properties
    E = 1e7
    ν = 0.34

    # Applied load
    fext = -100

    print("Calculating plane strain...")

    Ce = planestrain(E, ν)

    println("Done")
    print("Calculating the per element stiffness matrix...")

    Ke = zeros(8, 8, length(elem))
    Threads.@threads for i = 1:length(elem)
        Ke[:,:,i] = stiffmatrix(elem[i], Ce)
    end

    println("Done")

    # Calculate the global stiffness matrix
    # Workaround for Julia's 1-array indexing

    print("Assembling global stiffness matrix...")

    Threads.@threads for i = 1:length(dof)
        dof[i] = dof[i] .- 1
    end

    globstiff = zeros(2*length(xyz), 2*length(xyz))

    Threads.@threads for i = 1:8
        for j = 1:8
            for m = 1:length(elem)
                globstiff[dof[m][j],dof[m][i]] += Ke[j,i,m]
            end
        end
    end

    println("Done")

    display(globstiff)

    # # TODO not working yet
    # def = deformation(E, ν, fext, xyz, Ce, Ke)

    # eps11 = 0
    # eps22 = 0
    # gamma12 = 0

    # # TODO not finished yet
    # σ = stress(xyz, Ce, de)

end

main(10000)
