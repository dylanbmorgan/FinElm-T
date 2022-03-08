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
    factor2 = ydir / xdir

    # Work out number of grid points in x and y directions
    if factor1 > 1
        factor = factor1
        ypoints = round(sqrt(totpoints / factor))
        xpoints = round(factor * ypoints)
    else
        factor = factor2
        xpoints = round(sqrt(totpoints / factor))
        ypoints = round(factor * xpoints)
    end

    return xpoints, ypoints
end

"""
    elements(xpoints, ypoints, xstart, ystart, etop, eright)

Create arrays of meshpoints for each side of the shape
"""
function elements(xpoints, ypoints, xstart, ystart, etop, eright)
    # Generate points in x-direction
    bot = zeros(xpoints)
    for (i,num) in enumerate(LinRange(xstart, eright, xpoints))
        bot[i] = num
    end

    # Generate points in y-direction
    left = zeros(ypoints)
    for (j,num) in enumerate(LinRange(ystart, etop, ypoints))
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

    # Generate coords of nodes 'xyz'
    xyz = zeros(nnodes, 2)
    Threads.@threads for i in 1:nnodesy
        yl = left[i] - left[1]
        dy = right[i] - left[i]

        for j in 1:nnodesx
            xb = bot[j] - bot[1] 
            dx = top[j] - bot[j]

            xcoor = (dx*yl + xb*ly) / (ly - dx*dy / lx)
            ycoor = dy / lx*xcoor + yl

            xyz[j + (i-1)*nnodesx, 1] = xcoor + bot[1]
            xyz[j + (i-1)*nnodesx, 2] = ycoor + left[1]
        end
    end

    # Node nums for elements
    nel = nelx * nely
    con = zeros(nel, 4)
    Threads.@threads for i in 1:nely
        for j in 1:nelx
            con[j + (i-1)*nelx, :] = [(j-1) + (i-1)*nnodesx, j + (i-1)*nnodesx,
                                      (j-1) + i*nnodesx + 1, (j-1) + i*nnodesx]
        end
    end

    # Global DOF for each element (4-node (linear) quadrilateral element)
    dof = zeros(Int, nel, 2 * 4)
    Threads.@threads for i in 1:nel
        dof[i, :] = [con[i, 1]*2, con[i, 2] * 2 - 1, con[i, 2]*2, con[i, 2] * 
                     2 + 1, con[i, 3]*2, con[i, 3] * 2 + 1, con[i, 4]*2, 
                     con[i, 4] * 2 + 1]
    end

    return xyz, con, dof
end

# New module shapefunctions

# TODO change to plane strain
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
    corners = [[rc[1,1],rc[end,1],rc[1,1],rc[end,1]];;
               [rc[1,2],rc[1,2],rc[end,2],rc[end,2]]]

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

    return dsB
end

function jacobian(ξ, η)
    J = [[1, 0] [0.125 - 0.125*η, 0.375 - ξ*0.125]]
    invJ = inv(J)

    return invJ
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
        dsB = straindisp(realcoors, ξ, η)

        # Jacobian and determinant
        J = jacobian(ξ, η)
        detJ = det(J)

        # Calculate K
        dsBT = dsB'
        sumdsBT = sum(dsBT, dims=2)
        dot1 = dsBT * Ce
        dot2 = dot1 * dsB
        Ke = Ke + dot2 * detJ * w
    end

    return Ke
end

function globstiff()
    # TODO 2, 3
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
    display(npBCid)
    display(npKe)
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
using LinearAlgebra

function main(points::Int=1000)
    ### Generate mesh grid ###
    if points < 5
        println(points, " grid point have be specified.")
        println("This can not run with less than 30 grid points!")
        println("It is also recommended to use at least 100 grid points.")
        throw(DomainError(points, "points argument must be ≥ 30"))
    elseif points < 100
        println("It is recommended to use at least 100 grid points.")
    end

    xstart = 0
    ystart = 0
    etop = 5
    eright = 4

    # Assign points to each solid element
    xpoints, ypoints = pointsdiff(points, xstart, eright, ystart, etop)

    println(Int(xpoints * ypoints), " grid points used.")

    # Obtain mesh arrays for solid element sides
    bot, top, left, right = elements(Int(xpoints), Int(ypoints), xstart, ystart,
                                     etop, eright)

    xyz, con, dof = mesh(bot, top, left, right)

    ### Plotting ###
    meshgrid = plot(
        xyz[:,1], xyz[:,2],
        seriestype = :scatter,
        aspect_ratio = :equal,
        xlim = (-0.2, 4.2),
        ylim = (-0.2, 5.2),
        legend = :none,
        markersize = 1, markershape = :x, markercolor = :black,
        xlabel = L"$x_1$",
        ylabel = L"$x_2$",
        title = "Meshgrid"
    )

    display(meshgrid)

    ### Deformation ###
    # Material properties
    E = 1e7
    ν = 0.34

    # Applied load
    fext = -100

    Ce = planestrain(E, ν)
    Ke = stiffmatrix(xyz, Ce)
    # TODO not working yet
    def = deformation(E, ν, fext, xyz, Ce, Ke)

    eps11 = 0
    eps22 = 0
    gamma12 = 0

    # TODO not finished yet
    σ = stress(xyz, Ce, de)

end

main()
# TODO remember to set limit back to 30
