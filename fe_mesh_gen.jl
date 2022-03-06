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

function planestress(E, ν)
    # Calculate plane stress
    constant = E / (1 - (ν*ν))
    mat = Array{Float64, [[1, ν, 0] [ν, 1, 0] [0, 0, 0.5*(1 - ν)]]}
    C = constant * mat

    return C
end

function cmat(ξ, η)
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

# New module elemental stiffness

using LinearAlgebra

function strain_displacement(realcoors, ξ, η)
    # Natural coors of quadrilateral
    natcoord = Array[[-1, 1, 1, -1] [-1, -1, 1, 1]]

    # Derivatives of shape functions wrt natural coords
    dNdnat = zeros((2, 4))
    dNdnat[1,:]= 0.25 * natcoord[1,:] * (1+natcoord[2,:] * η)
    dNdnat[2,:]= 0.25 * natcoord[2,:] * (1+natcoord[1,:] * ξ)

    # Elemental Jacobian matrix
    Jmat = dot(dNdnat,realcoors)
    J = det(Jmat)

    JmatInv = inv(Jmat)
    dNdx = dot(JmatInv,dNdnat)

    # TODO Implement the following python code
    # dsB=np.zeros((3,8))
    # dsB[0,0::2]=dNdx[0,:]
    # dsB[1,1::2]=dNdx[1,:]
    # dsB[2,0::2]=dNdx[1,:]
    # dsB[2,1::2]=dNdx[0,:]

    # return dsB
end

using Plots; gr()
using LaTeXStrings

function main(points=20)
    if points < 10
        println(points, " grid point(s) have be specified.")
        println("This can not run with less than 10 grid points!")
        println("It is also recommended to use at least 100 grid points.")
        throw(DomainError(points, "points argument must be ≥ 10"))
    elseif points < 100
        println("It is recommended to use at least 100 grid points.")
        println(points, " points have been specified here.")
    end

    xstart = 0
    ystart = 0
    etop = 5
    eright = 4

    # Assign points to each solid element
    xpoints, ypoints = pointsdiff(points, xstart, eright, ystart, etop)

    if points >= 100
        println(Int(xpoints * ypoints), " grid points used.")
    end

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
        legend = :none,
        markersize = 1, markershape = :x, markercolor = :black,
        xlabel = L"$x_1$",
        ylabel = L"$x_2$",
        title = "Meshgrid"
    )

    display(meshgrid)
end

main()
