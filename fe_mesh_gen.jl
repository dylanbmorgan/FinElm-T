#!/usr/bin/env julia

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
    # x-direction
    bot = zeros(xpoints)
    for (i,num) in enumerate(LinRange(xstart, eright, xpoints))
        bot[i] = num
    end

    # y-direction
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

# TODO Split into modules in different files
using Plots; gr()

function main(points=1000)
    if points < 10
        println(points, " grid point(s) have be specified.")
        println("This can not run with less than 10 grid points!")
        println("It is also recommended to use at least 100 grid points.")
        throw(DomainError(points, "points argument must be ≥ 10"))
    elseif points < 100
        println("It is recommended to use at least 100 grid points.")
        println(points, " points have been used here.")
    end

    xstart1 = 0; xstart2 = 1.5
    ystart1 = 4; ystart2 = 0
    etop1 = 5; etop2 = 4
    eright1 = 4; eright2 = 2.5

    # Assign points to each solid element
    x1points, y1points = pointsdiff(points, xstart1, eright1, ystart1, etop1)
    x2points, y2points = pointsdiff(points, xstart2, eright2, ystart2, etop2)

    # Remove top row from lower element
    y2range = LinRange(xstart2, eright2, Int(x2points))
    reinterpret(Float64, y2range)
    y2pen = y2range[end] - y2range[end - 1]

    # Obtain mesh arrays for solid element sides
    bot1, top1, left1, right1 = elements(Int(x1points), Int(y1points), xstart1,
                                         ystart1, etop1, eright1)
    bot2, top2, left2, right2 = elements(Int(x2points), Int(y2points), xstart2,
                                         ystart2, etop2 - y2pen, eright2)

    xyz1, con1, dof1 = mesh(bot1, top1, left1, right1)
    xyz2, con2, dof2 = mesh(bot2, top2, left2, right2)

    ### Plotting ###
    meshgrid = plot(
        xyz1[:,1], xyz1[:,2],
        seriestype = :scatter,
        aspect_ratio = :equal,
        xlim = (-0.2, 4.2),
        legend = :none,
        markersize = 1, markershape = :circle, markercolor = :black,
        title = "Meshgrid"
    )

    plot!(meshgrid,
          xyz2[:,1], xyz2[:,2],
          seriestype = :scatter,
          markersize = 1, markershape = :circle, markercolor = :black
          )

    display(meshgrid)
end

# TODO Parse arguments from the REPL
main()
