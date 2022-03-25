####################################
#--------- Generate Mesh ----------#
####################################

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
