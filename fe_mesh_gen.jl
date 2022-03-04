using Plots

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
    for i in 1:nnodesy
        yl = left[i] - left[1]
        dy = right[i] - left[i] 

        for j in 1:nnodesx
            xb = bot[j] - bot[1]
            dx = top[j] - bot[j]

            xcoor = (dx*yl + xb*ly) / (ly - dx*dy / lx)
            ycoor = dy / lx*xcoor + yl

            xyz[j + (i-1)*nnodesx, 1] = xcoor + left[1]
            xyz[j + (i-1)*nnodesx, 2] = ycoor + bot[1]
        end
    end

    # Node nums for elements
    nel = nelx * nely
    con = zeros(nel, 4)
    for i in 1:nely
        for j in 1:nelx
            con[j + (i-1)*nelx, :] = [(j-1) + (i-1)*nnodesx, j + (i-1)*nnodesx,
                                      (j-1) + i*nnodesx + 1, (j-1) + i*nnodesx]
        end
    end

    # Global DOFF for each element (4-node (linear) quadrilateral element)
    dof = zeros(Int, nel, 2 * 4)
    for i in 1:nel
        dof[i, :] = [con[i, 1]*2, con[i, 2] * 2 - 1, con[i, 2]*2, con[i, 2] * 
                     2 + 1, con[i, 3]*2, con[i, 3] * 2 + 1, con[i, 4]*2, 
                     con[i, 4] * 2 + 1]
    end

    return xyz, con, dof
end

xpoints = 20; ypoints = 20

# Split num of points by area of each part of shape
x1points = round(xpoints / 2)
x2points = xpoints - x1points
y1points = round(ypoints / 2)
y2points = ypoints - y1points

bot1, top1, left1, right1 = elements(Int(x1points), Int(y1points), 0, 4, 5, 4)
bot2, top2, left2, right2 = elements(Int(x2points), Int(y2points), 1.5, 0, 4, 2.5)

xyz1, con1, dof1 = mesh(bot1, top1, left1, right1)
xyz2, con2, dof2 = mesh(bot2, top2, left2, right2)

display(scatter(
    xyz1[:, 1], xyz1[:, 2]
))

display(scatter(
    xyz2[:, 1], xyz2[:, 2]
))
