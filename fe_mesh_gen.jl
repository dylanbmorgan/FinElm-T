using Plots

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

    # Global DOF for each element (4-node (linear) quadrilateral element)
    dof = zeros(Int, nel, 2 * 4)
    for i in 1:nel
        dof[i, :] = [con[i, 1]*2, con[i, 2] * 2 - 1, con[i, 2]*2, con[i, 2] * 
                     2 + 1, con[i, 3]*2, con[i, 3] * 2 + 1, con[i, 4]*2, 
                     con[i, 4] * 2 + 1]
    end

    return xyz, con, dof
end

function elements(verts, npoints)
    
end

# Shape verticies 
#top = [-2 2]
#bottom = [-2 -0.5 0.5 2] 
#left = [0 4 5]
#right = [0 4 5]

# Number of nodes in each direction
#horiznodes = 10
#vertnodes = 5

bot = [0, 1, 2, 3, 4]
top = [0, 1, 2, 3, 4]
left = [0, 1, 2, 3, 4, 5]
right = [0, 1, 2, 3, 4, 5]


# TEST 8 elements
#bot = [0, 0.5, 1, 1.5, 2]
#top = [0, 0.5, 1, 1.5, 2]
#left = [0, 0.5, 1]
#right = [0.5, 0.75, 1]

#generate mesh
xyz, con, dof = mesh(bot,top,left,right)

display(scatter(
    xyz[:, 1], xyz[:, 2]
))
