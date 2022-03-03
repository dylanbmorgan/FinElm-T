module MeshGen
export mesh

"""
    mesh(bot, top, left, right)

Generate a finite element mesh grid for some 2D object
"""
function mesh(bot, top, left, right)
    # Number of nodes and elements in the domain 
    nnodesx = size(bot, 1)  # TODO Set type of size
    nnodesy = size(left, 1)
    nelx = nnodesx - 1
    nely = nnodesy - 1
    nnodes = nnodesx * nnodesy

    # Dimensions of the domain 
    lx = bot[nnodesx - 1] - bot[1]
    ly = left[nnodesy - 1] - left[1]

    # Generate coords of nodes 'xyz'
    xyz = zeros(nnodes, 2)
    display(xyz)
    for i in 1:nnodesy
        yl = left[i] - left[1]
        dy = right[i] - left[i] 

        for j in 1:nnodesx
            xb = bot[j] - bot[1]
            dx = top[j] - bot[j]

            xcoor = (dx*yl + xb * ly) / (ly - dx*dy / lx)
            ycoor = dy / lx*xcoor + yl

            println(j + i * nnodesx)
            xyz[j + i*nnodesx, 1] = ycoor + left[1]
            xyz[j + i*nnodesx, 2] = xcoor + bot[1]
        end
    end

    # Node nums for elements
    nel = nelx * nely
    con = zeros(nel, 4)
    for i in 1:nely
        for j in 1:nelx
            con[j + i*nelx, :] = [j + i*nnodesx, j + i*nnodesx + 1, j + (i + 1)*
                nnodesx + 1, j + (i + 1)*nnodesx]
        end
    end

    # Global DOF for each element (4-node (linear) quadrilateral element)
    # dof = zeros(Int, nel, 2 * 4)
    # for i in range(nel)
    #     dof[i, :] = [con[i, 1]*2, con[i, 2] * 2 - 1, con[i, 2]*2, con[i, 2] * 2
    #                  + 1, con[i, 3]*2, con[i, 3] * 2 + 1, con[i, 4]*2, con[i, 4]
    #                  * 2 + 1]
    # end

    return xyz, con, dof
end

end
