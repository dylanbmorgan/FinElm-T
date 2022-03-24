#!/usr/bin/env julia

###########################################
#--------- Module Generate Mesh ----------#
###########################################

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

########################################
#--------- Module Deformation ---------#
########################################

function planestrain(E, ν)
    constant = E / ((1 + ν) * (1 - (2*ν)))
    mat = [[1 - ν, ν, 0] [ν, 1 - ν, 0] [0, 0, 0.5*(1 - (2 * ν))]]
    C = constant * mat

    return C
end


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
        Ke = Ke + dot2 * detJ * w
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

function deform(xyz, D, elem, de)
    # Deform nodes
    # Convert xyz to match dims of D
    defxyz = zeros(length(xyz) * 2)
    Threads.@threads for (i, j) in collect(zip(1:2:length(xyz)*2, 1:length(xyz)))
        defxyz[i] = xyz[j][1]
        defxyz[i + 1] = xyz[j][2]
    end

    # Deform nodes
    defxyz .+= D

    # Deform elements
    # Convert elem to match dims of de
    flatelem = fill(Float64[], length(elem))
    Threads.@threads for i = 1:length(elem)
        tmpelem = zeros(8)

        for (j, k) = zip(1:4, 1:2:8)
            tmpelem[k] = elem[i][j][1]
            tmpelem[k + 1] = elem[i][j][2]
        end

        flatelem[i] = tmpelem
        flatelem[i] += de[i]
    end

    return defxyz, flatelem
end

##########################################
#--------- Module strain stress ---------#
##########################################

# Needs to link with deform module to get straindisp function

function strainstress(elem, de, C)
    # Location of Gauss points
    a = 1/sqrt(3)

    # Weights of functions
    w = 1

    # Matrix of Gauss points
    gauss = [[-a, a, a, -a] [-a, -a, a, a]]

    # Vector of elemental stresses
    σe = zeros(3)

    for i = 1:4
        # Natural coors
        ξ = gauss[i,1]
        η = gauss[i,2]

        # B matrix
        Be, J = straindisp(elem, ξ, η)

        # Jacobian determinant
        detJ = det(J)

        # Elemental IP strain
        ϵeIP = Be * de

        # Elemental IP stress
        σeIP = C * ϵeIP

        σe += σeIP * detJ * w
    end

    ϵe = inv(C) * σe

    return σe[1], σe[2], σe[3], ϵe[1], ϵe[2], ϵe[3]
end

#################################
#--------- Module plot ---------#
#################################

using Plots; gr()
using LaTeXStrings
using Colors
using Statistics

function graph(elem, xyz, σ)
    ### Nodes ###
    # X1 = [i[1] for i in xyz]
    # Y1 = [j[2] for j in xyz]

    X1 = zeros(Int(length(xyz) / 2))
    Y1 = zeros(Int(length(xyz) / 2))
    for (i, j) in zip(1:2:length(xyz), 1:Int(length(xyz) / 2))
        X1[j] = xyz[i]
        Y1[j] = xyz[i + 1]
    end

    grid = plot(
        X1, Y1,
        seriestype = :scatter,
        markersize = 1.5, markershape = :circle, markercolor = :blue,
        aspect_ratio = :equal,
        framestyle = :origin,
        # xlim = (-0.2, 4.2),
        # ylim = (-0.2, 5.2),
        legend = :none,
        xlabel = L"$x_1$",
        ylabel = L"$x_2$",
        title = "Elements"
    )
    ###############

    ### Elements ###
    X2 = fill(Float64[], length(elem))
    Y2 = copy(X2)

    for i = 1:length(elem)
        X2[i] = zeros(4)
        Y2[i] = zeros(4)
    end

    # TODO write comments in this function
    for i = 1:length(elem)
        for (j, k) in zip(1:2:8, 1:4)
            X2[i][k] = elem[i][j]
            Y2[i][k] = elem[i][j + 1]
        end
    end

    # test = range(HSV(0,1,1), stop=HSV(360,1,1), length=90)
    # display(σ)

    for (x, y) in (zip(X2, Y2))
        plot!(
            Shape([x[1], x[2], x[3], x[4]], [y[1], y[2], y[3], y[4]]),
            fill = nothing #test[count]
        )
    end
    ###############

    display(grid)


    # X3 = zeros(length(X1))
    # Y3 = zeros(length(Y1))
    # for i = 1:length(X1)
    #     X3[i] = mean(X1[i])
    #     Y3[i] = mean(Y1[i])
    # end

    # Z3 = σ[:,1]

    # Z3 = zeros(length(X3), length(Y3))
    # for i = 1:length(X3)
    #     for j = 1:length(Y3)
    #         Z3[i,j] = i*j
    #     end
    # end

end

##################################
#--------- MAIN PROGRAM ---------#
##################################

function main(points::Int=1000)
    ### Generate mesh grid ###
    if points < 20
        println(points, " grid point have be specified.")
        println("This can not run with less than 20 grid points!")
        println("It is also recommended to use at least 100 grid points.")
        throw(DomainError(points, "points argument must be ≥ 20"))
    elseif points < 100
        println("Warning! It is recommended to use at least 100 grid points.")
        println()
    end

    points = points / 2.6

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

    print("Calculating meshpoint arrays...")

    # Obtain mesh arrays for solid part edges
    bot1, top1, left1, right1 = mesharrays(Int(xpoints1), Int(ypoints1), xstart1, ystart1, etop1, eright1)
    bot2, top2, left2, right2 = mesharrays(Int(xpoints2), Int(ypoints2), xstart2, ystart2, etop2, eright2)

    println("Done")
    print("Calculating meshgrid coordinates and solid elements...")

    # Return grid of nodes for each solid part
    xyz1, elem1 = mesh(bot1, top1, left1, right1)
    xyz2, elem2 = mesh(bot2, top2, left2, right2)

    println("Done"); println()

    # Solid nodes concatenated with duplicates removed
    xyz = union(xyz1, xyz2)

    println(length(xyz), " nodes used"); println()

    # Combine elements from 2 solid shapes
    elem = vcat(elem1, elem2)

    print("Calculating connecting grid and elemental degrees of freedom...")

    dof = condof(elem, xyz)

    println("Done")

    ### Deformation ###
    # Material properties
    E = 1e7
    ν = 0.31

    # Applied load
    fext = -100

    print("Calculating plane strain...")

    C = planestrain(E, ν)

    println("Done")
    print("Calculating the per element stiffness matrix...")

    Ke = zeros(8, 8, length(elem))
    Threads.@threads for i = 1:length(elem)
        Ke[:,:,i] = stiffmatrix(elem[i], C)
    end

    println("Done")
    print("Assembling global stiffness matrix...")

    Threads.@threads for i = 1:length(dof)
        dof[i] = dof[i] .- 1
    end

    # Global stiffness matrix
    K = zeros(2*length(xyz), 2*length(xyz))

    Threads.@threads for i = 1:8
        for j = 1:8
            for m = 1:length(elem)
                K[dof[m][j],dof[m][i]] += Ke[j,i,m]
            end
        end
    end

    println("Done")
    print("Calculating the force vector...")

    # Applied force per unit area
    F = forcevec(xyz, etop2, eright2, fext)

    println("Done")
    print("Reducing the force vector and global stiffness matrix...")

    xyzcopy = copy(xyz)
    Fred, Kred, bc = reducemat!(xyzcopy, K, F)

    println("done")
    print("Solving for D linearly...")

    # Solve
    Dred = (Fred\Kred)'

    # Add back fixed nodes
    D = zeros(length(Dred) + length(bc))
    Threads.@threads for i = 1:length(bc)
        D[i] = 0
    end

    Threads.@threads for (i, j) = collect(zip(length(bc):length(D), 1:length(Dred)))
        D[i + 1] = Dred[j]
    end

    println("done")
    print("Deforming nodes and elements...")

    de = zeros(length(dof), 8)
    Threads.@threads for (i, elem) in collect(enumerate(dof))
        for (j, corn) in enumerate(elem)
            de[i,j] = D[corn]
        end
    end

    # convert de to vector of vectors
    de = [de[i,:] for i in 1:size(de, 1)]

    defxyz, flatelem = deform(xyz, D, elem, de)

    println("Done")

    ### Calculate stress and strain ###
    print("Calculating stresses and strains...")

    σ = zeros(length(elem), 3)
    ϵ = zeros(length(elem), 3)

    Threads.@threads for i = 1:length(elem)
        σ[i,1], σ[i,2], σ[i,3], ϵ[i,1], ϵ[i,2], ϵ[i,3] = strainstress(elem[i], de[i], C)
    end

    println("Done")

    ### Plotting ###
    print("Plotting...")

    graph(flatelem, defxyz, σ)

    println("Done")

end

main(500)
