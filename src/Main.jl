##########################
#--------- MAIN ---------#
##########################

include("MeshGen.jl")
include("Deformation.jl")
include("StrainStress.jl")
include("Graph.jl")

using LinearAlgebra
using Statistics
using PrettyTables


"""
    finelm(time, points, scalefac, interactive, writefile, filename, figpath)

Perform finite element analysis of a T.
"""
function finelm(
    time::Bool=false,
    points::Int=2000,
    scalefac::Float64=1e11,
    interactive::Int=1,
    writefile::Bool=false,
    filename::String="../output/run.out",
    figpath::String="../output/"
    )

    # Check supplied number of nodes
    if points < 20
        println(points, " grid point have be specified.")
        println("This can not run with less than 20 grid points!")
        println("It is also recommended to use at least 100 grid points.")
        throw(DomainError(points, "points argument must be ≥ 20"))
    elseif points < 100
        println("Warning! It is recommended to use at least 100 grid points.")
        println()
    end

    # Write output to file if arg is supplied
    if writefile == true
        print("Output is being redirected to a file...")
        open(filename, "w") do io
            redirect_stdout(io) do
                if time == true
                    @time finelm(time, points, scalefac, interactive, false, filename, figpath)
                else
                    finelm(time, points, scalefac, interactive, false, filename, figpath)
                end
            end
        end
        println("Done")
        println()
        print("Time taken:")
        return nothing
    end

    ### Generate mesh grid ###
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

    println("################")
    println(length(xyz), " nodes used")
    println("################")
    println()

    # Combine elements from 2 solid shapes
    elem = vcat(elem1, elem2)

    print("Calculating connecting grid and elemental degrees of freedom...")

    dof = condof(elem, xyz)

    println("Done")

    ### Deformation ###
    # Material properties
    E = 1e10
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

    # Solve linearly
    Dred = (Fred\inv(Kred))'

    # Add back fixed nodes
    D = zeros(length(Dred) + length(bc))
    Threads.@threads for i = 1:length(bc)
        D[i] = 0
    end

    Threads.@threads for (i, j) = collect(zip(length(bc):length(D), 1:length(Dred)))
        D[i+1] = Dred[j]
    end

    println("done")
    print("Deforming nodes...")

    defxyz = deform(xyz, D, scalefac)

    println("Done")

    ### Calculate stress and strain ###
    print("Calculating stresses and strains...")

    # Get elemental displacements
    de = zeros(length(dof), 8)
    Threads.@threads for (i, elem) in collect(enumerate(dof))
        for (j, corn) in enumerate(elem)
            de[i,j] = D[corn]
        end
    end

    # convert de to vector of vectors
    de = [de[i,:] for i in 1:size(de, 1)]

    σ = zeros(length(elem), 3)
    ϵ = zeros(length(elem), 3)

    Threads.@threads for i = 1:length(elem)
        σ[i,1], σ[i,2], σ[i,3], ϵ[i,1], ϵ[i,2], ϵ[i,3] = strainstress(elem[i], de[i], C)
    end

    # Calculate average, max, and min of stress, strain, and deformation
    tableD = []
    σ1 = []
    σ2 = []
    σ3 = []
    ϵ1 = []
    ϵ2 = []
    ϵ3 = []

    append!(tableD, round(mean(D), sigdigits=6))
    append!(tableD, round(findmax(D)[1], sigdigits=6))
    append!(tableD, round(findmin(D)[1], sigdigits=6))

    for i = 1:3
        append!(σ1, round(mean(σ[:,i]), sigdigits=6))
        append!(σ2, round(findmax(σ[:,i])[1], sigdigits=6))
        append!(σ3, round(findmin(σ[:,i])[1], sigdigits=6))
        append!(ϵ1, round(mean(ϵ[:,i]), sigdigits=6))
        append!(ϵ2, round(findmax(ϵ[:,i])[1], sigdigits=6))
        append!(ϵ3, round(findmin(ϵ[:,i])[1], sigdigits=6))
    end

    prepend!(tableD, ["D"])
    prepend!(σ1, ["σ_11"])
    prepend!(σ2, ["σ_12"])
    prepend!(σ3, ["σ_22"])
    prepend!(ϵ1, ["ϵ_11"])
    prepend!(ϵ2, ["ϵ_22"])
    prepend!(ϵ3, ["γ_12"])

    # Format data to be displayed in table
    tabledata = permutedims([tableD;; σ1;; σ2;; σ3;; ϵ1;; ϵ2;; ϵ3])

    println("Done")

    ### Plotting ###
    print("Plotting...")

    graph(xyz, defxyz, elem, σ, ϵ, interactive, figpath)

    println("Done")

    # Output data
    println()
    println("###########################")
    println("### Computed Quantities ###")
    println("###########################")
    pretty_table(
        tabledata;
        header = (["Quantity", "Mean", "Max", "Min"])
    )
    println()
    println("###########################################")
    println("### Computed Quantities in LaTeX format ###")
    println("###########################################")
    pretty_table(
        tabledata;
        header = (["Quantity", "Mean", "Max", "Min"]),
        backend = Val(:latex)
    )

    # Display some run information
   if time == true
        println()
        println("Run information:")
        print("Total runtime: ")
    end
end
