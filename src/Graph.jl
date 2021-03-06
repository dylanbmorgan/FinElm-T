
###########################
#--------- Graph ---------#
###########################

using Plots; gr()
using PyCall
@pyimport matplotlib.pyplot as plt
@pyimport matplotlib.patches as patches
using LaTeXStrings


function graph(xyz, defxyz, elem, σ, ϵ, interactive, figpath)
    ### Original nodes ###
    X1 = [i[1] for i in xyz]
    Y1 = [j[2] for j in xyz]

    grid1 = plot(
        X1, Y1,
        seriestype = :scatter,
        markershape = :circle, markercolor = :blue,
        markerstrokewidth = 0, markersize = 1.5,
        aspect_ratio = :equal,
        framestyle = :origin,
        xlim = (-0.2, 4.2),
        ylim = (-0.2, 5.2),
        legend = :none,
        xlabel = L"$x_1$",
        ylabel = L"$x_2$",
        title = "Undeformed Grid",
        titlefontsize = 10
    )
    ###############

    ### Deformed nodes ###
    X2 = zeros(length(xyz))
    Y2 = copy(X2)
    Threads.@threads for (i, j) in collect(zip(1:2:length(defxyz), 1:length(xyz)))
        X2[j] = defxyz[i]
        Y2[j] = defxyz[i + 1]
    end

    grid2 = plot(
        X2, Y2,
        seriestype = :scatter,
        markershape = :circle, markercolor = :red,
        markerstrokewidth=0, markersize = 1.5,
        aspect_ratio = :equal,
        framestyle = :origin,
        xlim = (-0.2, 4.2),
        ylim = (-0.2, 5.2),
        legend = :none,
        xlabel = L"$x_1$",
        ylabel = L"$x_2$",
        title = "Scaled Deformed Grid",
        titlefontsize = 10
    )

    # Plot subplots
    disp = plot(grid1, grid2)

    # Change whether fig is saved or displayed depending on arg
    if interactive == 1
        savefig(disp, "$(figpath)displacement.png")
    elseif interactive == 2
        display(disp)
    elseif interactive == 3
        savefig(disp, "$(figpath)displacement.png")
        display(disp)
    end

    ######################

    ### Stresses ###
    # Split important bits of elem array up into a vector
    flatelem = zeros(length(elem) * 4)
    Threads.@threads for (i, j) in collect(zip(1:length(elem), 1:4:length(flatelem)))
        flatelem[j] = elem[i][1][1]
        flatelem[j+1] = elem[i][2][1]
        flatelem[j+2] = elem[i][3][1]
        flatelem[j+3] = elem[i][4][1]
    end

    lowerindex = trunc(Int, (findfirst(isequal(0), flatelem) - 1) / 4)

    # Lower half of shape
    X3L = zeros(lowerindex)
    Y3L = copy(X3L)
    Threads.@threads for i = 1:lowerindex
        xval = 0
        yval = 0

        for j = 1:4
            xval += elem[i][j][1]
            yval += elem[i][j][2]
        end

        X3L[i] = sum(xval) / 4
        Y3L[i] = sum(yval) / 4
    end

    # Upper half of shape
    upperindex = lowerindex + 1
    updiff = length(elem) - upperindex + 1
    X3U = zeros(updiff)
    Y3U = copy(X3U)
    Threads.@threads for (i, j) in collect(zip(upperindex:length(elem), 1:updiff))
        xval = 0
        yval = 0

        for k = 1:4
            xval += elem[i][k][1]
            yval += elem[i][k][2]
        end

        X3U[j] = sum(xval) / 4
        Y3U[j] = sum(yval) / 4
    end

    X3 = vcat(X3L, X3U)
    Y3 = vcat(Y3L, Y3U)

    # Plot
    # Using matplotlib for this from a python call isn't ideal
    # but it's probably the best/easiest way to plot it
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))

    # Components of stress vector
    g1 = ax[1,1].tricontourf(X3, Y3, σ[:,1],
                             levels=LinRange(findmin(σ[:,1])[1], findmax(σ[:,1])[1], 20))
    patch1 = patches.Rectangle((0, 0), X3[1], Y3U[1], fc="white")
    patch2 = patches.Rectangle((findmax(X3L)[1], 0), 4, Y3U[1], fc="white")
    ax[1,1].add_patch(patch1)
    ax[1,1].add_patch(patch2)
    ax[1,1].set_title(L"σ_{11}", fontsize=15)
    fig.colorbar(g1, ax=ax[1,1], shrink=0.7)

    g2 = ax[1,2].tricontourf(X3, Y3, σ[:, 2],
                             levels=LinRange(findmin(σ[:,2])[1], findmax(σ[:,2])[1], 20))
    patch3 = patches.Rectangle((0, 0), X3[1], Y3U[1], fc="white")
    patch4 = patches.Rectangle((findmax(X3L)[1], 0), 4, Y3U[1], fc="white")
    ax[1,2].add_patch(patch3)
    ax[1,2].add_patch(patch4)
    ax[1,2].set_title(L"σ_{22}", fontsize=15)
    fig.colorbar(g2, ax=ax[1,2], shrink=0.7)

    g3 = ax[1,3].tricontourf(X3, Y3, σ[:, 3],
                             levels=LinRange(findmin(σ[:,3])[1], findmax(σ[:,3])[1], 20))
    patch5 = patches.Rectangle((0, 0), X3[1], Y3U[1], fc="white")
    patch6 = patches.Rectangle((findmax(X3L)[1], 0), 4, Y3U[1], fc="white")
    ax[1,3].add_patch(patch5)
    ax[1,3].add_patch(patch6)
    ax[1,3].set_title(L"σ_{12}", fontsize=15)
    fig.colorbar(g3, ax=ax[1,3], shrink=0.7)

    # Components of strain vector
    g4 = ax[2,1].tricontourf(X3, Y3, ϵ[:, 1],
                             levels=LinRange(findmin(ϵ[:,1])[1], findmax(ϵ[:,1])[1], 20))
    patch7 = patches.Rectangle((0, 0), X3[1], Y3U[1], fc="white")
    patch8 = patches.Rectangle((findmax(X3L)[1], 0), 4, Y3U[1], fc="white")
    ax[2,1].add_patch(patch7)
    ax[2,1].add_patch(patch8)
    ax[2,1].set_title(L"ϵ_{11}", fontsize=15)
    fig.colorbar(g4, ax=ax[2,1], shrink=0.7)

    g5 = ax[2,2].tricontourf(X3, Y3, ϵ[:, 2],
                             levels=LinRange(findmin(ϵ[:,2])[1], findmax(ϵ[:,2])[1], 20))
    patch9 = patches.Rectangle((0, 0), X3[1], Y3U[1], fc="white")
    patch10 = patches.Rectangle((findmax(X3L)[1], 0), 4, Y3U[1], fc="white")
    ax[2,2].add_patch(patch9)
    ax[2,2].add_patch(patch10)
    ax[2,2].set_title(L"ϵ_{22}", fontsize=15)
    fig.colorbar(g5, ax=ax[2,2], shrink=0.7)

    g6 = ax[2,3].tricontourf(X3, Y3, ϵ[:, 3],
                             levels=LinRange(findmin(ϵ[:,3])[1], findmax(ϵ[:,3])[1], 20))
    patch11 = patches.Rectangle((0, 0), X3[1], Y3U[1], fc="white")
    patch12 = patches.Rectangle((findmax(X3L)[1], 0), 4, Y3U[1], fc="white")
    ax[2,3].add_patch(patch11)
    ax[2,3].add_patch(patch12)
    ax[2,3].set_title(L"γ_{12}", fontsize=15)
    fig.colorbar(g6, ax=ax[2,3], shrink=0.7)

    # Adjust some parameters
    for i in ax
        i.set_xlabel(L"$x_1$", fontsize=15)
        i.set_ylabel(L"$x_2$", fontsize=15)
        i.set_xlim(0, 4)
        i.set_ylim(0, 5)
        i.set_aspect("equal")
    end

    fig.tight_layout()  # Better spacing between subplots

    # Change whether fig is saved or displayed depending on arg
    if interactive == 1
        plt.savefig("$(figpath)strainstress.png")
    elseif interactive == 2
        plt.show()
    elseif interactive == 3
        plt.savefig("$(figpath)strainstress.png")
        plt.show()
    end
end

function plottime()
    # Data from runs
    runtime = [2.49, 2.51, 2.57, 2.64, 3.97, 9.70, 30.74, 60.47, 110.09, 325.88, 778.40, 1214.02]
    nnodes = [160, 337, 576, 882, 2176, 4060, 6506, 8448, 10665, 15861, 20420, 23780]

    # Plot
    display(plot(
        nnodes, runtime,
        yaxis = :log,
        marker = :circle,
        markersize = 3,
        legend = :none,
        xlabel = "Number of Nodes",
        ylabel = "Total Run Time (log) / sec"
    ))
end

function plotconv()
    # Data from runs
    meanD = [5.37614e-13, 7.66974e-13, 1.02092e-12, 1.28569e-12, 1.9914e-12, 2.97906e-12, 3.54956e-12,
             3.9345e-12, 4.48195e-12, 5.64084e-12, 6.3608e-12, 6.68799e-12]
    maxD = [1.69578e-12, 2.48807e-12, 3.18746e-12, 4.02003e-12, 6.17745e-12, 9.65003e-12, 1.11053e-11,
            1.21609e-11, 1.39048e-11, 1.77602e-11, 2.00271e-11, 2.06977e-11]
    nnodes = [160, 337, 576, 882, 2176, 4060, 6506, 8448, 10665, 15861, 20420, 23780]

    # plot
    graph = plot(
        nnodes, meanD,
        marker = :circle,
        markersize = 3, markerstrokewidth = 0.5,
        legend = :right,
        label = "Average Displacement",
        ylim = (1e-13, 2.1e-11),
        yticks = (1e-13:2.5e-12:2.5e-11),
        xlabel = "Number of Nodes",
        ylabel = "Displacement / m",
    )

    plot!(
        nnodes, maxD,
        marker = :square,
        markersize = 2.4, markerstrokewidth = 0.5,
        label = "Maximum Displacement"
    )

    # display(graph)
    savefig("/home/dylanmorgan/Documents/warwick/chapter_1/px912/assignments/solids/report/figures/d_convergence.png")
end
