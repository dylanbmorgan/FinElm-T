
###########################
#--------- Graph ---------#
###########################

using Plots; gr()
using PyCall; @pyimport matplotlib.pyplot as plt
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
        title = "Deformed Grid",
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
    for (i, j) in zip(1:length(elem), 1:4:length(flatelem))
        flatelem[j] = elem[i][1][1]
        flatelem[j+1] = elem[i][2][1]
        flatelem[j+2] = elem[i][3][1]
        flatelem[j+3] = elem[i][4][1]
    end

    lowerindex = trunc(Int, (findfirst(isequal(0), flatelem) - 1) / 4)

    # Lower half of shape
    X3L = zeros(lowerindex)
    Y3L = copy(X3L)
    for i = 1:lowerindex
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
    for (i, j) = zip(upperindex:length(elem), 1:updiff)
        xval = 0
        yval = 0

        for k = 1:4
            xval += elem[i][k][1]
            yval += elem[i][k][2]
        end

        X3U[j] = sum(xval) / 4
        Y3U[j] = sum(yval) / 4
    end

    # Plot
    # Using matplotlib for this from a python call isn't ideal
    # but it's probably the best/easiest way to plot it
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))

    # Components of stress vector
    g1 = ax[1,1].tricontourf(X3L, Y3L, σ[1:lowerindex, 1],
                             levels=LinRange(findmin(σ[:,1])[1], findmax(σ[:,1])[1], 16))
    ax[1,1].tricontourf(X3U, Y3U, σ[upperindex:length(elem), 1],
                        levels=LinRange(findmin(σ[:,1])[1], findmax(σ[:,1])[1], 16))
    ax[1,1].set_title(L"σ_{11}", fontsize=20)
    fig.colorbar(g1, ax=ax[1,1], shrink=0.7)

    g2 = ax[1,2].tricontourf(X3L, Y3L, σ[1:lowerindex, 2],
                             levels=LinRange(findmin(σ[:,2])[1], findmax(σ[:,2])[1], 16))
    ax[1,2].tricontourf(X3U, Y3U, σ[upperindex:length(elem), 2],
                        levels=LinRange(findmin(σ[:,2])[1], findmax(σ[:,2])[1], 16))
    ax[1,2].set_title(L"σ_{22}", fontsize=20)
    fig.colorbar(g2, ax=ax[1,2], shrink=0.7)

    g3 = ax[1,3].tricontourf(X3L, Y3L, σ[1:lowerindex, 3],
                             levels=LinRange(findmin(σ[:,3])[1], findmax(σ[:,3])[1], 16))
    ax[1,3].tricontourf(X3U, Y3U, σ[upperindex:length(elem), 3],
                        levels=LinRange(findmin(σ[:,3])[1], findmax(σ[:,3])[1], 16))
    ax[1,3].set_title(L"σ_{12}", fontsize=20)
    fig.colorbar(g3, ax=ax[1,3], shrink=0.7)

    # Components of strain vector
    g4 = ax[2,1].tricontourf(X3L, Y3L, ϵ[1:lowerindex, 1],
                             levels=LinRange(findmin(ϵ[:,1])[1], findmax(ϵ[:,1])[1], 16))
    ax[2,1].tricontourf(X3U, Y3U, ϵ[upperindex:length(elem), 1],
                        levels=LinRange(findmin(ϵ[:,1])[1], findmax(ϵ[:,1])[1], 16))
    ax[2,1].set_title(L"ϵ_{11}", fontsize=20)
    fig.colorbar(g4, ax=ax[2,1], shrink=0.7)

    g5 = ax[2,2].tricontourf(X3L, Y3L, ϵ[1:lowerindex, 2],
                             levels=LinRange(findmin(ϵ[:,2])[1], findmax(ϵ[:,2])[1], 16))
    ax[2,2].tricontourf(X3U, Y3U, ϵ[upperindex:length(elem), 2],
                        levels=LinRange(findmin(ϵ[:,2])[1], findmax(ϵ[:,2])[1], 16))
    ax[2,2].set_title(L"ϵ_{22}", fontsize=20)
    fig.colorbar(g5, ax=ax[2,2], shrink=0.7)

    g6 = ax[2,3].tricontourf(X3L, Y3L, ϵ[1:lowerindex, 3],
                             levels=LinRange(findmin(ϵ[:,3])[1], findmax(ϵ[:,3])[1], 16))
    ax[2,3].tricontourf(X3U, Y3U, ϵ[upperindex:length(elem), 3],
                        levels=LinRange(findmin(ϵ[:,3])[1], findmax(ϵ[:,3])[1], 16))
    ax[2,3].set_title(L"γ_{12}", fontsize=20)
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
