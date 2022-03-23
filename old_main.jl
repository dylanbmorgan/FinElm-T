
function main(points=1000)
    if points < 10
        println(points, " grid point(s) have be specified.")
        println("This can not run with less than 10 grid points!")
        println("It is also recommended to use at least 100 grid points.")
        throw(DomainError(points, "points argument must be ≥ 10"))
    elseif points < 100
        println("It is recommended to use at least 100 grid points.")
        println(points, " points have been specified here.")
    end

    xstart1 = 0; xstart2 = 1.5
    ystart1 = 4; ystart2 = 0
    etop1 = 5; etop2 = 4
    eright1 = 4; eright2 = 2.5

    # Assign points to each solid element
    x1points, y1points = pointsdiff(points, xstart1, eright1, ystart1, etop1)
    x2points, y2points = pointsdiff(points, xstart2, eright2, ystart2, etop2)

    println(Int(x1points * y1points), " grid points used.")

    # Remove top row from lower element
    y2range = LinRange(xstart2, eright2, Int(x2points))
    reinterpret(Float64, y2range)
    y2pen = y2range[end] - y2range[end - 1]

    # Obtain mesh arrays for solid element sides
    bot, top, left, right = elements(Int(xpoints1), Int(ypoints1), xstart1,
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
        markersize = 1, markershape = :x, markercolor = :black,
        xlabel = L"$x_1$",
        ylabel = L"$x_2$",
        title = "Meshgrid"
    )

    plot!(meshgrid,
          xyz2[:,1], xyz2[:,2],
          seriestype = :scatter,
          markersize = 1, markershape = :circle, markercolor = :black
          )

    display(meshgrid)
end
