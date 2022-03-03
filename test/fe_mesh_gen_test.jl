include("../src/fe_mesh_gen.jl")

using MeshGen
using Plots

# 8 elements
bot = [0, 0.5, 1, 1.5, 2]
top = [0, 0.5, 1, 1.5, 2]
left = [0, 0.5, 1]
right = [0.5, 0.75, 1]

#generate mesh
xyz, con, dof = mesh(bot,top,left,right)

display(plot(
    xyz[:, 1], xyz[:, 2]
))
