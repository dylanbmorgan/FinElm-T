#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# FE mesh generation function
def mesh(bot,top,left,right):

    # returns [XYZ, CON, DOF] 
    # XYZ - array of nodal coordinates [number of elements x 2]
    # CON - array of node numbers for elements (linear QUADS [number of elements x 4])
    # DOF - array of element DOFs (4-node (linear) quadrilateral element => [number of elements x 8]); 
    
    # number of nodes and elements in the domain
    nnodesx = len(bot)                     # number of horizontal nodes
    nnodesy = len(left)                    # number of vertical nodes 
    nelx = nnodesx-1                       # number of horizontal elements
    nely = nnodesy-1                       # number of vertical elements
    nnodes = nnodesx*nnodesy               # total number of nodes

    # dimensions of the domain
    lx = bot[nnodesx-1] - bot[0]           # length of the domain in x-direction (horizontal)
    ly = left[nnodesy-1] - left[0]         # length of the domain in y-direction (vertical)

    # GENERATE COORDINATES OF NODES 'XYZ'
    XYZ = np.zeros((nnodes,2))            # two-column array [nnodes x 2] containing all nodal coordinates  
    for i in range(nnodesy):              # loop over all nodes on the vertical sides 
        yl = left[i] - left[0]            # distance between node 'i' and left-bottom node '1'
        dy = right[i] - left[i]           # distance between the corresponing nodes j on top and bottom

        for j in range(nnodesx):          # loop over all nodes on the horizontal sides
            xb = bot[j] - bot[0]          # distance between node 'j' and bottom-left node '1' 
            dx = top[j] - bot[j]          # distance between nodes 'j' on opposite sides (top and bottom)

            x = (dx*yl+xb*ly)/(ly-dx*dy/lx) # x-coordinate (horizontal) of a node in the interior of the domain
            y = dy/lx*x+yl                  # y-coordinate (vertical) of a node in the interior of the domain

            XYZ[j+i*nnodesx, 0] = x + bot[0]  # coordinate 'x' in the global coordinate system 
            XYZ[j+i*nnodesx, 1] = y + left[0] # coordinate 'y' in the global coordinate system

    # NODE NUMBERS FOR ELEMENTS 
    nel = nelx*nely                              # total number of elements in the domain
    CON = np.zeros((nel,4), dtype=int)           # [nel*4] array of node number for each element
    for i in range(nely):                        # loop over elements in the vertical direction 
        for j in range(nelx):                    # loop over elements in the horizontal direction 
            # element 'el' and corresponding node numbers
            CON[j+i*nelx, :] = [j+i*nnodesx, j+i*nnodesx+1,j+(i+1)*nnodesx+1, j+(i+1)*nnodesx] 

    # Global DOF for each element (4-node (linear) quadrilateral element)
    DOF = np.zeros((nel,2*4), dtype=int)
    for i in range(nel):
        # defines single row of DOF for each element 'i'
        #print(CON[i, 1]*2+1)
        DOF[i,:] = [CON[i,0]*2, CON[i,1]*2-1, CON[i,1]*2, CON[i,1]*2+1,CON[i,2]*2, CON[i,2]*2+1, CON[i,3]*2, CON[i,3]*2+1]
        #print(DOF[i,:])

    return XYZ, CON, DOF


# Mesh generation procedure
# coordinates defining relevant external edges of the model  

# 8 elements 
bot = [0, 0.5, 1, 1.5, 2]                     # x-coordinates of bottom side nodes
top = [0, 0.5, 1, 1.5, 2]                     # x-coordinates of top side nodes
left = [0, 0.5, 1]                            # y-coordinates of left-hand side nodes
right = [0.5, 0.75, 1]                        # y-coordinates nodes of right-hand side nodes

#generate mesh
XYZ, CON, DOF = mesh(bot,top,left,right)

# plot the mesh 
plt.plot(XYZ[:, 0], XYZ[:, 1], 'sk')

for i in range(len(CON)):
    plt.fill(XYZ[CON[i, :], 0], XYZ[CON[i, :], 1], edgecolor='k', fill=False)

for i in range(4):                             #loop over all nodes within an element
    for j in range(len(CON)):                  #loop over all elements
        sh=0.01
        plt.text(XYZ[CON[j,i],0]+sh,XYZ[CON[j,i],1]+sh, CON[j,i])

# Set chart title.
plt.title("Mesh", fontsize=19)
# Set x axis label.
plt.xlabel("$x_1$", fontsize=10)
# Set y axis label.
plt.ylabel("$x_2$", fontsize=10)
# Set size of tick labels.
plt.tick_params(axis='both', which='major', labelsize=9)

plt.show()
