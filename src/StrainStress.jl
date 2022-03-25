#######################################
#--------- Strain and Stress ---------#
#######################################

include("Deformation.jl")


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
