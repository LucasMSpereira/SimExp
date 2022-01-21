module num

using Printf, LinearAlgebra, CUDA, InOut, SuiteSparse, SparseArrays, Statistics
using Base: print_array
import mtx

### This module provides numerical routines

function sample(element, nip)

    # Returns the local coordinates and weighting coefficients
    # for numerical integration of an element (gauss quadrature)

    if element == "quadrilateral"

        if nip == 1

            points = [0 0]
        
            weights = [4]

        elseif nip == 4

            points = [sqrt(1/3) sqrt(1/3);
                sqrt(1/3) -sqrt(1/3);
                -sqrt(1/3) sqrt(1/3);
                -sqrt(1/3) -sqrt(1/3)]
            
            weights = [1 1 1 1]

        elseif nip == 9

            points = [sqrt(3/5) sqrt(3/5);
                sqrt(3/5) -sqrt(3/5);
                -sqrt(3/5) sqrt(3/5);
                -sqrt(3/5) -sqrt(3/5);
                sqrt(3/5) 0;
                -sqrt(3/5) 0;
                0 sqrt(3/5);
                0 -sqrt(3/5);
                0 0]
            
            weights = [25/81; 40/81; 40/81; 61/81]

        end

    elseif element == "triangle"

        if nip == 1

            points = [1/3 1/3]

            weights = [1/2]

        elseif nip == 3

            points = [1/2 1/2;
            1/2 0;
            0 1/2]

            weights = [1/6; 1/6; 1/6]

        end

    elseif element == "HEXA20"

                points = [0.577350  0.577350  0.577350;
        -0.577350 0.577350  0.577350;
        0.577350 -0.577350  0.577350;
        -0.577350 -0.577350  0.577350;
        0.577350 0.577350  -0.577350;
        -0.577350 0.577350  -0.577350;
        0.577350 -0.577350  -0.577350;
        -0.577350 -0.577350  -0.577350]

        weights = [1; 1; 1; 1; 1; 1; 1; 1]

    end

    return points, weights
        
end

function rect_km(deltaX, deltaY,e, v)
    # This subroutine forms the "analytical" stiffness matrix for
    # rectangular 4- or 8-node plane strain elements using nip=4.
    km = Array{Float64}(undef, (8, 8))
    t2=-1+2*v
    t3=e*t2
    t4=deltaX^2
    t6=-1+v
    t7=2*e*t6
    t8=deltaY^2
    t10=t3*t4+t7*t8
    t11=1/deltaX
    t15=1/(1+v)
    t17=1/t2
    t18=1/deltaY*t15*t17
    t20=t10*t11*t18/6
    t23=e*t15*t17/8
    t25=-e*t2*t4
    t31=(t25+e*t6*t8)*t11*t18/6
    t37=e*(4*v-1)*t15*t17/8
    t40=-t10*t11*t18/12
    t46=(-t25-4*e*t6*t8)*t11*t18/12
    t48=t3*t8
    t49=t7*t4+t48
    t52=t49*t11*t18/6
    t58=(-4*e*t6*t4+t48)*t11*t18/12
    t61=-t49*t11*t18/12
    t67=(e*t6*t4-t48)*t11*t18/6
    km[1,1]=t20
    km[1,2]=-t23
    km[1,3]=t31
    km[1,4]=t37
    km[1,5]=t40
    km[1,6]=t23
    km[1,7]=t46
    km[1,8]=-t37
    km[2,2]=t52
    km[2,3]=-t37
    km[2,4]=t58
    km[2,5]=t23
    km[2,6]=t61
    km[2,7]=t37
    km[2,8]=t67
    km[3,3]=t20
    km[3,4]=t23
    km[3,5]=t46
    km[3,6]=t37
    km[3,7]=t40
    km[3,8]=-t23
    km[4,4]=t52
    km[4,5]=-t37
    km[4,6]=t67
    km[4,7]=-t23
    km[4,8]=t61
    km[5,5]=t20
    km[5,6]=-t23
    km[5,7]=t31
    km[5,8]=t37
    km[6,6]=t52
    km[6,7]=-t37
    km[6,8]=t58
    km[7,7]=t20
    km[7,8]=t23
    km[8,8]=t52
    for i=1:size(km, 1)
        for j=i+1:size(km, 1)
            km[j,i]=km[i,j]
        end
    end
    return km
end

function deemat(e, v, dim, type)

    # Calculates stress-strain matrix dee

    if dim == 2 && type == "strain"

        dee = e*(1 - v)/((1 + v)*(1 - 2 * v))*
        [1 v/(1 - v) 0;
        v/(1 - v) 1 0;
        0 0 (1 - 2*v)/(2*(1 - v))]

    elseif dim == 2 && type == "stress"
        dee = e/(1-v^2)*[
            1 v 0
            v 1 0
            0 0 (1-v)/2
        ]
    elseif dim == 2 && type == "axisymmetric"

        dee = e*(1 - v)/((1 + v)*(1 - 2 * v))*
        [1 v/(1 - v) 0 v/(1 - v);
        v/(1 - v) 1 0 v/(1 - v);
        0 0 (1 - 2*v)/(2*(1 - v)) 0;
        v/(1 - v) v/(1 - v) 0 1]

    elseif dim == 3

        dee = e * (1 - v) / ((1 + v) * (1 - 2 * v)) *
        [1 v/(1 - v) v/(1 - v) 0 0 0;
        v/(1 - v) 1 v/(1 - v) 0 0 0;
        v/(1 - v) v/(1 - v) 1 0 0 0;
        0 0 0 (1 - 2*v)/(2*(1 - v)) 0 0;
        0 0 0 0 (1 - 2*v)/(2*(1 - v)) 0;
        0 0 0 0 0 (1 - 2*v)/(2*(1 - v))]

    end

    return dee
        
end

function shapeFun(points, nip, nod, element)

    # Return shape functions at the current integrating point
    # Ref.: "Introduction to Finite Element Analysis Using MATLAB® and Abaqus", CRC, pg 207-210

    if element == "triangle"

        isoX, isoY = points[nip, :]
        
        if nod == 3
            
            h1 = isoX
            h2 = 1 - isoX - isoY
            h3 = isoY

            fun = [h3; h2; h1]

        elseif nod == 6

            

        end

    elseif element == "quadrilateral"

        isoX, isoY = points[nip, :]

        if nod == 4

            h1 = 0.25*(1 - isoX)*(1 - isoY)
            h2 = 0.25*(1 + isoX)*(1 - isoY)
            h3 = 0.25*(1 + isoX)*(1 + isoY)
            h4 = 0.25*(1 - isoX)*(1 + isoY)

            fun = [h1; h2; h3; h4]

        elseif nod == 5

            h1 = 0.25*(1 + isoX)*(1 + isoY)
            h2 = 0.25*(1 - isoX)*(1 + isoY)
            h3 = 0.25*(1 - isoX)*(1 - isoY)
            h4 = 0.25*(1 + isoX)*(1 + isoY)
            h5 = 0.5*(1 - isoX^2)*(1 + isoY)

            fun = [h1 - 0.5*h5; h2 - 0.5*h5; h3; h4; h5]

        elseif nod == 6

            h1 = 0.25*(1 + isoX)*(1 + isoY)
            h2 = 0.25*(1 - isoX)*(1 + isoY)
            h3 = 0.25*(1 - isoX)*(1 - isoY)
            h4 = 0.25*(1 + isoX)*(1 + isoY)
            h5 = 0.5*(1 - isoX^2)*(1 + isoY)
            h6 = 0.5*(1 - isoY^2)*(1 - isoX)

            fun = [h1 - 0.5*h5; h2 - 0.5*h5 - 0.5*h6; h3 - 0.5*h6; h4; h5; h6]

        elseif nod == 7

            h1 = 0.25*(1 + isoX)*(1 + isoY)
            h2 = 0.25*(1 - isoX)*(1 + isoY)
            h3 = 0.25*(1 - isoX)*(1 - isoY)
            h4 = 0.25*(1 + isoX)*(1 + isoY)
            h5 = 0.5*(1 - isoX^2)*(1 + isoY)
            h6 = 0.5*(1 - isoY^2)*(1 - isoX)
            h7 = 0.5*(1 - isoX^2)*(1 - isoY)

            fun = [h1 - 0.5*h5; h2 - 0.5*h5 - 0.5*h6; h3 - 0.5*h6 - 0.5*h7; h4 - 0.5*h7; h5; h6; h7]

        elseif nod == 8

            h1 = 0.25*(1 + isoX)*(1 + isoY)
            h2 = 0.25*(1 - isoX)*(1 + isoY)
            h3 = 0.25*(1 - isoX)*(1 - isoY)
            h4 = 0.25*(1 + isoX)*(1 + isoY)
            h5 = 0.5*(1 - isoX^2)*(1 + isoY)
            h6 = 0.5*(1 - isoY^2)*(1 - isoX)
            h7 = 0.5*(1 - isoX^2)*(1 - isoY)
            h8 = 0.5*(1 - isoY^2)*(1 + isoX)

            fun = [h1 - 0.5*h5 - 0.5*h8; h2 - 0.5*h5 - 0.5*h6; h3 - 0.5*h6 - 0.5*h7;
            h4 - 0.5*h7 - 0.5*h8; h5; h6; h7; h8]

        elseif nod == 9

            h1 = 0.25*(1 + isoX)*(1 + isoY)
            h2 = 0.25*(1 - isoX)*(1 + isoY)
            h3 = 0.25*(1 - isoX)*(1 - isoY)
            h4 = 0.25*(1 + isoX)*(1 + isoY)
            h5 = 0.5*(1 - isoX^2)*(1 + isoY)
            h6 = 0.5*(1 - isoY^2)*(1 - isoX)
            h7 = 0.5*(1 - isoX^2)*(1 - isoY)
            h8 = 0.5*(1 - isoY^2)*(1 + isoX)
            h9 = (1 - isoX^2)*(1 - isoY^2)

            fun = [h1 - 0.5*h5 - 0.5*h8 - 0.25*h9; h2 - 0.5*h5 - 0.5*h6 - 0.25*h9; h3 - 0.5*h6 - 0.5*h7 - 0.25*h9;
            h4 - 0.5*h7 - 0.5*h8 - 0.25*h9; h5 - 0.5*h9; h6 - 0.5*h9; h7 - 0.5*h9; h8 - 0.5*h9; h9]
            
        end

    elseif element == "HEXA20"

        # page 617 of the FEA book

        isoX, isoY, isoZ = points[nip, :]

        fun = [
        (1 - isoX)*(1 - isoY)*(1 - isoZ)*(-isoX - isoY - isoZ - 2)/8;
        (1 - isoX)*(1 - isoY)*(1 - isoZ^2)/4;
        (1 - isoX)*(1 - isoY)*(1 + isoZ)*(-isoX - isoY + isoZ - 2)/8;
        (1 - isoX^2)*(1 - isoY)*(1 + isoZ)/4;
        (1 + isoX)*(1 - isoY)*(1 + isoZ)*(isoX - isoY + isoZ - 2)/8;
        (1 + isoX)*(1 - isoY)*(1 - isoZ^2)/4;
        (1 + isoX)*(1 - isoY)*(1 - isoZ)*(isoX - isoY - isoZ - 2)/8;
        (1 - isoX^2)*(1 - isoY)*(1 - isoZ)/4;
        (1 - isoX)*(1 - isoY^2)*(1 - isoZ)/4;
        (1 - isoX)*(1 - isoY^2)*(1 + isoZ)/4;
        (1 + isoX)*(1 - isoY^2)*(1 + isoZ)/4;
        (1 + isoX)*(1 - isoY^2)*(1 - isoZ)/4;
        (1 - isoX)*(1 + isoY)*(1 - isoZ)*(-isoX + isoY - isoZ - 2)/8;
        (1 - isoX)*(1 + isoY)*(1 - isoZ^2)/4;
        (1 - isoX)*(1 + isoY)*(1 + isoZ)*(-isoX + isoY + isoZ - 2)/8;
        (1 - isoX^2)*(1 + isoY)*(1 + isoZ)/4;
        (1 + isoX)*(1 + isoY)*(1 + isoZ)*(isoX + isoY + isoZ - 2)/8;
        (1 + isoX)*(1 + isoY)*(1 - isoZ^2)/4;
        (1 + isoX)*(1 + isoY)*(1 - isoZ)*(isoX + isoY - isoZ - 2)/8;
        (1 - isoX^2)*(1 + isoY)*(1 - isoZ)/4]

    end

    return fun

end

function shapeDer(points, nip, nod, element)

    # Returns ('n° of DOFs/node' x 'n° of nodes/element')  matrix with shape
    # functions derivatives in isoparametric coordinates

    if element == "triangle"

        isoX, isoY = points[nip, :]

        if nod == 3

            h1 = [1; 0]
            h2 = [-1; -1]
            h3 = [0; 1]

            der = hcat(h3, h2, h1)

        elseif nod == 4

            h1 = [0.25*(1 + isoY); 0.25*(1 + isoX)]
            h2 = [-0.25*(1 + isoY); 0.25*(1 - isoX)]
            h3 = [-0.25*(1 + isoY); -0.25*(1 + isoX)]
            h4 = [0.25*(1 + isoY); 0.25*(1 + isoX)]

            der = hcat(h1, h2, h3, h4)

        elseif nod == 5

            h1 = [0.25*(1 + isoY); 0.25*(1 + isoX)]
            h2 = [-0.25*(1 + isoY); 0.25*(1 - isoX)]
            h3 = [-0.25*(1 + isoY); -0.25*(1 + isoX)]
            h4 = [0.25*(1 + isoY); 0.25*(1 + isoX)]
            h5 = [-isoX*(1 + isoY); 0.5*(1 - isoX^2)]

            der = hcat(h1 - h5.*0.5, h2 - h5.*0.5, h3, h4, h5)

        elseif nod == 6

            h1 = [0.25*(1 + isoY); 0.25*(1 + isoX)]
            h2 = [-0.25*(1 + isoY); 0.25*(1 - isoX)]
            h3 = [-0.25*(1 + isoY); -0.25*(1 + isoX)]
            h4 = [0.25*(1 + isoY); 0.25*(1 + isoX)]
            h5 = [-isoX*(1 + isoY); 0.5*(1 - isoX^2)]
            h6 = [-0.5*(1 - isoY^2); -isoY*(1 - isoX)]

            der = hcat(h1 - h5.*0.5, h2 - h5.*0.5 - h6.*0.5, h3 - h6.*0.5, h4, h5, h6)

        elseif nod == 7

            h1 = [0.25*(1 + isoY); 0.25*(1 + isoX)]
            h2 = [-0.25*(1 + isoY); 0.25*(1 - isoX)]
            h3 = [-0.25*(1 + isoY); -0.25*(1 + isoX)]
            h4 = [0.25*(1 + isoY); 0.25*(1 + isoX)]
            h5 = [-isoX*(1 + isoY); 0.5*(1 - isoX^2)]
            h6 = [-0.5*(1 - isoY^2); -isoY*(1 - isoX)]
            h7 = [-isoX*(1 - isoY); -0.5*(1 - isoX^2)]

            der = hcat(h1 - h5.*0.5, h2 - h5.*0.5 - h6.*0.5, h3 - h6.*0.5 - h7.*0.5, h4 - h7.*0.5, h5, h6, h7)

        elseif nod == 8

            h1 = [0.25*(1 + isoY); 0.25*(1 + isoX)]
            h2 = [-0.25*(1 + isoY); 0.25*(1 - isoX)]
            h3 = [-0.25*(1 + isoY); -0.25*(1 + isoX)]
            h4 = [0.25*(1 + isoY); 0.25*(1 + isoX)]
            h5 = [-isoX*(1 + isoY); 0.5*(1 - isoX^2)]
            h6 = [-0.5*(1 - isoY^2); -isoY*(1 - isoX)]
            h7 = [-isoX*(1 - isoY); -0.5*(1 - isoX^2)]
            h8 = [0.5*(1 - isoY^2); -isoY*(1 + isoX)]

            der = hcat(h1 - h5.*0.5 - h8.*0.5, h2 - h5.*0.5 - h6.*0.5, h3 - h6.*0.5 - h7.*0.5,
            h4 - h7.*0.5 - h8.*0.5, h5, h6, h7, h8)

        elseif nod == 9

            h1 = [0.25*(1 + isoY); 0.25*(1 + isoX)]
            h2 = [-0.25*(1 + isoY); 0.25*(1 - isoX)]
            h3 = [-0.25*(1 + isoY); -0.25*(1 + isoX)]
            h4 = [0.25*(1 + isoY); 0.25*(1 + isoX)]
            h5 = [-isoX*(1 + isoY); 0.5*(1 - isoX^2)]
            h6 = [-0.5*(1 - isoY^2); -isoY*(1 - isoX)]
            h7 = [-isoX*(1 - isoY); -0.5*(1 - isoX^2)]
            h8 = [0.5*(1 - isoY^2); -isoY*(1 + isoX)]
            h9 = [-2*isoX*(1 - isoY^2); -2*isoY*(1 - isoX^2)]

            der = hcat(h1 - h5.*0.5 - h8.*0.5 - h9.*0.25, h2 - h5.*0.5 - h6.*0.5 - h9.*0.25, h3 - h6.*0.5 - h7.*0.5 - h9.*0.25,
            h4 - h7.*0.5 - h8.*0.5 - h9.*0.25, h5 - h9.*0.5, h6 - h9.*0.5, h7 - h9.*0.5, h8 - h9.*0.5, h9)
            
        end

    elseif element == "quadrilateral"
        isoX, isoY = points[nip, :]
        h1 = [-0.25*(1 - isoY); -0.25*(1 - isoX)]
        h2 = [0.25*(1 - isoY); -0.25*(1 + isoX)]
        h3 = [0.25*(1 + isoY); 0.25*(1 + isoX)]
        h4 = [-0.25*(1 + isoY); 0.25*(1 - isoX)]
        der = hcat(h1, h2, h3, h4)

    elseif element == "HEXA20"

        isoX, isoY, isoZ = points[nip, :]

        der = [
        (2*isoX + isoY + isoZ + 1)*(1 - isoY)*(1 - isoZ)/8	(isoX + 2*isoY + isoZ + 1)*(1 - isoX)*(1 - isoZ)/8	(isoX + isoY + 2*isoZ + 1)*(1 - isoX)*(1 - isoY)/8;
        -(1 - isoY)*(1 - isoZ^2)/4	-(1 - isoX)*(1 - isoZ^2)/4	-isoZ*(1 - isoX)*(1 - isoY)/2;
        (2*isoX + isoY - isoZ + 1)*(1 - isoY)*(1 + isoZ)/8	(isoX + 2*isoY - isoZ + 1)*(1 - isoX)*(1 + isoZ)/8	(-isoX - isoY + 2*isoZ - 1)*(1 - isoX)*(1 - isoY)/8;
        -isoX*(1 - isoY)*(1 + isoZ)/2	-(1 - isoX^2)*(1 + isoZ)/4	(1 - isoX^2)*(1 - isoY)/4;
        (2*isoX - isoY + isoZ - 1)*(1 - isoY)*(1 + isoZ)/8	(-isoX + 2*isoY - isoZ + 1)*(1 + isoX)*(1 + isoZ)/8	(isoX - isoY + 2*isoZ - 1)*(1 + isoX)*(1 - isoY)/8;
        (1 - isoY)*(1 - isoZ^2)/4	-(1 + isoX)*(1 - isoZ^2)/4	-isoZ*(1 + isoX)*(1 - isoY)/2;
        (2*isoX - isoY - isoZ - 1)*(1 - isoY)*(1 - isoZ)/8	(-isoX + 2*isoY + isoZ + 1)*(1 + isoX)*(1 - isoZ)/8	(-isoX + isoY + 2*isoZ + 1)*(1 + isoX)*(1 - isoY)/8;
        -isoX*(1 - isoY)*(1 - isoZ)/2	-(1 - isoX^2)*(1 - isoZ)/4	-(1 - isoX^2)*(1 - isoY)/4;
        -(1 - isoY^2)*(1 - isoZ)/4	-isoY*(1 - isoX)*(1 - isoZ)/2	-(1 - isoX)*(1 - isoY^2)/4;
        -(1 - isoY^2)*(1 + isoZ)/4	-isoY*(1 - isoX)*(1 + isoZ)/2	(1 - isoX)*(1 - isoY^2)/4;
        (1 - isoY^2)*(1 + isoZ)/4	-isoY*(1 + isoX)*(1 + isoZ)/2	(1 + isoX)*(1 - isoY^2)/4;
        (1 - isoY^2)*(1 - isoZ)/4	-isoY*(1 + isoX)*(1 - isoZ)/2	-(1 + isoX)*(1 - isoY^2)/4;
        (2*isoX - isoY + isoZ + 1)*(1 + isoY)*(1 - isoZ)/8	(-isoX + 2*isoY - isoZ - 1)*(1 - isoX)*(1 - isoZ)/8	(isoX - isoY + 2*isoZ + 1)*(1 - isoX)*(1 + isoY)/8;
        -(1 + isoY)*(1 - isoZ^2)/4	(1 - isoX)*(1 - isoZ^2)/4	-isoZ*(1 - isoX)*(1 + isoY)/2;
        (2*isoX - isoY - isoZ + 1)*(1 + isoY)*(1 + isoZ)/8	(-isoX + 2*isoY + isoZ - 1)*(1 - isoX)*(1 + isoZ)/8	(-isoX + isoY + 2*isoZ - 1)*(1 - isoX)*(1 + isoY)/8;
        -isoX*(1 + isoY)*(1 + isoZ)/2	(1 - isoX^2)*(1 + isoZ)/4	(1 - isoX^2)*(1 + isoY)/4;
        (2*isoX + isoY+ isoZ - 1)*(1 + isoY)*(1 + isoZ)/8	(isoX + 2*isoY+ isoZ - 1)*(1 + isoX)*(1 + isoZ)/8	(isoX + isoY+ 2*isoZ - 1)*(1 + isoX)*(1 + isoY)/8;
        (1 + isoY)*(1 - isoZ^2)/4	(1 + isoX)*(1 - isoZ^2)/4	-isoZ*(1 + isoX)*(1 + isoY)/2;
        (2*isoX + isoY - isoZ - 1)*(1 + isoY)*(1 - isoZ)/8	(isoX + 2*isoY - isoZ - 1)*(1 + isoX)*(1 - isoZ)/8	(-isoX - isoY + 2* isoZ + 1)*(1 + isoX)*(1 + isoY)/8;
        -isoX*(1 + isoY)*(1 - isoZ)/2	(1 - isoX^2)*(1 - isoZ)/4	-(1 - isoX^2)*(1 + isoY)/4
        ]'

    end

    return der
    
end

function beemat(der, coord, type_2d, nodof, nod)

    # Returns the derivatives of the shape functions in the cartesian
    # system (bee) and the determinant of the jacobian

    jac = der*coord
    determ = det(jac)
    deriv = inv(jac)*der

    if nodof == 2

        if type_2d == "stress" || type_2d == "strain"

            bee = fill(0.0, (3, nodof*nod))
            for m in 1:nod

                k = 2*m
                l = k - 1
                global bee[1, l] = deriv[1, m]
                global bee[3, k] = deriv[1, m]
                global bee[2, k] = deriv[2, m]
                global bee[3, l] = deriv[2, m]

            end
            
        elseif type_2d == "axisymmetric"

            bee = fill(0.0, (4, nodof*nod))
            for m in 1:nod
                
                k = 2*m
                l = k - 1
                global bee[1, l] = deriv[1, m]
                global bee[3, k] = deriv[1, m]
                global bee[2, k] = deriv[2, m]
                global bee[3, l] = deriv[2, m]

            end

        end

    elseif nodof == 3

        bee = fill(0.0, (6, nodof*nod))
        for m in 1:nod

            n = 3*m
            k = n - 1
            l = k - 1
            bee[1, l] = deriv[1, m]
            bee[4, k] = deriv[1, m]
            bee[6, n] = deriv[1, m]
            bee[2, k] = deriv[2, m]
            bee[4, l] = deriv[2, m]
            bee[5, n] = deriv[2, m]
            bee[3, n] = deriv[3, m]
            bee[5, k] = deriv[3, m]
            bee[6, l] = deriv[3, m]

        end

    end

    return bee, determ
    
end

function para_beemat!(bee_d, deriv_d)

    # Returns the derivatives of the shape functions in the cartesian
    # system (bee).
    # Parallelized version of beemat()

    if size(bee_d, 1) == 3

        m = threadIdx().x
        k = 2*m
        l = k - 1
        bee_d[1, l] = deriv_d[1, m]
        bee_d[3, k] = deriv_d[1, m]
        bee_d[2, k] = deriv_d[2, m]
        bee_d[3, l] = deriv_d[2, m]

    elseif size(bee_d, 1) == 4

        m = threadIdx().x
        k = 2*m
        l = k - 1
        bee_d[1, l] = deriv_d[1, m]
        bee_d[3, k] = deriv_d[1, m]
        bee_d[2, k] = deriv_d[2, m]
        bee_d[3, l] = deriv_d[2, m]

    elseif size(bee_d, 1) == 6

        m = threadIdx().x
        n = 3*m
        k = n - 1
        l = k - 1
        bee_d[1, l] = deriv_d[1, m]
        bee_d[4, k] = deriv_d[1, m]
        bee_d[6, n] = deriv_d[1, m]
        bee_d[2, k] = deriv_d[2, m]
        bee_d[4, l] = deriv_d[2, m]
        bee_d[5, n] = deriv_d[2, m]
        bee_d[3, n] = deriv_d[3, m]
        bee_d[5, k] = deriv_d[3, m]
        bee_d[6, l] = deriv_d[3, m]

    end
    
    return
    
end

function para_fsparv!(preKv_indices_d, preKv_d, km, g, kdiag, intermPos)

    # Assembles element matrices into a symmetric
    # global matrix stored in skyline form.
    for i in 1:CUDA.size(g, 1)

        if g[i, threadIdx().x] != 0

            for j in 1:CUDA.size(g, 1)

                if g[j, threadIdx().x] != 0 && (g[i, threadIdx().x] - g[j, threadIdx().x]) >= 0

                    intermPos[threadIdx().x] += 1
                    # Store position to be altered in kv_d
                    preKv_indices_d[intermPos[threadIdx().x], threadIdx().x] = kdiag[g[i, threadIdx().x]] - (g[i, threadIdx().x] - g[j, threadIdx().x])

                    # Store value to be incremented in the respective position of kv_d
                    preKv_d[intermPos[threadIdx().x], threadIdx().x] = km[i, j, threadIdx().x]

                end

            end

        end

    end
    
    return
    
end

function fsparv!(kv, km, g, kdiag)

    # Assembles element matrices into a symmetric global
    # matrix stored in skyline form.

    for i in 1:size(g, 1)
        
        if g[i] != 0

            for j in 1:size(g, 1)
                
                if g[j] != 0 && (g[i] - g[j]) >= 0
                    
                    kv[kdiag[g[i]] - (g[i] - g[j])] += km[i, j]
    
                end

            end

        end

    end

    return kv
    
end

function sparin!(kv, kdiag)

    # Returns the (Choleski) factorised vector kv stored as a skyline

    n = size(kdiag, 1)
    kv[1] = sqrt(kv[1])
    for i = 2:n
        # a partir do segundo elemento de kdiag (quase 4000)

        ki = kdiag[i] - i
        l = kdiag[i - 1] - ki + 1
        for j = l:i
            # de {kdiag[i - 1] - kdiag[i] - i + 1} à posição atual de kdiag

            global x = kv[ki + j]
            kj = kdiag[j] - j
            if j != 1

                ll = kdiag[j - 1] - kj + 1
                l > ll && (ll = l)
                if ll != j

                    m = j - 1
                    for k = ll:m

                        x = x - kv[ki + k]*kv[kj + k]

                    end

                end

            end
            
            kv[ki + j] = x/kv[kj + j]

        end

        kv[ki + i] = sqrt(x)

    end

    return kv

end

function spabac!(kv, kdiag, disp)

    # Cholesky forward and back substitution on symmetric skyline global matrix

    disp[1] = disp[1]/kv[1]
    for i = 2:size(kdiag, 1)

        ki = kdiag[i] - i
        l = kdiag[i - 1] - ki + 1
        x = disp[i]

        if l != i

            m = i - 1   
            for j = l:m

                x = x - kv[ki + j]*disp[j]

            end

        end
        disp[i] = x/kv[ki + i]

    end
    for it = 2:size(kdiag, 1)

        i = size(kdiag, 1) + 2 - it
        ki = kdiag[i] - i
        x = disp[i]/kv[ki + i]
        disp[i] = x
        l = kdiag[i - 1] - ki + 1
        if l != i

            m = i - 1
            for k = l:m

                disp[k] = disp[k] - x*kv[ki + k]

            end

        end

    end

    disp[1] = disp[1]/kv[1]
    return disp
    
end

function pcg!(par, diag_precon, storkm, output)

    # Preconditioned conjugate gradient (pcg) method for solving FEA system
    # "element-by-element" (without assembly of global stiffness matrix).
    # Section 3.5 of the FEA book
    
    for i in 1:size(diag_precon, 1)

        diag_precon[i] != 0 ? diag_precon[i] = 1/diag_precon[i] : false
    
    end
    cg_iters = 0
    x = zeros(par.neq)
    d = zeros(par.neq)
    xnew = zeros(par.neq)
    x = zeros(par.neq)
    d .= diag_precon.*par.disp[:, 1]
    p = d

    for i in 1:par.cg_limit

        cg_iters += 1
        u = zeros(par.neq)
        for iel in 1:par.nels
            
            km = storkm[:, :, iel]
            m = zeros((size(par.g_g, 1), 1))
            for b in 1:size(par.g_g, 1)

                par.g_g[b, iel] != 0 ? m[b] = p[par.g_g[b, iel]] : false
                
            end
            l = km*m
            for b in 1:size(par.g_g, 1)

                par.g_g[b, iel] != 0 ? u[par.g_g[b, iel]] += l[b] : false
                
            end

        end
        # fixed_freedoms != 0 ? u[no] = p[no]*store : false
        up = par.disp[:, i]'*d
        alpha = up/(p'*u)
        xnew = x + alpha*p
        par.disp[:, i + 1] = par.disp[:, i] - alpha*u
        d .= diag_precon.*par.disp[:, i + 1]
        beta = (par.disp[:, i + 1]'*d)/up
        p = d + beta*p
        println("\n$(i)     $(maximum(abs.(xnew - x))/maximum(abs.(xnew)))\n")
        maximum(abs.(xnew - x))/maximum(abs.(xnew)) <= par.cg_tol ? break : false
        x = xnew

    end

    if (maximum(abs.(xnew - x))/maximum(abs.(xnew))) <= par.cg_tol

        write(output, "\nNumber of PCG iterations to convergence was $(cg_iters).\n\n")
        println("\nNumber of PCG iterations to convergence was $(cg_iters).\n")

    else

        write(output, "\n *** WARNING: PCG didn't convergence in $(par.cg_limit) iterations. ***\n\n")
        println("\n *** WARNING: PCG didn't convergence in $(par.cg_limit) iterations. ***\n\n")
        
    end

    return x
    
end

function formku!(ku, km, g)
    # Assembles element stiffness matrices into symmetrical
    # global matrix (stored as an upper rectangle).
    for i = 1:size(km, 1)
        if g[i] != 0
            for j = 1:size(km, 1)
                if g[j] != 0
                    icd = g[j]-g[i]+1
                    icd >= 1 && (ku[g[i], icd] += km[i,j])
                end
            end
        end
    end
    return ku
end

function elmat(area, density, element)
    # Form lumped mass matrix for quadrilateral
    # 4- or 8-node plane strain elements.
    if element == "quadrilateral"
        mm = 0.25*area*density
    else
        println("num.elmat: Invalid type of element.")
    end
    return mm
end

function StandardEigen!(diag, nband, ku)
    # Reduce to standard eigenvalue problem
    rrmass = Array{Float64}(undef, length(diag))
    if minimum(diag) <= 0
        println("num.StandardEigen: vector diag (global mass matrix) has nonpositive component.")
    else
        for i = 1:length(diag)
            rrmass[i] = 1/sqrt(diag[i])
        end
    end
    for i = 1:length(diag)
        if i <= length(diag)-nband
            k = nband+1
        else
            k=length(diag)-i+1
        end
        for j = 1:k
            ku[i,j] *= rrmass[i]*rrmass[i+j-1]
        end
    end
    return rrmass, ku
end

function bisect!(diag,udiag,etol, modes, output)
    # This subroutine finds the eigenvalues of a tridiagonal matrix,
    # given with its diagonal elements in the array diag[n] and
    # its subdiagonal elements in the last n - 1 stores of the
    # array udiag[n], using ql transformations. The eigenvalues are
    # overwritten on the diagonal elements in the array diag in
    # ascending order. The subroutine will fail if any one
    # eigenvalue takes more than 30 iterations.
    if size(diag,1)!=1
        for i=2:size(diag,1)
            udiag[i-1]=udiag[i]
        end
    end
    udiag[size(diag,1)]=0
    b=0
    f=0
    for l=1:size(diag,1)
        j=0
        h=etol*(abs(diag[l])+abs(udiag[l]))
        (b<h) && (b=h)
        # look for small subdiagonal element
        for m=l:size(diag,1)
            global lastM = m
            (abs(udiag[m])<=b) && (break)
        end
        if l!=lastM
            while true
                if j==30
                    ifail=1
                    break
                end
                j += 1
                # form shift
                g = diag[l]
                h = diag[l+1] - g
                if abs(h)<abs(udiag[l])
                    p=h*0.5/udiag[l]
                    r=sqrt(p^2 + 1)
                    h=p+r
                    (p<0) && (h=p-r)
                    diag[l]=udiag[l]/h
                else
                    p=2*udiag[l]/h
                    r=sqrt(p^2 + 1)
                    diag[l]=udiag[l]*p/(1+r)
                end
                h=g-diag[l]
                i1=l+1
                if i1<=size(diag,1)
                    for i=i1:size(diag,1)
                        diag[i]-=h
                    end
                end
                f+=h
                # ql transformation
                p=diag[lastM]
                c=1
                s=0
                m1=lastM-1
                for ii=l:m1
                    i=m1-ii+l
                    g=c*udiag[i]
                    h=c*p
                    if abs(p)>=abs(udiag[i])
                        c= udiag[i]/p
                        r=sqrt(c^2 + 1)
                        udiag[i+1]=s*p*r
                        s=c/r
                        c=1/r
                    else
                        c=p/udiag[i]
                        r=sqrt(c^2 + 1)
                        udiag[i+1]=s*udiag[i]*r
                        s=1/r
                        c=c/r
                    end
                    p=c*diag[i]-s*g
                    diag[i+1]=h+s*(c*g+s*diag[i])
                end
                udiag[l]=s*p
                diag[l]=c*p
                abs(udiag[l])<=b && (break)
            end
        end
        p=diag[l]+f
        # order eigenvalue
        aux=0
        if l!=1
            for ii=2:l
                i=l-ii+2
                if p>=diag[i-1]
                    aux=1
                    break
                end
                diag[i]=diag[i-1]
            end
        end
        if aux==0
            i=1
        end
        diag[i] = p
        ifail=0
    end
    # Output eigenvalues
    InOut.EigenValOut(diag, modes, output)
    return diag
end

function sparin_gauss!(kdiag, kv)
    # This subroutine performs Gaussian factorisation of a skyline matrix.
    for j=1:(length(kdiag)-1)
        den=kv[kdiag[j]]
        ii=0                 
        for i = (j+1):length(kdiag)
            ii+=1
            l=kdiag[i]-ii
            if (l-kdiag[i-1])>0
                num=kv[l]
                fac=num/den
                kk=-1
                for k=i:length(kdiag)
                    kk+=1
                    l1=kdiag[i+kk]-kk
                    l2=l1-ii
                    l3=kdiag[i+kk-1]
                    (l2-l3>0) && (kv[l1]-=fac*kv[l2])
                end
            end
        end
    end
    return kv
end

function spabac_gauss!(kv,udiag,kdiag)
    # This subroutine performs Gaussian forward and
    # back-substitution on a skyline matrix
    for j=1:(length(kdiag)-1)
        den=kv[kdiag[j]]
        global ii=0
        for i=(j+1):length(kdiag)
            global ii+=1
            l=kdiag[i]-ii
            if l-kdiag[i-1]>0
                num=kv[l]
                fac=num/den
                udiag[i]-=fac*udiag[j]
            end
        end
    end
    udiag[length(kdiag)]/=kv[kdiag[length(kdiag)]]
    for i=(length(kdiag)-1):-1:1
        jj=0
        asum=0
        for j=(i+1):length(kdiag)
            jj+=1
            l1=kdiag[i+jj]-jj
            l2=kdiag[i+jj-1]
            (l1-l2>0) && (asum+=kv[l1]*udiag[j])
        end
        udiag[i]=(udiag[i]-asum)/kv[kdiag[i]]
    end
    return udiag
end

function CentroidStress(par, disp, output)
    
    if size(par.props, 2) == 1
        
        # Return elastic stress–strain dee matrix in 2D (plane strain) or
        # 3D. Uses Young’s modulus and Poisson’s ratio.
        dee = num.deemat(par.props[1, 1], par.props[2, 1], par.nodof, par.type_2d)
        
    end
    
    #InOut.HeaderStress(par.type_2d, par.nodof, output)
    for iel = 1:par.nels
        
        eld = zeros(size(par.g_g, 1))
        for l in 1:size(par.g_g, 1)
            
            if par.g_g[l, iel] != 0
                
                global eld[l] = disp[par.g_g[l, iel]]
                
            end
            
        end
        
        if size(par.props, 2) > 1
            
            # Return elastic stress–strain dee matrix in 2D (plane strain) or
            # 3D. Uses Young’s modulus and Poisson’s ratio.
            dee2 = num.deemat(par.props[1, [etype[iel]]], par.props[2, [etype[iel]]], par.nodof, par.type_2d)
            
        end
        
        
        # Remaining operations would usually be inside 
        # a loop in the integrating points
        # Shape function derivatives in normalized coordinate system
        der = num.shapeDer([0 0], 1, par.nod, par.element)
        # Shape function derivatives in global coordinate system
        # and jacobian
        bee, determ = num.beemat(der, par.coord[par.g_num[:, iel], :], par.type_2d,
                                par.nodof, par.nod)
        
        if par.type_2d == "axisymmetric"
            
            fun = num.shapeFun(points, i, nod, element)
            global gc = fun*coord
            bee[4, 1:2:nodof*nod-1] = fun[:]/gc[1]
            
        end
        
        if size(par.props, 2) > 1
            
            sigma = dee2*(bee*eld)
            
        else
            
            sigma = dee*(bee*eld)
            
        end
        
        #InOut.StressOut(iel, par.gc, sigma, 1, par.nip, output)
        num.vonMises(sigma, iel, 10, output)

    end
end

function centerPos(element, coord, g_num, iel)
    # Determine geometrical center of element
    if element=="quadrilateral"
        x=mean(coord[g_num[[1 2], iel],1])
        y=mean(coord[g_num[[1 3], iel], 2])
    end
    return x, y
end

function vonMises(sigma, iel, nxe, output)
    vm = sqrt(sigma'*V*sigma)
    InOut.vonMisesOut(vm, iel, nxe, output)
end

function EigenVec(nmodes, kh, kdiag, diag, penalty, udiag, rrmass, output)
    # This function determines the eigenvectors corresponding to
    # the eigenvalues calculated previously
    for i=1:nmodes
        kv=kh
        kv[kdiag].-=diag[i]
        kv[1]+=penalty
        udiag=zeros(Float64, length(udiag))
        udiag[1]=kv[1]
        kv = num.sparin_gauss!(kdiag, kv)
        udiag = num.spabac_gauss!(kv,udiag,kdiag)
        udiag.*=rrmass
        InOut.EigenVecOut(i, udiag, output)
    end
end

function globals!(kmGlobal, mmGlobal, mm, km, g)
    for i in 1:length(g)
        for j in 1:length(g)
            if g[i]*g[j] != 0
                kmGlobal[g[i], g[j]] += km[i, j]
                mmGlobal[g[i], g[j]] += mm[i, j]
            end
        end
    end
    return kmGlobal, mmGlobal
end

function globalKm!(kmGlobal, km, g)
    for i in 1:length(g)
        for j in 1:length(g)
            if g[i]*g[j] != 0
                kmGlobal[g[i], g[j]] += km[i, j]
            end
        end
    end
    return kmGlobal
end

function ecmat(nod, nodof, ndof, fun)
    # Forms consistent element matrix
    nt=zeros(Float64, (ndof, nodof))
    tn=zeros(Float64, (nodof, ndof))
    for i=1:nod 
        for j=1:nodof
            nt[(i-1)*nodof+j,j]=fun[i]
            tn[j,(i-1)*nodof+j]=fun[i]
        end
    end
    return nt*tn
end

function mass!(mm, ecm, rho, det, weights, nip, niptot)
    ecm *= weights*det
    for i in 1:size(ecm,1)
        for j in 1:size(ecm,1)
            mm[i, j] += ecm[i,j]
        end
    end
    nip == niptot && (mm *= rho)
    return mm
end

function simpleEigenProblem(kmGlobal, mmGlobal)
    div = inv(mmGlobal)*kmGlobal
    freqsHz = sqrt.(eigvals(div))./(2*pi)
    println(freqsHz[1:3])
    println(maximum(freqsHz))
    eigVec = eigvecs(div)
    println(size(eigVec))
end

function jacobiRots!(k, m, nmodes, s = 12)
    # Bathe section 11.3
    λold = zeros(Float64, size(k,1))
    𝚽 = I
    crit1 = true
    crit2 = true
    count = 0
    while crit1 || crit2
        count += 1
        println(count)
        P = zeros(Float64, size(k))
        for i in 1:size(k,1)
            for j in 1:size(k, 1)
                if i!=j
                    kBar_ii = -m[i, i]*k[i, j]
                    kBar_jj = -m[j, j]*k[i, j]
                    kBar = k[i, i]*m[j, j] - k[j, j]*m[i, i]
                    sign(kBar) == 0 ? bb = 1 : bb = sign(kBar)
                    x = kBar/2 + bb*sqrt((kBar/2)^2 + kBar_ii*kBar_jj)
                else
                    P[i,j] = 1.0
                end
                if j < i
                    P[i, j] = -kBar_ii/x
                elseif j > i
                    P[i, j] = kBar_jj/x
                end
            end
        end
        # Update matrices k, m and the product of all Pₗ
        println("P = $P")
        k = P'*k*P
        m = P'*m*P
        𝚽 *= P
        # Update eigenvalue criterion
        λnew = diag(k)./diag(m)
        crit1 = reduce(*, abs.(λnew .- λold)./λnew .< 10.0^-s)
        λold = λnew
        # Update k and m criterion
        for i in 1:size(k,1)
            for j in 1:size(k, 1)
                i >= j && continue
                crit2a = (sqrt(k[i, j]^2/(k[i, i]*k[j,j])) < 10.0^-s)
                crit2b = (sqrt(m[i, j]^2/(m[i, i]*m[j,j])) < 10.0^-s)
                crit2 *= crit2a*crit2b
            end
        end
    end
    𝚽 *= inv(sqrt(m))
    # Obtain ordered list of eigenvalues
    vals = zeros(Float64, (2, size(k, 1)))
    vals[1,:] = diag(k/m)
    vals[2, :] = 1:size(k, 1)
    
    vals = sortslices(vals, dims=2)
    # Order eigenvectors accordingly
    vecs = 𝚽[:,convert(Vector{Int32},vals[2,:])]
    return vals[1:nmodes], vecs[:,1:nmodes]
end

function finDiff(newObjective, oldObjective, changeSens)
    # Numerical aproximation of sensitivity
    return (newObjective-oldObjective)/changeSens
end

function compSens(p, dens, para, iel, k, disp)
    # Analytical calculation of compliance sensitivity (without self-weight)
    u = num.eleDisp(para.ndof, disp, para.g_g[:,iel])
    return -0.5*p*dens^(p-1)*u'*k*u
end

function compSensExt(p, dens, para, iel, k, disp)
    # Analytical calculation of compliance sensitivity (without self-weight)
    u = num.eleDisp(para.ndof, disp, para.g_g[:,iel])
    return -0.5*p*(para.props[1,1]-para.props[3,1])*dens^(p-1)*u'*k*u
end

function selfWeight!(nodof, nels, nod, nf, g_num, disp, spacing, rho, densities, etype, np_types)
    # Include self-weight loads in load vector
    # Assumes SI units
    if np_types > 1
        if nodof == 2
            for iel in 1:nels
                for node in 1:nod
                    nf[2,g_num[node,iel]]!=0&&(disp[nf[2,g_num[node,iel]]] -= 0.25*prod(spacing)*(rho[etype[iel]]*densities[iel])*9.81)
                end
            end
        elseif nodof == 3

        end
    else
        if nodof == 2
            for iel in 1:nels
                for node in 1:nod
                 nf[2,g_num[node,iel]]!=0&&(disp[nf[2,g_num[node,iel]]] -= 0.25*prod(spacing)*(rho[1]*densities[iel])*9.81)
                end
            end
        elseif nodof == 3

        end
    end
    return disp
end

function eleDisp(ndof, disp, g_g)
    # Gather element displacements
    u = zeros(Float64, ndof)
    for dof in 1:ndof
        g_g[dof] != 0 && (u[dof] = disp[g_g[dof]])
    end
    return u
end

function compSensWeight(p, dens, para, iel, k, spacing)
    # Analytical calculation of compliance sensitivity including self-weight

    u = num.eleDisp(para.ndof, para.disp, para.g_g[:,iel])
    # Build vector that distributes element weight accross its nodes: [0 1 0 1 0 1 0 1]'
    # (loads in y DOFs of nodes)
    if para.nodof == 2
        q = Float64.([i%2!=0 ? 0 : 1 for i in 1:2*para.nod])
    else

    end
    size(para.props, 2) > 1 ? rho=para.props[3,para.etype[iel]] : rho=para.props[3,1]
    return (-0.5*p*dens^(p-1)*u'*k*u)[1]-rho*prod(spacing)*9.81*dot(q,u)/4
end

function volFrac(dens, nels)
    #=
    Fraction of current volume relative to initial value (fully solid structure).
    Assumes elements have the same dimensions.
    element_volume*sum(dens)/(nels*element_volume)
    element_volume = prod(spacing)*1
    This assumption also results in:
    sensitivity = d(volFrac)/dx_i = x_i/nels
    =#
    return sum(dens)/nels
end

function eigenValSens(para, vecs, p, dens, k, eigenVal, m, iel, mGlobal)
    sens = zeros(para.nmodes)
    for mode in 1:para.nmodes
        u = num.eleDisp(para.ndof, vecs[:,mode], para.g_g[:,iel])
        sens[mode] = (u'*(p*dens^(p-1)*k-eigenVal[mode]*m)*u)/(vecs[:,mode]'*mGlobal*vecs[:,mode])
    end
    return sens
end

function check(nelx,nely,rmin,x,dc)
    # Apply sensitivity filter
    dcn=zeros(length(dc))
    for iel in 1:length(dc)
        # Line of current element
        i=floor((iel-1)/nelx)+1
        # Column of current element
        j=iel-floor((iel-1)/nelx)*nelx
        # Scan (square) neighborhood of current element and apply filter
        global summ=0
        # Move vertically across neighborhood of current element (lines)
        for k = maximum([i-floor(rmin),1]):minimum([i+floor(rmin),nely])
            # Move horizontally across neighborhood of current element (columns)
            for l = maximum([j-floor(rmin),1]):minimum([j+floor(rmin),nelx])
                # Weighting factor H_f
                fac = rmin-sqrt((i-k)^2+(j-l)^2)
                # Accumulate weight factor H_f for elements inside filter radius
                global summ += maximum([0,fac])
                dcn[iel] += maximum([0,fac])*x[convert(Int32,(k - 1)*nelx + l)]*dc[convert(Int32,(k - 1)*nelx + l)]
                if k == minimum([i+floor(rmin),nely]) && l == minimum([j+floor(rmin),nelx])
                    dcn[iel] /= (x[iel]*summ)
                end
            end
        end
    end
    return dcn
end

function OC(nelx,nely,x,volfrac,dc)
    # Optimality criteria update
    global l1 = 0
    global l2 = 100000
    global move = 0.2
    while (l2-l1 > 1e-4)
        lmid = (l1+l2)/2
        global xnew = max.(0.001,max.(x.-move,min.(1,min.(x.+move,x.*sqrt.(abs.(dc)/lmid)))))
        if (sum(sum(xnew)) - volfrac*nelx*nely) > 0
            l1 = lmid
        else
            l2 = lmid
        end
    end
    return xnew
end

function buildFilterMatrix(para,radius,spacing)
            
    # % Estimating the number of elements in filter matrix: elements within filter
    # % area multiplied by the number of elements.
    nfilterel = convert(Int32, para.nels*4*round(radius/spacing[1])^2)

    if nfilterel != 0
        # % Selecting coordinates of design nodes
        X_den = zeros(para.nels)
        Y_den = zeros(para.nels)

        k=0
        # % Element weights
        R = zeros(nfilterel)
        # % Index vectors
        ifilt = zeros(nfilterel)
        jfilt = zeros(nfilterel)
        

        for iel in 1:para.nels
            aa = num.centerPos(para.element, para.coord, para.g_num, iel)
            X_den[iel] = aa[1]
            Y_den[iel] = aa[2]
        end

        # % Looping at the elements inside the design domain
        for el = 1:para.nels

            # % Weight function (linear with radial distance)
            r = max.(0, 1.0 .- (((X_den .- X_den[el]).^2 + (Y_den .- Y_den[el]).^2).^(1/2) / radius))
            #= r = the maximum between 0 and [1 - (distances between reference element and others)/filter_radius]
            r[i] = 0 if distance equal to or greater than filter radius
            r[i] cases for element_side < filter radius < 2*element_side rectangular mesh:
                4 if element in corner
                6 if element in border but not in corner
                9 otherwise
            =#
            notzero = length(filter(!iszero,r))

            k += notzero

            # % Index vectors building
            global ifilt[k-notzero+1:k] = el*ones(Int32, notzero)
            global jfilt[k-notzero+1:k] = findnz(sparse(r))[1]

            # % Weights update
            R[k-notzero+1:k] = filter(x->x>0,r)/sum(r)

            # if isinteger(log(10,el))
            #     println("\nelement $el")
            #     println("nfilterel $nfilterel")
            #     num.megaPrint(r,"r",true)
            #     println("k $k")
            #     num.megaPrint(ifilt,"ifilt",true)
            #     num.megaPrint(jfilt,"jfilt",true)
            # end

        end

        # % H filter matrix

        return sparse(ifilt[1:k],jfilt[1:k],R[1:k])
    end
end

function stiffSerialAssembly(par, output)

    # Element stiffness calculations and assembly into global stiffness matrix.
    # No parallelization
        
    t_kdiag = 0
    t_shapeDer = 0
    t_beemat = 0
    t_kmm = 0
    t_fsparvv = 0
    t_sparinn = 0
    t_spabacc = 0

    # Returns maximum bandwidth kdiag for each
    # row of a skyline storage system from g
    t_k = @timed kdiag = mtx.fkdiag(par.g_g, par.neq, par.nels)
    t_kdiag += t_k.time

    # Global stiffness matrix in skyline storage
    kv = zeros(last(kdiag))

    # Return elastic stress–strain dee matrix in 2D (plane strain) or
    # 3D. Uses Young’s modulus and Poisson’s ratio.
    if size(par.props, 2) == 1
        
        dee = num.deemat(par.props[1, 1], par.props[2, 1], par.nodof, par.type_2d)
        
    end
        
    # Output number of equations and skyline storage
    write(output, "There are $(par.neq) equations and the skyline storage is $(kdiag[par.neq])\n")
     
    println("Start of stiffness and assembly procedures")

    for iel = 1:par.nels

        # loop for all elements

        if size(par.props, 2) > 1
                
            # Return elastic stress–strain dee matrix in 2D (plane strain) or
            # 3D. Uses Young’s modulus and Poisson’s ratio.
            dee2 = num.deemat(par.props[1, par.etype[iel]], par.props[2, par.etype[iel]], par.nodof, par.type_2d)

        end
            
        # element stiffness
        km = zeros(Float64, par.ndof, par.ndof)
            
        # loop for all integration points of the current element
                        
        for i = 1:par.nip
                
            # Returns the shape function derivatives der at
            # the ith integrating point (isoparametric coordinates).
            t_d = @timed der = num.shapeDer(par.points, i, par.nod, par.element)
            t_shapeDer += t_d.time

            # jacobian's determinant and bee matrix (shape functions derivatives in cartesian system)
            # Refer to section 3.7.3 of the FEA book
            t_b = @timed bee, determ = num.beemat(der, par.coord[par.g_num[:, iel], :], par.type_2d, par.nodof, par.nod)
            t_beemat += t_b.time

            if par.type_2d == "axisymmetric"
                    
                # Returns the shape functions fun at the ith integrating point.
                fun = num.shapeFun(par.points, i, par.nod, par.element)
                global par.gc = fun*par.coord
                bee[4, 1:2:par.nodof*par.nod - 1] = fun[:]/par.gc[1]
                
            end
                
            if size(par.props, 2) > 1
                t_km = @timed km += bee'*dee2*bee*determ*par.weights[i]*par.gc[1]
                t_kmm += t_km.time                    
            else
                t_km = @timed km += bee'*dee*bee*determ*par.weights[i]*par.gc[1]
                t_kmm += t_km.time
            end
                
        end

        # Returns lower triangular global stiffness matrix kv stored
        # as a vector in skyline form.
        t_fsparv = @timed num.fsparv!(kv, km, par.g_g[:, iel], kdiag)
        t_fsparvv += t_fsparv.time

        iel == trunc(par.nels/2) && println("Halfway through stiffness and assembly procedures")

    end

    println("Done with stiffness and assembly procedures")
        

        
    ##########################################################################
    #----------------------- EQUATION SOLUTION-------------------------------#
    ##########################################################################
        
        
        
    # Choleski factorization and skyline storage of kv
    t_sparin = @timed kv = num.sparin!(kv, kdiag)
    t_sparinn += t_sparin.time

    # Find displacements solution by forward and back substitution
    t_spabac = @timed disp = num.spabac!(kv, kdiag, par.loads)
    t_spabacc += t_spabac.time


    println("kdiag   $(t_kdiag)")
    println("shapeDer   $(t_shapeDer)")
    println("beemat   $(t_beemat)")
    println("km   $(t_kmm)")
    println("fsparv   $(t_fsparvv)")
    println("sparin!   $(t_sparinn)")
    println("spabac!   $(t_spabacc)")
    println("Solver runtime:")

    return disp

end

function stiffSerialPCG(par, output)
        
    # Element stiffness calculations using preconditioned gradient method/pcg
    # (no assembly of global stiffness matrix). Section 3.5 of the FEA book.
    # No parallelization

        
    write(output, "There are $(par.neq) equations\n")
    diag_precon = zeros(par.neq)
    storkm = Array{Float64}(undef, (par.ndof, par.ndof, par.nels))
        
    # Return elastic stress–strain dee matrix in 2D (plane strain) or
    # 3D. Uses Young’s modulus and Poisson’s ratio.
    if size(par.props, 2) == 1
        
        dee = num.deemat(par.young[1], par.poisson[1], par.nodof, par.type_2d)
        
    end
                
    for iel = 1:par.nels

        # loop for all elements
        if size(par.props, 2) > 1
                
            # Return elastic stress–strain dee matrix in 2D (plane strain) or
            # 3D. Uses Young’s modulus and Poisson’s ratio.
            dee2 = num.deemat(par.young[par.etype[iel]], par.poisson[par.etype[iel]], par.nodof, par.type_2d)

        end
        
        # element stiffness
        km = zeros(par.ndof, par.ndof)
            
        # loop for all integration points of the current element
                        
        for i = 1:par.nip
                
            # Returns the shape function derivatives der at
            # the ith integrating point (isoparametric coordinates).
            der = num.shapeDer(par.points, i, par.nod, par.element)
    println("element $iel nip $i")
            # jacobian's determinant and bee matrix (shape functions derivatives in cartesian system)
            # Refer to section 3.7.3 of the FEA book
            bee, determ = num.beemat(der, par.coord[par.g_num[:, iel], :], par.type_2d, par.nodof, par.nod)
    
            if par.type_2d == "axisymmetric"
                    
                # Returns the shape functions fun at the ith integrating point.
                fun = num.shapeFun(par.points, i, par.nod, par.element)
                global par.gc = fun*par.coord
                bee[4, 1:2:par.nodof*par.nod - 1] = fun[:]/par.gc[1]
                
            end
                
            if size(par.props, 2) > 1
                    
                km += bee'*dee2*bee*determ*par.weights[i]*par.gc[1]
                    
            else
                                        
                km += bee'*dee*bee*determ*par.weights[i]*par.gc[1]
                    
            end
                
        end

        storkm[:, :, iel] = km
        
        #=
        Build preconditioner based on elements that would
        be on the main diagonal of the global stifness matrix.
        Each one of these values will  be inverted inside the
        "num.pcg()" function below. See pages 68-70 of the FEA book
        =#
        for k in 1:par.ndof
                
            if par.g_g[k, iel] != 0
                    
                diag_precon[par.g_g[k, iel]] += km[k,k] 
                    
            end
                
        end
            
    end
       

        
  ##########################################################################
  #-----------------------PCG EQUATION SOLUTION----------------------------#
  ##########################################################################
        
  par.disp = pcg!(par, diag_precon, storkm, output)

end

function stiffParaAssembly(par, output)

    # Element stiffness calculations and assembly into global stiffness matrix.
    # Solution parallelized on the GPU.


    # Returns maximum bandwidth kdiag for each
    # row of a skyline storage system from g
    kdiag = mtx.fkdiag(par.g_g, par.neq, par.nels)

    # Global stiffness matrix
    kv = Array{Float64}(zeros(last(kdiag)))

    # Initialization of bee matrix
    if par.type_2d == "plane"

        bee_d = CUDA.zeros(3, par.nodof*par.nod)

    elseif par.type_2d == "axisymmetric"

        bee_d = CUDA.zeros(4, par.nodof*par.nod)

    elseif par.nodof == 3

        bee_d = CUDA.zeros(6, par.nodof*par.nod)

    else

        println("Invalid type_2d/nodof combination. (Beginning of stiffness and assembly procedures)")
        
    end

    if size(par.props, 2) == 1
        
        # Return elastic stress–strain dee matrix in 2D (plane strain) or
        # 3D. Uses Young’s modulus and Poisson’s ratio.
        dee_d = cu(num.deemat(par.young[1], par.poisson[1], par.nodof, par.type_2d))

    end
        
    # Output number of equations and skyline storage
    write(output, "There are $(par.neq) equations and the skyline storage is $(kdiag[par.neq])\n")
    
    println("Start of stiffness and assembly procedures")

    # element stiffness
    #km = zeros(Float64, par.ndof, par.ndof, par.nels)

    for iel = 1:par.nels

        # loop for all elements

        # element stiffness
        km_d = CUDA.zeros(Float64, par.ndof, par.ndof)

        if size(par.props, 2) > 1
                
            # Return elastic stress–strain dee matrix in 2D (plane strain) or
            # 3D. Uses Young’s modulus and Poisson’s ratio.
            dee2_d = num.para_deemat(par.young[par.etype[iel]], par.poisson[par.etype[iel]], par.nodof, par.type_2d)

        end
            
        # loop for all integration points of the current element
                        
        for i = 1:par.nip
                
            # Returns the shape function derivatives der at
            # the ith integrating point (isoparametric coordinates).
            der = num.shapeDer(par.points, i, par.nod, par.element)
    
            # bee matrix (shape functions derivatives in cartesian
            #system), jacobian and its determinant.
            # Section 3.7.3 of the FEA book
            jac = cu(der)*cu(par.coord[par.g_num[:, iel], :])
            determ = det(Array(jac))
            deriv_d = CuArray{Float64}(undef, (par.nodof, par.nod))
            deriv_d = cu(inv(Array(jac))*der)
            @cuda threads = par.nod num.para_beemat!(bee_d, deriv_d)
            if par.type_2d == "axisymmetric"
                    
                # Returns the shape functions fun at the ith integrating point.
                fun = num.shapeFun(par.points, i, par.nod, par.element)
                global par.gc = fun*par.coord
                bee_d[4, 1:2:par.nodof*par.nod - 1] = fun[:]/par.gc[1]
                
            end

            if size(par.props, 2) > 1

                km_d[:, :] += bee_d'*dee2_d*bee_d*cu(determ)*cu(par.weights[i])*cu(par.gc[1])
                # km[:, :, iel] += bee'*dee*bee*determ*par.weights[i]*par.gc[1]
 
            else

                km_d[:, :] += bee_d'dee_d*bee_d*cu(determ)*cu(par.weights[i])*cu(par.gc[1])
                # km[:, :, iel] += bee'*dee*bee*determ*par.weights[i]*par.gc[1]

            end
                
        end

        iel == trunc(par.nels/2) && println("Halfway through stiffness and assembly procedures")
        
        # Returns lower triangular global stiffness matrix kv stored
        # as a vector in skyline form.
        km = Array(km_d)
        num.fsparv!(kv, km, par.g_g[:, iel], kdiag)
        
    end

    #
        # for linha in 1:size(par.g_g, 1)

        #     linha == 1 && (global nulo = 0)
        #     for coluna in 1:size(par.g_g, 2)

        #         par.g_g[linha, coluna] == 0 && (nulo += 1)
                
        #     end
            
        # end

        # Inside the kernel num.para_fsparv!, this matrix will
        # store the values to be input later into kv
        #preKv = zeros(Int16, par.ndof^2 - par.ndof, par.nels)

        # Inside the kernel num.para_fsparv!, this matrix will
        # store the respective positions to be altered in kv
        #preKv_indices = zeros(Int16, (par.ndof^2 - par.ndof, par.nels))

        #intermPos = zeros(Int16, size(par.g_g, 2))

        # CUDA.@time @cuda threads = par.nels para_fsparv!(cu(preKv_indices), cu(preKv), cu(km), cu(par.g_g), cu(kdiag), cu(intermPos))
        # preKv_indices = Array(preKv_indices_d)
        # preKv = Array(preKv_d)

        #=
        preKv and preKv_indices store the information necessary
        to "build" kv, which wasn't done inside of num.para_fsparv!
        because of synchronization and control flow (i.e. "if"
        statements) issues. Moreover, since this task will be done
        in a serialized manner, it's better to be done in the CPU and not
        the GPU (according to item 5.2.1 "Application level" of the
        official CUDA C++ programming guide).
        =# 
        # for ele in 1:par.nels

        #     for DOF in 1:par.ndof

        #         if preKv_indices[DOF, ele] != 0

        #             kv[preKv_indices[DOF, ele]] += preKv[DOF, ele]

        #         end

        #     end

        # end
        # println("\n\n$(kv[1:10])\n\n")
    #


    println("Done with stiffness and assembly procedures")


        
    ##########################################################################
    #----------------------- EQUATION SOLUTION-------------------------------#
    ##########################################################################
        
        

    # Choleski factorization and skyline storage of kv
    kv = num.sparin!(kv, kdiag)
    # println("sparin!\n$(kv[1:100])")

    # Find displacements solution
    par.disp = num.spabac!(kv, kdiag, par.disp)
    # println("spabac!\n$(par.disp[1:100])")
    
    return par.disp
    
end

function eigenSimple(par)
    # Solver for eigenvalues/eigenvectors problem
    
    # Return elastic stress–strain dee matrix in 2D (plane strain) or
    # 3D. Uses Young’s modulus and Poisson’s ratio.
    if size(par.props, 2) == 1
        
        dee = num.deemat(par.props[1, 1], par.props[2, 1], par.nodof, par.type_2d)
        
    end
    
    # Output number of equations and skyline storage
    #InOut.neqStorageOut(par.neq, kdiag[end], output)
    
    println("Start of assembly procedures")
    
    # Global stiffness matrix
    kmGlobal = zeros(Float64, (par.neq, par.neq))
    # Global mass matrix
    mmGlobal = zeros(Float64, (par.neq, par.neq))
    # Loop elements for stiffness and mass information
    for iel = 1:par.nels
        
        # element stiffness
        
        if size(par.props, 2) > 1
            
            # Return elastic stress–strain dee matrix in 2D (plane strain) or
            # 3D. Uses Young’s modulus and Poisson’s ratio.
            dee2 = num.deemat(par.props[1, par.etype[iel]], par.props[2, par.etype[iel]], par.nodof, par.type_2d)
            
        end
        
        # element stiffness
        km = zeros(Float64, par.ndof, par.ndof)
        
        # element mass
        mm = zeros(Float64, par.ndof, par.ndof)
        
        # loop for all integration points of the current element
        for i = 1:par.nip
            
            # Returns the shape function derivatives der at
            # the ith integrating point (isoparametric coordinates).
            der = num.shapeDer(par.points, i, par.nod, par.element)
            fun = num.shapeFun(par.points, i, par.nod, par.element)
            
            # jacobian's determinant and bee matrix (shape functions derivatives in cartesian system)
            # Refer to section 3.7.3 of the FEA book
            bee, determ = num.beemat(der, par.coord[par.g_num[:, iel], :], par.type_2d, par.nodof, par.nod)
            
            if par.type_2d == "axisymmetric"
                
                # Returns the shape functions fun at the ith integrating point.
                fun = num.shapeFun(par.points, i, par.nod, par.element)
                global par.gc = fun*par.coord
                bee[4, 1:2:par.nodof*par.nod - 1] = fun[:]/par.gc[1]
                
            end
            
            if size(par.props, 2) > 1
                km += bee'*dee2*bee*determ*par.weights[i]*par.gc[1]
            else
                km += bee'*dee*bee*determ*par.weights[i]*par.gc[1]
            end

            if size(par.props, 2) > 1
                mm = num.mass!(mm, num.ecmat(par.nod, par.nodof, par.ndof, fun), par.props[3,par.etype[iel]], determ, par.weights[i], i, par.nip)
            else
                mm = num.mass!(mm, num.ecmat(par.nod, par.nodof, par.ndof, fun), par.props[3,1], determ, par.weights[i], i, par.nip)
            end
            
        end
        
        # Global matrices assembly
        kmGlobal, mmGlobal = num.globals!(kmGlobal, mmGlobal, mm, km, par.g_g[:, iel])
        iel == trunc(par.nels/2) && println("Halfway through stiffness and assembly procedures")
        
    end
    println("Done with stiffness and assembly procedures")

    # Extract eigenvalues and eigenvectors from generalized eigenproblem
    eigenval, eigenvecs = eigen(kmGlobal, mmGlobal)
    freqs = sqrt.(eigenval[1:par.nmodes])./=(2*pi)
    println(freqs)

    return eigenvecs[:,par.eigenvec]

end

function eigenSerialAssembly(par, output)
    # Solver for eigenvalues/eigenvectors problem

    # Returns maximum bandwidth kdiag for each
    # row of a skyline storage system from g
    kdiag = mtx.fkdiag(par.g_g, par.neq, par.nels)

    # Determine maximum bandwidth of elements
    nband = 0
    for iel in 1:par.nels
        (nband < mtx.bandwidth(par.g_g[:, iel])) && (nband = mtx.bandwidth(par.g_g[:, iel]))
    end
    InOut.bandwidthOut(nband, output)
        
    # Output number of equations and skyline storage
    InOut.neqStorageOut(par.neq, kdiag[end], output)

    # Return elastic stress–strain dee matrix in 2D (plane strain) or
    # 3D. Uses Young’s modulus and Poisson’s ratio.
    if size(par.props, 2) == 1
        
        dee = num.deemat(par.props[1, 1], par.props[2, 1], par.nodof, par.type_2d)
        
    end
     
    println("Start of stiffness and assembly procedures")

    # Lumped global mass matrix as a vector
    diag = zeros(Float64, par.neq)
    # Global stiffness matrix stored as upper triangle
    ku = zeros(Float64, (par.neq,nband+1))
    # Loop elements for stiffness and mass information
    for iel = 1:par.nels
        
        
        # Element area
        yUpLeft = par.coord[par.g_num[3, iel], 2]
        xDownRight = par.coord[par.g_num[2, iel], 1]
        xDownLeft = par.coord[par.g_num[1, iel], 1]
        yDownLeft = par.coord[par.g_num[1, iel], 2]
        area = (xDownRight - xDownLeft)*(yUpLeft - yDownLeft)
        
        # element stiffness
        
        if size(par.props, 2) > 1
            
            # Return elastic stress–strain dee matrix in 2D (plane strain) or
            # 3D. Uses Young’s modulus and Poisson’s ratio.
            dee2 = num.deemat(par.props[1, par.etype[iel]], par.props[2, par.etype[iel]], par.nodof, par.type_2d)
            
        end
        
        # element stiffness
        km = zeros(Float64, par.ndof, par.ndof)
        
        # loop for all integration points of the current element
        
        for i = 1:par.nip
            
            # Returns the shape function derivatives der at
            # the ith integrating point (isoparametric coordinates).
            der = num.shapeDer(par.points, i, par.nod, par.element)
            
            # jacobian's determinant and bee matrix (shape functions derivatives in cartesian system)
            # Refer to section 3.7.3 of the FEA book
            bee, determ = num.beemat(der, par.coord[par.g_num[:, iel], :], par.type_2d, par.nodof, par.nod)
            
            if par.type_2d == "axisymmetric"
                
                # Returns the shape functions fun at the ith integrating point.
                fun = num.shapeFun(par.points, i, par.nod, par.element)
                global par.gc = fun*par.coord
                bee[4, 1:2:par.nodof*par.nod - 1] = fun[:]/par.gc[1]
                
            end
            
            if size(par.props, 2) > 1
                km += bee'*dee2*bee*determ*par.weights[i]*par.gc[1]
            else
                km += bee'*dee*bee*determ*par.weights[i]*par.gc[1]
            end
            
        end
        
        # Element mass matrix
        if size(par.props, 2) > 1
            mm = num.elmat(area,par.props[3,par.etype[iel]], par.element)
        else
            mm = num.elmat(area,par.props[3,1], par.element)
        end
        # Global mass assembly
        mtx.formlump!(diag,mm,par.g_g[:, iel])
        # Global stiffness assembly
        ku = num.formku!(ku,km,par.g_g[:, iel])
        iel == trunc(par.nels/2) && println("Halfway through stiffness and assembly procedures")
        
    end
    println("Done with stiffness and assembly procedures")

    # Reduce to standard eigenvalue problem
    rrmass, ku = num.StandardEigen!(diag, nband, ku)
    # Convert to skyline form
    kh = mtx.EigenSkyline!(kdiag, ku)
    # Tridiagonalize stiffness matrix
    diag, udiag = mtx.bandred!(ku, par.neq)
    # Extract eigenvalues
    diag = num.bisect!(diag,udiag,par.etol,par.nmodes, output)
    vals = sqrt.(diag[1:par.nmodes])
    vals ./= (2*pi)
    println(vals)
    # Extract eigenvectors
    num.EigenVec(par.nmodes, kh, kdiag, diag, par.penalty, udiag, rrmass, output)

end

function staticSens(par, densities)
    
    # Solver for static cases
    
    # Return elastic stress–strain dee matrix in 2D (plane strain) or
    # 3D. Uses Young’s modulus and Poisson’s ratio.
    if size(par.props, 2) == 1
        
        dee = num.deemat(par.props[1, 1], par.props[2, 1], par.nodof, par.type_2d)
        
    end
    
    # Output number of equations and skyline storage
    # InOut.neqStorageOut(par.neq, kdiag[end], output)
    
    # Global stiffness
    kmGlobal = zeros(Float64, (par.neq, par.neq))

    # Loop elements for stiffness and mass information
    for iel = 1:par.nels
        
        # element stiffness
        
        if size(par.props, 2) > 1
            
            # Return elastic stress–strain dee matrix in 2D (plane strain) or
            # 3D. Uses Young’s modulus and Poisson’s ratio.
            dee2 = num.deemat(par.props[1, par.etype[iel]], par.props[2, par.etype[iel]], par.nodof, par.type_2d)
            
        end
        
        # element stiffness
        km = zeros(Float64, par.ndof, par.ndof)
        
        # loop for all integration points of the current element
        for i = 1:par.nip
            
            # Returns the shape function derivatives der at
            # the ith integrating point (isoparametric coordinates).
            der = num.shapeDer(par.points, i, par.nod, par.element)
            
            # jacobian's determinant and bee matrix (shape functions derivatives in cartesian system)
            # Refer to section 3.7.3 of the FEA book
            bee, determ = num.beemat(der, par.coord[par.g_num[:, iel], :], par.type_2d, par.nodof, par.nod)
            
            if par.type_2d == "axisymmetric"
                
                # Returns the shape functions fun at the ith integrating point.
                fun = num.shapeFun(par.points, i, par.nod, par.element)
                global par.gc = fun*par.coord
                bee[4, 1:2:par.nodof*par.nod - 1] = fun[:]/par.gc[1]
                
            end
            
            if size(par.props, 2) > 1
                km += bee'*dee2*bee*determ*par.weights[i]*par.gc[1]
            else
                km += bee'*dee*bee*determ*par.weights[i]*par.gc[1]
            end
            
        end
        
        # Store element stiffness matrices
        if iel == 1
            global elementsKm = collect([km])
        else
            push!(elementsKm, km)
        end
        # Global stiffness assembly
        kmGlobal = num.globalKm!(kmGlobal, km*densities[iel], par.g_g[:, iel])

    end
    
    # InOut.megaPrint(densities, "densities"; values=true)
    println(maximum(densities))
    println(minimum(densities))
    println()
    
    disp = ldlt(sparse(Symmetric(kmGlobal)))\par.loads
    
    return disp, elementsKm, 0.5*dot(par.loads,disp)

end

function staticSensExt(par, densities)
    
    # Solver for static cases
    
    # Return elastic stress–strain dee matrix in 2D (plane strain) or
    # 3D. Uses Young’s modulus and Poisson’s ratio.
    if size(par.props, 2) == 1

        dee = num.deemat(par.props[1, 1], par.props[2, 1], par.nodof, par.type_2d)

    end
    
    # Output number of equations and skyline storage
    # InOut.neqStorageOut(par.neq, kdiag[end], output)
    
    # Global stiffness
    kmGlobal = zeros(Float64, (par.neq, par.neq))

    # Loop elements for stiffness and mass information
    for iel = 1:par.nels
        
        # element stiffness
        
        if size(par.props, 2) > 1
            
            # Return elastic stress–strain dee matrix in 2D (plane strain) or
            # 3D. Uses Young’s modulus and Poisson’s ratio.
            dee2 = num.deemat(par.props[1, par.etype[iel]], par.props[2, par.etype[iel]], par.nodof, par.type_2d)
            
        end
        
        # element stiffness
        km = zeros(Float64, par.ndof, par.ndof)
        
        # loop for all integration points of the current element
        for i = 1:par.nip
            
            # Returns the shape function derivatives der at
            # the ith integrating point (isoparametric coordinates).
            der = num.shapeDer(par.points, i, par.nod, par.element)
            
            # jacobian's determinant and bee matrix (shape functions derivatives in cartesian system)
            # Refer to section 3.7.3 of the FEA book
            bee, determ = num.beemat(der, par.coord[par.g_num[:, iel], :], par.type_2d, par.nodof, par.nod)
            
            if par.type_2d == "axisymmetric"
                
                # Returns the shape functions fun at the ith integrating point.
                fun = num.shapeFun(par.points, i, par.nod, par.element)
                global par.gc = fun*par.coord
                bee[4, 1:2:par.nodof*par.nod - 1] = fun[:]/par.gc[1]
                
            end
            
            if size(par.props, 2) > 1
                km += bee'*dee2*bee*determ*par.weights[i]*par.gc[1]
            else
                km += bee'*dee*bee*determ*par.weights[i]*par.gc[1]
            end
            
        end
        
        # Store element stiffness matrices
        if iel == 1
            global elementsKm = collect([km])
        else
            push!(elementsKm, km)
        end
        
        # Global stiffness assembly
        kmGlobal = num.globalKm!(kmGlobal, km*(par.props[3,1]+densities[iel]*(par.props[1,1]-par.props[3,1])),
                        par.g_g[:, iel])

    end
    println(maximum(densities))
    # disp = ldlt(sparse(Symmetric(kmGlobal)))\par.loads
    # disp = sparse((kmGlobal+kmGlobal')/2)\par.loads
    disp = sparse(kmGlobal)\par.loads
    println(minimum(densities))
    println()
    
    return disp, elementsKm, 0.5*dot(par.loads,disp)

end

function eigenSens(par, densities, SIMPpenalty)
    # Solver for eigenvalues/eigenvectors problem
    
    # Return elastic stress–strain dee matrix in 2D (plane strain) or
    # 3D. Uses Young’s modulus and Poisson’s ratio.
    if size(par.props, 2) == 1
        
        dee = num.deemat(par.props[1, 1], par.props[2, 1], par.nodof, par.type_2d)
        
    end
    
    # Output number of equations and skyline storage
    #InOut.neqStorageOut(par.neq, kdiag[end], output)
    
    # Global stiffness matrix
    kmGlobal = zeros(Float64, (par.neq, par.neq))
    # Global mass matrix
    mmGlobal = zeros(Float64, (par.neq, par.neq))
    # Loop elements for stiffness and mass information
    for iel = 1:par.nels
        
        # element stiffness
        if size(par.props, 2) > 1
            
            # Return elastic stress–strain dee matrix in 2D (plane strain) or
            # 3D. Uses Young’s modulus and Poisson’s ratio.
            dee2 = num.deemat(par.props[1, par.etype[iel]], par.props[2, par.etype[iel]], par.nodof, par.type_2d)
            
        end

        # element stiffness
        km = zeros(Float64, par.ndof, par.ndof)
        
        # element mass
        mm = zeros(Float64, par.ndof, par.ndof)
        
        # loop for all integration points of the current element
        for i = 1:par.nip
            
            # Returns the shape function derivatives der at
            # the ith integrating point (isoparametric coordinates).
            der = num.shapeDer(par.points, i, par.nod, par.element)
            fun = num.shapeFun(par.points, i, par.nod, par.element)
            
            # jacobian's determinant and bee matrix (shape functions derivatives in cartesian system)
            # Refer to section 3.7.3 of the FEA book
            bee, determ = num.beemat(der, par.coord[par.g_num[:, iel], :], par.type_2d, par.nodof, par.nod)
            
            if par.type_2d == "axisymmetric"
                
                # Returns the shape functions fun at the ith integrating point.
                fun = num.shapeFun(par.points, i, par.nod, par.element)
                global par.gc = fun*par.coord
                bee[4, 1:2:par.nodof*par.nod - 1] = fun[:]/par.gc[1]
                
            end

            # Element stiffness
            if size(par.props, 2) > 1
                km += bee'*dee2*bee*determ*par.weights[i]*par.gc[1]
            else
                km += bee'*dee*bee*determ*par.weights[i]*par.gc[1]
            end

            # Element mass
            if size(par.props, 2) > 1
                mm = num.mass!(mm, num.ecmat(par.nod, par.nodof, par.ndof, fun), par.props[3,par.etype[iel]], determ, par.weights[i], i, par.nip)
            else
                mm = num.mass!(mm, num.ecmat(par.nod, par.nodof, par.ndof, fun), par.props[3,1], determ, par.weights[i], i, par.nip)
            end
            
        end

        # Global matrices assembly
        kmGlobal, mmGlobal = num.globals!(kmGlobal, mmGlobal, mm*densities[iel], km*(densities[iel]^SIMPpenalty), par.g_g[:, iel])

        # Store element matrices
        if iel == 1
            global elementsKm = collect([km])
            global elementsMm = collect([mm])
        else
            push!(elementsKm, km)
            push!(elementsMm, mm)
        end

    end

    # Extract eigenvalues and eigenvectors from generalized eigenproblem
    eigenval, eigenvecs = eigen(kmGlobal, mmGlobal)

    return eigenval[1:par.nmodes], eigenvecs[:,1:par.nmodes], elementsKm, elementsMm, mmGlobal

end


end