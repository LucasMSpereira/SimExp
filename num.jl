module num

using Printf, LinearAlgebra, CUDA
using Base: print_array
import mtx

### This module provides numerical routines

function sample(element, nip, process)

    # Returns the local coordinates and weighting coefficients
    # for numerical integration of an element (gauss quadrature)

    write(process, "\nnum.sample\n\n")

    if element == "quadrilateral"

        if nip == 1

            points = [0 0]
        
            weights = 4

        elseif nip == 4

            points = [sqrt(1/3) sqrt(1/3);
                sqrt(1/3) -sqrt(1/3);
                -sqrt(1/3) sqrt(1/3);
                -sqrt(1/3) -sqrt(1/3)]
            
            weights = 1

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

    write(process, "points:\n")
    for i in 1:size(points, 1)

        write(process, "$(points[i, :])\n")
  
    end
    
    write(process, "weights:\n")
    for i in 1:size(weights, 1)

        write(process, "$(weights[i, :])\n")
  
    end

    return points, weights
        
end

function deemat(e, v, dim, type, process)

    # Calculates stress-strain matrix dee

    write(process,  "\nnum.deemat\n\n")

    if dim == 2 && type == "plane"

        dee = e*(1 - v)/((1 + v)*(1 - 2 * v))*
        [1 v/(1 - v) 0;
        v/(1 - v) 1 0;
        0 0 (1 - 2*v)/(2*(1 - v))]

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

    else

        write(process, "Invalid combination of dimension and type of analysis. (function num.deemat)\n")

    end

    write(process, "dee:\n")
    for i in 1:size(dee, 1)

        write(process, "$(dee[i, :])\n")
        
    end

    return dee
        
end

function shapeFun(points, nip, nod, element, process)

    # Return shape functions at the current integrating point
    # Ref.: "Introduction to Finite Element Analysis Using MATLAB® and Abaqus", CRC, pg 207-210

    write(process, "\nnum.shapefun\n\n")

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

            h1 = 0.25*(1 + isoX)*(1 + isoY)
            h2 = 0.25*(1 - isoX)*(1 + isoY)
            h3 = 0.25*(1 - isoX)*(1 - isoY)
            h4 = 0.25*(1 + isoX)*(1 + isoY)

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

    write(process, "fun:\n")
    for i in 1:size(fun, 1)

        write(process, "$(fun[i, :])\n")
        
    end

    return fun

end

function shapeDer(points, nip, nod, element, process)

    # Returns ('n° of DOFs/node' x 'n° of nodes/element')  matrix with shape
    # functions derivatives in isoparametric coordinates

    write(process, "\nnum.shapeDer\n\n")

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

    write(process, "der:\n")
    for i in 1:size(der, 1)

        write(process, "$(der[i, :])\n")
        
    end

    return der
    
end

function beemat(der, coord, type_2d, nodof, nod, process)

    # Returns the derivatives of the shape functions in the cartesian
    # system (bee) and the determinant of the jacobian

    write(process, "\nnum.beemat\n\n")

    jac = der*coord
    determ = det(jac)
    deriv = inv(jac)*der

    if nodof == 2

        if type_2d == "plane"

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

        else

            write(process, "Invalid type of 2D analysis. (function num.beemat)\n")
            
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

    else

        write(process, "Invalid dimension (nodof) of analysis. (function num.beemat)\n")

    end

    write(process, "Jacobian:\n")
    for i in 1:size(jac, 1)

        write(process, "$(jac[i, :])\n")
        
    end

    write(process, "\nbee:\n")
    for i in 1:size(bee, 1)

        write(process, "$(bee[i, :])\n")
        
    end

    write(process, "\nJacobian determinant = $(determ)\n")

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

function fsparv!(kv, km, g, kdiag, process)

    # Assembles element matrices into a symmetric global matrix stored in skyline form.

    write(process, "\nnum.fsparv\n\n")

    for i in 1:size(g, 1)
        
        if g[i] != 0

            for j in 1:size(g, 1)
                
                if g[j] != 0 && (g[i] - g[j]) >= 0
                    
                    kv[kdiag[g[i]] - (g[i] - g[j])] += km[i, j]
    
                end

            end

        end

    end

    # write(process, "kv:\n")
    # for i in 1:size(kv, 1)

    #     write(process, "$(kv[i, :])\n")
        
    # end

    return kv
    
end

function sparin!(kv, kdiag, process)

    # Returns the (Choleski) factorised vector kv stored as a skyline

    write(process, "\nnum.sparin\n\n")

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

    # write(process, "kv:\n")
    # for i in 1:size(kv, 1)

    #     write(process, "$(kv[i, :])\n")
        
    # end

    return kv

end

function spabac!(kv, kdiag, disp, process)

    # Cholesky forward and back substitution on symmetric skyline global matrix

    write(process, "\nnum.spabac\n\n")

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
    write(process, "Displacements:\n")
    for i in 1:size(disp, 1)

        write(process, "$(disp[i, :])\n")
        
    end

    return disp
    
end

function pcg!(parameters, diag_precon, storkm, process, output)

    # Preconditioned conjugate gradient (pcg) method for solving FEA system
    # "element-by-element" (without assembly of global stiffness matrix).
    # Section 3.5 of the FEA book

    
    write(process, "\nnum.pcg\n")

    for i in 1:size(diag_precon, 1)

        diag_precon[i] != 0 ? diag_precon[i] = 1/diag_precon[i] : false
    
    end
    cg_iters = 0
    x = zeros(parameters.neq)
    d = zeros(parameters.neq)
    xnew = zeros(parameters.neq)
    x = zeros(parameters.neq)
    d .= diag_precon.*parameters.disp[:, 1]
    p = d

    write(process, "\nIteration     Convergence criterion\n")
    for i in 1:parameters.cg_limit

        cg_iters += 1
        u = zeros(parameters.neq)
        for iel in 1:parameters.nels
            
            km = storkm[:, :, iel]
            m = zeros((size(parameters.g_g, 1), 1))
            for b in 1:size(parameters.g_g, 1)

                parameters.g_g[b, iel] != 0 ? m[b] = p[parameters.g_g[b, iel]] : false
                
            end
            l = km*m
            for b in 1:size(parameters.g_g, 1)

                parameters.g_g[b, iel] != 0 ? u[parameters.g_g[b, iel]] += l[b] : false
                
            end

        end
        # fixed_freedoms != 0 ? u[no] = p[no]*store : false
        up = parameters.disp[:, i]'*d
        alpha = up/(p'*u)
        xnew = x + alpha*p
        parameters.disp[:, i + 1] = parameters.disp[:, i] - alpha*u
        d .= diag_precon.*parameters.disp[:, i + 1]
        beta = (parameters.disp[:, i + 1]'*d)/up
        p = d + beta*p
        @printf process "%i" i
        @printf process "\t\t\t%.4e\n" maximum(abs.(xnew - x))/maximum(abs.(xnew))
        println("\n$(i)     $(maximum(abs.(xnew - x))/maximum(abs.(xnew)))\n")
        maximum(abs.(xnew - x))/maximum(abs.(xnew)) <= parameters.cg_tol ? break : false
        x = xnew

    end

    if (maximum(abs.(xnew - x))/maximum(abs.(xnew))) <= parameters.cg_tol

        write(output, "\nNumber of PCG iterations to convergence was $(cg_iters).\n\n")
        println("\nNumber of PCG iterations to convergence was $(cg_iters).\n")

    else

        write(output, "\n *** WARNING: PCG didn't convergence in $(parameters.cg_limit) iterations. ***\n\n")
        println("\n *** WARNING: PCG didn't convergence in $(parameters.cg_limit) iterations. ***\n\n")
        
    end

    return x
    
end

function stiffSerialAssembly(parameters, process, output)

    # Element stiffness calculations and assembly into global stiffness matrix.
    # No parallelization
        
    write(process, "\nnum.stiffSerialAssembly\n")

    t_kdiag = 0
    t_shapeDer = 0
    t_beemat = 0
    t_kmm = 0
    t_fsparvv = 0
    t_sparinn = 0
    t_spabacc = 0

    # Returns maximum bandwidth kdiag for each
    # row of a skyline storage system from g
    t_k = @timed kdiag = mtx.fkdiag(parameters.g_g, parameters.neq, parameters.nels, process)
    t_kdiag += t_k.time

    # Global stiffness matrix
    kv = zeros(last(kdiag))

    # Return elastic stress–strain dee matrix in 2D (plane strain) or
    # 3D. Uses Young’s modulus and Poisson’s ratio.
    if parameters.np_types == 1
        
        dee = num.deemat(parameters.young[1], parameters.poisson[1], parameters.nodof, parameters.type_2d, process)
        
    end
        
    # Output number of equations and skyline storage
    write(output, "There are $(parameters.neq) equations and the skyline storage is $(kdiag[parameters.neq])\n")
     
    println("Start of stiffness and assembly procedures")

    for iel = 1:parameters.nels

        # loop for all elements

        if parameters.np_types > 1
                
            # Return elastic stress–strain dee matrix in 2D (plane strain) or
            # 3D. Uses Young’s modulus and Poisson’s ratio.
            dee2 = num.deemat(parameters.young[parameters.etype[iel]], parameters.poisson[parameters.etype[iel]], parameters.nodof, parameters.type_2d, process)

        end
            
        # element stiffness
        km = zeros(Float32, parameters.ndof, parameters.ndof)
            
        # loop for all integration points of the current element
            
        write(process, "\n/////////////////  SITFFNESS OF ELEMENT $(iel) //////////////\n")
            
        for i = 1:parameters.nip
                
            # Returns the shape function derivatives der at
            # the ith integrating point (isoparametric coordinates).
            t_d = @timed der = num.shapeDer(parameters.points, i, parameters.nod, parameters.element, process)
            t_shapeDer += t_d.time

            # jacobian and bee matrix (shape functions derivatives in cartesian system)
            # Refer to section 3.7.3 of the FEA book
            t_b = @timed bee, determ = num.beemat(der, parameters.coord[parameters.g_num[:, iel], :], parameters.type_2d, parameters.nodof, parameters.nod, process)
            t_beemat += t_b.time

            if parameters.type_2d == "axisymmetric"
                    
                # Returns the shape functions fun at the ith integrating point.
                fun = num.shapeFun(parameters.points, i, parameters.nod, parameters.element, process)
                global parameters.gc = fun*parameters.coord
                bee[4, 1:2:parameters.nodof*parameters.nod - 1] = fun[:]/parameters.gc[1]
                
            end
                
            if parameters.np_types > 1
                    
                t_km = @timed km += bee'*dee2*bee*determ*parameters.weights[i]*parameters.gc[1]
                t_kmm += t_km.time                    
   
            else
                                        
                t_km = @timed km += bee'*dee*bee*determ*parameters.weights[i]*parameters.gc[1]
                t_kmm += t_km.time                    
            end
                
        end

        # Returns lower triangular global stiffness matrix kv stored
        # as a vector in skyline form.
        t_fsparv = @timed num.fsparv!(kv, km, parameters.g_g[:, iel], kdiag, process)
        t_fsparvv += t_fsparv.time

        iel == trunc(parameters.nels/2) && println("Halfway through stiffnes and assembly procedures")

    end

    println("Done with stiffness and assembly procedures")
        

        
    ##########################################################################
    #----------------------- EQUATION SOLUTION-------------------------------#
    ##########################################################################
        
        
        
    # Choleski factorization and skyline storage of kv
    t_sparin = @timed kv = num.sparin!(kv, kdiag, process)
    t_sparinn += t_sparin.time

    # Find displacements solution
    t_spabac = @timed parameters.disp = num.spabac!(kv, kdiag, parameters.disp, process)
    t_spabacc += t_spabac.time

    println("kdiag   $(t_kdiag)")
    println("shapeDer   $(t_shapeDer)")
    println("beemat   $(t_beemat)")
    println("km   $(t_kmm)")
    println("fsparv   $(t_fsparvv)")
    println("sparin!   $(t_sparinn)")
    println("spabac!   $(t_spabacc)")
    println("Solver runtime:")

    return parameters.disp

end

function stiffSerialPCG(parameters, process, output)
        
    # Element stiffness calculations using preconditioned gradient method/pcg
    # (no assembly of global stiffness matrix). Section 3.5 of the FEA book.
    # No parallelization

    write(process, "\nnum.stiffSerialPCG\n")
        
    write(output, "There are $(parameters.neq) equations\n")
    diag_precon = zeros(parameters.neq)
    storkm = Array{Float32}(undef, (parameters.ndof, parameters.ndof, parameters.nels))
        
    # Return elastic stress–strain dee matrix in 2D (plane strain) or
    # 3D. Uses Young’s modulus and Poisson’s ratio.
    if parameters.np_types == 1
        
        dee = num.deemat(parameters.young[1], parameters.poisson[1], parameters.nodof, parameters.type_2d, process)
        
    end
                
    for iel = 1:parameters.nels

        # loop for all elements

        if parameters.np_types > 1
                
            # Return elastic stress–strain dee matrix in 2D (plane strain) or
            # 3D. Uses Young’s modulus and Poisson’s ratio.
            dee2 = num.deemat(parameters.young[parameters.etype[iel]], parameters.poisson[parameters.etype[iel]], parameters.nodof, parameters.type_2d, process)

        end
        
        # element stiffness
        km = zeros(parameters.ndof, parameters.ndof)
            
        # loop for all integration points of the current element
            
        write(process, "\n/////////////////  SITFFNESS OF ELEMENT $(iel) //////////////\n")
            
        for i = 1:parameters.nip
                
            # Returns the shape function derivatives der at
            # the ith integrating point (isoparametric coordinates).
            der = num.shapeDer(parameters.points, i, parameters.nod, parameters.element, process)
    println("element $iel nip $i")
            # jacobian and bee matrix (shape functions derivatives in cartesian system)
            # Refer to section 3.7.3 of the FEA book
            bee, determ = num.beemat(der, parameters.coord[parameters.g_num[:, iel], :], parameters.type_2d, parameters.nodof, parameters.nod, process)
    
            if parameters.type_2d == "axisymmetric"
                    
                # Returns the shape functions fun at the ith integrating point.
                fun = num.shapeFun(parameters.points, i, parameters.nod, parameters.element, process)
                global parameters.gc = fun*parameters.coord
                bee[4, 1:2:parameters.nodof*parameters.nod - 1] = fun[:]/parameters.gc[1]
                
            end
                
            if parameters.np_types > 1
                    
                km += bee'*dee2*bee*determ*parameters.weights[i]*parameters.gc[1]
                    
            else
                                        
                km += bee'*dee*bee*determ*parameters.weights[i]*parameters.gc[1]
                    
            end
                
        end

        iel == trunc(parameters.nels/2) && println("Halfway through stiffness and assembly procedures")
        storkm[:, :, iel] = km
        
        #=
        Build preconditioner based on elements that would
        be on the main diagonal of the global stifness matrix.
        Each one of these values will then be inverted inside the
        "num.pcg()" function below. See pages 68-70 of the FEA book
        =#
        for k in 1:parameters.ndof
                
            if parameters.g_g[k, iel] != 0
                    
                diag_precon[parameters.g_g[k, iel]] += km[k,k] 
                    
            end
                
        end
            
    end
       

        
  ##########################################################################
  #-----------------------PCG EQUATION SOLUTION----------------------------#
  ##########################################################################
        
  parameters.disp = pcg!(parameters, diag_precon, storkm, process, output)

end

function stiffParaAssembly(parameters, process, output)

    # Element stiffness calculations and assembly into global stiffness matrix.
    # Solution parallelized on the GPU.

    write(process, "\nnum.stiffParaAssembly\n")

    # Returns maximum bandwidth kdiag for each
    # row of a skyline storage system from g
    kdiag = mtx.fkdiag(parameters.g_g, parameters.neq, parameters.nels, process)

    # Global stiffness matrix
    kv = Array{Float32}(zeros(last(kdiag)))

    # Initialization of bee matrix
    if parameters.type_2d == "plane"

        bee_d = CUDA.zeros(3, parameters.nodof*parameters.nod)

    elseif parameters.type_2d == "axisymmetric"

        bee_d = CUDA.zeros(4, parameters.nodof*parameters.nod)

    elseif parameters.nodof == 3

        bee_d = CUDA.zeros(6, parameters.nodof*parameters.nod)

    else

        println("Invalid type_2d/nodof combination. (Beginning of stiffness and assembly procedures)")
        
    end

    if parameters.np_types == 1
        
        # Return elastic stress–strain dee matrix in 2D (plane strain) or
        # 3D. Uses Young’s modulus and Poisson’s ratio.
        dee_d = cu(num.deemat(parameters.young[1], parameters.poisson[1], parameters.nodof, parameters.type_2d, process))

    end
        
    # Output number of equations and skyline storage
    write(output, "There are $(parameters.neq) equations and the skyline storage is $(kdiag[parameters.neq])\n")
    
    println("Start of stiffness and assembly procedures")

    # element stiffness
    #km = zeros(Float32, parameters.ndof, parameters.ndof, parameters.nels)

    for iel = 1:parameters.nels

        # loop for all elements

        # element stiffness
        km_d = CUDA.zeros(Float32, parameters.ndof, parameters.ndof)

        if parameters.np_types > 1
                
            # Return elastic stress–strain dee matrix in 2D (plane strain) or
            # 3D. Uses Young’s modulus and Poisson’s ratio.
            dee2_d = num.para_deemat(parameters.young[parameters.etype[iel]], parameters.poisson[parameters.etype[iel]], parameters.nodof, parameters.type_2d, process)

        end
            
        # loop for all integration points of the current element
            
        write(process, "\n/////////////////  SITFFNESS OF ELEMENT $(iel) //////////////\n")
            
        for i = 1:parameters.nip
                
            # Returns the shape function derivatives der at
            # the ith integrating point (isoparametric coordinates).
            der = num.shapeDer(parameters.points, i, parameters.nod, parameters.element, process)
    
            # bee matrix (shape functions derivatives in cartesian
            #system), jacobian and its determinant.
            # Section 3.7.3 of the FEA book
            jac = cu(der)*cu(parameters.coord[parameters.g_num[:, iel], :])
            write(process, "Jacobian:\n")
            for i in 1:size(jac, 1)
        
                write(process, "$(jac[i, :])\n")
                
            end
            determ = det(Array(jac))
            deriv_d = CuArray{Float32}(undef, (parameters.nodof, parameters.nod))
            deriv_d = cu(inv(Array(jac))*der)
            write(process, "\nnum.para_beemat\n\n")
            @cuda threads = parameters.nod num.para_beemat!(bee_d, deriv_d)
            write(process, "\nbee:\n")
            for i in 1:size(bee_d, 1)
        
                write(process, "$(bee_d[i, :])\n")
                
            end
            write(process, "\nJacobian determinant = $(determ)\n")

            if parameters.type_2d == "axisymmetric"
                    
                # Returns the shape functions fun at the ith integrating point.
                fun = num.shapeFun(parameters.points, i, parameters.nod, parameters.element, process)
                global parameters.gc = fun*parameters.coord
                bee_d[4, 1:2:parameters.nodof*parameters.nod - 1] = fun[:]/parameters.gc[1]
                
            end

            if parameters.np_types > 1

                km_d[:, :] += bee_d'*dee2_d*bee_d*cu(determ)*cu(parameters.weights[i])*cu(parameters.gc[1])
                # km[:, :, iel] += bee'*dee*bee*determ*parameters.weights[i]*parameters.gc[1]
 
            else

                km_d[:, :] += bee_d'dee_d*bee_d*cu(determ)*cu(parameters.weights[i])*cu(parameters.gc[1])
                # km[:, :, iel] += bee'*dee*bee*determ*parameters.weights[i]*parameters.gc[1]

            end
                
        end

        iel == trunc(parameters.nels/2) && println("Halfway through stiffness and assembly procedures")
        
        # Returns lower triangular global stiffness matrix kv stored
        # as a vector in skyline form.
        km = Array(km_d)
        num.fsparv!(kv, km, parameters.g_g[:, iel], kdiag, process)
        
    end

    #
        # for linha in 1:size(parameters.g_g, 1)

        #     linha == 1 && (global nulo = 0)
        #     for coluna in 1:size(parameters.g_g, 2)

        #         parameters.g_g[linha, coluna] == 0 && (nulo += 1)
                
        #     end
            
        # end

        # Inside the kernel num.para_fsparv!, this matrix will
        # store the values to be input later into kv
        #preKv = zeros(Int16, parameters.ndof^2 - parameters.ndof, parameters.nels)

        # Inside the kernel num.para_fsparv!, this matrix will
        # store the respective positions to be altered in kv
        #preKv_indices = zeros(Int16, (parameters.ndof^2 - parameters.ndof, parameters.nels))

        #intermPos = zeros(Int16, size(parameters.g_g, 2))

        #write(process, "\nnum.para_fsparv\n")
        # CUDA.@time @cuda threads = parameters.nels para_fsparv!(cu(preKv_indices), cu(preKv), cu(km), cu(parameters.g_g), cu(kdiag), cu(intermPos))
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
        # for ele in 1:parameters.nels

        #     for DOF in 1:parameters.ndof

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
    kv = num.sparin!(kv, kdiag, process)
    # println("sparin!\n$(kv[1:100])")

    # Find displacements solution
    parameters.disp = num.spabac!(kv, kdiag, parameters.disp, process)
    # println("spabac!\n$(parameters.disp[1:100])")

    return parameters.disp

end


end