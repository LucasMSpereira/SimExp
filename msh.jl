module msh

using Meshes, InOut

### This module provides mesh routines

function triQuadMshData!(nxe, nye, nf, type, spacing, coords)
    
    # Generate coords (nodeID x coordinate), g_num (nodeID x element)
    # and g_g (DOF x element) matrices from mesh
    deltax = -convert(Float64, spacing[1]*nxe/2)
    deltay = -convert(Float64, spacing[2]*nye/2)
    recMesh = Meshes.CartesianGrid((nxe, nye), (deltax, deltay), spacing)
    
    if type == "quadrilateral"
        
        g_num = Array{Int64}(undef, (4, nelements(recMesh)))
        g_g = Array{Int64}(undef, (8, nelements(recMesh)))
        
    elseif type == "triangle"
        
        g_num = Array{Int64}(undef, (3, 2*nelements(recMesh)))
        g_g = Array{Int64}(undef, (6, 2*nelements(recMesh)))
        
    end
    
    for elem in 1:nelements(recMesh)
        
        if (elem - 1)/size(recMesh)[1] == floor(Int16, (elem - 1)/size(recMesh)[1])
            
            # if first element in this row
            
            node1 = floor(Int16, (elem - 1)/size(recMesh)[1]) + elem
            coords[node1, :] = coordinates(vertices(recMesh[elem])[1])
            coords[node1 + 1, :] = coordinates(vertices(recMesh[elem])[2])
            
            if elem >= (size(recMesh)[1]*(size(recMesh)[2] - 1) + 1)
                
                node1 = floor(Int16, (elem - 1)/size(recMesh)[1]) + elem
                node4 = node1 + size(recMesh)[1] + 1
                coords[node4 + 1, :] = coordinates(vertices(recMesh[elem])[3])
                coords[node4, :] = coordinates(vertices(recMesh[elem])[4])
                
            end
            
        else
            
            node1 = floor(Int16, (elem - 1)/size(recMesh)[1]) + elem
            coords[node1 + 1, :] = coordinates(vertices(recMesh[elem])[2])
            
            if elem >= (size(recMesh)[1]*(size(recMesh)[2] - 1) + 1)
                
                node1 = floor(Int16, (elem - 1)/size(recMesh)[1]) + elem
                node3 = node1 + size(recMesh)[1] + 2
                coords[node3, :] = coordinates(vertices(recMesh[elem])[3])
                
            end
            
        end
        
        if type == "quadrilateral"
            
            g_num[1, elem] = floor(Int16, (elem - 1)/size(recMesh)[1]) + elem
            g_num[2, elem] = g_num[1, elem] + 1
            g_num[3, elem] = g_num[1, elem] + size(recMesh)[1] + 2
            g_num[4, elem] = g_num[3, elem] - 1
            
            for node in 1:4
                
                g_g[2*node - 1, elem] = nf[1, g_num[node, elem]]
                g_g[2*node, elem] = nf[2, g_num[node, elem]]
                
            end
            
        elseif type == "triangle"
            
            # separate each "rectangle" in two adjacent triangles
            
            g_num[1, 2*elem - 1] = floor(Int16, (elem - 1)/size(recMesh)[1]) + elem
            g_num[2, 2*elem - 1] = g_num[1, 2*elem - 1] + 1
            g_num[3, 2*elem - 1] = g_num[1, 2*elem - 1] + size(recMesh)[1] + 2
            
            g_num[1, 2*elem] = floor(Int16, (elem - 1)/size(recMesh)[1]) + elem
            g_num[2, 2*elem] = g_num[1, 2*elem] + size(recMesh)[1] + 2
            g_num[3, 2*elem] = g_num[2, 2*elem] - 1
            
            for node in 1:3
                
                g_g[2*node - 1, 2*elem - 1] = nf[1, g_num[node, 2*elem - 1]]
                g_g[2*node, 2*elem - 1] = nf[2, g_num[node, 2*elem - 1]]
                
                g_g[2*node - 1, 2*elem] = nf[1, g_num[node, 2*elem]]
                g_g[2*node, 2*elem] = nf[2, g_num[node, 2*elem]]
                
            end
            
        end
        
    end

    return coords, g_num, g_g
    
end

function RctMshSize(nxe, nye, element)

    # 'Rectangular mesh size' returns nels (number of elements) and
    # nn (number of nodes) in a rectangular mesh

    if element == "quadrilateral"

        nels = nxe*nye

    else element == "triangle"

        nels = 2*nxe*nye

    end

    nn = (nxe + 1)*(nye + 1)

    return nels, nn
    
end

function TriQuadMsh(element, x_coords, y_coords, dir, k, nf)

    # Returns nodal coordinates, IDs and steering vectors for quadrilateral or triangular mesh
    # This routine is equivalent to "geom_rect" from the original FORTRAN code

    coords = Array{Float64}(undef, (length(x_coords) * length(y_coords), 2))
    
    if element == "quadrilateral"

        g_num = Array{Int64}(undef, (4, (length(x_coords) - 1)*(length(y_coords) - 1)))
        g_g = Array{Int64}(undef, (8, (length(x_coords) - 1)*(length(y_coords) - 1)))

    else

        g_num = Array{Int64}(undef, (3, 2*(length(x_coords) - 1)*(length(y_coords) - 1)))
        g_g = Array{Int64}(undef, (6, 2*(length(x_coords) - 1)*(length(y_coords) - 1)))
        
    end       
    
    if dir == "x"

        primary = x_coords
        secondary = y_coords

    elseif dir == "y"

        primary = y_coords
        secondary = x_coords

    else

    end

    for i in 1:(length(secondary))

        for j in 1:(length(primary))

            if i == 1 && j == 1

                n = 1
                
            end
            
            global coords[n, :] =  [primary[j], secondary[i]]
            global n += 1

        end
        
    end

    for i in 1:(length(secondary) - 1)

        for j in 1:(length(primary) - 1)
            
            if element == "quadrilateral" && dir == "x"
                
                global g_num[1, k] = (i - 1) * length(primary) + j
                global g_num[2, k] = g_num[1, k] + length(primary)
                global g_num[3, k] = g_num[2, k] + 1
                global g_num[4, k] = g_num[1, k] + 1
                
            elseif element == "quadrilateral" && dir == "y"

                global g_num[1, k] = (i - 1) * length(primary) + j
                global g_num[2, k] = g_num[1, k] + 1
                global g_num[3, k] = g_num[2, k] + 1
                global g_num[4, k] = g_num[1, k] + length(primary)

            elseif element == "triangle" && dir == "x"

                # Per k quadrilateral in mesh, define 2 adjacent triangles

                o = (i - 1) * length(primary) + j

                global g_num[1, 2 * k - 1] = o + 1
                global g_num[2, 2 * k - 1] = o + length(primary)
                global g_num[3, 2 * k - 1] = o

                global g_num[1, 2 * k] = o + length(primary)
                global g_num[2, 2 * k] = o + 1
                global g_num[3, 2 * k] = o + length(primary) + 1

            elseif element == "triangle" && dir == "y"

                # Per k quadrilateral in mesh, define 2 adjacent triangles

                o = (i - 1) * length(primary) + j

                global g_num[1, 2 * k] = o
                global g_num[2, 2 * k] = o + 1
                global g_num[3, 2 * k] = o + length(primary)
                
                global g_num[1, 2 * k - 1] = o + 1
                global g_num[2, 2 * k - 1] = o + length(primary) + 1
                global g_num[3, 2 * k - 1] = o + length(primary)
                
                for l in 1:3

                    g_g[2 * l - 1, 2 * k - 1] = 2 * g_num[l, 2 * k - 1] - 1
                    g_g[2 * l, 2 * k - 1] = 2 * g_num[l, 2 * k - 1]
                    
                    g_g[2 * l - 1, 2 * k] = 2 * g_num[l, 2 * k] - 1
                    g_g[2 * l, 2 * k] = 2 * g_num[l, 2 * k]
                
                end
                
            end

            k += 1

        end
        
    end
    
    for i in 1:size(g_num, 2)

        for j in 1:size(g_num, 1)
            
            global g_g[2*j - 1, i] = nf[1, g_num[j, i]]
            global g_g[2*j, i] = nf[2, g_num[j, i]]
            
        end
        
    end

    if element == "quadrilateral" && dir == "y"

        coords = circshift(coords, (0, -1))

    end

    return coords, g_g, g_num

end

function hexaSize(nxe, nye, nze)

    # calculate number of nodes in 20-node hexahedral mesh.
    # suppose ordered xyz_coord vectors

    x = (1 + nye)*(1 + nze)*nxe
    y = (1 + nxe)*(1 + nze)*nye
    z = (1 + nye)*(1 + nxe)*nze
    nn = (nxe + 1)*(nye + 1)*(nze + 1) + x + y + z
    nels = nxe*nye*nze

    return nn, nels
    
end

function hexaMesh(x_coords, y_coords, z_coords, nf, nod, nxe, nye, nze, nn)

    # Returns nodal coordinates, IDs and steering vectors for hexahedra ("bricks")

    if nod == 8
    end

    if nod == 14
    end

    if nod == 20
        
        g_num = Array{Int64}(undef, ((nod, ((length(x_coords) - 1)*(length(y_coords) - 1))*(length(z_coords) - 1))))
        g_g = Array{Int64}(undef, (3*nod, (((length(x_coords) - 1)*(length(y_coords) - 1))*(length(z_coords) - 1)))) 
        coords = fill(7777.123, (nn, 3))
        
        # placing all VERTICES in coords
        node_id = 0

        for i in 1:length(x_coords)
            
            i == 1 ? node_id = 1 : false
            i != 1 ? node_id = 2*i - 1 : false
            for j in 1:length(y_coords)
                
                for k in 1:length(z_coords)
                    
                    coords[node_id, :] = [x_coords[i] y_coords[j] z_coords[k]]
                    k != length(z_coords) ? node_id += 3*nxe + 2 : false
                    
                end
                node_id += 1 + 2*nxe + 1 + nxe + (1 + nxe)*nze
                
            end
            
        end
        
        # place in coords: intermidiate nodes that form planes perpendicular to the x axis
        
        for i in 2:(nn - 1)
            
            if coords[i - 1, 1] != 7777.123 && coords[i + 1, 1] != 7777.123
                
                coords[i, :] = (coords[i - 1, :] + coords[i + 1, :])/2
                
            end
            
        end
        
        # place in coords: intermidiate nodes that form planes perpendicular to the y or z axes
        
        for p in 1:(nye + 1)
            
            p == 1 ? node_id = 2*nxe + 2 : false
            for j in 1:nze
                
                for i in 1:(nxe + 1)
                    
                    coords[node_id, :] = [x_coords[i] y_coords[p] (z_coords[j] + z_coords[j + 1])/2]
                    i != (nxe + 1) ? node_id += 1 : false
                    
                end
                node_id += 2*nxe + 2
                
            end
            
            if p < nye + 1
                for j in 1:(nze + 1)
                    
                    for i in 1:(nxe + 1)
                        
                        coords[node_id, :] = [x_coords[i] (y_coords[p] + y_coords[p + 1])/2 z_coords[j]]
                        if i != (nxe + 1) || j != (nze + 1)
                            
                            node_id += 1
                            
                        end
                        
                    end
                    
                end
                node_id += 2*nxe + 2
            end
            
        end
        
        # construct g_num (nodes per element)
        
        for p in 1:nxe
            
            # loop for each nodal plane perpendicular to x
            # see figure 5.22, pg 194 FEA book
            
            n = 2*p - 1
            for i in 1:nye
                
                for j in 1:nze
                    
                    if p == 1 && i == 1 && j == 1
                        
                        global el = 1
                        
                    end
                    middleNodes = 1 + 2*nxe + (1 + nxe + 1 + 2*nxe)*(1 + nze - j) + (j - 1)*(nxe + nxe*(j - 1)) + n
                    backNodes = middleNodes + 1 + nxe + (nze - j + 1)*(1 + nxe) + (j - 1)*(1 + 2*nxe + 1 + nxe)
                    g_num[1, el] = n
                    g_num[2, el] = n + 1
                    g_num[3, el] = n + 2
                    g_num[4, el] = n + 4
                    g_num[5, el] = n + 7
                    g_num[6, el] = n + 6
                    g_num[7, el] = n + 5
                    g_num[8, el] = n + 3
                    g_num[9, el] = middleNodes
                    g_num[10, el] = middleNodes + 1
                    g_num[11, el] = middleNodes + 3
                    g_num[12, el] = middleNodes + 2
                    g_num[13, el] = backNodes
                    g_num[14, el] = backNodes + 1
                    g_num[15, el] = backNodes + 2
                    g_num[16, el] = backNodes + 4
                    g_num[17, el] = backNodes + 7
                    g_num[18, el] = backNodes + 6
                    g_num[19, el] = backNodes + 5
                    g_num[20, el] = backNodes + 3
                    j != nze ? n += 5 : n += 10 + 2*nze
                    el += 1
                    
                end
                
            end
            
        end
        
        # construct g_g (DOFs per element)
        
        for i in 1:size(g_num, 2)
            
            for j in 1:size(g_num, 1)
                
                if j == 1
                    
                    a = 1
                    
                else
                    
                    a = 3*(j - 1) + 1
                    
                end
                
                global g_g[a, i] = nf[1, g_num[j, i]]
                global g_g[a + 1, i] = nf[2, g_num[j, i]]
                global g_g[a + 2, i] = nf[3, g_num[j, i]]
                
            end
            
        end
        
        return coords, g_g, g_num
        
    end

end

end