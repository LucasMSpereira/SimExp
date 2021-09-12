# Started writing in 07/12/20
# Lucas Pereira, SP - São Paulo, Brazil

# Picking essencial steps from the FEA book program 5.1, pg 170
# "FEA book" refers to "Programming the Finite Element Method", I. M. Smith et al, 5 ed., Wiley

# Initial logic based on fortran code provided in the book's
# website (can be obtained without access to the book).

push!(LOAD_PATH, ".")
import InOut          # IO module
import msh            # Mesh routines module
import num            # FEA specific numerical routines
import mtx            # Matrix processing routines
import graph       # Graphics related routines
using Parameters
using CUDA

@time begin
output = open("results.res", "w")
process = open("process.res", "w")

# Gets data from input file as an array of strings
data = InOut.Extract("Input\\p51_1.dat", process)
for i in 1:size(data, 1)

  if data[i] == "type_2d"

    global type_2d = data[i + 1]            # type of strain to be considered

  elseif data[i] == "element"

    global element = data[i + 1]            # type of element

  elseif data[i] == "nod"
    
    global nod = parse(Int32, data[i + 1])           # number of nodes per element

    # 29: DOF/elem = node/elem * DOF/node
    global ndof = nod * nodof

  elseif data[i] == "nxe"

    # apparently has to be Int64 to work with Meshes.CartesianGrid()
    global nxe = parse(Int64, data[i + 1])           # number of elements in x direction

  elseif data[i] == "nye"

    # apparently has to be Int64 to work with Meshes.CartesianGrid()
    global nye = parse(Int64, data[i + 1])           # number of elements in y direction

  elseif data[i] == "nze"

    global nze = parse(Int32, data[i + 1])

  elseif data[i] == "nip"

    global nip = parse(Int32, data[i + 1])           # number of integrating points

  elseif data[i] == "np_types"

    global np_types = parse(Int32, data[i + 1])      # number of variations of physical properties

  elseif data[i] == "young"

    global young = Array{Float32}(undef, np_types)
    for ll in 1:np_types

      young[ll] = parse(Float32, data[i + ll])

    end

  elseif data[i] == "poisson"

    global poisson = Array{Float32}(undef, np_types)
    for ll in 1:np_types

      poisson[ll] = parse(Float32, data[i + ll])

    end

  elseif data[i] == "fixed_freedoms"

    global fixed_freedoms = parse(Int32, data[i + 1])              # number of fixed displacements
    
  elseif data[i] == "nodof"

    global nodof = parse(Int32, data[i + 1])                  # number os DOFs per node

  elseif data[i] == "nr"

    global nr = parse(Int32, data[i + 1])              # number of nodes with restrained DOF(s)

    if nodof == 2

      global nn = (nxe + 1)*(nye + 1)
      global nels = 2*nxe*nye

    elseif nodof == 3 && element == "HEXA20"

      global nn, nels = msh.hexaSize(nxe, nye, nze, process)

    else

      println("Invalid nodof")

    end

    # node freedom matrix (for example mechanical support *see picture in pg 173)
    global nf = ones(Int32, nodof, nn)
    a = b = j = 1

    if nodof == 2

      for c in 1:(nr*(nodof + 1))
      
        if b == 1

          j = parse(Int32, data[i + 2 + c])
          b += 1
          
        elseif b == 2
              
          nf[1, j] = parse(Int32, data[i + 2 + c])
          b += 1
          
        elseif b == 3
              
          nf[2, j] = parse(Int32, data[i + 2 + c])
          b += 1
          
        end
          
        if b == 4
              
          b = 1
              
        end
        c += 1
        continue

      end

    elseif nodof == 3

      for c in 1:(nr*(nodof + 1))
      
        if b == 1
              
          j = parse(Int32, data[i + 2 + c])
          b += 1
          
        elseif b == 2
              
          nf[1, j] = parse(Int32, data[i + 2 + c])
          b += 1
          
        elseif b == 3
              
          nf[2, j] = parse(Int32, data[i + 2 + c])
          b += 1

        elseif b == 4

          nf[3, j] = parse(Int32, data[i + 2 + c])
          b += 1
          
        end
          
        if b == 5
              
          b = 1
              
        end
        
        c += 1
        continue

      end
        
    else

      println("Invalid nodof input. Has to be either 2 or 3.")

    end
    nf = mtx.formnf(nf, process)

    # 43: number of equations = number of DOFs
    # without restriction/imposed displacement
    global neq = maximum(nf)
    write(process, "Number of equations = $(neq)\n")
    
  elseif data[i] == "etype"

    #=
    In case physical/mechanical properties vary among
    elements (that is, in the physical body being simulated), etype (nels x 1)
    will associate an ID to each element. This ID will later be
    used as an index to access the vector of the respective property.
    For example, if the 23rd component of vector etype is 4, it
    means that element 23 is of "type" 4. Then, when the algorithm
    needs to know a certain property of that element, say it's Young’s
    modulus, it'll do so by accessing the vector young (nels x 1)
    in the 4th component.
    =#

    if data[i + 1] != "0"

      global etype = Array{Int32}(undef, nels)        # element property type vector
      for ll in 1:nels

        etype[ll] = parse(Int32, data[i + ll])

      end

    else

      global etype = [0]               # constant properties among elements

    end

  elseif data[i] == "loaded_nodes"
    
    global loaded_nodes = parse(Int32, data[i + 1])           # number of loaded nodes
    global loadedNodes_pos = i

  elseif data[i] == "spacing"

    # Grid spacing in x and y axis for rectangular meshes
    
    global spacing = (parse(Float64, data[i + 1]), parse(Float64, data[i + 2]))

  elseif data[i] == "cg_tol"

    if data[i + 1] != "0"

      global cg_tol = parse(Float32, data[i + 1])        # preconditioned conjugate gradient tolerance
    
    else

      cg_tol = 0
    
    end

  elseif data[i] == "cg_limit"

    global cg_limit = parse(Int32, data[i + 1])        # preconditioned conjugate gradient iterations limit

  elseif data[i] == "GPU"

    global GPU = parse(Int32, data[i + 1])

  end

end



################################################
### LOOP ELEMENTS TO FIND GLOBAL ARRAY SIZES ###
################################################



# Build global mesh references (relation between elements, node IDs, DOFs)
# for later use.
if element == "HEXA20"

  coord, g_g, g_num = msh.hexaMesh(x_coords, y_coords, z_coords, nf, nod, nxe, nye, nze, nn, process)

elseif element == "triangle"

  coord = Array{Float32}(undef, (nxe + 1)*(nye + 1), 2)
  coord, g_num, g_g = msh.triQuadMshData!(nxe, nye, nf, element, spacing, coord, process)

end



##################################################
######### ELEMENT STIFFNESS INTEGRATION ##########
##################################################



write(process, "\n/////////////////////////////////////////////\n")
write(process, "\n     ELEMENTS STIFFNESS\n")
write(process, "\n/////////////////////////////////////////////\n")

disp = zeros(neq, cg_limit + 1)               # disp initially stores applied loads

for vv in 1:(nodof + 1):(loaded_nodes*(nodof + 1))
  
  for j in 1:size(nf, 1)
    
    if nf[j, parse(Int32, data[loadedNodes_pos + 2 + vv])] != 0
      
      global disp[nf[j, parse(Int32, data[loadedNodes_pos + 2 + vv])], 1] = parse(Float32, data[loadedNodes_pos + 2 + vv + j])
      
    end
    
  end
  
end

gc = ones(nodof)

# Coordinates and weighting coefficients for gauss quadrature numerical integration 
points, weights = num.sample(element, nip, process)

@with_kw mutable struct parameters

  gc::Array{Float32} = gc
  points::Array{Float32} = points   # normal coordinates of gauss quadrature integrating points
  weights::Array{Float32} = weights   # weights for gauss quadrature
  cg_limit::Int32 = cg_limit  # limit of pcg iterations
  g_g::Array{Int32} = g_g   # DOFs x element
  nels::Int32 = nels    # total number of elements
  neq::Int32 = neq    # number of equations
  ndof::Int32 = ndof   # number of DOFs per element
  nip::Int32 = nip  # number of integrating points per element
  element::String = element  # type of element
  # Int64 to avoid error in CuArray initialization (num.para_beemat)
  nod::Int64 = nod    # number of nodes per element  
  coord::Array{Float32} = coord   # list of cartesian coordinates of each node
  g_num::Array{Int32} = g_num   # node IDs x element
  type_2d::String = type_2d   # type of 2D analysis
  # Int64 to avoid error in CuArray initialization (num.para_beemat)
  nodof::Int64 = nodof      # number of DOFs per node
  young::Array{Float32} = young   # young's modulus
  poisson::Array{Float32} = poisson   # Poisson's coefficient
  etype::Array{Int32} = etype   # type of each element, in case physical properties vary
  disp::Array{Float32} = disp   # initialy nodal loads, then nodal displacements
  cg_tol::Float32 = cg_tol    # tolerance threshold for pcg convergence
  np_types::Int32 = np_types  # number of variations of physical properties

end

para = parameters()

# Pick solution method according to input.

if cg_limit == 0 && GPU == 0

  # Serial and with assembly

  @time disp = num.stiffSerialAssembly(para, process, output)

elseif cg_limit == 0 && GPU == 1

  # Parallelized and with assembly

  CUDA.allowscalar(false)
  println("disp")
  @time disp = num.stiffParaAssembly(para, process, output)

elseif cg_limit > 1 && GPU == 0

  # Serial and without assembly (uses PCG)

  @time disp = num.stiffSerialPCG(para, process, output)

elseif cg_limit > 1 && GPU == 1

  # Parallelized and without assembly (uses PCG)

  CUDA.allowscalar(false)
  @time disp = num.stiffParaPCG(para, process, output)

else

  println("Check input for assembly strategy and/or parallelization")

end

# 91-100: Include fixed displacements boundary conditions (in case there are any)
# if fixed_freedoms != 0

#   for i = 1:fixed_freedoms

#     # no: fixed DOF numbers vector, nf: nodal freedom matrix
#     no[i] = nf[sense[i],node[i]]

#   end

#   # 98: penalty strategy to fix DOF value, see section 3.6 of the FEA book
#   kv(kdiag(no))=kv(kdiag(no))+penalty 
#   loads(no)=kv(kdiag(no))*value
# end

# 110: write displacements line by line
# as in "[node id]" "[x- disp]" "[y- disp]"
InOut.DispOut(output, type_2d, disp, nf, nn, nodof, process)



##################################################
### RECOVER STRESSES AT nip INTEGRATING POINTS ###
##################################################



write(process, "\n/////////////////////////////////////////////\n")
write(process, "\nCalculating stresses at integrating points\n")
write(process, "\n/////////////////////////////////////////////\n")

if np_types == 1

  # Return elastic stress–strain dee matrix in 2D (plane strain) or
  # 3D. Uses Young’s modulus and Poisson’s ratio.
  dee = num.deemat(young[1], poisson[1], nodof, type_2d, process)

end

InOut.HeaderStress(type_2d, nodof, process, output)
for iel = 1:nels

  eld = zeros(size(g_g, 1))
  for l in 1:size(g_g, 1)

    if g_g[l, iel] != 0
    
      global eld[l] = disp[g_g[l, iel]]

    end
    
  end
  
  if np_types > 1

    # Return elastic stress–strain dee matrix in 2D (plane strain) or
    # 3D. Uses Young’s modulus and Poisson’s ratio.
    dee2 = num.deemat(young[etype[iel]], poisson[etype[iel]], nodof, type_2d, process)

  end

  for i = 1:nip

    write(process, "\nIntegrating point $(nip) of element $(iel)\n")

    der = num.shapeDer(points, i, nod, element, process)
    bee, determ = num.beemat(der, coord[g_num[:, iel], :], type_2d, nodof, nod, process)
    
    if type_2d == "axisymmetric"

      fun = num.shapeFun(points, i, nod, element, process)
      global gc = fun*coord
      bee[4, 1:2:nodof*nod-1] = fun[:]/gc[1]

    end

    if np_types > 1
      
      sigma = dee2*(bee*eld)

    else
    
      sigma = dee*(bee*eld)

    end

    InOut.StressOut(iel, gc, sigma, i, nip, process, output)

  end

end
println("Time until end of simulation and output:")
end


#####################################
######### GRAPHICAL OUTPUT ##########
#####################################

println("Beginning of graphics procedures")
steps = 2000
deformScale = 2.5e5
window, ProgramID = graph.openGLSetup(nxe, nye, coord, g_num, "Shaders\\SimExpVertex.vs", "Shaders\\SimExpFragment.fs", nodof, disp, nf)
graph.renderLoop(steps, ProgramID, window, deformScale)
println("End of graphics procedures")

close(output)
close(process)