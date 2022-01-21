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
import graph          # Graphics related routines
using Parameters, CUDA, Printf, Statistics, Plots,
  JuMP, Nonconvex, ChainRulesCore, FiniteDifferences, ChainRulesTestUtils

#@time begin
output = open("results.res", "w")

# Gets data from input file as an array of strings
data = InOut.Extract("Input\\Tarefa3A120x40.dat")
# Loop array of strings to define variables
for i in 1:size(data, 1)

  if data[i] == "type_2d"

    global type_2d = data[i + 1]            # type of strain to be considered

  elseif data[i] == "selfWeight"
    global selfWeight = parse(Bool, data[i+1])
  elseif data[i] == "changeSens"

    global changeSens = parse(Float64, data[i+1]) # Fixed densities update

  elseif data[i] == "SIMPpenalty"
    
    global SIMPpenalty = parse(Float64, data[i+1]) # penalty for SIMP method

  elseif data[i] == "nmodes"

    # Number of vibration modes to be extracted
    global nmodes = parse(Int64, data[i + 1])
    # Vibration mode to be plotted
    global eigenvec = parse(Int64, data[i + 3])

  elseif data[i] == "penalty"
    global penalty = parse(Float64, data[i + 1])
  elseif data[i] == "element"

    global element = data[i + 1]            # type of element

  elseif data[i] == "nod"
    
    global nod = parse(Int64, data[i + 1])           # number of nodes per element

    # 29: DOF/elem = node/elem * DOF/node
    global ndof = nod * nodof

  elseif data[i] == "nxe"

    # apparently has to be Int64 to work with Meshes.CartesianGrid()
    global nxe = parse(Int64, data[i + 1])           # number of elements in x direction

  elseif data[i] == "nye"

    # apparently has to be Int64 to work with Meshes.CartesianGrid()
    global nye = parse(Int64, data[i + 1])           # number of elements in y direction

  elseif data[i] == "nze"

    global nze = parse(Int64, data[i + 1])

  elseif data[i] == "nip"

    global nip = parse(Int64, data[i + 1])           # number of integrating points

  elseif data[i] == "fixed_freedoms"

    global fixed_freedoms = parse(Int64, data[i + 1])              # number of fixed displacements
    
  elseif data[i] == "nodof"

    global nodof = parse(Int64, data[i + 1])                  # number os DOFs per node

  elseif data[i] == "nr"

    global nr = parse(Int64, data[i + 1])              # number of nodes with restrained DOF(s)

    if nodof == 2
      global nels = nxe*nye
      global nn = (nxe + 1)*(nye + 1)
      element == "triangle" && (nels *= 2)
    elseif nodof == 3 && element == "HEXA20"

      global nn, nels = msh.hexaSize(nxe, nye, nze)

    else

      println("Invalid nodof")

    end

    # node freedom matrix (for example mechanical support *see picture in pg 173)
    global nf = ones(Int64, nodof, nn)
    a = b = j = 1

    if nodof == 2

      for c in 1:(nr*(nodof + 1))
      
        if b == 1

          j = parse(Int64, data[i + 2 + c])
          b += 1
          
        elseif b == 2
              
          nf[1, j] = parse(Int64, data[i + 2 + c])
          b += 1
          
        elseif b == 3
              
          nf[2, j] = parse(Int64, data[i + 2 + c])
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
              
          j = parse(Int64, data[i + 2 + c])
          b += 1
          
        elseif b == 2
              
          nf[1, j] = parse(Int64, data[i + 2 + c])
          b += 1
          
        elseif b == 3
              
          nf[2, j] = parse(Int64, data[i + 2 + c])
          b += 1

        elseif b == 4

          nf[3, j] = parse(Int64, data[i + 2 + c])
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
    nf = mtx.formnf(nf)

    # 43: number of equations = number of DOFs
    # without restriction/imposed displacement
    global neq = maximum(nf)

  elseif data[i] == "np_types"

    # Number of types of elements with regards to
    # combination of physical properties
    global np_types = parse(Int64, data[i + 1])
    # Number of physical/mechanical properties
    global nprops = parse(Int64, data[i + 3])
    
    global props = Array{Float64}(undef, (nprops, np_types))
    global n = 4
    for prop in 1:nprops
      for type in 1:np_types
        global props[prop, type] = parse(Float64, data[i + n])
        global n += 1
      end
    end

    #=
    In case physical/mechanical properties vary among
    elements, etype (nels x 1) will associate
    an ID to each element. This ID will later be
    used as an index to access the respective properties.
    For example, if the 23rd component of vector etype is 4, it
    means that element 23 is of "material type" 4. Then, when the algorithm
    needs to know a certain property of that element, say it's Young’s
    modulus, it'll do so by accessing the matrix props in the 4th
    column and the line that corresponds to Young's modulus.
    =#
    if np_types > 1

      global etype = Array{Int64}(undef, nels)
      for ll in 1:nels

        etype[ll] = parse(Int64, data[i + n + 1 + ll])

      end

    else

      global etype = [0]               # constant properties among elements

    end

  elseif data[i] == "loaded_nodes"
    
    global loads = InOut.fillInDisp(data, neq, cg_limit, nodof,
              parse(Int64, data[i + 1]), nf, i+1)

  elseif data[i] == "spacing"

    # Grid spacing in x and y axis for rectangular meshes
    
    global spacing = (parse(Float64, data[i + 1]), parse(Float64, data[i + 2]))

  elseif data[i] == "cg_tol"

    if data[i + 1] != "0"

      global cg_tol = parse(Float64, data[i + 1])        # preconditioned conjugate gradient tolerance
    
    else

      cg_tol = 0
    
    end

  elseif data[i] == "cg_limit"

    global cg_limit = parse(Int64, data[i + 1])        # preconditioned conjugate gradient iterations limit

  elseif data[i] == "GPU"

    global GPU = parse(Bool, data[i + 1])

  elseif data[i] == "probType"

    global probType = data[i + 1]

  elseif data[i] == "etol"

    global etol = parse(Float64, data[i + 1])

  # elseif data[i] == "optIgnore"

  #   # In case of topological optimization, define region(s)
  #   # in which material or void must be guaranteed
  #   #=
  #     optIgnore 3
  #     dist 0 3.5 40.2 5.0
  #     xy 1
  #     40.0 >
  #     -20.0 <
  #     square 0
  #     20.0 20.0
  #     30.0 20.0
  #     30.0 30.0
  #     20.0 30.0
  #   =#
    
  #   global optRules = parse(Int32, data[i+1])
  #   for rule in 1:optRules
  #     rule == 1 && s = 1
  #     if data[i+s+rules] == "dist"
  #       [data[i+s+rules+1]]
  #       s = 5
  #     end
  #   end
  elseif data[i] == "volfrac"
    global volfrac = parse(Float64, data[i+1])
  elseif data[i] == "rmin"
    global rmin = parse(Float64, data[i+1])
  elseif data[i] == "filter"
    global filter = data[i+1]
  end

end


##########################
### GENERATE MESH DATA ###
##########################


# Build global mesh references (relation between elements, node IDs, DOFs)
if element == "HEXA20"

  coord, g_g, g_num = msh.hexaMesh(x_coords, y_coords, z_coords, nf, nod, nxe, nye, nze, nn)

elseif element == "triangle" || element == "quadrilateral"

  coord = Array{Float64}(undef, (nxe + 1)*(nye + 1), 2)
  coord, g_num, g_g = msh.triQuadMshData!(nxe, nye, nf, element, spacing, coord)

end


# Define initial disp vector containing loads
if probType in ["static", "staticSens"]
  gc = ones(nodof)
else
  gc = [1.0]
end
# Include self-weight loads
if selfWeight==true
  loads = num.selfWeight!(nodof, nels, nod, nf, g_num, loads, spacing, props[3,:], ones(nels), etype, np_types)
end

# Coordinates and weighting coefficients for gauss quadrature numerical integration 
points, weights = num.sample(element, nip)

@with_kw mutable struct parameters

  gc::Array{Float64} = gc
  points::Array{Float64} = points   # normal coordinates of gauss quadrature integrating points
  weights::Array{Float64} = weights   # weights for gauss quadrature
  cg_limit::Int64 = cg_limit  # limit of pcg iterations
  g_g::Array{Int64} = g_g   # DOFs x element
  nels::Int64 = nels    # total number of elements
  neq::Int64 = neq    # number of equations
  ndof::Int64 = ndof   # number of DOFs per element
  nip::Int64 = nip  # number of integrating points per element
  element::String = element  # type of element
  # Int64 to avoid error in CuArray initialization (num.para_beemat)
  nod::Int64 = nod    # number of nodes per element  
  coord::Array{Float64} = coord   # list of cartesian coordinates of each node
  g_num::Array{Int64} = g_num   # node IDs x element
  type_2d::String = type_2d   # type of 2D analysis
  # Int64 to avoid error in CuArray initialization (num.para_beemat)
  nodof::Int64 = nodof      # number of DOFs per node
  etype::Array{Int64} = etype   # type of each element, in case physical properties vary
  loads::Array{Float64} = loads   # initialy nodal loads, then nodal displacements
  cg_tol::Float64 = cg_tol    # tolerance threshold for pcg convergence
  props::Array{Float64} = props  # physical properties of each material type of element
  nmodes::Int64 = nmodes # number of modes of vibration for eigenproblems
  etol::Float64 = etol # tolerance for iterative eigenvalue routine
  penalty::Float64 = penalty
  eigenvec::Int64 = eigenvec # vibration mode to be plotted

end

para = parameters()

################################################################
######### CALL SOLVER ACCORDING TO TYPE OF SIMULATION ##########
################################################################

# Static, structural, serialized and with assembly
[probType, iszero(cg_limit), GPU] == ["static", true, false] && (@time disp = num.stiffSerialAssembly(para, output))
# Static, structural, serialized and without assembly (uses PCG)
[probType, iszero(cg_limit), GPU] == ["static", false, false] && (@time disp = num.stiffSerialPCG(para, output))
# Static, structural, parallelized in the GPU and with assembly
if [probType, iszero(cg_limit), GPU] == ["static", true, true]
  CUDA.allowscalar(false)
  @time disp = num.stiffParaAssembly(para, output)
end
# Static, structural, parallelized and without assembly (uses PCG)
if [probType, iszero(cg_limit), GPU] == ["static", false, true]
  CUDA.allowscalar(false)
  @time disp = num.stiffParaPCG(para, output)
end
[probType, iszero(cg_limit), GPU] == ["static", false, true] && (@time disp = num.stiffSerialPCG(para, output))
# Free vibration analysis
[probType, iszero(cg_limit), GPU] == ["eigenSimple", true, false] && (@time disp = num.eigenSimple(para, output))



if false
  # Calculation of sensitivities for static analysis, with or without self-weight
  if probType == "staticSens"
    # Define elements which will be tested
    global eleRange = 1:2:para.nels
    # Define step sizes to be tested
    global stepRange = -4.0:-1:-5.0
    # Loop in step size values
    for changeSens in stepRange
      println("step=10^$changeSens=$(10^changeSens)")
      if changeSens==stepRange[1]
        # Prepare some variables
        global sensAnalitic = zeros(length(eleRange))
        global sensFinDiff = zeros(length(eleRange))
        global densities = fill(volfrac, para.nels)
        global error = zeros(Float64, length(eleRange), length(stepRange))
        # Get compliance for case of unit densities
        global compUni = 1
        para.disp, elementsKm, compUni = num.staticSens(para, densities.^SIMPpenalty)
        @printf "Compliance for unitary densities: %.4e\n" compUni
      end
      @printf "Test point\tElement\t\tCompliance\tAnalytic\tApprox.\t\tDifference (%%)\n"
      global point = 1
      # Loop all elements to be tested
      for iel in eleRange
        # Element density step
        global densities[iel] += 10^changeSens
        # Refill "disp" with initial loads (fills with zeros in case of only self-weight)
        global para.disp = InOut.fillInDisp(data, para.neq, para.cg_limit, para.nodof, loaded_nodes, nf, loadedNodes_pos)
        # Include self-weight loads, in case they should be considered
        if selfWeight==true
          para.disp = num.selfWeight!(nodof, nels, nod, nf, g_num, para.disp, spacing,
          para.props[3,:], densities, etype, np_types)
        end
        # Get compliance for case of element density step
        para.disp, elementsKm, compNew = num.staticSens(para, densities.^SIMPpenalty)
        # Calculate compliance sensitivities analitically
        if selfWeight==true
          # With self-weight
          global sensAnalitic[point] = num.compSensWeight(SIMPpenalty, densities[iel], para, iel, elementsKm[iel], spacing)
        else
          # Without self-weight
          global sensAnalitic[point] = num.compSens(SIMPpenalty, densities[iel], para, iel, elementsKm[iel])
        end
        # Verify compliance sensitivities by finite differences
        global sensFinDiff[point] = num.finDiff(compNew, compUni, 10^changeSens)
        # Error between analitical and numerical sensitivities
        global error[point,findfirst(x->x==changeSens,stepRange)]=(sensFinDiff[point]-sensAnalitic[point])/sensAnalitic[point]*100
        # Print data
        @printf "%i\t\t%i\t\t%.4e\t%.4e\t%.4e\t%0.2f\n" point iel compNew sensAnalitic[point] sensFinDiff[point] error[point,findfirst(x->x==changeSens,stepRange)]
        global densities = fill(volfrac, para.nels)
        global point += 1
      end
    end
    # # Output error to process.res
    # for i in 1:size(error, 1)
    #   for j in 1:size(error, 2)
    #     if j!=size(error, 2)
    #       @printf process "%0.4e " error[i,j]
    #     else
    #       @printf process "%0.4e\n" error[i,j]
    #     end
    #   end
    # end
    println("Best step: $(stepRange[findmin([mean(error[:,i]) for i in 1:length(stepRange)])[2]])")

  end
  
end

if false
  # Calculation of sensitivities for free vibration analysis
  if probType == "eigenSens"
    # Define elements which will be tested
    eleRange = 1:para.nels
    stepRange = -4.0:-0.5:-7.0
    # Loop in all step size values
    for changeSens in stepRange
      if changeSens==stepRange[1]
        # Prepare some variables
        global sensAnalitic = zeros(length(eleRange), para.nmodes)
        global sensFinDiff = zeros(length(eleRange), para.nmodes)
        global densities = ones(para.nels)
        global error = zeros(length(eleRange), length(stepRange), para.nmodes)
        # Get compliance for case of unit densities
        global valsUni=zeros(para.nmodes)
        global mGlobal = zeros(para.neq,para.neq)
        global vecs = zeros(para.neq,para.nmodes)
        valsUni, vecs, elementsKm, elementsMm, mGlobal = num.eigenSens(para, densities, SIMPpenalty, output)
        print("Eigenvalues for unitary densities: ")
        [@printf "%.4e " valsUni[mode] for mode in 1:para.nmodes]
      end
      println("\nstep=10^$changeSens=$(10^changeSens)")
      global point = 1
      # Loop in elements of domain to test sensitivities
      for iel in eleRange
        # Element density step
        densities[iel] += 10^changeSens
        # Get compliance for case of element density step
        valsStep, vecs, elementsKm, elementsMm, mGlobal = num.eigenSens(para, densities, SIMPpenalty, output);
        # Calculate compliance sensitivities analitically
        global sensAnalitic[point,:] = num.eigenValSens(para, vecs, SIMPpenalty, densities[iel], elementsKm[iel], valsStep, elementsMm[iel], iel, mGlobal)
        # Verify compliance sensitivities by finite differences
        global sensFinDiff[point, :] = num.finDiff(valsStep, valsUni, 10^changeSens)
        # Error between analitical and numerical sensitivities
        error[point,findfirst(x->x==changeSens,stepRange),:]=(sensFinDiff[point,:]-sensAnalitic[point,:])./sensAnalitic[point,:]*100
        # Print data
        if iel%2==0
          @printf "Test point: %i\tElement: %i\n" point iel
          @printf "Eigenvalues\t\tAnalytical\t\tApprox.\t\t\tError\n"
          @printf "%.4e\t\t%.4e\t\t%.4e\t\t%0.2f\n" valsStep[1] sensAnalitic[point,1] sensFinDiff[point, 1] error[point,findfirst(x->x==changeSens,stepRange),1]
        end
          densities = ones(Float64, para.nels)
        point += 1
      end
    end
    # Output errors to process.res
    # for k in 1:size(error, 3)
    #   for i in 1:size(error, 1)
    #     for j in 1:size(error, 2)
    #       j!=size(error, 2) ? (@printf process "%0.4e " error[i,j, k]) : (@printf process "%0.4e\n" error[i,j, k])
    #     end
    #   end
    #   @printf process "\n"
    # end
    [println("Best step for mode $mode: $(stepRange[findmin([mean(abs.(error[:,i, mode])) for i in 1:length(stepRange)])[2]])") for mode in 1:para.nmodes]
  end
end






#
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
  # InOut.DispOut(output, type_2d, para.disp, nf, nn, nodof)

#

##################################################
### RECOVER STRESSES AT nip INTEGRATING POINTS ###
##################################################

#num.CentroidStress(para, disp, output)

#println("Time for solver, post-processing and output:")
#end


#####################################
######### GRAPHICAL OUTPUT ##########
#####################################
if false
  println("Beginning of graphics procedures")
  # para.disp=vecs[:,3]
  # Duration of animation
  steps = 600
  # Scale factor for deformations
  deformScale = 0.02
  # Scale factor for the geometry
  # (Replace with correct "normalized device coordinates")
  coordScale = 0.4
  # General setup e shader compilation procedures
  window, ProgramID, vertexShader, fragmentShader = graph.glSetupDisp(element, nod, nxe, nye, coord, g_num, "Shaders\\vertexDisp.vs", "Shaders\\fragDisp.fs", nodof, para.disp, nf, coordScale)
  # Loop with OpenGL functions for drawing
  # and processing of mouse and keyboard inputs
  graph.renderDisp(steps, ProgramID, window, deformScale, nxe*nye, vertexShader, fragmentShader)
  println("End of graphics procedures")
end

close(output)