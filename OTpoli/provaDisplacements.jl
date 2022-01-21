




push!(LOAD_PATH, ".")
import graph

nxe = 50
nye = 15
g_num = convert.(Int64, g_num)'
coord[:,1] .-= abs(maximum(coord[:,1])-minimum(coord[:,1]))/2
coord[:,2] .-= abs(maximum(coord[:,2])-minimum(coord[:,2]))/2

nf = zeros(Int64, (2,(nxe+1)*(nye+1)))
for j in 1:size(nf,2)
  for i in 1:2
    if j + i == 2
      k = 0
    end
    global k += 1
    nf[i,j] = k
  end  
end

println("\ngraphics")
if true
steps = 1000
# Scale factor for deformations
deformScale = 40
# Scale factor for the geometry
# (Replace with correct "normalized device coordinates")
coordScale = 1e-2
# General setup e shader compilation procedures
window, ProgramID, vertexShader, fragmentShader = graph.glSetupDisp("quadrilateral", 4, nxe, nye, coord, g_num, "Shaders\\vertexDisp.vs", "Shaders\\fragDisp.fs", 2, disp, nf, coordScale)
# Loop with OpenGL functions for drawing
# and processing of mouse and keyboard inputs
graph.renderDisp(steps, ProgramID, window, deformScale, nxe*nye, vertexShader, fragmentShader)
end