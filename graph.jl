module graph

using Base: Float32
import GLFW
using ModernGL, LinearAlgebra, Printf, InOut, num, Statistics
include("utilities.jl")

    ### This module provides graphics processing routines

function deformData(nf, disp, steps)
    
    #=
        Create the equivalent of nf (dof/node), but replacing the DOF's IDs
        with their respectives displacements, which are expressed by ranges
        from 0 to their full values. This data will be used to graphically
        plot the displacements in a "linearized" animation.
    =#

    dynNodes = zeros(size(nf)[1], size(nf)[2], steps + 1)
    for node in 1:size(nf, 2)
        
        if nf[1, node] != 0

            if disp[nf[1, node]] != 0
                
                range = 0:disp[nf[1, node]]/steps:disp[nf[1, node]]
                dynNodes[1, node, 1:length(range)] = range

            end

        end
        if nf[2, node] != 0

            if disp[nf[2, node]] != 0
    
                range = 0:disp[nf[2, node]]/steps:disp[nf[2, node]]
                dynNodes[2, node, 1:length(range)] = 0:disp[nf[2, node]]/steps:disp[nf[2, node]]

            end
        end
        
    end
    
    return dynNodes
    
end

function projection(near, far, left, right, top, bottom)

    # Computes the projection matrix for perspective based on the bounds of the view frustum

    p = GLfloat[
        2*near/(right - left) 0 (right + left)/(right - left) 0
        0 2*near/(top - bottom) (top + bottom)/(top - bottom) 0
        0 0 -(far + near)/(far - near) -2*far*near/(far - near)
        0 0 -1.0 0
    ]

    return p
    
end

function getShader(shaderFile)
    
    #=
    This function gets the content of an external text file that
    defines the shader (in GLSL). It will alter
    this string by including on its first line
    (through pushfirst!()) the version of the GLSL to be considered.
    Finally, it will output this string.
    =#
    
    # GLSLversion = get_glsl_version_string() # get version of GLSL
    fileID = open(shaderFile, "r")
    content = readlines(fileID)
    # pushfirst!(content, GLSLversion)
    shader = content[1]
        
    for i in 2:(length(content) - 1)
            
        shader *= "$(content[i])\n"

    end

    shader *= "$(content[length(content)])"
        
    close(fileID)
        
    return shader
        
end

function ProgramSetup(vsh, fsh)
        
    # Compile shaders and check for errors
    vertexShader = createShader(vsh, GL_VERTEX_SHADER)
    # success = GLint[0]
    # glGetShaderiv(vertexShader, GL_COMPILE_STATUS, success)
    # if !success
    #     println("entrou")
    #     logSize = GLint[0]
    #     glGetShaderiv(shader, GL_INFO_LOG_LENGTH, logSize)
    #     infoLogg = Array{String}(undef, logSize)
    #     glGetShaderInfoLog(vertexShader, logSize, l, infoLogg)
    #     println(infoLogg)
    # end
    fragmentShader = createShader(fsh, GL_FRAGMENT_SHADER)

    # Create a "program" to connect both shaders
    program = createShaderProgram(vertexShader, fragmentShader)

    glUseProgram(program)
        
    return program, vertexShader, fragmentShader
        
end

function mMatrix(s, ??, t, i, axis)

    # define scale, rotation and tanslation operations.
    # inputs s, theta and t refer to arrays with ranges of each operation.
    # this funciton will be called every iteration to update the
    # transformation matrix and change the vertices positions accordingly

    scale = GLfloat[
        1.0 0.0 0.0 s[i]
        0.0 1.0 0.0 s[i]
        0.0 0.0 1.0 s[i]
        0.0 0.0 0.0 1.0
    ]

    if axis == "x"

        rotation = GLfloat[
            1.0 0.0 0.0 0.0
            0.0 cos(??[i]) -sin(??[i]) 0.0
            0.0 sin(??[i]) cos(??[i]) 0.0
            0.0 0.0 0.0 1.0
        ]

    elseif axis == "y"

        rotation = GLfloat[
            cos(??[i]) 0.0 sin(??[i]) 0.0
            0.0 1.0 0.0 0.0
            -sin(??[i]) 0.0 cos(??[i]) 0.0
            0.0 0.0 0.0 1.0
        ]

    elseif axis == "z"

        rotation = GLfloat[
            cos(??[i]) -sin(??[i]) 0.0 0.0
            sin(??[i]) cos(??[i]) 0.0 0.0
            0.0 0.0 1.0 0.0
            0.0 0.0 0.0 1.0
        ]

    end

    translation = GLfloat[
        1.0 0.0 0.0 t[i]
        0.0 1.0 0.0 0.0
        0.0 0.0 1.0 0.0
        0.0 0.0 0.0 1.0
    ]

    return scale*rotation*translation
    
end

function vMatrix()

    # Computes view matrix to map world coordinates onto view/camera space

    orientation = GLfloat[
        1.0 0.0 0.0 0.0
        0.0 1.0 0.0 0.0
        0.0 0.0 -1.0 0.0
        0.0 0.0 0.0 1.0
    ]

    translation = GLfloat[
        1.0 0.0 0.0 0.0
        0.0 1.0 0.0 0
        0.0 0.0 1.0 -1
        0.0 0.0 0.0 1.0
    ]

    return orientation*translation
    
end

function glSetupDisp(element, nod, nxe, nye, coords, g_num, vertexShader, fragmentShader, nodof, disp, nf, coordScale)
    
    # Many setup steps for OpenGL
    
    vertices = Array{GLfloat}(undef, (nxe + 1)*(nye + 1)*4)
    indices = Array{GLuint}(undef, 3*2*nxe*nye)
    coords *= coordScale
    for i in 0:(size(nf, 2) - 1)
        vertices[4*i + 1] = coords[i + 1, 1]
        vertices[4*i + 2] = coords[i + 1, 2]
        nf[1, i + 1] == 0 ? vertices[4*i + 3] = 0 : (vertices[4*i + 3] = disp[nf[1, i + 1]])
        nf[2, i + 1] == 0 ? vertices[4*i + 4] = 0 : (vertices[4*i + 4] = disp[nf[2, i + 1]])
    end
    maxDisp = 0
    for i in 0:(convert(Int32, length(vertices)/4) - 1)
        l = sqrt(vertices[4*i + 3]^2 + vertices[4*i + 4]^2)
        l > maxDisp && (maxDisp = l)
    end

    # Organize indices for vertices according
    # to type of element and number of nodes per element
    if element == "quadrilateral"
        if nod == 4
            for i in 1:size(g_num, 2)
                i == 1 && (a = 1)
                indices[a] = g_num[1, i] - 1
                indices[a+1] = g_num[2, i] - 1
                indices[a+2] = g_num[3, i] - 1
                indices[a+3] = g_num[1, i] - 1
                indices[a+4] = g_num[3, i] - 1
                indices[a+5] = g_num[4, i] - 1
                global a += 6
            end
        end
    elseif element == "triangle"
        if nod == 3
            for i in 1:size(g_num, 2)
                i == 1 && (a = 1)
                indices[a] = g_num[1, i] - 1
                indices[a+1] = g_num[2, i] - 1
                indices[a+2] = g_num[3, i] - 1
                a += 3
            end
        end
    end
    GLFW.Init()
    
    # Create a windowed mode window and its OpenGL context
    window = GLFW.CreateWindow(900, 900, "Displacements")
    # Make the window's context current
    GLFW.MakeContextCurrent(window)

    ### generate buffers and objects ###

    # element buffer
    ebo = glGenBuffer()

    # Generate a vertex array object and array buffer for the data
    vao = glGenVertexArray()
    glBindVertexArray(vao)

    # Create virtual buffer object so data can be copied to GPU
    vbo = glGenBuffer()

    # Bind virtual objects to respective buffers to activate them
    glBindBuffer(GL_ARRAY_BUFFER, vbo)
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo)

    # Copy information to respective buffer
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices,
    GL_STATIC_DRAW)

    glEnable(GL_DEPTH_TEST)

    # shaders access and setup
    vsh = graph.getShader(vertexShader)
    fsh = graph.getShader(fragmentShader)

    # creating OpenGL "program" and linking shaders
    ProgramID, vertexShader, fragmentShader = ProgramSetup(vsh, fsh)

    glUniform1f(glGetUniformLocation(ProgramID, "maxDisp"), maxDisp)

    # position attribute
    glVertexAttribPointer(0, nodof, GL_FLOAT, GL_FALSE, 2*nodof*sizeof(Float32), C_NULL)
    glEnableVertexAttribArray(0)

    # displacement attribute
    glVertexAttribPointer(1, nodof, GL_FLOAT, GL_FALSE, 2*nodof*sizeof(Float32), Ptr{Nothing}(nodof*sizeof(Float32)))
    glEnableVertexAttribArray(1)

    # wireframe mode: GL_LINE
    # fill mode: GL_FILL
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)

    glDisable(GL_LINE_SMOOTH)
    glLineWidth(3)

    return window, ProgramID, vertexShader, fragmentShader

end

function renderDisp(steps, ProgramID, window, deformScale, nels, vertexShader, fragmentShader)

    # Starts the render loop and then terminates GLFW

    glUniform1i(glGetUniformLocation(ProgramID, "stepsTot"), steps)
    glUniform1f(glGetUniformLocation(ProgramID, "deformScale"), deformScale)

    mMatrix = GLfloat[
        3 0 0 0
        0 3 0 0
        0 0 3 0
        0 0 0 1
    ]
    glUniformMatrix4fv(glGetUniformLocation(ProgramID, "model"), 1, GL_FALSE, mMatrix)
    # Starts the render loop and then terminates GLFW
    # while (!GLFW.WindowShouldClose(window))
    for s in 1:steps

        s%20 == 0 && println("$(floor(Int32, s/steps*100))%")
        s == 1 && (p = 0)
        # Set background color
        glClearColor(0.2, 0.2, 0.5, 1)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        glUseProgram(ProgramID)
        if s > 0.95*steps
            
        elseif s > 0.15*steps && s < 0.9*steps
            global p += 1
        end
        glUniform1i(glGetUniformLocation(ProgramID, "stepsCurrent"), p)
        #=
        In graphical procedures, mesh is considered as a subdivision of
        quadrilateral elements into two triangles each.
        (n?? of elements)*(2 triangles/element)*(3 vertices/triangle) = total number of vertices
        =#
        glDrawElements(GL_TRIANGLES, nels*2*3, GL_UNSIGNED_INT, C_NULL)
        # check and call events
        GLFW.PollEvents()
        # swap the buffers
        GLFW.SwapBuffers(window)

    end
    glDeleteShader(vertexShader)
    glDeleteShader(fragmentShader)
    GLFW.Terminate()

end

function renderDens(ProgramID, window, nels, vertexShader, fragmentShader)

    # Starts the render loop and then terminates GLFW
    mMatrix = GLfloat[
        3 0 0 0
        0 3 0 0
        0 0 3 0
        0 0 0 1
    ]
    glUniformMatrix4fv(glGetUniformLocation(ProgramID, "model"), 1, GL_FALSE, mMatrix)

    while (!GLFW.WindowShouldClose(window))

        # Set background color
        glClearColor(0.1, 0.2, 0.3, 1)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        glUseProgram(ProgramID)
        
        #=
        In graphical procedures, mesh is considered as a subdivision of
        quadrilateral elements into two triangles each.
        (n?? of elements)*(2 triangles/element)*(3 vertices/triangle) = total number of vertices
        =#
        glDrawArrays(GL_POINT, 0, nels)
        # check and call events
        GLFW.PollEvents()
        # swap the buffers
        GLFW.SwapBuffers(window)

    end
    glDeleteShader(vertexShader)
    glDeleteShader(fragmentShader)
    GLFW.Terminate()

end

function glSetupDens(nxe, nye, coord, vertexShader, fragmentShader, nodof, dens, g_num, element)
    
    # Many setup steps for OpenGL
    
    vertices = Array{GLfloat}(undef, nxe*nye*3)
    coord /= maximum(coord)

    for i in 0:(nxe*nye - 1)
        vertices[3*i + 1], vertices[3*i + 2] = num.centerPos(element, coord, g_num, i+1)
        vertices[3*i + 3] = convert(Float32, dens[i + 1])
    end
    GLFW.Init()
    
    size=900
    # Create a windowed mode window and its OpenGL context
    window = GLFW.CreateWindow(size, size, "Resulting densities")
    # Make the window's context current
    GLFW.MakeContextCurrent(window)

    ### generate buffers and objects ###

    # Generate a vertex array object and array buffer for the data
    vao = glGenVertexArray()

    # Create virtual buffer object so data can be copied to GPU
    vbo = glGenBuffer()

    # Bind virtual objects to respective buffers to activate them
    glBindBuffer(GL_ARRAY_BUFFER, vbo)
    glBindVertexArray(vao)

    # Copy information to respective buffer
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW)

    glEnable(GL_DEPTH_TEST)

    # shaders access and setup
    vsh = graph.getShader(vertexShader)
    fsh = graph.getShader(fragmentShader)

    # creating OpenGL "program" and linking shaders
    ProgramID, vertexShader, fragmentShader = ProgramSetup(vsh, fsh)

    # position attribute
    glVertexAttribPointer(0, nodof, GL_FLOAT, GL_FALSE, (nodof+1)*sizeof(Float32), C_NULL)
    glEnableVertexAttribArray(0)

    # density attribute
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(Float32), Ptr{Nothing}(nodof*sizeof(Float32)))
    glEnableVertexAttribArray(1)

    # glPointSize(floor(0.9*size/(nxe>nye ? nxe : nye)))
    glPointSize(20)

    return window, ProgramID, vertexShader, fragmentShader

end

end

################################ OLD STUFF FOR FUTURE REFERENCE ################################



    # printArray(coords, "coords")
    # printArray(g_num, "g_num")
    # printArray(g_g, "g_g")

    # Vertices coordinates

        # data = GLfloat[0.5, -0.5, 0.0,
        #               -0.5, -0.5, 0.0,
        #                0.0, 0.5, 0.0]

        # data = GLfloat[0.5, -0.5, 0.0,
        # -0.5, -0.5, 0.0,
        # 0.0, 0.5, 0.0]

        # data = GLfloat[
        # 0.5, 0.0, 0.0, 0.0, # 0a
        # 0.5, 0.5, 0.0, 0.0, # 1a
        # 0, 0.5, 0.0, 0.0, # 2a
        # 0, -0.5, 0.0, 0.0, # 3b
        # -0.5, -0.5, 0.0, 0.0, # 4b
        # -0.5, 0, 0.0, 0.0 # 5b
        # ]

        # indices = GLuint[  # note that we start from 0!
        # 0, 1, 2, # first triangle
        # 3, 4, 5 # second triangle
        # ]

        # data = GLfloat[
        #     -0.5, -0.5, -0.5,
        #      0.5, -0.5, -0.5,
        #      0.5,  0.5, -0.5,
        #      0.5,  0.5, -0.5,
        #     -0.5,  0.5, -0.5,
        #     -0.5, -0.5, -0.5,

        #     -0.5, -0.5,  0.5,
        #      0.5, -0.5,  0.5,
        #      0.5,  0.5,  0.5,
        #      0.5,  0.5,  0.5,
        #     -0.5,  0.5,  0.5,
        #     -0.5, -0.5,  0.5,

        #     -0.5,  0.5,  0.5,
        #     -0.5,  0.5, -0.5,
        #     -0.5, -0.5, -0.5,
        #     -0.5, -0.5, -0.5,
        #     -0.5, -0.5,  0.5,
        #     -0.5,  0.5,  0.5,

        #      0.5,  0.5,  0.5,
        #      0.5,  0.5, -0.5,
        #      0.5, -0.5, -0.5,
        #      0.5, -0.5, -0.5,
        #      0.5, -0.5,  0.5,
        #      0.5,  0.5,  0.5,

        #     -0.5, -0.5, -0.5,
        #      0.5, -0.5, -0.5,
        #      0.5, -0.5,  0.5,
        #      0.5, -0.5,  0.5,
        #     -0.5, -0.5,  0.5,
        #     -0.5, -0.5, -0.5,

        #     -0.5,  0.5, -0.5,
        #      0.5,  0.5, -0.5,
        #      0.5,  0.5,  0.5,
        #      0.5,  0.5,  0.5,
        #     -0.5,  0.5,  0.5,
        #     -0.5,  0.5, -0.5,
        # ]

        # data = GLfloat[
        # 0.5, 0.5, 0.2, 0.5, # top right
        # 0.5, -0.5, 0.4, -0.2, # bottom right
        # -0.5, -0.5, 0.3, -0.3, # bottom left
        # -0.5, 0.5, -0.4, 0.1 # top left
        # ]
    #

    # #scale
        # si  = 0.3
        # sf = 2
        # s = si:(sf - si)/steps:sf
        # # rotation
        # ??i = 0
        # ??f = pi
        # ?? = ??i:(??f - ??i)/steps:??f
        # i = j = 1
        # # view matrix construction
        # viewMatrix = vMatrix()
        # # view frustum for projection matrix
        # near = 1.0
        # far = 10
        # left = -0.5
        # right = 0.5
        # top = 0.5
        # bottom = -0.5
        # vecVaria = Array{GLfloat}(undef, 2)


    # transformations
    # modelMatrix = mMatrix(s, ??, t, i, "x")
    # projMatrix = projection(near, far, left, right, top, bottom)
    # printArray(projMatrix,"projMatrix")


    # transformation uniforms
    # modLoc = glGetUniformLocation(ProgramID, "model")
    # viewLoc = glGetUniformLocation(ProgramID, "view")
    # projLoc = glGetUniformLocation(ProgramID, "projection")
    # glUniformMatrix4fv(modLoc, 1, GL_FALSE, modelMatrix)
    # glUniformMatrix4fv(viewLoc, 1, GL_FALSE, viewMatrix)
    # glUniformMatrix4fv(projLoc, 1, GL_FALSE, projMatrix)


    ### data for transformations    
    # steps = 100
    # # translation
    # ti = 0
    # tf = 0.3
    # t = ti:(tf - ti)/steps:tf

    # for p in 1:3

    #     if p == 1

    #     modelMatrix = GLfloat[
    #         1.0 0.0 0.0 t[i]
    #         0.0 1.0 0.0 0.0
    #         0.0 0.0 1.0 0.0
    #         0.0 0.0 0.0 1.0
    #     ]

    #     elseif p == 2

    #         modelMatrix = GLfloat[
    #         1.0 0.0 0.0 0.0
    #         0.0 1.0 0.0 t[i]
    #         0.0 0.0 1.0 0.0
    #         0.0 0.0 0.0 1.0
    #     ]
        
    #     elseif p == 3

    #         modelMatrix = GLfloat[
    #         1.0 0.0 0.0 0.0
    #         0.0 1.0 0.0 0.0
    #         0.0 0.0 1.0 0.0
    #         0.0 0.0 0.0 1.0
    #     ]

    #     end

    #     modLoc = glGetUniformLocation(ProgramID, "model")
    #     glUniformMatrix4fv(modLoc, 1, GL_FALSE, modelMatrix)
    #     glDrawElements(GL_POINTS, 3, GL_UNSIGNED_INT, Ptr{Nothing}((21 + p)*sizeof(Float32)))
        
    # end

    # for node in 1:sizeMesh
                
        #         modelMatrix = GLfloat[
        #             1 0 0 dynNodes[1, node, s]
        #             0 1 0 dynNodes[2, node, s]
        #             0 0 1 0
        #             0 0 0 1
        #         ]
        #         modLoc = glGetUniformLocation(ProgramID, "model")
        #         glUniformMatrix4fv(modLoc, 1, GL_FALSE, modelMatrix)
        #         # glDrawElements(GL_LINES, 9, GL_UNSIGNED_INT, Ptr{Nothing}(2*(node - 1)*sizeof(Float32)))
    # end

        #@printf "%d\t%.4e\t%.4e\t%.4e\t%.4e\n" s dynNodes[1, 8, s] dynNodes[2, 8, s] dynNodes[1, 9, s] dynNodes[2, 9, s]
    

        # for node in 1:sizeMesh
                
        #     modelMatrix = GLfloat[
        #         1 0 0 dynNodes[1, node, s]
        #         0 1 0 dynNodes[2, node, s]
        #         0 0 1 0
        #         0 0 0 1
        #     ]
        #     modLoc = glGetUniformLocation(ProgramID, "model")
        #     glUniformMatrix4fv(modLoc, 1, GL_FALSE, modelMatrix)
        #     # glDrawElements(GL_LINES, 2, GL_UNSIGNED_INT, Ptr{Nothing}(2*(node - 1)*sizeof(Float32)))
        #     glDrawArrays(GL_POINTS, node - 1, 1)

        # end

        #@printf "%d\t%.4e\t%.4e\t%.4e\t%.4e\n" s dynNodes[1, 8, s] dynNodes[2, 8, s] dynNodes[1, 9, s] dynNodes[2, 9, s]


#
