module InOut

using Printf, Statistics


### This module provides IO routines

function Extract(file)

    # Returns single array containing strings with the data in input file

    data = String[]
    open(file) do input
        fileContent = readlines(input)
        for line in 1:length(fileContent)

            pieces = split(fileContent[line])
            for j in 1:length(pieces)
            
                push!(data, pieces[j])
    
            end

        end
        
    end
    println("InOut.Extract: length(data) = $(length(data))")
    return data
    
end

function neqStorageOut(neq, store, output)
    # Output number of equations and skyline storage
    write(output, "There are $neq equations and the skyline storage is $store.\n")
end

function DispOut(output, type_2d, disp, nf, nn, nodof)

    # Outputs displacements


    if nodof == 2

        if type_2d == "axisymmetric"

            @printf output "Node\tr-disp\t\t\tz-disp\n\n"
        
        else
            
            @printf output "Node\tx-disp\t\t\ty-disp\n\n"
        
        end

    elseif nodof == 3

        @printf output "\nNode\tx-disp\t\t\ty-disp\t\t\tz-disp\n\n"

    else

        write(results, "Invalid nodof (InOut.DispOut)")

    end

    for k = 1:nn 

        @printf output "%i" k
        for i in 1:nodof

            if nf[i, k] != 0

                @printf output "\t\t%.4e" disp[nf[i, k], end]
                i == nodof ? (@printf output "\n") : false

            else

                @printf output "\t\t0.0000    "
                i == nodof ? (@printf output "\n") : false
                
            end
            
        end
      
    end
    
end

function printArray(variable, varName)

    # Function mainly for debugging purposes

    println("$(varName):")
    for i in 1:size(variable, 1)

        println("$(i)   $(variable[i, :])")
    
    end
    println()
    
end

function StressOut(iel, gc, sigma, i, nip, output)

    # Outputs Stress values

    @printf output "\t\t%i\t\t\t\t%i\t\t%i" iel i gc[1]

    for mm in 1:size(sigma, 1)

        @printf output "\t\t%.4e" sigma[mm]
        mm == size(sigma, 1) && (@printf output "\n")
        
    end

    i == nip && (@printf output "\n")
    
end

function HeaderStress(type_2d, nodof, output)

    # Outputs header of stress table

    if nodof == 2

        if type_2d == "axisymmetric"

            @printf output "\nElement r-coord\t\tz-coord\t\tsig_r\t\tsig_z\t\ttau_rz\t\tsig_t\n\n"
            
        else
        
            @printf output "\nElement\t\tnip\t\tgc\t\t  sig_x\t\t\t   sig_y\t\t\t tau_xy\n\n"
        
        end

    else

        @printf output "\nElement\t\tnip\t\tgc\t\t  sig_x\t\t\t   sig_y\t\t   sig_z\t\t tau_xy\t\t   tau_yz\t\t   tau_zx\n\n"

    end
    
end

function bandwidthOut(nband, output)
    # Outputs bandwidth
    write(output, "\nThe half-bandwidth (including diagonal) is $(nband+1).\n")
end

function EigenValOut(diag, output)
    write(output, "The eigenvalues are: ")
    for i = 1:length(diag)
      @printf output "%.4e  " diag[i]
    end
    @printf output "\n"

end

function EigenVecOut(i, udiag, output)
    write(output, "Eigenvector number $i is: ")
    m = maximum(abs.(udiag))
    for j = 1:length(udiag)
        @printf output "%.4e  " udiag[j]/m
    end
    @printf output "\n"
end

function fillInDisp(data, neq, cg_limit, nodof, loaded_nodes, nf, loadedNodes_pos)
    # disp initially stores applied loads
    loads = zeros(neq, cg_limit + 1)
    
    for vv in 1:(nodof + 1):(loaded_nodes*(nodof + 1))
        
        for j in 1:size(nf, 1)
            
            if nf[j, parse(Int64, data[loadedNodes_pos + vv])] != 0
                
                global loads[nf[j, parse(Int64, data[loadedNodes_pos + vv])], 1] = parse(Float64, data[loadedNodes_pos + vv + j])
                
            end
            
        end
        
    end
    return loads
end

function vonMisesOut(vm, iel, nxe, output)
    if iel == 1
        write(output, "\nvon Mises stress at centroid of elements:\n")
        write(output, "Element\t\tStress\n")
    end
    @printf output "%0.4e " vm
    iel%nxe==0 && @printf output "\n"
end

function megaPrint(vec,vecName;values=false, percentage=5)
    #=
    Debuggin aid. Prints:
        a certain percentage of random elements of a vector
        the vector's standard deviation
        the vector's mean
    =#
    println(vecName)
    values == true && (println(vec[rand(1:length(vec),floor(Int32,length(vec)*percentage/100))]))
    println("mean: $(mean(vec))")
    println("std: $(std(vec))")
end

end