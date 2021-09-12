module InOut

using Printf



### This module provides IO routines

function Extract(file, process)

    # Returns single array containing strings with the data in input file

    write(process, "InOut.Extract\n")

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

    return data
    
end

function DispOut(output, type_2d, disp, nf, nn, nodof, process)

    # Outputs displacements

    write(process, "\nInOut.DispOut\n")

    if nodof == 2

        if type_2d == "axisymmetric"

            @printf output "Node\tr-disp\t\t\tz-disp\n\n"
        
        else
            
            @printf output "Node\tx-disp\t\t\ty-disp\n\n"
        
        end

    elseif nodof == 3

        @printf output "Node\tx-disp\t\t\ty-disp\t\t\tz-disp\n\n"

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

function StressOut(iel, gc, sigma, i, nip, process, output)

    # Outputs Stress values

    write(process, "\nInOut.StressOut\n")
    @printf output "\t\t%i\t\t\t\t%i\t\t%i" iel i gc[1]

    for mm in 1:size(sigma, 1)

        @printf output "\t\t%.4e" sigma[mm]
        mm == size(sigma, 1) ? (@printf output "\n") : false
        
    end

    i == nip && (@printf output "\n")
    
end

function HeaderStress(type_2d, nodof, process, output)

    # Outputs header of stress table

    write(process, "\nInOut.HeaderStress\n")    

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


end