module mtx

    ### This module provides matrix algebra and processing routines

function fkdiag(g_g, neq, nels, process)

    # Computes the skyline profile
    
    write(process, "\nmtx.fkdiag\n\n")
    
    kdiag = zeros(Int32, neq)
    for ll = 1:nels
    
        for i = 1:size(g_g, 1)
            # dof não nulo base do ele
    
            iwp1 = 1
            if g_g[i, ll] != 0
                
                for j = 1:size(g_g, 1)
                    # dof não nulo de comparação do ele
    
                    if g_g[j, ll] != 0
    
                        im = g_g[i, ll] - g_g[j, ll] + 1
                        im > iwp1 && (iwp1 = im) # iwp1 salva maior diferença de DOF + 1 p/ ele
    
                    end
    
                end

                # maior dif. de DOFs não nulos > kdiag[dof não nulo base] então
                # kdiag[dof não nulo base] = maior diferença de DOF + 1
                iwp1 > kdiag[g_g[i, ll]] && (kdiag[g_g[i, ll]] = iwp1)
    
            end
    
        end
    
    end
    
    for i = 2:neq
    
        kdiag[i] += kdiag[i - 1]
          
    end
    
    write(process, "kdiag = $(kdiag)\n")
    
    return kdiag
        
end

function para_fkdiag(g_g, neq, nels, process)

    # Computes the skyline profile
    
    write(process, "\nmtx.fkdiag\n\n")
    
    kdiag = zeros(Int32, neq)
    for ll = 1:nels
    
        for i = 1:size(g_g, 1)
            # dof não nulo base do ele
    
            iwp1 = 1
            if g_g[i, ll] != 0
                
                for j = 1:size(g_g, 1)
                    # dof não nulo de comparação do ele
    
                    if g_g[j, ll] != 0
    
                        im = g_g[i, ll] - g_g[j, ll] + 1
                        im > iwp1 && (iwp1 = im) # iwp1 salva maior diferença de DOF + 1 p/ ele
    
                    end
    
                end

                # maior dif. de DOFs não nulos > kdiag[dof não nulo base] então
                # kdiag[dof não nulo base] = maior diferença de DOF + 1
                iwp1 > kdiag[g_g[i, ll]] && (kdiag[g_g[i, ll]] = iwp1)
    
            end
    
        end
    
    end
    
    for i = 2:neq
    
        kdiag[i] += kdiag[i - 1]
          
    end
    
    write(process, "kdiag = $(kdiag)\n")
    
    return kdiag
        
end

function formnf(nf, process)

    # Forms the node freedom matrix nf
    
    write(process, "\nmtx.formnf\n\n")

    m = 0
    for j = 1:size(nf, 2)

        for i = 1:size(nf, 1)

            if nf[i, j] != 0

                m += 1
                nf[i, j] = m

            end

        end
    end

    write(process, "nf:\n")
    for i in 1:size(nf, 1)

        write(process, "$(nf[i, :])\n")
        
    end

    return nf
    
end

end