module mtx
import InOut
    ### This module provides matrix algebra and processing routines

function fkdiag(g_g, neq, nels)

    # Computes the skyline profile
    
    kdiag = zeros(Int64, neq)
    for ll = 1:nels
    
        for i = 1:size(g_g, 1)
    
            iwp1 = 1
            if g_g[i, ll] != 0
                
                for j = 1:size(g_g, 1)
    
                    if g_g[j, ll] != 0
    
                        im = g_g[i, ll] - g_g[j, ll] + 1
                        im > iwp1 && (iwp1 = im)
    
                    end
    
                end
                (iwp1 > kdiag[g_g[i, ll]]) && (kdiag[g_g[i, ll]] = iwp1)
            end
    
        end
    
    end
    
    for i = 2:neq
    
        kdiag[i] += kdiag[i - 1]
          
    end
    
    return kdiag
        
end

function para_fkdiag(g_g, neq, nels)

    # Computes the skyline profile
    
    kdiag = zeros(Int64, neq)
    for ll = 1:nels
    
        for i = 1:size(g_g, 1)
    
            iwp1 = 1
            if g_g[i, ll] != 0
                
                for j = 1:size(g_g, 1)
    
                    if g_g[j, ll] != 0
    
                        im = g_g[i, ll] - g_g[j, ll] + 1
                        im > iwp1 && (iwp1 = im) 
    
                    end
    
                end

                iwp1 > kdiag[g_g[i, ll]] && (kdiag[g_g[i, ll]] = iwp1)
    
            end
    
        end
    
    end
    
    for i = 2:neq
    
        kdiag[i] += kdiag[i - 1]
          
    end
    
    return kdiag
        
end

function formnf(nf)

    # Forms the node freedom matrix nf
    
    m = 0
    for j = 1:size(nf, 2)

        for i = 1:size(nf, 1)

            if nf[i, j] != 0

                m += 1
                nf[i, j] = m

            end

        end
    end

    return nf
    
end

function bandwidth(g)
    # Returns the bandwidth of element with DOFs g
    return maximum(g)-minimum(g)
end

function formlump!(diag, mm, g)
    # Assembles the lumped global mass matrix as a vector
    for i = 1:length(g)
        g[i] != 0 && (diag[g[i]] += mm)
    end
    return diag
end

function EigenSkyline!(kdiag, ku)
    # Convert to skyline form
    kh = Array{Float64}(undef, kdiag[end])
    kh[1]=ku[1,1]
    for i=2:length(kdiag)
        idiag=kdiag[i]-kdiag[i-1]
        for j = 1:idiag
            (j == 1) ? k = 1 : (k+=1)
            kh[k] = ku[i+j-idiag,1-j+idiag]
        end
    end
    return kh
end

function bandred!(ku, neq)
    # This subroutine transforms a real symmetric band matrix ku,
    # of order n and band width iw, to tridiagonal form by an appropriate
    # sequence of Jacobi rotations. During the transformation, the
    # property of the band matrix is maintained. The method yields
    # a tridiagonal matrix, the diagonal elements of which are in
    # diag[n] and off-diagonal elements in udiag[n].
    udiag = zeros(Float64, neq)
    if (size(ku,1)-2)>=1
        for k = 1:(size(ku,1)-2)
            maxr=size(ku,2)-1
            (size(ku,1)-k<(size(ku,2)-1)) && (maxr=size(ku,1)-k)
            for irr = 2:maxr
                ir=2+maxr-irr
                kr=k+ir
                for j = kr:size(ku,1):(size(ku,2)-1)
                    if j!=kr
                        abs(g)<1e-40 && break
                        jm=j-(size(ku,2)-1)
                        b=-ku[jm-1,size(ku,2)]/g
                        iugl=j-(size(ku,2)-1)
                    else
                        ku[k,ir+1]<1e-40 && break
                        b=-ku[k,ir]/ku[k,ir+1]
                        iugl=k
                    end
                    s=1/sqrt(1+b*b)
                    c=b*s
                    c2=c^2
                    s2=s^2
                    cs=c*s
                    u=c2*ku[j-1,1]-2*cs*ku[j-1,2]+s2*ku[j,1]
                    u1=s2*ku[j-1,1]+2*cs*ku[j-1,2]+c2*ku[j,1]
                    ku[j-1,2]=cs*(ku[j-1,1]-ku[j,1])+(c2-s2)*ku[j-1,2]
                    ku[j-1,1]=u
                    ku[j,1]=u1
                    j2=j-2
                    for l=iugl:j2
                        jl=j-l
                        u=c*ku[l,jl]-s*ku[l,jl+1]
                        ku[l,jl+1]=s*ku[l,jl]+c*ku[l,jl+1]
                        ku[l,jl]=u
                    end
                    jm=j-(size(ku,2)-1)
                    j!=kr && (ku[jm-1,size(ku,2)]=c*ku[jm-1,size(ku,2)]-s*g)
                    maxl=size(ku,2)-2
                    (size(ku,1)-j<(size(ku,2)-2)) && (maxl=size(ku,1)-j)
                    if maxl>0
                        for l=1:maxl
                            u=c*ku[j-1,l+2]-s*ku[j,l+1]
                            ku[j,l+1]=s*ku[j-1,l+2]+c*ku[j,l+1]
                            ku[j-1,l+2]=u
                        end
                    end 
                    if (j+size(ku,2)-1)<=size(ku,1)
                        g=-s*ku[j,end]
                        ku[j,end]*=c
                    end
                end
            end
        end
    end 
    udiag[1]=0
    if 2<=size(ku,1)
        for i=2:size(ku,1)
            udiag[i]=ku[i-1,2]
        end
    end 
    return ku[:,1], udiag
end

end