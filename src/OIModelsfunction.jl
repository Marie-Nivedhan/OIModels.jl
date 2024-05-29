##Functions for all Models
#Definition of imagingdata(M::Model,tabuv,pixsize,N)
""" 
    This fuction calculates the image matrix of the model 'M' using the fourier space of the data 'tabuv' observed. The result will be 
    sized N*N where each pixel pixsize sized 
    
    Parameters
    ----------
    M : Model
        Predefined model

    tabuv : Array{Quantity{Float64, NoDims, Unitful.FreeUnits{(rad^-1,), NoDims, nothing}}, 3}
            Detected fourier space coordinates (u,v) for 3 wavelengths 

    pixsize : Quantity{Int64, NoDims, Unitful.FreeUnits{(mas,), NoDims, nothing}}
            The size of one pixel 
    
    N: Int64
        Size of the matrix

    Returns
    -------
    result : Matrix{Float64}
    Matrix image
"""
function imagingdata(M::Model,tabuv,pixsize,N)

    
    u = tabuv[1,:,:] .|>u"mas^-1"
    v = tabuv[2,:,:] .|>u"mas^-1"
    uv = vcat(reshape(v,1,:),reshape(u,1,:)) .* pixsize .|> NoUnits
    nfftplan  = plan_nfft(uv,(N,N))
    img = nfftplan'*interferometry_fourier(M,u[:],v[:])./(N*N)
    img1 = real.(img)
    return img1
end

#Definition of imagingxy(M::Model,pixsize,N)
""" 
    This fuction calculates the image matrix of the model 'M' in a real space. The result will be 
    sized N*N where each pixel pixsize sized 
    
    Parameters
    ----------
    M : Model
        Predefined model

    pixsize : Quantity{Int64, NoDims, Unitful.FreeUnits{(mas,), NoDims, nothing}}
            The size of one pixel 

    N: Int64
        Size of the matrix

    Returns
    -------
    result : Matrix{Float64}
    Matrix image
"""
function imagingxy(M::Model,pixsize,N)

    iy = ((1:N).-N/2 .-1).*pixsize
    ix = ((1:N).-N/2 .-1).*pixsize
    y = iy * ones(N)'
    x = ones(N)*ix'
    img = interferometry_image(M,x,y; atol=pixsize/2)
    img1 = real.(img)
    return img1

end


#Definition of imaginguv(M::Model,pixsize,N)
""" 
    This fuction calculates the image matrix of the model 'M' in a fourier space. The result will be 
    sized N*N where each pixel pixsize sized 
    
    Parameters
    ----------
    M : Model
        Predefined model

    pixsize : Quantity{Int64, NoDims, Unitful.FreeUnits{(mas,), NoDims, nothing}}
            The size of one pixel 

    N: Int64
        Size of the matrix
        
    Returns
    -------
    result : Matrix{Float64}
    Matrix image
"""
function imaginguv(M::Model,pixsize,N)

    fov = pixsize*N
    iv =  ((1:N).-N/2 .-1)./fov
    iu =  ((1:N).-N/2 .-1)./fov
    u2 = iu * ones(N)'
    v2 = ones(N)*iv'
    uv2 = vcat(reshape(v2,1,:),reshape(u2,1,:)) .* pixsize .|> NoUnits
    nfftplan2  = plan_nfft(uv2,(N,N))
    img = nfftplan2'*interferometry_fourier(M,u2[:],v2[:])./(N*N)
    img1 = real.(img)
    return img1

end


function findobject(M::Star,donnees,intervalle,uv)
    v, re = destructure(M)
    for i in p1:p2
        v[1]=i*u"mas"
        for j in p1:p2    
            v[2]=j*u"mas"
            if length(v)==3
                for k in p1:p2
                    v[3]=abs(k*u"mas")
                    M=re(v)
                    if interferometry_fourier(M,uv)==donnees
                        return M
                    end
                end
            else
                M=re(v)
                if interferometry_fourier(M,uv)==donnees
                    return M
                end
            end
        end
    end
    return "L'objet ne se trouve pas dans la plage de position implémentée"
end
    
function finddisk(M::Disk,donnees,intervalle,uv)
    v, re = destructure(M)
    for i in p1:p2
        v[1]=i*u"mas"
        for j in p1:p2    
            v[2]=j*u"mas"
            M=re(v)
            if interferometry_fourier(M,uv)==donnees
                return M
            end
        end
    end
    return "L'objet ne se trouve pas dans la plage de position implémentée"
end
    