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


function findimg(M::Star,data,range,uv)
    v, re = destructure(M)
    I=length(range)
    img = zeros(Float64, I, I)
    T=Vector{Float64}()
        for x in range
            v[1]=x*u"mas"
            t=0
            for y in range   
                v[2]=y*u"mas"
                M=re(v)
                t=sum(abs2,interferometry_fourier(M,uv).-data)
                push!(T,t)
            end    
        end 
    for idx in eachindex(img)
        img[idx]=T[idx]
    end
    return img
end

function findobj(img)
    min=minimum(img)
    vect=[]
    x=length(img[:,1])
    y=length(img[1,:])
    for i in 1:x
        for j in 1:y
            if img[i,j]==min
                vect=[(j-y/2)*u"mas",(i-x/2)*u"mas"]
                return vect
            end
        end
    end
end

#= 
function findimg(M::Disk,donnees,intervalle,uv)
    v, re = destructure(M)
    I=length(intervalle)
    img = zeros(Float64, I, I)
    T=Vector{Float64}()
    for x in intervalle
        v[1]=x*u"mas"
        t=0s
        for y in intervalle   
            v[2]=y*u"mas"
                for R in 1:I/2
                v[3]=R*u"mas"
                M=re(v)
                t=sum(abs2,interferometry_fourier(M,uv).-donnees)
                push!(T,t)
            end
        end
    end
    for idx in eachindex(img)
        img[idx]=T[idx]
    end
    return img
end
 =#

function stripeunits(x::Quantity)
    ustrip(upreferred(x))
end

function stripeunits(x::AbstractArray)
    stripeunits.(x)
end
stripeunits(x) =x

stripeunits(x::Model)= fmap(stripeunits,x)

function findmodelmcmc(M::Star,tabuv,donnees,iteration)
    # Prepare the data for Turing
    d = hcat(real.(donnees), imag.(donnees))[:]

    # Calculate the standard deviation of the observed data
    obs_std = std(d)

    # Définir le modèle Turing

    @model function gdemo(obs,uv,obs_std)

        x ~ Uniform(-127.0, 128.0)  # Likelihood
        y ~ Uniform(-127.0, 128.0) #Normal(s, sigma) 
        # Simulez les données en utilisant la fonction interferometry_fourier
        model_data = interferometry_fourier(Star(x*1u"mas", y*1u"mas"), uv)
        model_data =hcat(real.(model_data),imag.(model_data))[:]
        
        # Vraisemblance 
        obs ~ MvNormal(model_data,obs_std^2 * I)

    end
    
    chain = sample(gdemo(d,tabuv,obs_std), MH(), iteration)#iteration=100000 c'est bon
    # Extraire les échantillons de x et y
    x_samples = chain[:x]
    y_samples = chain[:y]

    x_values = sample(x_samples)
    y_values = sample(y_samples)
    # Trouver les valeurs les plus probables (pics) de x et y
    x_peak = mode(x_values)
    y_peak = mode(y_values)

    M=Star(x_peak*u"mas",y_peak*u"mas")
    return M
end

function findmodelmcmc(M::Disk,tabuv,donnees,iteration)#iteration =2000000 ça fonctionne mais ~15 min
    # Prepare the data for Turing
    d = hcat(real.(donnees), imag.(donnees))[:]

    # Calculate the standard deviation of the observed data
    obs_std = std(d)

    # Définir le modèle Turing

    @model function gdemo(obs,uv,obs_std)

        x ~ Uniform(-127.0, 128.0)  # Likelihood
        y ~ Uniform(-127.0, 128.0) #Normal(s, sigma) 
        R ~ Uniform(0.0, 128.0)
        # Simulez les données en utilisant la fonction interferometry_fourier
        model_data = interferometry_fourier(Disk(x*1u"mas", y*1u"mas",R*1u"mas"), uv)
        model_data =hcat(real.(model_data),imag.(model_data))[:]
        
        # Vraisemblance 
        obs ~ MvNormal(model_data,obs_std^2 * I)

    end
    
    chain = sample(gdemo(d,tabuv,obs_std),MH(), iteration)
    
    # Extraire les échantillons de x et y
    x_samples = chain[:x]
    y_samples = chain[:y]
    R_samples = chain[:R]

    x_values = sample(x_samples)
    y_values = sample(y_samples)
    R_values = sample(R_samples)
    # Trouver les valeurs les plus probables (pics) de x et y
    x_peak = mode(x_values)
    y_peak = mode(y_values)
    R_peak = mode(R_values)

    M=Disk(x_peak*u"mas",y_peak*u"mas",R_peak*u"mas")
    return M
end

function findmodelmcmc(M::Gauss,tabuv,donnees,iteration)#iteration =2000000 ça fonctionne mais ~10 min
    # Prepare the data for Turing
    d = hcat(real.(donnees), imag.(donnees))[:]

    # Calculate the standard deviation of the observed data
    obs_std = std(d)

    # Définir le modèle Turing

    @model function gdemo(obs,uv,obs_std)

        x ~ Uniform(-127.0, 128.0)  # Likelihood
        y ~ Uniform(-127.0, 128.0) #Normal(s, sigma) 
        R ~ Uniform(0.0, 128.0)
        # Simulez les données en utilisant la fonction interferometry_fourier
        model_data = interferometry_fourier(Gauss(x*1u"mas", y*1u"mas",R*1u"mas"), uv)
        model_data =hcat(real.(model_data),imag.(model_data))[:]
        
        # Vraisemblance 
        obs ~ MvNormal(model_data,obs_std^2 * I)

    end
    
    chain = sample(gdemo(d,tabuv,obs_std),MH(), iteration)
    
    # Extraire les échantillons de x et y
    x_samples = chain[:x]
    y_samples = chain[:y]
    R_samples = chain[:R]

    x_values = sample(x_samples)
    y_values = sample(y_samples)
    R_values = sample(R_samples)
    # Trouver les valeurs les plus probables (pics) de x et y
    x_peak = mode(x_values)
    y_peak = mode(y_values)
    R_peak = mode(R_values)

    M=Gauss(x_peak*u"mas",y_peak*u"mas",R_peak*u"mas")
    return M
end

function findmodelmcmc(M::Ring,tabuv,donnees,iteration)#iteration =4000000 ça fonctionne mais ~40 min
    # Prepare the data for Turing
    d = hcat(real.(donnees), imag.(donnees))[:]

    # Calculate the standard deviation of the observed data
    obs_std = std(d)

    # Définir le modèle Turing

    @model function gdemo(obs,uv,obs_std)

        x ~ Uniform(-127.0, 128.0)  # Likelihood
        y ~ Uniform(-127.0, 128.0) #Normal(s, sigma) 
        R1 ~ Uniform(0.0, 128.0)
        R2 ~ Uniform(0.0, 128.0)
        # Simulez les données en utilisant la fonction interferometry_fourier
        model_data = interferometry_fourier(Ring(x*1u"mas", y*1u"mas",R1*1u"mas",R2*1u"mas"), uv)
        model_data =hcat(real.(model_data),imag.(model_data))[:]
        
        # Vraisemblance 
        obs ~ MvNormal(model_data,obs_std^2 * I)

    end
    
    chain = sample(gdemo(d,tabuv,obs_std),MH(), iteration)
    
    # Extraire les échantillons de x et y
    x_samples = chain[:x]
    y_samples = chain[:y]
    R1_samples = chain[:R1]
    R2_samples = chain[:R2]

    x_values = sample(x_samples)
    y_values = sample(y_samples)
    R1_values = sample(R1_samples)
    R2_values = sample(R2_samples)

    # Trouver les valeurs les plus probables (pics) de x et y
    x_peak = mode(x_values)
    y_peak = mode(y_values)
    R1_peak = mode(R1_values)
    R2_peak = mode(R2_values)

    M=Ring(x_peak*u"mas",y_peak*u"mas",R1_peak*u"mas",R2_peak*u"mas")
    return M
end