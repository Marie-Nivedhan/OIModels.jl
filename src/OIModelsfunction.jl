##Functions for all Models
#Definition of imagingdata(M::Model,tabuv,pixsize,N)
""" 
    This fuction calculates the image matrix of th dirty image of the model 'M' using the fourier space of the data 'tabuv' observed. The result will be 
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

    u = tabuv[1,:,:] .|>u"mas^-1"       #= Extracts and converts the u-coordinates of the fourier space datas tabuv 
                                           to milliarcseconds inverse=#
    v = tabuv[2,:,:] .|>u"mas^-1"       #= Extracts and converts the v-coordinates of the fourier space datas tabuv 
                                           to milliarcseconds inverse=#
    uv = vcat(reshape(v,1,:),reshape(u,1,:)) .* pixsize .|> NoUnits     #= Concatenates and reshapes the u and v coordinates,
                                                                           scales by the pixel size, and removes units for the next step. =#
    nfftplan  = plan_nfft(uv,(N,N))     #= Creates a square N size Non-uniform Fast Fourier Transform (NFFT) plan for the given coordinates matrix uv.  =#
    img = nfftplan'*interferometry_fourier(M,u[:],v[:])./(N*N)      #=Applies the inverse NFFT to the fourier transform data obtained from 
                                                                      the interferometry_fourier function for the model M, normalizing by the 
                                                                      total number of pixels.=#
    return real.(img)       #=The function returns the real part of the reconstructed image=#
    
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

    iy = ((1:N).-N/2 .-1).*pixsize      #= ix represent the x coordinates in the image plane, scaled by the pixel size and centered 
                                           around zero. =#
    ix = ((1:N).-N/2 .-1).*pixsize      #= iy represent the y coordinates in the image plane, scaled by the pixel size and centered
                                           around zero. =#
    y = iy * ones(N)'   #=Matrix representing the y coordinates for each pixel in the image.=#
    x = ones(N)*ix'     #=Matrix representing the x coordinates for each pixel in the image.=#
    img = interferometry_image(M,x,y; atol=pixsize/2)       #= Calls the interferometry_image function to reconstruct the image from the model M using the x and y coordinates.
                                                               The atol parameter is set to half the pixel size for accuracy.=#
    return real.(img)       #= The function returns the real part of the reconstructed image. =#

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

    fov = pixsize*N     #= Calculates the field of view by multiplying the pixel size by the number of pixels. =#
    iv =  ((1:N).-N/2 .-1)./fov     #= iv represent the v coordinates in the image plane, normalized by the field of view and 
                                       centered around zero. =#
    iu =  ((1:N).-N/2 .-1)./fov     #= iu represent the u coordinates in the image plane, normalized by the field of view and 
                                       centered around zero. =#
    u = iu * ones(N)'       #=Matrix representing the u coordinates =#
    v = ones(N)*iv'     #=Matrix representing the v coordinates =#
    uv = vcat(reshape(v,1,:),reshape(u,1,:)) .* pixsize .|> NoUnits     #= Concatenates and reshapes the u and v arrays, scales 
                                                                           by the pixel size, and removes units. =#
    nfftplan  = plan_nfft(uv,(N,N))     #= Creates a square N size Non-uniform Fast Fourier Transform (NFFT) plan for the given 
                                           coordinates matrix uv.  =#
    img = nfftplan'*interferometry_fourier(M,u[:],v[:])./(N*N)      #= Applies the inverse NFFT to the fourier transform data obtained from 
                                                                       the interferometry_fourier function for the model M, normalizing by 
                                                                       the total number of pixels.=#
    return real.(img)       #=The function returns the real part of the reconstructed image=#

end

#Definition of findimg(M::Star,data,range,uv)
""" 
    This fuction calculates the image matrix of the Star 'M' corresponding to the datas of 'data' in the fourier space 'uv'. The result 
    will be sized N*N where each pixel pixsize sized 
    
    Parameters
    ----------
    M : Model
        Predefined model

    data : Quantity{Int64, NoDims, Unitful.FreeUnits{(mas,), NoDims, nothing}}
           fourier space datas

    range: Int64
            range scope for finding the components of 'M' 
        
   Returns
   -------
   result : Matrix{Float64}
            Matrix image    
"""
function findimg(M::Star,data,range,uv)
    v, re = destructure(M)      #= Extracts components v and re from the M model. =#
    I=length(range)     #= Determines the length of the range, which will determine the size of the resulting image img. =#
    img = zeros(Float64, I, I)      #= Initializes an image matrix filled with zeros of size I x I. =#
    T=Vector{Float64}()     #=  Initializes an empty vector T to store intermediate results. =#
        for x in range      #=  Iterates over each element x in the range. =#
            v[1]=x*u"mas"       #= Sets the first component of v to x scaled by milliarcseconds =#
            t=0     #=  Initializes t to zero. =#
            for y in range      #= Iterates over each element y in the range. =#
                v[2]=y*u"mas"       #= Sets the second component of v to y scaled by milliarcseconds =#
                M=re(v)     #= Updates M using the function re with v as an argument. =#
                t=sum(abs2,interferometry_fourier(M,uv).-data)      #= Computes the sum of squared differences between the interferometry
                                                                       Fourier transform of M with uv coordinates and the given data. This
                                                                       quantifies the difference or error between the model and the actual 
                                                                       data. =#
                push!(T,t)      #= Appends the computed error t to the vector T =#
            end    
        end 
    for idx in eachindex(img)       #= Iterates over each index in the img matrix. =#
        img[idx]=T[idx]     #= Assigns the corresponding value from T to the img matrix. =#
    end
    return img      #= Returns the img matrix, which represents the optimized image =#
end

#Definition of findobj(img)
""" 
    This fuction calculates the vector position the Star 'M' in the result of findimg.
    
    Parameters
    ----------
    img : Matrix{Float64}
          Matrix image 

   Returns
   -------
   result :  
"""
function findobj(img)

    min=minimum(img)        #= Finds the minimum value in the entire img matrix. =#
    vect=[]     #= Initializes an empty vector vect to store the position of the minimum value. =#
    x=length(img[:,1])      #=  Gets the number of rows (x) in the image matrix. =#
    y=length(img[1,:])      #=  Gets the number of columns (y) in the image matrix. =#
    for i in 1:x        #= Iterates over rows. =#
        for j in 1:y        #= Iterates over columns. =#
            if img[i,j]==min        #= Checks if the current element img[i, j] equals min. =#
                vect=[(j-y/2)*u"mas",(i-x/2)*u"mas"]      #= Calculates the position (j - y/2) and (i - x/2) scaled by milliarcseconds  =#  
                return vect     #= Returns the vector vect containing the position of the minimum value in milliarcseconds relative to the 
                                   center of the image. =#
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

#Definition of stripeunits(x::Quantity)
""" 
    This fuction converts the Quantity 'x' to the international scientific units.
    
    Parameters
    ----------
    x : Quantity

   Returns
   -------
   result :  typeof(x)
"""
function stripeunits(x::Quantity)
    ustrip(upreferred(x))
end
function stripeunits(x::AbstractArray)
    stripeunits.(x)
end
stripeunits(x) =x
stripeunits(x::Model)= fmap(stripeunits,x)



function findgraph(M::Star,data,range,uv) # ~1 min
    f(x) = sum(abs2.(interferometry_fourier(Star([x[1],x[2]]), uv) .- data))
    # Génération de données pour le tracé
    x = y = range.*u"mas"
    z = [f([a, b]) for a in x, b in y]


    # Trouver les indices du minimum de z
    min_index = argmin(z)
    min_x, min_y = x[min_index[1]], y[min_index[2]]
    min_z = minimum(z)
    println(min_z)
    M=Star([min_x, min_y])
    return M
end

function findgraph(M::Disk,data,range,uv) # ~1 min
    f(x) = sum(abs2.(interferometry_fourier(Disk([x[1]],[x[2]]), uv) .- data))
    # Génération de données pour le tracé
    x = y = range.*u"mas"
    R=(0:step(range):last(range)).*u"mas"
    z = [f([[a, b],[c]]) for a in x, b in y,c in R]


    # Trouver les indices du minimum de z
    min_index = argmin(z)
    min_x, min_y, min_R= x[min_index[1]], y[min_index[2]],R[min_index[3]]
    min_z = minimum(z)
    println(min_z)
    M=Star([min_x, min_y],[min_R])
    return M
end

#Definition of findmodelmcmc(M::Model,tabuv,donnees,iteration)
""" 
    This fuction function uses Markov Chain Monte Carlo (MCMC) sampling to fit a model to interferometric data. 
    The function defines a probabilistic model based on the type of astronomical object (e.g., star, disk, Gaussian, ring),
    samples from the posterior distribution, and then finds the most probable parameters for the object.
    
    Parameters
    ----------
    M : Model
        Predefined model

    tabuv : Array{Quantity{Float64, NoDims, Unitful.FreeUnits{(rad^-1,), NoDims, nothing}}, 3}
            Detected fourier space coordinates (u,v) for 3 wavelengths 

    donnees : 
    
    iteration: Int64
                number of iteration made during the sampling

    Returns
    -------
    result : Model
             object that we try to find in 'donnees'
"""
function findmodelmcmc(M::Model,tabuv,donnees,iteration)
    
    d = hcat(real.(donnees), imag.(donnees))[:]     #= Concatenates the real and imaginary parts of the data. =#

    obs_std = std(d)        #= Calculates the standard deviation of the observed data. =#       

   
    if M isa Star       #= Checks if M is a Star. =#
        # Defines the Turing model
        @model function gdemostar(obs,uv,obs_std)

            # Prior distributions for the parameters   
            #= Defines uniform priors over the ranges [-127.0, 128.0] for both x and y =#     
            x ~ Uniform(-127.0, 128.0) 
            y ~ Uniform(-127.0, 128.0) 

            # Simulate the model data using interferometry_fourier      
            #= Inside the model, simulates the interferometric data for a Star model using interferometry_fourier =#
            model_data = interferometry_fourier(Star(x*1u"mas", y*1u"mas"), uv)     
            model_data =hcat(real.(model_data),imag.(model_data))[:]
            
            # Likelihood: observed data is normally distributed around model data
            obs ~ MvNormal(model_data,obs_std^2 * I)
    
        end

        # Perform MCMC sampling using Metropolis-Hasting
        chain = sample(gdemostar(d,tabuv,obs_std), MH(), iteration)     #= advised iteration=100000 =#

        # Extract samples
        x_samples = chain[:x]
        y_samples = chain[:y]
    
        x_values = sample(x_samples)
        y_values = sample(y_samples)

        # Find the most probable values
        x_peak = mode(x_values)
        y_peak = mode(y_values)
    
        M=Star(x_peak*u"mas",y_peak*u"mas")
        #Returns the estimated Star object M using the most probable x and y values scaled to milliarcseconds
        return M        

    elseif  M isa Disk      #= Checks if M is a Disk. =#
        # Defines the Turing model
        @model function gdemodisk(obs,uv,obs_std)
            # Prior distributions for the parameters        
            #= Defines uniform priors over the ranges [-127.0, 128.0] for both x and y and over the ranges [0.0, 128.0] for R =#
            x ~ Uniform(-127.0, 128.0)  
            y ~ Uniform(-127.0, 128.0) 
            R ~ Uniform(0.0, 128.0)

            # Simulate the model data using interferometry_fourier  
            #= Inside the model, simulates the interferometric data for a Disk model using interferometry_fourier =#
            model_data = interferometry_fourier(Disk(x*1u"mas", y*1u"mas",R*1u"mas"), uv)       
            model_data =hcat(real.(model_data),imag.(model_data))[:]
            
            # Likelihood: observed data is normally distributed around model data
            obs ~ MvNormal(model_data,obs_std^2 * I)
    
        end
        # Perform MCMC sampling using Metropolis-Hasting
        chain = sample(gdemodisk(d,tabuv,obs_std),MH(), iteration)      #= advised iteration=1000000 =#
        
        # Extract samples
        x_samples = chain[:x]
        y_samples = chain[:y]
        R_samples = chain[:R]
    
        x_values = sample(x_samples)
        y_values = sample(y_samples)
        R_values = sample(R_samples)
        
        # Find the most probable values
        x_peak = mode(x_values)
        y_peak = mode(y_values)
        R_peak = mode(R_values)
    
        M=Disk(x_peak*u"mas",y_peak*u"mas",R_peak*u"mas")
        #Returns the estimated Disk object M using the most probable x, y and R values scaled to milliarcseconds
        return M

    elseif M isa Gauss      #= Checks if M is a Gauss. =#
        # Defines the Turing model
        @model function gdemogauss(obs,uv,obs_std)
            # Prior distributions for the parameters        
            #= Defines uniform priors over the ranges [-127.0, 128.0] for both x and y and over the ranges [0.0, 128.0] for R =#
            x ~ Uniform(-127.0, 128.0)  
            y ~ Uniform(-127.0, 128.0)
            R ~ Uniform(0.0, 128.0)

            # Simulate the model data using interferometry_fourier  
            #= Inside the model, simulates the interferometric data for a Gauss model using interferometry_fourier =#
            model_data = interferometry_fourier(Gauss(x*1u"mas", y*1u"mas",R*1u"mas"), uv)
            model_data =hcat(real.(model_data),imag.(model_data))[:]
            
            # Likelihood: observed data is normally distributed around model data
            obs ~ MvNormal(model_data,obs_std^2 * I)
    
        end

        # Perform MCMC sampling using Metropolis-Hasting
        chain = sample(gdemogauss(d,tabuv,obs_std),MH(), iteration)     #=advised iteration =2000000=# 
        
        # Extract samples
        x_samples = chain[:x]
        y_samples = chain[:y]
        R_samples = chain[:R]
    
        x_values = sample(x_samples)
        y_values = sample(y_samples)
        R_values = sample(R_samples)
        
        # Find the most probable values
        x_peak = mode(x_values)
        y_peak = mode(y_values)
        R_peak = mode(R_values)
    
        M=Gauss(x_peak*u"mas",y_peak*u"mas",R_peak*u"mas")
        #Returns the estimated Gauss object M using the most probable x, y and R values scaled to milliarcseconds
        return M

    elseif M isa Ring       #= Checks if M is a Ring. =#
        # Defines the Turing model
        @model function gdemoring(obs,uv,obs_std)
            # Prior distributions for the parameters        
            #= Defines uniform priors over the ranges [-127.0, 128.0] for both x and y and over the ranges [0.0, 128.0] for both R1 and R2 =#
            x ~ Uniform(-127.0, 128.0)  
            y ~ Uniform(-127.0, 128.0) 
            R1 ~ Uniform(0.0, 128.0)
            R2 ~ Uniform(0.0, 128.0)

            # Simulate the model data using interferometry_fourier  
            #= Inside the model, simulates the interferometric data for a Ring model using interferometry_fourier =#
            model_data = interferometry_fourier(Ring(x*1u"mas", y*1u"mas",R1*1u"mas",R2*1u"mas"), uv)
            model_data =hcat(real.(model_data),imag.(model_data))[:]
            
            # Likelihood: observed data is normally distributed around model data
            obs ~ MvNormal(model_data,obs_std^2 * I)
    
        end

        # Perform MCMC sampling using Metropolis-Hasting
        chain = sample(gdemoring(d,tabuv,obs_std),MH(), iteration)      #=advised iteration =5000000=# 
        
        # Extract samples
        x_samples = chain[:x]
        y_samples = chain[:y]
        R1_samples = chain[:R1]
        R2_samples = chain[:R2]
    
        x_values = sample(x_samples)
        y_values = sample(y_samples)
        R1_values = sample(R1_samples)
        R2_values = sample(R2_samples)
    
        # Find the most probable values
        x_peak = mode(x_values)
        y_peak = mode(y_values)
        R1_peak = mode(R1_values)
        R2_peak = mode(R2_values)

        M=Ring(x_peak*u"mas",y_peak*u"mas",R1_peak*u"mas",R2_peak*u"mas")
        #Returns the estimated Ring object M using the most probable x, y, R1 and R2 values scaled to milliarcseconds
        return M
    end
end


#Definition of functions findmodeloptim
""" 
    These fuctions use limited memory Broyden–Fletcher–Goldfarb–Shanno (LBFGS)algorithm  optimizer to fit a model to interferometric data. 
    These functions define an objective function 'f' that will be optimized. Then using loops it will find the most probable parameters for 
    the object.
    
    Parameters
    ----------
    M : Model
        Predefined model

    tabuv : Array{Quantity{Float64, NoDims, Unitful.FreeUnits{(rad^-1,), NoDims, nothing}}, 3}
            Detected fourier space coordinates (u,v) for 3 wavelengths 

    data : 
    
    
    range: 
            

    Returns
    -------
    result : Model
             object that we try to find in 'data'
"""
function findmodeloptim(M::Star,range,data,tabuv) # advised range =-127.0:8.0:128.0
    # Objective function f(b) to minimize
    #= Defines f(b) to compute the sum of squared differences between observed data and model predictions (interferometry_fourier) 
    for a Star model with parameters b. =#
    function f(b, donnees,uv)
        x = b[1] * u"mas"
        y = b[2] * u"mas"
        sum(abs2, interferometry_fourier(Star([x, y]), uv) .- donnees)
    end

    # Loop over each optimization method
    #= Iterates over the range of initial values for x and y. For each combination, it initializes x0 and performs optimization
       using optimize with the LBFGS() method =#
    for i in range
        for j in range
                x0 = [i,j]
                # Optimization using LBFGS method
                opt_result = optimize(b -> f(b, data,tabuv), x0,LBFGS(); autodiff = :forward)

                # Retrieve results
                global_min = opt_result.minimum
                minimizer = opt_result.minimizer

                # Check if minimum value is zero (indicating a good fit)
                if isapprox(global_min,0.0; atol=5)
                    return Star(minimizer.*u"mas")
                end
        end
    end    
    error("Optimization did not find a satisfactory solution.")
end

function findmodeloptim(M::Gauss,range,data,tabuv) # advised range =-127.0:64.0:128.0
    # Objective function f(b) to minimize
    #= Defines f(b) to compute the sum of squared differences between observed data and model predictions (interferometry_fourier) 
    for a Gauss model with parameters b. =#
    function f(b, donnees,uv)
        x = b[1]u"mas"
        y = b[2]u"mas"
        R = b[3]u"mas"
        sum(abs2, interferometry_fourier(Gauss([x, y],[R]), uv) .- donnees)
    end

    # Loop over each optimization method
    #= Iterates over the range of initial values for x and y. For each combination, it initializes x0 and performs optimization
       using optimize with the LBFGS() method =#
    for i in range
        for j in range
            for k in 0.0:step(range)*2:last(range)
                x0 = [i,j,k]
                # Optimization using LBFGS method
                opt_result = optimize(b -> f(b, data,tabuv), x0,LBFGS(); autodiff = :forward)

                # Retrieve results
                global_min = opt_result.minimum
                minimizer = opt_result.minimizer

                # Uses isapprox(global_min, 0.0; atol=0.5) to determine if the optimization successfully minimized the objective function
                if isapprox(global_min,0.0; atol=0.5)
                    return Gauss(minimizer[1].*u"mas",minimizer[2].*u"mas",minimizer[3].*u"mas")
                end
            end
        end
    end   
    error("Optimization did not find a satisfactory solution.") 
end

function findmodeloptim(M::Ring,range,data,tabuv) # advised range =-127.0:8.0:128.0
    # Objective function f(b) to minimize
    #= Defines f(b) to compute the sum of squared differences between observed data and model predictions (interferometry_fourier) 
    for a Ring model with parameters b. =#
    function f(b, donnees,uv)
        x = b[1]u"mas"
        y = b[2]u"mas"
        R1 = b[3]u"mas"
        R2 = b[4]u"mas"
        sum(abs2, interferometry_fourier(Ring(x, y, R1, R2), uv) .- donnees)
    end        

    # Loop over each optimization method
    #= Iterates over the range of initial values for x and y. For each combination, it initializes x0 and performs optimization
       using optimize with the LBFGS() method =#
    for i in -127.0:16.0:128.0
        for j in -127.0:16.0:128.0
            for k in 0.0:32.0:128.0
                for l in 0.0:32.0:128.0
                    x0 = [i,j,k,l]
                    # Optimisation avec la méthode actuelle
                    opt_result = optimize(b -> f(b, data,tabuv), x0,LBFGS())

                    # Retrieve results
                    global_min = opt_result.minimum
                    minimizer = opt_result.minimizer

                    # Uses isapprox(global_min, 0.0; atol=0.5) to determine if the optimization successfully minimized the objective function
                    if isapprox(global_min,0.0; atol=0.5)
                        return Ring(minimizer[1]*u"mas",minimizer[2]*u"mas",minimizer[3]*u"mas",minimizer[4]*u"mas")
                    end
                end
            end
        end
    end    
    error("Optimization did not find a satisfactory solution.")
end

function findmodeloptim(M::Disk,range,data,tabuv) # advised range =-127.0:8.0:128.0
    
    # Objective function f(b) to minimize
    #= Defines f(b) to compute the sum of squared differences between observed data and model predictions (interferometry_fourier) 
    for a Disk model with parameters b. =#
    function f(b, donnees,uv)
        x = b[1]u"mas"
        y = b[2]u"mas"
        R = b[3]u"mas"
        sum(abs2, interferometry_fourier(Disk([x, y],[R]), uv) .- donnees)
    end
    # Loop over each optimization method
    #= Iterates over the range of initial values for x and y. For each combination, it initializes x0 and performs optimization
       using optimize with the LBFGS() method =#
    for i in range
        for j in range
            for k in 0.0:step(range)*2:last(range)
                x0 = [i,j,k]
                # Optimization using LBFGS method without autodiff
                opt_result = optimize(b -> f(b, data,tabuv), x0,LBFGS())

                # Retrieve results
                global_min = opt_result.minimum
                minimizer = opt_result.minimizer

                # Uses isapprox(global_min, 0.0; atol=0.5) to determine if the optimization successfully minimized the objective function
                if isapprox(global_min,0.0; atol=0.5)
                    return Disk(minimizer[1].*u"mas",minimizer[2].*u"mas",minimizer[3].*u"mas")
                end
            end
        end
    end    
    error("Optimization did not find a satisfactory solution.")
end
